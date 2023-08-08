"""
This module contains the classes and tools to handle a GRB object, i.e.
a collection of lightcurve bins for each of which is given an energy spectrum,
an explosion (trigger) time and physical characteristics.

Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""

import os
import sys
import pickle
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import cm

import astropy
import astropy.units         as u
from   astropy.table         import Table, QTable
from   astropy.io            import fits
from   astropy.time          import Time
from   astropy.coordinates   import SkyCoord
from   astropy.coordinates   import AltAz
from   astropy.visualization import quantity_support

import observatory as obs
from visibility import Visibility

from niceprint import warning, failure
from niceplot import single_legend
from utilities import get_filename

from gammapy.modeling.models import Models
from gammapy.modeling.models import PointSpatialModel, SkyModel
from gammapy.modeling.models import EBLAbsorptionNormSpectralModel
from gammapy.modeling.models import TemplateSpectralModel

# Transform warnings into errors - useful to find who is guilty !
import warnings
#warnings.filterwarnings('error')
warnings.filterwarnings('ignore')

__all__ = ["GammaRayBurst"]

###############################################################################
class GammaRayBurst(object):
    """
    A class to store GRB properties.

    The GRB information is read from files.
    A GRB is composed of several parameters among which, a name, a redshift,
    a position in the sky, and measurement points in time at which
    an energy spectrum is given.

    Each energy spectrum is given as a series of flux measurements but is then
    stored as an interpolated spectrum* (i.e. the energy bin width has a modest
    importance).

    The GRB properties like `z` and `Eiso` connect the afterglow and prompt
    emissions together.

    The original energy bins can be limited by a maximal maximal energy
    (Limitation in time can be obtained from the visibility computation).

    The spectra are usually modified for an EBL absorption, or the absorbed
    flux is given in the file in some cases.

    The afterglow flux can be adjusted by a multiplicative factor for
    various tests.

    If requested, the associated prompt spectrum is read and added to the
    afterglow spectra to give the total contribution.

    **Note:**

    Following a question on the `Slack gammapy` channel on November 27 :sup:`th` ,
    and the answer by Axel Donath:

    The following statement later in the code gave an error
    (:code:`dlist_onoff` is a collection of `Dataset`)
    :code:`dlist_onoff.write(datapath,prefix="cls",overwrite=True)`
    gives:

     ..  code-block:: python

        ..\\gammapy\\modeling\\models\\spectral.py",
        line 989, in to_dict "data": self.energy.data.tolist(),
        NotImplementedError: memoryview: unsupported format >f


    This error comes from the fact that the energy list as to be
    explicitely passed as a float as done below: :code:`cls.Eval.astype(float)`
    (A `Quantity` is passed as requested but the underlying numpy
    dtype is not supported by :code:`energy.data.tolist()`)

    """
    ignore = ["id","filename", "Eval","tval", "tstart", "tstop", "fluxval",
              "spec_afterglow", "prompt", "id90", "E_prompt", "flux_prompt",
              "spec_prompt", "models", "vis"]
    """ Ignore these members when dumping the class out."""

    ###------------------------------------------------------------------------
    def __init__(self):
        """
        This initializes a default GRB.

        Returns
        -------
        None.

        """

        self.id       = 'dummy'
        self.filename = None # File with the spectral and lightcurve data
        self.z        = 0 # Redshift
        self.eblmodel = None

        # GRB properties - Dummy default values
        self.radec     = SkyCoord(ra=100*u.deg, dec= -15*u.deg, frame="icrs")
        self.Eiso      = 0.*u.erg
        self.Liso      = 0.*u.Unit("erg/s")
        self.Epeak     = 0.*u.keV
        self.t90       = 0.*u.s
        self.G0H       = 0.*u.dimensionless_unscaled
        self.G0W       = 0.*u.dimensionless_unscaled
        self.Fpeak     = 0.*u.erg
        self.Fpeak_GBM = 0.*u.erg
        self.gamle     = 0.*u.dimensionless_unscaled
        self.gamhe     = 0.*u.dimensionless_unscaled

        # GRB explosion time
        self.t_trig   = Time('2020-01-01T02:00:00', format="isot", scale='utc')

        #------------
        ### Afterglow
        #------------

        # Afterglow Flux table - Spectra at a series of points in time
        self.Eval           = [0]*u.GeV # Energy values
        self.tval           = [0]*u.s   # Time values
        self.tstart         = 0*u.s     # Start time
        self.tstop          = 0*u.s     # Stop time

        # Raw flux arrays
        self.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")

        # Afterglow interpolated arrays
        self.spec_afterglow = []

        #-------------------
        ### Prompt component
        #-------------------
        # Default : no prompt
        self.prompt      = False

        # Afterglow Slice id at which the prompt stops
        self.id90        = -1

        # Prompt energy bins
        self.E_prompt    = [0]*u.GeV

        # Prompt raw flux value
        self.flux_prompt = 1*u.Unit("1 / (cm2 GeV s)")

        # One Interpolated, non attenuated E-spectrum
        self.spec_prompt = None

        #-------------
        ### Total flux
        #-------------
        # If the prompt flux exist, this the sum of both flux until the
        # last time bin of the afterglow flux below t90 of the prompt.

        # Gammapy Skymodel models (one per t slice)
        self.models = []

        #-------------
        ### Visibility
        #-------------
        # Visibility - dictionnary with site entries
        self.vis  = dict()

    ###------------------------------------------------------------------------
    @classmethod
    def from_fits(cls,
                  filename,
                  prompt  = None,
                  ebl     = None,
                  Emax    = None,
                  dt      = 0.0,
                  magnify = 1,
                  vis_from_file = False,
                  debug = False):

        """
        Read the GRB data from a `fits` file.
        Fluxes are given for a series of (t,E) values
        So far no flux exists beyond the last point (no extrapolation).

        The spectrum is stored as a table, :obj:`TemplateSpectralModel`, that takes
        as an input a series of flux values as a function of the energy.
        The model will return values interpolated in log-space with
        :func:`scipy.interpolate.interp1d`, returning zero for energies
        outside of the limits of the provided energy array.

        In this function, the trigger time (time of the GRB explosion), is
        expected to come from the input file, with the possibility to add a
        constant time shift in days, or is set randomly over the default year.
        Alternatively, the dates can be overwritten by dates stored in
        an external file in the `SoHAPPy.py` main process using
        :func:`grb.GammaRayBurst.set_visibility` function.

        Parameters
        ----------
        cls : GammaRayBurst class
            GammaRayBurst class instance.
        filename : Path
            GRB input file name.
        prompt: Boolean
            If True, read information from the associated prompt component.
        ebl : String, optional
            The EBL absorption model considered. If ebl is `built-in`, uses
            the absorbed specrum from the data if available. Can be `None`.
        Emax : Astropy Quantity, optionnal
           Maximal energy considered in the data. The default is `None`.
        dt: float, optionnal
            A number of Julian days to be added to the present source trigger
            time. The default is 0.0.
        magnify : float, optional
            Flux multiplicative factor to the afterglow model flux for tests.
            Default is 1.
        vis_from_file: Boolean, optional
            If True the visibility is read from the file.
            The default is `False`.
        debug: Boolean, optional
            Debugging flag. The default is False.


        Returns
        -------
        A GammaRayBurst instance.

        """

        ### -----------------------------------------------------
        ### Open file, get header, keys, and data fill the class members
        ### -----------------------------------------------------
        hdul   = fits.open(get_filename(filename))
        hdr    = hdul[0].header
        keys_0 = list(hdul[0].header.keys())

        cls = GammaRayBurst() # Default constructor

        cls.filename = filename
        cls.z        = hdr['Z']
        cls.eblmodel = ebl

        # Get identifier by removing all extensions
        cls.id = str(filename.name).rstrip(''.join(filename.suffixes))

        cls.radec    = SkyCoord(ra = hdr['RA']*u.deg, dec = hdr['DEC']*u.deg,
                                frame="icrs")

        cls.Eiso     = hdr['EISO']*u.erg

        if "LISO" in keys_0:
            cls.Liso = hdr["LISO"]*u.Unit("erg/s") # New large prod. files
        else:
            cls.Liso = 0*u.Unit("erg/s")

        cls.Epeak    = hdr['EPEAK']*u.keV
        cls.t90      = hdr['Duration']*u.s
        cls.G0H      = hdr['G0H']
        if "G0W" in keys_0: # Not in SHORTFITS
            cls.G0W      = hdr['G0W']

        if "PHFLUX" in keys_0:
            cls.Fpeak = hdr['PHFLUX']*u.Unit("cm-2.s-1")
        elif "PHFX" in keys_0:
            cls.Fpeak = hdr['PHFX']*u.Unit("cm-2.s-1")

        cls.gamle = hdr['LOWSP']
        cls.gamhe = hdr['HIGHSP']

        if "PHFX_GBM" in keys_0:
            cls.Fpeak_GBM = hdr['PHFX_GBM']*u.Unit("cm-2.s-1")
        else:
            cls.Fpeak_GBM = 0*u.Unit("cm-2.s-1")

        ###--------------------------
        ### GRB trigger time
        ###--------------------------
        # If the input file does not contain a trigger date, then set a date
        # by default. It is very likely that this source belongs to a set
        # for which a list of explosion (trigger) times should be generated.

        if "GRBJD" in keys_0: # Not in SHORTFITS
            cls.t_trig = Time(hdr['GRBJD']*u.day + dt*u.day,
                              format="jd",scale="utc")
        elif "GRBTIME" in keys_0: # Add start date
            cls.t_trig = Time(hdr['GRBTIME'] + dt*u.day,
                              format="jd",scale="utc")
        else:
            warning(f"{__name__:}.py: Trigger time absent from file, using random value")
            import random
            cls.t_trig += random.uniform(0, 365)*u.day

        ###--------------------------
        ### Time intervals - so far common to afterglow and prompt
        ###--------------------------
        tval  = QTable.read(hdul["TIMES (AFTERGLOW)"])

        # Temporary : in short GRB fits file, unit is omitted
        # The flux value is given on an interval ("col0", "col1°°.
        # The default is to consider that the flux is valid at the end of the
        # interval.

        if isinstance(tval[0][0], astropy.units.quantity.Quantity):
            cls.tval = np.array(tval[tval.colnames[0]].value)*tval[0][0].unit
        else: # In SHORT GRB, the time bin is given, [t1, t2].
            cls.tval = np.array(tval["col1"])*u.s

        cls.tstart   = cls.t_trig
        cls.tstop    = cls.t_trig + cls.tval[-1]

        ###--------------------------
        ### Afterglow Energies - Limited to Emax if defined
        ###--------------------------
        tab_key = "Energies (afterglow)"
        col_key = Table.read(hdul[tab_key]).colnames[0]
        cls.Eval  = Table.read(hdul[tab_key])[col_key].quantity
        cls.Eval = np.array(cls.Eval)*cls.Eval[0].unit

        if Emax is not None and Emax<=cls.Eval[-1]:
            warning(" Data up to {:5.3f} restricted to {:5.3f}"
                    .format(cls.Eval[-1],Emax))
            cls.Eval = cls.Eval[cls.Eval<= Emax]

        ###--------------------------
        ### Afterglow flux
        ###--------------------------

        # Get flux, possibly already absorbed
        if ebl == "in-file": # Read default absorbed model
            flux = QTable.read(hdul["EBL-ABS. SPECTRA (AFTERGLOW)"])
        else:
            flux = QTable.read(hdul["SPECTRA (AFTERGLOW)"])

        flux_unit    = u.Unit(flux.meta["UNITS"])
        if str(flux_unit).find("ph") > -1:
            flux_unit = flux_unit/u.Unit("ph") # Removes ph

        # Store the flux. Note the transposition
        itmax = len(cls.tval)-1
        jEmax = len(cls.Eval) - 1
        cls.fluxval = np.zeros( (itmax+1,jEmax+1) )*flux_unit

        for i in range(0,itmax+1):
            for j in range(0,jEmax+1):
                cls.fluxval[i][j] = magnify* flux[j][i]*flux_unit # transp!

        # Build time series of interpolated spectra - limited to dtmax
        for i in range(len(cls.tval)):
            glow = TemplateSpectralModel(energy = cls.Eval.astype(float),
                                         values = cls.fluxval[i],
                                         interp_kwargs={"values_scale": "log"})
            cls.spec_afterglow.append(glow)

        ###--------------------------
        ### Prompt - a unique energy spectrum
        ###--------------------------

        # Get the prompt if potentially visible and if requested
        cls.prompt = False # No prompt component was found
        if prompt is not None:
            if Path(prompt).is_dir():
                # if cls.vis["North"].vis_prompt or cls.vis["South"].vis_prompt:
                # Deduce prompt folder from GRB path name
                cls.spec_prompt = cls.get_prompt(folder = prompt,
                                                 debug  = debug)
                if cls.spec_prompt is not None:
                    cls.prompt = True
            else:
                sys.exit(" {} is not a valid data folder :".format(prompt))

        ###--------------------------
        ### Total attenuated spectra
        ###--------------------------
        for i,t in enumerate(cls.tval):
            spec_tot = cls.spec_afterglow[i]
            if cls.prompt:
                if i<=cls.id90: spec_tot += cls.spec_prompt[i]
            m = SkyModel(spectral_model= cls.EBLabsorbed(spec_tot,ebl),
                         spatial_model = PointSpatialModel(lon_0=cls.radec.ra,
                                                           lat_0=cls.radec.dec,
                                                           frame="icrs"),
                             name="model_"+str(i))
            cls.models.append(m)

        ### Close fits file
        hdul.close()

        return cls

    ###------------------------------------------------------------------------
    @classmethod
    def historical_from_yaml(cls, item, ebl=None, magnify=1):

        """
        Read the characteristics of a parameterised GRB from a `yaml` file.
        It is assumed that the energy and time decays are not correlated and
        follow two independent power laws.
        The function associates spectra to a list of time intervals
        in order to comply with the more general case in this class.
        Note that the spectra are considered to be from the afterglow only.

        Returns
        -------
        A `GammaRayBurst` instance.

        """

        cls = GammaRayBurst() # This calls the constructor

        cls.filename  = Path(Path(__file__).absolute().parent,
                         "data/historical/GRB_"+item+".yml")

        with open(cls.filename) as f:
            import yaml
            from yaml.loader import SafeLoader
            data = yaml.load(f, Loader=SafeLoader)


        cls.z        = data["z"]
        cls.eblmodel = ebl
        cls.id       = data["name"]

        cls.radec  = SkyCoord(data["ra"], data["dec"], frame='icrs')
        cls.Eiso   = u.Quantity(data["Eiso"])
        # self.Liso  = 0.*u.Unit("erg/s")

        cls.Epeak  = u.Quantity(data['Epeak'])
        cls.t90    = u.Quantity(data['t90'])
        cls.G0H    = data['G0H']
        cls.G0W    = data['G0W']
        cls.Fpeak  = u.Quantity(data['Fluxpeak'])
        # self.Fpeak_GBM = 0.*u.erg
        cls.gamle  = data['gamma_le']
        cls.gamhe  = data['gamma_he']

        ###--------------------------
        ### GRB trigger time
        ###--------------------------
        cls.t_trig = Time(data["t_trig"],format="datetime",scale="utc")

        ###--------------------------
        ### Time intervals
        ###--------------------------
        tmin     = max(u.Quantity(data["tmin"]),1*u.s) # Cannot be zero
        tmax     = u.Quantity(data["tmax"])
        ntbin    =  data["ntbin"]

        if ntbin != 1:
            cls.tval = np.logspace(np.log10(tmin.to(u.s).value),
                                   np.log10(tmax.to(u.s).value),
                                   ntbin)*u.s
        else: # A single time window
            cls.tval = np.array([tmin.value,tmax.to(tmin.unit).value])*tmin.unit

        cls.tstart   = cls.t_trig
        cls.tstop    = cls.t_trig + cls.tval[-1]

        ###--------------------------
        ### Afterglow Energies - Limited to Emin, Emax
        ###--------------------------
        Emin     = u.Quantity(data["Emin"])
        Emax     = u.Quantity(data["Emax"])

        cls.Eval = np.asarray([Emin.value,
                               Emax.to(Emin.unit).value])*Emin.unit

        ###--------------------------
        ### Afterglow flux
        ###--------------------------

        flux_unit = u.Unit("1/(cm2 GeV s)")
        cls.fluxval = np.zeros( (len(cls.tval),len(cls.Eval)) )*flux_unit

        for i,t in enumerate(cls.tval):
            for j,E in enumerate(cls.Eval):
                dnde = (u.Quantity(data["K"])*(E/data["E0"])**-data["gamma"]
                                 *(t/data["t0"])**-data["beta"])
                #print(i,j,dnde)
                cls.fluxval[i][j] = magnify* dnde.to(flux_unit)

        ###--------------------------
        ### Prompt - No prompt foreseen in this case
        ###--------------------------
        cls.prompt = False

        ###--------------------------
        ### Total attenuated spectra
        ###--------------------------
        # See note in from_fits for explanatins on some parameters
        for i,t in enumerate(cls.tval):

            tab = TemplateSpectralModel(energy = cls.Eval.astype(float),
                                        values = cls.fluxval[i],
                                        interp_kwargs={"values_scale": "log"})

            cls.spec_afterglow.append(tab) # Needed for display


            m = SkyModel(spectral_model= cls.EBLabsorbed(tab),
                         spatial_model = PointSpatialModel(lon_0=cls.radec.ra,
                                                           lat_0=cls.radec.dec,
                                                           frame="icrs"),
                         name="model_"+str(i))

            cls.models.append(m)

        return cls

    ###------------------------------------------------------------------------
    @classmethod
    def prompt(cls, filename, glowname= None, ebl= None,
                             z=0*u.dimensionless_unscaled, magnify=1):
        """
        This function is for tests using time-resolved spectra.
        It reads prompt data from a file and associate the prompt to the
        afterglow if requested (or keep the default from the constructor
        otherwise, with possibility to supersede the redshift and the flux
        multiplicative factor). It does not add the two spectra as this would
        require to work on the energy and time binnings that are different and
        this is not implemented yet.

        Parameters
        ----------
        filename : Path
            Prompt file path.
        glowname : Path, optional
            Afterglow file path. The default is None.
        ebl : string, optional
            Name of the EBL model. The default is None.
        z : Quantity, optional
            redshift. The default is 0*u.dimensionless_unscaled.
        magnify : float, optional
            Multiplicative factor of the flux. The default is 1.

        Returns
        -------
        GammaRayBurst
            New instance.

        """

        if not filename.exists():
            sys.exit(f"File {filename:} not found")

        # Read prompt data -supersede defaults
        cls = GammaRayBurst()
        cls. z       = z
        cls.magnifiy = magnify
        cls.eblmodel = ebl

        # Copy missing data from the corresponding afterglow
        if glowname is not None:
            glow = GammaRayBurst.from_fits(glowname)
            # print(glow)
            cls.z         = glow.z
            cls.radec     = glow.radec
            cls.Eiso      = glow.Eiso
            cls.Liso      = glow.Liso
            cls.Epeak     = glow.Epeak
            cls.Fpeak     = glow.Fpeak
            cls.t90       = glow.t90
            cls.G0H       = glow.G0H
            cls.G0W       = glow.G0W
            cls.gamle     = glow.gamle
            cls.gamhe     = glow.gamhe
            cls.Fpeak_GBM = glow.Fpeak_GBM
            cls.t_trig    = glow.t_trig

        # Reads prompt data
        hdul = fits.open(filename)
        cls.id = Path(filename.name).stem

        cls.Eval     = Table.read(hdul,hdu=1)["energy"].quantity*u.TeV

        cls.tval     = Table.read(hdul,hdu=2)["time"].quantity*u.s
        cls.tstart   = cls.t_trig
        cls.tstop    = cls.t_trig + cls.tval[-1]

        flux         = Table.read(hdul,hdu=3)
        flux_unit = u.Unit("1/(cm2 TeV s)")
        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Note the transposition from flux to fluxval
        cls.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                f = flux[j][i]
                if f>1:
                    # print(i,j,f)
                    f =0 # Correct a bug in event #172 - to be removed
                cls.fluxval[i][j] = magnify*f*flux_unit # transp!

        for i in range(len(cls.tval)):
            # Note that TableModel makes an interpolation
            prpt = TemplateSpectralModel(energy  = cls.Eval.astype(float),
                                         values = cls.fluxval[i],
                                         interp_kwargs={"values_scale": "log"})
            # Use the afterglow place holder to store the prompt
            cls.spec_afterglow.append(prpt)

            m = SkyModel(spectral_model= cls.EBLabsorbed(prpt,ebl),
                         spatial_model = PointSpatialModel(lon_0=cls.radec.ra,
                                                           lat_0=cls.radec.dec,
                                                           frame="icrs"),
                                                            name="model_"+str(i))
            cls.models.append(m)

        # Not adding the avraged on time prompt
        cls.prompt = False

        hdul.close()

        return cls

    ###------------------------------------------------------------------------
    def EBLabsorbed(self, tab, debug=False):

        """
        Returns the EBL-absorbed model of the current instance.
        Absorption data are either obtained from the Gammapy datasets or from
        external proprietary files. Data are unchanged if the GRB has
        already attenuated spectra or if no absorption is considered.

        Parameters
        ----------
        tab : Gammapy TemplateSpectralModel
            A flux versus energy as an interpolated table.
        model : String
            An EBL model name among those available.
        debug : Boolean, optional
            If True let's talk a bit. The default is False.

        Returns
        -------
        attflux : TemplateSpectralModel
            An attenuated flux versus energy as an interpolated table.

        """

        attflux = tab # Initialise to the unabsorbed flux
        if self.eblmodel is None or self.eblmodel == 'in-file':
            return attflux

        if self.eblmodel != "gilmore":
            eblabs = EBLAbsorptionNormSpectralModel.read_builtin(self.eblmodel,
                                                              redshift=self.z)
            attflux = tab*eblabs
        else:
            from ebl import EBL_from_file
            eblabs = EBL_from_file("data/ebl/others/ebl_gilmore12-10GeV.dat")
            # Change flux for absorption, recreate model
            # I di not fine a clever way
            attenuated = []
            for E in self.Eval:
                attenuated.append(tab(E).value*eblabs(E,self.z))
            attenuated *= tab(np.mean(self.Eval)).unit
            attflux = TemplateSpectralModel(energy = self.Eval.astype(float),
                                            values = attenuated,
                                            interp_kwargs={"values_scale": "log"})

        return attflux

    ###------------------------------------------------------------------------
    def set_visibility(self, item, loc, info=None, n_night=None, n_skip=None,
                             status="", dbg=False):
        """
        Attach a visibility to a GRB instance.
        Either recompute it if a keyword has been given and a dictionnary
        retrieved or read it from the specified folder or dictionnary.

        * If a dictionnary is given, it can be:

            1. From the :obj:`visibility.yaml` file - the visibility is
               computed on the fly.
            2. Not from the :obj:`visibility.yaml` file, the dictionnary
               contains a computed visibility from a `json` file.

        * The keyword is among those:

            * **built-in**, the visibility is read from the fits file.
            * **forced**, the visibility is build assuming one infinite
              night.
            * **permanent**, the visibility is permanent, i.e. the night is
              infinite and the GRB is above the horizon. In that case it
              can be useful to force the zenith angle to be fixed at a
              certain value.

        """

        ### Update the default instance
        self.vis[loc] = Visibility(pos    = self.radec,
                                   site   = obs.xyz["CTA"][loc],
                                   window = [self.tstart, self.tstop],
                                   name   = self.id+"_"+loc,
                                   status = status)

        # At this stage the visibility is maximal
        if info == "permanent":
            self.status = "permanent"
            return self

        if isinstance(info,dict): # A dictionnary not a keyword

            if "altmin" in info.keys(): # Compute from dictionnary

                # Supersede max. number of nights and nights to be skipped
                if n_night is not None:
                    info["depth"] = n_night
                if n_skip is not None:
                    info["skip"] = n_skip

                # Compute from dictionnary elements
                self.vis[loc] = self.vis[loc].compute(param = info,debug = dbg)

            else: # A dictionnary of all visibilities (from a .json file)
                self.vis[loc]  = Visibility.from_dict(info[str(item)+"_"+loc])

        # Not a dictionnary
        elif info == "built-in":
            # sys.exit("{}.py : built-in vis. reading to be reimplemented"
            #          .format(__name__))
            # infolder is not passed
            self.vis[loc]  = Visibility.from_fits(self, loc=loc)

        else: # Special or from disk

            if info == "forced": # Infinitite nights
                self.vis[loc] = self.vis[loc].force_night()

            else: # Binary - Obsolete - Archived bin files might be not valid
                import pickle
                with open(Path(info,self.id+"_"+loc+"_vis.bin"),"r") as f:
                    self.vis[loc] =  pickle.load(f)

###------------------------------------------------------------------------
    def altaz(self,loc="",dt=0*u.s):
        """
        Get altitude azimuth for the GRB at a given site at GRB time t (s)

        Parameters
        ----------
        location : string
            Either `North` or `South`.
        dt : Quantity (time)
            The time elapsed since the trigger time

        Returns
        -------
        altaz : astropy.coordinates.AltAz
            The GRB poistion in the sky at the given time and site.

        """
        if type(dt) is not astropy.units.quantity.Quantity:
            sys.exit("Time has to be a quantity (weird behaviour otherwise")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            altaz = self.radec.transform_to(AltAz(obstime  = dt + self.t_trig,
                                                  location = self.vis[loc].site))
        return altaz

    ###------------------------------------------------------------------------
    def get_prompt(self, folder = None, strict=False, debug = False):
        """
        Get the time averaged prompt component associated to the afterglow.
        The prompt spectra are produced so that they correspond to the values
        of the physical parameters from the population (Lorentz Factor,
        peak energy, photon flux, and redshift).
        A single set of spectra is given spanning over the `T90` GRB duration.
        The flux is an "average" prompt spectrum over `T90`.

        Parameters
        ----------
        folder : String, optional
            Where to find the prompt data below the 'lightcurves/prompt'
            folder. The default is None.
        strict : Boolean, optional
            If True, requires that the density profile agrees between the
            prompt and the afterglow.
        debug : Boolean, optional
            If True, let's talk a bit. The default is False.

        Returns
        -------
        List of Models
            List of `Gammapy` models for the Prompt component.

        """

        # Find back the source id number from the id (assumes only the id has digits)
        gid = ""
        for s in self.id:
            if s.isdigit():
                gid += s

        if strict:
            if gid in [6, 30 ,191]:
                failure(f" GRB prompt {gid:} simulated with a Wind profile")
                failure(" It cannot be associated to the current afterglow !")
                return -1, None

        ###----------------------------
        ### Compute weihts to "Mask" the time interval beyond t90
        ###----------------------------

        # Find the closest time bin at t90 in the Afterglow
        # id90 is the index of the time following the t90 value
        # t90 is therefore between the index id90-1 and id90
        self.id90 = np.where(self.tval >= self.t90)[0][0]
        if self.id90 == 0:
            warning(f"{self.id:} : t90 = {self.t90:} "\
                    "is before the first afterglow time bin")
            return None

        # The given flux is per unit of time.
        # In order to take into account that the last time bin is not completely
        # covered, it is corrected by the fraction in time in this bin.
        dtprompt = self.t90 - self.tval[self.id90-1] # Prompt is active
        dtslice  = self.tval[self.id90]-self.tval[self.id90-1] # Total slice duration
        fraction = dtprompt/dtslice

        # All slices have weight 1 before the slice containing t90, the
        # fraction above for the slice containing t90, and zero beyond
        weight = np.concatenate(  (np.ones(self.id90),
                                  [fraction],
                                   np.zeros(len(self.tval)-self.id90-1)))

        ###----------------------------
        ### Read prompt data
        ###----------------------------
        # Open file - read the lines
        filename = Path(folder,gid+"-spec.dat")
        if not filename.exists():
            failure("{} does no exist".format(filename))
            return None

        file  = open(filename,"r")
        lines = file.readlines()

        # get Gamma and z
        data         = lines[0].split()
        gamma_prompt = float(data[2])
        redshift     = float(data[5])

        # Consistency check - Zeljka data are rounded at 2 digits
        if np.abs(self.z - redshift)>0.01 or np.abs(self.G0H - gamma_prompt)>0.01:
            failure(f" {self.id:10s}: "\
                    f"Afterglow / prompt : z= {self.z:4.2f} / {redshift:4.2f}"\
                    f"  G0H/gamma= {self.G0H:5.2f} / {gamma_prompt:5.2f}"\
                    f" (G0W= {self.G0W:5.2f})")

        # Get units - omit "ph" in flux unit
        data     = lines[1].split()
        unit_e   = data[0][1:-1]
        unit_flx = ''.join([ x + " " for x in data[2:]])[1:-2]

        # Get flux points - assign units to data, use afterglow units
        data = lines
        Eval = []
        fluxval = []
        for line in lines[4:]:
            data = line.split()
            Eval.append(float(data[0]))
            fluxval.append(float(data[1]))

        self.E_prompt    = Eval*u.Unit(unit_e)
        self.flux_prompt = fluxval*u.Unit(unit_flx)

        if unit_e != self.Eval.unit:
            if debug:
                warning(f"Converting energy units from {unit_e:}"\
                        f" to {self.Eval[0].unit:}")
            self.E_prompt = self.E_prompt.to(self.Eval[0].unit)

        if unit_flx != self.fluxval[0].unit:
            if debug:
                warning("Converting flux units from {} to {}"
                        .format(unit_flx, self.fluxval[0][0].unit))
            self.flux_prompt = self.flux_prompt.to(self.fluxval[0][0].unit)

        ###----------------------------
        ### Create a list of weighted models for non-zero weights
        ###----------------------------
        models = []
        for i in range(self.id90+1):
            flux_w = self.flux_prompt*weight[i]
            models.append(TemplateSpectralModel(energy = self.E_prompt,
                                        values = flux_w,
                                        interp_kwargs={"values_scale": "log"}))

        if debug:
            print("Prompt associated to ",self.id)
            print("Gamma = ",gamma_prompt," redshift =",redshift)
            print("  {:8} : {:10}".format(unit_e, unit_flx))
            for E, flux in zip(self.E_prompt,self.flux_prompt):
                print("  {:8.2f} : {:10.2e} "
                      .format(E.value,flux.value))
            islice=0
            print(" {:>8} - {:>5} ".format("Time","Weight"),end="")
            for t, w in zip(self.tval, weight):
                print(" {:8.2f} - {:5.2f} ".format(t,w),end="")
                print("*") if islice== self.id90 else print("")
                islice+=1

            print("Prompt is between: ",self.tval[self.id90-1],self.tval[self.id90])

        return models

    ###------------------------------------------------------------------------
    ### INPUT/OUPUT
    ###------------------------------------------------------------------------
    def __str__(self):
        """
        Printout the GRB properties (Not the visibilities).
        """

        txt=""
        if self.filename is not None:
            txt += '  Read from      : {}\n'.format(self.filename)
        txt += '  RA, DEC        : {:6.2f} {:6.2f}\n' \
            .format(self.radec.ra.value, self.radec.dec.value)
        txt += '  Redshift       : {:6.2f}\n'.format(self.z)
        txt += '  EBL model      : "{:}"\n'.format(self.eblmodel)
        txt += '  Eiso           : {:6.2e}\n'.format(self.Eiso)
        txt += '  Epeak          : {:6.2f}\n'.format(self.Epeak)
        txt += '  t90            : {:6.2f}\n'.format(self.t90)
        txt += '  G0H / G0W      : {:6.2f} / {:6.2f}\n' \
        .format(self.G0H,self.G0W)
        txt += '  Flux peak      : {:6.2f}\n'.format(self.Fpeak)
        txt += '  Flux peak (GBM): {:6.2f}\n'.format(self.Fpeak_GBM)
        txt += '  gamma LE / HE  : {:6.2f} / {:6.2f}\n' \
        .format(self.gamle,self.gamhe)

        txt += '  t_trig         : {}\n'.format(Time(self.t_trig,format="iso"))
        txt += '  Duration       : {:6.2f} {:10.2f}\n' \
            .format( (self.tval[-1]-self.tval[0]).to(u.d),
                     (self.tval[-1]-self.tval[0]).to(u.s))
        # txt += '  Bins : E, t, dt    : {:4d} {:4d} {:4d} \n' \
        # .format(len(self.Eval), len(self.tval), len(self.time_interval))
        txt += '  Bins : E, t    : {:4d} {:4d} \n' \
        .format(len(self.Eval), len(self.tval))
        if self.prompt:
            txt += ' Prompt component  : \n'
            txt += '  Up to slice    : {:3d}\n'.format(self.id90)
            txt += "  Bins : E       : {:3d}\n".format(len(self.E_prompt))
        else:
            txt += ' Prompt component not considered\n'
        if self.vis is None:
            txt += ' Visibility not available\n'

        return  txt

    ###------------------------------------------------------------------------
    def write_to_bin(self, folder, debug=True):
        """
        Write current instance to a binary file for further
        reuse.

        Parameters
        ----------
        folder : String
            Folder name.
        debug : Boolean, optional
            If True, talks  bit. The default is True.

        Returns
        -------
        None.

        """

        filename = Path(folder,self.id+".bin")
        outfile  = open(filename,"wb")
        pickle.dump(self,outfile)
        outfile.close()

        if debug:
            print(" GRB saved to : {}".format(filename))

    ###------------------------------------------------------------------------
    def model_to_yaml(self, output="yaml"):
        """
        Dump the time values for each energy specturm into a text file.
        Dump all spectra associated to the source into `yaml` files, tar and
        compresss the files, delete the originals.

        Parameters
        ----------
        output : string, optional
            Output folder name. The default is "yaml".

        Returns
        -------
        None.

        """

        # Create dedicated folder
        folder    =Path(output,grb.name+"/")
        if not folder.is_dir():
            import os
            os.makedirs(folder)

        print(" Now dumping ",grb.name,"to ",folder)

        # Dump measurement times
        name= Path(folder,grb.name+"_times.txt")
        out = open(name,"w")
        for i,t in enumerate(grb.tval):
            print("{:3d} {:>8.2f}".format(i,t),file=out)
        out.close()

        # Dump spectra and add to tar file
        from gammapy.modeling.models import Models
        import tarfile
        tar = tarfile.open(Path(folder,grb.name+".tar.gz"), "w:gz")

        for i, spec in enumerate(grb.models):
            specname     = "spec_{:3s}.yaml".format(str(i+1).zfill(3))
            specfilename = Path(folder,specname)

            # Dump data in yaml format
            out          = open(specfilename,"w")
            print(Models([spec]).to_yaml(),file=out)
            out.close()

            # Put spectrum in tar file
            tar.add(specfilename,arcname=specname)

            # Delete file
            Path.unlink(specfilename)

        tar.close()

    ###------------------------------------------------------------------------
    def models_to_fits(self, output="model2fits", debug=False):

        """
        Function created to provide input to the Data Challenge, 2023.
        Store the spectra as yaml records inside an astropy.table Table
        associated to the valid time interval.

        Parameters
        ----------
        output : string, optional
            The output folder name. The default is "fits".
        debug : Boolean, optional
            If True, say something. The default is False.

        Returns
        -------
        Path.
            Output file Path.

        """

        # Create dedicated folder
        folder = Path(output)
        if not folder.is_dir():
            os.makedirs(folder)

        # Generate all yaml models and get the maximum length
        yaml_list = [ Models([spec]).to_yaml() for spec in self.models]
        maxlength = max([len(m) for m in yaml_list])
        if debug:
            print("   Max. yaml model length is ",maxlength)

        # Fits filename
        output_name = Path(folder,"DC_"+self.id+".fits")
        print(" Now dumping ",grb.name,"to ",output_name," as a fits table")

        # Create an astropy table with 3 columns for spectra (same number of rows)
        table = Table()
        table["time_start"] = np.append([0e0],grb.tval[:-1].value)*grb.tval.unit
        table["time_stop"]  = grb.tval
        table["mod"] = maxlength*" "

        for i, mod in enumerate(yaml_list):
            table["mod"][i] = mod

        # Add this stage, the data can be dumped into a fits file with
        # table.write(output_name, overwrite=True)
        # However, in order to modify the header we add a step
        hdu_spectra = fits.BinTableHDU(data=table,name="spectra")

        # Add explosion time in the header
        header = hdu_spectra.header
        header["trigger"]= self.t_trig.isot

        # Add visibility
        table = Table()
        if len(self.vis["North"].t_true[0]) !=0:
            table["tstart"] = [ t[0].isot   for t in self.vis["North"].t_true]
            table["tstop"]  = [ t[1].isot   for t in self.vis["North"].t_true]
        else:
            table["tstart"]  = [-1]
            table["tstop"]   = [-1]

        hdu_visN = fits.BinTableHDU(data=table,name="vis_N")

        table = Table()
        if len(self.vis["South"].t_true[0]) !=0:
            table["tstart"] = [ t[0].isot   for t in self.vis["South"].t_true]
            table["tstop"]  = [ t[1].isot   for t in self.vis["South"].t_true]
        else:
            table["tstart"]  = [-1]
            table["tstop"]   = [-1]
        hdu_visS = fits.BinTableHDU(data=table,name="vis_S")

        # Create the hdulist, write to output
        hdul = fits.HDUList()
        hdul.append(hdu_spectra)
        hdul.append(hdu_visN)
        hdul.append(hdu_visS)

        # Write to file
        hdul.writeto(output_name,overwrite=True)
        hdul.close()

        return output_name

    ###------------------------------------------------------------------------
    @staticmethod
    def check_models_from_fits(filename):
        """
        To be checked and re-implemented or deleted.

        Parameters
        ----------
        filename : string
            Input file name.

        Returns
        -------
        None.

        """

        hdul = fits.open(filename)
        print(hdul.info())

        # Get trigger time
        hdu = hdul[1]
        print("Trigger/explosion time :",hdu.header["TRIGGER"])

        # Get spectral data
        data = hdu.data
        print(data.columns)
        data["time_start"]
        data["time_stop"]
        for m in data["mod"]:
            print(Models.from_yaml(m))

        # Get visibility periods
        print("Visibility windows in North")
        print(hdul["VIS_N"].data["tstart"])
        print(hdul["VIS_N"].data["tstop"])

        print("Visibility windows in South")
        print(hdul["VIS_S"].data["tstart"])
        print(hdul["VIS_S"].data["tstop"])

        hdul.close()

    ###------------------------------------------------------------------------
    ### PLOTS
    ###------------------------------------------------------------------------

    def plot(self, pdf=None):
        """
        A series of plot of the present instance. The figures can also be
        printed to a pdf file.

        Parameters
        ----------
        pdf : PdfPages instance, optional
            pdf output destination. The default is None.

        Returns
        -------
        None.

        """

        #plt.style.use('seaborn-talk')
        plt.style.use('seaborn-poster') # Bug with normal x marker !!!

        fig_ts = self.plot_time_spectra()
        fig_es = self.plot_energy_spectra()
        fig_vn = self.vis["North"].plot(self)
        fig_vs = self.vis["South"].plot(self)

        if pdf is not None:
            pdf.savefig(fig_ts)
            pdf.savefig(fig_es)
            pdf.savefig(fig_vn)
            pdf.savefig(fig_vs)

    ###------------------------------------------------------------------------
    def plot_time_spectra(self, n_E_2disp = 6,
                                e_min = 10*u.GeV,
                                e_max = 10*u.TeV,
                                ymin=1.e-20, ax=None):
        """
        Plot a series of lightcurves for some energies.

        Parameters
        ----------
        n_E_2disp : integer, optional
            Number of energies to be sampled. The default is 6.
        e_min : Astropy Qauntity, optional
            Minimum energy. The default is 10*u.GeV.
        e_max : Astropy Qauntity, optional
            Maximal energy. The default is 10*u.TeV.
        ymin : float, optional
            Minimal ordinate value. The default is 1.e-20.
        ax : Matplotlib Axes, optional
            Current axis. The default is None.

        Returns
        -------
        fig : Matplotlib figure
            Current figure.

        """
        if ax is None:
            fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,6))
        else:
            fig = ax.get_figure()

        # Energy sampling
        Emin = max( e_min.value,
                    min(self.Eval).to(e_min.unit).value )
        Emax = min( e_max.to(e_min.unit).value,
                    max(self.Eval).to(e_min.unit).value)
        E = np.logspace(np.log10(Emin), np.log10(Emax), n_E_2disp)*e_min.unit

        ### -----------------------------------
        ### Compute the partial and total fluxes
        ### -----------------------------------

        # Afterglow is always here
        unit     = self.spec_afterglow[0](100*u.GeV).unit
        flx_glow = np.array([f(E).value for f in self.spec_afterglow])
        flx_glow = flx_glow*unit

        if self.prompt:
            # Prompt
            unit = self.spec_prompt[0](100*u.GeV).unit
            flx_prompt = np.array([f(E) for f in self.spec_prompt])
            flx_prompt = flx_prompt*unit

            # Prompt + afterglow
            flx_tot = flx_prompt.value + flx_glow[:(self.id90+1),:].value
            flx_tot = np.concatenate( (flx_tot,flx_glow[self.id90+1:,:].value) )*unit
        else:
            flx_tot = flx_glow

        max_flx_tot = np.max(flx_tot)

        # Attenuated flux
        unit = self.models[0].spectral_model(100*u.GeV).unit
        flx_tot_att = np.array([f.spectral_model(E).value for f in self.models])
        flx_tot_att = flx_tot_att*unit

        ### -----------------------------------
        ### Plot the various fluxes for the E series
        ### -----------------------------------
        with quantity_support():

            for i in range(len(E)):
                color = cm.Dark2(i/len(E))  # cm.cool, cm.rainbow...

                ax.plot(self.tval,flx_tot[:,i],
                        color=color,ls="-", alpha=0.8, lw=1,
                        label="Total")

                if self.prompt:
                    ax.plot(self.tval,flx_glow[:,i],
                            color= color,alpha=0.7,ls=":",
                            label="afterglow")

                    ax.plot(self.tval[:(self.id90+1)],
                            flx_prompt[:,i],
                            color=color,ls="--",alpha=0.5,
                            label="prompt")

                ax.plot(self.tval,flx_tot_att[:,i],
                        color=color,ls="-",marker=".",
                        label=str(round(E[i].value)*E[i].unit))

            if self.prompt:
                ax.axvline(self.t90,color="grey",ls=":",alpha=0.2)
                ax.text(x=self.t90, y= 1.2*ymin,
                         s= "$t_{90}$="+str(round(self.t90.value,1)*self.t90.unit),
                         rotation=0, fontsize=14)
            else:
                ax.axvline(30*u.s,color="black", ls="--", alpha=0.5)
                ax.text(x=35*u.s, y=1.25*ymin, s="30 s",
                        rotation=0, fontsize=14)

                if max(self.tval) >= 1*u.d:
                    ax.axvline(1*u.d,color="black", ls="--", alpha=0.5)
                    ax.text(x=1.05*u.d, y=1.25*ymin, s="1 d",
                            rotation=0, fontsize=14)

                if max(self.tval) >= 2*u.d:
                    ax.axvline(2*u.d,color="black", ls="--", alpha=0.5)
                    ax.text(x=2.05*u.d, y=1.25*ymin, s="2 d",
                            rotation=0, fontsize=14)
        ax.set_xlabel("Time (s)")
        title = "{}: {:>2d} E points".format(self.id,len(self.Eval))
        ax.set_title(title,fontsize=12)
        ax.set_ylim(ymin = ymin, ymax = 2*max_flx_tot)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.grid(which="both",ls="-",color="lightgrey",alpha=0.5)
        single_legend(ax,loc='upper left',fontsize=11,bbox_to_anchor=[1.02,0.5])

        plt.tight_layout()

        return fig
    ###------------------------------------------------------------------------
    def plot_energy_spectra(self, n_t_2disp = 5, ymin = 1e-16,
                            e_min = 1*u.GeV, eindex=2,
                            ax=None):
        """
        Energy spectra for various measurement points in time.

        Note: the plot methods handle the labels, should not be overwritten

        Parameters
        ----------
        n_E_2disp : integer, optional
            Number of curves to show in energy spectra. The default is 5.
        ymin : float, optional
            Minimal flux value in current unit. The default is 1e-16.
        e_min : Quantity Energy, optional
            Minimal energy in the plot. The default is 1*u.GeV.
        eindex : integer, optionnal
            Flux is multiplied by Energy to eindex. The default is 2.
        ax : matplotlib axis, optional
            Figure axis. The default is None.

        Returns
        -------
        fig : Matplotlib figure
            Current figure.

        """

        import matplotlib.pyplot as plt
        from astropy.visualization import quantity_support
        from niceprint import t_str

        if ax is None:
            fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,6))
        else:
            fig = ax.get_figure()

        # Compute number of spectra to be shown
        nspectra = len(self.tval)-1
        if nspectra > n_t_2disp:
            dnt = int(round(nspectra/n_t_2disp))
            tidx = list(range(0,nspectra,dnt))
            # tidx = list(range(0,self.id90))
            # tidx = [12, 13, 14, 15] + tidx # Typical Prompt times
        else:
            dnt=1
            tidx = [1]


        # Plot in color some of the energy spectra and the initial data points
        with quantity_support():
            e_unit    = self.Eval[0].unit
            flux_unit = self.fluxval[0,0].unit
            # In the present version, there is only one prompt spectrum up to t90
            if self.prompt: # if no prompt id90 =-1
                ax.plot(self.E_prompt.to(e_unit),
                        (self.E_prompt**eindex)*self.flux_prompt,
                        alpha=0.5, ls="--", color = "black",
                        label=r"$Prompt \ 0-t_{90}$")

            for i in tidx:
                t    = self.tval[i]
                flux = self.fluxval[i,:] # Afterglow

                # Plot the GRB interpolated spectra
                self.models[i].spectral_model.plot([min(self.Eval),max(self.Eval)],
                                    ax,
                                    energy_unit=e_unit,  energy_power=eindex,
                                    flux_unit = flux_unit,
                                    label="t= {:>s}".format(t_str(t)))

                # Plot the initial data points (should fit)
                c = ax.lines[-1].get_color() # Last color
                ax.plot(self.Eval,
                        self.Eval**eindex*flux,
                        ls='--', lw=1., marker=".", alpha=0.5, color=c)

            # Plot the rest faded out
            for i in range(0,nspectra):
                t = self.tval[i]
                self.models[i].spectral_model.plot([min(self.Eval),max(self.Eval)],ax,
                                    energy_unit = e_unit, energy_power=eindex,
                                    flux_unit = flux_unit,
                                    alpha=0.2, color="grey",lw=1)

            ax.axvline(10*u.GeV, color="grey",alpha=0.2)
            ax.text(x = 10*u.GeV, y=ymin*50, s="10 GeV",
                      rotation=270, fontsize=14,va="bottom")
            title = "{}: z={:3.1f} - {} EBL -{:>2d} Flux points"\
            .format(self.id,self.z, self.eblmodel, len(self.tval))

        if n_t_2disp <= 15:
            ax.legend(fontsize = 12) # Too many t slices
        ax.set_title(title,fontsize = 12)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(ymin=ymin)
        ax.set_xlim(xmin=e_min.to(e_unit).value)

        plt.tight_layout()

        return fig

    ###------------------------------------------------------------------------
    def energy_and_time_2d(self):
        """
        Two-dimensionnal plot of flux versus energy and time.

        Returns
        -------
        None.

        """

        with quantity_support():

            fig, (a1) = plt.subplots(nrows=1,ncols=1,figsize=(12,10))

            h = a1.contourf(self.tval,self.Eval,np.log10(self.fluxval.value.T) )

            a1.set_xscale('log')
            a1.set_yscale('log')
            a1.set_xlabel("Time (s)")
            a1.set_ylabel("E (GeV)")
            cbar = fig.colorbar(h, ax=a1)
            cbar.set_label(str((self.fluxval[0][0]).unit), rotation=270)

    ###------------------------------------------------------------------------
    def energy_over_timeslices(self):
        """
        Plot energy spectra along the time bins.
        The number of Time bin is variable.
        """

        # Compute grid size
        idxmax = len(self.tval)
        ncols=6
        nrows = int(idxmax/ncols)+1
        if idxmax%ncols == 0:
            nrows=nrows-1

        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15,15),
                               sharex=True,
                               sharey=True)

        for icol in range(0,ncols):
            for irow in range(0,nrows):

                # print(irow,icol,icol+ncols*irow)
                idx = icol+ncols*irow
                a = ax[irow][icol]
                if idx < idxmax:
                    a.plot(np.log10(self.Eval.value),
                           np.log10(self.fluxval[idx].value),
                           marker='.',
                           markerfacecolor='r',
                           markeredgecolor='r',
                           label="t= {:06.2f}".format(self.tval[idx]))
                    a.grid("both")
                    a.legend()
                a.set_xlim(-0.5,4.5)
                a.set_ylim(-22,-4)
                # This has been thought ot be necessary - strange it works without
                #if (irow != nrows-1): a.set_xticklabels([])
                # if (icol != 0):       a.set_yticklabels([])

        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none',
                        top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel("$log(E (GeV))$",size=20)
        plt.ylabel("$log(Flux (Gev^{-1} cm^{-2} s^{-1}) )$",size=20)

        fig.suptitle(self.id,size=20,y=0.9)
        plt.subplots_adjust(hspace = .001, wspace=0.001)

    ###------------------------------------------------------------------------
    def plot_altitude_and_flux(self, times,
                          site=None,
                          tshift=0*u.s,
                          Eref = 100*u.GeV,
                          ax=None):
        """
        Plot the altitude and flux of the GRB at the given energy, with
        a ccertain time shift.

        Parameters
        ----------
        times : list of Astropy Time
            list of times, defined elsewhere.
        site : string, optional
            A string defining the site. The default is None.
        tshift : Astropy time, optional
            A time shift for the plot. The default is 0*u.s.
        Eref : Astropy Quantity, optional
            The reference energy. The default is 100*u.GeV.
        ax : Matplotlib Axes, optional
            Current axis. The default is None.

        Returns
        -------
        ax : Matplotlib Axes
            Current axis.

        """

        ax = plt.gca() if ax is None else ax

        ### Altitude sampling for plots - absolute time
        altaz = self.radec.transform_to(AltAz(obstime  = Time(times,format="mjd"),
                                               location = site))

        ### Change time reference if requested
        tref =  self.t_trig - tshift
        times = times - tshift # Change reference

        ### GRB altitude and minimum
        ax.plot(times.datetime,
                altaz.alt.value,
                color="darkblue",alpha=0.5, marker="+",label="Altitude")

        ### FLux points -  Plot typical spectrum from the current source
        axx  = ax.twinx()

        flux = [g.spectral_model(Eref).value \
                for g in self.models]*self.models[0].spectral_model(Eref).unit
        axx.plot((Time(self.t_trig,format="mjd") + self.tval -tshift).datetime,
                  flux,
                  marker=".",
                  ls = "--",
                  lw = 1,
                  color="tab:purple",
                  label = r"$E \sim {}$".format(Eref))
        axx.set_yscale("log")
        axx.legend(loc='center left', bbox_to_anchor=(1.1, 0.9),fontsize=12)


        ### Trigger (dot and vertical line)
        alttrig = self.radec.transform_to(AltAz(obstime  = Time(self.t_trig,format="mjd"),
                                           location = site)).alt.value
        ax.plot(tref.datetime,alttrig,label="Trigger",marker="o",markersize=10,
                color="tab:orange")
        ax.axvline(tref.datetime,ls=":",color="tab:orange")

        ### Trigger + 1 day (dot and vertical line)
        altend = self.radec.transform_to(AltAz(obstime  = Time(self.t_trig,format="mjd") + 1*u.day,
                                          location = site)).alt.value
        ax.plot((tref+1*u.day).datetime,altend,
                label="Trigger + 1 day",marker="o",markersize=10,color="black")
        ax.axvline((tref+1*u.day).datetime,ls=":",color="grey")

        ax.set_ylim(ymin=0,ymax=1.2*max(altaz.alt.value))
        ax.set_ylabel("altitude (°)")

        return ax

    ###------------------------------------------------------------------------
    def time_over_energyband(self):
        """
        Plot time spectra along the energy bins.
        The number of energy bins is fixed.

        To be rearranged.
        """

        # coln = "blue"
        # cols = "red"
        # tb_n1 = self.vis["North"].t_true[0]
        # tb_n2 = self.vis["North"].t_true[1]
        # tb_s1 = self.vis["South"].t_true[0]
        # tb_s2 = self.vis["South"].t_true[1]

        # dt_n = (tb_n2 - tb_n1)/3600
        # dt_s = (tb_s2 - tb_s1)/3600
        # print(tbound_n, tbound_s)

        fig, ax = plt.subplots(nrows=8, ncols=5, figsize=(20,20))
        for irow in range(0,8):
            for icol in range(0,5):

                #print(irow,icol,icol+8*irow)
                idx = icol+5*irow
                a = ax[irow][icol]
                label = "E= {:6.2f}".format(self.Eval[idx])
                #print(idx,self.Eval[idx])
                # Artificially adds 1 s to avoid log(0)
                a.plot(np.log10(self.tval.value+1),np.log10(self.fluxval[:,idx].value),
                       marker='.',markerfacecolor='grey',markeredgecolor='grey',
                       label=label)
                # Display visibility boundaties - Add articifially 1 second  to
                # avoid log(0)
                # if (dt_n):
                #     a.axvline(x=np.log10(tb_n1+1),linestyle =":",color=coln)
                #     a.axvline(x=np.log10(tb_n2+1),linestyle =":",color=coln)
                # if (dt_s):
                #     a.axvline(x=np.log10(tb_s1+1),linestyle =":",color=cols)
                #     a.axvline(x=np.log10(tb_s2+1),linestyle =":",color=cols)

                a.set_xlim(-0.5,5.5)
                a.set_ylim(-22,-4)
                if irow != 7:
                    a.set_xticklabels([])
                if icol != 0:
                    a.set_yticklabels([])
                a.grid("both")
                a.legend()
                #a.set_xscale("log")
                #a.set_yscale("log")

        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none',
                        top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel("$log(t + 1$ $(s))$",size=20)
        plt.ylabel("$log(Flux$ $(Gev^{-1} cm^{-2} s^{-1}) )$",size=20)

        # title = self.id \
        # + " Obs. : North = {:6.2f}h - South = {:6.2f}h".format(dt_n,dt_s)

        # fig.suptitle(title,size=16,y=0.9)
        plt.subplots_adjust(hspace = .001, wspace=0.001)
        plt.show(block=False)

    ###------------------------------------------------------------------------
    def animated_spectra(self, emin = 20 * u.GeV,
                               emax = 10 * u.TeV,
                               savefig = False,
                               outdir='./out/'):

        """
        Create a gif animation of time slices. To be finalised.


        Parameters
        ----------
        emin : Astropy quantity, optional
            Minimal energy. The default is 0.02 * u.TeV.
        emax : Astropy quantity, optional
            Maximal energy. The default is 10 * u.TeV.
        savefig : Boolean, optional
            If True, save the animation to a file. The default is False.
        outdir : string, optional
            Output folder name. The default is './out/'.

        Returns
        -------
        None.

        """

        sys.exit(" animated_spectra to be sorted out - but should work !")

        from matplotlib.animation import FuncAnimation

        def get_model(i):
            return self.models[i].spectral_model

        # Initialise plot
        fig_kw  = dict(num=self.id + '_models')
        fig, ax = plt.subplots(**fig_kw)

        model_init = get_model(0)

        fmin, fmax = model_init(energy=emax), model_init(energy=emin)

        x = np.linspace(emin.value, emax.value, 100) * u.TeV
        y = model_init(energy=x)

        ax.set(xlim=(emin.value,emax.value), ylim=(fmin.value, fmax.value))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy [TeV]')
        ax.set_ylabel('Flux [1 / (cm2 s TeV)]')
        ax.set_ylim([1.e-20, 1.e-6])
        ax.grid(which='both')
        line = ax.plot(x,y, color='k', lw=2)[0]

        ### --------------------------------------------------
        def animate(i):
            model = get_model(i)
            y = model(energy=x)
            line.set_ydata(y)
            ax.set_title(self.id + '; z={:.2f}; dt={:6.2f} s'
                         .format(self.z,self.tval[i].value))
         ### --------------------------------------------------

        anim = FuncAnimation(fig, animate, interval=500,
                             frames=len(self.tval))
        if savefig is True:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            filename = Path(outdir, self.id + '_animate.gif')
            anim.save(str(filename), writer='ffmpeg')
            # anim.save(filename, writer='imagemagick')
        else:
            plt.draw()

###############################################################################
if __name__ == "__main__":

    """
    A standalone function to read a GRB and make various tests

    """

    from configuration import Configuration
    from niceprint import heading

    # This is required to have the EBL models read from gammapy
    os.environ['GAMMAPY_DATA'] = str(Path(Path(__file__).absolute().parent,
                                          "data"))

    ###---------------------------
    ### Input source data - copy of first SoHAPPy steps
    ###---------------------------
    cf = Configuration.command_line()

    cf.ifirst      = 1
    cf.nsrc        = 1 # 250
    cf.prompt_dir  = None
    cf.save_grb    = False
    grblist        = cf.source_ids()

    data_path   = Path(cf.infolder,cf.data_dir) # Input data
    cf.create_output_folder()

    ###---------------------------
    ### visibility parameters
    ###---------------------------
    case = 1

    # Computed from a keyword
    if case==1: cf.visibility = "strictmoonveto"

    # Obtained from a directory as json files
    if case==2: cf.visibility = "visibility/long/vis_301_400_strictmoonveto.json"

    visinfo = cf.decode_visibility_keyword()

    ###---------------------------
    ### Loop over sources
    ###---------------------------
    for i, item in enumerate(grblist):

        fname = cf.prefix+str(item)+cf.suffix

        grb = GammaRayBurst.from_fits(Path(data_path,fname),
                                      prompt  = cf.prompt_dir,
                                      ebl     = cf.ebl_model)

        for loc in ["North","South"]:
            grb.set_visibility(item,loc,info=visinfo)

        heading(grb.id)
        print(grb)

        # Save GRB if requested
        if cf.save_grb : grb.write(Path("./"))

        if grb.vis.keys():
            for loc in ["North","South"]:
                grb.vis[loc].print()

        grb.plot()
        grb.energy_and_time_2d()
        # grb.animated_spectra(savefig=True)
        grb.energy_over_timeslices()
        grb.time_over_energyband()

    # ### -------------------
    # ### Actions for the Data Challenge
    # ### -------------------
    # grb_print = True
    # grb_plot  = True
    # save_grb  = False # (False) GRB saved to disk -> use grb.py main
    # dump2yaml = False # Obsolete
    # dump2fits = False
    # debug     = True # dump2fits


    #     if dump2yaml: grb.model_to_yaml() # Dump GRB model to yaml files

    #     if dump2fits:
    #         filename = grb.models_to_fits(debug=debug) # Dump GRB model to fits files
    #         check_models_from_fits(filename) # Not a method
