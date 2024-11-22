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
import tarfile
import random
import warnings

from pathlib import Path
import numpy as np

import yaml
from yaml.loader import SafeLoader

import matplotlib.pyplot as plt
from matplotlib import cm

import astropy
import astropy.units as u
from astropy.table import Table, QTable
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.visualization import quantity_support

import observatory as obs
from visibility import Visibility

from ebl import EBL_from_file

from niceprint import warning, failure, t_fmt, t_str, heading
from niceplot import single_legend

from gammapy.modeling.models import Models
from gammapy.modeling.models import PointSpatialModel, SkyModel
from gammapy.modeling.models import EBLAbsorptionNormSpectralModel
from gammapy.modeling.models import TemplateSpectralModel

# Transform warnings into errors - useful to find who is guilty !
# warnings.filterwarnings('error')
warnings.filterwarnings('ignore')

__all__ = ["GammaRayBurst"]


# #############################################################################
class GammaRayBurst():
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

    Following a question on the `Slack gammapy` channel on November
    27 :sup:`th` , and the answer by Axel Donath:

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
    ignore = ["id", "filename", "Eval", "tval", "tstart", "tstop", "fluxval",
              "spec_afterglow", "prompt", "id90", "E_prompt", "flux_prompt",
              "spec_prompt", "models", "vis"]
    """ Ignore these members when dumping the class out."""

    e_range = [10*u.GeV, 10*u.TeV]
    """ Default energy range for the energy spectra"""

    t_range = [1*u.s, 5*u.d]
    """ Default time range (after trigger) for the lightcurve"""

    # ##-----------------------------------------------------------------------
    def __init__(self):
        """
        This initializes a default GRB.

        Returns
        -------
        None.

        """

        self.id = 'dummy'
        self.filename = None  # File with the spectral and lightcurve data
        self.z = 0  # Redshift
        self.eblmodel = None

        # GRB properties - Dummy default values
        self.radec = SkyCoord(ra=100*u.deg, dec=-15*u.deg, frame="icrs")
        self.Eiso = 0.*u.erg
        self.Liso = 0.*u.Unit("erg/s")
        self.Epeak = 0.*u.keV
        self.t90 = 0.*u.s
        self.G0H = 0.*u.dimensionless_unscaled
        self.G0W = 0.*u.dimensionless_unscaled
        self.Fpeak = 0.*u.erg
        self.Fpeak_GBM = 0.*u.erg
        self.gamle = 0.*u.dimensionless_unscaled
        self.gamhe = 0.*u.dimensionless_unscaled

        # GRB explosion time and observation
        self.t_trig = Time('2020-01-01T02:00:00', format="isot", scale='utc')

        # ------------
        # ## Afterglow
        # ------------

        # Afterglow Flux table - Spectra at a series of points in time
        self.Eval = [0]*u.GeV  # Energy values
        self.tval = [0]*u.s  # Time values
        self.tstart = self.t_trig  # Start time in seconds
        self.tstop = self.t_trig   # Stop time in seconds

        # Raw flux arrays
        self.fluxval = [0]*u.Unit("1 / (cm2 GeV s)")

        # Afterglow interpolated arrays
        self.spec_afterglow = []

        # -------------------
        # ## Prompt component
        # -------------------
        # Default : no prompt
        self.prompt = False

        # Afterglow Slice id at which the prompt stops
        self.id90 = -1

        # Prompt energy bins
        self.E_prompt = [0]*u.GeV

        # Prompt raw flux value
        self.flux_prompt = 1*u.Unit("1 / (cm2 GeV s)")

        # One Interpolated, non attenuated E-spectrum
        self.spec_prompt = None

        # -------------
        # ## Total flux
        # -------------
        # If the prompt flux exist, this the sum of both flux until the
        # last time bin of the afterglow flux below t90 of the prompt.

        # Gammapy Skymodel models (one per t slice)
        self.models = []

        # -------------
        # ## Visibility
        # -------------
        # Visibility - dictionnary with site entries
        self.vis = dict.fromkeys(["North", "South"])

    # ## ----------------------------------------------------------------------
    @classmethod
    def from_fits(cls,
                  file_name,
                  prompt=None,
                  ebl=None,
                  elimit=None,
                  tlimit=None,
                  dt=0.0,
                  magnify=1,
                  debug=False):
        """
        Read the GRB data from a `fits` file.
        Fluxes are given for a series of (t,E) values
        So far no flux exists beyond the last point (no extrapolation).

        The spectrum is stored as a table, :obj:`TemplateSpectralModel`, that
        takes as an input a series of flux values as a function of the energy.
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
        file_name : Path
            GRB input file name.
        prompt: Boolean
            If True, read information from the associated prompt component.
        ebl : String, optional
            The EBL absorption model considered. If ebl is `built-in`, uses
            the absorbed specrum from the data if available. Can be `None`.
        elimit : Astropy Quantity
            Upper limit to the energies
        tlimit: Astropy Time Quantity
            Upper limit to the time sampling.
        dt: float, optionnal
            A number of Julian days to be added to the present source trigger
            time. The default is 0.0.
        magnify : float, optional
            Flux multiplicative factor to the afterglow model flux for tests.
            Default is 1.
        debug: Boolean, optional
            Debugging flag. The default is False.


        Returns
        -------
        A GammaRayBurst instance.

        """

        # ## -----------------------------------------------------
        # ## Open file, get header, keys, and data - Fill the class members
        # ## -----------------------------------------------------
        # hdul = fits.open(get_filename(filename))
        hdul = fits.open(file_name)
        hdr = hdul[0].header
        keys_0 = list(hdul[0].header.keys())

        cls = GammaRayBurst()  # Default constructor

        cls.filename = file_name
        cls.z = hdr['Z']
        cls.eblmodel = ebl

        # Get identifier by removing all extensions
        cls.id = str(file_name.name).rstrip(''.join(file_name.suffixes))

        cls.radec = SkyCoord(ra=hdr['RA']*u.deg, dec=hdr['DEC']*u.deg,
                             frame="icrs")

        cls.Eiso = hdr['EISO']*u.erg

        if "LISO" in keys_0:
            cls.Liso = hdr["LISO"]*u.Unit("erg/s")  # New large prod. files
        else:
            cls.Liso = 0*u.Unit("erg/s")

        cls.Epeak = hdr['EPEAK']*u.keV
        cls.t90 = hdr['Duration']*u.s
        cls.G0H = hdr['G0H']
        if "G0W" in keys_0:  # Not in SHORTFITS
            cls.G0W = hdr['G0W']

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

        # ##--------------------------
        # ## GRB trigger time
        # ##--------------------------
        # If the input file does not contain a trigger date, then set a date
        # by default. It is very likely that this source belongs to a set
        # for which a list of explosion (trigger) times should be generated.

        if "GRBJD" in keys_0:  # Not in SHORTFITS
            cls.t_trig = Time(hdr['GRBJD']*u.day + dt*u.day,
                              format="jd", scale="utc")
        elif "GRBTIME" in keys_0:  # Add start date
            cls.t_trig = Time(hdr['GRBTIME'] + dt*u.day,
                              format="jd", scale="utc")
        else:
            warning(f"{__name__:}.py: Trigger time absent from file,"
                    f" using random value")
            cls.t_trig += random.uniform(0, 365)*u.day

        # ##--------------------------
        # ## Time intervals - so far common to afterglow and prompt
        # ##--------------------------
        tval = QTable.read(hdul["TIMES (AFTERGLOW)"])

        # Temporary : in short GRB fits file, unit is omitted
        # The flux value is given on an interval ("col0", "col1°°.
        # The default is to consider that the flux is valid at the end of the
        # interval.

        if isinstance(tval[0][0], astropy.units.quantity.Quantity):
            cls.tval = np.array(tval[tval.colnames[0]].value)*tval[0][0].unit
        else:  # In SHORT GRB, the time bin is given, [t1, t2].
            cls.tval = np.array(tval["col1"])*u.s

        # Select point up to tmax. Force last point to be tmax.
        if tlimit is not None and tlimit <= cls.tval[-1]:
            ninit = len(cls.tval)
            cls.tval = cls.tval[(tlimit - cls.tval) > 0]
            cls.tval[-1] = tlimit.to(u.s)
            warning(f" Data up to {cls.tval[-1]:5.3f} "
                    f"restricted to {tlimit:5.3f}")
            print(f" GRB {cls.id:}: slices {ninit:} -> {len(cls.tval):}")

        cls.tstart = cls.t_trig
        cls.tstop = cls.t_trig + cls.tval[-1]
        cls.to_min = cls.tval[0]
        cls.to_max = cls.tval[-1]

        # ##--------------------------
        # ## Afterglow Energies - Limited to elimit if defined
        # ##--------------------------
        tab_key = "Energies (afterglow)"
        col_key = Table.read(hdul[tab_key]).colnames[0]
        cls.Eval = Table.read(hdul[tab_key])[col_key].quantity
        cls.Eval = np.array(cls.Eval)*cls.Eval[0].unit

        if elimit is not None and elimit <= cls.Eval[-1]:
            warning(f" Data up to {cls.Eval[-1]:5.3f} "
                    f"restricted to {elimit:5.3f}")
            cls.Eval = cls.Eval[cls.Eval <= elimit]

        # ##--------------------------
        # ## Afterglow flux
        # ##--------------------------

        # Get flux, possibly already absorbed
        if ebl == "in-file":  # Read default absorbed model
            flux = QTable.read(hdul["EBL-ABS. SPECTRA (AFTERGLOW)"])
        else:
            flux = QTable.read(hdul["SPECTRA (AFTERGLOW)"])

        flux_unit = u.Unit(flux.meta["UNITS"])
        if str(flux_unit).find("ph") > -1:
            flux_unit = flux_unit/u.Unit("ph")  # Removes ph

        # Store the flux. Note the transposition
        itmax = len(cls.tval)-1
        jEmax = len(cls.Eval) - 1
        cls.fluxval = np.zeros((itmax+1, jEmax+1))*flux_unit

        for i in range(0, itmax+1):
            for j in range(0, jEmax+1):
                cls.fluxval[i][j] = magnify * flux[j][i] * flux_unit  # transp!

        # Build time series of interpolated spectra - limited to dtmax
        for i in range(len(cls.tval)):
            glow = TemplateSpectralModel(energy=cls.Eval.astype(float),
                                         values=cls.fluxval[i],
                                         interp_kwargs={"values_scale": "log"})
            cls.spec_afterglow.append(glow)

        # ##--------------------------
        # ## Prompt - a unique energy spectrum
        # ##--------------------------

        # Get the prompt if potentially visible and if requested
        cls.prompt = False  # No prompt component was found
        if prompt is not None:
            if Path(prompt).is_dir():
                # Deduce prompt folder from GRB path name
                cls.spec_prompt = cls.get_prompt(folder=prompt, debug=debug)
                if cls.spec_prompt is not None:
                    cls.prompt = True
            else:
                sys.exit(f" {prompt:} is not a valid data folder :")

        # ##--------------------------
        # ## Total attenuated spectra
        # ##--------------------------
        for i, _ in enumerate(cls.tval):
            spec_tot = cls.spec_afterglow[i]
            if cls.prompt:
                if i <= cls.id90:
                    spec_tot += cls.spec_prompt[i]
            m = SkyModel(spectral_model=cls.EBLabsorbed(spec_tot),
                         spatial_model=PointSpatialModel(lon_0=cls.radec.ra,
                                                         lat_0=cls.radec.dec,
                                                         frame="icrs"),
                         name="model_" + str(i))
            cls.models.append(m)

        # ## Close fits file
        hdul.close()

        return cls

    # ##-----------------------------------------------------------------------
    @classmethod
    def historical_from_yaml(cls,
                             filename,
                             ebl=None,
                             magnify=1,
                             dtmax=0.5*u.h,
                             nebin=25):
        """
        Read the characteristics of a parameterised GRB from a `yaml` file.
        It is assumed that the energy and time decays are not correlated and
        follow two independent power laws.
        The function associates spectra to a list of time intervals
        in order to comply with the more general case in this class.
        Note that the spectra are considered to be from the afterglow only.

        Parameters
        ----------
        item : string
            Source identifier related to the input yaml file
        ebl : String, optional
            The EBL absorption model considered. If ebl is `built-in`, uses
            the absorbed specrum from the data if available. Can be `None`.
        magnify : float, optional
            Flux multiplicative factor to the afterglow model flux for tests.
        elimit : Astropy Quantity
            Upper limit to the energies
        tlimit: Astropy Time Quantity
            Upper limit to the time sampling.
        dtmax : Astropy Quantity, optional
            Time bin maximal length. The default is 0.5*u.h.
        debug: Boolean, optional
            Debugging flag. The default is False.

        Returns
        -------
        A `GammaRayBurst` instance.
        """

        cls = GammaRayBurst()  # This calls the constructor

        cls.filename = filename

        with open(cls.filename, encoding="utf-8") as f:
            data = yaml.load(f, Loader=SafeLoader)

        cls.z = data["z"]
        cls.eblmodel = ebl
        cls.id = data["name"]

        cls.radec = SkyCoord(data["ra"], data["dec"], frame='icrs')
        cls.Eiso = u.Quantity(data["Eiso"])

        cls.Epeak = u.Quantity(data['Epeak'])
        cls.t90 = u.Quantity(data['t90'])
        cls.G0H = data['G0H']
        cls.G0W = data['G0W']
        cls.Fpeak = u.Quantity(data['Fluxpeak'])
        cls.gamle = data['gamma_le']
        cls.gamhe = data['gamma_he']

        # ##--------------------------
        # ## GRB trigger time
        # ##--------------------------
        cls.t_trig = Time(data["t_trig"], format="datetime", scale="utc")

        # ##--------------------------
        # ## Time intervals for the simulation
        # ##--------------------------
        ntbin = data["ntbin"]

        # Implement a time binning in log space, never exceeding a certain
        # duration
        if ntbin != 1:

            # First generate in log
            tsamp = np.logspace(np.log10(cls.t_range[0].to(u.s).value),
                                np.log10(cls.t_range[1].to(u.s).value),
                                ntbin)

            # Finf when differences exceed the maximal allowed value
            idmax = np.where(np.diff(tsamp) > dtmax.to(u.s).value)[0][0]

            # Generate flat from that value
            tflat = np.arange(tsamp[idmax],
                              cls.t_range[1].to(u.s).value,
                              dtmax.to(u.s).value)

            # Concatenate the two arrays
            cls.tval = np.concatenate((tsamp[:idmax], tflat))*u.s

        else:  # A single time window
            cls.tval = np.array([cls.t_range[0].to(u.s).value,
                                 cls.t_range[1].to(u.s).value])*u.s

        cls.tstart = cls.t_trig + cls.tval[0]
        cls.tstop = cls.t_trig + cls.tval[-1]

        # ##--------------------------
        # ## Afterglow Energies - Default Energu ranges
        # ##--------------------------
        eo_min = 10*u.GeV
        eo_max = (10*u.TeV).to(eo_min.unit)

        cls.Eval = np.logspace(np.log10(eo_min.value),
                               np.log10(eo_max.value), nebin)*eo_min.unit

        # ##--------------------------
        # ## Afterglow flux
        # ##--------------------------

        flux_unit = u.Unit("1/(cm2 GeV s)")
        cls.fluxval = np.zeros((len(cls.tval), len(cls.Eval)))*flux_unit

        for i, t in enumerate(cls.tval):
            for j, E in enumerate(cls.Eval):
                dnde = (u.Quantity(data["K"])*(E/data["E0"])**-data["gamma"]
                        * (t/data["t0"])**-data["beta"])
                dnde = dnde.to(1/(u.cm**2)/u.s/eo_min.unit)
                # print(f"t={cls.tval[i]:8.1f}    E={cls.Eval[j]:8.1f} "
                #       f"   ->  dN/dE = {dnde:8.2e}")
                cls.fluxval[i][j] = magnify * dnde.to(flux_unit)

        # ##--------------------------
        # ## Prompt - No prompt foreseen in this case
        # ##--------------------------
        cls.prompt = False

        # ##--------------------------
        # ## Total attenuated spectra
        # ##--------------------------
        # See note in from_fits for explanatins on some parameters
        for i, t in enumerate(cls.tval):

            tab = TemplateSpectralModel(energy=cls.Eval.astype(float),
                                        values=cls.fluxval[i],
                                        interp_kwargs={"values_scale": "log"})

            cls.spec_afterglow.append(tab)  # Needed for display

            m = SkyModel(spectral_model=cls.EBLabsorbed(tab),
                         spatial_model=PointSpatialModel(lon_0=cls.radec.ra,
                                                         lat_0=cls.radec.dec,
                                                         frame="icrs"),
                         name="model_"+str(i))

            cls.models.append(m)

        return cls

    # ##-----------------------------------------------------------------------
    @classmethod
    def time_resolved_prompt(cls, filename, glowname=None,
                             ebl=None,
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
        cls. z = z
        cls.magnifiy = magnify
        cls.eblmodel = ebl

        # Copy missing data from the corresponding afterglow
        if glowname is not None:
            glow = GammaRayBurst.from_fits(glowname)
            # print(glow)
            cls.z = glow.z
            cls.radec = glow.radec
            cls.Eiso = glow.Eiso
            cls.Liso = glow.Liso
            cls.Epeak = glow.Epeak
            cls.Fpeak = glow.Fpeak
            cls.t90 = glow.t90
            cls.G0H = glow.G0H
            cls.G0W = glow.G0W
            cls.gamle = glow.gamle
            cls.gamhe = glow.gamhe
            cls.Fpeak_GBM = glow.Fpeak_GBM
            cls.t_trig = glow.t_trig

        # Reads prompt data
        hdul = fits.open(filename)
        cls.id = Path(filename.name).stem

        cls.Eval = Table.read(hdul, hdu=1)["energy"].quantity*u.TeV

        cls.tval = Table.read(hdul, hdu=2)["time"].quantity*u.s
        cls.tstart = cls.t_trig
        cls.tstop = cls.t_trig + cls.tval[-1]

        flux = Table.read(hdul, hdu=3)
        flux_unit = u.Unit("1/(cm2 TeV s)")
        icol_t = len(flux.colnames)  # column number - time
        jrow_E = len(flux[flux.colnames[0]])  # row number

        # Note the transposition from flux to fluxval
        cls.fluxval = np.zeros((icol_t, jrow_E))*flux_unit
        for i in range(0, icol_t):
            for j in range(0, jrow_E):
                f = flux[j][i]
                if f > 1:
                    # print(i,j,f)
                    f = 0  # Correct a bug in event #172 - to be removed
                cls.fluxval[i][j] = magnify*f*flux_unit  # transp!

        for i in range(len(cls.tval)):
            # Note that TableModel makes an interpolation
            prpt = TemplateSpectralModel(energy=cls.Eval.astype(float),
                                         values=cls.fluxval[i],
                                         interp_kwargs={"values_scale": "log"})
            # Use the afterglow place holder to store the prompt
            cls.spec_afterglow.append(prpt)

            m = SkyModel(spectral_model=cls.EBLabsorbed(prpt, ebl),
                         spatial_model=PointSpatialModel(lon_0=cls.radec.ra,
                                                         lat_0=cls.radec.dec,
                                                         frame="icrs"),
                         name="model_"+str(i))
            cls.models.append(m)

        # Not adding the averaged on time prompt
        cls.prompt = False

        hdul.close()

        return cls

    # ##-----------------------------------------------------------------------
    def EBLabsorbed(self, tab):
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

        attflux = tab  # Initialise to the unabsorbed flux
        if self.eblmodel is None or self.eblmodel == 'in-file':
            return attflux

        if self.eblmodel != "gilmore":
            eblabs = \
             EBLAbsorptionNormSpectralModel.read_builtin(self.eblmodel,
                                                         redshift=self.z)
            attflux = tab*eblabs
        else:
            eblabs = EBL_from_file("data/ebl/others/ebl_gilmore12-10GeV.dat")
            # Change flux for absorption, recreate model
            # I di not fine a clever way
            attenuated = []
            for E in self.Eval:
                attenuated.append(tab(E).value*eblabs(E, self.z))
            attenuated *= tab(np.mean(self.Eval)).unit
            attflux = TemplateSpectralModel(energy=self.Eval.astype(float),
                                            values=attenuated,
                                            interp_kwargs={"values_scale":
                                                           "log"})

        return attflux

    # ##-----------------------------------------------------------------------
    def set_visibility(self,
                       item, loc, observatory="CTAO",
                       tmin=None, tmax=None,
                       info=None, n_night=None, n_skip=None,
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

        Parameters
        ----------
        item : TYPE
            DESCRIPTION.
        loc : TYPE
            DESCRIPTION.
        observatory : TYPE, optional
            DESCRIPTION. The default is "CTAO".
        tmin : TYPE, optional
            DESCRIPTION. The default is None.
        tmax : TYPE, optional
            DESCRIPTION. The default is None.
        info : TYPE, optional
            DESCRIPTION. The default is None.
        n_night : TYPE, optional
            DESCRIPTION. The default is None.
        n_skip : TYPE, optional
            DESCRIPTION. The default is None.
        status : TYPE, optional
            DESCRIPTION. The default is "".
        dbg : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        # Define observation window
        tmin = self.tstart if tmin is None\
            else self.tstart + tmin
        tmax = self.tstop if tmax is None\
            else self.tstart + tmax

        # ## Update the default - At this stage the visibility is maximal
        self.vis[loc] = Visibility(pos=self.radec,
                                   site=obs.xyz[observatory][loc],
                                   tmin=tmin, tmax=tmax,
                                   name=self.id+"_"+loc,
                                   status=status)

        if info == "permanent":
            self.vis[loc].status = info  # Keep track of that special case
            return self

        if isinstance(info, dict):  # A dictionnary not a keyword

            if "altmin" in info.keys():  # Compute from dictionnary

                # Supersede max. number of nights and nights to be skipped
                if n_night is not None:
                    info["depth"] = n_night
                if n_skip is not None:
                    info["skip"] = n_skip

                # Compute from dictionnary elements
                self.vis[loc] = self.vis[loc].compute(param=info, debug=dbg)

            else:  # A dictionnary of all visibilities (from a .json file)
                self.vis[loc] = Visibility.from_dict(info[str(item) + "_"+loc])

        # Not a dictionnary
        elif info == "built-in":
            # sys.exit("{}.py : built-in vis. reading to be reimplemented"
            #          .format(__name__))
            # infolder is not passed
            self.vis[loc] = Visibility.from_fits(self, loc=loc)

        else:  # Special or from disk

            if info == "forced":  # Infinitite nights
                self.vis[loc] = self.vis[loc].force_night()

            else:  # Binary - Obsolete - Archived bin files might be not valid
                with open(Path(info, self.id+"_"+loc+"_vis.bin"), "rb") as f:
                    self.vis[loc] = pickle.load(f)

# ##-----------------------------------------------------------------------
    def altaz(self, loc="", dt=0*u.s):
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
        if not isinstance(dt, astropy.units.quantity.Quantity):
            sys.exit("Time has to be a quantity (weird behaviour otherwise")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            altaz = self.radec.transform_to(AltAz(obstime=dt + self.t_trig,
                                                  location=self.vis[loc].site))
        return altaz

    # ##-----------------------------------------------------------------------
    def get_prompt(self, folder=None, strict=False, debug=False):
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

        # Find back the source id number from the id
        # (assumes only the id has digits)
        gid = ""
        for s in self.id:
            if s.isdigit():
                gid += s

        if strict:
            if gid in [6, 30, 191]:
                failure(f" GRB prompt {gid:} simulated with a Wind profile")
                failure(" It cannot be associated to the current afterglow !")
                return -1, None

        # ##----------------------------
        # ## Compute weihts to "Mask" the time interval beyond t90
        # ##----------------------------

        # Find the closest time bin at t90 in the Afterglow
        # id90 is the index of the time following the t90 value
        # t90 is therefore between the index id90-1 and id90
        self.id90 = np.where(self.tval >= self.t90)[0][0]
        if self.id90 == 0:
            warning(f"{self.id:} : t90 = {self.t90:} "
                    F"is before the first afterglow time bin")
            return None

        # The given flux is per unit of time.
        # In order to take into account that the last time bin is not
        # completely covered, it is corrected by the fraction in time in this
        # bin.
        dtprompt = self.t90 - self.tval[self.id90 - 1]  # Prompt is active
        dtslice = self.tval[self.id90] - self.tval[self.id90 - 1]  # Total
        fraction = dtprompt/dtslice

        # All slices have weight 1 before the slice containing t90, the
        # fraction above for the slice containing t90, and zero beyond
        weight = np.concatenate((np.ones(self.id90),
                                 [fraction],
                                 np.zeros(len(self.tval) - self.id90-1)))

        # ##----------------------------
        # ## Read prompt data
        # ##----------------------------
        # Open file - read the lines
        filename = Path(folder, gid+"-spec.dat")
        if not filename.exists():
            failure("{filename:} does no exist")
            return None

        with open(filename, "r", encoding="utf-8") as file:
            lines = file.readlines()

        # get Gamma and z
        data = lines[0].split()
        gamma_prompt = float(data[2])
        redshift = float(data[5])

        # Consistency check - Zeljka data are rounded at 2 digits
        if np.abs(self.z - redshift) > 0.01 \
           or np.abs(self.G0H - gamma_prompt) > 0.01:
            failure(f" {self.id:10s}: "
                    f"Afterglow / prompt : z= {self.z:4.2f} / {redshift:4.2f}"
                    f"  G0H/gamma= {self.G0H:5.2f} / {gamma_prompt:5.2f}"
                    f" (G0W= {self.G0W:5.2f})")

        # Get units - omit "ph" in flux unit
        data = lines[1].split()
        unit_e = data[0][1:-1]
        unit_flx = ''.join([x + " " for x in data[2:]])[1:-2]

        # Get flux points - assign units to data, use afterglow units
        data = lines
        Eval = []
        fluxval = []
        for line in lines[4:]:
            data = line.split()
            Eval.append(float(data[0]))
            fluxval.append(float(data[1]))

        self.E_prompt = Eval*u.Unit(unit_e)
        self.flux_prompt = fluxval*u.Unit(unit_flx)

        if unit_e != self.Eval.unit:
            if debug:
                warning(f"Converting energy units from {unit_e:}"
                        f" to {self.Eval[0].unit:}")
            self.E_prompt = self.E_prompt.to(self.Eval[0].unit)

        if unit_flx != self.fluxval[0].unit:
            if debug:
                warning(f"Converting flux units from {unit_flx:} "
                        f"to {self.fluxval[0][0].unit:}")
            self.flux_prompt = self.flux_prompt.to(self.fluxval[0][0].unit)

        # ##----------------------------
        # ## Create a list of weighted models for non-zero weights
        # ##----------------------------
        models = []
        for i in range(self.id90+1):
            flux_w = self.flux_prompt*weight[i]
            models.append(TemplateSpectralModel(energy=self.E_prompt,
                                                values=flux_w,
                                                interp_kwargs={"values_scale":
                                                               "log"}))

        if debug:
            print("Prompt associated to ", self.id)
            print(f"Gamma = {gamma_prompt:}  redshift = {redshift:}")
            print(f"  {unit_e:8} : {unit_flx:10}")
            for E, flux in zip(self.E_prompt, self.flux_prompt):
                print(f"  {E.value:8.2f} : {flux.value:10.2e} ")
            islice = 0
            print(f" {'Time':>8} - {'Weight':>5} ", end="")
            for t, w in zip(self.tval, weight):
                print(f" {t:8.2f} - {w:5.2f} ", end="")
                print("*") if islice == self.id90 else print("")
                islice += 1

            print("Prompt is between: ", self.tval[self.id90-1],
                  self.tval[self.id90])

        return models

    # ##------------------------------------------------------------------------
    # ## INPUT/OUPUT
    # ##------------------------------------------------------------------------
    def __str__(self):
        """
        Printout the GRB properties (Not the visibilities).
        """

        txt = ""
        if self.filename is not None:
            txt += f'  Read from      : {self.filename}\n'
        txt += f'  RA, DEC        : {self.radec.ra.value:>8.2f} '\
               f'{self.radec.ra.value:>8.2f}\n'
        txt += f'  Redshift       : {self.z:>8.2f}\n'
        txt += f'  EBL model      : "{self.eblmodel:}"\n'
        txt += f'  Eiso           : {self.Eiso:6.2e}\n'
        txt += f'  Epeak          : {self.Epeak:>8.2f}\n'
        txt += f'  t90            : {self.t90:>8.2f}\n'
        txt += f'  G0H / G0W      : {self.G0H:>8.2f} / {self.G0W:>8.2f}\n'
        txt += f'  Flux peak      : {self.Fpeak:>8.2f}\n'
        txt += f'  Flux peak (GBM): {self.Fpeak_GBM:>8.2f}\n'
        txt += f'  gamma LE / HE  : {self.gamle:>8.2f} / {self.gamhe:>8.2f}\n'

        txt += f'  t_trig         : {Time(self.t_trig, format="iso"):}\n'
        txt += f'  Obs. window    : {t_fmt(min(self.tval)):>8.2f} '\
               f'{t_fmt(max(self.tval)):>8.2f}\n'
        txt += f'  Energy window  : {min(self.Eval):>8.2f} '\
               f'{ max(self.Eval):>8.2f}\n'
        txt += f'  Eff. duration  : '\
               f'{t_fmt(self.tval[-1]-self.tval[0]):>8.2f}\n'
        txt += f'  Bins : E, t    : {len(self.Eval):4d} '\
               f'{len(self.tval):4d} \n'
        if self.prompt:
            txt += ' Prompt component  : \n'
            txt += '  Up to slice    : {self.id90:3d}\n'
            txt += "  Bins : E       : {len(self.E_prompt):3d}\n"
        else:
            txt += '  + Prompt component not considered\n'
        if self.vis is None:
            txt += '  + Visibility not available\n'

        return txt

    # ##-----------------------------------------------------------------------
    def write_to_bin(self, folder, debug=True):
        """
        Write current instance to a binary file for further
        reuse.

        Parameters
        ----------
        folder : String
            Folder name.
        debug : Boolean, optional
            If True, talks a bit. The default is True.

        Returns
        -------
        None.

        """

        filename = Path(folder, self.id + ".bin")
        with open(filename, "wb") as outfile:
            pickle.dump(self, outfile)

        if debug:
            print(" GRB saved to : {filename:}")

    # ##-----------------------------------------------------------------------
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
        folder = Path(output, self.filename + "/")
        if not folder.is_dir():
            os.makedirs(folder)

        print(" Now dumping ", self.filename, "to ", folder)

        # Dump measurement times
        name = Path(folder, self.filename + "_times.txt")
        with open(name, "w", encoding="utf-8") as out:
            for i, t in enumerate(self.tval):
                print(f"{i:3d} {t:>8.2f}", file=out)

        # Dump spectra and add to tar file
        with tarfile.open(Path(folder,
                               self.filename+".tar.gz"), "w:gz") as tar:

            for i, spec in enumerate(self.models):
                specname = f"spec_{str(i+1).zfill(3):3s}.yaml"
                specfilename = Path(folder, specname)

                # Dump data in yaml format
                with open(specfilename, "w", encoding="utf-8") as out:
                    print(Models([spec]).to_yaml(), file=out)

                # Put spectrum in tar file
                tar.add(specfilename, arcname=specname)

                # Delete file
                Path.unlink(specfilename)

    # ##-----------------------------------------------------------------------
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
        yaml_list = [Models([spec]).to_yaml() for spec in self.models]
        maxlength = max([len(m) for m in yaml_list])
        if debug:
            print("   Max. yaml model length is ", maxlength)

        # Fits filename
        output_name = Path(folder, "DC_" + self.id + ".fits")
        print(" Now dumping ", grb.name, "to ", output_name,
              " as a FITS table")

        # Create an astropy table with 3 columns for spectra (same # of rows)
        table = Table()
        table["time_start"] = np.append([0e0],
                                        grb.tval[:-1].value)*grb.tval.unit
        table["time_stop"] = grb.tval
        table["mod"] = maxlength*" "

        for i, mod in enumerate(yaml_list):
            table["mod"][i] = mod

        # Add this stage, the data can be dumped into a fits file with
        # table.write(output_name, overwrite=True)
        # However, in order to modify the header we add a step
        hdu_spectra = fits.BinTableHDU(data=table, name="spectra")

        # Add explosion time in the header
        header = hdu_spectra.header
        header["trigger"] = self.t_trig.isot

        # Add visibility
        table = Table()
        if len(self.vis["North"].t_true[0]) != 0:
            table["tstart"] = [t[0].isot for t in self.vis["North"].t_true]
            table["tstop"] = [t[1].isot for t in self.vis["North"].t_true]
        else:
            table["tstart"] = [-1]
            table["tstop"] = [-1]

        hdu_visN = fits.BinTableHDU(data=table, name="vis_N")

        table = Table()
        if len(self.vis["South"].t_true[0]) != 0:
            table["tstart"] = [t[0].isot for t in self.vis["South"].t_true]
            table["tstop"] = [t[1].isot for t in self.vis["South"].t_true]
        else:
            table["tstart"] = [-1]
            table["tstop"] = [-1]
        hdu_visS = fits.BinTableHDU(data=table, name="vis_S")

        # Create the hdulist, write to output
        hdul = fits.HDUList()
        hdul.append(hdu_spectra)
        hdul.append(hdu_visN)
        hdul.append(hdu_visS)

        # Write to file
        hdul.writeto(output_name, overwrite=True)
        hdul.close()

        return output_name

    # ##-----------------------------------------------------------------------
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
        print("Trigger/explosion time :", hdu.header["TRIGGER"])

        # Get spectral data
        data = hdu.data
        print(data.columns)
        # data["time_start"]
        # data["time_stop"]
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

    # ##-----------------------------------------------------------------------
    # ## PLOTS
    # ##-----------------------------------------------------------------------
    def plot(self, e_unit="GeV", f_unit="1/ (GeV cm2 s)", pdf=None):
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

        fig_ts = self.plot_time_spectra(e_unit=e_unit,
                                        f_unit=f_unit).get_figure()
        fig_es = self.plot_energy_spectra(e_unit=e_unit,
                                          f_unit=f_unit).get_figure()

        # Visibility might not have been computed yet
        fig_vn = self.vis["North"].plot(self) \
            if self.vis["North"] is not None else None
        fig_vs = self.vis["South"].plot(self) \
            if self.vis["South"] is not None else None

        if pdf is not None:
            pdf.savefig(fig_ts)
            pdf.savefig(fig_es)
            if fig_vn is not None:
                pdf.savefig(fig_vn)
            if fig_vs is not None:
                pdf.savefig(fig_vs)

    # ##-----------------------------------------------------------------------
    def plot_time_spectra(self,
                          n_E_2disp=6,
                          t_min=0*u.s,
                          t_max=1*u.yr,
                          e_min=10*u.GeV,
                          e_max=10*u.TeV,
                          e_index=2,
                          e_unit="GeV",
                          f_unit="1/ (GeV cm2 s)",
                          f_min=1.e-20*u.Unit("1/ (GeV cm2 s)"),
                          ax=None):
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
        e_unit : string, optional
            Astropy unit for energies in the plot. The default is GeV.
        f_unit: string, optional
            Astropy units for flux. The default is "1/ (GeV cm2 s)".
        f_min: Astropy Quantity, optional
            Minimal flux to be displayed. The default is
            1.e-20*u.Unit("1/ (GeV cm2 s)")
        ax : Matplotlib Axes, optional
            Current axis. The default is None.

        Returns
        -------
        ax : matplotlib axes
            Current axis.

        """

        if ax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
        else:
            fig = ax.get_figure()

        # Time limits
        t_min = max(t_min, min(self.tval))
        t_max = min(t_max, max(self.tval))

        # Compute the mask corresponding to this window
        t_unit = self.tval[0].unit
        time_mask = ((self.tval.value >= t_min.to(t_unit).value)
                     & (self.tval.value <= t_max.to(t_unit).value))
        times = self.tval

        # Energy sampling - convert units and restricty to available data
        e_min = e_min.to(e_unit)
        e_max = e_max.to(e_unit)
        f_min = f_min.to(f_unit)

        e_min = max(e_min.value, min(self.Eval).to(e_unit).value)
        e_max = min(e_max.value, max(self.Eval).to(e_unit).value)
        E = np.logspace(np.log10(e_min), np.log10(e_max),
                        n_E_2disp)*u.Unit(e_unit)

        # ## -----------------------------------
        # ## Compute the partial and total fluxes
        # ## -----------------------------------

        # Afterglow is always here
        # unit     = self.spec_afterglow[0](100*u.GeV).unit
        flx_glow = np.array([f(E).to(f_unit).value
                             for f in self.spec_afterglow])
        flx_glow = flx_glow*u.Unit(f_unit)

        if self.prompt:
            # Prompt
            unit = self.spec_prompt[0](100*u.GeV).unit
            flx_prompt = np.array([f(E).to(f_unit) for f in self.spec_prompt])
            flx_prompt = flx_prompt*u.Unit(f_unit)

            # Prompt + afterglow
            glow_before = flx_glow[:(self.id90+1), :].value
            flx_tot = flx_prompt.value + glow_before
            glow_after = flx_glow[self.id90+1:, :].value
            flx_tot = np.concatenate((flx_tot, glow_after)) * unit.to(f_unit)
        else:
            flx_tot = flx_glow

        max_flx_tot = np.max(flx_tot[time_mask])  # In the time window

        # Attenuated flux
        unit = self.models[0].spectral_model(100*u.GeV).unit
        flx_tot_att = np.array([f.spectral_model(E).value
                                for f in self.models])
        flx_tot_att = flx_tot_att*unit.to(f_unit)

        # ## -----------------------------------
        # ## Plot the various fluxes for the E series
        # ## -----------------------------------

        with quantity_support():

            for i, _ in enumerate(E):
                color = cm.Dark2(i/len(E))  # cm.cool, cm.rainbow...

                ax.plot(times, flx_tot[:, i],
                        color=color, ls="-", alpha=0.8, lw=1, label="Total")

                if self.prompt:
                    ax.plot(times, flx_glow[:, i],
                            color=color, alpha=0.7, ls=":", label="afterglow")

                    ax.plot(times[:(self.id90+1)],
                            flx_prompt[:, i],
                            color=color, ls="--", alpha=0.5, label="prompt")

                ax.plot(times, flx_tot_att[:, i],
                        color=color, ls="-", marker=".",
                        label=str(round(E[i].value)*E[i].unit))

            if self.prompt:
                ax.axvline(self.t90, color="grey", ls=":", alpha=0.2)
                ax.text(x=self.t90, y=1.2*f_min,
                        s="$t_{90}$="+str(round(self.t90.value,
                                                1)*self.t90.unit),
                        rotation=0, fontsize=14)
            else:
                if t_min <= 30*u.s <= t_max:
                    ax.axvline(30*u.s, color="black", ls="--", alpha=0.5)
                    ax.text(x=35*u.s, y=1.25*f_min, s="30 s",
                            rotation=0, fontsize=14)

                if t_min <= 1*u.d <= t_max:
                    ax.axvline(1*u.d, color="black", ls="--", alpha=0.5)
                    ax.text(x=1.05*u.d, y=1.25*f_min, s="1 d",
                            rotation=0, fontsize=14)

                if t_min <= 2*u.d <= t_max:
                    ax.axvline(2*u.d, color="black", ls="--", alpha=0.5)
                    ax.text(x=2.05*u.d, y=1.25*f_min, s="2 d",
                            rotation=0, fontsize=14)

            ax.set_ylim(ymin=f_min, ymax=2*max_flx_tot)
            ax.set_xlim(xmin=0.95*t_min, xmax=1.05*t_max)

        ax.set_xlabel("Time (s)")
        title = f"{self.id:}: {len(self.Eval):>2d} E points"
        ax.set_title(title, fontsize=12)
        ax.set_yscale("log")
        ax.set_xscale("log")

        ax.grid(which="both", ls="-", color="lightgrey", alpha=0.5)

        single_legend(fig, loc='upper left', fontsize=11,
                      bbox_to_anchor=[1.02, 0.5])

        plt.tight_layout()

        return ax

    # ##------------------------------------------------------------------------
    def plot_energy_spectra(self, n_t_2disp=5,
                            e_min=1*u.GeV,
                            e_index=2,
                            e_unit="GeV",
                            f_unit="1/ (GeV cm2 s)",
                            f_min=1.e-20*u.Unit("1/ (GeV cm2 s)"),
                            ax=None):
        """
        Energy spectra for various measurement points in time.

        Note: the plot methods handle the labels, should not be overwritten

        Parameters
        ----------
        n_E_2disp : integer, optional
            Number of curves to show in energy spectra. The default is 5.
        e_min : Quantity Energy, optional
            Minimal energy in the plot. The default is 1*u.GeV.
        e_index : integer, optionnal
            Flux is multiplied by Energy to e_index. The default is 2.
        e_unit : string, optional
            Astropy unit for energies in the plot. The default is GeV.
        f_unit: string, optional
            Astropy units for flux. The default is "1/ (GeV cm2 s)".
        f_min: Astropy Quantity, optional
            Minimal flux to be displayed. The default is
            1.e-20*u.Unit("1/ (GeV cm2 s)")
        ax : matplotlib axis, optional
            Figure axis. The default is None.

        Returns
        -------
        ax : matplotlib axes
            Current axis.

        """

        if ax is None:
            _, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

        # Compute number of spectra to be shown
        nspectra = len(self.tval) - 1
        if nspectra > n_t_2disp:
            dnt = int(round(nspectra/n_t_2disp))
            tidx = list(range(0, nspectra, dnt))
            # tidx = list(range(0,self.id90))
            # tidx = [12, 13, 14, 15] + tidx # Typical Prompt times
        else:
            dnt = 1
            tidx = [1]

        # Convert energies and flux
        e_min = e_min.to(e_unit)
        f_min = f_min.to(f_unit)

        # Plot in color some of the energy spectra and the initial data points
        with quantity_support():

            # In the present version, there is only one prompt spectrum
            # up to t90
            if self.prompt:  # if no prompt id90 =-1
                ax.plot(self.E_prompt.to(e_unit),
                        (self.E_prompt**e_index)*self.flux_prompt,
                        alpha=0.5, ls="--", color="black",
                        label=r"$Prompt \ 0-t_{90}$")

            for i in tidx:
                t = self.tval[i]
                flux = self.fluxval[i, :]  # Afterglow

                # Plot the GRB interpolated spectra
                self.models[i].spectral_model\
                    .plot([min(self.Eval), max(self.Eval)],
                          ax,
                          energy_power=e_index,
                          # flux_unit = f_unit,
                          label=f"t= {t_str(t):>s}")

                # Plot the initial data points (should fit)
                c = ax.lines[-1].get_color()  # Last color
                ax.plot(self.Eval,
                        self.Eval**e_index*flux,
                        ls='--', lw=1., marker=".", alpha=0.5, color=c)

            # Plot the rest faded out
            for i in range(0, nspectra):
                t = self.tval[i]
                self.models[i].spectral_model\
                    .plot([min(self.Eval), max(self.Eval)],
                          ax,
                          energy_power=e_index,
                          # flux_unit = f_unit,
                          alpha=0.2, color="grey", lw=1)

            ax.axvline(10*u.GeV, color="grey", alpha=0.2)
        #       ax.text(x = 10*u.GeV, y=ymin*50, s="10 GeV",
        #                   rotation=270, fontsize=14,va="bottom")
            title = f"{self.id}: z={self.z:3.1f} - "\
                    f"{self.eblmodel} EBL -{len(self.tval):>2d} Flux points"

        if n_t_2disp <= 15:
            ax.legend(fontsize=12)  # Too many t slices
        ax.set_title(title, fontsize=12)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(ymin=(e_min**e_index*f_min).value)
        ax.set_xlim(xmin=e_min.to(e_unit).value)

        plt.tight_layout()

        return ax

    # ##------------------------------------------------------------------------
    def energy_and_time_2d(self):
        """
        Two-dimensionnal plot of flux versus energy and time.

        Returns
        -------
        None.

        """

        with quantity_support():

            fig, (a1) = plt.subplots(nrows=1, ncols=1, figsize=(12, 10))

            h = a1.contourf(self.tval, self.Eval,
                            np.log10(self.fluxval.value.T))

            a1.set_xscale('log')
            a1.set_yscale('log')
            a1.set_xlabel("Time (s)")
            a1.set_ylabel("E (GeV)")
            cbar = fig.colorbar(h, ax=a1)
            cbar.set_label(str((self.fluxval[0][0]).unit), rotation=270)

    # ##-----------------------------------------------------------------------
    def energy_over_timeslices(self, ncols=6):
        """
        Plot energy spectra along the time bins.
        The number of Time bin is variable.
        """

        # Compute grid size
        idxmax = len(self.tval)
        nrows = int(idxmax/ncols)+1
        if idxmax % ncols == 0:
            nrows = nrows-1

        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 15),
                               sharex=True,
                               sharey=True)

        for icol in range(0, ncols):
            for irow in range(0, nrows):

                # print(irow,icol,icol+ncols*irow)
                idx = icol+ncols*irow
                a = ax[irow][icol]
                if idx < idxmax:
                    a.plot(np.log10(self.Eval.value),
                           np.log10(self.fluxval[idx].value),
                           marker='.',
                           markerfacecolor='r',
                           markeredgecolor='r',
                           label=f"t= {self.tval[idx]:6.2f}")
                    a.grid("both")
                    a.legend()
                a.set_xlim(-0.5, 4.5)
                a.set_ylim(-22, -4)
                # This has been thought necessary - strange it works without
                # if (irow != nrows-1): a.set_xticklabels([])
                # if (icol != 0):       a.set_yticklabels([])

        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none',
                        top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel("$log(E (GeV))$", size=20)
        plt.ylabel("$log(Flux (Gev^{-1} cm^{-2} s^{-1}) )$", size=20)

        fig.suptitle(self.id, size=20, y=0.9)
        plt.subplots_adjust(hspace=.001, wspace=0.001)

    # ##-----------------------------------------------------------------------
    def plot_altitude(self, times, site=None, ax=None, tshift=0*u.s):
        """
        Plot the altitude of the GRB with a certain time shift.

        Parameters
        ----------
        times : list of Astropy Time
            list of times, defined elsewhere.
        site : string, optional
            A string defining the site. The default is None.
        tshift : Astropy time, optional
            A time shift for the plot. The default is 0*u.s.
        ax : Matplotlib Axes, optional
            Current axis. The default is None.

        Returns
        -------
        ax : Matplotlib Axes
            Current axis.
        """

        ax = plt.gca() if ax is None else ax

        # ## Altitude sampling for plots - absolute time
        altaz = self.radec.transform_to(AltAz(obstime=Time(times,
                                                           format="mjd"),
                                              location=site))

        # ## Change time reference if requested
        tref = self.t_trig - tshift
        times = times - tshift  # Change reference

        # ## GRB altitude and minimum
        ax.plot(times.datetime,
                altaz.alt.value,
                color="darkblue", alpha=0.5,
                marker="", label="Altitude")

        # ## Trigger (dot and vertical line)
        alttrig = self.radec.transform_to(AltAz(obstime=Time(self.t_trig,
                                                             format="mjd"),
                                                location=site)).alt.value
        ax.plot(tref.datetime, alttrig, label="Trigger",
                marker="o", markersize=10,
                color="tab:orange")
        ax.axvline(tref.datetime, ls=":", color="tab:orange")

        # ## Trigger + 1 day (dot and vertical line)
        trigtime = Time(self.t_trig, format="mjd") + 1*u.day
        altend = self.radec.transform_to(AltAz(obstime=trigtime,
                                               location=site)).alt.value
        ax.plot((tref + 1*u.day).datetime, altend,
                label="Trigger + 1 day",
                marker="o", markersize=10, color="black")
        ax.axvline((tref+1*u.day).datetime, ls=":", color="grey")

        ax.set_ylim(ymin=0, ymax=1.2*max(altaz.alt.value))
        ax.set_ylabel("Altitude (°)")

        return ax

    # ##-----------------------------------------------------------------------
    def plot_flux(self, ax=None, Eref=100*u.GeV, tshift=0*u.s):
        """
        Plot the flux of the GRB at the given energy, with
        a certain time shift.

        Parameters
        ----------
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

        flux = [g.spectral_model(Eref).value
                for g in self.models]*self.models[0].spectral_model(Eref).unit
        ax.plot((Time(self.t_trig, format="mjd")
                 + self.tval - tshift).datetime,
                flux,
                marker=".", ls="--", lw=1, color="tab:purple",
                label=rf"$E \sim {Eref}$")
        ax.set_yscale("log")

        ax.legend()

        return ax

    # ##-----------------------------------------------------------------------
    def time_over_energyband(self, ncols=5):
        """
        Plot time spectra along the energy bins.
        The number of energy bins is fixed.

        Parameters
        ----------
        ncols : integer, optional
            Number of columns in the plot. The default is 5.

        Returns
        -------
        None.

        """

        ncols = 5
        nrows = int(len(self.Eval)/ncols)
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20))
        for irow in range(0, nrows):
            for icol in range(0, ncols):

                # print(irow,icol,icol+8*irow)
                idx = icol+5*irow
                a = ax[irow][icol]
                label = f"E= {self.Eval[idx]:6.2f}"

                # Artificially adds 1 s to avoid log(0)
                a.plot(np.log10(self.tval.value),
                       np.log10(self.fluxval[:, idx].value),
                       marker='.',
                       markerfacecolor='grey', markeredgecolor='grey',
                       label=label)

                a.set_xlim(-0.5, 5.5)
                a.set_ylim(-22, -4)
                if irow != nrows-1:
                    a.set_xticklabels([])
                if icol != 0:
                    a.set_yticklabels([])
                a.grid("both")
                a.legend()

        fig.add_subplot(111, frameon=False)

        # hide tick and tick label of the big axis
        plt.tick_params(labelcolor='none',
                        top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel("$log(t $ $(s))$", size=20)
        plt.ylabel("$log(Flux$ $(Gev^{-1} cm^{-2} s^{-1}) )$", size=20)

        title = self.id
        # + " Obs. : North = {:6.2f}h - South = {:6.2f}h".format(dt_n,dt_s)

        fig.suptitle(title, size=16, y=0.9)
        plt.subplots_adjust(hspace=.001, wspace=0.001)
        plt.show(block=False)

    # ##------------------------------------------------------------------------
    def animated_spectra(self,
                         emin=20 * u.GeV,
                         emax=10 * u.TeV,
                         savefig=False,
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
        fig_kw = {"num": self.id + '_models'}
        fig, ax = plt.subplots(**fig_kw)

        model_init = get_model(0)

        fmin, fmax = model_init(energy=emax), model_init(energy=emin)

        x = np.linspace(emin.value, emax.value, 100) * u.TeV
        y = model_init(energy=x)

        ax.set(xlim=(emin.value, emax.value), ylim=(fmin.value, fmax.value))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy [TeV]')
        ax.set_ylabel('Flux [1 / (cm2 s TeV)]')
        ax.set_ylim([1.e-20, 1.e-6])
        ax.grid(which='both')
        line = ax.plot(x, y, color='k', lw=2)[0]

        # ## --------------------------------------------------
        def animate(i):
            model = get_model(i)
            y = model(energy=x)
            line.set_ydata(y)
            ax.set_title(f'{self.id:} ; '
                         f'z={self.z:.2f}; dt={self.tval[i].value:6.2f} s')
        # ## --------------------------------------------------

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

    # A standalone function to read a GRB and make various tests

    import seaborn as sns
    from configuration import Configuration

    # Bigger texts and labels
    sns.set_context("paper")  # poster, talk, notebook, paper

    os.environ["HAPPY_IN"] = r"D:\\CTAO\SoHAPPy\input"
    os.environ["HAPPY_OUT"] = r"D:\\CTAO\SoHAPPy\output"

    # This is required to have the EBL models read from gammapy
    os.environ['GAMMAPY_DATA'] = str(Path(Path(__file__).absolute().parent,
                                          "data"))

    # ##---------------------------
    # ## Input source data - copy of first SoHAPPy steps
    # ##---------------------------
    sys.argv = ["", "-c", "data/config_ref_omega.yaml"]
    cf = Configuration.command_line()

    cf.ifirst = [343]
    cf.nsrc = 1
    cf.prompt_dir = None
    cf.save_grb = False
    cf.emax = u.Quantity("10 TeV")
    cf.tmax = u.Quantity("10 h")
    grblist = cf.source_ids(os.environ["HAPPY_IN"])

    # ##---------------------------
    # ## visibility parameters
    # ##---------------------------
    case = 2

    # Computed from a keyword
    if case == 1:
        cf.visibility = "strictmoonveto"

    # Obtained from a directory as json files
    if case == 2:
        cf.visibility = "visibility/long/vis_301_400_strictmoonveto.json"

    visinfo = cf.decode_visibility_keyword()

    # ##---------------------------
    # ## Loop over sources
    # ##---------------------------
    for item, fname in enumerate(grblist):

        grb = GammaRayBurst.from_fits(Path(fname),
                                      prompt=cf.prompt_dir,
                                      ebl=cf.ebl_model,
                                      elimit=cf.elimit,
                                      tlimit=cf.tlimit)

        for loc in ["North", "South"]:
            grb.set_visibility(item, loc, info=visinfo)

        heading(grb.id)
        print(grb)

        # Save GRB if requested
        if cf.save_grb:
            grb.write_to_bin(Path("./"))

        if grb.vis.keys():
            for loc in ["North", "South"]:
                grb.vis[loc].print()

        grb.plot()
        grb.energy_and_time_2d()
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
    # Dump GRB model to fits files
    #         filename = grb.models_to_fits(debug=debug)
    #         check_models_from_fits(filename) # Not a method
