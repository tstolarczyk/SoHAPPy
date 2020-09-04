# -*- coding: utf-8 -*-
"""
This module contains the classes and tools to handle a GRB object, i.e.
a collection of lightcurve bins for each of which is given an energy spectrum.

Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import yaml
from yaml import CLoader as Loader
#from glob import glob

import sys
import warnings

from pathlib import Path

import numpy as np

import astropy
import astropy.units         as u
from   astropy.table         import Table
from   astropy.io            import fits
from   astropy.time          import Time
from   astropy.coordinates   import SkyCoord
from   astropy.coordinates   import EarthLocation
from   astropy.coordinates   import AltAz
#from   astropy.coordinates   import Distance

import gammapy
# Suppress warnings and in particular the deprecation warning
warnings.filterwarnings("ignore")

if (gammapy.__version__ == "0.12"):
    from gammapy.spectrum.models import Absorption, AbsorbedSpectralModel
    from gammapy.spectrum.models import TableModel
    from gammapy.image.models    import SkyPointSource

if (gammapy.__version__ == "0.16"):
    from gammapy.modeling.models import Absorption, AbsorbedSpectralModel
    from gammapy.modeling.models import TemplateSpectralModel
    from gammapy.modeling.models import PointSpatialModel

# Prevent downloading Earth position correction from outside
# Results in a warning : https://docs.astropy.org/en/stable/utils/iers.html
from astropy.utils import iers
iers.conf.auto_download = False

# iers.conf.auto_max_age = None # to be checked
warnings.filterwarnings("error") # Help to identify the warning locations

# The following information can be obtained from EarthLocation.of_site
# but this requires to have an internet connection
# Get all available site names :
# sitelist = EarthLocation.get_site_names()

# Get coordinattes of Paranal and LaPalma
# Note : with the Gammapy-0.12 environment,this requires to
#        downgrade numpy to version 1.16.2 (instead on >1.18)
#        See discussion here :
#        https://github.com/astropy/astropy/issues/9306

# xyz_north = EarthLocation.of_site('Roque de los Muchachos')
# xyz_south = EarthLocation.of_site('Paranal Observatory')

xyz_north = EarthLocation.from_geocentric( 5327448.9957829,
                                          -1718665.73869569,
                                           3051566.90295403,
                                           unit="m")
xyz_south = EarthLocation.from_geocentric( 1946404.34103884,
                                          -5467644.29079852,
                                          -2642728.20144425,
                                           unit="m")

# Values take from Maria Grazia Bernardini - See Readme, July 28, 2020
# xyz_north = EarthLocation.from_geocentric( 5327285.09211954,
#                                           -1718777.11250295,
#                                           3051786.7327476,
#                                           unit="m")
# xyz_south = EarthLocation.from_geocentric(1946635.7979987,
#                                           -5467633.94561753,
#                                           -2642498.5212285,
#                                           unit="m")

# xyz_north = EarthLocation.from_geodetic('342.1184', '28.7606', 2326. * u.meter)
# xyz_south = EarthLocation.from_geodetic('289.5972',  '-24.6253', 2635. * u.meter)

__all__ = ["GammaRayBurst","GRBDummy"]

###############################################################################
class GRBDummy():
    """
    A list of obaservation times and an integrated flux at each point
    Useful for slice merging tests or prototype for e.g. Fermi extrapolated
    spectra.

    """

    def __init__(self,tmin = 6*u.h, tmax = 54*u.h, nbin = 7, beta= 1):
        self.tmeas = tmeas = np.linspace(tmin,tmax,7)
        self.flux  = u.Unit("1/(GeV cm2 s)")*(tmeas/(1*u.h))**(-beta)

    def __str__(self):
        with np.printoptions(precision=3, suppress=False):
            print("Points measured at  : ",self.tmeas)
            print("Flux                : ",self.flux)
        return ""

    def plot(self,ax):
        ax.plot(self.tmeas,self.flux,marker="*")

###############################################################################
class GammaRayBurst(object):
    """
    Class to store GRB properties

    The GRB information is read from files.
    A GRB is composed of several parameters among which, a name, a redshift,
    a position in the sky, a visibility window, measurement points at which
    an energy spectrum is given.
    Each energy spectrum is given as a series of flux measurements.
    The GRB properties like z, Eiso are generated in the population synthesis
    (G. Ghirlanda) and they connect the afterglow and prompt emission together
    (Other parameters like the magnetic field are independent in the two
    models).

    """

    ###########################################################################
    def __init__(self, ebl_model=None,z=0):
        """
        This initializes a default GRB.
        The default visibility has been put by hand for backward compatibility
        with the old file format and definitions where it was assumed that the
        GRB was seen without any delay. The dummy GRB is visible in North and
        South (inspired form GRB 289 in Lara's 1000 GRB files)

        Parameters
        ----------
        ebl_model : string, optional
            A normalized name of an EBL abosrtpion model (e.gh. domibnguez")
            from the list available in gammapy. The default is None.
        z : Dimensionnles quantity, optional
            DESCRIPTION. The default is 0*u.dimensionless_unscaled.

        Returns
        -------
        None.

        """

        self.site_keys = ["North","South"] # Put it somewhere else !
        self.pos_site  = {"North":xyz_north, "South":xyz_south}
        self.name      = 'dummy'
        self.reduction_factor = 1

        # GRB properties - theory
        if (ebl_model !=None):
            self.eblabs = Absorption.read_builtin(ebl_model) # Initialise absorption
        else:
            self.eblabs = None
            print(" EBL absorption not defined")
        self.z        = z
        self.radec    = SkyCoord(ra=100*u.degree, dec=-15*u.degree,
                                 frame='icrs')
        self.Eiso     = 0.*u.erg
        self.Epeak    = 0.*u.keV
        self.t90      = 0.*u.s
        self.G0H      = 0.*u.dimensionless_unscaled
        self.G0W      = 0.*u.dimensionless_unscaled
        self.FluxPeak = 0.*u.erg
        self.gamma_le = 0.*u.dimensionless_unscaled
        self.gamma_he = 0.*u.dimensionless_unscaled

        # GRB alert received
        self.altmin   = 0*u.deg
        self.t_trig   = Time('2000-01-01 02:00:00', scale='utc')
        self.site     = {"North": 'Roque de los Muchachos',
                         "South": 'Paranal'}

        # Visible any moment of the year from the site
        self.vis         = {"North":True,"South":True}
        # Visible any moment within the 24h after the trigger from the site
        self.vis_tonight = {"North":True,"South":True}
        # Visible at the moment of the GRB trigger
        self.vis_prompt  = {"North":True,"South":True}

        # The start and stop dates are searched during 24h after the trigger
        # and correspond to the first visibility interval.
        # Does not report a second visibility interval during 24h,
        # that should be possible only if the target is promptly visible
        # (about 15% of the targets)
        # By default seen in the North.
        # Default times are defined for the first GRB sample with no
        # visibikity given. The North visibility is chosen to have the
        # trigger therein. The South visibility starts after the trigger.
        # This allows having the two conditions explored.

        # Visible above min. alt.
        self.t_true = { "North":[
                                 [Time('2000-01-01 00:00:00', scale='utc'),
                                  Time('2000-01-01 08:00:00', scale='utc')]
                                ],
                        "South":[
                                 [Time('2000-01-01 03:00:00', scale='utc'),
                                  Time('2000-01-01 11:00:00', scale='utc')]
                                ]}

        # Astronomical twilight
        self.t_twilight = {"North":[
                                     [Time('2000-01-01 00:00:00', scale='utc'),
                                      Time('2000-01-01 8:00:00', scale='utc')]
                                   ],
                           "South":[
                                     [Time('2000-01-01 03:00:00', scale='utc'),
                                      Time('2000-01-01 11:00:00', scale='utc')]
                                   ]}

        # Rise and set of the target
        self.t_event    = {"North":[
                                    [Time('2000-01-01 00:00:00', scale='utc'),
                                     Time('2000-01-01 08:00:00', scale='utc')]
                                    ],
                           "South":[
                                    [Time('2000-01-01 03:00:00', scale='utc'),
                                     Time('2000-01-01 11:00:00', scale='utc')]
                                    ]}

        # Flux table - Flux at a series of points
        self.Eval           = [0]*u.GeV
        self.tval           = [0]*u.s
        self.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")

        # GRB spectral and spatial model
        self.spectra = [] # Gammapy models (one per t slice)
        self.pointing = None

        if (gammapy.__version__=="0.12"):
            self.spatial  =  SkyPointSource(lon_0=0*u.deg,lat_0=0*u.deg)
        if (gammapy.__version__=="0.16"):
            self.spatial = PointSpatialModel(lon_0=0*u.deg,lat_0=0*u.deg)

        # Stacked obs
        # useless self.stack_obs = None

        return

    ###########################################################################
    @classmethod
    def read(cls, filename, ebl= None):

        """
        Fluxes are given for a series of (t,E) values
        So far no flux is given beyond the last point.

        The spectrum is stored as a table, TableModel, that takes as an input
        a series of flux values as a function of the energy.
        The model will return values interpolated in log-space with
        math::`scipy.interpolate.interp1d`, returning zero for energies
        outside of the limits of the provided energy array.

        Parameters
        ----------
        cls : GammaRayBurst class
            GammaRayBurst class instance
        filename : STRING
            GRB input file name
        ebl : STRING, optional
            The EBL absorption model considered. There is no default set.

        Returns
        -------
        A GammaRayBurst instance.

        """

        grb = cls(ebl_model=ebl) # This calls the constructor

        hdul = fits.open(filename)

        grb.name = Path(filename.name).stem

        ### GRB properties ---
        hdr = hdul[0].header
        grb.radec    = SkyCoord(ra=hdr['RA']*u.degree,
                                 dec=hdr['DEC']*u.degree,
                                 frame='icrs')
        grb.z        = hdr['Z']
        grb.Eiso     = hdr['EISO']*u.erg
        grb.Epeak    = hdr['EPEAK']*u.keV
        grb.t90      = hdr['Duration']*u.s
        grb.G0H      = hdr['G0H']
        grb.G0W      = hdr['G0W']
        grb.Fluxpeak = hdr['EISO']*u.erg
        grb.gamma_le = hdr['LOWSP']
        grb.gamma_he = hdr['HIGHSP']

        ### Visibilities ---
        grbvis = Table.read(hdul,hdu=1)

        # Visibility has been computed with this minimum altitude
        grb.altmin      = grbvis.meta["MIN_ALT"]*u.deg # Minimum altitude

        # GRB alert received
        grb.t_trig      = Time(hdr['GRBJD']*u.day,format="jd",scale="utc")
        grb.vis         = {"North":hdr['V_N_ANYT'],"South":hdr['V_S_ANYT']}
        grb.vis_tonight = {"North":hdr['V_N_TNG'], "South":hdr['V_S_TNG']}
        grb.vis_prompt  = {"North":hdr['V_N_PR'],  "South":hdr['V_S_PR']}

        # Decode default time windows
        #-----------
        def f(t):
            t = Time(t.data,format="jd",scale="utc")
            return {"North": [ [ t[0:2][0], t[0:2][1] ] ],
                    "South": [ [ t[2:4][0], t[2:4][1] ] ]}
        #-----------
        grb.t_true     = f(grbvis["True"])
        grb.t_twilight = f(grbvis["Twilight"])
        grb.t_event    = f(grbvis["Event"])

        ### Energies, times and fluxes ---
        grb.Eval     = Table.read(hdul,hdu=2)["Energies (afterglow)"].quantity
        grb.tval     = Table.read(hdul,hdu=3)["Times (afterglow)"].quantity
        flux         = Table.read(hdul,hdu=4)
        #flux_unit    = u.Unit(flux.meta["UNITS"])/u.Unit("ph") # Removes ph
        # flux_unit    = u.Unit(flux.meta["UNITS"]) # Removes ph
        flux_unit = u.Unit("1/(cm2 GeV s)")

        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Note the transposition from flux to fluxval
        grb.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                grb.fluxval[i][j] = flux[j][i]*flux_unit # transp!

        for i,t in enumerate(grb.tval):
            # Note that TableModel makes an interpolation
            if (gammapy.__version__ == "0.12"):
                tab = TableModel(energy   = grb.Eval,
                             values       = grb.fluxval[i],
                             norm         = 1.,
                             values_scale = 'lin')
            if (gammapy.__version__ == "0.16"):
                tab = TemplateSpectralModel(energy       = grb.Eval,
                             values       = grb.fluxval[i],
                             norm         = 1.,
                             interp_kwargs={"values_scale": "linear"})

            if (grb.eblabs != None):
                absmodel = AbsorbedSpectralModel(spectral_model = tab,
                                                 absorption    = grb.eblabs,
                                                 parameter     = grb.z)
                grb.spectra.append(absmodel)


        return grb

    ###########################################################################
    @classmethod
    def read_prompt(cls, filename, glow= None, ebl= None,
                    z=0*u.dimensionless_unscaled):
        """
        Read prompt data froma file and associate it to the afterglow
        information if requested (or keep the default from the constructor
        otherwise)

        Parameters
        ----------
        cls : GammaRayBurst
            Class intance
        filename : String
            DESCRIPTION.
        glow : boolean, optional
            If defined, the Prompt characteritics are obatined from the
            afterglow with same id. The default is None.
        ebl : string, optional
            The name of an EBL absoprtion model in the Gammapy data.
            The default is None.
        z : float, optional
            GRB redshift. The default is 0*u.dimensionless_unscaled.

        Returns
        -------
        grb : GammaRayBurst
            A GammaRayBurst instance

        """

        if (glow != None):
            grb = glow # Copy afterglow parameters
            #grb.z = 0
            # Flux table - Flux at a series of points
            grb.Eval           = [0]*u.GeV
            grb.tval           = [0]*u.s
            grb.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")
            grb.spectra = [] # Gammapy models (one per t slice)
        else:
            grb = cls(ebl_model = ebl,z=z) # This calls the constructor

        hdul = fits.open(filename)

        grb.name = Path(filename.name).stem

        # Energies, times and fluxes
        grb.Eval     = Table.read(hdul,hdu=1)["energy"].quantity*u.TeV
        grb.tval     = Table.read(hdul,hdu=2)["time"].quantity*u.s
        flux         = Table.read(hdul,hdu=3)
        flux_unit = u.Unit("1/(cm2 TeV s)")

        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Note the transposition from flux to fluxval
        grb.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                f = flux[j][i]
                if (f>1):
                    # print(i,j,f)
                    f =0 # Correct a bug in event #172 - to be removed
                grb.fluxval[i][j] = f*flux_unit # transp!

        for i,t in enumerate(grb.tval):
            # Note that TableModel makes an interpolation
            if (gammapy.__version__ == "0.12"):
                tab = TableModel(energy   = grb.Eval,
                             values       = grb.fluxval[i],
                             norm         = 1.,
                             values_scale = 'lin')
            if (gammapy.__version__ == "0.16"):
                tab = TemplateSpectralModel(energy       = grb.Eval,
                             values       = grb.fluxval[i],
                             norm         = 1.,
                             interp_kwargs={"values_scale": "linear"})

            if (grb.eblabs != None):
                absmodel = AbsorbedSpectralModel(spectral_model = tab,
                                                      absorption    = grb.eblabs,
                                                      parameter     = grb.z)
                grb.spectra.append(absmodel)

        return grb

    ###########################################################################
    @classmethod
    def read_old(cls, path, ebl = None,reduction_factor=1,
             pointing = None):

        """
        Reads one of the 10 old GRB test files.
        In this temporary files :
        - The GRB position is not given
        - The GRB occurence time is not given

        The data consist in time intervals for which an energy specturm is
        given. The time intervals are set by hand because of a file format
        that do not permit to retrieve them.

        In order to be consistent with the new file format the tval values
        are recomputed from the time interval boundaries (start dates).
        A flux is added at the last time point which has the value of the
        flux in the last interval. This is then coherent with the new file
        format and in particular.


        Parameters
        ----------
        cls : GammaRayBurst class
            A class instance.
        path : STRING
            The folder where to find the GRB data file.
        ebl : STRING, optional
            The EBL absorption model. The default is 'dominguez'.

        Returns
        -------
        A GammaRayBurst instance.

        """

        ###------------------------------
        def get_grb_properties(filename):
            """Get GRB properties"""
            with open(filename,'r') as stream:
                data = yaml.load(stream,Loader=Loader)
            delta_t = data['time intervals'].split()
            intervals = []
            for idx in range(0,24,2):
                # in seconds
                intervals.append([float(delta_t[idx]), float(delta_t[idx + 1])] * u.s)
            data['time intervals'] = intervals

            return data
        ###------------------------------

        grb = cls(ebl_model = ebl) # This calls the constructor

        data = grb.get_grb_properties(path.joinpath('grb_properties.txt'))

        grb.name          = data['name']
        grb.z             = data['redshift']

        # HACK, waiting for Lara
        #time_interval     = data['time intervals']
        dt = [[30., 50.], [50., 80.], [80., 125.], [125., 200.],
              [200., 315.], [315., 500.], [500., 800.], [800., 1250.],
              [1250., 2000.], [2000., 3150.], [3150., 5000.],
              [5000., 10000.]]

        # Reads spectral shape from each time interval
        grb.spectra = []
        fluxlist = []
        tvalues = []

        for t in dt:

            filename = '{}_{:.0f}-{:.0f}.txt'.format(data["name"],t[0],t[1])
            table = Table.read(path.joinpath(filename), format='ascii')

            energy = table['col1'] * u.TeV
            flux = (table['col2'] / reduction_factor) * u.Unit('1 / (cm2 s TeV)')
            fluxlist.append(flux) # Each flux contains len(energy) elements
            tvalues.append(t[0])

        # Duplicate the n-1 information to the last point
        fluxlist.append(fluxlist[-1]) # Add the last flux to the last time
        tvalues.append(dt[-1][1])

        # Fill canonical members
        grb.tval = tvalues * u.s
        grb.Eval     = energy # yes the last one - but they are all athe same
        grb.fluxval  = np.asarray(fluxlist) \
                       *fluxlist[0][0].unit # asarray remove units !

        for i,t in enumerate(grb.tval):
            if (gammapy.__version__ == "0.12"):
                tab = TableModel(energy   = grb.Eval,
                              values      = grb.fluxval[i],
                              norm         = 1.,
                              values_scale = 'lin')
            if (gammapy.__version__ == "0.16"):
                tab = TemplateSpectralModel(energy = grb.Eval,
                              values               = grb.fluxval[i],
                              norm                 = 1.,
                              interp_kwargs={"values_scale": "linear"})

            absmodel = AbsorbedSpectralModel(spectral_model = tab,
                                      absorption     = grb.eblabs,
                                      parameter      = grb.z)
            grb.spectra.append(absmodel)

        grb.pointing = pointing

        return grb

    ###------------------------------------------------------------------------
    def altaz(self,loc="",dt=0*u.s):
        """
        Get altitude azimuth for the GRB at a given site at GRB time t (s)

        Parameters
        ----------
        location : string
            Either North or South
        time : Quantity (time)
            The time elapsed since the trigger time

        Returns
        -------
        altaz : astropy.coordinates.AltAz
            The GRB poistion in the sky at the given time and site

        """
        if (type(dt) is not astropy.units.quantity.Quantity):
            sys.exit("Time has to be a quantity (weird behaviour otherwise")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            altaz = self.radec.transform_to(AltAz(obstime  = dt + self.t_trig,
                                                  location = self.pos_site[loc]))

        return altaz

    ###------------------------------------------------------------------------
    def update_visibility(self,vis,debug=False):

        self.altmin               = vis.altmin
        self.vis[vis.loc]         = vis.vis
        self.vis_prompt[vis.loc]  = vis.prompt
        self.vis_tonight[vis.loc] = vis.tonight

        self.t_event[vis.loc]    = vis.t_above
        self.t_twilight[vis.loc] = vis.t_night

        if (vis.tonight):
            self.t_true[vis.loc]     = vis.t_vis
        else:
            self.t_true[vis.loc]     = [-9]

        if (debug):
            print(" Visibility updated for altmin = ",self.altmin)

        return
    ###------------------------------------------------------------------------
    def __str__(self):
        """ Printout the GRB properties """

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ttrig = self.t_trig.datetime

        # After glow duration
        txt = '\n'
        txt += "==================================================================\n"
        txt += "===                             {:10s}                     ===\n"\
        .format(self.name)
        txt += "==================================================================\n"
        txt += '  RA, DEC            : {:6.2f} {:6.2f}\n'\
        .format(self.radec.ra, self.radec.dec)
        txt += '  Redshift           : {:6.2f}\n'.format(self.z)
        #ètxt += '  Distance           : {:6.2f}\n'.format(Distance(z=self.z,unit=u.Glyr))
        txt += '  Eiso               : {:6.2e}\n'.format(self.Eiso)
        txt += '  Epeak              : {:6.2f}\n'.format(self.Epeak)
        txt += '  t90                : {:6.2f}\n'.format(self.t90)
        txt += '  G0H / G0W          : {:6.2f} / {:6.2f}\n' \
        .format(self.G0H,self.G0W)
        txt += '  Flux peak          : {:6.2f}\n'.format(self.FluxPeak)
        txt += '  gamma LE / HE      : {:6.2f} / {:6.2f}\n' \
        .format(self.gamma_le,self.gamma_he)

        txt += '  t_trig             : {}\n'.format(ttrig)
        txt += '  Duration           : {:6.2f} {:10.2f}\n' \
            .format( (self.tval[-1]-self.tval[0]).to(u.d),
                     (self.tval[-1]-self.tval[0]).to(u.s))
        # txt += '  Bins : E, t, dt    : {:4d} {:4d} {:4d} \n' \
        # .format(len(self.Eval), len(self.tval), len(self.time_interval))
        txt += '  Bins : E, t   : {:4d} {:4d} \n' \
        .format(len(self.Eval), len(self.tval))

        # for t in self.tval: print(t)

        return txt

   ###------------------------------------------------------------------------
    def print_visibility(self,loc=None,log=None, alt=None):
        """
        Print out the GRB night, above-the-horizon periods, and default
        visibility window.

        Parameters
        ----------
        loc : String, optional
            Either "North" or "South". The default is None.
        log : TextIO, optional
            Log file. The default is None.
        alt : Visibility instance
            Get the values from the Visibility object (used for comparisons)

        Returns
        -------
        None.

        """

        if (loc == None):
            self.show_visibility(loc="North",log=log)
            self.show_visibility(loc="South",log=log)
        else:

            if (alt==None):
                label       = "Orig. vis"
                vis         = self.vis[loc]
                vis_tonight = self.vis_tonight[loc]
                vis_prompt  = self.vis_prompt[loc]
                altmin      = self.altmin
                t_event     = self.t_event[loc]
                t_twilight  = self.t_twilight[loc]
                t_true      = self.t_true[loc]
            else :
                label       = "Alt.  vis"
                vis         = alt.vis
                vis_tonight = alt.tonight
                vis_prompt  = alt.prompt
                altmin      = alt.altmin
                t_event     = alt.t_above
                t_twilight  = alt.t_night
                t_true      = alt.t_vis

            log.prt("================= {}   {:10s} {:6s}  ================"
              .format(label,self.name,loc))
            log.prt(' Visible    : {} - tonight, prompt : {}, {} > {:5.1f}'
                  .format(vis, vis_tonight,vis_prompt,altmin))
            #log.prt(" Trigger: {}".format(self.t_trig.datetime))

            if (vis_tonight): # Seen within 24hr after the trigger
                log.prt("+----------------------------------------------------------------+")

                ###------------------
                def show(t,case="Dummy"):
                    for elt in t:
                        t1 = elt[0]
                        t2 = elt[1]
                        log.prt(" {:6s} : {} * {}".format(case,
                                                          t1.datetime,
                                                          t2.datetime))
                        t1  = (t1-self.t_trig).sec*u.s
                        t2  = (t2-self.t_trig).sec*u.s
                        # log.prt("        : {:7.2f} {:6.2f} * {:7.2f} {:6.2f}"
                        #       .format(t1,self.altaz(loc=loc,dt=t1).alt,
                        #               t2,self.altaz(loc=loc,dt=t2).alt))
                    return
                #-------------------

                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    show(t_event,case="Event") # Event - above altmin
                    show(t_twilight,case="Twil.")  # Twilight - Dark time
                    show(t_true,case="True")  # True : dark time + altmin + triggered
                log.prt("+----------------------------------------------------------------+")

        return

#------------------------------------------------------------------------------
if __name__ == "__main__":
    """
    A standalone functionto read a GRB and make various tests
    """
    import os
    from utilities import Log
    os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
    import pickle

    import grb_plot as gplt
    import gammapy
    from SoHAPPy import get_grb_fromfile, init

    if (gammapy.__version__ == "0.12"):
        from   gammapy.spectrum.models import Absorption
    if (gammapy.__version__ == "0.16"):
        from gammapy.modeling.models import Absorption

    import ana_config as cf # Steering parameters

    cf.dbg  = 0
    cf.day_after  = 0
    #cf.altmin = 10.2*u.deg
    cf.ngrb = 1 # 250
    cf.niter = 10# 110
    cf.ifirst = [85]


    log_filename    = cf.res_dir  + cf.logfile
    log = Log(name  = log_filename, talk=not cf.silent)

     # GRB list to be analysed
    if type(cf.ifirst)!=list:
        grblist = list(range(cf.ifirst,cf.ifirst+cf.ngrb))
    else:
        grblist = cf.ifirst


    out = None
    if (cf.altmin == 10*u.deg): out = open("deltas.txt", 'w')

    # Loop over GRB list accordng to the config. file
    for i in grblist:

        # Get GRBs
        init("")
        cf.niter = 10 # to avoid debugging when checking visibility
        grb = get_grb_fromfile(i,log=log)
        if (cf.newvis):
            grb.update_visibility(altmin=cf.altmin, debug=True)
        print(grb)
        grb.show_visibility(log=log)

    if out: out.close()


