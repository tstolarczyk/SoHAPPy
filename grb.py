# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import yaml
from yaml import CLoader as Loader
#from glob import glob

from pathlib import Path

import numpy as np

import astropy.units         as u
from   astropy.table         import Table
from   astropy.io            import fits
from   astropy.time          import Time
from   astropy.coordinates   import SkyCoord, AltAz, EarthLocation

from gammapy.maps import WcsGeom, MapAxis
from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel
from gammapy.image.models import SkyPointSource
from gammapy.spectrum import SpectrumDatasetOnOffStacker
from gammapy.data import Observations
from gammapy.stats.poisson import significance_on_off

# Prevent downloading Earth position correction from outside
# Results in warning
# https://docs.astropy.org/en/stable/utils/iers.html
from astropy.utils import iers
iers.conf.auto_download = False

__all__=['GammaRayBurst','GRBTarget']

###############################################################################
class GammaRayBurst(object):
    """
    Class to store GRB properties

    This class includes the calls for simulated observations.
    The GRB information is read from file.
    A GRB is composed of :
        - a name
        - a redshift
        - a position in the sky (temporarily given outside the file)
        - a list of time-intervals with correspoding :
            - energy spectrum.
            - observations
            - simulations

    Simulations are done with appropriate functions
    ~run_OnOffsimulation
    ~run_3Dsimulation
    """

    ###########################################################################
    def __init__(self):
        """
        The default visibility has been put by hand for backward compatibility
        with the olf file format and defintions where it was assumed that the
        GRB was seen without any delay. The dummy GRB is visible in North and
        South (inspired form GRB 289)
        """
        # Path of the GRB properties/model
        self.filepath = 'dummy'
        self.name     = 'dummy'
        self.reduction_factor = 1

        # GRB properties - theory
        self.z        = 0*u.dimensionless_unscaled
        self.radec    = SkyCoord(ra=100*u.degree,
                                 dec=-20*u.degree,
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
        self.altmin   = 90*u.deg
        self.t_trig   = Time('2000-01-01 00:00:00', scale='utc')
        self.site       = {"North": 'Roque de los Muchachos',
                           "South": 'Paranal'}

        # Visible any moment of the year from the site
        self.vis         = {"North":True,"South":True}
        # Visible any moment within the 24h after the trigger from the site
        self.vis_tonight = {"North":True,"South":True}
        # Visible at the moment of the GRB trigger
        self.vis_prompt  = {"North":True,"South":True}

        # The start and stop time are searched during 24h after the trigger
        # and correspond to the first visibility interval.
        # Does not report a second visibility interval during 24h,
        # that should be possible only if the target is promptly visible
        # (about 15% of the targets)
        # By default seen in the North

        # Visible above min. alt.
        self.t_true     = {"North":[Time('2000-01-01 00:01:00', scale='utc'),
                                    Time('2000-01-01 12:00:00', scale='utc')],
                           "South":[Time('2000-01-01 03:00:00', scale='utc'),
                                    Time('2000-01-01 12:00:00', scale='utc')]}

        # Astronomical twilight
        self.t_twilight = {"North":[Time('2000-01-01 00:01:00', scale='utc'),
                                    Time('2000-01-01 12:00:00', scale='utc')],
                           "South":[Time('2000-01-01 03:00:00', scale='utc'),
                                    Time('2000-01-01 12:00:00', scale='utc')]}

        # Rise and set of the target
        self.t_event    = {"North":[Time('2000-01-01 00:01:00', scale='utc'),
                                    Time('2000-01-01 12:00:00', scale='utc')],
                           "South":[Time('2000-01-01 03:00:00', scale='utc'),
                                    Time('2000-01-01 12:00:00', scale='utc')]}

        # Flux table - Flux at a series of points
        self.Eval           = [0]*u.GeV
        self.tval           = [0]*u.s
        self.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")

        # Valid time intervals, and observing times (end of intervals)
        self.time_interval  = [0]*u.s
        self.t_s            = [0]*u.s

        # GRB spectral and spatial model
        self.spectra = [] # Gammapy models (one per t slice)
        self.pointing = None
        self.spatial  =  SkyPointSource(lon_0=0*u.deg,lat_0=0*u.deg)

        # Stacked obs
        self.stack_obs = None

        return

    ###########################################################################
    @classmethod
    def read(cls, filename, ebl= None):

        """
        Fluxes are given for a series of (t,E) values
        From this we create a list of time intervals omitting the last
        point in t. An alternative is to add and other t point close to
        infinite, with a very low flux which is unlikely to give
        significant detection.

        Once this is done, it is considered that the flux at is valid
        for all t values until the end of the interval. An alternative
        would be to compute an evolutive flux along the time interval, or
        to assign a mean value from both edges to the whole interval.

        The observing time is considered to be the end of the interval.

        - number of values in t  : $N_{val}$
        - number of intervals    : $N_{val} - 1
        - number of observations : $N_{val} - 1

        Parameters
        ----------
        cls : GammaRayBurst class
            GammaRayBurst class instance
        filename : STRING
            GRB input file name
        ebl : STRING, optional
            The EBL absorption model considered. The default is 'dominguez'.

        Returns
        -------
        A GammaRayBurst filled instantiation.

        """

        grb = cls() # This calls the constructor

        hdul = fits.open(filename)

        grb.name = Path(filename.name).stem

        # GRB properties
        hdr = hdul[0].header
        grb.radec     = SkyCoord(ra=hdr['RA']*u.degree,
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

        # GRB alert received
        grb.t_trig      = Time(hdr['GRBJD']*u.day,format="jd")
        grb.vis         = {"North":hdr['V_N_ANYT'],"South":hdr['V_S_ANYT']}
        grb.vis_tonight = {"North":hdr['V_N_TNG'], "South":hdr['V_S_TNG']}
        grb.vis_prompt  = {"North":hdr['V_N_PR'],  "South":hdr['V_S_PR']}

        grbvis = Table.read(hdul,hdu=1)

        grb.altmin = grbvis.meta["MIN_ALT"]*u.deg # Minimum altitude

        # Start/stop times
        t  = Time(grbvis["True"].data,format="jd")
        grb.t_true     = {"North":t[0:2],"South":t[2:4]}

        t  = Time(grbvis["Twilight"].data,format="jd")
        grb.t_twilight = {"North":t[0:2],"South":t[2:4]}

        t  = Time(grbvis["Event"].data,format="jd")
        grb.t_event    = {"North":t[0:2],"South":t[2:4]}

        # Energies, times and fluxes
        grb.Eval     = Table.read(hdul,hdu=2)["Energies (afterglow)"].quantity
        grb.tval     = Table.read(hdul,hdu=3)["Times (afterglow)"].quantity
        flux         = Table.read(hdul,hdu=4)
        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Warning : the units here should comply with
        # hdulist[4].header["units"]
        # Note the transpoisition from flux to fluxval
        grb.fluxval = np.zeros( (icol_t,jrow_E) ) /(u.cm)**2/(u.s)/(u.GeV)
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                grb.fluxval[i][j] = flux[j][i]/((u.cm)**2*u.s*u.GeV) # transp!

        # Compute time intervals
        grb.time_interval = []
        for i,t in enumerate(grb.tval[:-1]): # Forget last point
            grb.time_interval = grb.time_interval + [[t,grb.tval[i+1]]]

        # Set observing time to end of interval
        grb.t_s = grb.tval[1:]

        # Create spectral model for each time interval
        for i,t in enumerate(grb.time_interval):
            # We decide that the valid flux is the one at the start of the
            # time interval
            tab = TableModel(energy       = grb.Eval,
                             values       = grb.fluxval[i],
                             norm         = 1.,
                             values_scale = 'lin')

            grb.spectra.append(AbsorbedSpectralModel(spectral_model = tab,
                                                      absorption    = ebl,
                                                      parameter     = grb.z))

        return grb

    ###########################################################################
    @classmethod
    def read_old(cls, path, ebl = None,reduction_factor=1,
             pointing = None, dbg=False):

        """

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
        A GammaRayBurst filled instantiation.

        """

        """
        This method reads one of the 10 old GRB test files.
        In this temporary files :
         - The GRB position is not given
         - The GRB occurence time is not given

        The data consist in time intervals for which an energy specturm is
        given. The time intervals are set by hand because of a file format
        that do not permit to retrieve them.

        The observing time is set to the ends of the time intervals.

        In order to be consistent with the new file format the tval values
        are recomputed from the time interval boundaries.
        A flux is added at the last time point which has the value of the
        flux in the last interval.
        This is then coherent with the new file format and in particular :

        - number of values in t  : $N_{val}$
        - number of intervals    : $N_{val} - 1
        - number of observations : $N_{val} - 1

        """

        grb = cls() # This calls the constructor

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
        fluxlist = [] # new

        ## grb.fluxval = np.zeros( (icol_t,jrow_E) ) /(u.cm)**2/(u.s)/(u.GeV)
        #HERE
        for t in dt:
            filename = '{}_{:.0f}-{:.0f}.txt'.format(data["name"],t[0],t[1])
            table = Table.read(path.joinpath(filename), format='ascii')

            energy = table['col1'] * u.TeV
            flux = (table['col2'] / reduction_factor) * u.Unit('1 / (cm2 s TeV)')
            fluxlist.append(flux) # Each flux contains len(energy) elts

            table_model = TableModel(energy      = energy,
                                     values      = flux,
                                     norm        = 1.,
                                     values_scale= 'lin')
            grb.spectra.append(
                AbsorbedSpectralModel(spectral_model = table_model,
                                      absorption     = ebl,
                                      parameter      = grb.z))

        fluxlist.append(fluxlist[-1]) # Add the last flux to the last time

        # The observing times are considered to be the end of the interval
        grb.t_s = [x[1] for x in dt] *u.s

        # Re-create t_val values
        grb.tval     = ( [dt[0][0]] + ([x[1] for x in dt] ) ) * u.s
        grb.Eval     = energy
        grb.fluxval  = np.asarray(fluxlist) \
                       *fluxlist[0][0].unit # asarray remove units !
        grb.time_interval = dt*u.s

        grb.pointing = pointing
        return grb

    ###########################################################################
    def plot(self):
        """
        Display GRB characteristics
        """
        print("READY TO PLOT")

        return

    ###########################################################################
    def __str__(self):
        """ Printout the GRB properties """
        txt = '\n'
        txt += "=============================================================\n"
        txt += "===                             GRB summary               ===\n"
        txt += "=============================================================\n"
        txt += '  Name           : {}\n'.format(self.name)
        txt += '  RA             : {}\n'.format(self.radec.ra)
        txt += '  DEC            : {}\n'.format(self.radec.dec)
        txt += '  Redshift       : {}\n'.format(self.z)
        txt += '  Eiso           : {}\n'.format(self.Eiso)
        txt += '  Epeak          : {}\n'.format(self.Epeak)
        txt += '  t90            : {}\n'.format(self.t90)
        txt += '  G0H / G0W      : {} / {}\n'.format(self.G0H,self.G0W)
        txt += '  Flux peak      : {}\n'.format(self.FluxPeak)
        txt += '  gamma LE / HE  : {} / {}\n'.format(self.gamma_le,self.gamma_he)
        txt += '  t_trig         : {}\n'.format(self.t_trig.datetime)

        txt += '  Energy bins    : {}\n'.format(len(self.Eval))
        txt += '  Time bins      : {}\n'.format(len(self.tval))
        txt += '  Time intervals : {}\n'.format(len(self.time_interval))

        return txt

    ###########################################################################
    def add_stack_OnOffobs(self,simulations):
        """Stack observations"""
        stack = SpectrumDatasetOnOffStacker(Observations(simulations))
        #stack.run()
        return stack.stacked_obs

    ###########################################################################
    def get_cumulative_stats(self,simulations):
        """Get cumulative statistics"""

        # Init vectors
        nsim = len(simulations)
        tot_time  = np.zeros(nsim)
        tot_n_on  = np.zeros(nsim)
        tot_n_off = np.zeros(nsim)
        tot_alpha = np.zeros(nsim)
        delta_t   = np.zeros(nsim)

        # Loop on observations
        for idx, obs in enumerate(simulations):

            alpha = obs.alpha

            if idx == 0:
                tot_time[idx]  = obs.total_stats_safe_range.livetime.value
                tot_n_on[idx]  = obs.total_stats_safe_range.n_on
                tot_n_off[idx] = obs.total_stats_safe_range.n_off
                tot_alpha[idx] = obs.total_stats_safe_range.alpha
            else:
                tot_time[idx]  += tot_time[idx-1] + obs.total_stats_safe_range.livetime.value
                tot_n_on[idx]  += tot_n_on[idx-1] + obs.total_stats_safe_range.n_on
                tot_n_off[idx] += tot_n_off[idx-1] + obs.total_stats_safe_range.n_off
                tot_alpha[idx]  = obs.total_stats_safe_range.alpha

            delta_t[idx] =  obs.total_stats_safe_range.livetime.value
            tot_excess = tot_n_on - alpha * tot_n_off
            tot_bkg = alpha * tot_n_off
            tot_significance = significance_on_off(tot_n_on,
                                                   tot_n_off,
                                                   tot_alpha)
#        print(" cumulated stat : On={:10.2f} Off={:10.2f} bck={:10.2f} dt={:10.2f} "
#           .format(tot_n_on[nsim-1], 
#                   tot_bkg[nsim-1],
#                   tot_n_off[nsim-1],
#                   tot_time[nsim-1])) 
#        print(" Cumulated On   = ",tot_n_on)
#        print(" Cumulated bkg  = ",tot_bkg)
#        print(" Cumulated Off  = ",tot_n_off)
#        print(" Cumulated time = ",tot_time) 

        return dict(livetime = tot_time,
                    excess   = tot_excess,
                    bkg      = tot_bkg,
                    sigma    = tot_significance,
                    n_on     = tot_n_on,
                    n_off    = tot_n_off,
                    alpha    = tot_alpha,
                    delta_t  = delta_t)

    # ###########################################################################
    # def quicklook(self,simulations,plot=False):

    #     for idx, obs in enumerate(simulations):
    #         print("GRB ",self.name," - Observation: ",idx)
    #         #obs.peek() # n_on, alpha*n_off Espec. eff.area,, Ematrix and stat
    #         print("    - lifetime:",obs.total_stats_safe_range.livetime.value)
    #         print("    - excess vector:",obs.excess_vector.data)
    #         if (plot):
    #             obs.excess_vector.peek()
    #         plt.show()

    #     return

    @staticmethod
    ###########################################################################
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

###############################################################################
#
###############################################################################
class GRBTarget(object):
    """
    Observation target information.

    Parameters
    ----------
    name : `str`
        Name of the source
    model : `~gammapy.spectrum.models.SpectralModel`
        Model of the source
    """

    ###########################################################################
    def __init__(self, name=None, model=None,pointing = None):
        self.name = name
        self.model = model
        self.pointing = pointing
    ###########################################################################
    def __str__(self):
        """Target report (`str`)."""
        ss = '*** Target parameters ***\n'
        ss += 'Name={}\n'.format(self.name)
        for par in self.model.parameters.parameters:
            ss += '{}={} {}\n'.format(par.name, str(par.value), par.unit)
        return ss
    ###########################################################################
    def from_fermi_lat_catalogue(name):
        raise NotImplementedError

###############################################################################
#
###############################################################################
class GRBObservationParameters(object):
    """
    Container for observation parameters.

    Parameters:

    alpha : `~astropy.units.Quantity`
        Normalisation between ON and OFF regions

    livetime :  `~astropy.units.Quantity`
        Observation time

    emin :  `~astropy.units.Quantity`
        Minimal energy for simulation

    emax : `~astropy.units.Quantity`
        Maximal energy for simulation
    """

    ###########################################################################
    def __init__(self, alpha=None, livetime=None,
                 emin=None, emax=None, pointing=None, fov=None, binsize=None):
        self.alpha = alpha
        self.livetime = livetime
        self.emin = emin
        self.emax = emax
        self.fov = fov
        self.binsize = binsize

        # Create sky map frame - Emin, Emax, and the bin size could be extracted from the IRF
        nedges = 10
        logEedges = np.logspace(np.log10(self.emin.value),
                                np.log10(self.emax.value), nedges)
        if (pointing):
            self.axis = MapAxis.from_edges(logEedges, unit="TeV",
                                           name="energy", interp="log")
            self.geom = WcsGeom.create(skydir   = pointing,
                                       binsz    = binsize,
                                       width    = (fov, fov),
                                       coordsys ="GAL",
                                       axes=[self.axis]) # Width in degrees

#        print(" - Number of Energy log. bins : ",nedges-1," from ",self.emin," to ",self.emax)
#        print(" - Number of spatial bins : ",self.fov/self.binsize," x ",self.fov/self.binsize)
#        print(" *** Sky map")

    ###########################################################################
    def __str__(self):
        """Observation summary report (`str`)."""
        ss = '*** Observation parameters summary ***\n'
        ss += 'alpha={} [{}]\n'.format(self.alpha.value, self.alpha.unit)
        ss += 'livetime={} [{}]\n'.format(self.livetime.value,
                                          self.livetime.unit)
        ss += 'emin={} [{}]\n'.format(self.emin.value, self.emin.unit)
        ss += 'emax={} [{}]\n'.format(self.emax.value, self.emax.unit)
        return ss

#------------------------------------------------------------------------------
if __name__ == "__main__":

    import grbplot as gplt
    from gammapy.spectrum.models import Absorption
    import ana_config as cf # Steering parameters

    cf.old_file   = False # Use old GRB file
    cf.dbg_level  = 0
    absorption = Absorption.read_builtin(cf.EBLmodel)
    for i in range(cf.ifirst,cf.ngrb+cf.ifirst):
        if (cf.old_file):
            loc = Path(cf.grb_oldfolder + "LGRB"+str(i))
            grb = GammaRayBurst.read_old(loc,ebl=absorption)
        else:
            loc = Path(cf.grb_folder + "/Event"+str(i)+".fits")
            grb = GammaRayBurst.read(loc,ebl=absorption)

        print(grb)

        tmin = grb.t_trig # GRB trigger time
        tmax = grb.t_event["North"][1] # End of visibility
        duration = tmax - tmin
        print(duration)

        # Visibility points
        tvis = np.linspace(0,duration.value,100)
        sitecoord = EarthLocation.of_site(grb.site["North"])
        altazvis  = grb.radec.transform_to( AltAz(obstime= tvis + grb.t_trig,
                                        location=sitecoord))

        # GRB points with the visibility
        tgrb = grb.t_s
        tgrb = [ t.value for t in tgrb if t < duration]*tgrb[0].unit
        altazgrb = grb.radec.transform_to( AltAz(obstime= grb.t_trig+tgrb,
                                        location=sitecoord))
#        # Measurement points
#
        from   astropy.visualization import quantity_support
        import matplotlib.pyplot as plt
        with quantity_support():
            plt.plot(tvis,altazvis.alt,label="Visibility")
            plt.plot(tgrb,altazgrb.alt,marker='X',ms=10.0,ls='',label="Obs. point")
            plt.legend()
            plt.show()
        #gplt.visibility(grb,opt='plot')
        #gplt.spectra(grb,opt="2D")
        #gplt.spectra(grb,opt="Packed")
        #gplt.spectra(grb,opt="Time")
        #gplt.spectra(grb,opt="Energy")



