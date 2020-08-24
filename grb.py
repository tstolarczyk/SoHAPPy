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
from   astropy.coordinates   import SkyCoord

from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel
from gammapy.image.models import SkyPointSource

# Prevent downloading Earth position correction from outside
# Results in a warning : https://docs.astropy.org/en/stable/utils/iers.html
from astropy.utils import iers
iers.conf.auto_download = False

__all__=['GammaRayBurst']

###############################################################################
class GammaRayBurst(object):
    """
    Class to store GRB properties

    The GRB information is read from files.
    A GRB is composed of several parameters among which:
        - A name
        - A redshift
        - A position in the sky
        - A visibility window
    It is characterised by a list of time-intervals with a given energy 
    spectrum.
    
    The energy spectrum is given as a series of points given the flux at that 
    point.
    
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
        self.radec    = SkyCoord(ra=100*u.degree, dec=-20*u.degree, 
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
        self.site       = {"North": 'Roque de los Muchachos',
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
        # By default seen in the North
        # Default times are defined for the firts GRB sample with no 
        # visibikity given. The North visibility is chosen to have the 
        # trigger therein. The South visibility starts after the trigger. 
        # This allows having the two conditions explored.

        # Visible above min. alt.
        self.t_true     = {"North":[Time('2000-01-01 00:00:00', scale='utc'),
                                    Time('2000-01-01 08:00:00', scale='utc')],
                           "South":[Time('2000-01-01 03:00:00', scale='utc'),
                                    Time('2000-01-01 11:00:00', scale='utc')]}

        # Astronomical twilight
        self.t_twilight = {"North":[Time('2000-01-01 00:00:00', scale='utc'),
                                    Time('2000-01-01 8:00:00', scale='utc')],
                           "South":[Time('2000-01-01 03:00:00', scale='utc'),
                                    Time('2000-01-01 11:00:00', scale='utc')]}

        # Rise and set of the target
        self.t_event    = {"North":[Time('2000-01-01 00:00:00', scale='utc'),
                                    Time('2000-01-01 08:00:00', scale='utc')],
                           "South":[Time('2000-01-01 03:00:00', scale='utc'),
                                    Time('2000-01-01 11:00:00', scale='utc')]}

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
        From this we create a list of time intervals omitting the first
        point in t. The flux does not exist beyond the last point. 
        
        Once this is done, it is considered that the flux at is valid
        for all t values until the end of the interval. An alternative
        would be to compute an evolutive flux along the time interval, or
        to assign a mean value from both edges to the whole interval.

        The observing time is considered to be the end of the interval.
        
        The spectrum is stored as a table, TableModel, that take as an input
        a series of flux values as a function of the energy.
        The model will return values interpolated in log-space with 
        math::`scipy.interpolate.interp1d`, returning zero for energies 
        outside of the limits of the provided energy array.

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
            The EBL absorption model considered. There is no default set.

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
        flux_unit    = u.Unit(flux.meta["UNITS"])/u.Unit("ph") # Removes ph
        
        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Note the transpoisition from flux to fluxval
        grb.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                grb.fluxval[i][j] = flux[j][i]*flux_unit # transp!

        # Compute time intervals
        grb.time_interval = []
        for i,t in enumerate(grb.tval[:-1]): # Forget last point
            grb.time_interval = grb.time_interval + [[t,grb.tval[i+1]]]

        # Set observing time to end of intervals
        grb.t_s = grb.tval[1:]

        # Create spectral model for each time interval
        for i,t in enumerate(grb.time_interval):
            # We decide that the valid flux is the one at the start of the
            # time interval
            # Note that TableModel makes an interpolation
            tab = TableModel(energy       = grb.Eval,
                             values       = grb.fluxval[i],
                             norm         = 1.,
                             values_scale = 'lin')

            if (ebl != None):
                grb.spectra.append(AbsorbedSpectralModel(spectral_model = tab,
                                                      absorption    = ebl,
                                                      parameter     = grb.z))
            else:
                print(" EBL absorption not defined")
                
        return grb

    ###########################################################################
    @classmethod
    def read_old(cls, path, ebl = None,reduction_factor=1,
             pointing = None):

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
    def __str__(self):
        """ Printout the GRB properties """
        
        # After glow duration
        dt_afterglow = self.time_interval[-1][1]-self.time_interval[0][0]
        txt = '\n'
        txt += "==================================================================\n"
        txt += "===                             {:10s}                     ===\n"\
        .format(self.name)
        txt += "==================================================================\n"
        txt += '  RA, DEC            : {:6.2f} {:6.2f}\n'\
        .format(self.radec.ra, self.radec.dec)
        txt += '  Redshift           : {:6.2f}\n'.format(self.z)
        txt += '  Eiso               : {:6.2e}\n'.format(self.Eiso)
        txt += '  Epeak              : {:6.2f}\n'.format(self.Epeak)
        txt += '  t90                : {:6.2f}\n'.format(self.t90)
        txt += '  G0H / G0W          : {:6.2f} / {:6.2f}\n' \
        .format(self.G0H,self.G0W)
        txt += '  Flux peak          : {:6.2f}\n'.format(self.FluxPeak)
        txt += '  gamma LE / HE      : {:6.2f} / {:6.2f}\n' \
        .format(self.gamma_le,self.gamma_he)
        txt += '  t_trig             : {}\n'.format(self.t_trig.datetime)
        txt += '  Afterglow duration : {:6.2f}\n'.format(dt_afterglow.to(u.d))
        txt += '  Bins : E, t, dt    : {:4d} {:4d} {:4d} \n' \
        .format(len(self.Eval), len(self.tval), len(self.time_interval))
        

        return txt
###########################################################################
    def visibility(self,opt="None"):
    
        site  = ["North","South"]
        vis   = {"North":False,"South":False}
        
        print("  Visibility (Minimum altitude : {})".format(self.altmin))
    
        for loc in site:
            print("+----------------------------------------------------------------+")
            print('{}  Visible    : {} - tonight, prompt : {}, {}'
                  .format(loc,
                          self.vis[loc],
                          self.vis_tonight[loc],
                          self.vis_prompt[loc]))
            if (self.vis_tonight[loc]): # Seen within 24hr after the trigger
                vis[loc] = True
                print(" Event  : {} * {}"
                      .format(self.t_event[loc][0].datetime,
                              self.t_event[loc][1].datetime))
    
                print(" Twil.  : {} * {}"
                      .format(self.t_twilight[loc][0].datetime,
                              self.t_twilight[loc][1].datetime))
    
                print(" True   : {} * {}"
                      .format(self.t_true[loc][0].datetime,
                              self.t_true[loc][1].datetime))
        print("+----------------------------------------------------------------+")
        
        return (vis["North"] or vis["South"])
    
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
#
################################################################################
##
################################################################################
#class GRBTarget(object):
#    """
#    Observation target information.
#
#    Parameters
#    ----------
#    name : `str`
#        Name of the source
#    model : `~gammapy.spectrum.models.SpectralModel`
#        Model of the source
#    """
#
#    ###########################################################################
#    def __init__(self, name=None, model=None,pointing = None):
#        self.name = name
#        self.model = model
#        self.pointing = pointing
#    ###########################################################################
#    def __str__(self):
#        """Target report (`str`)."""
#        ss = '*** Target parameters ***\n'
#        ss += 'Name={}\n'.format(self.name)
#        for par in self.model.parameters.parameters:
#            ss += '{}={} {}\n'.format(par.name, str(par.value), par.unit)
#        return ss
#    ###########################################################################
#    def from_fermi_lat_catalogue(name):
#        raise NotImplementedError

#------------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
    
    import grb_plot as gplt
    from gammapy.spectrum.models import Absorption
    import ana_config as cf # Steering parameters

    cf.old_file   = False # Use old GRB file
    cf.dbg_level  = 0
    absorption = Absorption.read_builtin(cf.EBLmodel)
    
    for i in range(cf.ifirst,cf.ngrb+cf.ifirst):

        # Get GRBs
        if (cf.old_file):
            loc = Path(cf.grb_olddir + "LGRB"+str(i))
            grb = GammaRayBurst.read_old(loc,ebl=absorption)
        else:
            loc = Path(cf.grb_dir + "/Event"+str(i)+".fits")
            grb = GammaRayBurst.read(loc,ebl=absorption)

        print(grb)
        gplt.spectra(grb,opt="Packed")

        grb.visibility()
        gplt.visibility_chart(grb)
        gplt.visibility_alt_plot(grb)

        ###################################
        # Investigate next day visibility
        ###################################
        # If the GRB was not visible, no chance it will be the day after
        # -9 values a unchanged
        for loc in ["North","South"]:
            for i in range(2):
                
                # The time of visiility is when the GRB is above the horizon 
                # (at min altitude the same) theday after. However, contrarily to 
                # the first day, the trigger has occured.
                if (grb.t_true[loc][i] != -9): 
                    
                    grb.t_true[loc][i] += 1*u.d
                
                # The dark time is the same + approx. one day
                if (grb.t_twilight[loc][i] != -9): 
                    grb.t_twilight[loc][i] += 1*u.d
                
                # The time above the horizon is the same + approx. one day
                if (grb.t_event[loc][i] != -9): 
                    grb.t_event[loc][i] += 1*u.d   
                
        loc = "North"
        print("the day after :")
        print(" Visible : unchanged -> ",grb.vis[loc])
        print(" Tonight : unchanged -> ",grb.vis_tonight[loc])
        print(" Prompt : unchanged -> "
             ,grb.vis_tonight[loc])

#        from   astropy.visualization import quantity_support
#        import matplotlib.pyplot as plt
#        with quantity_support():
#            plt.plot(tvis,altazvis.alt,label="Visibility")
#            plt.plot(tgrb,altazgrb.alt,marker='X',ms=10.0,ls='',label="Obs. point")
#            plt.legend()
#            plt.show()
        #gplt.spectra(grb,opt="2D")
        #gplt.spectra(grb,opt="Packed")
        #gplt.spectra(grb,opt="Time")
        #gplt.spectra(grb,opt="Energy")



