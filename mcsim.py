# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""

import sys
import numpy as np
import os
import time

import concurrent.futures

import astropy.units as u
from   astropy.coordinates   import AltAz

#from scipy.stats import chi2
#from scipy.stats import norm

from irf import IRF
from irf_onaxis import CTAPerf_onaxis

from gammapy.irf import load_cta_irfs
from gammapy.spectrum.models import PowerLaw
from gammapy.maps import WcsNDMap

from gammapy.cube import MapDataset
from gammapy.cube.models import SkyModel, BackgroundModel
from gammapy.spectrum import SpectrumDatasetOnOffStacker
from gammapy.data import Observations

from fit_onoff import mc_onoff, cumulative_OnOffstats
from fit_3d import fit_3d, DataSetShow

import mcsim_config as mcf
import ana_config as cf

from utilities import banner, warning

__all__ = ['MonteCarlo']
###############################################################################
class MonteCarlo():
    """
    This class handles a GRB simulation, either on one of the two sites, and
    in a next version combining both sites.
    It uses several additionnal modules :
        - mcsim_config : contains variables and flags with some useful options 
        that are not intended to be changed by the user (suggested values are 
        mentionned between brackets):
            - det_level (0.9): 
                Fraction declaring a detection above the 3,5 sigma threshold
            - alpha (0.2) : related to the number of off N zones in the 
            on-off analysis, alpha = 1/N
            - fov (5.*u.deg) : Full Field-of-view for the sky map
            - pxsize (0.125*u.deg) : square pixel length for the sky map
            - dtslew (30*u.s) : slewing time (supposed identival for all 
            telescope types)  
            - fixslew (True) : if True teh slewing time if fixed to dtslew,
            or randomly generated between 0 and dtslew if False.
            - n_lima_min = 10 # Count number below which Li&Ma cannot be trusted anymore
write_slices = True # Store detailed information on slices if True


from   astropy.coordinates   import EarthLocation
# The following information can be obtained from EarthLocation.of_site
# but this requires to have an internet connection
# Get all available site names :
# sitelist = EarthLocation.get_site_names()

# Get coordinattes of Paranal and LaPalma

# xyz_north = EarthLocation.of_site('Roque de los Muchachos')
xyz_north = EarthLocation.from_geocentric( 5327448.9957829,
                                          -1718665.73869569,
                                           3051566.90295403,
                                           unit="m")

# xyz_south = EarthLocation.of_site('Paranal Observatory')
xyz_south = EarthLocation.from_geocentric( 1946404.34103884,
                                          -5467644.29079852,
                                          -2642728.20144425,
                                          unit="m")

pos_site ={"North":xyz_north, "South":xyz_south} 
        - mcsim_plot : generate plot from the simulation results.
    """
    ###########################################################################
    def __init__(self,grb,
                 niter  =1,
                 where  = "North",
                 ndof   = 0,
                 debug  = 0):
        """
        Initialize class members to default values

        """
        zdeg = 0*u.deg
        
        self.dbg          = debug
        
        # Input parameters and objects
        self.grb          = grb    # A GRB instantiation
        self.niter        = niter     # Number of trials
        self.ndof         = ndof   # Deg. of freedom in likelikood, 0 if on-off
        self.where        = where  # Site ("North" or "South")
        
        # Create performance for each Observation time slice
        self.ntbins       = 0  # Number of time slices in the simulation
        self.tstart       = 0 *u.s # Observation start after trigger
        self.tstop        = 0 *u.s # Observation stop after trigger
        self.pos_start    = AltAz(alt=zdeg,az=zdeg) # Alt-az at tstart
        self.pos_stop     = AltAz(alt=zdeg,az=zdeg) # Att-az at tstop
        self.tobs         = [] # Observation time slices from grb intervals
        self.altaz        = [] # Position in the sky at start/end of interval
        self.perf         = [] # List of performances for each time slice
        self.Emin         = [] # Min energy from the IRF
        self.Emax         = [] # Max energy from the IRF

        # Initialise mode associated to the MC
        self.rnd_seed     = 'random-seed'
        self.source_model = self.source_spectrum()
        self.null_model   = self.source_spectrum(norm=0/(u.cm)**2/u.s/u.TeV)

        # list of simulations (one per slice)
        self.simulations = None

        # Significance over simulations
        self.id_smax_list  = [] # Slice number to get back the time/ altaz
        self.smax_list     = [] # List of max. significances along trials
        self.nex_smax_list = [] # List of excess counts at max. signif. 
        self.nb_smax_list  = [] # List of background counts at max. signif.

        self.id_3s_list    = [] # Slice number ot get back the time/altaz
        self.nex_3s_list   = [] # List of excess
        self.nb_3s_list    = [] # List of background
        self.detect_3s     = 0  # Number of trials 3 sigma was reached

        self.id_5s_list    = [] # Slice number to get back the time/altaz
        self.nex_5s_list   = [] # List of excess
        self.nb_5s_list    = [] # List of background
        self.detect_5s     = 0  # Number of trials 5 sigma was reached

        # Mean sigma versus time - one value per time slice
        self.sigma_mean = []
        self.sigma_std  = []
        
        # Slice number with error or warning
        self.err_slice  = [] 

        self.mctime = 0.00
        self.err    = -999 # error code : default, simulation is not completed 
            
        return

###############################################################################
#
# Monte Carlo simulation (loop over slices)
#
###############################################################################
    def run(self):
        """
        - Run simulations of the current grb, for a given hemisphere.
          A simulation correspond to a series of observations corresponding
          to the GRB time slices.

        - Compute and store :
            - Significance values for all slices in the trials;
            - Maximum significance reached for the GRB;
            - Time of maximum siginificance, number of excess and background 
            events
            - Time to reach 3 or 5 sigma, number of excess and background 
            events
            
        """
                 
        print()
        banner("+==============+=================================================+")
        banner("| SIMULATION   | {:>10s} - {:<10s}                         |"
                  .format(self.grb.name,self.where))
        
        if (self.grb.vis_tonight[self.where] != True):
            return # mc variables remain at their default (see Constructor)

        if (self.observation_slices() == 0):  # Create observation time slices
            warning(" No observatioon slices found")
            return
        
        self.create_perf(folder=cf.irf_dir) # Get performance for each time slice

        sigma_sum   = 0           # Compute mean and std of sig for slices
        sigma2_sum  = 0           #     "          "
        self.mctime = time.time() # Starts chronometer
        if (cf.do_accelerate == False):
            abort_test  = True   # Abortion test supposed already performed
        else:
            abort_test  = False  # Abortion test not yet performed
        self.err    = self.niter # GRB is visible, simul. ought to be complete
        
        if (cf.write_slices): # File for problem tracking - deleted if none
            sdump = cf.res_dir +"/"+ self.grb.name+"-"+self.where+"_slices.txt"
            fslice  = open(sdump,"w")
            print("   ---",sdump," opened")
            to_be_dumped = False        
            self.dump_slices(header=True, file=fslice)
            
        ###############################################
        ### Monte Carlo iterations
        ###############################################
        print(" *** Iterating")
        iMC=1
        
        while(iMC <= self.niter):
            if (iMC <= 10) or (np.mod(iMC,10) == 0): print("#",iMC," ",end="")
            
            self.loop_slices_onoff() # Loop on GRB time slices - get the signal

            obs_stat = cumulative_OnOffstats(self.simulations, 
                                             n_min = mcf.n_lima_min,
                                             alpha = mcf.alpha)
            sigma    = obs_stat['sigma']
            non      = obs_stat['n_on']
            noff     = obs_stat['n_off']
            nex      = non - mcf.alpha*noff
            nb       = mcf.alpha*noff
            
            if (cf.write_slices): # dump slices to track problems
                status = self.dump_slices(iMC=iMC,
                                          non=non,noff=noff,sigma=sigma,
                                          file=fslice)
                if (to_be_dumped == False): to_be_dumped = status
                
            if (self.dbg > 2): 
                from mcsim_plot import onetrial
                onetrial(self)

#            if (self.ndof != 0): ################################ 3D
#                fcnsb = np.array([x[0] for x in self.simulations]) # lik_fit
#                fcnb  = np.array([x[1] for x in self.simulations]) # lik_bonly
#                nex   = np.array([x[2] for x in self.simulations])
#                nb    = np.array([x[3] for x in self.simulations])
#                # ns_bonly, nb_bonly - not used
#                TS =  fcnb - fcnsb
#                TS[TS < 1e-3] = 0 # Avoid sligtly negative values
#
#                # pvalue  = 1 - chi2.cdf(TS, ndof)
#                # sigma   = norm.ppf(1-pvalue)  # Integral from -Inf, can be<0
#                # print(" Warning sqrt TS")
#                sigma = np.sqrt(TS)
#                sigma[sigma < 0 ] = 0
#                sigma[sigma > 99] = 99

            # Acummulate sum and sum**2 for mean / error in each slice
            sigma_sum  += sigma
            sigma2_sum += sigma**2

            self.fill_stat(sigma, nex, nb) # Update stat list

            # In case 3 signma is not reached in the first 10% of the trials
            # then the 90% CL can not be reached.
            if (abort_test==False):
                if (iMC/self.niter > 1 - mcf.det_level):
                    abort_test = True
                    if (self.detect_3s ==0):
                        print("\n *** Will never reach ",
                              100*mcf.det_level,"% of CL",end="")
                        print(" ===> Simulation stopped")
                        self.err = iMC
                        break

            iMC+=1 # End of MC loop
            
        if (cf.write_slices) :
            print()
            fslice.close()
            print("   ---",sdump," closed")
            if (not to_be_dumped): 
                print("   ---",sdump," deleted")
                os.remove(sdump)                 
        
        self.mctime = (time.time() - self.mctime)/self.niter

        ### Compute mean sigma and rms for all trials
        self.sigma_mean = sigma_sum/self.niter
        sigma2_mean     = sigma2_sum/self.niter
        self.sigma_std  = np.sqrt(sigma2_mean-self.sigma_mean**2)
        
        return
    ##########################################################################
    def fill_stat(self,sigma, nex, nb):
        """
        Get the statistics and hnadle the eceptions of the current MC 
        simulation
        """
 
        ### Find maximum - cannot be a slice with non or noff below limit
        sigmax = np.nanmax(sigma) # Returns Nan only if all are Nan
        if (sigmax == np.nan):
            print(" All sigma values are Nan !!! ")
            maxidx = -1
            sigmax = -999
            nexmax = -1
            nbmax  = -1
        else:
            maxidx = np.nanargmax(sigma) # If sigma is nan, it would be the max !
            nexmax = nex[maxidx]
            nbmax  = nb[maxidx]
            
        self.id_smax_list.append(maxidx)
        self.smax_list.append(sigmax)
        self.nex_smax_list.append(nexmax)
        self.nb_smax_list.append(nbmax)

        # Find where 3 sigma is reached
        nex_3s = -1
        nb_3s  = -1
        t_3s   = -1
        alt_3s = -1
        id_3s = np.where(sigma>=3)[0] # This ignore Nan values 
        if (np.size(id_3s) != 0):
            id_3s  = id_3s[0] # First one
            nex_3s = nex[id_3s]
            nb_3s  = nb[id_3s]
            t_3s   = self.tobs[id_3s][1] # end of slice
            alt_3s = self.altaz[id_3s][1].alt # end of slice
            self.id_3s_list.append(id_3s)
            self.detect_3s += 1

        self.nex_3s_list.append(nex_3s)
        self.nb_3s_list.append(nb_3s)

        # Find where 5 sigma is reached
        nex_5s = -1
        nb_5s  = -1
        t_5s   = -1
        alt_5s = -1      
        id_5s = np.where(sigma>=5)[0]  # This ignore Nan values 
        if (np.size(id_5s) != 0):
            id_5s  = id_5s[0] # First one
            nex_5s = nex[id_5s]
            nb_5s  = nb[id_5s]
            t_5s   = self.tobs[id_5s][1] # end of slice
            alt_5s = self.altaz[id_5s][1].alt # end of slice
            self.id_5s_list.append(id_5s)
            self.detect_5s += 1

        self.nex_5s_list.append(nex_5s)
        self.nb_5s_list.append(nb_5s)

        if (self.dbg>2):
            # Bug corrected on April 27th, 2020
            # The convention is that the observation time is the end of the 
            # time slice, i.e. "tobs[maxidx][0] -> tobs[maxidx][1]
            # also corrected in mcsim_res
            print(" >>> t_smax = {:5.2f}, ({:5.2f} sigma at slice {:2d})"
                  .format(self.tobs[maxidx][0],sigma[maxidx],maxidx+1))
            print("     t_3s / t_5s  at {:5.2f}/{:5.2f} <- alt = {:6.2f}/{:6.2f}"
                  .format(t_3s,t_5s,alt_3s,alt_5s))

        return

###########################################################################
#
# Initialisations
#
###########################################################################
    def observation_slices(self):
        """
        Create the time intervals for observation as the result of the
        GRB time interval and the visibility slot.
        Get the boundaries of the observation window (tstart, tstop).
        Get the sky position at the boundaries.
        
        2020/01/23 : OK
            - check the t00 trick, don't remeber why it is here
        """
        print(" *** Create observations slices...",end="\n")

        on  = False
        off = False
        t1  = self.grb.t_true[self.where][0] # Visibilty window for the site
        t2  = self.grb.t_true[self.where][1]
        if (cf.day_after != 0):
            # Note : The next day(s)), by defintion the GRB is visible 
            # (i.e. the trigger date is outside the detection window)
            t1 = max(self.grb.t_twilight[self.where][0] + cf.day_after*u.d, 
                     self.grb.t_event[self.where][0]    + cf.day_after*u.d)    
            t2 = t2 + cf.day_after*u.d
        
        # Define visibility start and stop in GRB trigger reference system
        # tstart can only be positive since t_true start at latest when
        # the GRB is seen
        # It is possible that the visibility falls into the first 
        # GRB time slice from 0 to t1.
        # But 0-t1 is not in the list.
        # Therefore the first avalaible time is at most the first GRB point 
        tstart = (t1 - self.grb.t_trig).sec * u.s
        t00 = self.grb.time_interval[0][0]
        if tstart < t00:
            tstart = t00
            # print("tstart =",self.tstart," reshifted to t00 = ",t00)
        # Add slewing time 
        if (cf.fixslew == True):
            tstart = tstart + cf.dtslew
        else:
            tstart = tstart + cf.dtslew*np.random.random()
        
        # End of observation is end of visibility
        tstop  = (t2 - self.grb.t_trig).sec * u.s
 
        nslice = 0
        for i,t in enumerate(self.grb.time_interval):
            #ÿ print(i+1," ",end="")
            
            if (tstart >= t[0] and tstart <= t[1] and on == False):
#                print("Slice {:3d} at {:9.2f}, visible in   [{:9.2f}, {:9.2f}]"
#                      .format(i,tstart,t[0],t[1]))
                on = True
                self.tobs = [ [tstart, t[1]] ] # Create first tobs interval
                nslice += 1
                First = True
                # Check if, by chance tstop is in the same interval
                if (tstop <= t[1]):
#                    print(" Burst visibility stops in same interval")
                    self.tobs =  [ [tstart, tstop ] ]
                    on = False

            if (tstop >= t[0] and tstop <= t[1] and on == True):
#                print("Slice {:3d} at {:9.2f}, invisible in [{:9.2f}, {:9.2f}]"
#                      .format(i,tstop,t[0],t[1]))
                self.tobs.append([t[0],tstop]) # Create last tobs interval
                nslice +=1
                off = True

            if (on == True and off==False and First==False):
                self.tobs.append(t) # Append next intervals except if first
                nslice += 1

            First = False
        # end of loop over GRB time intervals

        self.tstart = tstart
        self.tstop  = tstop
        
        site     = mcf.pos_site[self.where]
        t = tstart + self.grb.t_trig
        self.pos_start = self.grb.radec.transform_to(AltAz(obstime = t,
                                                    location = site))
        t = tstop + self.grb.t_trig
        self.pos_stop = self.grb.radec.transform_to(AltAz(obstime = t,
                                                    location = site))
        self.ntbins = nslice
        if (self.dbg>0):
            print("                                    ===> {} slices found"
              .format(self.ntbins),end="\n")
        return nslice

    ##########################################################################
    def compare_time_interval(self):
        """
        Print the grb time intervals and the observation time intervals
        as a visual check that the observation intervalls were correctly
        set.
        2020/01/23 : OK
        """
        j = 0
        for i,t in enumerate(self.grb.time_interval):
            print("{:3d} : [{:9.2f}, {:9.2f}]"
                  .format(i,t[0],t[1]),end="")
            if (t[1] >= self.tobs[0][0] and j < len(self.tobs) ):
                print(" ***  [{:9.2f}, {:9.2f}] -> Alt =  [{:6.2f}, {:6.2f}]"
                      .format(self.tobs[j][0],self.tobs[j][1],
                              self.altaz[j][0].alt,self.altaz[j][1].alt) )
                j = j+1
            else:
                print()
        return

    ##########################################################################
    def create_perf(self,folder="."):
        """
        For each observation interval get the IRF file from the
        current position of the GRB in the sky and the cumulated observation
        time.
        - Observation point : end of observing interval (counts cumulated)
        - Position in the sky : beginning of intervall (higher flux)

        """
        print(" *** Create IRF set...",end="\n")

        obstime = 0
        site     = mcf.pos_site[self.where]

        for i, dt in enumerate(self.tobs):
            # print(i+1," ",end="")

            # Define observation window 
            obstime = dt[1] - dt[0] # Duration (not cumulated !)
            
            # - Date and alt-az boundaries
            treal = [dt[0] + self.grb.t_trig, dt[1] + self.grb.t_trig] 
            altaz = self.grb.radec.transform_to(AltAz(obstime  = treal,
                                                        location = site))
            self.altaz.append(altaz)

            # Get Irf from zenith, azimuth at beginning of slice, obstime
            if (self.ndof == 0):

                irf_file = IRF.get_irf_file(folder  = folder,
                                            theta   = 90*u.degree-altaz[0].alt,
                                            phi     = altaz[0].az,
                                            obstime = obstime,
                                            khem    = self.where)

                perfi        = CTAPerf_onaxis.read(irf_file)
                reco_energy = perfi.bkg.energy.edges
                # self.perf.peek()
            else:
                sys.exit("Likelihood not reimplemented")
                break
                perfi = load_cta_irfs(irf_file) # Load IRF
                reco_energy = perfi["bkg"].data.axis("energy").edges

            self.perf.append(perfi)
            self.Emin.append(min(reco_energy))
            self.Emax.append(max(reco_energy))
        
        if (self.dbg>0):
            print("                       ===> {} IRF found"
              .format(len(self.Emin)),end="\n")

        return
    ##########################################################################
    def source_spectrum(self,
                        index = 2,
                        norm  = 1e-11/(u.cm)**2/u.s/u.TeV,
                        Eref ="1 TeV"):
        """
        Likelihood analysis
        """

        spectral_fitmodel = PowerLaw(index     = index,
                                     amplitude = norm,
                                     reference = Eref)

        model = SkyModel(spatial_model  = self.grb.spatial,
                         spectral_model = spectral_fitmodel,
                         name           = "tobefitted")

        model.parameters["amplitude"].min = 0. # prevent negative signal

        # Fit spectral index
        if (self.ndof>1):
            model.parameters["index"].frozen     = False
        else:
            model.parameters["index"].frozen     = True

        # Fit source location
        if (self.ndof > 2):
            model.parameters["lon_0"].frozen     = False
            model.parameters["lat_0"].frozen     = False
        else:
            model.parameters["lon_0"].frozen     = True
            model.parameters["lat_0"].frozen     = True

        if (norm.value ==0):
            model.parameters["amplitude"].frozen = True
        else:
            model.parameters["amplitude"].frozen = False

#    #    ### Spectral model
#        if (self.dbg>1):
#            fig = plt.figure(figsize=(4,4))
#            a1 = fig.add_subplot(111)
#            a1.set_ylim(ymin=1e-15,ymax =1e-7)
#            model.plot(energy_range=[0.1*u.TeV,100*u.TeV],ax=a1)
#            plt.show()

        return model

    ##########################################################################
    def theory_prediction(self, irf, spectrum):
        """
        Likelihood analysis
        """

        expo, bckg, psf, edisp = irf
        ### Signal and background theory (to generate a MC trial)
        bckg_theory = BackgroundModel(bckg)
        grb_theory  = SkyModel(spatial_model  = self.grb.spatial,
                               spectral_model = spectrum,
                               name           = "Theory")
        datatheory  = MapDataset(model            = grb_theory,
                                 exposure         = expo,
                                 background_model = bckg_theory,
                                 psf              = psf,
                                 edisp            = edisp)

        if (self.dbg>1): DataSetShow(datatheory,label = "THEORY  ",plot = False)
        
        return datatheory.npred() # Countmap

###########################################################################
#
# Slices (observation) in the Monte Carlo GRB simulation
#
###########################################################################
    def loop_slices_onoff(self):
        """
        Get counts for each time slice.
        Stack observations.

        Loop over time interval and create a target, a set of observation parameters.
        Make the Monte Carlo simulations
        
        """
        self.simulations = []
        for idx, dt in enumerate(self.tobs):         # Loop over Time slices

            # if (debug_level>1): print("*===> Simulating interval ",idx)
            simu = mc_onoff(mc=self,
                            alpha = mcf.alpha, 
                            obstime =dt[1] - dt[0],
                            idx = idx)                
 
            self.simulations.append(simu)
            self.stack_obs = SpectrumDatasetOnOffStacker(Observations(simu))
            #stack.run()

        return

    def loop_slices_deprecated(self):
        """
        Get counts for each time slice.
        Stack observations.

        Loop over time interval and create a target, a set of observation parameters.
        Make the Monte Carlo simulations
        
        """
        self.simulations = []
        for idx, dt in enumerate(self.tobs):         # Loop over Time slices

            # if (debug_level>1): print("*===> Simulating interval ",idx)
            if (self.ndof == 0):
                simu = mc_onoff(mc=self,
                                alpha = mcf.alpha, 
                                obstime =dt[1] - dt[0],
                                idx = idx)                
            else:
                print("Not implemented")
                break
#            cumulated_npred  = None
#            obs_param = GRBObservationParameters(alpha    = self.alpha,
#                                                 livetime = livetime,
#                                                 emin     = self.Emin[idx],
#                                                 emax     = self.Emax[idx],
#                                                 pointing = self.grb.pointing,
#                                                 fov      = self.fov,
#                                                 binsize  = self.binsize)
#                livetime  = dt[1] - dt[0] # do not Cumulate livetime
#
                ### Get IRF for the current observation
#                irf = self.get_observation_irf(obs_param,
#                                               self.grb.pointing,
#                                               self.perf[idx],
#                                               debug_level)#
#                ### Get count map
#                npred = self.theory_prediction(irf,
#                                               self.grb.spectra[idx])
#                if (cumulated_npred == None):
#                    cumulated_npred = npred
#                else:
#                    cumulated_npred += npred
#                ### Generate MC trial
#                #rng = np.random.RandomState(seed=51) # same events at each trial
#                rng    = np.random.RandomState()
#                counts = rng.poisson(cumulated_npred.data) # Compute counts
#                counts_map = WcsNDMap(obs_param.geom, counts)
#
#                simu  =  fit_3d(self, counts_map, irf, debug = debug_level)

            
            # End of loop overs time slices
            self.simulations.append(simu)

            if (self.ndof == 0): 
                self.stack_obs = SpectrumDatasetOnOffStacker(Observations(simu))
                #stack.run()


        # This goes back to the Monte Carlo module to treat the
        # information collected in this simulation
        return
    ###########################################################################
    def dump_slices(self,iMC=0,non=0,noff=0,sigma=0,file=None,header=False):
        """
        Print ou the list of on and off counts and the corresponding
        Li & Ma value. If non and off are less than a minimul value
        return a True status flag for further actions
        """
        status = False
        
        if (header):
            print("     ",file=file)
            for i, dt in enumerate(self.tobs):
                print("{:8.1f}".format(dt[1]-dt[0]), end="    ",file = file)
            print(file=file)

        else:
            print("{:4d}".format(iMC),end="",file = file)
            for i in range(len(non)): 
                if(non[i]<=mcf.n_lima_min): 
                    sep="*   "
                    status = True
                    self.err_slice += [i]
                else:
                    sep="    "
                print("{:8.1f}".format(non[i]), end=sep,file = file)
            print(file = file)
    
            print("{:4d}".format(iMC),end="",file = file)
            for i in range(len(noff)): 
                if(noff[i]<=mcf.n_lima_min): 
                    sep="*   "
                    status = True
                    self.err_slice += [i]
                else:
                    sep="    "
                print("{:8.1f}".format(noff[i]), end=sep,file = file)
            print(file = file)
            
            print("{:4d}".format(iMC),end="",file = file)
            for i in range(len(sigma)):
                if (sigma[i]>=3):
                    sep = "*** "
                else:
                    sep= "    "
                print("{:8.1f}".format(sigma[i]), end=sep,file = file)
            print(file = file)

        return status
    