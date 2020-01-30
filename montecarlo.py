# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import time

from   pathlib import Path

import astropy.units as u
from   astropy.visualization import quantity_support
from   astropy.coordinates   import AltAz, EarthLocation

from scipy.stats import chi2
#from scipy.stats import norm

from fit_3d import fit_3d, DataSetShow
from grb import GRBObservationParameters
from irf import IRF
from cta_irf_onaxis import CTAPerf_onaxis
from gammapy.spectrum.models import PowerLaw
from gammapy.maps import WcsNDMap

from gammapy.cube import MapDataset
from gammapy.cube.models import SkyModel, BackgroundModel

from fit_onoff import mc_onoff
sys.path.append("../../../utilities_ths/")  
sys.path.append("analysis/") 
from utilities_ths import MyLabel
from utilities_ths import stamp

__all__ = ['MonteCarlo']
###############################################################################
class MonteCarlo():
    """
    2020/01/23 : OK
    - Maybe find the Site location from 'where' here to avoid 
    requesting it at several places in the code ?
    """

    def __init__(self,grb,niter=1,
                 alpha  = 0.2, 
                 fov    = [5.*u.deg, 0.125*u.deg],
                 dtslew = [30*u.s, True],
                 where  = "North",
                 fixed  = False,
                 ndof   = 0,
                 clevel = 0.9):
        """
        Initialize class members to default values

        """
        zdeg = 0*u.deg
        # Input parameters and objects
        self.niter        = niter     # Number of trials
        self.fov          = fov[0]    # Field-of-view width 
        self.binsize      = fov[1]    # Field-of-view bin size
        self.dtslew       = dtslew[0] # Slewing time
        self.fix_slew     = dtslew[1] # Fixed (True, or random from zero)
        self.fixed        = fixed  # Source alt-az changes if False
        self.alpha        = alpha  # Th on/off ratio in the on-off analysis
        self.ndof         = ndof   # Deg. of freedom in likelikood, 0 if on-off
        self.clevel       = clevel # Confidence level
        self.grb          = grb    # A GRB instantiation
        self.where        = where  # Site (North or South)

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

        # list of simulations
        self.simulations = None

        # Significance over simulations
        self.id_smax_list  = [] # Slice number to get back the time/ altaz
        self.smax_list     = [] # List of max. significances along trials
        self.nex_smax_list = [] # List of excess counts at max. signif. 
        self.nb_smax_list  = [] # List of background counts at max. signif.

        self.id_3s_list    = [] # Slice number ot get back the time/altaz
        self.nex_3s_list   = [] # List of excess
        self.nb_3s_list    = [] # List of background

        self.id_5s_list    = [] # Slice number to get back the time/altaz
        self.nex_5s_list   = [] # List of excess
        self.nb_5s_list    = [] # List of background

        # Related to time and time bins
        self.detect_3s = 0  # Fraction of trials 3 sigma was reached
        self.detect_5s = 0  # Fraction of trials 5 sigma was reached

        # Mean sigma versus time - one value per time slice
        self.sigma_mean = []
        self.sigma_rms  = []

        self.mctime = 0
        self.abort  = -1 # Event simulation status :
                         #  abort = -1 : never started (GRB not visible)
                         #  abort =  0 : completed (all trials simulated) 
                         #  abort >  0 : aborted at 'abort' slice number 

        return

    ##########################################################################
    def observation_slices(self):
        """
        Create the time intervals for observation as the result of the
        GRB time interval and the visibility slot.
        Get the boundaries of the observation window (tstart, tstop).
        Get the sky position at the boundaries.
        
        2020/01/23 : OK
            - check the t00 trick, don't remeber why it is here
        """

        on  = False
        off = False
        tvis   = self.grb.t_true[self.where] # Visibilty window for the site

        # Define visibility start and stop in GRB trigger reference system
        # tstart can only be positive since t_true start at latest when
        # the GRB is seen
        tstart = (tvis[0] - self.grb.t_trig).sec * u.s
        t00 = self.grb.time_interval[0][0]
 
        # Patch - check that
        if tstart < t00:
            tstart = t00
            print("tstart =",self.tstart," reshifted to t00 = ",t00)
 
        if (self.fix_slew):
            tstart = tstart + self.dtslew
        else:
            tstart = tstart + self.dtslew*np.random.random()
        tstop  = (tvis[1] - self.grb.t_trig).sec * u.s
 
        ibins = 0
        for i,t in enumerate(self.grb.time_interval):

            if (tstart >= t[0] and tstart <= t[1] and on == False):
#                print("Slice {:3d} at {:9.2f}, visible in   [{:9.2f}, {:9.2f}]"
#                      .format(i,tstart,t[0],t[1]))
                on = True
                self.tobs = [ [tstart, t[1]] ] # Create first tobs interval
                ibins += 1
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
                off = True

            if (on == True and off==False and First==False):
                self.tobs.append(t) # Append next intervals except if first
                ibins += 1

            First = False
        # end of loop over GRB time intervals

        self.tstart = tstart
        self.tstop  = tstop
        
        site     = EarthLocation.of_site(self.grb.site[self.where])
        t = tstart + self.grb.t_trig
        self.pos_start = self.grb.radec.transform_to(AltAz(obstime = t,
                                                    location = site))
        t = tstop + self.grb.t_trig
        self.pos_stop = self.grb.radec.transform_to(AltAz(obstime = t,
                                                    location = site))
        self.ntbins = ibins

        return

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
    def create_perf(self):
        """
        For each observation interval get the IRF file from the
        current position of the GRB in the sky and the cumulated observation
        time.
        - Observation point : end of observing interval (counts cumulated)
        - Position in the sky : beginning of intervall (higher flux)

        """

        obstime = 0
        site     = EarthLocation.of_site(self.grb.site[self.where])

        for i, dt in enumerate(self.tobs):

            # Define observation window 
            obstime = dt[1] - dt[0] # Duration (not cumulated !)
            
            # - Date and alt-az boundaries
            dtreal = [dt[0] + self.grb.t_trig, dt[1] + self.grb.t_trig] 
            altaz = self.grb.radec.transform_to(AltAz(obstime  = dtreal,
                                                        location = site))
            self.altaz.append(altaz)

#            print(" Where = ",self.where)
#            print(" treal = ",dtreal[0].datetime,dtreal[1].datetime)
#            print(" alt   = ",altaz.alt,"-> zen",90*u.deg-altaz.alt )
#            print(" az    = ",altaz.az)
#            print(" observed since ",obstime)

            # Get Irf from zenith, azimuth at beginning of slice, obstime
            if (self.ndof == 0):

                irf_file = IRF.get_irf_file(theta   = 90*u.degree-altaz[0].alt,
                                            phi     = altaz[0].az,
                                            obstime = obstime,
                                            khem    = self.where)

                perfi        = CTAPerf_onaxis.read(irf_file)
                reco_energy = perfi.bkg.energy.edges
                # self.perf.peek()
            else:
                sys.exit("Likelihood not reimplemented")
                break
                perfi = IRF.load_cta_irfs(irf_file) # Load IRF
                reco_energy = perfi["bkg"].data.axis("energy").edges

            self.perf.append(perfi)
            self.Emin.append(min(reco_energy))
            self.Emax.append(max(reco_energy))

        return
    ##########################################################################
    def source_spectrum(self,
                        index = 2,
                        norm  = 1e-11/(u.cm)**2/u.s/u.TeV,
                        reference ="1 TeV",
                        debug=True):
        """
        Likelihood analysis
        """

        spectral_fitmodel = PowerLaw(index     = index,
                                     amplitude = norm,
                                     reference = reference)

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

    #    ### Spectral model

#        if (debug):
#            print("==========================================================")
#            fig = plt.figure(figsize=(4,4))
#            a1 = fig.add_subplot(111)
#            a1.set_ylim(ymin=1e-15,ymax =1e-7)
#            model.plot(energy_range=[0.1*u.TeV,100*u.TeV],ax=a1)
#            plt.show()

        return model

    ##########################################################################
    def theory_prediction(mc, irf, spectrum, debug=True):
        """
        Likelihood analysis
        """

        expo, bckg, psf, edisp = irf
        ### Signal and background theory (to generate a MC trial)
        bckg_theory = BackgroundModel(bckg)
        grb_theory  = SkyModel(spatial_model  = mc.grb.spatial,
                               spectral_model = spectrum,
                               name           = "Theory")
        datatheory  = MapDataset(model            = grb_theory,
                                 exposure         = expo,
                                 background_model = bckg_theory,
                                 psf              = psf,
                                 edisp            = edisp)

        if (debug): DataSetShow(datatheory,label = "THEORY  ",plot = False)
        return datatheory.npred() # Countmap


    ##########################################################################
    def run(self, dbg_level):
        """
        - Run simulations of the current grb - Likelihood method

        - Compute :
            - significance on the full observation,
            - significance max on intervals,
            - time of siginificance max,
            - time of 3 or 5 sigma reached
        """

        # RUN only if vsisible
        if (self.grb.vis_tonight[self.where]):
            self.observation_slices() # Create observation time slices
            self.create_perf()        # Get performance for each time slice
        else:
            return # mc variables remain at their default (see Constructor)

        n_3s = 0 # Number of time 3 sigma was reached
        n_5s = 0 # Number of time 5 sigma was reached
        sigma_sum = 0
        sigma2_sum = 0

        self.mctime = time.time() # Starts chronometer
        abort_test  = False        # Abortion test not performed
        self.abort  = 0            # GRB is visible
        
        ###############################################
        ### Monte Carlo iterations
        ###############################################
        iMC=1
        
        while(iMC <= self.niter):
            if (self.niter <= 10) or (np.mod(iMC,10) == 0):
                print("#",iMC," ",end="")

            ### Loop on GRB time slices - get the signal
            self.loop_slices(dbg_level)

            ### Analyze the current simulation
            if (self.ndof != 0):
                ################################ 3D
                fcnsb = np.array([x[0] for x in self.simulations]) # lik_fit
                fcnb  = np.array([x[1] for x in self.simulations]) # lik_bonly
                nex   = np.array([x[2] for x in self.simulations])
                nb    = np.array([x[3] for x in self.simulations])
                # ns_bonly - not used
                # nb_bonly - not used

                TS =  fcnb - fcnsb
                TS[TS < 1e-3] = 0 # Avoid sligtly negative values

                # pvalue  = 1 - chi2.cdf(TS, ndof)
                # sigma   = norm.ppf(1-pvalue)  # Integral from -Inf, can be<0
                # print(" Warning sqrt TS")
                sigma = np.sqrt(TS)
                sigma[sigma < 0 ] = 0
                sigma[sigma > 99] = 99

                if (dbg_level>1): self.plot_trial(sigma,nex,nb)

            if (self.ndof == 0):
                ################################ On-Off
                ### Get cumulated statistics
                obs_stat = self.grb.get_cumulative_stats(self.simulations)
                sigma    = obs_stat['sigma']
                non      = obs_stat['n_on']
                noff     = obs_stat['n_off']

                nex = non -self.alpha*noff
                nb  = self.alpha*noff

            # End of test on likelihood case

            # Acummulate sum and sum**2 for mean / error in each sl
            sigma_sum  += sigma
            sigma2_sum += sigma**2

            # Update stat list
            i3s, i5s = self.fill_stat(sigma, nex, nb, dbg_level)
            n_3s += i3s
            n_5s += i5s

            # Test whether it is worth to continue
            # In case 3 signma is not reached in the first 10% of the trials
            # then the 90% CL can not be reached.
            if (abort_test==False):
#                print("Trial {:4d} /{:5d} : {:4.2f} (s = {:5.2f}, n3s = {:5d})"
#                      .format(iMC,self.niter, iMC/self.niter,max(sigma),n_3s))

                if (iMC/self.niter > 1 - self.clevel):
                    abort_test = True
                    if (n_3s ==0):
                        print()
                        print("No chance to reach ",100*self.clevel,"% of CL")
                        print(" Simulation stopped")
                        self.abort = iMC
                        break

            iMC+=1 # End of MC loop

        self.mctime = (time.time() - self.mctime)/self.niter


        ### Compute mean sigma and rms for all trials
        self.sigma_mean = sigma_sum/self.niter
        sigma2_mean     = sigma2_sum/self.niter
        self.sigma_rms  = np.sqrt(sigma2_mean-self.sigma_mean**2)

        ### Compute detection level
        self.detect_3s = n_3s/self.niter
        self.detect_5s = n_5s/self.niter

        return


    ##########################################################################
    def fill_stat(self,sigma, nex, nb, dbg_level):

        ### Find maximum
        maxidx = np.nanargmax(sigma) # If sigma is nan, it is the max !
        self.id_smax_list.append(maxidx)
        self.smax_list.append(sigma[maxidx])
        self.nex_smax_list.append(nex[maxidx])
        self.nb_smax_list.append(nb[maxidx])

        # Find when 3 sigma is reached
        nex_3s = -1
        nb_3s  = -1
        t_3s   = -1
        alt_3s = -1
        i3s    = 0

        id_3s = np.where(sigma>=3)[0]
        if (np.size(id_3s) != 0):
            id_3s  = id_3s[0]
            nex_3s = nex[id_3s]
            nb_3s  = nb[id_3s]
            t_3s   = self.tobs[id_3s][1] # end of slice
            alt_3s = self.altaz[id_3s][1].alt # end of slice
            self.id_3s_list.append(id_3s)
            i3s = 1

        self.nex_3s_list.append(nex_3s)
        self.nb_3s_list.append(nb_3s)

        # Find when 5 sigma is reached
        nex_5s = -1
        nb_5s  = -1
        t_5s   = -1
        alt_5s = -1
        i5s    = 0

        id_5s = np.where(sigma>=5)[0]
        if (np.size(id_5s) != 0):
            id_5s  = id_5s[0]
            nex_5s = nex[id_5s]
            nb_5s  = nb[id_5s]
            t_5s   = self.tobs[id_5s][1] # end of slice
            alt_5s = self.altaz[id_5s][1].alt # end of slice
            self.id_5s_list.append(id_5s)
            i5s = 1

        self.nex_5s_list.append(nex_5s)
        self.nb_5s_list.append(nb_5s)

        if (dbg_level):
            print(" >>> t_smax = {:5.2f}, ({:5.2f} sigma at slice {:2d})"
                  .format(self.tobs[maxidx][0],sigma[maxidx],maxidx+1))
            print("     t_3s / t_5s  at {:5.2f}/{:5.2f} <- alt = {:6.2f}/{:6.2f}"
                  .format(t_3s,t_5s,alt_3s,alt_5s))

        return i3s, i5s

    ###########################################################################
    def loop_slices(self,debug_level):
        """
        Compute counts for each time slice

        Loop over time interval and create a target, a set of observation parameters.
        Make the Monte Carlo simulations
        
        
        
        """
        self.simulations = []

        # Loop over Time slices
        for idx, dt in enumerate(self.tobs):

            # if (debug_level>1): print("*===> Simulating interval ",idx)

            if (self.ndof != 0):
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
#                obs_param = GRBObservationParameters(alpha    = self.alpha,
#                                                 livetime = livetime,
#                                                 emin     = self.Emin[idx],
#                                                 emax     = self.Emax[idx],
#                                                 pointing = self.grb.pointing,
#                                                 fov      = self.fov,
#                                                 binsize  = self.binsize)
#
                ### Get IRF for the current observation
#                irf = self.get_observation_irf(obs_param,
#                                               self.grb.pointing,
#                                               self.perf[idx],
#                                               debug_level)
#
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

            else:
                livetime = dt[1] - dt[0] # Current slice livetime - not cumulated !

                param = GRBObservationParameters(alpha    = self.alpha,
                                                 livetime = livetime,
                                                 emin     = self.Emin[idx],
                                                 emax     = self.Emax[idx],
                                                 pointing = self.grb.pointing,
                                                 fov      = self.fov,
                                                 binsize  = self.binsize)
                
                simu = mc_onoff(mc = self, 
                                idx = idx, 
                                obs_param = param,
                                debug = debug_level)

            # End of loop overs time slices
            self.simulations.append(simu)

            if (self.ndof == 0): 
                self.stack_obs = self.grb.add_stack_OnOffobs(simu)

        # This goes back to the Monte Carlo module to treat the
        # information collected in this simulation
        return


#    @staticmethod
#    def plot(simu, target):
#        """Plot some simulation results"""
#
#        import matplotlib.pyplot as plt
#        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
#                                       figsize=(10, 5))
#
#        # Spectrum plot
#        energy_range = [0.01 * u.TeV, 100 * u.TeV]
#        target.model.plot(ax=ax1, energy_range=energy_range,
#                          label='Model')
#        plt.text(0.55, 0.65, target.__str__(),
#                 style='italic', transform=ax1.transAxes, fontsize=7,
#                 bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
#        ax1.set_xlim([energy_range[0].value, energy_range[1].value])
#        ax1.set_ylim(1.e-17, 1.e-5)
#        ax1.grid(which='both')
#        ax1.legend(loc=0)
#
#        # Counts plot
#        on_off = simu.on_vector.data.data.value
#        off = 1. / simu.off_vector.backscal * simu.off_vector.data.data.value
#        excess = on_off - off
#        bins = simu.on_vector.energy.lo.value
#        x = simu.on_vector.energy.nodes.value
#        ax2.hist(x, bins=bins, weights=on_off,
#                 facecolor='blue', alpha=1, label='ON')
#        ax2.hist(x, bins=bins, weights=off,
#                 facecolor='green', alpha=1, label='OFF')
#        ax2.hist(x, bins=bins, weights=excess,
#                 facecolor='red', alpha=1, label='EXCESS')
#        ax2.legend(loc='best')
#        ax2.set_xscale('log')
#        ax2.set_xlabel('Energy [TeV]')
#        ax2.set_ylabel('Expected counts')
#        ax2.set_xlim([energy_range[0].value, energy_range[1].value])
#        ax2.set_ylim([0.0001, on_off.max() * (1 + 0.05)])
#        ax2.vlines(simu.lo_threshold.value, 0, 1.1 * on_off.max(),
#                   linestyles='dashed')
#        ax2.grid(which='both')
#        ax2.text(0.55, 0.05, simu.__str__(),
#                 style='italic', transform=ax2.transAxes, fontsize=7,
#                 bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
#        plt.tight_layout()
#
#

    ###########################################################################
    def result(self,debug = 0, out=sys.stdout,popfile=None,write_header="False"):
        """Print out results"""


        print("\n========================== RESULTS",
              "site =",self.where," ===",file=out,)
        print("  GRB name             : ",self.grb.name,file=out)

        print("  Redshift           z : {:6.2f} "
              .format(self.grb.z),file=out)

        print("  Fixed alt-az         : ",self.fixed,file=out)

        print("  Slewing time (fixed) : ",self.dtslew,file=out)
        
        print("  Observation window   : ")
        print("     - Time       :  {:6.2f} - {:6.2f}"
              .format(self.tstart,self.tstop,file=out))
        print("     - Alt. (zen) :  {:6.2f} - {:6.2f} ({:6.2f} - {:6.2f})"
              .format(self.pos_start.alt,
                      self.pos_stop.alt,
                      90*u.deg-self.pos_start.alt,
                      90*u.deg-self.pos_stop.alt,
                      file=out))
        print("     - Azimuth    :  {:6.2f} to {:6.2f}"
              .format(self.pos_start.az,
                      self.pos_stop.az,
                      file=out))
        print("----------------------------------------------------",file=out)
        print("STATUS ",end="",file=out)
        if (self.abort != 0):
            print("===> Simulation was aborted at trial ",self.abort,file=out)
        else:
            if (self.detect_3s >= self.clevel):
                print(" ===> 3 sigma detected",file=out)
                if (self.detect_5s >= self.clevel):
                    print(" ===> 5 sigma detected",file=out)
                else:
                    print(" ===> NOT detected at 5 sigma",file=out)
            else:
                print(" ===> NOT detected at 3 sigma",file=out)

        # If simulation went well

        t_unit = u.s
        alt_unit = u.deg

        smax_mean    = -1
        err_smax     = -1
        t_smax       = -1 *t_unit
        err_t_smax   = -1 *t_unit
        alt_smax     = -1 *alt_unit
        err_alt_smax = -1 *alt_unit
        az_smax      = -1 *alt_unit
        err_az_smax  = -1 *alt_unit
        nex_smax     = -1
        nb_smax      = -1

        nex_3s      = -1
        nb_3s       = -1
        t_3s        = -1 *t_unit
        err_t_3s    = -1 *t_unit
        alt_3s      = -1 *alt_unit
        err_alt_3s  = -1 *alt_unit
        az_3s       = -1 *alt_unit
        err_az_3s   = -1 *alt_unit
        
        nex_5s      = -1
        nb_5s       = -1
        t_5s        = -1 *t_unit
        err_t_5s    = -1 *t_unit
        alt_5s      = -1 *alt_unit
        err_alt_5s  = -1 *alt_unit
        az_5s       = -1 *alt_unit
        err_az_5s   = -1 *alt_unit
        
        if (self.abort ==0):

            t_unit   = self.tobs[0][0].unit  # Get time unit, useful in what follows
            alt_unit = self.altaz[0][0].alt.unit  # Get time unit, useful in what follows

            # Max significance is reached
            smax_mean  = np.mean(self.smax_list)
            err_smax   = np.std(self.smax_list)
            nex_smax   = np.mean(self.nex_smax_list)
            nb_smax    = np.mean(self.nb_smax_list)

            t = [self.tobs[i][0].value
                 for i in self.id_smax_list
                 if self.tobs[i][0].value>=0]
            t_smax     = np.mean(t) * t_unit
            err_t_smax = np.std(t) * t_unit
            
            alt = [self.altaz[i][0].alt.value
                                 for i in self.id_smax_list
                                 if self.altaz[i][0].alt.value>=0]
            alt_smax   = np.mean(alt)* alt_unit
            err_alt_smax = np.std(alt)* alt_unit
            
            az = [self.altaz[i][0].az.value
                                 for i in self.id_smax_list
                                 if self.altaz[i][0].az.value>=0]
            az_smax   = np.mean(az)* alt_unit
            err_az_smax = np.std(az)* alt_unit
            
            ### 3 sigma
            if (len(self.id_3s_list)):
                nex = [n for n in self.nex_3s_list if n >=0]
                nb  = [n for n in self.nb_3s_list  if n >=0]
                nex_3s    = np.mean(nex)
                nb_3s     = np.mean(nb)
                
                t = [self.tobs[i][0].value
                     for i in self.id_3s_list
                     if self.tobs[i][0].value>=0]
                t_3s      = np.mean(t) *t_unit
                err_t_3s  = np.std(t) *t_unit
                
                alt =  [self.altaz[i][0].alt.value
                                 for i in self.id_3s_list
                                 if self.altaz[i][0].alt.value>=0]
                alt_3s    = np.mean(alt) * alt_unit
                err_alt_3s = np.std(alt)* alt_unit
                
                az =  [self.altaz[i][0].az.value
                                 for i in self.id_3s_list
                                 if self.altaz[i][0].az.value>=0]
                az_3s    = np.mean(az) * alt_unit
                err_az_3s = np.std(az)* alt_unit

            ### 5 sigma
            if (len(self.id_5s_list)):
                nex = [n for n in self.nex_5s_list if n >=0]
                nb  = [n for n in self.nb_5s_list  if n >=0]
                nex_5s    = np.mean(nex)
                nb_5s     = np.mean(nb)
                
                t = [self.tobs[i][0].value
                     for i in self.id_5s_list
                     if self.tobs[i][0].value>=0]
                t_5s      = np.mean(t) *t_unit
                err_t_5s  = np.std(t) *t_unit
                
                alt =  [self.altaz[i][0].alt.value
                                 for i in self.id_5s_list
                                 if self.altaz[i][0].alt.value>=0]
                alt_5s    = np.mean(alt) * alt_unit
                err_alt_5s = np.std(alt)* alt_unit
                
                az =  [self.altaz[i][0].az.value
                                 for i in self.id_5s_list
                                 if self.altaz[i][0].az.value>=0]
                az_5s    = np.mean(az) * alt_unit
                err_az_5s = np.std(az)* alt_unit

            print("----------------------------------------------------",file=out)

            print("  E range              : {:6.2f} {:6.2f}"
                  .format(self.Emin[0], self.Emax[0]),file=out)

            print("  Observation slices   : {:3d}"
                  .format(len(self.tobs)),end="",file=out)

            if ( len(self.tobs) > 1) :
                print(" from  {:6.2f} to {:6.2f}"
                      .format(self.tobs[0][0],self.tobs[-1][1]),file=out)
            else:
                print(" from  {:6.2f} to {:6.2f}"
                      .format(self.tobs[0][0],self.tobs[0][1]),file=out)

            print("====================================================",file=out)

            print("  Max. significance       : {:5.1f}   +/- {:5.1f}"
                  .format(smax_mean,err_smax,file=out))

            print("                     time : {:5.1f} +/- {:5.1f}"
                   .format(t_smax,err_t_smax,file=out))
            print("                     alt. : {:5.1f} +/- {:5.1f}"
                   .format(alt_smax,err_alt_smax,file=out))
            print("                      az. : {:5.1f} +/- {:5.1f}"
                   .format(az_smax,err_az_smax,file=out))
            print("              S, B counts : {:8.1f} {:0.1f}"
                  .format(nex_smax,nb_smax,file=out))
            print()

            print("  3 sigma reached at time : {:5.1f} +/- {:5.1f}"
                  .format(t_3s,err_t_3s,file=out))

            print("                     alt. : {:5.1f} +/- {:5.1f}"
                   .format(alt_3s,err_alt_3s,file=out))
            print("                      az. : {:5.1f} +/- {:5.1f}"
                   .format(az_3s,err_az_3s,file=out))
            
            print("              S, B counts : {:8.1f}, {:8.1f}"
                  .format(nex_3s,nb_3s,file=out))

            print("                  Fraction:",100*self.detect_3s,"%",file=out)
            print()


            print("  5 sigma reached at time : {:5.1f} +/- {:5.1f}"
                  .format(t_5s,err_t_5s,file=out))
            print("                     alt. : {:5.1f} +/- {:5.1f}"
                   .format(alt_5s,err_alt_5s,file=out))
            print("                      az. : {:5.1f} +/- {:5.1f}"
                   .format(az_5s,err_az_5s,file=out))


            print("              S, B counts : {:8.1f} {:8.1f}"
                  .format(nex_5s,nb_5s,file=out))

            print("                  Fraction:",100*self.detect_5s,"%",file=out)

            print("====================================================\n",file=out)


        # Dump result into population file - even if simulation was not
        # performed or aborted
        if (write_header==True):
            txt = "{:>10} {:>9} {:>5} " # GRB, Eiso, z
            txt+= "{:>6} {:>6} {:>6} "  # site, Ra, Dec
            txt+= "{:>15} " # Trigger time
            txt+= "{:>10} {:>10} " # tstart, tstop,  (ref:trigger date)
            txt+= "{:>10} {:>10} " # alt-start, alt-stop
            txt+= "{:>10} {:>10} " # az-start, az-stop, 
            txt+= "{:>3} {:>4} " # nbins, ndof
            txt+= "{:>8} {:>10} {:>10} " # sigmax nex_max nb_max 
            txt+= "{:>10} {:>10} {:>10} " # tmax  alt_max az_max 
            txt+= "{:>10} {:>10} {:>10} {:>6} {:>6} " # nex_3s, nb_3s, t_3s, alt_3s, az_3s
            txt+= "{:>5} " # detect_3s
            txt+= "{:>10} {:>10} {:>10} {:>6} {:>6} " # nex_5s, nb_5s, t_5s, alt_5s, az_5s
            txt+= "{:>5} " #  detect_5s
            txt+= "{:>6} {:>5s}" # computing time for this GRB
            print(txt.format("name", "eiso",  "z",
                             "site","ra","dec",
                             "ttrig",
                             "t1","t2",
                             "alt1","alt2",
                             "az1","az2",
                             "nt", "ndof",
                             "sigmax", "nex_max","nb_max",
                             "t_max","altmax","azmax",
                             "nex_3s", "nb_3s",  "t3s", "alt3s",  "az3s",
                             "d3s",
                             "nex_5s", "nb_5s",  "t5s", "alt5s",  "az5s",
                             "d5s",
                             "mct","abort"),
                              file = popfile)

        txt = "{:>10} {:9.2e} {:5.2f} " # GRB, Eiso, z
        txt+= "{:>6} {:>6.2f} {:>6.2f} " # IRF : site, Ra, Dec
        txt+= "{:>15.8f} " # Trigger time
        txt+= "{:>10.2f} {:>10.2f} " # tstart, tstop, (ref:trigger date)
        txt+= "{:>10.2f} {:>10.2f} " # alt-start, alt-stop
        txt+= "{:>10.2f} {:>10.2f} " # az-start, az-stop, 
        txt+= "{:>3d} {:>4d} " # nbins, ndof
        txt+= "{:>8.2f} {:>10.2f} {:>10.2f} " # sigmax nex_max nb_max 
        txt+= "{:>10.2f} {:>10.2f} {:>10.2f} " # tmax  alt_max az_max 
        txt+= "{:>10.2f} {:>10.2f} {:>10.2f} {:>6.2f} {:>6.2f} "
        txt+= "{:>5.2f} "
        txt+= "{:>10.2f} {:>10.2f} {:>10.2f} {:>6.2f} {:>6.2f} "
        txt+= "{:>5.2f} "
        txt+= "{:>6.2f} {:>5d}"

        print(txt.format(
              self.grb.name, self.grb.Eiso.value, self.grb.z,
              self.where, self.grb.radec.ra.value, self.grb.radec.dec.value,
              self.grb.t_trig.mjd,
              self.tstart.value, self.tstop.value,
              self.pos_start.alt.value, self.pos_stop.alt.value,
              self.pos_start.az.value, self.pos_stop.az.value,
              self.ntbins,self.ndof,
              smax_mean, nex_smax, nb_smax, 
              t_smax.value, alt_smax.value, az_smax.value,
              nex_3s, nb_3s, t_3s.value, alt_3s.value, az_3s.value, 
              self.detect_3s,
              nex_5s, nb_5s, t_5s.value, alt_5s.value, az_5s.value,
              self.detect_5s,
              self.mctime, self.abort),
              file = popfile)

        if (debug>1 and self.abort != -1) :
            self.compare_time_interval()

        return

 ###############################################################################
    def plot2(self, saveplots, outfile,ref="VIS"):

        # Plot sigma versus time on the full
        # GRB lifetime together with visibility

        coord = EarthLocation.of_site(self.grb.site[self.where])
        dtplot = 500*u.s # Add some space for the plot

        # Visibility altitude curve
        tvis1 = self.grb.t_true[self.where][0] # End of visibility
        tvis2 = self.grb.t_true[self.where][1] # End of visibility
        dtvis = (tvis2-tvis1).to(u.s)
        tvis = np.linspace(0,dtvis.value,100)*dtvis.unit

        # Set the reference time
        event_name = self.grb.name+" - "+self.where
        if (ref=="VIS"):
            tref = tvis1
            text_ref =event_name + " (Ref: vis. start)"
        else:
            tref = self.grb.t_trig
            text_ref = event_name+ " (Ref: GRB trigger)"

        t_trig_rel = self.grb.t_trig - tref


        # Compute alt-az for visibility time sampling
        tvis_abs = tvis1 + tvis
        tvis_rel = tvis_abs - tref
        altazvis = self.grb.radec.transform_to( AltAz(obstime= tvis_abs,
                                                      location=coord))

        # Compute GRB measurement points from the given intervals
        tgrb_abs   = self.grb.t_trig + self.grb.t_s # Absolute time
        tgrb_rel   = tgrb_abs - tref
        altazgrb   = self.grb.radec.transform_to( AltAz(obstime= tgrb_abs,
                                                        location=coord))

        # Compute GRB observation points from the given intervals
        tobs = [t[1].value for t in self.tobs ]*self.tobs[0][0].unit
        tgrb_obs_abs   = self.grb.t_trig + tobs # Absolute time
        tgrb_obs_rel   = tgrb_obs_abs - tref
        altazobs   = self.grb.radec.transform_to( AltAz(obstime= tgrb_obs_abs,
                                                        location=coord))

        with quantity_support():
            fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(20,8))

            # Relative to reference
            ax1.plot(tvis_rel.sec, altazvis.alt,
                     linestyle="--",
                     color="b",
                     label="Altitude")

            ax1.plot(tgrb_rel.sec,altazgrb.alt,
                     linestyle="",
                     marker="*",
                     markerfacecolor="b",
                     label="GRB end of intervals")

            ax1.plot(tgrb_obs_rel.sec,altazobs.alt,
                     linestyle="",
                     marker="o",
                     markersize=10,
                     markerfacecolor="none",
                     markeredgewidth = 1.5,
                     markeredgecolor ="r",
                     label="GRB measurements")

            ax1.axvline(x=t_trig_rel.sec,
                        ls=":",
                        color="grey",
                        label="Trigger")

            ax1.axvline(x=(tvis1 - tref).sec,
                        ls=":",
                        color="green",
                        label="Start")
            ax1.axvline(x=(tvis2 - tref).sec,
                        ls=":",
                        color="red",
                        label="Stop")
            ax1.axhline(y=10*u.deg,
                        ls="dashdot",
                        color="grey")

            ax1.set_title(text_ref,fontsize=12)
            ax1.set_xlim(xmin=-dtplot,xmax= (tvis2 -tref + dtplot ).sec)
            ax1.set_ylim(ymin=0*u.deg)

            if (ref=="VIS"):
                ax1.set_xlabel("Elapsed time since visible (s)")
            else:
                ax1.set_xlabel("Elapsed time since Trigger (s)")
            ax1.set_ylabel("Altitude")

            #ax1.set_xscale("log")
            # ax1.legend(loc=8) # This could overlab the points/curves

            # Absolute time
            ax2.plot(tvis_abs.datetime,altazvis.alt,
                     linestyle="--",
                     color="b",
                     label="Altitude")
            ax2.plot(tgrb_abs.datetime,altazgrb.alt,
                     linestyle="",
                     marker="*",
                     markerfacecolor="b",
                     label="GRB end of intervals")
            ax2.plot(tgrb_obs_abs.datetime,altazobs.alt,
                     linestyle="",
                     marker="o",
                     markersize=10,
                     markerfacecolor="none",
                     markeredgewidth = 1.5,
                     markeredgecolor ="r",
                     label="GRB measurements")

            ax2.axvline(x=self.grb.t_trig.datetime ,
                        ls=":",
                        color="grey",
                        label = "Trigger")

            ax2.axvline(x=tvis1.datetime,
                        ls=":",
                        color="green",
                        label="Start")
            ax2.axvline(x=tvis2.datetime,
                        ls=":",
                        color="red",
                        label="Stop")

            ax2.axhline(y=10*u.deg,
                        ls="dashdot",
                        color="grey")

            ax2.tick_params(axis='x', rotation=70)
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            ax2.set_title(text_ref,fontsize=12)
            ax2.set_xlim(xmin=(tref-dtplot).datetime,xmax=(tvis2 + dtplot).datetime)
            ax2.set_ylim(ymin=0*u.deg)
            ax2.legend(bbox_to_anchor=(1.05, 0.5))

            if (ref=="VIS"):
                ax2.set_xlabel("Date since visible (hh:mm)")
            else:
                ax2.set_xlabel("Date since Trigger (hh:mm)")
            ax2.set_ylabel("Altitude")


            plt.show()

#            plt.errorbar(dtgrb+self.grb.t_s,self.sigma_mean,yerr=self.sigma_rms,fmt='o')
#            plt.show()
#        dt = np.linspace(0,duration.value,100)
#        t =  tmin + dt
#
#        with quantity_support():
#
#         plt.xticks(rotation=70)
#
#         where = EarthLocation.of_site(self.grb.site[loc])
#
#             if (grb.vis_tonight[loc]):
#                 altaz =  grb.radec.transform_to( AltAz(obstime=t,
#                                                  location=where ))
#                 plt.plot(t.datetime,altaz.alt,
#                          linestyle=":",
#                          color=color[loc])
#
#                 # For the given site, above 10°
#                 tstart = grb.t_true[loc][0]
#                 tstop  = grb.t_true[loc][1]
#                 duration = tstop - tstart
#                 t10 = tstart +  np.linspace(0,duration.value,100)
#                 altaz =  grb.radec.transform_to( AltAz(obstime=t10,
#                                                     location=where ))
#                 plt.plot(t10.datetime,altaz.alt, color=color[loc],label=loc)
#
#             else:
#                 print("Not visible in :",loc," within 24h after trigger")
#
#         if (grb.vis_tonight[site[0]] or grb.vis_tonight[site[1]]):
#             plt.hlines(y=10*u.deg,xmin=tmin.datetime,xmax=tmax.datetime,
#                        linestyles='--',colors="grey",label="Min alt.")
#             plt.title(grb.name)
#             plt.xlabel("Date (MM-DD-HH) UTC")
#             plt.ylim(ymin=0*u.deg,ymax=90*u.deg)
#             plt.legend()
#             plt.show()


        return
    ###############################################################################
    def plot(self, saveplots, outfile):
        """Plot simulations of the current grb"""

        if (saveplots == 0): return
        if (self.ndof != 0):
            sig = '$\sigma_{\sqrt{TS}}$'
        else:
            sig = '$\sigma_{Li&Ma}$'


        with quantity_support():

            
            # Compute relevant quantities
            eventid = self.grb.name +  " " + self.where
            t_s = [t[1] for t in self.tobs] # observation at end of interval
            smax = np.max(self.smax_list)
            err_smax = np.std(self.smax_list)
            t3s = [self.tobs[i][0].value
                     for i in self.id_3s_list
                     if self.tobs[i][0].value>=0]
            t5s = [self.tobs[i][0].value
                   for i in self.id_5s_list
                   if self.tobs[i][0].value>=0]
            
            fig = plt.figure(figsize=(15,15))
            fig.subplots_adjust(left   = None,
                                bottom = None,
                                right  = None,
                                #top    = None,
                                wspace = 0.4,
                                hspace = 0.3)

            ### Mean significance at each slice
            a1 = plt.subplot(221)
            a1.set_xlabel('Observation duration (s)')
            a1.set_ylabel(sig)
            a1.set_title(sig+' for '+str(self.niter)+' realisations')
            a1.set_xscale("log", nonposx='clip')
            a1.grid(which='both')

            label = " {}".format(eventid)
            a1.errorbar(t_s,
                        self.sigma_mean,
                        yerr=self.sigma_rms,
                        fmt='o',
                        label = label)
            a1.legend(fontsize=14,loc="upper left")
            
            ### Max significance in the run
            a2 = plt.subplot(222)
            a2.set_xlabel(sig)
            a2.set_ylabel('#Realisation')
            a2.set_title('Max. '+sig+' = '
                         + str( round(np.mean(self.smax_list),1))
                         + ' +/- '
                         + str( round(np.std(self.smax_list),1)))
            a2.grid(which='both')

            a2.hist(self.smax_list,
                    bins=int(self.niter/2),
                    range=[smax-3*err_smax,smax+3*err_smax],
                    color="grey",
                    label = eventid)
            
            a2.legend(fontsize=14)


            ### Plot S, B, S+B counts at  max significance
            a3 = plt.subplot(223)
            a3.set_xlabel('Log. of Counts')
            a3.set_ylabel('#Realisations')
            a3.set_title('S & B at max. '+sig)
            a3.grid(which='both')
            
            ntot = np.add(self.nex_smax_list,self.nb_smax_list)
            
            n_min = min(min(self.nex_smax_list),
                       min(self.nb_smax_list))
            n_max = max(ntot)
            ln_min = np.log10(n_min)
            ln_max = np.log10(n_max)
            
            n, bins, _ = a3.hist(np.log10(ntot),
                                 range=[ln_min,ln_max],
                                 bins=int(self.niter/2),
                                 color="tab:blue",
                                 alpha=0.5,
                                 label= MyLabel(ntot,
                                                eventid+"\nExcess+Bckgd"))
            
            a3.hist(np.log10(self.nb_smax_list), 
                    bins=bins,
                    color="tab:green",
                    alpha=0.5,
                    label= MyLabel(self.nb_smax_list,"Bckgd"))
            a3.hist(np.log10(self.nex_smax_list),
                    bins=bins,
                    color="tab:orange",
                    alpha=0.5,
                    label=MyLabel(self.nex_smax_list,
                                  "Excess"))
            plt.legend()

            ### Plot time to get 3/ 5 sigma for detected cases

            a5 = plt.subplot(224)
            a5.set_xlabel('Delay from trigger (s)')
            a5.set_ylabel('#Realisation')
#            title = "$T_{3\sigma}$ " \
#                  + "=  {:5.1f} +/- {:5.1f} ({:5.1f}%)" \
#                  .format(np.mean(t3s),np.std(t3s),
#                                 100*self.detect_3s)
#            a5.set_title(title)
            a5.grid(which='both')
            n, bins,_ = a5.hist(t3s,bins=10,alpha=0.5,label= MyLabel(t3s,"$t_{3\sigma}$"))
            a5.hist(t5s,bins=bins,alpha=0.5,label= MyLabel(t5s,"$t_{5\sigma}$"))
            a5.legend()

            plt.show()

            if (saveplots==2 or saveplots==3):
                fig.savefig(Path(outfile).with_suffix('.pdf'))
                # unsupported fig.savefig(outfile.with_suffix('.jpg'))
        return

    ###############################################################################
    def plot_trial(self,sigma,nex,nb):

        pvalue  = 1 - chi2.cdf(sigma**2, self.ndof)

        fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(12,15))
        ax1 = ax[0,0]
        ax2 = ax[0,1]
        ax3 = ax[1,0]
        ax4 = ax[1,1]
        ax1.plot(self.grb.t_s,pvalue,linestyle='-',marker='.')
        ax1.set_title("pvalue (the smaller the better)")
        ax1.set_xscale('log')

        ax2.plot(self.grb.t_s,sigma,linestyle='-',marker='.')
        ax2.set_title("sigma")
        ax2.set_xscale('log')

        ax3.plot(self.grb.t_s,nex,linestyle='-',marker='.')
        ax3.plot(self.grb.t_s,[x[2] for x in self.simulations],'o')
        ax3.set_title("nex")
        ax3.set_xscale('log')

        ax4.plot(self.grb.t_s,nb,linestyle='-',marker='.')
        ax4.plot(self.grb.t_s,[x[3] for x in self.simulations],'o')
        ax4.set_title("nb")
        ax4.set_xscale('log')

        plt.show()

        return
