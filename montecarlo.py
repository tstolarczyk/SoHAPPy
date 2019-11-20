# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""

import numpy as np
import sys

import astropy.units as u

#import numpy as np

from scipy.stats import chi2
#from scipy.stats import norm

import matplotlib.pyplot as plt
from fit_3d import fit_3d, DataSetShow
from fit_onoff import simulate_onoff
from grb import GRBObservationParameters

from gammapy.spectrum.models import PowerLaw
from gammapy.maps import WcsNDMap

from gammapy.cube import MapDataset
from gammapy.cube.models import SkyModel, BackgroundModel
from gammapy.cube import make_map_exposure_true_energy, make_map_background_irf
from gammapy.cube import PSFKernel

__all__ = ['MonteCarlo']
###############################################################################
class MonteCarlo():

    def __init__(self,niter,fov,binsize,alpha,Emin,Emax,perf,grb,ndof=1):

        # Input parameters and object
        self.niter    = niter
        self.perf     = perf
        self.grb      = grb
        self.Emin     = Emin
        self.Emax     = Emax
        self.binsize  = binsize
        self.fov      = fov
        self.alpha    = alpha
        self.rnd_seed = 'random-seed'
        self.ndof     = ndof
        self.source_model = self.source_spectrum(ndof=ndof)
        self.null_model   = self.source_spectrum(ndof=ndof,
                                                 norm=0/(u.cm)**2/u.s/u.TeV)

        # list of simulations
        self.simulations = None

        # Max siginificance
        self.sigmax_list     = []
        self.nex_sigmax_list = []
        self.nb_sigmax_list  = []
        self.t_sigmax_list   = []
        # Mean values
        self.sigmax       = 0
        self.err_sigmax   = 0
        self.nex_sigmax   = -1
        self.nb_sigmax    = -1
        self.t_sigmax     = 0
        self.err_t_sigmax = 0

        ### 3 sigma detection
        self.nex_3sigma_list = []
        self.nb_3sigma_list  = []
        self.t_3sigma_list   = []
        # Mean values
        self.nex_3sigma   = -1
        self.nb_3sigma    = -1
        self.t_3sigma     = -1e9
        self.err_t_3sigma = -1e9

        ### 5 sigma detection
        self.nex_5sigma_list = []
        self.nb_5sigma_list  = []
        self.t_5sigma_list   = []
        # Mean values
        self.nex_5sigma    = -1
        self.nb_5sigma     = -1
        self.t_5sigma      = -1
        self.err_t_5sigma  = -1e9

        # Related to time and time bins
        self.detect_3sigma = 0  # Fraction of trials 3 sigma was reached
        self.detect_5sigma = 0  # Fraction of trials 5 sigma was reached

        # Mean sigma versus time plot
        self.sigma_mean = []
        self.sigma_rms  = []

        return

    ##########################################################################
    def source_spectrum(self,
                        ndof,
                        index = 2,
                        norm  = 1e-11/(u.cm)**2/u.s/u.TeV,
                        reference ="1 TeV",
                        debug=True):

        spectral_fitmodel = PowerLaw(index     = index,
                                     amplitude = norm,
                                     reference = reference)

        model = SkyModel(spatial_model  = self.grb.spatial_model,
                         spectral_model = spectral_fitmodel,
                         name           = "tobefitted")

        model.parameters["amplitude"].min = 0. # prevent negative signal

        # Fit spectral index
        if (ndof>1):
            model.parameters["index"].frozen     = False
        else:
            model.parameters["index"].frozen     = True

        # Fit source location
        if (ndof > 2):
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

        expo, bckg, psf, edisp = irf
        ### Signal and background theory (to generate a MC trial)
        bckg_theory = BackgroundModel(bckg)
        grb_theory  = SkyModel(spatial_model  = mc.grb.spatial_model,
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
    def run(self, ndof, dbg_level):
        """
        - Run simulations of the current grb - Likelihood method

        - Compute :
            - significance on the full observation,
            - significance max on intervals,
            - time of siginificance max,
            - time of 3 or 5 sigma reached
        """

        iMC=1
        prt_frequency = 10

        n_3sigma = 0 # Number of time 3 sigma was reached
        n_5sigma = 0 # Number of time 5 sigma was reached

        sigma_sum = 0
        sigma2_sum = 0

        while(iMC <= self.niter):
            if (self.niter <= 10) or (np.mod(iMC,prt_frequency) == 0):
                print("#",iMC," ",end="")

            ### Loop on GRB time slices
            self.loop_slices(ndof, dbg_level)

            if (ndof != 0):
                ################################ 3D
                fcnsb = np.array([x[0] for x in self.simulations]) # lik_fit
                fcnb  = np.array([x[1] for x in self.simulations]) # lik_bonly
                nex   = np.array([x[2] for x in self.simulations])
                nb    = np.array([x[3] for x in self.simulations])
                # ns_bonly - not used
                # nb_bonly - not used

                TS =  fcnb - fcnsb
                TS[TS < 1e-3] = 0 # Avoid sligtly negative value

                # pvalue  = 1 - chi2.cdf(TS, ndof)
                # sigma   = norm.ppf(1-pvalue)  # Integral from -Inf, can be<0
                # print(" Warning sqrt TS")
                sigma = np.sqrt(TS)
                sigma[sigma < 0 ] = 0
                sigma[sigma > 99] = 99

                if (dbg_level>1): self.plot_trial(sigma,nex,nb,ndof)

            if (ndof == 0):
                ################################ On-Off
                ### Get cumulated statistics
                obs_stat = self.grb.get_cumulative_stats(self.simulations)
                sigma = obs_stat['sigma']
                non        = obs_stat['n_on']
                noff       = obs_stat['n_off']

                nex = non -self.alpha*noff
                nb  = self.alpha*noff

            # End of test on likelihood case

            ### Find maximum
            maxidx     = np.argmax(sigma)
            sigmax     = sigma[maxidx]
            nex_sigmax = nex[maxidx]
            nb_sigmax  = nb[maxidx]
            t_sigmax   = self.grb.t_s[maxidx]

            self.t_sigmax_list.append(t_sigmax)
            self.nex_sigmax_list.append(nex_sigmax)
            self.nb_sigmax_list.append(nb_sigmax)
            self.sigmax_list.append(sigmax)

            # Find when 3 sigma is reached
            ts3sid = np.where(sigma>=3)
            if (np.size(ts3sid) != 0):
                ts3sid = ts3sid[0][0]
                t_3sigma   = self.grb.t_s[ts3sid]
                nex_3sigma = nex[ts3sid]
                nb_3sigma  = nb[ts3sid]
                n_3sigma  +=1
            else:
                t_3sigma = -1e9
                nex_3sigma = -1
                nb_3sigma = -1

            self.t_3sigma_list.append(t_3sigma)
            self.nex_3sigma_list.append(nex_3sigma)
            self.nb_3sigma_list.append(nb_3sigma)

            # Find when 5 sigma is reached
            ts5sid = np.where(sigma>=5)
            if (np.size(ts5sid) != 0):
                ts5sid = ts5sid[0][0]
                t_5sigma   = self.grb.t_s[ts5sid]
                nex_5sigma = nex[ts5sid]
                nb_5sigma  = nb[ts5sid]
                n_5sigma  += 1
            else:
                t_5sigma = -1e9
                nex_5sigma = -1
                nb_5sigma = -1

            self.t_5sigma_list.append(t_5sigma)
            self.nex_5sigma_list.append(nex_5sigma)
            self.nb_5sigma_list.append(nb_5sigma)

            # Acummulate sum and sum**2 for mean / error
            sigma_sum  += sigma
            sigma2_sum += sigma**2

            if (dbg_level):
                print(" t(sigmax) = {:5.2f}, sig= {:5.2f} at slice {:2d}"
                      .format(t_sigmax,sigmax,maxidx+1))
                print(" t(3s)     = {:5.2f}".format(t_3sigma))
                print(" (t5s)     = {:5.2f}".format(t_5sigma))

            iMC+=1
            # End of MC loop

        ### Compute mean sigma and rms for all trials
        self.sigma_mean = sigma_sum/self.niter
        sigma2_mean     = sigma2_sum/self.niter
        self.sigma_rms  = np.sqrt(sigma2_mean-self.sigma_mean**2)
        #print("Mean sigma =",self.sigma_mean,"rms=",self.sigma_rms)

        ### Compute mean values and erros
        self.sigmax         = round(np.mean(self.sigmax_list),3)
        self.err_sigmax     = round(np.std( self.sigmax_list),3)
        self.nex_sigmax     = round(np.std( self.nex_sigmax_list),3)
        self.nb_sigmax      = round(np.std( self.nb_sigmax_list),3)
        self.t_sigmax       = round(np.mean(self.t_sigmax_list),3)
        self.err_t_sigmax   = round(np.std( self.t_sigmax_list),3)


        t_detected   = [time for time in self.t_3sigma_list   if time >= 0]
        nex_detected = [nex  for nex  in self.nex_3sigma_list if nex>=0]
        nb_detected  = [nb   for nb   in self.nb_3sigma_list  if nb>=0]
        if (len(t_detected)):
            self.nex_3sigma     = round(np.std(nex_detected),3)
            self.nb_3sigma      = round(np.std(nb_detected),3)
            self.t_3sigma       = np.mean(t_detected)
            self.err_t_3sigma   = np.std( t_detected)

        t_detected =   [time for time in self.t_5sigma_list if time >= 0]
        nex_detected = [nex  for nex  in self.nex_5sigma_list if nex>=0]
        nb_detected  = [nb   for nb   in self.nb_5sigma_list  if nb>=0]


        if (len(t_detected)):
            self.nex_5sigma     = round(np.std(self.nex_5sigma_list),3)
            self.nb_5sigma      = round(np.std(self.nb_5sigma_list),3)
            self.t_5sigma       = round(np.mean(t_detected),3)
            self.err_t_5sigma   = round(np.std( t_detected),3)

        ### Compute detection level
        self.detect_3sigma = n_3sigma/self.niter
        self.detect_5sigma = n_5sigma/self.niter


        return

    ###########################################################################
    #
    ###########################################################################
    def loop_slices(self,ndof,debug_level):
        """
        Run simulation - OnOff analysis

        Loop over time interval and create a target, a set of observation parameters.
        Run the simulation (includes the On-Off analysis)
        """
        cumulated_npred = None
        livetime = 0
        self.simulations = []
        for idx, interval in enumerate(self.grb.time_interval):


            if (debug_level>1): print("*===> Simulating interval ",idx)

            livetime += interval[1] - interval[0] # Cumulate livetime

            obs_param = GRBObservationParameters(alpha    = self.alpha,
                                                 livetime = livetime,
                                                 emin     = self.Emin,
                                                 emax     = self.Emax,
                                                 pointing = self.grb.pointing,
                                                 fov      = self.fov,
                                                 binsize  = self.binsize)

            if (ndof != 0):
                ### Get IRF for the current observation
                irf = self.get_observation_irf(obs_param,
                                               self.grb.pointing,
                                               self.perf,
                                               debug_level)

                ### Get count map
                npred = self.theory_prediction(irf,
                                               self.grb.spectral_model[idx])
                if (cumulated_npred == None):
                    cumulated_npred = npred
                else:
                    cumulated_npred += npred
                ### Generate MC trial
                #rng = np.random.RandomState(seed=51) # same events at each trial
                rng    = np.random.RandomState()
                counts = rng.poisson(cumulated_npred.data) # Compute counts
                counts_map = WcsNDMap(obs_param.geom, counts)

                simu  =  fit_3d(self,
                                counts_map,
                                irf,
                                debug        = debug_level)

            else:
                simu = simulate_onoff(self,
                                      spectrum     = self.grb.spectral_model[idx],
                                      pointing     = self.grb.pointing,
                                      obs_param    = obs_param,
                                      debug        = debug_level)

            self.simulations.append(simu)

            if (ndof == 0): self.stack_obs = self.grb.add_stack_OnOffobs(simu)

        # This goes back to the Monte Carlo module to treat the
        # information collected in this simulation
        return

    ###########################################################################
    #
    ###########################################################################
    def result(self,out=sys.stdout):
        """Print out results"""

        print("\n========================== RESULTS",file=out)
        print("  GRB name               : ",self.grb.name,file=out)
        print("  Redshift             z : ",self.grb.z,file=out)
        print("==========================  ",file=out)
        print("  Max. significance       : ",round(self.sigmax,1),
              "+/-",round(self.err_sigmax,1),file=out)
        print("                  at time : ",round(self.t_sigmax,1),
              "+/-",round(self.err_t_sigmax,1),file=out)
        print("  3 sigma reached at time : ",round(self.t_3sigma,1),
              "+/-",round(self.err_t_3sigma,1),file=out)
        print("                  Fraction:",100*self.detect_3sigma,"%",file=out)
        print("  5 sigma reached at time : ",round(self.t_5sigma,1),
              "+/-",round(self.err_t_5sigma,1),file=out)
        print("                  Fraction:",100*self.detect_5sigma,"%",file=out)
        print("==========================\n",file=out)

        return

    ###############################################################################
    #
    ###############################################################################
    def plot(self, saveplots, out_filename):
        """Plot simulations of the current grb"""

        if (saveplots == 0): return
        if (self.ndof != 0):
            sig = '$\sigma_{\sqrt{TS}}$'
        else:
            sig = '$\sigma_{Li&Ma}$'

        # Prepare plots
        fig = plt.figure(figsize=(12,18))
        fig.subplots_adjust(left   = None,
                            bottom = None,
                            right  = None,
                            top    = None,
                            wspace = 0.3,
                            hspace = 0.3)

        ### Mean significance at each slice
        a1 = plt.subplot(321)
        a1.set_xlabel('Observation duration (s)')
        a1.set_ylabel(sig)
        a1.set_title(sig+' for '+str(self.niter)+' realisations')
        a1.set_xscale("log", nonposx='clip')
        a1.grid(which='both')
        a1.errorbar(self.grb.t_s,self.sigma_mean,yerr=self.sigma_rms,fmt='o')

        ### Max significance in the run
        a2 = plt.subplot(322)
        a2.set_xlabel(sig)
        a2.set_ylabel('#Realisation')
        a2.set_title('Max. '+sig+' = '
                     + str( round(self.sigmax,1))
                     + ' +/- '
                     + str( round(self.err_sigmax,1)))
        a2.grid(which='both')
        smax = int(max(self.sigmax_list))+1
        a2.hist(self.sigmax_list,bins=25,range=[0,smax],color="grey")
        #y,binEdges = np.histogram(sigmax_list)
        #bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        #plt.bar(bincenters, y, yerr=np.sqrt(y))

        ### Plot S, B, S+B counts at  max significance
        a3 = plt.subplot(323)
        a3.set_xlabel('Log. of Counts')
        a3.set_ylabel('#Realisation')
        a3.set_title('S & B at max. '+sig)
        a3.grid(which='both')
        ltot = np.log10(np.add(self.nex_sigmax_list,self.nb_sigmax_list))
        n, bins, p_ = a3.hist(ltot,
                              range = [0,max(ltot)+1],
                              bins=25,
                              color="grey",
                              alpha=0.5,
                              label="Excess+Bckgd")
        a3.hist(np.log10(self.nb_sigmax_list), bins=bins,
                color="pink",
                alpha=0.5,
                label="Bckgd")
        a3.hist(np.log10(self.nex_sigmax_list),bins=bins,
                color="green",
                alpha=0.5,
                label="Excess")
        plt.legend()
        
        # Plot correlation between max significance and time 
        a4 = plt.subplot(324)
        a4.set_ylabel('Observation duration')
        a4.set_xlabel('$\sigma_{Max}$')
        a4.set_title('Duration versus $\sigma$ correlation')
        a4.grid(which='both')
        a4.set_yscale("log")
        a4.set_ylim([10,1e4])
        a4.scatter(self.sigmax_list,self.t_sigmax_list,marker='+')
    

        ### Plot time to get 3 sigma for detected cases
        a5 = plt.subplot(325)
        a5.set_xlabel('Observation duration (s)')
        a5.set_ylabel('#Realisation')
        a5.set_title('$T_{3\sigma}$ = '
                     + str( round(self.t_3sigma,1))
                     + ' +/- '
                     + str( round(self.err_t_3sigma,1))
                     + ' ('
                     + str(round(100*self.detect_3sigma,1))
                     +'%)')
        a5.grid(which='both')
        t_detected = [time for time in self.t_3sigma_list if time >= 0]
        a5.hist(t_detected, color="grey")

        ### Plot time to get 5 sigma for detected cases
        a6 = plt.subplot(326)
        a6.set_xlabel('Observation duration (s)')
        a6.set_ylabel('#Realisation')
        a6.set_title('$T_{5\sigma}$ = '
                     + str( round(self.t_5sigma,1))
                     + ' +/- '
                     + str( round(self.err_t_5sigma,1))
                     + ' ('
                     + str(round(100*self.detect_5sigma,1))
                     +'%)')
        a6.grid(which='both')
        t_detected = [time for time in self.t_5sigma_list if time >= 0]
        a6.hist(t_detected, color="grey")



        if (saveplots==1 or saveplots==2): plt.show()
        if (saveplots==2 or saveplots==3):
            #fig.savefig(out_filename+'.eps',transparent=True)
            fig.savefig(out_filename+'.jpg')
        return

    ###############################################################################
    def plot_trial(self,sigma,nex,nb,ndof):


        pvalue  = 1 - chi2.cdf(sigma**2, ndof)

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

    ###############################################################################
    def get_observation_irf(self,obs_param, pointing, irfs, debug=False):

        """Get the IRf components"""

        axis = obs_param.axis
        duration = obs_param.livetime
        geom = obs_param.geom
        fov = obs_param.fov

        offset = 0*u.deg

        ### Compute the exposure map
        exposure = make_map_exposure_true_energy(
                    pointing = pointing.galactic, # does not accept SkyCoord altaz
                    livetime = duration,
                    aeff     = irfs["aeff"],
                    geom     = geom)

        ### Compute the background map
        background = make_map_background_irf(
                        pointing = pointing.galactic, # does not accept SkyCoord altaz
                        ontime   = duration,
                        bkg=irfs["bkg"],
                        geom=geom
                        )

        ### Point spread function
        psf = irfs["psf"].to_energy_dependent_table_psf(theta=offset)
        psf_kernel = PSFKernel.from_table_psf(psf, geom, max_radius=fov/2) # Why this max_radius ?

        ### Energy dispersion
        edisp = irfs["edisp"].to_energy_dispersion(offset,
                    e_reco=axis.edges,
                    e_true=axis.edges)

        if (debug>2):
            islice = 3

            exposure.slice_by_idx({"energy": islice-1}).plot(add_cbar=True);
            # fig, ax = plt.subplots(ncols=nedges-1,nrows=1,figsize=((nedges-1)*2,4))
            # islice = 1
            # while (islice < nlogedges):
            #     exposure.slice_by_idx({"energy": islice-1}).plot(ax=ax[islice-1],add_cbar=True);
            #     islice += 1
            plt.title("Exposure map")

            background.slice_by_idx({"energy": islice}).plot(add_cbar=True);
            # fig, ax = plt.subplots(ncols=nedges-1,nrows=1,figsize=((nedges-1)*2,4))
            # islice = 1
            # while (islice < nedges):
            #     background.slice_by_idx({"energy": islice-1}).plot(ax=ax[islice-1],add_cbar=True);
            #     islice += 1
            plt.title("Background map")

            psf_kernel.psf_kernel_map.sum_over_axes().plot(stretch="log");
            plt.title("point spread function")

            edisp.plot_matrix();
            plt.title("Energy dispersion matrix")
            plt.show()

        return exposure, background, psf_kernel, edisp
