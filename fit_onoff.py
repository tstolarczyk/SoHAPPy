# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import sys
import os
import numpy as np
import ana_config as cf

# # Avoid deprecation Astropy warnings in gammapy.maps
import warnings
with warnings.catch_warnings():
    import gammapy
    if (gammapy.__version__=="0.12"):
        from gammapy.spectrum import SpectrumDatasetOnOff
        from gammapy.spectrum import SpectrumEvaluator
        from gammapy.spectrum.core import PHACountsSpectrum
        from gammapy.utils.random import get_random_state
        from gammapy.stats.poisson import significance_on_off


###########################################################################
def mc_onoff(aslice, mc = None, alpha   = 0.2, debug=False):
    
    """
    Simulate one observation (slice idx) with given parameters.

    Compute predicted number of signal and background counts in an observation 
    region, from the livetime.
    Fluctuate the count numbers (except if requested to not to do so).    
    
    Accept to specific flags:
        - :math:`set_signal_to_zero` will simulate background only 
        (signal stricltly force to zero, no fluctuations)
        - math:`do_fluctuate` will apply the Poisson fluctuation
        
    Return the on-region and off-region respective counts as an observation

    Parameters

    perf : `~gammapy.scripts.CTAPerf`
        CTA performance

    target : `~gammapy.scripts.Target`
        Source

    obs_param : `~gammapy.scripts.ObservationParameters`
        Observation parameters

    obs_id : `int`, optional
        Observation Id

    random_state : {int, 'random-seed', 'global-rng', `~numpy.random.RandomState`}
        Defines random number generator initialisation.
        Passed to `~gammapy.utils.random.get_random_state`.

    """
    if (gammapy.__version__ == "0.16"): return
    
    #debug=True
        
    obstime      = aslice.ts2()-aslice.ts1()
    idx          = aslice.idt()
    spectrum     = mc.slot.grb.spectra[aslice.fid()]
    #spectrum     = mc.slot.grb.spectra[idx] # This has been like this :-()
    random_state = get_random_state(mc.rnd_seed) # Randomise counts later
    
    # Note : the number of on/off counts versus energy is strictly related to
    # background energy bins (=energy egdes -1)
    perf_list   = aslice.perf()
    reco_energy = perf_list[0].bkg.energy.edges
    nEbins      = len(reco_energy) -1    
    
    tot_on_counts  = np.zeros(nEbins)
    tot_off_counts = np.zeros(nEbins)
    
    if (debug):
        if (idx==0):
            print("\n")
            print("{:>2s} {:>8s} {:>4s} {:>4s} {:>10s} {:>10s} {:>10s} {:>10s} {:>s}"
                  .format("","dt (s)","perf","flux","on","off","tot_on","tot_off","Alt")  )      

    for ip,perf in enumerate(perf_list):
        
        if (len(perf.bkg.energy.edges)-1 != nEbins):
            print(" Perf E-bins assumed is ",nEbins,
                  "actual is ",len(perf.bkg.energy.edges))
            sys.exit("Fatal error in the IRF E binning")
    
        ### BACKGROUND / OFF COUNTS
        bkg_mean_values = perf.bkg.data.data * obstime.to('s')
        if (cf.do_fluctuate):
            bkg_counts = random_state.poisson(bkg_mean_values.value)
            off_counts = random_state.poisson(bkg_mean_values.value / alpha)
        else:
            bkg_counts = bkg_mean_values.value
            off_counts = bkg_counts / alpha
                
        ### SIGNAL / EXCESS / ON COUNTS

        if (cf.signal_to_zero == False):
        # Compute signal
        
            predicted_counts = SpectrumEvaluator(model    = spectrum,
                                                 aeff     = perf.aeff,
                                                 livetime = obstime,
                                                 edisp    = perf.rmf) 
            npred = predicted_counts.compute_npred()

            # set negative values to zero (interpolation issues)
            npred.data.data[ np.where(npred.data.data < 0.) ] = 0
        
            if (cf.do_fluctuate) : 
                on_counts  = random_state.poisson(npred.data.data.value)
            else:
                on_counts = npred.data.data.value
        else:
        # Force signal to zero
            on_counts = 0

        #on_counts += bkg_counts  # counts in ON region
        on_counts = on_counts + bkg_counts  # counts in ON region
        
        # CUMULATE the counts of both sites if it is the case (perf)
        # tot_on_counts += on_counts
        # tot_off_counts += off_counts
        tot_on_counts  = tot_on_counts  + on_counts
        tot_off_counts = tot_off_counts + off_counts
        
        if(debug):        
            if (ip == 0):  
                print("{:2d} {:8.2f} ".format(idx,obstime.value),end="")
            else: 
                print("{:2s} {:8s} ".format("",""),end="")
            if (mc.slot.site != "Both"):
                alt = mc.slot.grb.altaz(loc=mc.slot.site,dt=aslice.tobs()).alt
            else: 
                alt = -99
            print("{:4d} {:3d} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:5.2f}  {}"
                  .format(ip,
                          aslice.fid(),
                          sum(on_counts),
                          sum(off_counts),
                          sum(tot_on_counts),
                          sum(tot_off_counts),
                          alt,
                          os.path.basename(perf.name)))
                 
        # End of loop over performance - cumulate on/off count is the slice
        
    on_vector = PHACountsSpectrum(data      = tot_on_counts,
                                  backscal  = 1,
                                  energy_lo = reco_energy[:-1],
                                  energy_hi = reco_energy[1:],
                                  )
    on_vector.livetime = obstime
 
    off_vector = PHACountsSpectrum(energy_lo=reco_energy[:-1],
                                    energy_hi=reco_energy[1:],
                                    data=tot_off_counts,
                                    backscal=1. / alpha,
                                    is_bkg=True,
                                    )
    off_vector.livetime = obstime
 
    # We take as perf the firts one - this is wrong when more than site
    # will have to be investigated
    obs = SpectrumDatasetOnOff(counts      = on_vector,
                                counts_off = off_vector,
                                aeff       = perf_list[0].aeff,
                                edisp      = perf_list[0].rmf,
                                livetime   = obstime )
 
    #obs.obs_id = idx

    # Set threshold according to the closest energy reco from bkg bins
    emin = min(reco_energy)
    emax = max(reco_energy)
    idx_min = np.abs(reco_energy - emin).argmin()
    idx_max = np.abs(reco_energy - emax).argmin()
    obs.lo_threshold = reco_energy[idx_min]
    obs.hi_threshold = reco_energy[idx_max]
    
    #debug=False
    if (debug):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8))
        on_vector.plot( ax=ax,show_energy=emin,label="On "+str(idx))
        off_vector.plot(ax=ax,show_energy=emin,label="Off "+str(idx))
        ax.legend()
        
        
    return obs
    
###########################################################################
def cumulative_OnOffstats(simulations,n_min=10,alpha=0.2,debug=False):
    """
    Get cumulative statistics of the current MC simulation, composed of 
    observations corresponding to the GRB time slices.
    If math:`non` or math:`noff` are below the chosen limit, the significance is not 
    trustable and set to math:`nan` (Not a number).
    """

    # Init vectors
    nsim      = len(simulations)
    tot_time  = np.zeros(nsim)
    tot_n_on  = np.zeros(nsim)
    tot_n_off = np.zeros(nsim)
    tot_alpha = np.zeros(nsim)
    delta_t   = np.zeros(nsim)

    # Loop on observations
    for idx, obs in enumerate(simulations):

        alpha = obs.alpha
        n_on  = obs.total_stats_safe_range.n_on
        n_off = obs.total_stats_safe_range.n_off
        # Possible only if a model is set
        # obs.plot_counts()
        
        if idx == 0:
            tot_time[idx]  = obs.total_stats_safe_range.livetime.value
            tot_n_on[idx]  = n_on
            tot_n_off[idx] = n_off
            tot_alpha[idx] = obs.total_stats_safe_range.alpha
        else:
            tot_time[idx]  += tot_time[idx-1] + obs.total_stats_safe_range.livetime.value
            tot_n_on[idx]  += tot_n_on[idx-1] + n_on
            tot_n_off[idx] += tot_n_off[idx-1] + n_off
            tot_alpha[idx]  = obs.total_stats_safe_range.alpha

            delta_t[idx] =  obs.total_stats_safe_range.livetime.value
            
    # Find problematic slices if any
    idon  = np.where(tot_n_on  <= n_min)[0] # Indices with problem
    idoff = np.where(tot_n_off <= n_min)[0] # Indices with problem
    idrm  = np.hstack([idon,idoff])        # Concatenate
    idrm  = np.asarray(list(set(idrm)))    # Sort, simplify
    
    # Compute excess,background and significance from the cumulated counts
    tot_excess       = tot_n_on - alpha * tot_n_off
    tot_bkg          = alpha * tot_n_off
    tot_significance = significance_on_off(tot_n_on,tot_n_off,tot_alpha)
    
    # Set significiance to nan if not trustable
    if (len(idrm) != 0): tot_significance[idrm] = np.nan
#            print("Li&Ma not applicable at slices :",idrm)
#            print("    - Non  = ",tot_n_on[idrm])
#            print("    - Noff = ",tot_n_off[idrm])
#            print("    - Sig  = ",tot_significance[idrm])

    return dict(livetime = tot_time,
                excess   = tot_excess,
                bkg      = tot_bkg,
                sigma    = tot_significance,
                n_on     = tot_n_on,
                n_off    = tot_n_off,
                alpha    = tot_alpha,
                delta_t  = delta_t)    