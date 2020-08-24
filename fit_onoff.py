# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import numpy as np
import astropy.units as u

import ana_config as cf

from gammapy.spectrum import SpectrumDatasetOnOff
from gammapy.spectrum import SpectrumEvaluator
from gammapy.spectrum.core import PHACountsSpectrum
from gammapy.utils.random import get_random_state
from gammapy.stats.poisson import significance_on_off


###########################################################################
def mc_onoff(mc      = None,
             alpha   = 0.2,
             obstime = 0 *u.s,
             idx     = 0):
    
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
    reco_energy = mc.perf[idx].bkg.energy.edges # This is already in the MC object
    emin     = mc.Emin[idx]
    emax     = mc.Emax[idx]
    
    # Randomise counts later
    random_state = get_random_state(mc.rnd_seed)
 
    # Compute expected counts (not a rate!)
    bkg_mean_values = mc.perf[idx].bkg.data.data * obstime.to('s')
 
    if (cf.signal_to_zero == False):
        
        predicted_counts = SpectrumEvaluator(model    = mc.grb.spectra[idx],
                                             aeff     = mc.perf[idx].aeff,
                                             livetime = obstime,
                                             edisp    = mc.perf[idx].rmf)
 
        npred = predicted_counts.compute_npred()
 
        # set negative values to zero (interpolation issues)
        npred.data.data[ np.where(npred.data.data < 0.) ] = 0
        
        if (cf.do_fluctuate) : 
            on_counts  = random_state.poisson(npred.data.data.value)  # excess
        else:
            on_counts = npred.data.data.value
            
    else:
        on_counts = 0

    if (cf.do_fluctuate):
        bkg_counts = random_state.poisson(bkg_mean_values.value)
        off_counts = random_state.poisson(bkg_mean_values.value / alpha)
    else:
        bkg_counts = bkg_mean_values.value
        off_counts = bkg_counts / alpha
    
    on_counts += bkg_counts  # counts in ON region

#     print(" mc_onoff {:3d} : On={:10.2f} bkg={:10.2f} off={:10.2f} dt={:10.2f}"
#           .format(idx,sum(on_counts), sum(bkg_counts),sum(off_counts),livetime))    
     
    on_vector = PHACountsSpectrum(data      = on_counts,
                                  backscal  =1,
                                  energy_lo = reco_energy[:-1],
                                  energy_hi = reco_energy[1:],
                                  )
    on_vector.livetime = obstime
#   on_vector.plot()
 
    off_vector = PHACountsSpectrum(energy_lo=reco_energy[:-1],
                                   energy_hi=reco_energy[1:],
                                   data=off_counts,
                                   backscal=1. / alpha,
                                   is_bkg=True,
                                   )
    off_vector.livetime = obstime
 
    obs = SpectrumDatasetOnOff(counts     = on_vector,
                               counts_off = off_vector,
                               aeff       = mc.perf[idx].aeff,
                               edisp      = mc.perf[idx].rmf,
                               livetime   = obstime )
 
    # obs.obs_id = idx
 
    # Set threshold according to the closest energy reco from bkg bins
    idx_min = np.abs(reco_energy - emin).argmin()
    idx_max = np.abs(reco_energy - emax).argmin()
    obs.lo_threshold = reco_energy[idx_min]
    obs.hi_threshold = reco_energy[idx_max]
 
    return obs
    
###########################################################################
def cumulative_OnOffstats(simulations,n_min=10,alpha=0.2):
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