# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import numpy as np

from gammapy.spectrum import SpectrumDatasetOnOff
from gammapy.spectrum import SpectrumEvaluator
from gammapy.spectrum.core import PHACountsSpectrum
from gammapy.utils.random import get_random_state

###########################################################################
def mc_onoff(mc = None,idx = 0, obs_param = None,debug=0):
     """
     Simulate observation with given parameters.
 
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
     livetime    = obs_param.livetime
     alpha       = obs_param.alpha
     emin        = obs_param.emin  # Check what is the diff with mc.Emin ...
     emax        = obs_param.emax
 
     # Compute expected counts
     bkg_mean_values = mc.perf[idx].bkg.data.data * livetime.to('s') # ThS - This is not a rate
 
     predicted_counts = SpectrumEvaluator(model    = mc.grb.spectra[idx],
                                          aeff     = mc.perf[idx].aeff,
                                          livetime = livetime,
                                          edisp    = mc.perf[idx].rmf)
 
     npred = predicted_counts.compute_npred()
 
     # set negative values to zero (interpolation issues)
     npred.data.data[ np.where(npred.data.data < 0.) ] = 0
 
     # Randomise counts
     random_state = get_random_state(mc.rnd_seed)
 
     # Default selection
     on_counts  = random_state.poisson(npred.data.data.value)  # excess
     bkg_counts = random_state.poisson(bkg_mean_values.value)  # bkg in ON region
     off_counts = random_state.poisson(bkg_mean_values.value / alpha)  # bkg in OFF region
      
     on_counts += bkg_counts  # evts in ON region

#     print(" mc_onoff {:3d} : On={:10.2f} bkg={:10.2f} off={:10.2f} dt={:10.2f}"
#           .format(idx,sum(on_counts), sum(bkg_counts),sum(off_counts),livetime))    
     

     on_vector = PHACountsSpectrum(data      = on_counts,
                                   backscal  =1,
                                   energy_lo = reco_energy[:-1],
                                   energy_hi = reco_energy[1:],
                                   )
     on_vector.livetime = livetime
 #        on_vector.plot()
 
     off_vector = PHACountsSpectrum(energy_lo=reco_energy[:-1],
                                    energy_hi=reco_energy[1:],
                                    data=off_counts,
                                    backscal=1. / alpha,
                                    is_bkg=True,
                                    )
     off_vector.livetime = livetime
 
     obs = SpectrumDatasetOnOff(counts    = on_vector,
                               counts_off = off_vector,
                               aeff       = mc.perf[idx].aeff,
                               edisp      = mc.perf[idx].rmf,
                               livetime   = livetime )
 
     # Gpy 0.9 - cannot be initialised here now, it seems.
     # Was it used anyhow ?
     # Forget it in the new version
     # obs.obs_id = obs_id
 
     # Set threshold according to the closest energy reco from bkg bins
     idx_min = np.abs(reco_energy - emin).argmin()
     idx_max = np.abs(reco_energy - emax).argmin()
     obs.lo_threshold = reco_energy[idx_min]
     obs.hi_threshold = reco_energy[idx_max]
 
     return obs