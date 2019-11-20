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
def simulate_onoff(mc,spectrum, pointing, obs_param,debug=0):
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

    livetime = obs_param.livetime
    alpha    = obs_param.alpha
    emin     = obs_param.emin
    emax     = obs_param.emax

    # Compute expected counts
    reco_energy = mc.perf.bkg.energy.edges
    bkg_mean_values = mc.perf.bkg.data.data * livetime.to('s') # ThS - This is not a rate

    predicted_counts = SpectrumEvaluator(model    = spectrum,
                                         aeff     = mc.perf.aeff,
                                         livetime = livetime,
                                         edisp    = mc.perf.rmf)

    npred = predicted_counts.compute_npred()

    # set negative values to zero (interpolation issue)
    idx = np.where(npred.data.data < 0.)
    npred.data.data[idx] = 0

    # Randomise counts
    random_state = get_random_state(mc.rnd_seed)

    # Default selection
    on_counts  = random_state.poisson(npred.data.data.value)  # excess
    bkg_counts = random_state.poisson(bkg_mean_values.value)  # bkg in ON region
    off_counts = random_state.poisson(bkg_mean_values.value / alpha)  # bkg in OFF region

# Try 95% containment radius (default radius in EventDisplay is 68%):
# The signal is scale by 0.95/0.68, background is scaled by b=r95**2/r68**2
# The radius r is given by sigma*np.sqrt(-2*log(1-alpha)) where alpha is
# the containment. Thereore b is scaled by log(1-0.95)/log(1-0.68) = 2.63
#        on_counts = random_state.poisson((0.95/0.68)*npred.data.data.value)  # excess
#        bkg_counts = random_state.poisson(2.63*bkg_mean_values.value)  # bkg in ON region
#        off_counts = random_state.poisson(2.63*bkg_mean_values.value / alpha)  # bkg in OFF region

    on_counts += bkg_counts  # evts in ON region
    on_vector = PHACountsSpectrum(data=on_counts,
                                  backscal=1,
                                  energy_lo=reco_energy[:-1],
                                  energy_hi=reco_energy[1:],
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
                              aeff       = mc.perf.aeff,
                              edisp      = mc.perf.rmf,
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
