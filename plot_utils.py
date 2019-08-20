import os
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

from gammapy.stats.poisson import (significance_on_off,
                                   excess_error,
                                   background_error)

class PlotGammaRayBurst(object):
    """
    Class to plot GRB simulations result
    """

    @staticmethod
    def plot_stats_detection(grb, ax_excess=None, ax_relative=None, ax_bkg=None, ax_sigma=None, savefig=False, outdir='./out/'):
        stats = grb.get_cumulative_stats()

        import matplotlib.pyplot as plt
        fig = plt.figure(num=grb.name + ' stats', figsize=(18, 6))

        if ax_excess is None:
            ax_excess = plt.subplot2grid((1, 4), (0, 0))
        yerr = excess_error(
            stats['n_on'],
            stats['n_off'],
            stats['alpha']
        )  
        ax_excess.errorbar(stats['livetime'], stats['excess'],
                           yerr=yerr,
                           color='black', fmt='o')
        ax_excess.set_xlabel('Livetime [s]', fontweight='bold')
        ax_excess.set_ylabel('#Evts', fontweight='bold')
        ax_excess.set_title('Cumulated excess', fontweight='bold')
        ax_excess.set_ylim(0., (stats['excess'] + yerr).max() * (1.1))
        ax_excess.grid(which='both')
        ax_excess.set_xscale('log')
        
        if ax_relative is None:
            ax_relative = plt.subplot2grid((1, 4), (0, 1))
        yerr = excess_error(
            stats['n_on'],
            stats['n_off'],
            stats['alpha']
        )  
        ax_relative.errorbar(stats['livetime'], yerr/stats['excess'],
                           yerr=0,
                           color='black', fmt='o')
        ax_relative.set_xlabel('Livetime [s]', fontweight='bold')
        ax_relative.set_ylabel('Excess relative error', fontweight='bold')
        ax_relative.set_title('Relative error evolution', fontweight='bold')
        ax_relative.grid(which='both')
        ax_relative.set_xscale('log')        
        
        if ax_bkg is None:
            ax_bkg = plt.subplot2grid((1, 4), (0, 2))
        yerr = background_error(
            stats['n_off'],
            stats['alpha']
        )
        ax_bkg.errorbar(stats['livetime'], stats['bkg'],
                        yerr=background_error(
                            stats['n_off'],
                            stats['alpha']
                        ),
                        color='black', fmt='o')

        ax_bkg.set_xlabel('Livetime [s]', fontweight='bold')
        ax_bkg.set_ylabel('#Evts', fontweight='bold')
        ax_bkg.set_title('Cumulated background', fontweight='bold')
        ax_bkg.set_ylim(stats['bkg'].min() * 0.9, (stats['bkg'] + yerr).max() * (1.1))
        ax_bkg.grid(which='both')
        ax_bkg.set_xscale('log')
        ax_bkg.set_yscale('log')
        
        if ax_sigma is None:
            ax_sigma = plt.subplot2grid((1, 4), (0, 3))
        ax_sigma.errorbar(stats['livetime'], stats['sigma'],
                          color='black', fmt='o')
        ax_sigma.set_xlabel('Livetime [s]', fontweight='bold')
        ax_sigma.set_ylabel('Significance', fontweight='bold')
        ax_sigma.set_title('Significance (Li & Ma)', fontweight='bold')
        ax_sigma.set_ylim(0., stats['sigma'].max() * (1.1))
        ax_sigma.grid(which='both')
        ax_sigma.set_xscale('log') # CHANGED
        
        plt.tight_layout()
        if savefig == True:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            plt.savefig(outdir + '/' + grb.name + '.png')
        return ax_excess, ax_relative, ax_bkg, ax_sigma


    @staticmethod
    def make_gif_from_models(grb, emin=0.02 * u.TeV, emax=10 * u.TeV, savefig=False, outdir='./out/'):
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        def get_model(i):
            return grb.spectral_model[i]

        def get_time(i):
            return grb.time_interval[i]
        
        # Initialise plot
        fig_kw = dict(num=grb.name + ' models')
        fig, ax = plt.subplots(**fig_kw)
        model_init = get_model(1)
        fmin, fmax = model_init(energy=emax), model_init(energy=emin)
        x = np.linspace(emin.value, emax.value, 100) * u.TeV
        y = model_init(energy=x)
        ax.set(xlim=(emin.value,emax.value), ylim=(fmin.value, fmax.value))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy [TeV]')
        ax.set_ylabel('Flux [1 / (cm2 s TeV)]')
        ax.set_ylim([1.e-20, 1.e-6])
        ax.grid(which='both')
        line = ax.plot(x,y, color='k', lw=2)[0]

        def animate(i):
            model = get_model(i)
            y = model(energy=x)
            line.set_ydata(y)
            ax.set_title(grb.name + '; z={:.2f}; dt={}--{} s'.format(grb.z,
                                                                get_time(i)[0].value,
                                                                get_time(i)[1].value))
        anim = FuncAnimation(fig, animate, interval=500,
                             frames=len(grb.spectral_model))
        plt.draw()
        if savefig == True:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            anim.save(outdir + '/' + grb.name + '_animate.gif', writer='imagemagick')
            

