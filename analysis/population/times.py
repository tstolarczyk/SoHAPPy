# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 17:11:21 2022.

Characterizes the time needed to get a detection or the mean maximal
significance. Plots in particular for...

@author: Stolar
"""
import sys

import numpy as np
import astropy.units as u
from astropy.visualization import quantity_support

import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from niceplot import MyLabel, stamp, col_size, vals_legend, \
                     show_noticeable_times
from niceprint import heading, t_fmt

from population import Pop
from pop_io import get_data


# Bigger texts and labels
sns.set_context("notebook")  # poster, talk, notebook, paper

__all__ = ["sigmx_after_delay", "time_to_sigma",
           "time_to_detection", "sigmx_versus_tmx",
           "sig_vs_time_errors",
           "delay_correlation_both_sites"]

codefolder = "../../"
sys.path.append(codefolder)


# ##----------------------------------------------------------------------------
def sigmx_after_delay(pop, sigmx_min=5, tmin=1*u.d):
    """
    Display mean maximal significances reached after a certain time.

    Significances are above a minimal value and time ae counted from the
    explosion. It also displays the fraction of trials reaching the
    thershold.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    sigmx_min : float, optional
        Mean siginificance value to be considered. The default is 5.
    tmin : Astropy Time, optional
        Mean duration to be considered. The default is 1*u.d.

    Returns
    -------
    None.

    """
    heading(f"Max significance after {t_fmt(tmin):}")

    print(f" Events above {sigmx_min:} (meanvalue), "
          f"{t_fmt(tmin):} after the trigger")

    for g, txt in zip([pop.g_n, pop.g_s, pop.g_b, pop.g_tot],
                      ["North", "South", "Both", "All"]):

        mask = g.sigmx > sigmx_min
        late = g[mask][g[mask].tmx > tmin.to(u.s).value]

        print(f"\n----- {txt:5s}: {len(late):4d}")

        print(f"{'Source':<25s}: {'Date':>10s}   {'sigma':>5s}   {'%':>5s}")

        for s, t, sig, d5s in zip(late.name,
                                  late.tmx.values*u.s,
                                  late.sigmx.values,
                                  late.d5s.values):

            print(f"{s:25s}: {t_fmt(t):>8.2f}   {sig:>5.1f}   {d5s:>5.2f}",
                  end="")

            if d5s >= pop.eff_lvl:
                print(" ***")
            elif sig >= 3:
                print(" *")
            else:
                print()


# ##----------------------------------------------------------------------------
def time_to_sigma(gpop, color="grey", tag="Unknown",
                  sigmx=True, ax=None, sigmin=5, alpha=0.3, lvl=90):
    """
    Plot time for reaching 5 sigma or sigma max above a minimum value.

    Parameters
    ----------
    gpop : pandas data frame
        Current population
    sigmx : boolean, optional
        If True, plots for sigma max instead of 5 sigma. The default is True.
    alpha : float, optional
        Plot transparency. The default is 0.3.

    Returns
    -------
    None.

    """
    if ax is None:
        fig, ax = plt.subplots()
        first = True
    else:
        first = False

    if sigmx:
        var = gpop[gpop.sigmx >= sigmin].tmx
        xlabel = (rf"$Mean \ Time \ to \ reach \ mean \ \sigma_{{max}} \ "
                  rf"\geq {{{sigmin:.0f}}} $")
    else:
        var = gpop[gpop.d5s > lvl].t5s
        xlabel = (rf"$Mean \ Time \ to \ reach \ 5 \sigma \ at  "
                  rf"\ {lvl:.0f} \% \ CL $")

    # Log bins, specific times and labels
    bins = np.logspace(2, np.log10(168*3600), 50)

    with quantity_support():
        ax.hist(var,  bins=bins, color=color, alpha=alpha, label=tag)

        axx = ax.twinx()
        hist, edges = np.histogram(var, bins=bins)
        xcenter = 0.5*(edges[1:] + edges[:-1])
        vals = np.cumsum(hist)/np.sum(hist)
        axx.plot(xcenter, vals, marker="", color=color,
                 alpha=min([1, 2*alpha]))
        axx.grid(axis="y")
        axx.set_ylabel("Fraction detected")
        # ax.set_yscale("log")

        if first:
            first = False
            ax.set_xscale("log")
            ax.set_xlabel(xlabel + r"$\ (s)$")
        else:
            ax.legend(loc="upper center")
            show_noticeable_times(axx, vpos=1.02*axx.get_ylim()[1])
            stamp(pop.tag[0], where="left")
    return ax


# ##----------------------------------------------------------------------------
def time_to_detection(sub_pop, binw=1, yscale="log", title="generic",
                      **kwargs):
    """
    Display the distribution of times for the 5 sigma maximal significance.

    Shows the delay starting from the burst or the start of the detection.
    All times are in hours.
    The plots says how log ot is necessary to wait for such detections (median
    values)

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    binw : float, optional
        Histogram bin width in hour(s). The default is 1.
    yscale : string, optional
        Log or Linear for y-scale. The default is "log".
    **kwargs : Dictionnary
        Extra parameters.

    Returns
    -------
    None.

    """
    varlist = [sub_pop.t5s/3600,
               sub_pop.tmx/3600]  # Go to hours
    masks = [(sub_pop.d5s > pop.eff_lvl),
             (sub_pop.tmx >= 0) & (sub_pop.sigmx >= 5)]
    tags = [r"$5\sigma$",
            r"$\sigma_{max} \geq 5$"]

    tobs = sub_pop.t1/3600  # That's the observation start

    # Interesting but useless here
    if masks is None:
        masks = np.ones(len(varlist)).astype(bool)

    # Define one hour binning from the min and max
    tmax = int(max([max(v[m]) for v, m in zip(varlist, masks)]))
    bins = range(0, tmax + 1, binw)
    # print(" Max time :",tmax, " bins:",bins)

    fig, ax = plt.subplots(nrows=len(varlist), ncols=1,
                           figsize=(12, 3.5*len(varlist)), sharex=True)

    first = True
    for ax0, tag, var, mask in zip(ax, tags, varlist, masks):

        _, bins, _ = ax0.hist(var[mask], bins=bins, alpha=0.5,
                              label=MyLabel(var[mask],
                                            label=tag+" (Burst)",
                                            stat="med"), **kwargs)
        if len(tobs) != 0:
            ax0.hist(var[mask]-tobs[mask], bins=bins, alpha=0.5,
                     label=MyLabel(var[mask]-tobs[mask],
                                   label=tag+" (Start)",
                                   stat="med"), **kwargs)

        ax0.set_yscale(yscale)

        ax0.legend()
        ax0.grid(which='both')
        minor_ticks = bins
        ax0.set_xticks(minor_ticks, minor=True)

        ax0.grid(which="major", ls=":")
        ax0.grid(which="minor", ls=":", alpha=0.5)
        ax0.set_ylabel("$h^{-1}$")
        if first:
            ax0.set_title(title)
            first = False

        axx = inset_axes(ax0, width="50%", height=1.2, loc="upper center")
        axx.hist(var[mask][var[mask] < 1]*60,
                 bins=range(0, 60, 1), alpha=0.5)
        axx.axvline(107/60, ls=":", label=" Alert + CTAO delays")
        axx.set_xlabel(r"$\Delta t$ (min)")
        axx.set_ylabel(r"$min^{-1}$")
        axx.legend()

    ax[-1].set_xlabel("Mean detection time (h)")  # Last one
    stamp(pop.tag[0], axis=fig)
    plt.tight_layout()


# ##---------------------------------------------------------------------------
def sigmx_versus_tmx(gpop, tag="unknown", axis="sigmx"):
    """
    Display mean maximal significance versus the mean time to reach it.

    Parameters
    ----------
    gpop : Pop instance
        Data on disk
    tag : String, optional
        Data descriptor. The default is "".

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))

    xlabel = r"$Mean \ Time \ to \ reach \ mean \ \sigma_{max} \ ($"

    with quantity_support():

        show_noticeable_times(ax)
        ax.axhline(3, ls=":", color="black", lw="2")
        ax.axhline(5, color="black", lw="2")
        ax.text(100, 3*0.50, r"$3\sigma$", ha="right")
        ax.text(100, 5*1.5, r"$5\sigma$", ha="right")

        # Main plots
        if axis == "sigmx":
            mask = gpop.sigmx > -np.Inf
            colors, sizes = col_size(gpop[mask].sigmx)
            title = r"$\sigma_{max}$"
            patches = vals_legend()

        elif axis == "z":
            mask = gpop.z > 0
            colors, sizes = col_size(gpop[mask].z,
                                     var_min=0.1,
                                     var_max=5,
                                     scale=100, log=False,
                                     colormap="cool")
            title = r"$redshift$"
            patches = vals_legend(vals=[0.1, 1, 2, 3, 4],
                                  var_max=5, colormap="cool")

        else:
            sys.exit("times:sigmx_vs_tmx : not implemented")

        ax.scatter(gpop[mask].tmx*u.s, gpop[mask].sigmx,
                   marker="o", color=colors, alpha=0.7, s=sizes)
        ax.scatter(gpop[mask].tmx*u.s, gpop[mask].sigmx,
                   marker=".", alpha=0.5, s=5, color="grey")
        ax.set_ylim(ymin=0.1)

        # Decorations
        ax.set_xlabel(xlabel + ax.get_xlabel() + "$)$")
        ax.set_ylabel(r"$ mean \ \sigma_{max}$")
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_title(tag)
    ax.grid(axis="y", ls=":")

    fig.legend(title=title, handles=patches,
               bbox_to_anchor=(1.1, 0.98))

    stamp(pop.tag[0], where="right")

    plt.tight_layout()


# ##---------------------------------------------------------------------------------------------
def sig_vs_time_errors(pop, mask, ax=None, xscale="log",
                       yscale="log", title=""):
    """
    Display the mean maximal significnace versus the mean time to reach it.

    Inludes error bars.

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    mask : Boolean sequence
        Data subselection.
    ax : Matplotlib axes, optional
        Current axis. The default is None.
    xscale : String, optional
        "Log" or "linear". The default is "log".
    yscale : String, optional
        "Log" or "linear". The default is "log".
    title : String, optional
        Descriptor. The default is "".

    Returns
    -------
    None.

    """
    ax = plt.gca() if ax is None else ax
    ax.errorbar(pop[mask].tmx/3600,
                pop[mask].sigmx,
                xerr=pop[mask].etmx/3600,
                yerr=pop[mask].esigmx,
                ls="", marker="o", ecolor="tab:blue", color="red", alpha=0.5,
                label=MyLabel(pop[mask].sigmx, stat=None))
    # ax.axhline(y=3,ls=":",label="$3\sigma$",color="lightgrey")
    # ax.axhline(y=5,ls=":",label="$5\sigma$",color="lightgrey")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel("Mean Time of max detection (h)")
    ax.set_ylabel("Mean max sigificance")
    ax.set_title(title)
    ax.legend()


# ##---------------------------------------------------------------------------------------------
def delay_correlation_both_sites(pop, var="tmx", sigmin=5, dtmax=26):
    """
    Find N and S detection for all both detected GRB.

    Get name of GRB detected on both sites and retrieve the independent N and S
    detections.

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    var: string
        Variable name to plot.
    sigmin: float
        Minimal mean maximal significance in the data. Default is 5.
    dtmax: flota
        Maximal time in hours on both axes.

    Returns
    -------
    None.

    """
    if var not in ["tmx", "t5s", "t3s", "t1"]:
        print(" Variable not allowed")
        return

    # ## Get North and South population seen onboth site - reindex the frames
    gbs = pop.g_s[(pop.g_s.name.isin(pop.g_b.name))].copy()
    gbn = pop.g_n[(pop.g_n.name.isin(pop.g_b.name))].copy()
    gb = pop.g_b.copy()

    gbs.reset_index(inplace=True)
    gbn.reset_index(inplace=True)
    gb.reset_index(inplace=True)

    # Apply cuts
    gbs = gbs[(gb.sigmx >= sigmin)]
    gbn = gbn[(gb.sigmx >= sigmin)]
    gb = gb[(gb.sigmx >= sigmin)]

    # Now that all frames have same index and size, do plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 7))
    color, size = col_size(gb.sigmx)

    ax.scatter(gbn[var]/3600, gbs[var]/3600, alpha=0.5, c=color, s=size,
               label=MyLabel(gbn[var], r"$\sigma_{max}$"))
    ax.set_title(r"Seen on both sites - $\sigma_{max} \geq "+str(sigmin)+"$")
    ax.set_xlabel("Mean detection time (h) - " + var + " North")
    ax.set_ylabel("Mean detection time (h) - " + var + " South")
    xmax = dtmax
    ymax = dtmax
    ax.set_xlim(xmax=xmax)
    ax.set_ylim(ymax=ymax)

    # ZOOM
    xmin_z = 0
    xmax_z = 1
    ymin_z = 0
    ymax_z = 10
    ax.axvspan(xmin_z - 1, xmax_z,
               ymin=ymin_z - 1, ymax=(ymax_z + 1)/ymax,
               color="grey", ls="-", alpha=0.5, lw=4, ec="black", fill=False)

    axx = inset_axes(ax, width="20%", height=4, loc="upper right")
    axx.scatter(gbn.t5s/3600, gbs.t5s/3600,  alpha=0.5, c=color, s=size,
                label=MyLabel(gbn.t5s, r"$\sigma_{max}$"))
    axx.set_xlim(xmin=-0.25, xmax=xmax_z)
    axx.set_ylim(ymin=-1, ymax=ymax_z)

    stamp(pop.tag[0], axis=fig, where="bottom")
    patches = vals_legend()
    fig.legend(title=r"$\sigma_{max}$", handles=patches,
               bbox_to_anchor=(1.01, 0.88))


# #############################################################################
if __name__ == "__main__":

    # nyears, files, tag = get_data(parpath=None, debug=True)
    parpath = "parameter_100k_ISM_alpha.yaml"
    nyears, files, tag = get_data(parpath=parpath, debug=False)

    pop = Pop(files, tag=tag, nyrs=nyears, fpeakmin=1)

    # ------------------------
    # Print some statistics
    # ------------------------
    # Print max significance reached after some time, above 3/5 sigma
    sigmx_after_delay(pop, sigmx_min=5, tmin=4*u.d)

    # ## ------------------------
    # ## Plots
    # ## ------------------------

    # Plot time for reaching 5 sigma or sigma max.

    # ax = time_to_sigma(pop.g_s, color=pop.col_s, tag="South")
    # ax = time_to_sigma(pop.g_n, color=pop.col_n, tag="North", ax=ax)

    # ax = time_to_sigma(pop.g_s, sigmin=10, color=pop.col_s, tag="South")
    # time_to_sigma(pop.g_n, sigmin=10, color=pop.col_n, tag="North", ax=ax)

    ax = time_to_sigma(pop.g_s, sigmx=False, color=pop.col_s, tag="South")
    time_to_sigma(pop.g_n, sigmx=False, color=pop.col_n, tag="North", ax=ax)


    # Time to get 5 sigma or sigmax in various conditions
    # time_to_detection(pop.g_tot, title="All")
    # time_to_detection(pop.g_n, title="North")
    # time_to_detection(pop.g_s, title="South")

    # # Scatter plot significance versus time to reach it
    # sigmx_versus_tmx(pop.g_n, tag="North")
    # sigmx_versus_tmx(pop.g_s, tag="South")
    # sigmx_versus_tmx(pop.g_tot, tag="Combined")

    # sigmx_versus_tmx(pop.g_tot, tag="Combined", axis="z")

    # # For GRBS detected on both sites, detection delays
    # delay_correlation_both_sites(pop, var="t3s", sigmin=5, dtmax=75)
    # delay_correlation_both_sites(pop, var="t5s", sigmin=5, dtmax=75)
    # delay_correlation_both_sites(pop, var="tmx", sigmin=5, dtmax=75)

    # # Scatter plot of sigma versus time with errors

    # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 8))

    # title = "All - " + r"$\sigma_{max} \geq 20$"
    # sig_vs_time_errors(pop.g_tot, (pop.g_tot.sigmx >= 20), yscale="log",
    #                    ax=ax1, title=title)
    # title = "All - " + r"$5 < \sigma_{max} < 20$"
    # sig_vs_time_errors(pop.g_tot, (pop.g_tot.sigmx >= 5)
    #                    & (pop.g_tot.sigmx < 20),
    #                    xscale="linear", yscale="linear", ax=ax2, title=title)
    # stamp(pop.tag[0], axis=fig)
    # plt.tight_layout()
