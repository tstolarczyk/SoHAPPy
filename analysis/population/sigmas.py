# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 17:43:24 2022.

Detection rates as a function of the mean maximal significance thresholds,
evolution of rates with this threshold, distribution of mean maximal
significances.

@author: Stolar
"""
import sys
import numpy as np

import seaborn as sns

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from niceprint import heading
from niceplot import MyLabel, stamp, single_legend, old_single_legend

from population import Pop
from pop_io import get_data

codefolder = "../../"
sys.path.append(codefolder)

__all__ = ["detection_level", "plot_sigmax_all", "plot_sigmax", "high_sigma",
           "plot_high_sigma_all", "sigmax_above_threshold",
           "rate_above_threshold"]


# ##---------------------------------------------------------------------------
def detection_level_vs_sigmx(pop, glist, bins=25, **kwargs):
    """
    Display the detection level for a list of subpopulation.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    gname : list of String
        Subpopulation names.
    nbin : integer, optional
        Histogram bin number. The default is 25.
    **kwargs : Dictionnary
        Extra arguments.

    Returns
    -------
    None.

    """
    _, ax = plt.subplots(nrows=1, ncols=len(glist),
                         figsize=(7*len(glist), 7),
                         sharey=True)
    if len(glist) == 1:
        ax = [ax]

    first = True
    for gname, axi in zip(glist, ax):

        subpop = pop.__dict__[gname]
        mask = subpop.sigmx > 0
        axi.hist2d(100*subpop[mask].d5s/pop.niter,
                   np.log10(subpop[mask].sigmx),
                   cmap="PuBu", bins=100, norm=mpl.colors.LogNorm())
        axi.grid(which="both", axis="both", ls=":", alpha=0.5)
        axi.axhline(np.log10(3),
                    color="tab:blue", ls="--", label=r"$3 \sigma$")
        axi.axhline(np.log10(5),
                    color="tab:green", ls="--", label=r"$5 \sigma$")
        axi.axhline(np.log10(20),
                    color="tab:orange", ls="--", label=r"$20 \sigma$")
        axi.axhline(np.log10(50),
                    color="red", ls="--", label=r"$50 \sigma$")
        axi.set_xlabel("Detection level (%)")
        axi.set_ylabel(r"$log_{10} \sigma_{max}$")
        axi.legend()
    plt.tight_layout()

    _, ax = plt.subplots(nrows=1, ncols=len(glist),
                         figsize=(7*len(glist), 7),
                         sharey=True)
    if len(glist) == 1:
        ax = [ax]

    for gname, axi in zip(glist, ax):
        subpop = pop.__dict__[gname]
        first = True
        for cut, col in zip([3, 5, 7, 10],
                            ["tab:blue", "tab:green", "tab:orange", "red"]):
            if first:
                bins = 100
                first = False
            mask = subpop.sigmx > cut
            _, bins, _ = axi.hist(100*subpop[mask].d5s/pop.niter, bins=bins,
                                  edgecolor=col, facecolor=None,
                                  histtype="step", label="> "+str(cut))
            axi.grid(which="both", axis="both", ls=":", alpha=0.5)
            axi.set_xlabel("Detection level (%)")
            axi.set_yscale("log")
            axi.legend()

    plt.tight_layout()


# ##---------------------------------------------------------------------------
def detection_level(pop, glist, bins=25, **kwargs):
    """
    Display the detection level for a list of subpopulation.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    gname : list of String
        Subpopulation names.
    nbin : integer, optional
        Histogram bin number. The default is 25.
    **kwargs : Dictionnary
        Extra arguments.

    Returns
    -------
    None.

    """
    _, ax = plt.subplots(nrows=1, ncols=len(glist),
                         figsize=(7*len(glist), 7),
                         sharey=True)
    first = True
    for gname, axi in zip(glist, ax):

        subpop = pop.__dict__[gname]

        for sigma, var, col in zip([3, 5],
                                   [subpop.d3s, subpop.d5s],
                                   [Pop.col_3, Pop.col_5]):

            var = 100*var/pop.niter
            label = str(sigma)+r"$\sigma$"

            axi.hist(var, bins, range=[0, 100],
                     label=MyLabel(var, label), color=col, alpha=0.5, **kwargs)

            axi.set_xlabel("Confidence level (%)")
            axi.set_yscale("log")
            axi.set_title(gname)
            axi.axvline(x=pop.eff_lvl,
                        color="red", ls=":", label="min. C.L. (%)")
            axi.grid("both")

            if first:
                axi.set_ylabel("Event counts")
                first = False

        old_single_legend(axi)

    plt.tight_layout()


# ##---------------------------------------------------------------------------
def plot_sigmax_all(pop, nbin=25, logx=True):
    """
    Plot the mean maximal significance distributions.

    The plots are done for North and North only,
    South and South only, and GRBs detected on both site. Check that detections
    with both sites give higher significances.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    nbin : TYPE, optional
        Number of bins in histograms. The default is 25.

    Returns
    -------
    None.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3,
                                        figsize=(18, 6), sharey=True)

    # North and North only
    ax, bins = plot_sigmax(pop.g_n, tag="North",
                           ax=ax1, bins=nbin, logx=logx,
                           color=Pop.col_n, alpha=0.5)
    plot_sigmax(pop.g_n0, tag="North only",
                ax=ax1, logx=logx, bins=bins,
                color=Pop.col_n, alpha=0.5)

    # South and South only
    ax, bins = plot_sigmax(pop.g_s, tag="South",
                           ax=ax2, bins=nbin, logx=logx,
                           color=Pop.col_s, alpha=0.5)
    plot_sigmax(pop.g_s0, tag="South only",
                ax=ax, logx=logx, bins=bins,
                color=Pop.col_s, alpha=0.5)

    # Both, all
    ax, bins = plot_sigmax(pop.g_b, tag="Both",
                           logx=logx, color=Pop.col_b,
                           ax=ax3, bins=nbin, alpha=0.5)
    plot_sigmax(pop.g_tot, tag="N+S",
                logx=logx, facecolor="none", edgecolor="black",
                ax=ax3, bins=bins, alpha=1,
                weight=len(pop.g_b)/len(pop.g_tot))

    stamp(pop.tag[0], axis=fig, where="bottom")
    plt.tight_layout()


# ##---------------------------------------------------------------------------
def plot_sigmax(g, ax=None, logx=True, bins=25,
                tag="", weight=1, xmax=None, **kwargs):
    """
    Plot mean maximal significance for a given population.

    Parameters
    ----------
    g : Pandas table
        Current subpopulation.
    ax : Matplotlib axes, optional
        Current axis. The default is None.
    logx : Boolean, optional
        Abscissa in log scale if True. The default is True.
    bins : integer, optional
        number of bisn in histrogram. The default is 25.
    tag : string, optional
        Descriptor of the data. The default is "".
    weight : float, optional
        Histrogram reweighting. The default is 1.
    xmax : Float, optional
        Abscissa maximal value. The default is None.
    **kwargs : Dictionnary
        Extra matplotlib arguments.

    Returns
    -------
    ax : Matplotlib axes
        Current axis.
    bins : TYPE
        DESCRIPTION.

    """
    if ax is None:
        _, ax = plt.subplots()

    if xmax is not None:
        mask = g.sigmx < xmax
    else:
        mask = (np.ones(len(g.sigmx), dtype=bool)) & (g.sigmx != -1)

    if logx:
        mask = (g.sigmx > 0) & (mask)
        x = np.log10(g[mask].sigmx)

        ax.set_title("Positive values only")
        ax.set_xlabel(r"$log_{10}(\sigma_{max})$")
        ax.axvline(x=np.log10(3), color="tab:orange", ls=":",
                   label=r"$3\sigma$")
        ax.axvline(x=np.log10(5), color="tab:green", ls=":",
                   label=r"$5\sigma$")
    else:
        x = g[mask].sigmx
        ax.set_xlabel(r"$\sigma_{max}$")
        ax.axvline(x=3., color="tab:orange", ls=":", label=r"$3\sigma$")
        ax.axvline(x=5., color="tab:green", ls=":", label=r"$5\sigma$")

    _, bins, _ = ax.hist(x, bins=bins, label=MyLabel(x, tag, stat="med"),
                         weights=np.ones(len(x))*weight, **kwargs)

    old_single_legend(ax)
    return (ax, bins)


# ##---------------------------------------------------------------------------
def high_sigma(gpop, ax=None, inset=True, sigmin=5,
               sig2show=100, tag="", weight=1, **kwargs):
    """
    Display the high values of the mean maximal significance.

    Tag the highest value with the event number.

    Parameters
    ----------
    gpop : Pandas table
        Current subpopulation.
    ax : Matplotlib axis, optional
        Current axis. The default is None.
    inset : Boolean, optional
        If True, add a log-scale display in an inset box. The default is True.
    sigmin : float, optional
        Minimal significance to consider. The default is 5.
    sig2show : float, optional
        Minimal significance to be tagged. They willnot be shown if negative.
        The default is 100.
    tag : String, optional
        Data descriptor. The default is "".
    **kwargs : Dictionnary
        Extra parameters.

    Returns
    -------
    ax : Matplotlib axes
        Current axis.
    bins : sequence
        Histogram binning.

    """
    var = gpop.sigmx[gpop.sigmx >= sigmin]
    weights = np.ones(len(var))*weight

    _, bins, _ = ax.hist(var, label=MyLabel(var*weights, tag, stat=None),
                         weights=weights, **kwargs)
    ax.set_xlabel(rf"$\sigma_{{max}} \geq {sigmin:}$")
    ax.set_ylabel("Event rate")

    # Inset in log
    if inset:
        axx = inset_axes(ax, width="75%", height=1.2, loc="upper right")
        _, _, _ = axx.hist(np.log10(gpop[gpop.sigmx > 0].sigmx),
                           bins=25, color="grey", edgecolor="black", alpha=0.5)
        axx.axvline(np.log10(3), ls=":", color="red", lw="2",
                    label=r"$3\sigma$")
        axx.axvline(np.log10(5), color="red", lw="2", label=r"$5\sigma$")
        axx.set_xlabel(r"Log $\sigma_{max}$")
        axx.legend()

        ax.legend(loc="upper left")
    else:
        ax.legend()

    # Outliers
    if sig2show > 0:
        outlier = zip(gpop[gpop.sigmx > sig2show].name,
                      gpop[gpop.sigmx > sig2show].sigmx)
        i = 1
        for g in outlier:
            y = (0 + (1.5*i + 0.1))*weight
            ax.text(x=g[1], y=y, s=g[0][5:], rotation=45)
            i += 1

    return ax, bins


# ##---------------------------------------------------------------------------
def plot_high_sigma_all(pop, nbin=100, sigmin=5, inset=False):
    """
    Plot the high value tail of the mean maximal sigificance values.

    Tag the data with the event number

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    nbin: integer
        Histrogram number of bins. Default is 100
    sigmin: float
        Minimal significance values in the plots. The default is 5
    inset : boolean, optional
        If True, display the same plot in log-scale in an inset box.
        The default is False.

    Returns
    -------
    None.

    """
    weight = 1/pop.nyears

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(12, 9))

    # North and North only
    ax0, bins = high_sigma(pop.g_n, sigmin=sigmin,
                           bins=nbin, color=Pop.col_n, ax=ax[0],
                           tag="N", alpha=0.5, weight=weight, inset=inset)
    high_sigma(pop.g_n0,  sigmin=sigmin,
               bins=bins, color=Pop.col_n, ax=ax0, tag="N only",
               weight=weight, inset=inset)
    ax0.set_xlabel(None)
    ax0.set_title(r" Mean $\sigma_{max}$ for "
                  + str(pop.nyears)
                  + " years of observation")
    ax0.grid("both", ls=":")

    # South and South only
    ax0, bins = high_sigma(pop.g_s, sigmin=sigmin,
                           bins=nbin, color=Pop.col_s, ax=ax[1], tag="S",
                           alpha=0.5, weight=weight, inset=inset)
    high_sigma(pop.g_s0, sigmin=sigmin,
               bins=bins, color=Pop.col_s, ax=ax0, tag="S only",
               weight=weight, inset=inset)
    ax0.set_xlabel(None)

    # Both sites
    ax0, bins = high_sigma(pop.g_b, sigmin=sigmin,
                           bins=nbin, color=Pop.col_b, ax=ax[2],
                           tag="Both", alpha=0.5, weight=weight, inset=inset)
    ax0.grid("both", ls=":")

    stamp(pop.tag[0], axis=fig, where="bottom")
    plt.tight_layout()


# ##----------------------------------------------------------------------------
def sigmax_above_threshold(pop, sigmin=0):
    """
    Print event number with mean maximal significance above a given threshold.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    sigmin : float, optional
        Signiifcnace threshold. The default is 0.

    Returns
    -------
    None.

    """
    heading(f"Source ids with mean maximal sigma above {sigmin:5.1f}")

    poplist = [pop.g_n, pop.g_n0, pop.g_s, pop.g_s0, pop.g_b]
    taglist = ["North", "N only", "South", "S only", "Both"]

    # Mean sigmax value of events above 5 sigma at 90%CL and sigmin
    for gsub, tag in zip(poplist, taglist):
        print("---- ", tag)
        gsub = gsub[gsub.d5s >= pop.eff_lvl]
        for sig, name in zip(gsub.sigmx, gsub["name"]):
            if sig >= sigmin:
                print(" - ", name, " :", sig)


# ##---------------------------------------------------------------------------
def rate_above_threshold(pop, tobs=1, plot=True):
    """
    Print detection rates for various mean maximal significance thresholds.

    Parameters
    ----------
    pop : TYPE
        DESCRIPTION.
    duration : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    poplist = [pop.g_n, pop.g_s, pop.g_n0, pop.g_s0, pop.g_b]
    taglist = ["North", "South", "N only", "S only", "Both"]
    siglist = [1, 3, 5, 10, 20, 50, 100, 200]

    # Compute values for each population, for a duration.
    vals = {}
    vtot = {}
    for gpop, tag in zip(poplist, taglist):
        # mask = (gpop.z >= 0.3) & (gpop.z < 0.7)
        mask = (gpop.z >= 0.7) & (gpop.z < 5)
        # mask = gpop.z > 0
        vals[tag] = [len(gpop[(gpop.sigmx >= sig) & mask])*tobs/pop.nyears
                     for sig in siglist]

    # masktot = (pop.g_tot.z >= 0.3) & (pop.g_tot.z < 0.7)
    masktot = (pop.g_tot.z >= 0.7) & (pop.g_tot.z < 5)
    vtot = [len(pop.g_tot[(pop.g_tot.sigmx >= sig) & masktot])*tobs/pop.nyears
            for sig in siglist]

    # Create and display table with numerical results
    heading(f"Detection rate with mean sigmax above sig for {tobs:} years")

    print(f"{'Sig min':>10s}", end="")
    for tag in taglist:
        print(f"{tag:>8s}", end="")
    print(f"{'Tot':>8s}")
    print(60*"-")

    for isig, sig in enumerate(siglist):

        # One line per min sigma
        print(f"{sig:>10.1f}", end="")
        for site in vals.keys():
            print(f"{vals[site][isig]:>8.2f}", end="")

        print(f"{vtot[isig]:>8.2f}")  # Total statistics

    print(60*"-")

    # Display results in a plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    if plot:

        for gpop, tag in zip(poplist, taglist):
            plt.plot(siglist, vals[tag], label=tag, marker="o")
        plt.plot(siglist, vtot, label="Total", marker="o")

        width = 5
        ax.axvline(50, color="grey", alpha=0.5, ls="--", lw=width)
        ax.text(50*1.05, 0.05*tobs, r"$50\sigma$", size=12)

        ax.axvline(20, color="grey", alpha=0.5, ls="--", lw=width)
        ax.text(20*1.05, 0.05*tobs, r"$20\sigma$", size=12)

        ax.axvline(5, color="grey", alpha=0.5, ls=":", lw=width)
        ax.text(5*1.05, 0.05*tobs, r"$5\sigma$", size=12)

        # Indicate rates
        if tobs == 10:
            ax.text(min(siglist), 1.1, "One per decade", alpha=0.7, size=12)
            ax.text(min(siglist), 11., "One per year", alpha=0.7, size=12)

        ax.set_xlabel(r"$\sigma_{max}$")
        ax.set_ylabel(r"Rate above $\sigma_{max}$ for "
                      + str(tobs) + "  yr of operation")
        ax.set_yscale("log")
        ax.set_xscale("log")

        from niceplot import old_single_legend
        old_single_legend(ax, bbox_to_anchor=[1.0, 1.0])
        ax.grid(which="minor", ls=":")
        ax.grid(which="major", ls="-")

        stamp(pop.tag[0], axis=fig, where="bottom")
        plt.tight_layout()


# ##---------------------------------------------------------------------------
def sigma_correlation_both_sites(pop, sigmin=1):
    """Plot sigma performance South versus North."""
    # ## Get North and South population seen onboth site - reindex the frames
    gbs = pop.g_s[(pop.g_s.name.isin(pop.g_b.name))]
    gbn = pop.g_n[(pop.g_n.name.isin(pop.g_b.name))]
    gb = pop.g_b

    gbs.reset_index(inplace=True)
    gbn.reset_index(inplace=True)
    gb.reset_index(inplace=True)

    mask = (gbs.sigmx >= sigmin) & (gbn.sigmx >= sigmin)

    # Sigma South versus North
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

    ax1.scatter(gbs[mask].sigmx, gbn[mask].sigmx, marker=".",
                label=MyLabel(gb[mask].sigmx, "Combined"))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"$\sigma_{max} \ South$")
    ax1.set_ylabel(r"$\sigma_{max} \ North$")
    # ax1.grid(axis="both", which="minor", ls="-", color="lightgrey")
    ax1.grid(axis="both", which="major", ls="-")
    ax1.set_xlim([1, 1000])
    ax1.set_ylim([1, 1000])
    ax1.plot([0.001, 1e4], [0.001, 1e4], color="red", ls="--")
    ax1.legend()

    # Difference hsitorgram
    fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(15, 8))
    ax2.hist(gbs[mask].sigmx - gbn[mask].sigmx, bins=100,
             label=MyLabel(gbs[mask].sigmx - gbn[mask].sigmx, "South - North"))
    ax2.set_yscale("log")
    ax2.set_xlim([-500, 500])
    ax2.set_xlabel(r"$\Delta\sigma_{max} \ (South-North)$")
    ax2.legend()

    # Diference versus redhsift
    fig, ax3 = plt.subplots(nrows=1, ncols=1, figsize=(15, 8))

    ymin = -5
    ymax = 5
    xmin = 0
    xmax = 4
    nbin = 10

    ax3.scatter(gbs[mask].z, (gbs[mask].sigmx-gbn[mask].sigmx)/gbs[mask].sigmx,
                marker=".",
                label=MyLabel(gbs[mask].sigmx / gbn[mask].sigmx,
                              "South / North"))
    ax3.set_ylim([ymin, ymax])
    ax3.set_xlim([xmin, xmax])
    ax3.set_ylabel(r"$ratio \sigma_{max} \ (South/North)$")
    ax3.set_xlabel("Redshift")
    axx3 = ax3.twinx()
    HH, xe, ye = np.histogram2d(gbs[mask].z,
                                (gbs[mask].sigmx - gbn[mask].sigmx)
                                / gbs[mask].sigmx,
                                bins=nbin,
                                range=[[xmin, xmax], [ymin, ymax]])

    # Do it smarter !
    aver = []
    for item, x in enumerate(0.5*(xe[1:]+xe[:-1])):  # along z
        # print("z =",x)
        # print(" Averaged mean = ",
        #       np.average(0.5*(ye[1:]+ye[:-1]), weights= HH[item]))
        aver.append(np.average(0.5*(ye[1:]+ye[:-1]), weights=HH[item]))
    axx3.plot(0.5*(xe[1:] + xe[:-1]), aver, color="red")
    axx3.grid(axis="y", which="major", ls=":")
    axx3.axhline(np.mean(aver), ls="--", color="red")
    # ax3.legend()


###############################################################################
if __name__ == "__main__":

    # Bigger texts and labels
    sns.set_context("talk")  # poster, talk, notebook, paper

    nyears, files, tag = get_data(parpath=None,debug=True)
    # parpath = "parameter_100k_ISM_alpha.yaml"
    # parpath = "parameter_100k_ISM_omega.yaml"
    # nyears, files, tag = get_data(parpath=parpath, debug=False)

    pop = Pop(files, tag=tag, nyrs=nyears, fpeakmin=1)

    sigma_correlation_both_sites(pop)

    # ## ------------------
    # ## Detection rate versus minimal significance normalised to duration
    # ## ------------------
    rate_above_threshold(pop, tobs=10)
    rate_above_threshold(pop, tobs=1)

    # # Mean maximal significance above threshold
    # sigmax_above_threshold(pop, sigmin=100)

    # # ## ------------------
    # # ## Plots
    # # ## ------------------

    # # Detection level distributions
    # detection_level(pop, ["g_tot", "g_n", "g_s"])
    # detection_level(pop, ["g_n", "g_s"])
    # detection_level_vs_sigmx(pop, ["g_tot"])
    # # Mean max significance distribution - compare N/Nonly etc.
    # plot_sigmax_all(pop)

    # # Highest significance values
    # plot_high_sigma_all(pop, inset=False)
