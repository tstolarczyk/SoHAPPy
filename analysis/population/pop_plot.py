# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 11:58:24 2021.

Show distributions and population coverage as functions of `Eiso`, `Epeak`
and the redshift `z`.

@author: Stolar
"""
import sys

import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

import seaborn as sns

from niceplot import MyLabel, single_legend, stamp, projected_scatter, \
                     col_size, vals_legend, lower_limit
from historical import plot_historical, historical
from population import Pop
from pop_io import get_data

# Bigger texts and labels
sns.set_context("notebook")  # poster, talk, notebook, paper

codefolder = "../../"
sys.path.append(codefolder)

__all__ = ["distri", "coverage"]


# ##----------------------------------------------------------------------------
def plot_generated(pop):
    """
    Plot charcteristics of the generated sources.

    Parameters
    ----------
    pop : pandas table
        The population in hand.

    Returns
    -------
    None.

    """
    # Use the North classified population to address the GRBs once
    gpop = pop.g_n

    # Define the plot framework

    fig, ax, axh, axv = projected_scatter()

    # ## -------------
    # ## Projections
    # ## -------------
    nbin = 100
    axh.hist(gpop.z, bins=nbin)
    axh.set_yscale("log")
    axh.grid("both", ls=":", color="black")

    axv.hist(np.log10(gpop.Eiso), bins=nbin,  orientation="horizontal")
    axv.set_xscale("log")
    axv.grid("both", ls=":", color="black")

    # ## -------------
    # ## Population coverage
    # ## -------------

    # Scatter plot is nnot very effective for large populations
    # ax.scatter(gpop.z, np.log10(gpop.Eiso), marker=".", color="grey")

    # 2D histrogram with bin height in log scale
    _, _, _, img = ax.hist2d(gpop.z, np.log10(gpop.Eiso),
                             norm=mpl.colors.LogNorm(), cmap="PuBu", bins=100)

    # Put colorbar at the right of the projected axis
    # cbar = fig.colorbar(img, ax=axv)
    # cbar.set_label('Counts')

    # Superimpose contrours
    levels = [50, 100, 300, 1000]
    HH, xe, ye = np.histogram2d(gpop.z, np.log10(gpop.Eiso), bins=25)
    grid = HH.transpose()
    midpoints = (xe[1:] + xe[:-1])/2, (ye[1:] + ye[:-1])/2
    cntr = ax.contour(*midpoints, grid, colors="red", levels=levels)
    # clbl = cntr.clabel()  # Labels

    # # Remove contours with "0" label as it is chaotic
    # cntr.set_paths(cntr.get_paths()[1:])  # Remove last (chaotic) contour

    # for ic, c in enumerate(clbl):
    #     if c._text == "0":
    #         c.set_text("")

    # Decorate
    ax.grid("both", ls=":", color="black")
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$log_{10} (E_{iso} (erg))$")

    return ax, axh, axv


# ##----------------------------------------------------------------------------
def distri(pop, grbs,
           var="z", var_range=[0, 4], var_log=False, var_name="", tag="",
           weight=1, ax=None, nbin=20, color="tab:blue", color2="red",
           **kwargs):
    """
    Plot the given variable wrt to the reference population and the ratio.

    Parameters
    ----------
    pop : Pandas table
        A population.
    grbs : pandas table
        A subselection of the population.
    var : String, optional
        A column name in the Pandas table. The default is "z".
    range : List, optional
        The variable minimal and maximal values allowed in the plot.
        The default is [0, 4].
    var_log : Boolean, optional
        Plot logarithm10 of the variable if true. The default is False.
    var_name : String, optional
        Explicit name of the plotted variable. The default is "".
    weight: float
        A weight applied to th epopuelation data in case it would sum-up
        several ones. The default is one.
    tag: String, optional
        Describes the data in the plot. The default is "".
    ax : matplotlib axes, optional
        Current matplotlib axes. The default is None.
    nbin : integer, optional
        Number of bins in the histrogram. The default is 20.
    color : string, optional
        Selected subpopulation color. The default is "tab:blue".
    color2 : String, optional
        Ratio curve color. The default is "red".
    **kwargs : dictionnary
        matplotlib supplementary parameters.

    Returns
    -------
    ax : matplotlib axis
        Population distribution axis.
    axx : matplotlib axis
        Dsitribtution ratio secondary axis.

    """
    ax = plt.gca() if ax is None else ax

    # If the variable name is not explcitely given, use the column name
    if var_name == "":
        var_name = var

    # Plot the reference GRB population
    mask = (var_range[0] <= pop.ref[var]) & (pop.ref[var] <= var_range[1])

    if var_log:
        x = [np.log10(v) if v > 0 else 0 for v in pop.ref[mask][var]]
    else:
        x = pop.ref[mask][var]

    nref, bins, _ = ax.hist(x, bins=nbin,
                            facecolor="none", edgecolor="black")

    # Plot the requested population - can be weighted if the sum of several
    mask = (var_range[0] <= grbs[var]) & (grbs[var] <= var_range[1])

    if var_log:
        x = [np.log10(v) if v > 0 else 0 for v in grbs[mask][var]]
    else:
        x = grbs[mask][var]

    n, bins, _ = ax.hist(x,
                         bins=bins, color=color,
                         label=MyLabel(x, label=tag),
                         weights=weight*np.ones(len(x)),
                         **kwargs)
    ax.legend()
    ax.set_xlabel(var_name)
    ax.set_yscale("log")

    # Plot the ratio
    ratio = [(x/xtot if xtot != 0 else 0) for x, xtot in zip(n, nref)]
    axx = ax.twinx()
    axx.plot(bins[:-1] + 0.5*(bins[1:] - bins[:-1]), ratio,
             alpha=1, color=color2, label="Ratio")

    axx.grid(ls=":")
    axx.legend(loc="lower right")

    return ax, axx


# ##---------------------------------------------------------------------------
def coverage(varx, vary,
             ref=None, pop=None, mask=None,
             xrange=(None, None), yrange=(None, None),
             lblx="", lbly="", title="dummy",
             xscale="log", yscale="log",
             nbin=25,
             alpha_ref=0.5, color_ref="grey", alpha=0.5,
             projs=True, tag=""):
    """
    Plot population coverage in two variable space.

    Includes a plot for the reference pop.

    Parameters
    ----------
    varx : String
        Abscissa variable.
    vary : String
        ordinate variable.
    ref : pandas dataframe, optional
        Reference population. The default is None.
    pop : pandas dataframe, optional
        Studied population. The default is None.
    mask : array of boolean, optional
        mask applied to the studied population. The default is None.
    xrange : list of two elements, optional
        Abscissa limits. The default is [None,None].
    yrange : List of two elemenst, optional
        Ordinate limits. The default is [None, None].
    lblx : string, optional
        Abscissa label specific tag. The default is "".
    lbly : String, optional
        ordinate label specific tag. The default is "".
    xscale : String, optional
        Logarithmic or linear abscissa scale. The default is "log".
    yscale : String, optional
        logarithmic or linear ordinate scale. The default is "log".
    title : String, optional
        Plot title. The default is "dummy".
    nbin : integer, optional
        Nimber of bins in projected histrograms. The default is 25.
    alpha : float, optional
        Studied population transparency. The default is 0.5.
    alpha_ref : float, optional
        reference population transparency. The default is 0.5.

    Returns
    -------
    None.

    """
    if projs:
        fig, ax, axh, axv = projected_scatter()
    else:
        fig, ax = plt.subplots(figsize=(10, 6))

    # Create dummy mask if undefined
    if mask is None:
        mask = (pop[varx] == pop[varx])

    # change TUple into list ofr assignments.
    # Create dummy range limits if undefined
    xrange = list(xrange)
    yrange = list(yrange)

    if xrange[0] is None:
        xrange[0] = np.min(ref[varx])
    if xrange[1] is None:
        xrange[1] = np.max(ref[varx])
    if yrange[0] is None:
        yrange[0] = np.min(ref[vary])
    if yrange[1] is None:
        yrange[1] = np.max(ref[vary])

    # ## --------------------
    # ## Central scatter plot
    # ## --------------------
    # reference population

    ax.scatter(np.log10(ref[varx]) if xscale == "log" else ref[varx],
               np.log10(ref[vary]) if yscale == "log" else ref[vary],
               marker=".", color=color_ref, s=10, alpha=alpha_ref, label="All")

    xdet = np.log10(pop[mask][varx]) if xscale == "log" else pop[mask][varx]
    ydet = np.log10(pop[mask][vary]) if yscale == "log" else pop[mask][vary]

    colors, sizes = col_size(pop[mask].sigmx)
    ax.scatter(xdet, ydet, marker="o", s=sizes, c=colors, alpha=alpha)

    plot_historical(ax, historical(), obs=["H.E.S.S"])

    # Decoration
    if xscale == "log":
        ax.set_xlim(xmin=np.log10(xrange[0]), xmax=np.log10(xrange[1]))
        ax.set_xlabel(r"$log_{10} \ $" + lblx)
    else:
        ax.set_xlim(xmin=xrange[0], xmax=xrange[1])
        ax.set_xlabel(lblx)

    if yscale == "log":
        ax.set_ylim(ymin=np.log10(yrange[0]), ymax=np.log10(yrange[1]))
        ax.set_ylabel(r"$log_{10} \ $" + lbly)
    else:
        ax.set_ylim(yrange)
        ax.set_ylabel(lbly)

    ax.grid(axis="both", ls="--")

    patches = vals_legend(alpha=alpha)

    if projs:
        fig.legend(title=r"$\sigma_{max}$", handles=patches,
                   bbox_to_anchor=[1.05, 1.01], ncol=2)
    else:
        fig.legend(title=r"$\sigma_{max}$ - " + title, handles=patches,
                   bbox_to_anchor=[0.9, 0.85], ncol=2)

    if projs:
        # ## --------------------
        # ## horizontal data
        # ## --------------------
        hist_mask = (ref[varx] >= xrange[0]) & (ref[varx] <= xrange[1])

        xref = np.log10(ref[hist_mask][varx]) if xscale == "log" \
            else ref[hist_mask][varx]
        _, bins, _ = axh.hist(xref, bins=nbin, facecolor="none",
                              edgecolor="black")

        x = np.log10(pop[mask][varx]) if xscale == "log" else pop[mask][varx]
        axh.hist(x, bins=bins, color="purple", alpha=0.3)

        if xscale == "log":
            axh.set_xlim(xmin=np.log10(xrange[0]),
                         xmax=np.log10(xrange[1]))
        else:
            axh.set_xlim(xmin=xrange[0],
                         xmax=xrange[1])

        axh.set_title(tag+" - " + title)
        axh.set_yscale("log")

        # ## --------------------
        # ## Vertical data
        # ## --------------------
        hmask = (ref[vary] >= yrange[0]) & (ref[vary] <= yrange[1])

        yref = np.log10(ref[hmask][vary]) if yscale == "log" \
            else ref[hmask][vary]
        _, bins, _ = axv.hist(yref, bins=nbin,
                              facecolor="none", edgecolor="black",
                              orientation="horizontal",
                              label=MyLabel(yref, label="All"))

        y = np.log10(pop[mask][vary]) if yscale == "log" else pop[mask][vary]
        axv.hist(y, bins=bins,
                 color="purple", alpha=0.3, orientation="horizontal",
                 label=MyLabel(y, label="Detected"))
        axv.set_xscale(yscale)

    return ax  # Central plot


# ##---------------------------------------------------------------------------
def pop_lower_limits(varx, vary, ref_file=None,
                     fit_deg=3, bins=25, plot=True):
    """
    Compute lower limit in coverage for a reference population.

    Parameters
    ----------
    varx : String
        Abscissa
    vary : String
        Ordinate
    ref_file : pathlib Path, optional
        A parameter file with the reference population. The default is None.
    fit_deg : integer
        the order of the polynomial fit.
    plot : boolean, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """
    if ref_file is None:
        ref_file = "../../data/samples/max_detection_parameter.yaml"

    print("-> Prepare lower limit curve from reference: ", ref_file)
    print("   This can be long!")

    nyears, files, tag = get_data(parpath=ref_file, debug=False)
    ref_pop = Pop(files, tag=tag, nyrs=nyears)

    fig, ax = plt.subplots()

    mask = ref_pop.g_tot["d3s"] > 0
    x_det = ref_pop.g_tot[mask][varx]
    y_det = ref_pop.g_tot[mask][vary]
    logx = False
    logy = True

    if logx:
        x_det = np.log10(x_det)
    if logy:
        y_det = np.log10(y_det)

    pfit = lower_limit(x_det, y_det, ax=ax, fitdeg=fit_deg, bins=bins)

    return pfit


# #############################################################################
if __name__ == "__main__":

    # Bigger texts and labels
    sns.set_context("talk")  # poster, talk, notebook, paper

    # ---------------------------------------------
    # Standard stricmoonveto demo dataset (1000 GRBs)
    # ---------------------------------------------
    # parpath = None

    # ---------------------------------------------
    # Max detection new sample(>10000 GRBs)
    # ---------------------------------------------
    # parpath = "max_detection_parameter_100k.yaml"

    # parpath = "parameter.yaml"  # Select a particular data set
    parpath = "parameter_100k_ISM_alpha.yaml"  # Select a particular data set
    # parpath = "parameter_100k_ISM_omega.yaml"  # Select a particular data set

    # Read population
    nyears, files, tag = get_data(parpath=parpath, debug=False)

    # Read population - and compute population seen in N or S (x2 in number)
    pop = Pop(files, tag=tag, nyrs=nyears, fpeakmin=1)

    # Population seen in North or South (x2 in number)
    popNS = pop.grb[(pop.grb.loca == "North") | (pop.grb.loca == "South")]

    # ## ----------------------------------
    # ## Plot one-dim variable coverage
    # ## ----------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(21, 15))

    # Loop over the reference population and superimpose selected data
    for icol, grbs, weight, tag in zip(range(2),
                                       [popNS, pop.g_tot],
                                       [0.5, 1],
                                       [r"North & South 5$\sigma$",
                                        r"Combined 5$\sigma$"]):
        print(icol)

        # Selection
        grbs = grbs[grbs.d5s >= pop.eff_lvl]

        distri(pop, grbs, var="z", tag=tag,
               var_range=[0, 4],
               ax=ax[0][icol],
               var_name="Redshift (z)",
               weight=weight)

        distri(pop, grbs, var="Eiso", tag=tag,
               var_range=[5.5e50, 5.5e55], var_log=True,
               ax=ax[1][icol],
               var_name=r"$log_{}10 E_{iso}$",
               weight=weight)

        distri(pop, grbs, var="Epeak", tag=tag,
               var_range=[10, 1e5], var_log=True,
               ax=ax[2][icol],
               weight=weight)

        plt.tight_layout()

    # ## ------------------------------------------
    # ## Plot Eiso versus z coverage in 2-dim space
    # ## -------------------------------------------

    poplist = [pop.g_n, pop.g_s, pop.g_tot]
    taglist = ["North", "South", "Combined"]
    masklist = [pop.g_n.d5s >= pop.eff_lvl,
                pop.g_s.d5s >= pop.eff_lvl,
                pop.g_tot.d5s >= pop.eff_lvl]
    # masklist = [pop.g_n.d5s   >= 0,
    #             pop.g_s.d5s   >= 0,
    #             pop.g_tot.d5s >= 0]
    if pop.niter == 1:
        title = r"Detected at $5\sigma$ - No count fluctuation"
    else:
        title = rf"Detected at $5\sigma - ({pop.eff_lvl:2.0f}\% \ C.L.)$"

    varx = "z"
    vary = "Eiso"
    lblx = "Redshift $z$"
    lbly = r"$E_{iso} \ (erg)$"
    xrange = [0, 5]
    yrange = [5e50, 1e56]
    # xrange = [None, None]
    # yrange = [None, None]

    # Prepare lower limits from the max. detection population
    pfit = None
    # pfit = pop_lower_limits(varx, vary, ref_file="parameter_100k_ISM-max.yaml",
    #                         fit_deg=4, bins=20)

    # # Plot coverage
    for grbs, mask, tag in zip(poplist, masklist, taglist):

        # With projections
        ax = coverage(varx, vary, ref=pop.ref, pop=grbs,
                      xrange=xrange, yrange=yrange,
                      lblx=lblx, lbly=lbly,
                      xscale="linear",
                      mask=mask,
                      title=title,
                      alpha=0.5,
                      alpha_ref=0.1,
                      tag=tag)

        # Plot lower limits from max detection
        if pfit is not None:
            xlim = ax.get_xlim()
            xsample = np.linspace(xlim[0], xlim[1], 20)
            ax.plot(xsample, np.polyval(pfit, xsample),
                    ls="--", lw=4, alpha=0.5, color="red")

        # Identify the data set
        stamp(pop.tag[0], axis=fig, x=1.04, y=0.5, rotation=270)

    #     # # Without projections
    #     # ax = coverage(varx,vary,ref=pop.ref, pop=grbs,
    #     #           xrange=xrange,yrange=yrange,
    #     #           lblx = lblx, lbly = lbly,
    #     #           xscale="linear",
    #     #           mask=mask,
    #     #           title=tag,
    #     #           projs=False)

        stamp(pop.tag[0], axis=plt.gcf(), x=1.04, y=0.5, rotation=270)

    # ## ----------------------------------------------
    # ## Plot Eiso versus Epeak coverage in 2-dim space
    # ## ----------------------------------------------
    # Eiso versus Epeak
    varx = "Fpeak"
    vary = "Eiso"
    lblx = r"$F_{peak} \ (keV)$"
    lbly = r"$E_{iso} \ (erg)$"
    xrange = [0.01, 1e5]
    yrange = [5e50, 1e56]

    for grbs, mask, tag in zip(poplist, masklist, taglist):

        coverage(varx, vary, ref=pop.ref, pop=grbs,
                 xrange=xrange, yrange=yrange,
                 lblx=lblx, lbly=lbly,
                 mask=mask,
                 title=title,
                 tag=tag)
        stamp(pop.tag[0], axis=fig, x=1.04, y=0.5, rotation=270)
