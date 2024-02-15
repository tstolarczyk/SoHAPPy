# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 11:58:24 2021

Show distributions and population coverage as functions of `Eiso`, `Epeak`
and the redshift `z`.

@author: Stolar
"""
import sys

import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

from niceplot import MyLabel, single_legend, stamp, projected_scatter, col_size, vals_legend
from historical import plot_historical, historical
from population import Pop
from pop_io import get_data

# Bigger texts and labels
sns.set_context("notebook") # poster, talk, notebook, paper

codefolder = "../../"
sys.path.append(codefolder)

__all__ = ["distri","coverage"]

###-----------------------------------------------------------------------------
def distri(pop, grbs,
           var="z", var_range=[0,4], var_log = False, var_name="", tag="",
           weight=1, ax=None, nbin=20, color="tab:blue", color2="red",
           **kwargs):
    """
    Plot the given variable compared to the reference population and the ratio
    of the two.

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
    mask = (var_range[0] <= pop.ref[var]) & (pop.ref[var]<= var_range[1])

    if var_log:
        x = [np.log10(v) if v>0 else 0 for v in pop.ref[mask][var]]
    else:
        x = pop.ref[mask][var]

    nref, bins, _  = ax.hist(x, bins=nbin,
                             facecolor="none",edgecolor="black")

    # Plot the requested population - can be weighted if the sum of several
    mask = (var_range[0] <= grbs[var]) & (grbs[var] <= var_range[1])

    if var_log:
        x  =[ np.log10(v) if v>0 else 0 for v in grbs[mask][var]]
    else:
        x = grbs[mask][var]

    n, bins, _ = ax.hist(x,
                         bins=bins, color=color,
                         label= MyLabel(x, label = tag),
                         weights = weight*np.ones(len(x)),
                         **kwargs)
    ax.legend()
    ax.set_xlabel(var_name)
    ax.set_yscale("log")

    # Plot the ratio
    ratio = [ (x/xtot if xtot!=0 else 0) for x, xtot in zip(n, nref) ]
    axx = ax.twinx()
    axx.plot(bins[:-1] + 0.5*(bins[1:]-bins[:-1]),ratio,
             alpha=1,color=color2,label="Ratio")

    axx.grid(ls=":")
    axx.legend(loc="lower right")

    return ax, axx

###----------------------------------------------------------------------------
def coverage(varx, vary,
             ref = None, pop = None, mask = None,
             xrange = (None,None), yrange = (None, None),
             lblx   = "",          lbly   ="",
             xscale = "log",       yscale ="log",
             title  = "dummy",
             nbin=25, alpha=0.5):
    """


    Parameters
    ----------
    varx : TYPE
        DESCRIPTION.
    vary : TYPE
        DESCRIPTION.
    ref : TYPE, optional
        DESCRIPTION. The default is None.
    pop : TYPE, optional
        DESCRIPTION. The default is None.
    mask : TYPE, optional
        DESCRIPTION. The default is None.
    xrange : TYPE, optional
        DESCRIPTION. The default is [None,None].
    yrange : TYPE, optional
        DESCRIPTION. The default is [None, None].
    lblx : TYPE, optional
        DESCRIPTION. The default is "".
    lbly : TYPE, optional
        DESCRIPTION. The default is "".
    xscale : TYPE, optional
        DESCRIPTION. The default is "log".
    yscale : TYPE, optional
        DESCRIPTION. The default is "log".
    title : TYPE, optional
        DESCRIPTION. The default is "dummy".
    nbin : TYPE, optional
        DESCRIPTION. The default is 25.
    alpha : TYPE, optional
        DESCRIPTION. The default is 0.5.

    Returns
    -------
    None.

    """

    fig, ax, axh, axv = projected_scatter()

    ### --------------------
    ### Central scatter plot
    ### --------------------
    # Generated population

    ax.scatter(np.log10(ref[varx]) if xscale=="log" else ref[varx],
               np.log10(ref[vary]) if yscale=="log" else ref[vary],
               marker=".", color="black",s=10, alpha=0.5, label="All")

    # Detected population
    colors, sizes = col_size(pop[mask].sigmx)
    ax.scatter(np.log10(pop[mask][varx]) if xscale=="log" else pop[mask][varx],
               np.log10(pop[mask][vary]) if yscale=="log" else pop[mask][vary],
               marker="o", s= sizes, c=colors, alpha = alpha)



    plot_historical(ax,historical(), obs=["H.E.S.S"])

    # Decoration
    if xscale=="log":
        ax.set_xlim(xmin=np.log10(xrange[0]), xmax=np.log10(xrange[1]))
        ax.set_xlabel("$log_{10} \ $"+lblx)
    else:
        ax.set_xlim(xmin=xrange[0],xmax=xrange[1])
        ax.set_xlabel(lblx)

    if yscale == "log":
        ax.set_ylim(ymin=np.log10(yrange[0]),ymax=np.log10(yrange[1]))
        ax.set_ylabel("$log_{10} / $" +lbly)
    else:
        ax.set_ylim(yrange)
        ax.set_ylabel(lbly)

    ax.grid("both",ls="--")

    patches = vals_legend(alpha = alpha)
    fig.legend(title="$\sigma_{max}$",handles=patches, bbox_to_anchor=[0.98, 1.01],ncol=2)

    single_legend(ax)

    ### --------------------
    ### horizontal data
    ### --------------------
    hist_mask = (ref[varx] >= xrange[0]) & (ref[varx] <= xrange[1])

    xref = np.log10(ref[hist_mask][varx]) if xscale=="log" else ref[hist_mask][varx]
    _, bins, _ = axh.hist(xref, bins = nbin, facecolor="none",edgecolor="black")

    x = np.log10(pop[mask][varx]) if xscale=="log" else pop[mask][varx]
    axh.hist(x, bins = bins, color="purple",alpha=0.3)

    if xscale == "log":
        axh.set_xlim(xmin = np.log10(xrange[0]), xmax=np.log10(xrange[1]))
    else:
        axh.set_xlim(xmin = xrange[0], xmax=xrange[1])

    axh.set_title(tag+" - " + title)
    axh.set_yscale("log")

    ### --------------------
    ### Vertical data
    ### --------------------
    hist_mask = (ref[vary] >= yrange[0]) & (ref[vary] <= yrange[1])

    yref = np.log10(ref[hist_mask][vary]) if yscale=="log" else ref[hist_mask][vary]
    _, bins, _ = axv.hist(yref, bins = nbin,
                          facecolor="none",edgecolor="black",orientation="horizontal",
                          label=MyLabel(yref ,label="All"))

    y = np.log10(pop[mask][vary]) if yscale=="log" else pop[mask][vary]
    axv.hist(y,bins = bins,
             color="purple",alpha=0.3,orientation="horizontal",
             label = MyLabel(y ,label="Detected"))
    axv.set_xscale(yscale)

###############################################################################
if __name__ == "__main__":

    nyears, files, tag = get_data(parpath=None,debug=True)
    # nyears, files, tag = get_data(parpath="parameter.yaml",debug=False)

    # Read population
    pop   = Pop(files, tag=tag, nyrs= nyears)

    # Population seen in North or South (x2 in number)
    popNS = pop.grb[(pop.grb.loca == "North") | (pop.grb.loca == "South")]

    # Plot one-dim variabale coverage
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(15,10))

    # Loop over the reference population and superimpose selected data
    for icol, grbs, weight, tag in zip(range(2),
                               [popNS, pop.g_tot],
                               [0.5, 1],
                               [r"North & South 5$\sigma$","Combined 5$\sigma$"]):
        print(icol)

        # Selection
        grbs = grbs[grbs.d5s >= pop.eff_lvl]

        distri(pop, grbs, var="z", tag=tag,
               var_range=[0, 4],
               ax= ax[0][icol],
               var_name="Redshift (z)",
               weight = weight)

        distri(pop, grbs ,var="Eiso", tag=tag,
               var_range = [5.5e50, 5.5e55], var_log=True,
               ax=ax[1][icol],
               var_name=r"$log_{}10 E_{iso}$",
               weight = weight)

        distri(pop, grbs ,var="Epeak", tag=tag,
               var_range=[10, 1e5], var_log=True,
               ax=ax[2][icol],
               weight = weight)

        plt.tight_layout()

    ### Plot coverage in 2-dim space
    poplist  = [pop.g_n,pop.g_s,pop.g_tot]
    taglist  = ["North", "South", "Combined"]
    masklist = [pop.g_n.d5s   >= pop.eff_lvl,
                pop.g_s.d5s   >= pop.eff_lvl,
                pop.g_tot.d5s >= pop.eff_lvl]
    title = r"Detected at $5\sigma \ (90\% \ C.L.)$"

    # Eiso versus z
    varx    = "z"
    vary    = "Eiso"
    lblx    = "$z$"
    lbly    = "$E_{iso} \ (erg)$"
    xrange  = [0,5]
    yrange  = [1e50, 5e55]

    alpha = 0.5

    for grbs, mask, tag in zip(poplist, masklist, taglist):
        coverage(varx,vary,ref=pop.ref, pop=grbs,
                      xrange=xrange,yrange=yrange,
                      lblx = lblx, lbly = lbly,
                      xscale="linear",
                      mask=mask,
                      title=title)
        stamp(pop.tag[0], axis=fig, x=1.04, y = 0.5, rotation=270)

    # Eiso versus Epeak
    varx    = "Epeak"
    vary    = "Eiso"
    lblx    = r"$E_{peak} \ (keV)$"
    lbly    = r"$E_{iso} \ (erg)$"
    xrange  = [10,1e5]
    yrange  = [1e50, 5e55]

    for grbs, mask, tag in zip(poplist, masklist, taglist):

        coverage(varx,vary,ref=pop.ref, pop=grbs,
                      xrange=xrange,yrange=yrange,
                      lblx = lblx, lbly = lbly,
                      mask=mask,
                      title=title)
        stamp(pop.tag[0], axis=fig, x=1.04, y = 0.5, rotation=270)
