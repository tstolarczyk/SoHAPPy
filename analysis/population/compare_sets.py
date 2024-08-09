# -*- coding: utf-8 -*-
"""
Compare two (or several sets)
Created on Tue Jan 10 14:04:54 2023

@author: Stolar
"""
import sys
import collections
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt

from niceplot import MyLabel
from niceprint import heading

from times import detection_fraction_versus_time

from population import Pop
from pop_io import get_data
codefolder = "../../"
sys.path.append(codefolder)

__all__ = ["find_5sigma_diff", "check_visibilities", "var_plot",
           "var_scatter_all", "var_scatter"]


# ##---------------------------------------------------------------------------
def find_5sigma_diff(pop1, pop2):
    """
    Find and display sources being detected at 5 sigma in one
    simulation/analysis, not the other.

    Parameters
    ----------
    pop1 : Pop instance
        First population data.
    pop2 : Pop instance
        Second population data.

    Returns
    -------
    None.

    """

    grb1 = pop1.g_tot
    grb2 = pop2.g_tot

    first1 = True
    first2 = True

    # ----------------------
    def show():
        print(f"{g1['name']:8s} {g1.loca:5s} "
              + f"1: sigmx= {g1.sigmx:5.1f} +/-{g1.esigmx:5.1f}"
              + f" * CL={g1.d5s:<5d} * vis={g1.prpt:2d} *  err={g1.err:5d}")
        print(f" {'':8s} {'':5s}"
              + f"2:        {g2.sigmx:5.1f} +/-{g2.esigmx:5.1f}"
              + f"      {g2.d5s:<5d} *     {g2.prpt:2d} *      {g2.err:5d}")
    # ----------------------

    for (_, g1), (_, g2) in zip(grb1.iterrows(), grb2.iterrows()):

        if (g1.d5s >= pop1.eff_lvl and g2.d5s < pop2.eff_lvl):
            if first1:
                print("Error : detected in 1 not in 2")
                first1 = False
            show()

        if (g1.d5s < pop1.eff_lvl and g2.d5s >= pop2.eff_lvl):
            if first2:
                print("Error : detected in 2 not in 1")
                first2 = False
            show()

    print("All done!")


# ##----------------------------------------------------------------------------
def check_visibilities(pop1, pop2):
    """
    Check sources visible in one population and not in the other.

    Parameters
    ----------
    pop1 : Pop instance
        First population data.
    pop2 : Pop instance
        Second population data.

    Returns
    -------
    None.

    """
    grb1 = pop1.grb
    grb2 = pop2.grb

    namelist = set(grb1["name"].values)
    count = {"North": 0, "South": 0, "Both": 0}

    #     print("check : ",name)
    for loc in ["North", "South", "Both"]:
        print(f"{24*'-':24s}")
        print(f"{loc:10s} err1 err2")
        for name in namelist:
            mask1 = (grb1["name"] == name) & (grb1.loca == loc)
            mask2 = (grb2["name"] == name) & (grb2.loca == loc)
            if grb1[mask1].err.values[0] != grb2[mask2].err.values[0]:
                print(f"{name:10s} {grb1[mask1].err.values[0]:4d}"
                      + f" {grb2[mask2].err.values[0]:4d}")
                count[loc] += 1
    print(" Differences between the 2 populations :", count)


# ##---------------------------------------------------------------------------
def var_plot(var1, var2, tag=("o", "o"), varmin=None, varmax=None,
             nbin=100, xscale="log", yscale="log", xlabel=""):
    """
    Show distributions of a given variable of two different populations.

    Parameters
    ----------
    var1 : Sequence
        First population variable data sequence.
    var2 : Sequence
        Second population variable data sequence.
    tag : List of strings, optional
        Tags to appear on the x and y axis respectively.
        The default is ["o","o"].
    varmin : float, optional
        Minimal allowed value for the variable. The default is None.
    varmax : float, optional
        Maximal allowed value for the variable . The default is None.
    nbin : integer, optional
        Number of bins for both axes. The default is 100.
    xscale : string, optional
        "log" or "linear" for the  abscissa. The default is "log".
    yscale : string, optional
        "log" or "linear" for the ordinate. The default is "log".
    xlabel : string, optional
        Abscissa label describing the variable. The default is "".

    Returns
    -------
    None.

    """

    if varmin is None:
        varmin = min(min(var1), min(var2))
    if varmax is None:
        varmax = max(max(var1), max(var2))
    print("min = ", varmin, " max=", varmax)

    mask1 = (var1 <= varmax) & (var1 >= varmin)
    mask2 = (var2 <= varmax) & (var2 >= varmin)

    _, (axa, axb) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

    # First population
    _, bins, _ = axa.hist(var1[mask1], bins=nbin, alpha=0.5,
                          label=MyLabel(var1[mask1], tag[0]))

    # Second population
    axa.hist(var2[mask2], bins=bins, alpha=0.5,
             label=MyLabel(var2[mask2], tag[1]))

    # Decoration
    axa.set_xscale(xscale)
    axa.set_yscale(yscale)
    axa.set_xlabel(xlabel)
    axa.legend()

    # Ratio betwenn the two plots
    axb.hist(var2[mask2]/var1[mask1], bins=nbin, alpha=0.5,
             label=MyLabel(var2[mask2]/var1[mask1], tag[1]+"/"+tag[0]))
    axb.set_xlabel(xlabel+" ratio")
    axb.set_yscale(yscale)
    axb.legend()


# ##----------------------------------------------------------------------------
def var_scatter_all(var, pop1, pop2, title=""):
    """
    For this study, the two sets should have the same size. We start from
    the original popualtion set (even if not visible)

    Parameters
    ----------
    var : String
        Column name to be plotted.
    pop1 : Pop instance
        First popualtion.
    pop2 : Pop instance
        Second population.

    Returns
    -------
    None.

    """

    gn1 = pop1.grb[pop1.grb.loca == "North"]
    gs1 = pop1.grb[pop1.grb.loca == "South"]
    gb1 = pop1.grb[pop1.grb.loca == "Both"]
    gn2 = pop2.grb[pop2.grb.loca == "North"]
    gs2 = pop2.grb[pop2.grb.loca == "South"]
    gb2 = pop2.grb[pop2.grb.loca == "Both"]

    _, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 6), sharey=True)

    first = True
    for ax0, g1, g2, title in zip(ax, [gn1, gs1, gb1],
                                      [gn2, gs2, gb2],
                                      ["North - " + title,
                                       "South - " + title,
                                       "Both - " + title]):
        # var1 = g1[nmin:nmax][var]
        # var2 = g2[nmin:nmax][var]
        var1 = g1[var]
        var2 = g2[var]

        var_scatter(var1, var2, tag=[pop1.tag, pop2.tag], title=title, ax=ax0)

        if first:
            first = False
        else:
            ax0.set_ylabel(None)

    plt.tight_layout(h_pad=0)


# ##---------------------------------------------------------------------------
def var_scatter(var1, var2, tag=("o", "o"), title="",
                varmin=None, varmax=None,
                logscale=True, cut=0, ax=None):
    """
    The two variable sets should have the same size

    Parameters
    ----------
    var1 : List
        First variable.
    var2 : List
        Second variable.
    tag : List of string, optional
        DESCRIPTION. The default is ["o","o"].
    varmin : Float, optional
        Minimal value to be plotted. The default is None.
    varmax : Float, optional
        Maximal value to be plotted. The default is None.
    nbin : Integer, optional
        Bin number. The default is 100.
    logscale : Boolean, optional
        If True plot the logarithm of the variables. The default is True.
    ax : Matplotlib axes, optional
        Current axis. The default is None.

    Returns
    -------
    ax : Matplotlib axes
       Current axis.

    """

    if len(var1) != len(var2):
        print(" Set 1 : ", tag[0], len(var1))
        print(" Set 2 : ", tag[0], len(var2))
        print(" Samples should have the same size !")
        return ax

    if varmin is None:
        varmin = min(min(var1), min(var2))
    if varmax is None:
        varmax = max(max(var1), max(var2))
    # print("min = ",varmin," max=",varmax)

    if logscale:
        var1 = [max(cut, np.log10(x)) if x > 0 else cut for x in var1]
        var2 = [max(cut, np.log10(x)) if x > 0 else cut for x in var2]

        varmin = max(cut, np.log10(varmin)) if varmin > 0 else cut
        varmax = np.log10(varmax) if varmax > 0 else cut

    sig3 = np.log10(3) if logscale else 3
    sig5 = np.log10(5) if logscale else 5

    if ax is None:
        plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

    ax.scatter(var1, var2, marker="+")

    # Diagonal
    ax.plot([varmin, varmax], [varmin, varmax], ls=":", color="red")

    # 3 and 5 sigma lines
    ax.axvline(sig3, label=r"$ 3 \sigma$", color="tab:orange", ls=":")
    ax.axhline(sig3, label=r"$ 3 \sigma$", color="tab:orange", ls=":")
    ax.axvline(sig5, label=r"$ 5 \sigma$", color="tab:green", ls=":")
    ax.axhline(sig5, label=r"$ 5 \sigma$", color="tab:green", ls=":")

    ax.set_xlabel(tag[0])
    ax.set_ylabel(tag[1])

    ax.set_xlim(varmin, varmax)
    ax.set_ylim(varmin, varmax)

    ax.set_title(title)
    ax.grid(which="minor", axis="both", ls=":", alpha=0.5)
    ax.grid(which="major", axis="both", ls="-", alpha=0.5)

    handles, labels = ax.get_legend_handles_labels()
    by_label = collections.OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

    return ax


##############################################################################
if __name__ == "__main__":


    import yaml
    from yaml.loader import SafeLoader

    import os
    os.environ["HAPPY_IN"] = r"D:\\CTAO\SoHAPPy\input"
    os.environ["HAPPY_OUT"] = r"D:\\CTAO\SoHAPPy\output"

    # parpath = Path(codefolder, "data/samples/comparepop_parameter.yaml")
    parpath = Path("comparepop_parameter.yaml")

    # Separate the file list as they are indepenedent data
    xdict = yaml.load(open(parpath.as_posix()), Loader=SafeLoader)
    filelist = xdict["outfolders"]
    nyears = xdict["duration"]
    tags = xdict["tags"]

    nyears, filelist, tags = get_data(parpath=parpath, debug=False)

    poplist = []
    for file, tag in zip(filelist, tags):
        pop = Pop(file, tag=tag, nyrs=nyears)
        pop.sanity_check()
        poplist.append(pop)

    [pop1, pop2] = [poplist[0], poplist[1]]

    # Find 5 sigma detection differences in the two population
    heading("Compare 5 sigma detection")
    find_5sigma_diff(pop1, pop2)

    # # Check differences in visibility
    heading("Check visibilities")
    check_visibilities(pop1, pop2)

    # ###------------------------
    # ### Check rates
    # ###------------------------
    import rate
    pop1.print()
    pop2.print()
    rate.compute(pop1)
    rate.compute(pop2)

    # ###------------------------
    # ### Check some variables
    # ###------------------------

    nmin = 0
    nmax = 3000
    sigmin = 0

    gs1 = pop1.g_s
    gs2 = pop2.g_s

    # Check tmax relative error
    var1 = gs1[nmin:nmax][gs1.etmx > 0].etmx/gs1[nmin:nmax].tmx
    var2 = gs2[nmin:nmax][gs2.etmx > 0].etmx/gs2[nmin:nmax].tmx
    var_plot(var1, var2, tag=[pop1.tag, pop2.tag],
             xscale="linear", yscale="log",
             xlabel="$t_{max}$ relative error")
    plt.tight_layout()

    # Mean maximal significance
    var1 = pop1.g_tot[nmin:nmax][pop1.g_tot.sigmx > 0].sigmx
    var2 = pop2.g_tot[nmin:nmax][pop2.g_tot.sigmx > 0].sigmx
    var_plot(var1, var2, tag=[pop1.tag, pop2.tag],
             xscale="linear", yscale="log",
             xlabel=r"$\sigma_{max}$ ")
    plt.tight_layout()

    # Mean maximal excess number
    var1 = pop1.g_tot[nmin:nmax][pop1.g_tot.nexmx > 0].nexmx
    var2 = pop2.g_tot[nmin:nmax][pop2.g_tot.nexmx > 0].nexmx
    var_plot(var1, var2, tag=[pop1.tag, pop2.tag], nbin=250,
             xscale="linear", yscale="log",
             xlabel=r"$N_{ex}$ at max.")
    plt.tight_layout()

    # Mean maximal background number
    var1 = pop1.g_tot[nmin:nmax][pop1.g_tot.nbmx > 0].nexmx
    var2 = pop2.g_tot[nmin:nmax][pop2.g_tot.nbmx > 0].nexmx
    var_plot(var1, var2, tag=[pop1.tag, pop2.tag], nbin=250,
             xscale="linear", yscale="log",
             xlabel=r"$N_{b}$ at max.")

    plt.tight_layout()

    # Mean time of maximal significance
    var1 = gs1[nmin:nmax][gs1.sigmx > 0].tmx
    var2 = gs2[nmin:nmax][gs2.sigmx > 0].tmx

    var_plot(var1, var2, tag=[pop1.tag, pop2.tag],
             xscale="linear", yscale="log",
             xlabel="$t_{max}$ ")
    print(" Set 1 : ", len(var1))
    print(" Set 2 : ", len(var2))
    plt.tight_layout()

    # ##------------------------
    # ## Compare in scatter plots
    # ##------------------------

    # heading("Compare sigmax")
    var_scatter_all("sigmx", pop1, pop2, title=r"$\sigma_{max}$")

    # heading("Compare excess counts")
    var_scatter_all("nexmx", pop1, pop2, title="Excess counts")

    # heading("Compare background counts")
    var_scatter_all("nbmx", pop1, pop2, title="background counts")

    # ##------------------------
    # ## Cumulative 5 sigma detectionversus time
    # ##------------------------
    axs = detection_fraction_versus_time(pop1, pop1.g_tot,
                                         color="tab:blue",
                                         label=pop1.tag,
                                         title=" All")
    detection_fraction_versus_time(pop2, pop2.g_tot,
                                   color="tab:orange",
                                   label=pop2.tag,
                                   title="All", axs=axs)
