# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:26:59 2022.

Compute detection rates for various site-related population normalised to
the data taking duration to get mean yearly rates above 3 and 5 standard
deviations of detection significance.

@author: Stolar
"""
import sys
import numpy as np
from pop_io import get_data
from population import Pop

CODE = "../../"
sys.path.append(CODE)

__all__ = ["compute", "separator", "stat_mean", "prt_mean", "stat_line"]

pminus = "+-"

# ##---------------------------------------------------------------------------
def compute_to_output(pop, summary=False):
    """
    Compute and display the detection rates of the various subpopulations.

    Parameters
    ----------
    pop : Pop instance
        Data read from disk.
    summary : boolean, optional
        True for less information is printed. The default is False.

    Returns
    -------
    None.

    """

    # Check that N only an S only tags exist
    suppinfo = ("N" in pop.grb) and ("S" in pop.grb) and ("B" in pop.grb)
    if not suppinfo:
        sys.exit(" North/South combinatory does not exist")

    # Header
    if pop.nyears:
        print(f" Normalized to {pop.nyears:} year", end="")
        print("s") if pop.nyears > 1 else print()

    print(f" Rate      : {'N':>18} {'S':>19}"
          f" {'Nonly':>19} {'Sonly':>19} {'Both':>19} {'Total':>19}")

    if not summary:
        stat_line(pop, tag="Vis.", nyrs=pop.nyears, sep=False)

    # # # Analysed
    # poplist = [g[g.err == niter] for g in poplist]
    # if not summary: stat_line(poplist,tag="Ana.",ny=nyears)

    # 3 sigma mean detection
    stat_mean(pop, tag="3s", nyrs=pop.nyears, sep=False)

    # 5 sigma mean detection
    stat_mean(pop, tag="5s", nyrs=pop.nyears, sep=False)

    # 3 sigma 90%CL detected
    if not summary:
        stat_line(pop, tag="3s 90%CL", nyrs=pop.nyears, sep=False)

    # 5 sigma 90%CL detected
    stat_line(pop, tag="5s 90%CL", nyrs=pop.nyears, sep=False)


# ##---------------------------------------------------------------------------
def compute(pop, summary=False):
    """
    Compute and display the detection rates of the various subpopulations.

    Parameters
    ----------
    pop : Pop instance
        Data read from disk.
    summary : boolean, optional
        True for less information is printed. The default is False.

    Returns
    -------
    None.

    """
    # Check that N only an S only tags exist
    suppinfo = ("N" in pop.grb) and ("S" in pop.grb) and ("B" in pop.grb)
    if not suppinfo:
        sys.exit(" North/South combinatory does not exist")

    # Header
    print()
    print("" + (132 - len(str(pop.tag)) - 1)*"-" + " " + str(pop.tag))
    if pop.nyears:
        print(f" Normalized to {pop.nyears:} year", end="")
        print("s") if pop.nyears > 1 else print()

    print("", 132*"-")
    print(f" Rate      : {'N':>18} {'S':>19}"
          f" {'Nonly':>19} {'Sonly':>19} {'Both':>19} {'Total':>19}")

    if not summary:
        stat_line(pop, tag="Vis.", nyrs=pop.nyears)

    # # # Analysed
    # poplist = [g[g.err == niter] for g in poplist]
    # if not summary: stat_line(poplist,tag="Ana.",ny=nyears)

    # 3 sigma mean detection
    stat_mean(pop, tag="3s", nyrs=pop.nyears)

    # 5 sigma mean detection
    stat_mean(pop, tag="5s", nyrs=pop.nyears)

    # 3 sigma 90%CL detected
    if not summary:
        stat_line(pop, tag="3s 90%CL", nyrs=pop.nyears)

    # 5 sigma 90%CL detected
    stat_line(pop, tag="5s 90%CL", nyrs=pop.nyears)

    print("", 132*"-")


# ##---------------------------------------------------------------------------
def separator():
    """Just a separator line."""
    print(f" ----------- {18*'-':>19} {18*'-':>19}", end="")
    print(f"{18*'-':>20} {18*'-':>19} {18*'-':>19} {18*'-':>19}")


# ##---------------------------------------------------------------------------
def stat_mean(pop, tag="", nyrs=1, sep=True):
    """
    Print mean detection above a given significance for a population set.

    Parameters
    ----------
    glist : population list
        list of populations to be analysed.
    tag : String, optional
        "3s" or "5s" for 3 and 5 sigma respectively. The default is "".
    nyrs : integer, optional
        Number of years, normalisation. The default is 1.

    Returns
    -------
    None.

    """
    if sep:
        separator()
    tag2 = tag
    if nyrs != 1:
        prt_mean(pop, nyrs=1, sig=tag, tag=tag)
        print()
        tag2 = ""
    prt_mean(pop, nyrs=nyrs, sig=tag, tag=tag2)
    print()


# ##---------------------------------------------------------------------------
def prt_mean(pop, nyrs=1, sig=None, tag=None):
    """Mean detection above a certain significance."""
    if nyrs != 1:
        print(f" {'yr-1':>9s} :", end="")
    else:
        print(f" {tag:9s} :", end="")

    glist = [pop.g_n, pop.g_s, pop.g_n0, pop.g_s0, pop.g_b,
             [pop.g_n0, pop.g_s0, pop.g_b]]

    for grb in glist:
        if sig == "3s":
            if isinstance(grb, list):
                var = sum([sum(x.d3s) for x in grb])
            else:
                var = sum(grb.d3s)
        elif sig == "5s":
            if isinstance(grb, list):
                var = sum([sum(x.d5s) for x in grb])
            else:
                var = sum(grb.d5s)
        else:
            sys.exit(" Should be '3s' or '5s'")

        nmean = var/pop.niter
        print(f" {nmean/nyrs:9.1f} {pminus:2s} "
              f"{np.sqrt(nmean)/nyrs:6.1f}", end="")


# ##---------------------------------------------------------------------------
def stat_line(pop, tag="", nyrs=1, sep=True):
    """Print statistics in various conditions."""
    glist = [pop.g_n, pop.g_s, pop.g_n0, pop.g_s0, pop.g_b]
    if tag.find("3s") != -1:
        glist = [g[g.d3s >= pop.eff_lvl] for g in glist]
    if tag.find("5s") != -1:
        glist = [g[g.d5s >= pop.eff_lvl] for g in glist]

    glist = glist + [glist[2]+glist[3]+glist[4]]

    # --------------------------------------------------------
    # Stat line internal functions
    # --------------------------------------------------------
    def prt_line(nyrs=1, tag="dummy"):
        if nyrs == 1:
            print(f" {tag:9s} :", end="")
        else:
            print(f" {'yr-1':>9s} :", end="")
        for grb in glist:
            print(f" {len(grb)/nyrs:9.1f} {pminus:2s} "
                  f"{np.sqrt(len(grb))/nyrs:6.1f}", end="")

    # --------------------------------------------------------
    def prt_vis_prompt():
        print(f" {'@trig':>9s} :", end="")
        for grb in glist:
            nvis = len(grb[grb.prpt == 1])
            ratio = 100*nvis/len(grb) if len(grb) else 0
            print(f" {nvis:9d}  {ratio:7.1f}%", end="")
    # --------------------------------------------------------

    if nyrs != 1:
        if sep:
            separator()
        prt_line(nyrs=1, tag=tag)
        print()
        prt_line(nyrs=nyrs, tag=" ")
        print()
        prt_vis_prompt()
    else:
        if sep:
            separator()
        prt_line(nyrs=nyrs, tag=tag)
    print()


# ##---------------------------------------------------------------------------
def dump_alldata(pop, filename):
    with open(filename, 'a') as f:
        print(" Writing data out")
        f.write(pop.grb.to_string(header=True, index=False))


# #############################################################################
if __name__ == "__main__":

    # Use default as a demo
    # papath = None
    # parpath = "parameter_100k_ISM-max.yaml"
    parpath = "parameter_100k_ISM_alpha.yaml"
    # parpath = "parameter_long1000.yaml"
    # parpath = "parameter_100k_ISM_omega.yaml"

    nyears, files, tag = get_data(parpath=parpath, debug=True)

    pop = Pop(files, tag=tag,  nyrs=nyears, fpeakmin=1.0, debug=False)

    pop.print()
    compute(pop, summary=False)

    pminus = "  "
    compute_to_output(pop, summary=False)

    # Dump all data to a single file - this can be quite long!
    # dump_alldata(pop, "output_all_data100000.txt")
    print("...completed!")
