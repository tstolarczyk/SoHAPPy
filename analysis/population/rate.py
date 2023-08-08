# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:26:59 2022

@author: Stolar
"""
import numpy as np

__all__ = ["compute","separator","stat_mean","prt_mean",
           "stat_line"]
###----------------------------------------------------------------------------
def compute(pop, summary=False):
    """
    Compute and display the detection rates of the various subpopulations

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
    print("",107*"-")
    if pop.nyears:
        print(f" Normalized to {pop.nyears:} year",end="")
        print("s") if pop.nyears>1 else print()

    print("",107*"-")
    print(f" Rate : {'N':>15} {'S':>15}"\
          f"{'Nonly':>16} {'Sonly':>15} {'Both':>15} {'Total':>15}")

    if not summary:
        stat_line(pop,tag="Vis.",nyrs=pop.nyears)

    # # # Analysed
    # poplist = [g[g.err == niter] for g in poplist]
    # if not summary: stat_line(poplist,tag="Ana.",ny=nyears)

    # 3 sigma mean detection
    stat_mean(pop,tag="3s",nyrs=pop.nyears)

    # 5 sigma mean detection
    stat_mean(pop,tag="5s",nyrs=pop.nyears)

    # 3 sigma 90%CL detected
    if not summary:
        stat_line(pop,tag="3s 90%CL",nyrs=pop.nyears)

    # 5 sigma 90%CL detected
    stat_line(pop,tag="5s 90%CL",nyrs=pop.nyears)

    print("",107*"-")

###----------------------------------------------------------------------------
def separator():
    """
    Just a separator line
    """
    print(f" ----------- {14*'-':>15} {14*'-':>15}",end="")
    print(f"{14*'-':>16} {14*'-':>15} {14*'-':>15} {14*'-':>15}")

###----------------------------------------------------------------------------
def stat_mean(pop,tag="",nyrs=1):
    """
    Print mean detection above a given significance for a list of
    subpopulations

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
    separator()
    tag2 = tag
    if nyrs !=1:
        prt_mean(pop, nyrs=1,sig=tag, tag=tag)
        print()
        tag2 =""
    prt_mean(pop, nyrs=nyrs, sig=tag, tag=tag2)
    print()

###----------------------------------------------------------------------------
def prt_mean(pop, nyrs=1, sig=None, tag=None):
    """
    Mean detection above a certain significance
    """

    if nyrs!=1 :
        print(f" {'yr-1':>9s} :",end="") ### !!!
    else:
        print(f" {tag:9s} :",end="") ### !!!

    glist = [pop.g_n, pop.g_s, pop.g_n0, pop.g_s0, pop.g_b,
             [pop.g_n0,pop.g_s0,pop.g_b]]

    for grb in glist:
        if sig == "3s":
            if isinstance(grb,list):
                var = sum([sum(x.d3s) for x in grb])
            else:
                var = sum(grb.d3s)
        elif sig == "5s":
            if isinstance(grb,list):
                var = sum([sum(x.d5s) for x in grb])
            else:
                var= sum(grb.d5s)
        else:
            sys.exit(" Should be '3s' or '5s'")

        nmean = var/pop.niter
        print(f" {nmean/nyrs:7.1f} +- "\
              f"{np.sqrt(nmean)/nyrs:4.1f}",end="")

#--------------------------------------------------------
def stat_line(pop,tag="",nyrs=1):
    """
    Print statistics in various conditions
    """

    glist = [pop.g_n, pop.g_s, pop.g_n0, pop.g_s0, pop.g_b]
    if tag.find("3s") != -1:
        glist = [g[g.d3s >= pop.eff_lvl] for g in glist]
    if tag.find("5s") != -1:
        glist = [g[g.d5s >= pop.eff_lvl] for g in glist]

    glist = glist + [glist[2]+glist[3]+glist[4]]

    #--------------------------------------------------------
    # Stat line internal functions
    #--------------------------------------------------------
    def prt_line(nyrs=1,tag="dummy"):
        if nyrs==1:
            print(f" {tag:9s} :",end="") ### !!!
        else:
            print(f" {'yr-1':>9s} :",end="") ### !!!
        for grb in glist:
            print(f" {len(grb)/nyrs:7.1f} "\
                  f"+- {np.sqrt(len(grb))/nyrs:4.1f}",end="" )
    #--------------------------------------------------------
    def prt_vis():
        print(f" {'@trig':>9s} :",end="")
        for grb in glist:
            nvis = len(grb[grb.prpt==1])
            ratio = 100*nvis/len(grb) if len(grb) else 0
            print(f" {nvis:7d}  {ratio:5.1f}%",end="")
    #--------------------------------------------------------

    if nyrs !=1:
        separator()
        prt_line(nyrs=1,tag = tag)
        print()
        prt_line(nyrs=nyrs, tag=" ")
        print()
        prt_vis()
    else:
        separator()
        prt_line(nyrs=nyrs, tag=tag)
    print()

###############################################################################
if __name__ == "__main__":

    from pop_io import create_csv
    from population import Pop
    import sys
    CODE = "../../"
    sys.path.append(CODE)

    nyears, file, _ = create_csv(file="parameter.yaml",debug=True)
    pop = Pop(filename=file, nyrs= nyears)
    compute(pop,summary=False)
