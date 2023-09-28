# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:12:11 2023
 Compare analysed data obtained from with the EBLmodel available in `SoHAPPy`.

@author: Stolar
"""

import sys
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt
from niceplot import MyLabel

from population import Pop
from pop_io import get_data


codefolder = "../../"
sys.path.append(codefolder)


__all__ = ["plot_var"]
###----------------------------------------------------------------------------
def plot_var(var, poplist, taglist, logx=False, bins=20):
    """


    Parameters
    ----------
    var : TYPE
        DESCRIPTION.
    poplist : TYPE
        DESCRIPTION.
    taglist : TYPE
        DESCRIPTION.
    logx : TYPE, optional
        DESCRIPTION. The default is False.
    bins : TYPE, optional
        DESCRIPTION. The default is 20.

    Returns
    -------
    None.

    """


    fig, axlist = plt.subplots(nrows=1, ncols=len(poplist),
                           figsize=(5*len(poplist),6),sharey=True)

    first = True
    for pop, tag, ax0 in zip(poplist, taglist, axlist):

        x = pop.g_tot[(pop.g_tot.d5s>=pop.eff_lvl)][var]
        if logx: x = np.log10(x)

        hdata = ax0.hist(x, bins=bins,
                              alpha=0.7, color="tab:orange", edgecolor="grey",lw=2,
                              label=MyLabel(x,label=tag))
        if first:
            plot_ref = hdata
            bins = hdata[1]
            first = False
        else:
            n, bins, bars = plot_ref
            ax0.bar(bins[:-1]+0.5*(bins[1:]-bins[:-1]),
                    n,
                    width=(bins[1:]-bins[:-1]),
                    facecolor="none",edgecolor="grey",lw=2,alpha=1)

        ax0.legend()
        ax0.set_xlabel(var)
        ax0.grid(ls=":",alpha=0.5)

    plt.tight_layout()

##############################################################################
if __name__ == "__main__":

    parpath = Path(codefolder,"data/samples/compare_ebl.yaml")
    nyears, file, tags = get_data(parpath=parpath,debug=True)

    # Read popualtion files
    poplist = []
    for f in file:
        pop = Pop(f, tag=tags, nyrs= nyears)
        poplist.append(pop)
    print(" ====> Done !")

    # Plots
    plot_var("z",poplist,tags)
    plot_var("Eiso",poplist, tags, logx=True)