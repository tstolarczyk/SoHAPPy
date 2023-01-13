# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 13:54:31 2023

@author: Stolar
"""
import sys
codefolder = "../../"
sys.path.append(codefolder)

import numpy as np

import matplotlib.pyplot as plt
from niceplot import MyLabel

###----------------------------------------------------------------------------
def plot_var(var, poplist, taglist, logx=False, bins=20):
    
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

    # A standalone function to read a GRB and make various tests

    from population import Pop    
    from pop_io import create_csv

    nyears, file, tags = create_csv(file="compare_EBL.yaml",debug=True)
    
    # Read popualtion files
    poplist = []
    for f in file: 
        pop = Pop(filename=f, nyrs= nyears)
        poplist.append(pop)
    print(" ====> Done !")
    
    # Plots
    plot_var("z",poplist,tags)
    plot_var("Eiso",poplist, tags, logx=True)