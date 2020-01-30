# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 09:23:08 2020

@author: Stolar
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from astropy.table import Table, vstack
from IPython.display import display

import sys
sys.path.append("../../../../utilities_ths/")  # path is where this module is
from utilities_ths import MyLabel
from utilities_ths import stamp

plt.style.use('seaborn-talk') # Make the labels readable
#plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!

__all__ = ['stat','computing_time']
###----------------------------------------------------------------------------
def stat(grb, g_vis, g_ana,  g_abrt, g_3s, g_5s, i3s, i5s, north, south):
    
    nevt = len(grb)
    print(" This is a file with ",nevt," events, i.e. ",nevt/2," GRBs")
    print("")
    print("           : {:>7} {:>5} {:>5}".format("N or S","N","S"))
    print("---------- : {:>7} {:>5} {:>5}".format("-------","-----","-----"))
    print("Total      : {:>7d} {:>5d} {:>5d}".format(nevt,len(grb[north]),len(grb[south])))
    print("Visible    : {:>7d} {:>5d} {:>5d}".format(len(g_vis),len(g_vis[north]),len(g_vis[south])))
    print("Aborted    : {:>7d} {:>5d} {:>5d}".format(len(g_abrt),len(g_abrt[north]),len(g_abrt[south])))
    print("Full simul : {:>7d} {:>5d} {:>5d}".format(len(g_ana),len(g_ana[north]),len(g_ana[south])))
    print("---------- : {:>7} {:>5} {:>5}".format("-------","-----","-----"))
    print("3s 90%CL   : {:>7d} {:>5d} {:>5d}".format(len(g_3s),len(g_3s[north]),len(g_3s[south])))
    print("5s 90%CL   : {:>7d} {:>5d} {:>5d}".format(len(g_5s),len(g_5s[north]),len(g_5s[south])))
    
    return

###----------------------------------------------------------------------------
def computing_time(file,gpop,i3s, i5s,north,south,niter=100):
    nbin=25
    niter = niter
    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(15,6))
    
    n, bins, _ = ax1.hist(niter*gpop.mct,bins=nbin,
                          facecolor="none",
                          edgecolor="black",
                          label=MyLabel(niter*gpop.mct,stat="med"))
    ax1.hist(niter*gpop[north].mct,bins=bins,
             alpha=0.5,
             label=MyLabel(niter*gpop[north].mct,"North",stat="med"))
    ax1.hist(niter*gpop[south].mct,bins=bins,
             alpha=0.5,
             label=MyLabel(niter*gpop[south].mct,"South",stat="med"))
    ax1.set_title(f"Simulation duration N/S ({niter:3d} iter.)")
    ax1.set_xlabel("Time (s)")
    ax1.legend(fontsize=10)
    stamp(ax1,file)
    ax1.legend()
    
    n, bins, _ = ax2.hist(niter*gpop.mct,bins=nbin,
                          facecolor="none",
                          edgecolor="black",
                          label=MyLabel(niter*gpop.mct,stat="med"))
    ax2.hist(niter*gpop[i3s].mct,bins=bins,
             alpha=0.5,
             label=MyLabel(niter*gpop[i3s].mct,"$3\sigma$ det.",stat="med"))
    ax2.hist(niter*gpop[i5s].mct,bins=bins,
             alpha=0.5,
             label=MyLabel(niter*gpop[i5s].mct,"$5\sigma$ det.",stat="med"))
    ax2.set_title(f"Simulation duration $3, 5\sigma$ ({niter:3d} iter.)")
    ax2.set_xlabel("Time (s)")
    ax2.set_yscale('log')
    stamp(ax2,file)
    ax2.legend(fontsize=10)

    plt.show()

    return