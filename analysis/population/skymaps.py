# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 16:04:07 2023

Show skymap coverage in ra-dec and in altitude-azimuth

@author: Stolar
"""
import sys


import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

import astropy.units as u
from astropy.visualization import quantity_support
from astropy.coordinates import Angle

from niceplot import col_size, stamp, MyLabel, single_legend

from population import Pop
from pop_io import get_data

codefolder = "../../"
sys.path.append(codefolder)

# Bigger texts and labels
sns.set_context("notebook") # poster, talk, notebook, paper

__all__ = ["radec","detection_alt_az" ]
###----------------------------------------------------------------------------
def radec(ginit, gpop,
           ax = None,
           title = "No title",
           label ="No label"):
    """


    Parameters
    ----------
    ginit : TYPE
        DESCRIPTION.
    gpop : TYPE
        DESCRIPTION.
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is "No title".
    label : TYPE, optional
        DESCRIPTION. The default is "No label".

    Returns
    -------
    None.

    """

    with quantity_support():
        # Reference population
        ra  =  [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value for x in ginit.ra]
        dec =  [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value for x in ginit.dec]
        ax.scatter(ra, dec, facecolor = "grey", edgecolor = "black",
                    marker = '.', alpha = 0.2, s = 10, label     = 'All')

        # Detected population
        ra  =  [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value for x in gpop.ra]
        dec =  [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value for x in gpop.dec]
        colors, sizes = col_size(gpop.sigmx)
        ax.scatter(ra, dec, alpha = 0.5, c= colors, s = sizes, label = label)
    ax.set_xlabel("ra (°)")
    ax.set_ylabel("dec (°)")
    ax.set_title(title)
    ax.legend(loc="lower right")
    stamp(pop.tag[0], axis=fig,where="bottom")

###--------------------------------------------------------------------------------------
def detection_alt_az(sub_pop, tag="", ax = None, start=False, track=False, arrow=False):
    """


    Parameters
    ----------
    sub_pop : TYPE
        DESCRIPTION.
    tag : TYPE, optional
        DESCRIPTION. The default is "".
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    start : TYPE, optional
        DESCRIPTION. The default is False.
    track : TYPE, optional
        DESCRIPTION. The default is False.
    arrow : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """


    # Detection altitude-azumith
    if ax is None:
        ax = plt.gca()

    alt2 =  [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value for x in sub_pop.alt2]
    az2  =  [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value for x in sub_pop.az2]
    colors, sizes = col_size(sub_pop.sigmx)
    _ = ax.scatter(az2, alt2, alpha = 0.6, c= colors, s = sizes,label = MyLabel(az2,label=tag))

    if start: # Display strating altitude
        alt1 =  [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value for x in sub_pop.alt1]
        az1  =  [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value for x in sub_pop.az1]
        ax.scatter(az1, alt1,
                   alpha = 0.2, marker="o", color="black", s = sizes/3,
                   label = MyLabel(az1))

    if track: # Display line between start and stop
        for i in range(0,len(az1)):
            l = mlines.Line2D([az1[i],az2[i]],
                              [alt1[i],alt2[i]],
                              ls="--",lw=1.0,color=colors[i],alpha=0.2)
            ax.add_line(l)
    if arrow: # Add direction arrow
        for i in range(0,len(az1)):
            ax.arrow(az2[i], alt2[i], 3, 3, shape='full', lw=1,
                     length_includes_head=True, head_width=.05)

    # patches=vals_legend(ax)
    single_legend(ax,loc="lower right")
    ax.set_xlabel("Azimuth (°)")
    ax.set_ylabel("Altitude (°)")
    ax.set_xlim(xmin=-180,xmax=180)
    ax.set_ylim(ymin=0,ymax=90)


################################################################################################
if __name__ == "__main__":

    nyears, files, tag = get_data(parpath=None,debug=True)
    # nyears, files, tag = get_data(parpath="parameter.yaml",debug=False)

    pop = Pop(files, tag=tag, nyrs= nyears)
    popNS = pop.grb[(pop.grb.loca=="North") | (pop.grb.loca=="South")]

    # Sky coverage - ra-dec
    taglist = ["North", "South", "North only","South only","Both", "All"]
    poplist = [pop.g_n, pop.g_s, pop.g_n0, pop.g_s0, pop.g_b, pop.g_tot]

    for g, tag in zip(poplist,taglist):

        fig = plt.figure(figsize=(15,6))
        ax  = fig.add_subplot(111,projection='aitoff')

        radec(popNS, g, ax= ax, label="$\sigma_{max}$", title=tag)

    # Sky coverage, altitude - azimuth
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20,10), sharey=True)
    detection_alt_az(pop.g_n,tag="North",ax=ax1)
    detection_alt_az(pop.g_s,tag="South",ax=ax2)
    plt.tight_layout()
