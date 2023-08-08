# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:47:57 2022

@author: Stolar
"""

import matplotlib.pyplot as plt
from astropy.visualization import quantity_support
from   astroplan import moon_phase_angle, moon_illumination

__all__ = ["moonphase_plot","moon_alt_plot",  "moonlight_plot",
           "moon_dist_plot" ]

###---------------------------------------------------------------------------
def moonphase_plot(times, ax=None, color="red"):

    if ax is None:
        fig, ax = plt.subplots(figsize=(21,5))

    with quantity_support():
        ax.plot(times.datetime, moon_phase_angle(times),
                color=color,label="Brightness")
    ax.set_ylabel("Brightness")
    ax.set_ylim(ymin=0,ymax=1.)
    ax.legend(loc="upper right")

    return ax

###---------------------------------------------------------------------------
def moon_alt_plot(times, alt, ax=None, alpha=1, color="darkblue"):

    if ax is None:
        fig, ax = plt.subplots(figsize=(21,5))

    with quantity_support():
        ax.plot(times.datetime,alt,
                color=color,alpha=alpha,label="Moon altitude")

    ax.set_ylabel("Alt.(°)")
    ax.legend()

    return ax

###---------------------------------------------------------------------------
def moonlight_plot(times, ax=None, color="tab:orange"):

    import matplotlib.pyplot as plt
    from astropy.visualization import quantity_support

    if (ax==None): fig, ax = plt.subplots(figsize=(21,5))

    with quantity_support():
        ax.plot(times.datetime, moon_illumination(times),
                color=color,label="Illumination")

    ax.set_ylabel("Illumination")
    ax.set_ylim(ymin=0,ymax=1.)
    ax.legend(loc="lower right")

    return ax

###---------------------------------------------------------------------------
def moon_dist_plot(radec,times, moon_radec, site="None",
                   ax=None, alpha=1, color="purple"):

    if (ax==None): fig, ax = plt.subplots(figsize=(21,5))

    dist = radec.separation(moon_radec)

    with quantity_support():
        ax.plot(times.datetime,dist,
                color=color,alpha=alpha, label="Moon distance")

    ax.set_ylabel("Dist.")
    ax.legend()

    return ax