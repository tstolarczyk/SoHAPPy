# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 14:47:57 2022.

@author: Stolar
"""

import matplotlib.pyplot as plt
from astropy.visualization import quantity_support
from astroplan import moon_phase_angle, moon_illumination

__all__ = ["moonphase_plot", "moon_alt_plot", "moonlight_plot",
           "moon_dist_plot"]


# ##---------------------------------------------------------------------------
def moonphase_plot(times, ax=None, color="red"):
    """."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(21, 5))

    with quantity_support():
        ax.plot(times.datetime, moon_phase_angle(times),
                color=color, label="Brightness")
    ax.set_ylabel("Brightness")
    ax.set_ylim(ymin=0, ymax=1.)
    ax.legend(loc="upper right")

    return ax


# ##---------------------------------------------------------------------------
def moon_alt_period_plot(times, alt, ax=None, alpha=1, color="orange"):
    """."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(21, 5))

    with quantity_support():
        ax.bar(times.datetime, alt, width=(times[1] - times[0]).value,
               color=color, alpha=alpha, label="Moon altitude")
    ax.set_ylabel("Altitude (°)")
    ax.legend()
    return ax


# ##---------------------------------------------------------------------------
def moon_alt_plot(times, alt, ax=None, alpha=1, color="darkblue"):
    """."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(21, 5))

    with quantity_support():
        ax.plot(times.datetime, alt,
                color=color, alpha=alpha, label="Moon altitude")

    ax.set_ylabel("Alt.(°)")
    ax.legend()

    return ax


# ##---------------------------------------------------------------------------
def moonlight_plot(times, ax=None, color="tab:orange", norm=False, tag=None):
    """."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(21, 5))

    illumination = moon_illumination(times)  # Fraction of Moon illiminated

    # Normalize illumination to current axis range
    if norm is True:
        ymin, ymax = ax.get_ylim()
        renorm = ymin + (ymax - ymin)*illumination
    else:
        renorm = 1.

    if tag is None:
        tag = "Phase"

    with quantity_support():
        ax.plot(times.datetime, renorm*illumination,
                color=color, label=tag)

    ax.set_ylabel("Illumination")  # Can be superseded later
    ax.legend()

    return ax


# ##---------------------------------------------------------------------------
def moon_dist_plot(radec, times, moon_radec, site="None",
                   ax=None, alpha=1, color="purple", tag=None):
    """."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(21, 5))

    dist = radec.separation(moon_radec)

    if tag is None:
        tag = "Distance"
    with quantity_support():
        ax.plot(times.datetime, dist,
                color=color, alpha=alpha, label=tag)

    ax.set_ylabel("Distance")
    ax.legend()

    return ax
