# -*- coding: utf-8 -*-
"""
Explore population from output files.

The output population files contain by default all the GRB information
available in the original fits file.

Created on Mon Jan  9 14:08:58 2023

@author: Stolar
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import astropy.units as u
from   astropy.time import Time, TimeDelta
from astropy.visualization import quantity_support
from  astropy.coordinates import Angle, SkyCoord

from population import Pop
from pop_io import get_data

from niceplot import MyLabel, single_legend
from niceprint import heading


# Bigger texts and labels
sns.set_context("notebook") # poster, talk, notebook, paper

CODE = "../../"
sys.path.append(CODE)

__all__ = ["galactic_coverage", "var_coverage"]

###----------------------------------------------------------------------------
def galactic_coverage(grb, latmax= 2*u.deg):

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(15,7))

    with quantity_support():

        # Wide
        ax1.hist(grb.lat, bins=1000,label=MyLabel(grb.lat))
        ax1.set_xlabel(None)
        # ax1.legend(bbox_to_anchor=[1,1])

        # Zoomed
        mask = abs(grb.lat)<latmax.value
        latzoom  = grb.lat[mask]
        namelist = grb[mask].name

        ax2.plot(latzoom, np.ones(len(latzoom)),marker="o",ls="",label=MyLabel(grb.lat))
        i=0
        for l, txt in zip(latzoom, namelist):
            ax2.text(l,1.01+(i%5)*0.015,s=txt[5:],rotation=90)
            ax2.set_ylim(ymax=1.1)
            i+=1

        ax2.set_xlim([-1.05*latmax,1.05*latmax])
        ax2.axvline( latmax.value,color="red",ls="--",label="Lat. limit")
        ax2.axvline(-latmax.value,color="red",ls="--",label="Lat. limit")
        ax2.set_xlabel("Latitude (deg)")
        single_legend(ax1)

    plt.tight_layout(h_pad=0)

###----------------------------------------------------------------------------
def var_coverage(var,grb,logz=False):

    fig = plt.figure(figsize=(20,8))

    gs  = gridspec.GridSpec(1, 2, width_ratios=[2, 1])

    if logz: varz = np.array([np.log10(x) for x in grb[var]])
    else:    varz = grb[var]
    sizes = 20+ (50*  (varz - min(varz))/(max(varz)-min(varz)))**2
    ax1 = fig.add_subplot(gs[0], projection="mollweide")
    sc = ax1.scatter(grb.ra*np.pi/180, grb.dec*np.pi/180,marker='o', c=varz , s=sizes, cmap='rainbow',alpha=0.5)

    ax1.set_title("GRB "+var+" distribution")
    cbar = fig.colorbar(sc,  orientation='vertical')
    cbar.set_label(var, rotation=0)

    ax2 = fig.add_subplot(gs[1])
    ax2.hist(grb[var],alpha=0.5)

    plt.tight_layout()
###############################################################################
if __name__ == "__main__":

    nyears, files, tag = get_data(parpath=None,debug=True)
    # nyears, files, tag = get_data(parpath="parameter.yaml",debug=False)

    pop = Pop(files, tag=tag, nyrs= nyears)

    # Wrap ra and dec for further use, define galactic coordinates - note that units handling is not clear
    pop.ref.ra  = np.array([Angle(x*u.deg).wrap_at(180*u.deg).value for x in pop.ref.ra ])
    pop.ref.dec = np.array([Angle(x*u.deg).wrap_at(180*u.deg).value for x in pop.ref.dec])
    c = SkyCoord(ra=pop.ref.ra*u.deg, dec=pop.ref.dec*u.deg, frame='icrs')
    pop.ref["lat"] = c.galactic.b # float

    # Dates covered
    heading("Dates and durations") #======================

    dmin = Time(np.min(pop.ref.t_trig)*u.day,format="mjd",scale="utc")
    dmax = Time(np.max(pop.ref.t_trig)*u.day,format="mjd",scale="utc")
    duration = dmax - dmin
    print(" Simulation period from ",dmin.datetime," to ",dmax.datetime)
    print("   - Duration : ",TimeDelta(duration).to_value(u.yr)," years")

    # Sky coverage - Note possible GRBs in the galactic plane
    galactic_coverage(pop.ref)

    # Some variable covreage - Eiso and z
    var_coverage("z",pop.ref)
    var_coverage("Eiso",pop.ref,logz=True)