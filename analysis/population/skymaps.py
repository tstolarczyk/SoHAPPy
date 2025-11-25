# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 16:04:07 2023.

Show skymap coverage in ra-dec and in altitude-azimuth

@author: Stolar
"""
import sys
import numpy as np

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
sns.set_context("notebook")  # poster, talk, notebook, paper

__all__ = ["radec", "detection_alt_az"]


# ##---------------------------------------------------------------------------
def radec(ginit, gpop,
          ax=None,
          title="No title",
          label="No label",
          scatter=False,
          nbins=100):
    """
    Plot the Ra-Dec distribution of a reference and a selected population.

    Can be a scatter plot or a 2D hsitrogram. In the later case, the reference
    popualtion is not plotted.

    Parameters
    ----------
    ginit : pandas dataframe
        Reference popualtion.
    gpop : pandas dataframe
        Selected population.
    ax : Matplolib axis, optional
        Current axis. The default is None.
    title : String, optional
        Plot title. The default is "No title".
    label : String, optional
        Sub-population label. The default is "No label".

    Returns
    -------
    None.

    """
    # Detected population
    ra_d = [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value
            for x in gpop.ra]
    dec_d = [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value
             for x in gpop.dec]

    if scatter:
        # Reference population
        ra_o = [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value
                for x in ginit.ra]
        dec_o = [Angle(x*u.deg).wrap_at(180*u.deg).to("radian").value
                 for x in ginit.dec]

        with quantity_support():

            ax.scatter(ra_o, dec_o, facecolor="grey", edgecolor="black",
                       marker='.', alpha=0.2, s=10, label='All')

            colors, sizes = col_size(gpop.sigmx)
            ax.scatter(ra_d, dec_d, alpha=0.5, c=colors, s=sizes, label=label)
    else:
        # create bin edges
        ra_edges = np.linspace(-np.pi, np.pi, nbins + 1)
        dec_edges = np.linspace(-np.pi/2., np.pi/2, nbins + 1)

        # calculate 2D histogram, the shape of hist is (bin_number, bin_number)
        hist, ra_edges, dec_edges = np.histogram2d(
            ra_d, dec_d, bins=[ra_edges, dec_edges], density=True)

        img = ax.pcolor(ra_edges[:-1], dec_edges[:-1], hist.T, cmap="PuBu")

        # Using hist2d directly is simpler but lead to an error on trying to
        # set limits and a recommandation to use cartopy
        # _, _, _, img = ax.hist2d(ra_d, dec_d)
        ax.set_title(tag)
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label('Counts')

    ax.set_xlabel("ra (°)")
    ax.set_ylabel("dec (°)")
    ax.set_title(title)
    ax.legend(loc="lower right")
    stamp(pop.tag[0], axis=fig, where="bottom")


# ##-------------------------------------------------------------------------------------
def detection_alt_az(sub_pop, tag="", scatter=True, ax=None,
                     start=False, track=False, arrow=False, **kwargs):
    """
    Plot the Alt-Az distribution of a population.

    Can use a scatter plot or an historgram.

    Parameters
    ----------
    sub_pop : pandas dataframe
        A given subpopulation.
    tag : String, optional
        Descriptor. The default is "".
    ax : matplotlib axis, optional
        Current axis. The default is None.
    start : Boolean, optional
        If True display start position (scatter only). The default is False.
    track : Boolean, optional
        If True dispaly a line showing the source direction in the sky
        (scatter only). The default is False.
    arrow : Boolean, optional
        If True add an arrow to th eline to shwo the direction.
        The default is False.
    **kwargs;
        Additional parameters

    Returns
    -------
    None.

    """
    if ax is None:
        ax = plt.gca()

    alt2 = [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value
            for x in sub_pop.alt2]
    az2 = [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value
           for x in sub_pop.az2]

    if scatter:  # oNly if the number of point is not too large
        colors, sizes = col_size(sub_pop.sigmx)
        ax.scatter(az2, alt2, alpha=0.6, c=colors, s=sizes,
                   label=MyLabel(az2, label=tag))

        if start:  # Display starting altitude
            alt1 = [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value
                    for x in sub_pop.alt1]
            az1 = [Angle(x*u.deg).wrap_at(180*u.deg).to("degree").value
                   for x in sub_pop.az1]
            ax.scatter(az1, alt1,
                       alpha=0.2, marker="o", color="black", s=sizes/3,
                       label=MyLabel(az1))

        if track:  # Display line between start and stop
            for i in range(0, len(az1)):
                lin = mlines.Line2D([az1[i], az2[i]],
                                    [alt1[i], alt2[i]],
                                    ls="--", lw=1.0,
                                    color=colors[i], alpha=0.2)
                ax.add_line(lin)
        if arrow:  # Add direction arrow
            for i in range(0, len(az1)):
                ax.arrow(az2[i], alt2[i], 3, 3, shape='full', lw=1,
                         length_includes_head=True, head_width=.05)
        single_legend(ax.get_figure(), loc="lower right")

    else:  # For large populations
        _, _, _, img = ax.hist2d(az2, alt2,
                                 label=MyLabel(az2, label=tag), **kwargs)
        ax.set_title(tag)
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label('Counts')

    ax.set_xlabel("Azimuth (°)")
    ax.set_ylabel("Altitude (°)")
    ax.set_xlim(xmin=-180, xmax=180)
    ax.set_ylim(ymin=0, ymax=90)


# ##############################################################################################
if __name__ == "__main__":

    # Bigger texts and labels
    sns.set_context("talk")  # poster, talk, notebook, paper

    parpath = "parameter_100k_ISM_alpha.yaml"

    nyears, files, tag = get_data(parpath=parpath, debug=True)
    # nyears, files, tag = get_data(parpath="parameter.yaml", debug=False)

    pop = Pop(files, tag=tag, nyrs=nyears)
    popNS = pop.grb[(pop.grb.loca == "North") | (pop.grb.loca == "South")]

    # Sky coverage - ra-dec
    taglist = ["North", "South",
               "North only", "South only",
               "Both", "All"]
    poplist = [pop.g_n, pop.g_s,
               pop.g_n0, pop.g_s0,
               pop.g_b, pop.g_tot]

    # Choose between histrogram and scatter plot depending on the pop. size
    scatter = False if len(pop.grb) > 1e4 else True

    for g, tag in zip(poplist, taglist):

        fig = plt.figure(figsize=(15, 6))
        ax = fig.add_subplot(111, projection='aitoff')
        # ax = fig.add_subplot(111, projection='mollweide')

        radec(popNS, g, ax=ax,
              label=r"$\sigma_{max}$", title=tag, scatter=scatter, nbins=100)

    # Sky coverage, altitude - azimuth

    # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
    #                                figsize=(20, 10), sharey=True)
    # detection_alt_az(pop.g_n, tag="North", ax=ax1, scatter=scatter, bins=100)
    # detection_alt_az(pop.g_s, tag="South", ax=ax2, scatter=scatter, bins=100)
    # plt.tight_layout()
