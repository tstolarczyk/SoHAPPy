# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

This module gather functions for plotting the `dataset` contents.

@author: Stolar
"""

import sys

import itertools
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from astropy.time import Time
from astropy.visualization import quantity_support
from niceprint import t_str

sys.path.append("../")
sys.path.append("../../../SoHAPPy")

__all__ = ['panels', 'windows', 'counts_from_plot']


# -----------------------------------------------------------------------------
def windows(dsets, nbin=15, ysize=5, unit="s", plot=True):
    """
    Roughly display the existing time wndows in the Dataset list.
    In a stacked dataset, GTI are merged (only) if they are contiguous.

    Parameters
    ----------
    dsets : Dataset list
        Current dataset list.
    nbin : integer, optional
        Number of bins. The default is 25.
    ysize : float, optional
        Plot width. The default is 5.
    unit : astropy.unit string, optional
        Time unit. The default is "s".
    plot : Boolean, optional
        If True, plot the data. The default is True.

    Returns
    -------
    None.

    """

    print(f"{'Id.':<3s} | {'Sum':>10s} | "
          f"{'Start':^23s} {'Stop':^23s} | "
          f"{'dt':<3s}")

    for iset, ds in enumerate(dsets):
        t0 = ds.gti.time_start[0]
        t1 = ds.gti.time_stop[0]
        dt = [t_str(t) for t in ds.gti.time_delta]
        ttot = t_str(ds.gti.time_sum)

        print(f"{iset:<3d} | {ttot:>10s} | "
              f"{t0.isot:23s} {t1.isot:23s} | {dt:}")

    if plot is True:

        with quantity_support():

            fig, (ax, ax1) = plt.subplots(nrows=2, ncols=1,
                                          figsize=(15, 2*ysize))

            ax.hist([(ds.gti.time_delta.to(unit).value) for ds in dsets],
                    bins=nbin)
            dt_short = np.array([ds.gti.time_delta.to(unit).value
                                 for ds in dsets])
            ax.set_xlabel("Start times (s)")
            ax.set_ylabel("Occurences")

            ax1.set_xlabel("Observation periods $[t_{start}, t_{stop}]$")

            tmin = np.concatenate(np.array([(ds.gti.time_start).jd
                                            for ds in dsets]))
            tmax = np.concatenate(np.array([(ds.gti.time_stop).jd
                                            for ds in dsets]))

            it = 0
            for t1, t2 in zip(tmin, tmax):
                # print(t1,t2)
                color = cm.cool(it/len(tmin))
                ax1.axvspan(xmin=t1, xmax=t2, ymin=0, ymax=1-it/len(tmin),
                            color=color, alpha=0.5, label=str(it))
                it += 1
            if len(tmin) < 6:
                ax1.legend()

            dt_short = np.concatenate(dt_short)
            dt_short = dt_short[dt_short < 500]
            if len(dt_short):
                axx = inset_axes(ax, width="50%", height=1.2,
                                 loc="upper right")
                axx.hist(dt_short, bins=nbin)
                axx.set_xlim(xmax=200)
                axx.set_xlabel("Time (s)")

        plt.tight_layout(h_pad=0)
    return fig


# -----------------------------------------------------------------------------
def panels(dsets,
           func=None,
           nmaxcol=4, xsize=5, ysize=5,
           xscale="log", yscale="log",
           fixrange=True,
           tref=Time('2000-01-01T00:00:00', format='isot', scale='utc'),
           **kwargs):
    """
    Display the plots obtained from a function applied to each individual
    dataset in an optimised panel.

    Parameters
    ----------
    dsets : Datasets object
        A dataset collection.
    func : Function, optional
        The pointer ot a function with first argument a dataset object, and
        a matplotib axis. The rest of the arguments are passed through the
        keyword arguments (`kwargs`). The default is None.
    nmaxcol : Integer, optional
        The maximum number of columns in the panel. The default is 5.
    xsize : Float, optional
        The wifth of the columns in the panel. The default is 5.
    ysize : Float, optional
        The height of the rows in the panel. The default is 7.
    fixrange : Boolean, optional
        If True, all plot y axis are identical. The default is False.
    kwargs : Keyword arguments
        Extra arguments for the function (func)

    Returns
    -------
    None.

    """

    # Panel geometry
    nplots = len(dsets)
    ncols = min(nmaxcol, nplots)  # If nplots < nmaxcol, take nplots
    nrows = int(nplots/ncols) + 1*(nplots % ncols != 0)

    fig, ax = plt.subplots(ncols=ncols, nrows=nrows,
                           figsize=(xsize*ncols, ysize*nrows),
                           sharex=True, sharey=fixrange)

    with quantity_support():
        iplot = 0

        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            if nplots != 1:
                ax0 = ax[jrow][icol] if (nrows > 1) else ax[icol]
            else:
                ax0 = ax

            # If nothing to plot, blank space instead of empty plot
            if iplot >= nplots:
                ax0.axis('off')
                continue  # Next plot

            # Subset of datasets
            ds = dsets[iplot]

            # Function to be plotted - creates thelabel related to the
            # elapsed time
            dt = ds.gti.time_start[0] - tref

            label = "$t$ = " + t_str(dt) + "\n" + r"$\Delta$t = "\
                + t_str(ds.gti.time_sum)

            func(ds, ax=ax0, tag=label, **kwargs)

            # Compactify
            if jrow+1 != nrows:
                ax0.set_xlabel(None)
            if icol != 0:
                ax0.set_ylabel(None)
            ax0.tick_params(which='major', length=10, width=2, direction='in')
            ax0.tick_params(which='minor', length=5, width=2, direction='in')

            iplot += 1

            ax0.set_xscale(xscale)
            ax0.set_yscale(yscale)

        plt.subplots_adjust(right=0.8, hspace=0.05)

    return fig


# -----------------------------------------------------------------------------
def counts_from_plot(ax, dset, plot_line=False):
    """
    This a utility to get the counts form the plot and check the count_str
    displays the content of the dataset

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The Gammapy dataset count plot to be analysed
    dest: Dataset object
        The original Dataset object to be compared with
    plot_line : Boolean, optional
        If True, plot vertical and horizontal lines to the extracted points.
        The default is False.

    """

    Emu = ax.lines[0].get_xdata()
    ymu = ax.lines[0].get_ydata()
    # Edata = ax.lines[1].get_xdata()
    ydata = ax.lines[1].get_ydata()
    i = 0
    print("--- Counts from plot ", 35*"-")
    print(f"{'bin':3s} - {'E':>8s}     {' Model':>8s} {'Excess':>8s}")
    for _, y, yd in zip(Emu, ymu, ydata):
        print("{i:3d} - {E:8.2f} {y:8.2f} {yd:8.2f}")
        if plot_line:
            ax.axhline(y, ls=":", color="grey")
            ax.axhline(yd, ls=":", color="red")
        i += 1

    print(" Total  excess counts = ", ydata.sum())
    print(" Masked excess counts = ", ydata[3:-1].sum(),
          "from dset ", dset.excess.data[dset.mask_safe].sum())
    print(" Total counts         = ", ymu.sum())
    print(" Masked counts        = ", ymu[3:-1].sum())
