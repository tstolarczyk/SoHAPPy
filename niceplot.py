# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 09:56:58 2022.

A bunch of functions to manipulate objects related to plots.

@author: Stolar
"""
import sys
from collections import OrderedDict
import gammapy

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm


__all__ = ["pause", "draw_contours", "MyLabel", "single_legend", "vals_legend",
           "stamp", "projected_scatter", "ColorMap", "col_size", "draw_sphere"]


# ##---------------------------------------------------------------------------
def pause():
    """
    Pause plot display in interactive mode on a shell script.

    In the
    abscence of a call to that function figures wil stack on the screen during
    the run and all disappear at the end of the run.
    Using this, figures will be stacked on screen and displayed for each event
    until they are closed by the user.

    Returns
    -------
    None.

    """
    plt.show(block=True)


# -----------------------------------------------------------------------------
def lower_limit(xval, yval,
                threshold=0.2,
                norm="sum",
                fitdeg=3,
                nbins=25,
                ax=None,
                plot_lim=True,
                color="red",
                cmap="magma_r",
                cbar=True,
                debug=False):
    """
    From a 2D histrogram obtain the bins limiting the population.

    Finf the lowest empty bins and fit a polynomila function trhough the
    2D-bins found.

    Parameters
    ----------
    xval: numpy array
        x values (float)
    yval: numpy array
        y values (float)
    threshold: float, optional
        Define the acceptabale limit as a relative fraction with respect to
        the a computed value (see "norm") in each column of xval. The default
        is 0.1
    norm: String
        On which quantity the threshold is applied, either on the sum of the
        values in the column, or the maximal value. The default is "sum".
    fitdeg : Interger, optional
        Degree of the polynom used for fitting. The default is 3.
    nbins: integer
        Number of bins in x and y? The default is 25
    ax : matplotlib axis, optional
        Current axis. If undefined, no plot. The default is None.
    prt : boolean, optional
        If True, print count number in each bin. The default is False.
    plot_lim : boolean, optional
        If True, plot the limiting curve.
        The default is True.
    color : string, optional
        limit curve color. The default is "red".
    cmap : string, optional
        Colormap used for the 2D histogram. The default is "magma_r".

    Returns
    -------
    pfit : numpy array
        Parameters (float) of the polynomial function.

    """
    plot = False if ax is None else True

    if plot:
        H, xe, ye, img = ax.hist2d(xval, yval, bins=nbins,
                                   cmap=cmap, norm=mpl.colors.LogNorm())
    else:
        H, xe, ye = np.histogram2d(xval, yval, bins=nbins)

    # Bin centre and heights
    xctr = (xe[1:] + xe[:-1])/2
    yctr = (ye[1:] + ye[:-1])/2
    dy = ye[1:] - ye[:-1]

    ylow = []
    xlow = []

    for ibin in range(len(xctr)):

        # Display counts on each cell
        if debug:
            for jbin in range(len(yctr)):
                ax.text(xctr[ibin], yctr[jbin], str(int(H[ibin, jbin])))

        # Find position of minimal value at bottom of each column
        ycol = H[ibin, :]
        if norm == "sum":
            yref = np.sum(H[ibin, :])
        elif norm == "max":
            yref = np.max(H[ibin, :])
        else:
            sys.exit(f" norm={norm:} not implemented")
        ipos = np.where(ycol/yref > threshold)[0]

        # Print found positions on screen and draw a circle on the plot
        if len(ipos) != 0:
            xlow.append(xctr[ibin])
            ylow.append(yctr[ipos[0]] - 1.5*dy[ipos[0]])
            if debug:
                print("bin  #", ibin, " : ", ycol, ipos, yctr[ipos])
                ax.text(xctr[ibin], yctr[ipos[0]], "O", size=20)
        else:
            if debug:
                print("bin  # No data")

    # Fit low value points with a polynomial
    pfit = np.polyfit(xlow, ylow, fitdeg)

    # from scipy.interpolate import UnivariateSpline
    # spl = UnivariateSpline(xlow, ylow, k=3)

    # ## Plot results if required
    if plot:
        if plot_lim:
            # ax.plot(xlow, ylow, color="black",
            #         lw=1, ls=":", label="lower limits")
            ax.plot(xlow, np.polyval(pfit, xlow),
                    color=color, lw=2, label=" Fit order=" + str(fitdeg))

        ax.legend()

        if cbar:
            fig = plt.gcf()
            cbar = fig.colorbar(img, ax=ax)
            cbar.set_label('Counts')

    return pfit


# -----------------------------------------------------------------------------
def draw_contours(xval, yval, nbins=25, ax=None, **kwargs):
    """
    Plot a contour from the x and y values.

    Example of matplotlib contour parameters:
    levels=10, colors="black", alpha=0.2, ls=":", zorder=100, norm="linear"
    See documentation for details

    Parameters
    ----------
    xval: numpy array
        x values (float)
    yval: numpy array
        y values (float)
    nbins: integer
        Number of bins in x and y? The default is 25
    ax : Matplotlib axes
        Current axis. The default is None
    **kwargs : pointer
        Additionnal parameters to be passed to the matplotlib contour function.

    Returns
    -------
    None.

    """
    if ax is None:
        ax = plt.gca()

    HH, xe, ye = np.histogram2d(xval, yval, bins=nbins)
    grid = HH.transpose()

    midpoints = (xe[1:] + xe[:-1])/2, (ye[1:] + ye[:-1])/2
    cntr = ax.contour(*midpoints, grid, **kwargs)

    return cntr


# ##---------------------------------------------------------------------------
def MyLabel(var, label="", stat="std"):
    """
    Create label for plots and histograms with statistics.

    Add extra statistical information (counts, dispersions) to a classical
    matplotlib label.

    Parameters
    ----------
    var : A list or numpy array
        The plotted variable.
    label : String, optional
        The classical label text. The default is "".
    stat : String, optional
        A keyword defining the extra information to be displayed (on top
        of counts) : "std" for standard deviation, "med" for median, None to
        have counts only. The default is "std".

    Returns
    -------
    legend : String
        The modfied label text.

    """
    if len(label) != 0:
        label = label+"\n"

    legend = "bad option"
    if stat is None:
        return label + f"$n$ : {len(var):d}"

    if stat.find("std") != -1:
        legend = label  \
                + f"$n$ : {len(var):d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$\sigma$ : "+" {:5.3f}".format(np.std(var))
    elif stat.find("med") != -1:
        legend = label  \
                + "$n$ : {:d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$Med.$ : "+" {:5.3f}".format(np.median(var))
    else:  # Display only the sum
        legend = label + ": {:d}".format(len(var))

    return legend


# ##---------------------------------------------------------------------------
def single_legend(fig, debug=False, **kwargs):
    """
    Remove duplicated labels in legend.

    (e.g. a vertical and an horizontal
    lines defining and intersection having the same labels) and group all
    labels in the same text box.
    Works also on a figure with multiple plots.

    Parameters
    ----------
    fig : matplotlib figure
        Figure to be treated.
    debug : boolean, optional
        If True, describe the figure axes. The default is False.
    **kwargs : list of arguments
        Any list of arguments accepted by the matplotlib legend function

    Returns
    -------
    None.

    """
    # Collect all legend labels from the figure
    all_handles = []
    all_labels = []
    for axis in fig.axes:
        handles, labels = axis.get_legend_handles_labels()
        for handle in handles:
            all_handles.append(handle)
        for label in labels:
            all_labels.append(label)

    # Create a dictionnary and remove duplicated labels
    by_label = OrderedDict(zip(all_labels, all_handles))

    if debug:
        print(all_handles)
        print(all_labels)
        print(by_label)

    fig.legend(by_label.values(), by_label.keys(), **kwargs)

    # Remove all individual labels
    for axis in fig.axes:
        lgd = axis.get_legend()
        if lgd is not None:
            lgd.remove()


# ##---------------------------------------------------------------------------
def old_single_legend(ax, **kwargs):
    """
    Remove duplicated labels in legend.

    (e.g. a vertical and an horizontal
    lines defining and intersection having the same labels).
    Replace by a new fucntion. Temporarily kept for backward compatibility.

    Parameters
    ----------
    ax : matplotlib.axes
        Current axis.

    Returns
    -------
    None.

    """
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), **kwargs)


# ##---------------------------------------------------------------------------
def vals_legend(vals=None, alpha=0.5, var_max=1000, **kwargs):
    """
    Create Matplotlib patches with colored circles from a list of values.

    Parameters
    ----------
    vals: list of float
        List of values in legend
    alpha : float, optional
        Alpha value. The default is 0.5.
    var_max : float, optional
        Maximal value. The default is 1000.
    **kwargs : Dictionnary
        Extra parameters.

    Returns
    -------
    patches : legend elements
        Legend associated to the current axis.

    """
    if vals is None:
        vals = [5, 10, 20, 50, 100, 500]

    labels = [str(x) for x in vals]
    # symbol = Line2D(range(1), range(1), color="white", marker='o',
    # markerfacecolor="red")

    colors, sizes = col_size(vals, var_max=var_max)
    patches = [plt.plot([], [], marker="o", alpha=alpha,
               ms=7 + sizes[i]/50,
               ls="", mec=None,
               color=colors[i],
               label="{:s}"
               .format(labels[i]))[0]for i in range(len(labels), **kwargs)]

    return patches


# ##----------------------------------------------------------------------------
def stamp(text, axis=None,
          where="right", x=None, y=None, rotation=0,
          **kwargs):
    """
    Annotate the side of any plot referred from the axis.

    including the gammapy version.

    Parameters
    ----------
    text : String
       Text to be displayed
    where : String, optional
       Position of the text with respect to the axis.
       The default is `right`.
    x : float, optional
       Text x position, supersedes the `where` variable.
       The default is `None`.
    y : float, optional
       Text x position, supersedes the `where` variable.
       The default is `None`.
    rotation : float, optional
       Text rotation - If not given use default.
       The default is `None`.
    axis : matplotlib.axes, optional
       Current plot axis. The default is `None`.
    **kwargs :
       Any additionnal arguments for matplotlib.axes.text

    Returns
    -------
    None.

    """
    text = text + " - " + gammapy.__version__

    if axis is None:
        axis = plt.gcf()  # Currentr figure

    if x is None or y is None:
        if where == "right":
            (x, y) = (1, 0.5)
            rotation = 270
        elif where == "left":
            (x, y) = (0, 0.5)
            rotation = 90
        elif where == "top":
            (x, y) = (0.5,   1)
            rotation = 0
        elif where == "bottom":
            (x, y) = (0.5,   0)
            rotation = 0

    axis.text(x=x, y=y, s=text,
              horizontalalignment='center',
              verticalalignment="center",
              rotation=rotation, **kwargs)


# ##---------------------------------------------------------------------------
def projected_scatter(xsize=12,  ysize=8,
                      left=0.1, width=0.7,
                      bottom=0.1, height=0.7,
                      spacing=0.02):
    """
    Matplotlib template to display a scatter plot and projections.

    Adapted from :
    https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html#sphx-glr-gallery-lines-bars-and-markers-scatter-hist-py

    Parameters
    ----------
    xsize : float, optional
        Figure horizontal size. The default is 12.
    ysize : float, optional
        Figure vertical size. The default is 8.
    left : float, optional
        Scatter plot left position. The default is 0.1.
    width : float, optional
        Height of the histogram handling the y-projection. The default is 0.7.
    bottom : float, optional
        Scatter plot bottom position. The default is 0.1.
    height : float, optional
        Height of the histogram handling the x-projection. The default is 0.7.
    spacing : float, optional
        Space between the subplots. The default is 0.02.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    ax : matplotlib.axes
        Scatter plot axis.
    axh : matplotlib.axes
        Horizontal projection axis.
    axv : matplotlib.axes
        Vertical projection axis.

    """
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    fig = mpl.pyplot.figure(figsize=(xsize, ysize))
    ax = fig.add_axes(rect_scatter)
    axh = fig.add_axes(rect_histx, sharex=ax)
    axv = fig.add_axes(rect_histy, sharey=ax)

    axh.tick_params(axis="x", labelbottom=False)
    axv.tick_params(axis="y", labelleft=False)

    return fig, ax, axh, axv


# ##---------------------------------------------------------------------------
def ColorMap(threshold, maxval):
    """
    Create a colormap based on a threshold and a maximal value.

    Parameters
    ----------
    threshold : float
        Starting values.
    maxval : float
        End value.

    Returns
    -------
    mymap : matplotlib.colors.LinearSegmentedColormap
        A matplotlib color map.

    """
    color_cut = int(256*threshold/maxval)
    print("color_cut=", color_cut)
    colors1 = plt.cm.Greys(np.linspace(0, 1, color_cut))
    colors2 = plt.cm.plasma(np.linspace(0, 1, 256-color_cut))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    return mymap


# ##---------------------------------------------------------------------------
def col_size(var, var_min=1.1, var_max=1000):
    """
    Get a color and a size from a numerical value.

    Parameters
    ----------
    var : float
        Input number.
    var_min : float, optional
        Minimal value. The default is 1.1.
    var_max : float, optional
        Maximal value. The default is 1000.

    Returns
    -------
    color : matplotlib color
        Color.
    size : float
        matplotlib size.

    """
    # Limit values in case they are not yet limited
    var = np.clip(var, var_min, None)

    color = cm.cool(np.log10(var)/np.log10(var_max))
    size = 100*np.log10(var)**2

    return color, size


# ##---------------------------------------------------------------------------
def draw_sphere(radius=1, colormap=plt.cm.viridis, ax=None, **kwargs):
    """
    Draw a sphere in matplotlib.

    Parameters
    ----------
    radius : float, optional
        radius. The default is 1.
    colormap : matplolib color map, optional
        Color map. The default is plt.cm.viridis.
    ax : matplolib axes, optional
        Current axis. The default is None.
    **kwargs :
        Extra arguments.

    Returns
    -------
    None.

    """
    if ax is None:
        fig = plt.figure(figsize=(8, 8), dpi=300)
        ax = fig.add_subplot(111, projection='3d')

    if mpl.__version__ > "3.5":
        ax.set_box_aspect(aspect=(1, 1, 1))

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

    ls = mcolors.LightSource(azdeg=0, altdeg=65)
    rgb = ls.shade(z, colormap)

    ax.plot_surface(x, y, z,  rstride=1, cstride=1,
                    facecolors=rgb, linewidth=0, **kwargs)
