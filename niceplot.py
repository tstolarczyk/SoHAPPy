# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 09:56:58 2022

A bunch of functions to manipulate objects related to plots.

@author: Stolar
"""

from collections import OrderedDict
import gammapy

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm

__all__ = ["pause", "MyLabel","single_legend","vals_legend","stamp",
           "projected_scatter","ColorMap","col_size","draw_sphere"]
###----------------------------------------------------------------------------
def pause():
    """
    Used to pause plot display in interactive mode on a shell script. In the
    abscence of a call to that function figures wil stack on the screen during
    the run and all disappear at the end of the run.
    Using this, figures will be stacked on screen and displayed for each event
    until they are closed by the user.

    Returns
    -------
    None.

    """

    plt.show(block=True)

###----------------------------------------------------------------------------
def MyLabel(var, label="", stat="std"):

    """
    A label for plots and histograms with statistics.
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

    if len(label) != 0 :
        label = label+"\n"

    legend="bad option"
    if stat is None:
        return label + f"$n$ : {len(var):d}"

    if stat.find("std") != -1:
        legend = label  \
                + f"$n$ : {len(var):d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$\sigma$ : "+" {:5.3f}".format(np.std(var))
    elif stat.find("med") !=-1 :
        legend = label  \
                + "$n$ : {:d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$Med.$ : "+" {:5.3f}".format(np.median(var))

    return legend

###----------------------------------------------------------------------------
def single_legend(ax, **kwargs):

    """
    Remove duplicated labels in legend (e.g. a vertical and an horizontal
    lines defining and intersection having the same labels).

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
    ax.legend(by_label.values(), by_label.keys(),**kwargs)

###----------------------------------------------------------------------------
def vals_legend(vals = None, alpha = 0.5, var_max = 1000, **kwargs):

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

    print("niceplot.vals_legend: Investigate why ax is not used")
    if vals is None:
        vals = [5 , 10, 20, 50, 100, 500]

    labels = [str(x) for x in vals]
    # symbol = Line2D(range(1), range(1), color="white", marker='o',
    # markerfacecolor="red")

    colors, sizes = col_size(vals,var_max=var_max)
    patches = [ plt.plot([],[], marker="o", alpha=alpha,
                          ms = 7+sizes[i]/50,
                          ls="", mec=None,
                          color=colors[i],
                          label="{:s}".format(labels[i]) )[0]  for i in range(len(labels),
                          **kwargs) ]

    return patches

###----------------------------------------------------------------------------
def stamp(text, axis = None,
          where = "right", x = None, y = None, rotation = 0,
          **kwargs):
    """
    Annotation on the side of any plot referred from the axis, including
    the gammapy version.


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

    if x is None or y is None:
        if   where =="right":
            (x,y) = (  1, 0.5)
            rotation = 270
        elif where =="left":
            (x,y) = (  0, 0.5)
            rotation = 90
        elif where =="top":
            (x,y) = (0.5,   1)
            rotation = 0
        elif where =="bottom":
            (x,y) = (0.5,   0)
            rotation = 0

    axis.text(x=x,y=y,s=text,
              horizontalalignment='center',
              verticalalignment="center",
              rotation=rotation, **kwargs)

###----------------------------------------------------------------------------
def projected_scatter(xsize   = 12,  ysize  = 8,
                      left    = 0.1, width  = 0.7,
                      bottom  = 0.1, height = 0.7,
                      spacing = 0.02):
    """
    Matplotlib template to display a scatter plot and the horizontal and
    vertical projections of the data. Adapted from :
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
    spacing : TYPE, optional
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
    rect_histx   = [left, bottom + height + spacing, width, 0.2]
    rect_histy   = [left + width + spacing, bottom, 0.2, height]

    fig = matplotlib.pyplot.figure(figsize=(xsize, ysize))
    ax  = fig.add_axes(rect_scatter)
    axh = fig.add_axes(rect_histx, sharex=ax)
    axv = fig.add_axes(rect_histy, sharey=ax)

    axh.tick_params(axis="x", labelbottom=False)
    axv.tick_params(axis="y", labelleft=False)

    return fig, ax, axh, axv

###----------------------------------------------------------------------------
def ColorMap(threshold,maxval):
    """
    Create a colormap based on a threshold and a maximal value

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
    print("color_cut=",color_cut)
    colors1  = plt.cm.Greys(np.linspace(0, 1, color_cut))
    colors2  = plt.cm.plasma(np.linspace(0, 1, 256-color_cut))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    return mymap

###----------------------------------------------------------------------------
def col_size(var,var_min=1.1,var_max = 1000):

    """
    get a color and a size from a numerical value

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
    var   = np.clip(var, var_min, None)

    color = cm.cool(np.log10(var)/np.log10(var_max))
    size  = 100*np.log10(var)**2

    return color, size

###----------------------------------------------------------------------------
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
        fig = plt.figure(figsize=(8,8), dpi=300)
        ax = fig.add_subplot(111, projection='3d')

    if matplotlib.__version__ > "3.5":
        ax.set_box_aspect(aspect = (1,1,1))

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = radius* np.outer(np.cos(u), np.sin(v))
    y = radius* np.outer(np.sin(u), np.sin(v))
    z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

    ls = mcolors.LightSource(azdeg=0, altdeg=65)
    rgb = ls.shade(z, colormap)

    ax.plot_surface(x, y, z,  rstride=1, cstride=1,
#                     color=color,
                    facecolors=rgb, linewidth=0, **kwargs)
