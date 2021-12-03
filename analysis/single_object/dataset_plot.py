# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

@author: Stolar
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

import astropy # for type
from astropy.time import Time
from astropy.visualization import quantity_support

from utilities import t_str

__all__ = ['panels', 'windows']

#------------------------------------------------------------------------------
def windows(dsets, nbin=15, ysize=5, unit="s", plot=True):
    """
    Roughly display the existing time wndows in the Dataset list.
    In a stacked dataset, GTI are merged (only) if they are contiguous
   

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

    for i, ds in enumerate(dsets):
        dt = [t_str(t) for t in ds.gti.time_delta]
        print(" Set {:3d}: sum:{:>10} | start={}, stop={}, dt={}"
              .format(i, t_str(ds.gti.time_sum),
                         ds.gti.time_start, ds.gti.time_stop, dt))
    if plot:
        import matplotlib.cm as cm

        with quantity_support():

            fig, (ax,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(15,2*ysize))
            
            ax.hist([(ds.gti.time_delta.to(unit).value) for ds in dsets],
                    bins=nbin)
            dt_short = np.array([(ds.gti.time_delta.to(unit).value) for ds in dsets])
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Occurences")
            ax1.set_xlabel("Observation periods")
            tmin = np.concatenate(np.array([(ds.gti.time_start).jd for ds in dsets]))
            tmax = np.concatenate(np.array([(ds.gti.time_stop).jd for ds in dsets]))
            it = 0
            for t1,t2 in zip(tmin,tmax):
#                     print(t1,t2)
                color = cm.cool(it/len(tmin))
                ax1.axvspan(xmin=t1, xmax=t2, ymin=0, ymax=1-it/len(tmin),color=color,alpha=0.5,label=str(it))
                it+=1
            if len(tmin)<6:
                ax1.legend()

            dt_short = np.concatenate(dt_short)
            dt_short = dt_short[dt_short<500]
            if (len(dt_short)):
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                axx = inset_axes(ax, width="50%", height=1.2,loc="upper right")
                axx.hist(dt_short,bins=nbin)
                axx.set_xlim(xmax=200)
                axx.set_xlabel("Time (s)")
        plt.tight_layout()
    return fig

#------------------------------------------------------------------------------
def panels(dsets,func = None,
           nmaxcol  = 4, xsize = 5, ysize = 5,
           xscale   = "log", yscale ="log",
           max_margin = 1.1,
           fixrange = True,
           tref     = Time('2000-01-01T00:00:00',format='isot', scale='utc'),
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
        keyword arguments (**kwargs). The default is None.
    nmaxcol : Integer, optional
        The maximum number of columns in the panel. The default is 5.
    xsize : Float, optional
        The wifth of the columns in the panel. The default is 5.
    ysize : Float, optional
        The height of the rows in the panel. The default is 7.
    fixrange : Boolean, optional
        If True, all plot y axis are identical. The default is False.
    **kwargs : Keyword arguments
        Extra arguments for the function (func)

    Returns
    -------
    None.

    """

    # Panel geometry
    nplots = len(dsets)
    ncols  = min(nmaxcol,nplots) # If nplots < nmaxcol, take nplots
    nrows  = int(nplots/ncols)+ 1*(nplots%ncols != 0)
    
    # Create figure
    # if xsize*ncols>17:
    #     xsize = 17/ncols
    #     print("Plot automatically resized")
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows,
                           figsize=(xsize*ncols,ysize*nrows),
                           sharex=True, sharey=fixrange)
    # Plot common y limit - if fixrange os True
    ymax = -np.Inf
    ymin =  np.Inf

    with quantity_support():
        iplot = 0
        import itertools

        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            if nplots !=1:
                ax0 = ax[jrow][icol] if (nrows>1) else ax[icol]
            else:
                ax0=ax

            # If nothing to plot, blank space instead of empty plot
            if (iplot >= nplots):
                ax0.axis('off')
                continue # Next plot

            # Plot function here
            t_since_trigger = dsets[iplot].gti.time_start[0] - tref
            (y1, y2) = func(dsets[iplot], ax=ax0, elapsed=t_since_trigger, **kwargs)
            if y1 != y2:
                if (y2 > ymax): ymax = y2
                if (y1 < ymin): ymin = y1

            # Compactify
            if (jrow+1 != nrows): ax0.set_xlabel(None)
            if (icol !=0): ax0.set_ylabel(None)
            ax0.tick_params(which='major', length=10, width=2, direction='in')
            ax0.tick_params(which='minor', length=5, width=2, direction='in')

            iplot+=1

            ax0.set_xscale(xscale)
            ax0.set_yscale(yscale)
            
        if (fixrange):
            plt.subplots_adjust(right=0.8,hspace=0.05,wspace=0.05)
            for ax in fig.get_axes():
                ymax = max_margin*ymax # If scale is log, it is applied to the log value
                if type(ymax) == astropy.units.quantity.Quantity:
                    ymax = ymax.value
                if type(ymin) == astropy.units.quantity.Quantity:
                     ymin = ymin.value                   
                ax.set_ylim(bottom=ymin, top=ymax)
        else:
            plt.subplots_adjust(right=0.8,hspace=0.05)
    
        # use first label handles as the last plot can be empty
    #     h0, lbl0 = fig.axes[0].get_legend_handles_labels()
    #     fig.legend(h0, lbl0, loc="center right", borderaxespad=0.1)
        # if (fixrange): fig.tight_layout(h_pad=0, w_pad=0)
        # else: fig.tight_layout(h_pad=0, w_pad=0)

    return fig

#------------------------------------------------------------------------------
def counts_from_plot(ax, dset, plot_line=False):
    """
    This a utility to get the counts form the plot and check the count_str
    displays the content of the dataset

    Parameters
    ----------
    ax : Matplotlib axis
        The Gammapy dataset count plot to be analysed
    dest: Dataset object
        The original Dataset object to be compared with
    plot_line : Boolean, optional
        If True, plot vertical and horizontal lines to the extracted points.
        The default is False.

    Returns
    -------
    None.

    """

    Emu = ax.lines[0].get_xdata()
    ymu = ax.lines[0].get_ydata()
    #Edata = ax.lines[1].get_xdata()
    ydata = ax.lines[1].get_ydata()
    i=0
    print("--- Counts from plot ", 35*"-")
    print("{:3s} - {:>8s}     {:>8s} {:>8s}".format("bin","E"," Model", "Excess"))
    for E,y,yd in zip(Emu,ymu,ydata):
        print("{:3d} - {:8.2f} {:8.2f} {:8.2f}".format(i, E,y,yd))
        if (plot_line):
            ax.axhline(y,ls=":",color="grey")
            ax.axhline(yd,ls=":",color="red")
        i+=1

    print(" Total  excess counts = ",ydata.sum())
    print(" Masked excess counts = ",ydata[3:-1].sum(),
          "from dset ",dset.excess.data[dset.mask_safe].sum())
    print(" Total counts         = ",ymu.sum())
    print(" Masked counts        = ",ymu[3:-1].sum())

    return

############################################################################### 
if __name__ == "__main__":
    """
    Code example 
    """
    import pickle
    from pathlib import Path
    from utilities import stamp
    from dataset_tools import createonoff_from_simulation
    from dataset_flux import extract_spectrum
    from gammapy.utils.random import get_random_state
    import yaml
    from yaml.loader import SafeLoader
    
    #plt.style.use('seaborn-talk') # Make the labels readable
    plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!   
    
    ###----------------------    
    ### Get data
    ###----------------------    
    base = Path("D:/000_Today/SoHAPPy_tests/single_sources/")
    
    array   = "omega" # "omega", "alpha"
    GRB_id  = "980"
    site    = "South" # "North", "South"
    merging = False
    
    # Read parameters
    par = (yaml.load(open("parameter.yaml"), Loader=SafeLoader))
    GRB = GRB_id+site[0] # The key to find back the parameters
    
    # Openfile
    file    = Path(base,array,"Event"+GRB_id+"-"+site+"_sim.bin")
    infile  = open(file,"rb")
    mc      = pickle.load(infile)
    infile.close()
    filetag = "GRB "+GRB_id+" - South ($\\"+array +"$)"
    dset_init = createonoff_from_simulation(mc,
                                            random_state=get_random_state(2021))
    # check_datasets(dset_init,deeper=True,masked=True)
    fit_tag = par[GRB]["fit_tag"] # None
    fig = panels(dset_init[22:], func=extract_spectrum, stacked=merging,
                xsize     = 6, yscale="log",
                nmaxcol   = par[GRB]["nmaxcol"], 
                e_ref     = u.Quantity(par[GRB]["e_ref"]), 
                tref      = mc.slot.grb.t_trig,
                e_unit    = "TeV", 
                flux_unit = "TeV-1 cm-2 s-1",
                flux_min  = 1.e-14, 
                flux_max  = float(par[GRB]["flux_max"]), 
                fit_tag   = fit_tag,
                tag = "",
                debug     = True)
    stamp(GRB, axis=fig)
    