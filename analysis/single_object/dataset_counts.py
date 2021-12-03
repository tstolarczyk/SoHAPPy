# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

@author: Stolar
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

from utilities import t_str

__all__ = ['excess_counts','residuals',
           'excess_versus_time',
           'stacked_versus_time',
           'excess_versus_E_and_time',
           'lightcurve']

#------------------------------------------------------------------------------
def excess_versus_time(dsets, rate=False, unmasked=False, ax=None,
                       xscale="log", yscale="log",
                       color1 = "tab:blue",color2 = "tab:orange",
                       debug  = False):
    """
    Plot both masked and unmasked (if not maked yet) count evolutions. Can 
    display rates if required.

    Parameters
    ----------
    dsets : list of Dataset objects
        The list of Dataset to be plotted.
    rate : Boolean, optional
        If True display the rate, otherwise the counts. The default is True.
    unmasked : boolean, optional
        True if the dataset were masked beforehand. The default is False.
    ax : matplotlib.axes, optional
        Current axis. The default is None.
    xscale : String, optional
        Abscissa scale. The default is "linear".
    yscale : String, optional
        y-scale. The default is "linear".
    color1 : TYPE, optional
        Color of masked data. The default is "tab:blue".
    color2 : String, optional
        Color of unmasked data. The default is "tab:orange".
    debug : TYPE, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    None.

    """

    ax = plt.gca() if ax is None else ax

    with quantity_support():

        t0 = dsets[0].gti.time_start[0]
        tmax = 0*u.d

        for i,ds in enumerate(dsets):

            time      = (ds.gti.time_start[0]-t0).sec*u.s + ds.gti.time_sum/2
            errtime   = 0.5*ds.gti.time_sum
            npred     = ds.npred_signal().data.sum()
            npred_msk = ds.npred_signal().data[ds.mask_safe].sum()
            tmax      = max(tmax,time)
            xs        = ds.excess.data.sum()
            xs_msk    = ds.excess.data[ds.mask_safe].sum()

            norm      = 1/2/errtime if rate else 1*u.dimensionless_unscaled

            if (debug):
                print("{:2} : t = [{:10.2f} {:10.2f}] - y={:8.2f} ym={:8.2f} / th={:8.2f} thm={:8.2f}"
                      .format(i,(time-errtime).value, time+errtime,
                              xs*norm,xs_msk*norm,npred*norm,npred_msk*norm))
                ax.text(time,
                        max(1.2*abs(xs_msk)*norm,10*norm),
                        str(i))

            # Handling of quantities by Errorbar is strange
            # - requires numpy arrays - quantities are not scalars !
            ax.errorbar(x = [time.value]*time.unit,
                        y = [xs_msk*norm.value]*norm.unit,
                        xerr = [errtime.value]*errtime.unit,
                        yerr = [np.sqrt(xs_msk)*norm.value]*norm.unit,
                        label="Excess counts",color=color1)
            ax.bar( time,xs_msk*norm,width=2*errtime,alpha=0.2,color=color1)

            # Before masking
            if (unmasked):
                ax.errorbar(x     = [time.value]*time.unit,
                            y     = [xs*norm.value]*norm.unit,
                            xerr  = [errtime.value]*errtime.unit,
                            yerr  = [np.sqrt(xs)*norm.value]*norm.unit,
                            label ="Total excess counts",color=color2)
                ax.bar(time,npred*norm,width=2*errtime,alpha=0.2,color=color2)

        ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
        if rate:
            ax.set_ylabel("Excess Count rate (" + ax.yaxis.get_label_text()+ ")")
        else:
            ax.set_ylabel("Excess counts")
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if tmax>1*u.d: ax.axvline(x=1*u.d,ls=":",
                                      color="grey",label="One day")

        import collections
        handles, labels = ax.get_legend_handles_labels()
        by_label = collections.OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    return
#------------------------------------------------------------------------------
def plot3D(x, y, data,
                 tag    = "Unknown",
                 ax3d   = None,
                 zscale = "linear",
                 cmap   = "plasma",
                 debug  = False):
    """
    Create a 3D plot with x is energy, y is time, z is data

    Parameters
    ----------
    x : 1D Numpy array or list
        x values
    y : 1D Numpy array or list
        y values
    data : 2D Numpy array
        The data
    tag : String, optional
        Used as a title. The default is "Unknown".
    ax3d : Matplotlib 3D axis, optional
        Axis used for plotting. The default is None.
    zscale : String, optional
        z axis scale, either "liner" or "log". The default is "linear".
    cmap : String, optional
        A colormap name. The default is "plasma".
    debug : Boolean, optional
        If True print some information. The default is False.

    Returns
    -------
    None.

    """

    if (ax3d == None):
        fig = plt.figure()
        ax3d = fig.add_subplot(111, projection='3d')

    yE, xt = np.meshgrid(y,x)

    zlbl = "Counts"
    if (zscale=="log"):
        data = np.log10(data)
        zlbl = zlbl + "(log)"

    ax3d.plot_surface(xt, yE, data,cmap=cmap, alpha=1.0,edgecolor='none')

    ax3d.set_title(tag)
    ax3d.set_ylabel('log(E) ('+str(y.unit)+')')
    ax3d.yaxis.labelpad=20
    ax3d.set_xlabel('Elapsed time ('+str(x.unit)+')')
    ax3d.xaxis.labelpad=20
    ax3d.set_zlabel(zlbl)
    ax3d.zaxis.labelpad=20
    ax3d.view_init(30, 45)

    return
#------------------------------------------------------------------------------
def excess_versus_E_and_time(dsets):

    """
    Excess versus E and time in 3D

    Parameters
    ----------
    dsets : Datasets object
        A collection of Gammapy datasets

    Returns
    -------
    None.

    """
    from dataset_tools import get_axis
    dt, E = get_axis(dsets, tunit=u.h, Eunit=u.GeV)

    data_pred   = np.asarray([ds.npred_signal().data.flatten() for ds in dsets])
    data_excess = np.asarray([ds.excess.data.flatten() for ds in dsets])

    fig = plt.figure(figsize=(20,8))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    plot3D(dt, E, data_pred,  tag="Prediction",ax3d = ax1)
    plot3D(dt, E, data_excess, tag="Excess counts", ax3d = ax2)

    plt.tight_layout()
    plt.show()

    return

#------------------------------------------------------------------------------
def lightcurve(dsets,
               t_unit = u.s, tmin=0*u.s, tmax=1000*u.s, binwidth=10*u.s,
               tag="excess",
               ax=None, xscale = "linear", yscale="linear",
               style="bar", marker="o",
               color="black", fillcolor=None,
               debug=False, **kwargs):
    """
    Reproduce (generate) a lightcurve from a list of dataset with defined 
    observation slices. Use interpolation along the time slices.

    Parameters
    ----------
    dsets : Dataset List
        The list od Dataset.
    t_unit : astropy.unit, optional
        The unit time of the plot. The default is u.s.
    tmin : astropy.time, optional
        Plot minimal time. The default is 0*u.s.
    tmax : astropy.time, optional
        Plot maximal time. The default is 1000*u.s.
    binwidth : astropy.time, optional
        Plot time bin width. The default is 10*u.s.
    tag : String, optional
        A tag defining the data to plot among excess and prediction.
        The default is "excess".
    ax : matplotlib.axes, optional
        Current plot axis. The default is None.
    xscale : String, optional
        Abscissa matplotlib scale. The default is "linear".
    yscale : String, optional
        Ordinate matplotlib scale. The default is "linear".
    style : String, optional
        Matplotlib style, either "line" or "bar". The default is "bar".
    marker : String, optional
        matplotlib marker style. The default is "o".
    color : String, optional
        matplotlib color. The default is "black".
    fillcolor : String, optional
        matplotlib fillcolor. The default is None.
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.
    **kwargs : Arbitrary Keyword Arguments
        Arbitrary Keyword Arguments
    Returns
    -------
    ax : matplotlib.axes
        Current matplotlib axis.

    """
    
    from scipy.interpolate import interp1d

    tmin     = tmin.to(t_unit)
    tmax     = tmax.to(t_unit)
    binwidth = binwidth.to(t_unit)

    # Create the time interpolation of the counting rate
    tsamp, dndt, tslice = [], [], []
    t0 = dsets[0].gti.time_start[0]
    duration = 0*t_unit

    for i,ds in enumerate(dsets):

        livetime = ds.gti.time_sum

        tsamp.append( ((ds.gti.time_start[0]-t0).sec*u.s
                           + livetime/2).to(t_unit) )
        tslice.append(duration)
        duration = duration + livetime.to(t_unit)

        if (tag == "excess"):
            dndt.append(ds.excess.data[ds.mask_safe].sum()/livetime.to(t_unit))
        elif (tag=="prediction"):
            dndt.append(ds.npred_signal().data[ds.mask_safe].sum()/livetime.to(t_unit))
        else:
            print("light_curve: tag not implemented")

    # From now on, all variables have values corresponding to t_unit
    tsamp  = np.asarray([x.value for x in tsamp])
    dndt   = np.asarray([x.value for x in dndt])
    tslice = np.asarray([x.value for x in tslice])

    tmin = tmin.value
    tmax = tmax.value
    binwidth = binwidth.value

    # Warning : f is not a quantity - Therefore all time units should be
    # the same before the interpolation - f should have unit 1/t_unit like dndt
    f = interp1d(tsamp,dndt,fill_value="extrapolate")
    #print(tsamp,dndt,f(tsamp))

    # Generate a number of counts for each bin
    nbin    = int((tmax-tmin)/binwidth)
    t_edges = np.linspace(tmin,tmax,nbin+1)
    t_bins  = t_edges[:-1] + 0.5*(t_edges[1:]-t_edges[:-1])
    # print(t_bins)

    # Compute count number in each bin = rate * bin width"
    n_random = []
    for t in t_bins:
        counts = f(t)*binwidth
        # print(t,"rate=",f(t),"counts=",counts)
        if (tag != "prediction"):
            # If excess is negative, do not fluxtuate
            if (counts>=0):
                n_random.append(np.random.poisson(int(counts)))
            else:
                print(" Negative counts : ",counts," not fluctuated")
                n_random.append(int(counts))
        else:
            n_random.append(counts)
    #print(n_random)

    ax = plt.gca() if ax is None else ax
    
    if (style=="bar"):
        ax.bar(t_bins, n_random, width=binwidth, edgecolor=color, color=fillcolor,label=tag, **kwargs)
    if (style=="line"):
        ax.plot(t_bins, n_random, color=color,marker=marker,lw=1,label=tag, **kwargs)
    if (debug):
        ax.plot(t_bins,f(t_bins)*binwidth, ls=":",color="red")
        ax.plot(tsamp[tsamp<tmax],dndt*binwidth,marker="o",ls="",color="green",label="Initial data")

    ax.set_xlabel("Observation duration ("+str(t_unit)+")")
    ax.set_ylabel("Counts per bin (" + str(binwidth) + str(t_unit)+")")
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)

    for xl in tslice[tslice<tmax]:
        ax.axvline(x=xl,ls="--",lw=1,alpha=0.5, color="grey")
    ax.axhline(y=0)

    ax.legend()

    return ax
#------------------------------------------------------------------------------
def stacked_versus_time(dstacked, dsets=None,
                        rate=True, ax=None, debug=False,
                        xscale="linear", yscale="linear",
                        color="tab:blue"):
    """
    DRAFT - Display stacked statictics and model statistics if original datasets is 
    given (since model stacking is irrelevant). 

    Parameters
    ----------
    dstacked : TYPE
        DESCRIPTION.
    dsets : TYPE, optional
        DESCRIPTION. The default is None.
    rate : TYPE, optional
        DESCRIPTION. The default is True.
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    debug : TYPE, optional
        DESCRIPTION. The default is False.
    xscale : TYPE, optional
        DESCRIPTION. The default is "linear".
    yscale : TYPE, optional
        DESCRIPTION. The default is "linear".
    color : TYPE, optional
        DESCRIPTION. The default is "tab:blue".

    Returns
    -------
    None.

    """

    print("to be finalised !")
    return

    if ax == None : fig,ax = plt.subplots(figsize=(7,7))

    duration = 0*u.s
    npred = 0
    with quantity_support():
        t0 = dstacked[0].gti.time_start[0]
        i = 0
        for dst, ds in zip(dstacked,dsets):

            livetime= dst.gti.time_delta[i]
            time       = (dst.gti.time_start[i]-t0).sec*u.s + livetime/2

            print(i," >>> ",livetime)
            errtime    = 0.5*(livetime-duration)
            i+=1

            xs     = dst.excess.data.sum()
            npred += ds.npred_signal().data[ds.mask_safe].sum()

            norm   = 1*u.dimensionless_unscaled if rate is False else 1/2/errtime

            if (debug):
                print(" At {:8.2f} +/-{:8.2f} : nxs = {:8.2f} +/-{:8.2f} Ratio={:5.2f}"
                      .format(time.value,errtime,xs*norm.value,np.sqrt(xs)*norm,npred/xs))

            # Handling of quantities by Errorbar is strange - requires numpy arrays - quantities are not scalars !

            ax.errorbar(x = [time.value]*time.unit,
                        y = [xs*norm.value]*norm.unit,
                        xerr = [errtime.value]*errtime.unit,
                        yerr = [np.sqrt(xs)*norm.value]*norm.unit,
                        label="Excess counts",color=color)
            ax.bar(time,npred*norm,width=2*errtime,alpha=0.2,color=color,label="Theory")
            duration += livetime

        ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
        ax.set_ylabel("Counts (" + ax.yaxis.get_label_text()+ ")")
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        import collections
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = collections.OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    return
#------------------------------------------------------------------------------
def excess_counts(ds, elapsed = 0*u.s,
                  ax   = None, stacked = False, model_bar = True,
                  n_min   = 0.1,      n_max  = None, 
                  emin    = 10*u.GeV, emax   = 20*u.TeV,
                  color   = "tab:orange",tag = None,
                  model_label = "$\mu_{src}$", alpha_model   = 0.2, 
                  lgd_col = 2, lgd_size = 12,
                  debug = False,
                  **kwargs):
    """
    This is a modified version of dataset plot_counts method.
    It does not explicitely mask the energy bins    

    Parameters
    ----------
    ds : Dataset
        A gammapy Dataset.
    elapsed : astropy.time, optional
        The time of start of this particular dataset. The default is 0*u.s.
    ax : matplotlib.axes, optional
        Current axis. The default is None.
    stacked : Boolean, optional
        Indicate if the current Dataset was previously stacked. 
        The default is False.
    model_bar : Boolean, optional
        If True, plot the model as bars. The default is True.
    n_min : float, optional
        Minimum number of counts on the plot. The default is 0.1.
    n_max : float, optional
        Maximal number of counts on the plot. The default is None.
    emin : astropy.Quantity, optional
        Minimum energy of the plot. The default is 10*u.GeV.
    emax : stropy.Quantity, optional
        Maximal energy of the plot. The default is 20*u.TeV.
    color : String, optional
        matplotlib color of the plot. The default is "tab:orange".
    tag : String, optional
        Extra information on top of the label. The default is None.
    model_label : String, optional
        The model plot label. The default is "$\mu_{src}$".
    alpha_model : TYPE, optional
        trabsparency of the model (if plotted). The default is 0.2.
    lgd_col : integer, optional
        Number of columns of the legend. The default is 2.
    lgd_size : float, optional
        Plot legned font size. The default is 12.
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.
    **kwargs : Arbitrary Keyword Arguments
        Arbitrary Keyword Arguments.

    Returns
    -------
    n_min : float
        Minimum counts - useful for the panelling
    max_counts : float
        Maximal counts - useful for the panelling.

    """
    
    ax = plt.gca() if ax is None else ax
    
    Eedges  = ds.background.geom.axes[0].edges.flatten() # Need the unit 
    
    ###-----------------------------------------
    ### Plot excess count using gammapy function
    ###-----------------------------------------
    
    # Define label
    label = "$t$ = " + t_str(elapsed) \
          + "\n$\Delta$t = "+t_str(ds.gti.time_sum)
    if tag != None: label = tag+"\n"+label
    
    # If data are stacked, no model are plotted - draw a line
    ls="-" if stacked else ""

    # Plot, retain the max. counts for the panel function
    ds.excess.plot(ax=ax, label=label,ls=ls,color=color, **kwargs)
    max_counts = np.nanmax(ds.counts.data).item()
    
    ###-----------------------------------------
    ### Plot theory and energy range if possible (not stacked)
    ###-----------------------------------------
    if (not stacked): 
        
        ### Get predicted counts -> Modify max. counts for later
        npred      = ds.npred_signal().data.flatten()
        max_counts = np.max([max_counts, np.nanmax(npred).item()])

        ### use default plotting - quite unreadable with data points
        if (not model_bar): # Use default function
            ds.npred_signal().plot(ax=ax, label=model_label)

            
        ### Use bar to have the model easily readable
        else: 
            # in log scale center is not center, widths are assymetric
            Ecenter = ds.background.geom.axes[0].center.flatten()
            width_r = Eedges[1:]-Ecenter
            width_l = Ecenter - Eedges[:-1]

            ax.bar(Ecenter, npred, width = -width_l,
                   alpha = alpha_model, align="edge", color=color, label=model_label)
            ax.bar(Ecenter, npred, width = width_r,
                   alpha = alpha_model, align="edge", color=color)
            
        ax.axvline(ds.energy_range[0],ls="--",color="grey")
        ax.axvline(ds.energy_range[1],ls="--",color="grey")
        # ds._plot_energy_range(ax=ax, label="") # Show the dataset energy masking
        

    ax.set_xlim(emin.to(Eedges[0].unit).value,emax.to(Eedges[0].unit).value)
    ax.grid("both",ls="--",alpha=0.5)

    ax.set_ylabel("Excess count number")
    ax.set_ylim(ymin=n_min)    
    if (n_max != None):  
        ax.set_ylim(ymax=n_max)
        max_counts = n_max
    else: n_max = max_counts
    
    ax.set_title("")

    # ax.legend(ncol=lgd_col,fontsize=lgd_size)
    from utilities import single_legend
    single_legend(ax, ncol=lgd_col,fontsize=lgd_size)

    return (n_min, max_counts)
#------------------------------------------------------------------------------
def residuals(ds,ax=None,**kwargs):
    """
    Plot residuals - Just use the gammapy function

    Parameters
    ----------
    ds : Dataset
        Current dataset.
    ax : matplotlib.axes, optional
        Current axis. The default is None.
    **kwargs : Arbitrary Keyword Arguments
        Arbitrary Keyword Arguments.

    Returns
    -------
    None.

    """
    
    ds.plot_residuals(ax=ax,**kwargs)
    return
#------------------------------------------------------------------------------

def comparison_plot(GRB_id,site="South", merging=True, 
                    arrays = ["omega","alpha"], debug=False):
    """
    Plot counts for several configurations (arrays).

    Parameters
    ----------
    GRB_id : String
        GRB identifier (xyz for Eventxyz.fits).
    site : String, optional
        "South", "North". The default is "South".
    merging : Boolean, optional
        Indicated if the dataset has been stacked. The default is True.
    arrays : List of String, optional
        The arrays or configuration names to be compared. 
        The default is ["omega","alpha"].
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    None.

    """
    
    import pickle
    from pathlib import Path
    import astropy.units as u
    from dataset_tools import createonoff_from_simulation, check_datasets, sigmax , compactify
    from gammapy.utils.random import get_random_state
    import yaml
    from yaml.loader import SafeLoader

    #plt.style.use('seaborn-talk') # Make the labels readable
    plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!       
    ###----------------------    
    ### Get data
    ###----------------------    
    # base = Path("D:/000_Today/SoHAPPy_tests/single_sources/")
    base = Path("D:/000_Today/SoHAPPy_tests/single_sources_denser/")
    
   
    # Read parameters
    par = (yaml.load(open("parameter.yaml"), Loader=SafeLoader))
    GRB = GRB_id+site[0] # The key to find back the parameters
    
    # Open Omega file
    
    array   = arrays[0]
    file    = Path(base,array,"Event"+GRB_id+"-"+site+"_sim.bin")
    infile  = open(file,"rb")
    mc_1    = pickle.load(infile)
    infile.close()
    # filetag_1 = "GRB "+GRB_id+" - South ($\\"+array +"$)"
    dset_init_1 = createonoff_from_simulation(mc_1,
                                              random_state=get_random_state(2021))

    # Open Alpha file
    array=arrays[1]
    file    = Path(base,array,"Event"+GRB_id+"-"+site+"_sim.bin")
    infile  = open(file,"rb")
    mc_2    = pickle.load(infile)
    infile.close()
    # filetag_2 = "GRB "+GRB_id+" - South ($\\"+array +"$)"
    dset_init_2 = createonoff_from_simulation(mc_2,
                                              random_state=get_random_state(2021))
    

    ###----------------------    
    ### Merge slices if requested 
    ###----------------------  
    if merging: # Stacking is also possible, Cf. dataset_tools
        dset_list_1 = compactify(dset_init_1,
                                 dtmin=u.Quantity(par[GRB]["dtmin"]))
        dset_list_2 = compactify(dset_init_2,
                                 dtmin=u.Quantity(par[GRB]["dtmin"]))
    else:
        dset_list_1 = dset_init_1.copy()
        dset_list_2 = dset_init_2.copy()
        
    ## Check datasets
    if debug:
        check_datasets(dset_init_1,deeper=True,masked=False)
        check_datasets(dset_list_1,deeper=True,masked=False)
        check_datasets(dset_init_2,deeper=True,masked=False)
        check_datasets(dset_list_2,deeper=True,masked=False)    
    print(" Max. Significance before merging")
    print("   - {:6s} = {:5.1f}".format(arrays[0],sigmax(dset_init_1,stacked=False)))
    print("   - {:6s} = {:5.1f}".format(arrays[1],sigmax(dset_init_2,stacked=False)))
    print(" Max. Significance after (maybe) merging")
    print("   - {:6s} = {:5.1f}"
          .format(arrays[0],sigmax(dset_list_1, stacked=merging)))
    print("   - {:6s} = {:5.1f}"
                    .format(arrays[1],sigmax(dset_list_2, stacked=merging)))
 
    nslots = len(dset_list_1)
    if nslots != len(dset_list_2):
        print(" Problem : the two sets are different")
        
    ###----------------------    
    ### Plot excess counts
    ###----------------------          
    fig, ax = plt.subplots(nrows=1, ncols=nslots, figsize=(nslots*5.5,6),sharey=True)

    # If model not shown, use steps for the darta points
    ds = "steps-mid" if merging else "default"
    
    iplot = 0
    
    emin  = 20*u.GeV
    emax  = 10*u.TeV
    lw = 2 if merging else 2
    
    for ds_1, ds_2, ax  in zip(dset_list_1, dset_list_2, ax):

        elapsed = ds_1.gti.time_start[0] - mc_1.slot.grb.t_trig
        excess_counts(ds_1, elapsed = elapsed, 
                      ax=ax, stacked = merging, model_bar = True,
                      n_min = 1., n_max  = float(par[GRB]["n_max"]), 
                      emin  = emin, emax = emax,
                      color ="tab:blue", alpha = 0.5, fillstyle="full",
                      model_label=None,
                      drawstyle=ds, lw=lw,
                      tag = arrays[0],
                      debug = False)    
        elapsed = ds_2.gti.time_start[0] - mc_2.slot.grb.t_trig
        excess_counts(ds_2, elapsed = elapsed, 
                      ax=ax, stacked = merging, model_bar = True,
                      n_min = 1., n_max = float(par[GRB]["n_max"]),
                      emin = emin, emax = emax,
                      color = "tab:orange", alpha = 0.5, fillstyle="full",
                      model_label=None,
                      drawstyle=ds, lw = lw,
                      tag=arrays[1],
                      debug=False)
        if (iplot): ax.set_ylabel(None)
        iplot+=1
    fig.suptitle(GRB, fontsize=16)
    plt.tight_layout(w_pad=0)  
    
###############################################################################
if __name__ == "__main__":
    """
    Code example 
    """
    # Transform warnings into errors - useful to find who is guilty !
    import warnings
    #warnings.filterwarnings('error')
    warnings.filterwarnings('ignore') 
    
    glist = [343, 980, 465, 676, 785]
    for gid in glist:
        print(" **************** ",gid," *****************")
        comparison_plot(str(gid), site="South", merging=True)
        plt.show()
 
    