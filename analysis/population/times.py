# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 17:11:21 2022

@author: Stolar
"""
import sys
codefolder = "../../"
sys.path.append(codefolder)

import numpy as np
import astropy.units as u

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import setup as stp
from niceplot import MyLabel, stamp, single_legend, col_size,  vals_legend
from niceprint import heading, t_fmt

###-----------------------------------------------------------------------------
def sigmax_after_delay(pop, sigmx_min=5, tmin=1*u.d):

    """
    Find and display mean maximal significnaces, above a threshold, reached 
    after a certain time following the explosion, and the fraction of time
    this threshold is reached.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    sigmx_min : float, optional
        Mean siginificance value to be considered. The default is 5.
    tmin : Astropy Time, optional
        Mean duration to be considered. The default is 1*u.d.

    Returns
    -------
    None.

    """

    heading("Max significance the day after")

    print(" Events above ",sigmx_min," (meanvalue), ",t_fmt(tmin)," after the trigger")

    for g,txt in zip([pop.g_n, pop.g_s, pop.g_b, pop.g_tot],
                     ["North","South","Both","All"]):
        
        mask = g.sigmx > sigmx_min
        late = g[mask][ g[mask].tmx > tmin.to(u.s).value]

        print(f"\n----- {txt:5s}: {len(late):4d}")

        print(f"{'Source':<10s}: {'Date':>10s}   {'sigma':>5s}   {'%':>5s}")

        for s, t, sig, d5s in zip(late.name,late.tmx*u.s, late.sigmx, late.d5s):

            print(f"{s:10s}: {t_fmt(t):>8.2f}   {sig:>5.1f}   {d5s:>5.2f}"
                  ,end="")

            if d5s >= pop.eff_lvl:
                print(" ***")
            elif sig >= 3:
                print(" *")
            else:
                print()
                
###-----------------------------------------------------------------------------
def time_to_5sigma_fullrange(pop, logx=True, nbin=25):
    """
    Display mean times to get 5 sigma siginificance in North and South.

    Parameters
    ----------
    pop : Pop instane
        Data from disk.
    logx : Boolean, optional
        If True x axis in Log otherwise linear. The default is True.
    nbin : integer, optional
        Histogram bin numbers. The default is 10.

    Returns
    -------
    None.

    """

    fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(12,6))

    var = np.log10(pop.g_n[pop.g_n.d5s>pop.eff_lvl].t5s/3600)
    n, bins, _ = ax1.hist(var,
                          bins=nbin, color=stp.col_n, alpha=0.3,
                          label = MyLabel(var,"North"))

    var = np.log10(pop.g_s[pop.g_s.d5s>pop.eff_lvl].t5s/3600)
    n, bins, _ = ax1.hist(var,
                          bins=nbin, color=stp.col_s, alpha=0.3,
                          label = MyLabel(var,"South"))

    ax1.set_xlabel("log Time in hours")
    ax1.axvline(np.log10(24),  ls="-", color="grey",label="1 day")
    ax1.axvline(np.log10(1),   ls="--",color="grey",label="1 h")
    ax1.axvline(np.log10(1/60),ls=":", color="grey",label="1'")

    ax1.legend()
    ax1.set_title("Mean Time to reach 5 $\sigma$")
    stamp(pop.tag,axis=fig,where="bottom")
    plt.tight_layout()
    
###-----------------------------------------------------------------------------
def time_to_detection(sub_pop, binw=1, yscale="log", title="generic", **kwargs):
    """
    Display the distribution of the 5 sigma and the mean maximal significance detection times.
    Shows the delay starting from the burst or the start of the detection.
    All times are in hours.
    The plots says how log ot is necessary to wait for such detections (median values)
    
    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    binw : float, optional
        Histogram bin width in hour(s). The default is 1.
    yscale : string, optional
        Log or Linear for y-scale. The default is "log".
    **kwargs : Dictionnary
        Extra parameters.

    Returns
    -------
    None.

    """
    
    varlist  = [sub_pop.t5s/3600, 
                sub_pop.tmx/3600] # Go to hours
    masks    = [(sub_pop.d5s>pop.eff_lvl), 
                (sub_pop.tmx>=0) & (sub_pop.sigmx>=5)]
    tags     = ["$5\sigma$",
                "$\sigma_{max} \geq 5$"]
         
    tobs = sub_pop.t1/3600 # That's the observation start 

    # Interesting but useless here
    if masks is None: 
         masks = np.ones(len(varlist)).astype(bool)
         
    # Define one hour binning from the min and max
    tmax = int(max([max(v[m]) for v,m in zip(varlist,masks)]))
    bins = range(0, tmax+ 1, binw)
    # print(" Max time :",tmax, " bins:",bins)
    
    fig, ax = plt.subplots(nrows=len(varlist), ncols=1,
                           figsize=(12,3.5*len(varlist)),sharex=True)
    
    first = True
    for ax0, tag, var, mask in zip(ax,tags,varlist,masks):
        
        n, bins,_ = ax0.hist(var[mask], bins=bins, alpha=0.5,
                             label=MyLabel(var[mask],label=tag+" (Burst)",stat="med"),**kwargs)
        if len(tobs) != 0:
            ax0.hist(var[mask]-tobs[mask],bins=bins,alpha=0.5,
                     label=MyLabel(var[mask]-tobs[mask],label=tag+" (Start)",stat="med"),**kwargs)
            
        ax0.set_yscale(yscale)
        
        ax0.legend()
        ax0.grid(which='both')
        minor_ticks= bins
        ax0.set_xticks(minor_ticks, minor=True)
        
        ax0.grid(which="major",ls=":")
        ax0.grid(which="minor",ls=":",alpha=0.5)
        ax0.set_ylabel("$h^{-1}$") 
        if first: 
            ax0.set_title(title)
            first = False

        axx = inset_axes(ax0, width="50%", height=1.2,loc="upper center")
        nxx,binsxx,_ = axx.hist(var[mask][var[mask]<1]*60,bins=range(0,60,1),alpha=0.5)
        axx.axvline(107/60,ls=":",label=" Swift + LST delays")
        axx.set_xlabel("$\Delta t$ (min)")
        axx.set_ylabel("$min^{-1}$") 
        axx.legend()

    ax0.set_xlabel("Mean detection time (h)") # Last one
    stamp(pop.tag,axis=fig)    
    plt.tight_layout()

###----------------------------------------------------------------------------------------------
def detection_fraction_versus_time(pop, sub_pop, title="unknown", 
                                   bins=50, density=True, 
                                   axs = None, **kwargs):
    
    varlist  = [np.log10(sub_pop.t5s), np.log10(sub_pop.tmx)] # Go to hours
    masklist = [(sub_pop.d5s>pop.eff_lvl), (sub_pop.tmx>=0) & (sub_pop.sigmx>=5)]
    taglist  = ["$5 \ \sigma$","$\sigma_{max}$"]
    if axs is None:
        fig, axs = plt.subplots(ncols=len(varlist), nrows=1,
                                figsize=(7*len(varlist), 5),sharey=True)
    else:
        if len(axs): 
            fig = axs[0].get_figure()
        else:
            fig = axs.get_figure()
            
    first = True
    for ax0, var, mask, tag in zip(axs, varlist, masklist, taglist):
        if mask.all() is None:
            mask = np.ones(len(var))        
        n, bins,_ = ax0.hist(var[mask],bins,
                             histtype="step",lw=2,
                             cumulative=True,density=density, **kwargs)
       
        for t0, col in zip([107, 600, 3600, 3600*24],
                           ["grey","tab:green","tab:orange","tab:blue"]):
            ax0.axvline(np.log10(t0), label=str(t_fmt(t0*u.s)), color=col,ls="--",alpha=0.5)
            
        if first:
            ax0.set_ylabel("Fraction") 
            first = False
        ax0.grid(which='both')
        ax0.grid(which="major",ls=":")
        ax0.set_xlabel("$log_{10} \ t$")
        ax0.set_title(title+" - "+tag)
    
    single_legend(ax0,bbox_to_anchor=[1,1])
    stamp(pop.tag,axis=fig,where="left")
    plt.tight_layout()
    
    return axs

###----------------------------------------------------------------------------------------------
def plot_t_sigma(pop,sub_pop,title="unknown"):
    """
    Display mean maximal significance as the fucntion of the mean time to reach the maximal
    significance.

    Parameters
    ----------
    pop : Pop instance
        Data on disk
    sub_pop : pandas frame
        Sub poulation.
    tag : String, optional
        Data descriptor. The default is "".

    Returns
    -------
    None.

    """
    
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(12,6),sharex=True)

    mask = sub_pop.sigmx>0
    colors, sizes = col_size(sub_pop[mask].sigmx)
    ax.scatter(np.log10(sub_pop[mask].tmx),np.log10(sub_pop[mask].sigmx),
               marker="o",color=colors,alpha=0.7,s=sizes)
    ax.scatter(np.log10(sub_pop[mask].tmx),np.log10(sub_pop[mask].sigmx),
               marker="x",alpha=0.5,s=5,color="grey")

    ax.set_xlabel("$log_{10}(\Delta t_{max})$")
    ax.set_ylabel("$log_{10}(\sigma_{max})$")
    ax.axvline(x=np.log10(3600),ls=":",color="green",label="1 hour")
    ax.axvline(x=np.log10(12*3600),ls=":",color="blue",label="12 hour")
    ax.axvline(x=np.log10(24*3600),ls=":",color="red",label="One day")
    ax.axhline(np.log10(3),ls=":",color="black",lw="2",label="$3\sigma$")
    ax.axhline(np.log10(5),color="black",lw="2",label="$5\sigma$")
    
    ax.legend(loc="lower left")
    ax.set_title(title)
    
    patches = vals_legend(ax)
    fig.legend(title="$\sigma_{max}$",handles=patches, bbox_to_anchor=(1.1, 0.98))
    stamp(pop.tag,axis=fig,where="right")    
    plt.tight_layout()

###----------------------------------------------------------------------------------------------
def sig_vs_time_errors(pop, mask, ax=None, xscale="log", yscale="log",title=""):
    """
    Display the mean maximal significnace versus the mean time to have it reached,
    with error bars

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    mask : Boolean sequence
        Data subselection.
    ax : Matplotlib axes, optional
        Current axis. The default is None.
    xscale : String, optional
        "Log" or "linear". The default is "log".
    yscale : String, optional
        "Log" or "linear". The default is "log".
    title : String, optional
        Descriptor. The default is "".

    Returns
    -------
    None.

    """
    
    ax = plt.gca() if ax is None else ax
    ax.errorbar(pop[mask].tmx/3600,
                pop[mask].sigmx,
                xerr=pop[mask].etmx/3600,
                yerr=pop[mask].esigmx,
                ls="", marker="o",ecolor="tab:blue",color="red",alpha=0.5,
                label=MyLabel(pop[mask].sigmx,stat=None))
    # ax.axhline(y=3,ls=":",label="$3\sigma$",color="lightgrey")
    # ax.axhline(y=5,ls=":",label="$5\sigma$",color="lightgrey")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel("Mean Time of max detection (h)")
    ax.set_ylabel("Mean max sigificance")
    ax.set_title(title)
    ax.legend()
    
###----------------------------------------------------------------------------------------------
def delay_correlation_both_sites(pop, sigmin=5):
    """
    Find N and S detection for all both detected GRB
    Get name of GRB detected on both sites and retrieve the independent N and S detections

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    sigmin: float
        Minimal mean maximal significance in the data. Default is 5.
        
    Returns
    -------
    None.

    """
    
    ### Get North and South population seen onboth site - reindex the frames
    gbs = pop.g_s[(pop.g_s.name.isin(pop.g_b.name))]
    gbn = pop.g_n[(pop.g_n.name.isin(pop.g_b.name))]
    gb = pop.g_b    
    
    gbs.reset_index(inplace=True)
    gbn.reset_index(inplace=True)
    gb.reset_index(inplace=True)    
    
    # Apply cuts
    gbs = gbs[(gb.sigmx >= sigmin)]
    gbn = gbn[(gb.sigmx >= sigmin)]
    gb  = gb[(gb.sigmx >= sigmin)]    
    
    # Now that all frames have same index and size, do plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,7))
    color, size = col_size(gb.sigmx)
    
    ax.scatter(gbn.t1/3600,gbs.t1/3600,alpha=0.5,c= color, s=size, label = MyLabel(gbn.t1,"$\sigma_{max}$"))
    ax.set_title("Seen on both sites - $\sigma_{max} \geq "+str(sigmin)+"$")
    ax.set_xlabel("Mean detection time (h) - North")
    ax.set_ylabel("Mean detection time (h) - South")
    xmax = 26
    ymax = 26
    ax.set_xlim(xmax=xmax)
    ax.set_ylim(ymax=ymax)
    
    # ZOOM
    xmin_z = 0
    xmax_z = 2
    ymin_z = 0
    ymax_z = 5
    ax.axvspan(xmin_z-1, xmax_z, ymin=ymin_z-1, ymax=(ymax_z+1)/ymax, color="grey",ls=":",alpha=0.5,lw=2,ec="black",fill=False)
    
    axx = inset_axes(ax, width="30%", height=3,loc="upper right")
    axx.scatter(gbn.t1/3600, gbs.t1/3600,  alpha=0.5,c= color, s=size, label = MyLabel(gbn.t1,"$\sigma_{max}$"))
    axx.set_xlim(xmax=xmax_z)
    axx.set_ylim(ymax=ymax_z)
    
    stamp(pop.tag, axis=fig,where="bottom")
    patches = vals_legend(ax)
    fig.legend(title="$\sigma_{max}$",handles=patches,bbox_to_anchor=(0.95, 0.88))
    
################################################################################################
if __name__ == "__main__":

    # A standalone function to read a GRB and make various tests

    from population import Pop
    from pop_io import create_csv

    import sys
    codefolder = "../../"
    sys.path.append(codefolder)

    nyears, file, _ = create_csv(file="parameter.yaml",debug=True)
    pop = Pop(filename=file, nyrs= nyears)

    ### ------------------------
    ### Print some statistics
    ### ------------------------
    # Print max significance reached after some time, above 3/5 sigma
    sigmax_after_delay(pop, sigmx_min=5, tmin=1*u.d)
    
    ### ------------------------
    ### Plots
    ### ------------------------

    # Mean Time distribution to get 5 sigma in N and S
    time_to_5sigma_fullrange(pop)
    
    # Time to get 5 sigma or sigmax in various conditions
    time_to_detection(pop.g_tot,title="All")
    time_to_detection(pop.g_n,title="North")
    time_to_detection(pop.g_s,title="South")
    
    # Fraction of 5 sigma detection as a function fo time
    detection_fraction_versus_time(pop, pop.g_tot,title="All")    
    detection_fraction_versus_time(pop, pop.g_n,title="North")    
    detection_fraction_versus_time(pop, pop.g_s,title="South")
    
    # Scatter plot significance versus time to reach it
    plot_t_sigma(pop,pop.g_tot,title="All")   
    
    # Scatter plot of sigma versus time with errors
    fig,(ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(15,8))
    
    title = "All - " + "$\sigma_{max} \geq 20$"
    sig_vs_time_errors(pop.g_tot,(pop.g_tot.sigmx>=20), yscale="log", ax=ax1, title=title)
    title = "All - " + "$5 < \sigma_{max} < 20$"
    sig_vs_time_errors(pop.g_tot,(pop.g_tot.sigmx>=5) & (pop.g_tot.sigmx<20),
                       xscale="linear", yscale="linear", ax=ax2, title=title)
    stamp(pop.tag,axis=fig)
    plt.tight_layout()
    
    # Time differences between N and S for sources detected on both sites
    delay_correlation_both_sites(pop)
   