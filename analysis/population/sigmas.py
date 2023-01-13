# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 17:43:24 2022

@author: Stolar
"""
import sys
codefolder = "../../"
sys.path.append(codefolder)

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import setup as stp
from niceprint import heading
from niceplot import MyLabel, stamp, single_legend

plt.style.use('seaborn-talk') # Make the labels readable
#plt.style.use('seaborn-poster')

###-------------------------------------------------------------
def detection_level(pop, glist, bins=25, **kwargs):
    """
    Display the detection level for a list of subpopulation.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    gname : list of String
        Subpopulation names.
    nbin : integer, optional
        Histogram bin number. The default is 25.
    **kwargs : Dictionnary
        Extra arguments.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(nrows=1, ncols=len(glist), 
                           figsize=(6*len(glist),5),
                           sharey = True)
    first = True
    for gname, axi  in zip(glist,ax): 
        
        subpop = pop.__dict__[gname]
        
        for sigma, var, col in zip([3,5], 
                                   [subpop.d3s, subpop.d5s],
                                   [stp.col_3, stp.col_5]):
            
            var = 100*var/pop.niter
            label = str(sigma)+r"$\sigma$"

            axi.hist(var, bins, range=[0,100], 
                     label=MyLabel(var,label), color=col, alpha=0.5, **kwargs)
            
            axi.set_xlabel("Confidence level (%)")
            axi.set_yscale("log")
            axi.set_title(gname)
            axi.axvline(x=pop.eff_lvl,
                       color="red",ls=":",label="min. C.L. (%)")
            axi.grid("both")
            
            if first:
                axi.set_ylabel("Event counts")
                first = False
            
        single_legend(axi)

    plt.tight_layout()

###-----------------------------------------------------------------------------
def plot_sigmax_all(pop, nbin=25, logx=True):
    """
    Plot the mean maximal significance distributions for North and North only, 
    South and South only, and GRBs detected on both site. Check that detections
    with both sites give higher significances.
    
    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    nbin : TYPE, optional
        Number of bins in histograms. The default is 25.
        
    Returns
    -------
    None.
    """

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3,
                                        figsize=(18,6),sharey=True)
    
    # North and North only
    ax, bins = plot_sigmax(pop.g_n,tag="North",
                           ax=ax1,bins=nbin,logx=logx,color=stp.col_n,alpha=0.5)
    plot_sigmax(pop.g_n0,tag="North only",
                ax=ax1,logx=logx,bins=bins,color=stp.col_n,alpha=0.5)
    
    # South and South only
    ax, bins = plot_sigmax(pop.g_s,tag="South",
                        ax=ax2,bins=nbin,logx=logx,color=stp.col_s,alpha=0.5)
    plot_sigmax(pop.g_s0,tag="South only",
                ax=ax,logx=logx,bins=bins,color=stp.col_s,alpha=0.5)
    
    # Both, all
    ax, bins = plot_sigmax(pop.g_b,tag="Both",
                           logx=logx,color=stp.col_b,ax=ax3,bins=nbin,alpha=0.5)
    plot_sigmax(pop.g_tot,tag="N+S",
                logx=logx,facecolor="none",edgecolor="black",ax=ax3,bins=bins,alpha=1,
                weight=len(pop.g_b)/len(pop.g_tot))
    
    stamp(pop.tag,axis=fig,where="bottom")
    plt.tight_layout()    

###----------------------------------------------------------------------------
def plot_sigmax(g, 
                ax=None,logx=True,bins=25,tag="",weight=1,xmax=None,**kwargs):
    """
    Plot mean maximal significance for a given population

    Parameters
    ----------
    g : Pandas table
        Current subpopulation.
    ax : Matplotlib axes, optional
        Current axis. The default is None.
    logx : Boolean, optional
        Abscissa in log scale if True. The default is True.
    bins : integer, optional
        number of bisn in histrogram. The default is 25.
    tag : string, optional
        Descriptor of the data. The default is "".
    weight : float, optional
        Histrogram reweighting. The default is 1.
    xmax : Float, optional
        Abscissa maximal value. The default is None.
    **kwargs : Dictionnary
        Extra matplotlib arguments.

    Returns
    -------
    ax : Matplotlib axes
        Current axis.
    bins : TYPE
        DESCRIPTION.

    """
    if (ax==None): fig,ax = plt.subplots()
        
    if (xmax != None): mask = g.sigmx < xmax
    else: mask= (np.ones(len(g.sigmx),dtype=bool)) & (g.sigmx!= -1)
        
    if logx:
        mask = (g.sigmx>0) & (mask) 
        x = np.log10(g[mask].sigmx)
        
        ax.set_title("Positive values only")
        ax.set_xlabel("$log_{10}(\sigma_{max})$")   
        ax.axvline(x=np.log10(3),color="tab:orange",ls=":",label="$3\sigma$")
        ax.axvline(x=np.log10(5),color="tab:green",ls=":",label="$5\sigma$") 
    else:
        x = g[mask].sigmx
        ax.set_xlabel("$\sigma_{max}$")
        ax.axvline(x=3.,color="tab:orange",ls=":",label="$3\sigma$")
        ax.axvline(x=5,color="tab:green",ls=":"  ,label="$5\sigma$")   
        
    n, bins,_ = ax.hist(x, bins=bins, label=MyLabel(x,tag,stat="med"),
                        weights=np.ones(len(x))*weight,**kwargs) 
    
    single_legend(ax)
    return (ax, bins)        

###----------------------------------------------------------------------------
def high_sigma(gpop, ax=None, inset=True, sigmin= 5, 
               sig2show=100, 
               tag="", weight=1, **kwargs):
    """
    Display the high values of the mena maximal significance and tag the 
    highest value with the event number.

    Parameters
    ----------
    gpop : Pandas table
        Current subpopulation.
    ax : Matplotlib axis, optional
        Current axis. The default is None.
    inset : Boolean, optional
        If True, add a log-scale display in an inset box. The default is True.
    sigmin : float, optional
        Minimal significance to consider. The default is 5.
    sig2show : float, optional
        Minimal significance to be tagged. They willnot be shown if negative. 
        The default is 100.
    tag : String, optional
        Data descriptor. The default is "".
    **kwargs : Dictionnary
        Extra parameters.

    Returns
    -------
    ax : Matplotlib axes
        Current axis.
    bins : sequence
        Histogram binning.

    """

    var = gpop.sigmx[gpop.sigmx>=sigmin]
    weights = np.ones(len(var))*weight
    
    n,bins,_ = ax.hist(var,label=MyLabel(var*weights,tag,stat=None), 
                       weights = weights,**kwargs)
    ax.set_xlabel("$\sigma_{{max}} \geq {}$".format(sigmin))
    ax.set_ylabel("Event rate")

    # Inset in log
    if inset:
        axx = inset_axes(ax, width="75%", height=1.2,loc="upper right")
        nxx,binsxx,_ = axx.hist(np.log10(gpop[gpop.sigmx>0].sigmx),
                           bins=25,color="grey",edgecolor="black", alpha=0.5)
        axx.axvline(np.log10(3),ls=":",color="red",lw="2",label="$3\sigma$")
        axx.axvline(np.log10(5),color="red",lw="2",label="$5\sigma$")
        axx.set_xlabel("Log $\sigma_{max}$")
        axx.legend()
        
        ax.legend(loc="upper left")
    else:
        ax.legend()

    # Outliers
    if sig2show>0:
        outlier=zip(gpop[gpop.sigmx>sig2show].name,
                    gpop[gpop.sigmx>sig2show].sigmx)
        i=1
        for g in outlier:
            y = (0+ (1.5*i + 0.1))*weight
            ax.text(x=g[1],y=y,s=g[0][5:],rotation=45)
            i+=1
            
    return ax, bins

###----------------------------------------------------------------------------
def plot_high_sigma_all(pop, nbin=100, sigmin=5, inset=False):
    """
    Plot the high value tail of the mean maximal sigificance values and
    tag the data with the even number

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    nbin: integer
        Histrogram number of bins. Default is 100
    sigmin: float
        Minimal significance values in the plots. The default is 5
    inset : boolean, optional
        If True, display the same plot in log-scale in an inset box. 
        The default is False.

    Returns
    -------
    None.

    """
    weight = 1/pop.nyears
    
    fig,ax = plt.subplots(nrows=3, ncols=1,figsize=(10,9))
        
    # North and North only
    ax0, bins = high_sigma(pop.g_n, sigmin = sigmin,
                           bins = nbin, color=stp.col_n,ax=ax[0],
                           tag="N", alpha=0.5, weight=weight, inset=inset)
    high_sigma(pop.g_n0,  sigmin=sigmin,
               bins = bins, color=stp.col_n,ax=ax0,tag="N only", 
               weight=weight,inset=inset)
    ax0.set_xlabel(None)
    ax0.set_title(" Mean $\sigma_{max}$ for "+str(pop.nyears)+" years of observation")
    ax0.grid("both",ls=":")
    
    # South and South only
    ax0, bins = high_sigma(pop.g_s, sigmin=sigmin,
                           bins = nbin, color=stp.col_s,ax=ax[1],tag="S", 
                           alpha=0.5, weight=weight,inset=inset)
    high_sigma(pop.g_s0, sigmin=sigmin,
               bins = bins, color=stp.col_s,ax=ax0,tag="S only",
               weight=weight,inset=inset)
    ax0.set_xlabel(None)
    
    # Both sites
    ax0, bins = high_sigma(pop.g_b, sigmin=sigmin,
                           bins=nbin, color=stp.col_b,ax=ax[2],
                           tag="Both", alpha=0.5,weight=weight, inset=inset)
    ax0.grid("both",ls=":")
    
    stamp(pop.tag,axis=fig,where="bottom")
    plt.tight_layout()   

###-----------------------------------------------------------------------------
def sigmax_above_threshold(pop, sigmin=0):
    """
    Print event number with mean maximal significance above a given threshold.

    Parameters
    ----------
    pop : Pop instance
        Data from disk.
    sigmin : float, optional
        Signiifcnace threshold. The default is 0.

    Returns
    -------
    None.

    """
    heading(f"Source with mean maximal sigma above {sigmin:5.1f}")

    poplist = [pop.g_n,pop.g_n0,pop.g_s,pop.g_s0,pop.g_b]
    taglist = ["North","N only","South","S only","Both"]  
    
    # Mean sigmax value of events above 5 sigma at 90%CL and sigmin
    for gsub, tag in zip(poplist, taglist):
        print("---- ",tag)
        gsub = gsub[gsub.d5s>=pop.eff_lvl]
        for sig, name in zip(gsub.sigmx,gsub["name"]):
            if sig>= sigmin:
                print(" - ",name," :",sig)
                
###----------------------------------------------------------------------------
def rate_above_threshold(pop, duration=None):
    """
    Print detection rates for various mean maximal significance thresholds.

    Parameters
    ----------
    pop : TYPE
        DESCRIPTION.
    duration : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    if duration is None: duration = pop.nyears
    
    heading(f"Detection rate with mean sigmax above sig for {duration:} years")

    poplist = [pop.g_n,pop.g_n0,pop.g_s,pop.g_s0,pop.g_b]
    taglist = ["North","N only","South","S only","Both"]  
    siglist = [-0.5, 3, 5, 10, 20, 50, 100, 200]
    
    print("{:>10s}".format("Sig min"),end="")
    for tag in taglist:
        print("{:>8s}".format(tag),end="")
    print("{:>8s}".format("Tot"))
    print(60*"-")
    
    for sig in siglist:
        print("{:>10.1f}".format(sig),end="")
        ntot = len(pop.g_tot[pop.g_tot.sigmx>=sig])
        for p  in poplist:
            n = len(p[p.sigmx>=sig])
            print("{:>8.2f}".format(duration*n/pop.nyears),end="")
        print("{:>8.2f}".format(duration*ntot/pop.nyears))

    print(60*"-")
    
###----------------------------------------------------------------------------
def plot_sigmax_cumulative(pop, bins=10):
    """
    Show the detection rates as a function of the mean maximal significance 
    threshold.

    Parameters
    ----------
    pop : Pop instance
        Data on disk.
    nbin : TYPE, optional
        DESCRIPTION. The default is 10.

    Returns
    -------
    None.

    """
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9,6))

    for g, tag, w in zip([pop.g_tot, 
                          pop.g_b, 
                          pd.concat([pop.g_n,pop.g_s],axis=0)],
                          ["All", "Both", "0.5*(N+S)"],
                          [1, 1, 0.5]):
        
        mask = g.sigmx > 0
        var = np.log10(g[mask].sigmx)
        
        n, bins = np.histogram(var, bins, weights=w*np.ones(len(var)))
        
        ncum = np.cumsum(n)
        ax.plot(bins[:-1] + 0.5*(bins[1:]-bins[:-1]),
                (max(ncum)-ncum)*10/pop.nyears,alpha=1,
                marker="o",label=tag)
    

        
        # Display some sigma values
        ax.axvline(np.log10(50),color="grey",alpha=0.5,ls="--")
        ax.text(np.log10(50)*1.05,0.2,"$50\sigma$",size=12)
        
        ax.axvline(np.log10(20),color="grey",alpha=0.5, ls="--")
        ax.text(np.log10(20)*1.05,0.2,"$20\sigma$",size=12)
                
        ax.axvline(np.log10(5),color="grey",alpha=0.5,ls=":")
        ax.text(np.log10(5)*1.05,0.2,"$5\sigma$",size=12)
    
        ax.set_xlabel("$Log_{10}$ $\sigma_{max}$")
        ax.set_ylabel("Rate above $\sigma_{max}$ for 10 yr of operation")
        ax.set_yscale("log")
        single_legend(ax,bbox_to_anchor=[1.0,1.0])
        ax.grid("both",ls=":")       
        
    stamp(pop.tag, axis=fig, where="bottom")
    
###############################################################################
if __name__ == "__main__":

    from population import Pop
    from pop_io import create_csv

    import sys
    codefolder = "../../"
    sys.path.append(codefolder)

    nyears, file, _ = create_csv(file="parameter.yaml",debug=True)
    
    print(nyears, file)
    pop = Pop(filename=file, nyrs= nyears)

    ### ------------------    
    ### Print statistics 
    ### ------------------    
    
    # Detection rate versus minimal significance normalised to duration
    rate_above_threshold(pop, duration=1)
    
    # Mean maximal siginificance above threshold
    sigmax_above_threshold(pop, sigmin=100)
    
    ### ------------------    
    ### Plots  
    ### ------------------    
    
    # Detection level distributions
    detection_level(pop ,["g_tot","g_n", "g_s"])

    # Mean max significance distribution - compare N/Nonly etc.
    plot_sigmax_all(pop)
    
    # Highest significance values
    plot_high_sigma_all(pop, inset=False)
    
    # Rate versus mean maximal siginificance threshold
    plot_sigmax_cumulative(pop)

