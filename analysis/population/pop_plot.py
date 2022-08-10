# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 11:58:24 2021

@author: Stolar
"""
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from utilities import MyLabel, single_legend, stamp
from setup import dtpop


history = {"190114C":
                {"Observatory": "MAGIC",    
                 "z": 0.4245, 
                 "Eiso": 3.00E+53*u.erg, 
                 "t90": 25*u.s,
                 "marker":"^",
                 "col":"black"},
               
           "201216C":
                {"Observatory": "MAGIC",    
                 "z": 1.1, 
                 "Eiso": 5.00E+53*u.erg, 
                 "t90": 30*u.s,
                 "marker":">",
                 "col":"black"},               
               
           "180720B":
                {"Observatory": "H.E.S.S.", 
                 "z":  0.654, 
                 "Eiso": 6.00E+53*u.erg, 
                 "t90": 49*u.s,
                 "marker":"v",
                 "col":"red"},
               
           "190829A":
                {"Observatory": "H.E.S.S.", 
                 "z": 0.0785, 
                 "Eiso": 2.00E+50*u.erg, 
                 "t90": 63*u.s,
                 "marker":"<",
                 "col":"red"},
               
           "080916C":
                {"Observatory":"Fermi/LAT", 
                 "z" :4.3,    
                 "Eiso": 8.80E+54*u.erg, 
                 "t90": -1*u.s,
                 "marker":"*",
                 "col":"tab:green"},
           
            "090902B":
                {"Observatory":"Fermi/LAT", 
                 "z": 1.822,  
                 "Eiso": 2.20E+52*u.erg, 
                 "t90": -1*u.s,
                 "marker":"*"  ,
                 "col":"tab:green"},
            
            "130427A":
                {"Observatory":"Fermi/LAT", 
                 "z": 0.34,   
                 "Eiso": 9.60E+53*u.erg, 
                 "t90": -1*u.s,
                "marker":"*" ,
                 "col":"tab:green"}           
          }
###-----------------------------------------------------------------------------
def historical():
    return history

###-----------------------------------------------------------------------------
def plot_historical(ax, ddict, obs=[]):
    """
    Add z, Eiso from already detected GRBs 

    Parameters
    ----------
    ax : matplotlib axis
        Current axis
    ddict : Dictionnay
        A dictionnary with predefined values 
    obs : String, optional
        Observatories to be shown from the dictionnary. The default is [].

    Returns
    -------
    ax : matpotlib axis
        Current axis

    """    
           
    ax = plt.gca() if ax is None else ax
    
    for target, value in ddict.items():
        data = ddict[target]
        if data["Observatory"] in obs:
            print(target)
            ax.scatter(data["z"],
                       np.log10(data["Eiso"].value),
                       marker = data["marker"],
                       color  = data["col"],
                       label  = target)
    return ax

###-----------------------------------------------------------------------------
def redshift(g, gref, g0=None, tag="unknown", reference=True,
             color="grey", color2="tab:green", 
             ax=None, zmax=4, det_level=90, nbin=20, w=1, **kwargs):

    ax = plt.gca() if ax is None else ax

    
    if reference == True:
    # First plot the reference "1000" GRB population
        mask = (gref.z<=zmax)
    
        ntot, bins, _  = ax.hist(gref[(gref.z<=zmax)].z, bins=nbin, 
                                 weights= np.ones(len(gref[(gref.z<=zmax)].z))*10/dtpop,
                                 facecolor="none",edgecolor="black")
    else:
        bins = nbin
        ntot = 0
        
    # Plot the first sub-population - apply weight for summed sample
    mask = (g.d5s>det_level) & (g.z<=zmax)
    n, bins, _  = ax.hist(g[mask].z,bins=bins, weights= np.ones(len(g[mask].z))*w*10/dtpop,
                          color=color,label=MyLabel(g[mask].z,label=tag),**kwargs)


    # Plot the second sub-population  - apply weight for summed sample
    if g0 != None:
        mask = (g0.d5s>det_level) & (g0.z<=zmax)
        n0, bins, _    = ax.hist(g0[mask].z,bins=bins,weights= np.ones(len(g0[mask].z))*w*10/dtpop,
                                 color=color,*kwargs)

    # Decoration
    ax.legend()
    ax.set_xlabel("Redshift (z)")
    ax.set_yscale("log")

    # Plot raio
    ratio = [ (x/xtot if xtot!=0 else 0) for x, xtot in zip(n, ntot) ]
    axx=ax.twinx()
    axx.plot(bins[:-1] + 0.5*(bins[1:]-bins[:-1]),ratio,alpha=1,color=color2,label="Ratio")
    axx.grid(ls=":")
    axx.legend(loc="lower right")
#         if tag != "Combined":

    return ax, axx
###-----------------------------------------------------------------------------
def eiso(g, gref, g0=None, tag="unknown", color="grey", color2="tab:green",ax=None, 
         Eisomin=5.5e50, Eisomax=5.5e55,  det_level=90, nbin=20, w=1, **kwargs):

    ax = plt.gca() if ax is None else ax

    # First plot the reference "1000" GRB population
    mask = (gref.Eiso<=Eisomax) & (gref.Eiso>=Eisomin)
    ntot, bins, _  = ax.hist(np.log10(gref[mask].Eiso), bins=nbin,weights= np.ones(len(np.log10(gref[mask].Eiso)))*10/dtpop,
                             facecolor="none",edgecolor="black")

    # Plot the first sub-population - apply weight for summed sample
    mask = (g.d5s>det_level) 
    n, bins, _     = ax.hist(np.log10(g[mask].Eiso),bins=bins, weights= np.ones(len(g[mask].Eiso))*w*10/dtpop,
                             color=color,label=MyLabel(np.log10(g[mask].Eiso),label=tag),**kwargs)


    # Plot the second sub-population  - apply weight for summed sample
    if g0 != None:
        mask = (g0.d5s>det_level) 
        n0, bins, _    = ax.hist(np.log10(g0[mask].Eiso),bins=bins,weights= np.ones(len(g0[mask].Eiso))*w*10/dtpop,
                                 color=color,**kwargs)

    # Decoration
    ax.legend(loc="upper left")
    ax.set_xlabel("$log_{10} \  E_{iso}$")
    ax.set_yscale("log")
    
    # Plot raio
    ratio = [ (x/xtot if xtot!=0 else 0) for x, xtot in zip(n, ntot) ]
    axx=ax.twinx()
    axx.plot(bins[:-1] + 0.5*(bins[1:]-bins[:-1]),ratio,alpha=1,color=color2,label="Ratio")
    axx.grid(ls=":")
    axx.legend(loc="lower right")
#         if tag != "Combined":

    return ax, axx
###-----------------------------------------------------------------------------
def epeak(g, gref, g0=None, tag="unknown", color="grey", color2="tab:green",ax=None, 
         Epeakmin=10, Epeakmax=100000,  det_level=90, nbin=20, w=1, **kwargs):

    ax = plt.gca() if ax is None else ax

    # First plot the reference "1000" GRB population
    mask = (gref.Epeak<=Epeakmax) & (gref.Epeak>=Epeakmin)
    ntot, bins, _  = ax.hist(np.log10(gref[mask].Epeak), bins=nbin,weights= np.ones(len(np.log10(gref[mask].Epeak)))*10/dtpop,
                             facecolor="none",edgecolor="black")

    # Plot the first sub-population - apply weight for summed sample
    mask = (g.d5s>det_level) 
    n, bins, _     = ax.hist(np.log10(g[mask].Epeak),bins=bins, weights= np.ones(len(g[mask].Epeak))*w*10/dtpop,
                             color=color,label=MyLabel(np.log10(g[mask].Epeak),label=tag),**kwargs)


    # Plot the second sub-population  - apply weight for summed sample
    if g0 != None:
        mask = (g0.d5s>det_level) 
        n0, bins, _    = ax.hist(np.log10(g0[mask].Epeak),bins=bins,weights= np.ones(len(g0[mask].Epeak))*w*10/dtpop,
                                 color=color,**kwargs)

    # Decoration
    ax.legend(loc="upper left")
    ax.set_xlabel("$log_{10} \  E_{peak} (keV)$")
    ax.set_yscale("log")
    
    # Plot raio
    ratio = [ (x/xtot if xtot!=0 else 0) for x, xtot in zip(n, ntot) ]
    axx=ax.twinx()
    axx.plot(bins[:-1] + 0.5*(bins[1:]-bins[:-1]),ratio,alpha=1,color=color2,label="Ratio")
    axx.grid(ls=":")
    axx.legend(loc="lower right")
#         if tag != "Combined":

    return ax, axx
###-----------------------------------------------------------------------------
def sig_plot_cumulative(var, mask=None, title=None,
                        ax=None, xlabel="",nbins=50, density=True,
                        **kwargs):
    
    if mask.all() == None: 
        mask = np.ones(len(var))
    
    ax = plt.gca() if ax is None else ax
    
    n, bins,_ = ax.hist(var[mask],bins=nbins,
                         histtype="step",lw=2,
                         cumulative=True,density=density, **kwargs)

    ax.axvline(np.log10(107),    label  ="Ref. delays",
               color="grey",ls="--",alpha=0.5)
    ax.axvline(np.log10(600),    label  ="10 min",     
               color="tab:green",ls="--",alpha=0.5)
    ax.axvline(np.log10(3600),   label ="1 h",        
               color="tab:orange",ls="--",alpha=0.5)
    ax.axvline(np.log10(3600*24),label="1 day",    
               color="tab:blue",ls="--",alpha=0.5)
    ax.set_title(title)
    ax.grid(which='both')
    ax.grid(which="major",ls=":")
    ax.set_xlabel(xlabel)
    
    return 

###-----------------------------------------------------------------------------
def col_size(var,var_min=1.1,var_max = 1000):
    
    import matplotlib.cm as cm
    # Limit values in case they are not yet limited
    var   = np.clip(var, var_min, None)
    
    color = cm.cool(np.log10(var)/np.log10(var_max))
    size  = 100*np.log10(var)**2

    return color, size

###-----------------------------------------------------------------------------
def sig_legend(ax,alpha=0.5, var_max=1000, **kwargs):
    
    siglist =  [5 , 10, 20, 50, 100, 500]
    labels = [str(x) for x in siglist]
    # symbol = Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="red")
    
    colors, sizes = col_size(siglist,var_max=var_max)
    patches = [ plt.plot([],[], marker="o", alpha=alpha,
                         ms = 7+sizes[i]/50, 
                         ls="", mec=None, 
                         color=colors[i], 
                         label="{:s}".format(labels[i]) )[0]  for i in range(len(labels)) ]
    # ax.legend(title="$\sigma_{max}$",handles=patches,**kwargs)
    
    return patches