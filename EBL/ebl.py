# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 15:24:21 2021

@author: Stolar
"""

import sys
codefolder = "../" 
sys.path.append(codefolder) 
import numpy as np
import astropy.units as u
from scipy.interpolate import interp2d
from utilities import single_legend
#------------------------------------------------------------------------------------------
def EBL_from_file(file, debug=False):
    data = np.loadtxt(file)

    # Load data from file
    if file.find("gilmore") != -1:
        z_dat = data[0,1:]
        E_dat = (data[1:,0]/1000)*u.GeV
        EBLabs = np.array(data[1:,1:]).T
        EBLabs = np.exp(-EBLabs)
        if debug: print(" Gilmore from Lara")
        
    elif file.find("dominguez") != -1:
        E_dat = data[0,1:]*u.GeV
        z_dat = data[1:,0] 
        EBLabs = data[1:,1:]  
        if debug: print(" Dominguez from Renaud")

    attenuation =interp2d(E_dat,z_dat, EBLabs)        
    return attenuation

#------------------------------------------------------------------------------------------
from itertools import cycle
import matplotlib.cm as cm
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt
import numpy

def EBL_plot(att_list, zlist=0.1, elist=100*u.GeV, tags="",
             ax = None, xsize=8, ysize=8, 
             color=None, ls = None,
             debug=False, **kwargs):
    
    # If not a list but a single element
    if not isinstance(att_list, list): att_list = [att_list]
    if not isinstance(tags, list): tags = [tags]
    if not isinstance(zlist,numpy.ndarray): zlist = [zlist]
    if not isinstance(elist,numpy.ndarray): elist = [elist]
        
    lines = ["-","--","-.",":"] # A line style per attenuation model
    linecycler = cycle(lines)    
    
    set_color = False
    set_line = False
    if ax==None: fig,ax = plt.subplots(figsize=(xsize,ysize))
    if color == None: set_color = True
    if ls == None : set_line = True 
    
    with quantity_support():

        for att,tag in zip(att_list,tags):
            if set_line: ls = next(linecycler) 
            for i, z in enumerate(zlist):
                if set_color : color = cm.cool(i/len(zlist)) # One color per z
                ax.plot(elist,att(elist,z),
                        label=tag+" z="+str(round(z,2)),
                        ls=ls, color=color,**kwargs)
        ax.set_xscale("log")
        ax.set_ylim(ymin=1e-2,ymax=2)
        
        single_legend(ax)
    return ax

#------------------------------------------------------------------------------
if __name__ == "__main__":
        
    ### Read absorptions
    eblfilename = "data/"+"ebl_gilmore12.dat"
    gilmore = EBL_from_file(eblfilename)

    eblfilename = "data/EBL_abs_RBelmont-dominguez-20170425.dat"
    dominguez = EBL_from_file(eblfilename)
    
    ### Plot absorptions
    zmin  = 1
    zmax  = 5
    nzbin = 4
    emin = 10*u.GeV
    emax = 10*u.TeV 
    nebin=100
    zlist   = np.append([0], np.logspace(np.log10(zmin),np.log10(zmax),nzbin))    
    Elist   = np.logspace(np.log10(emin.value),np.log10(emax.to(emin.unit).value),nebin)*emin.unit
    taglist = ["gilmore", "dominguez"]
    ax = EBL_plot([gilmore, dominguez], tags=taglist, zlist=zlist, elist=Elist)
    ax.set_yscale("log")
  
    