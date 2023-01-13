# -*- coding: utf-8 -*-
"""
Read a set of Fits files to access the source information in the form of a 
GammaRayBusrt SoHAPPy class. As this is time consuming, the classes can be 
stored on disk for further use.

Created on Mon Jan  9 14:33:09 2023

@author: Stolar
"""
import os


import pickle
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.visualization import quantity_support

from grb import GammaRayBurst
import configuration
from niceplot import MyLabel, single_legend

###----------------------------------------------------------------------------
def get_population(ifirst=1, ilast=10, folder="."):
    """
    Create a list of GammaRayBurst classes from the input fits file

    Parameters
    ----------
    ifirst : integer, optional
        first source number. The default is 1.
    ilast : integer, optional
        Last source folde. The default is 10.
    folder : string, optional
        Source input folder name. The default is ".".

    Returns
    -------
    grbpop : list
        List pf GammaRayBurst instances.

    """

    cf = configuration.Configuration()
    cf.ifirst     = ifirst
    cf.nsrc       = ilast
    grblist    = cf.source_ids()    
    
    grbpop=[]
    for i, item in enumerate(grblist):
        if ((cf.nsrc <= 10) or (np.mod(i,40) == 0)): print("#",i+1," ",end="")
        fname = Path(folder,"Event"+str(item)+".fits.gz" )
        grb = GammaRayBurst.from_fits(fname, ebl     = "dominguez")
        grbpop.append(grb)
    grbpop = np.array(grbpop)
    print("Population is created")
    return grbpop

###----------------------------------------------------------------------------
def check_binning(grbpop):
    nEval = []
    ntval = []
    for g in grbpop:
        nEval.append(len(g.Eval))
        ntval.append(len(g.tval))
    fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    ax1.hist(nEval,label=MyLabel(nEval))
    ax1.set_xlabel("Energy array length")
    ax1.legend()
    ax2.hist(ntval,label=MyLabel(ntval))
    ax2.set_xlabel("Time array length")
    ax2.legend()
    plt.show()
    
###----------------------------------------------------------------------------
def flux_peak_time(grb,Eref=100*u.GeV):
    
    idxref  = np.abs(grb.Eval-Eref).argmin() # Get flux close to Eref
    fluxref = grb.fluxval[:,idxref] # Obtain time spectrum for that energy
    idmax   = fluxref.argmax() # Get time index of max value  
    
    return idmax, idxref

###----------------------------------------------------------------------------
def flux_peak_time_all(energies, nbin=25):
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=len(energies), 
                                   figsize=(5*len(energies),10),
                                   sharey=True)

    for Eref, ax01, ax02 in zip(energies, ax1, ax2):
        tpeak = []
        Fpeak = []        
        for grb in grbpop:
            t_id, E_id = flux_peak_time(grb,Eref=Eref)
            tmax = grb.tval[t_id]
            Emax = grb.Eval[E_id]
            flux = grb.fluxval[t_id][E_id]
            # print(f" Max flux = {flux:8.2e} at t={tmax:8.2f} (E={Emax:5.2f})")
            tpeak.append(tmax.value)
            Fpeak.append(flux.value)
        Fpeak = np.array(Fpeak)*flux.unit
        tpeak = np.array(tpeak)*tmax.unit
        
        with quantity_support():
            ax01.hist(tpeak, bins=nbin,label=MyLabel(tpeak))
            ax01.set_title(f"{Emax.value:5.2f} {Emax.unit:s}")
            ax01.set_xlabel("$t_{peak}$ (" + ax01.get_xlabel()+")")
            ax01.legend()
            
            ax02.hist(Fpeak,bins=nbin)
            ax02.set_xlabel("$F_{peak}$ (" + ax02.get_xlabel()+")")
        
        plt.tight_layout()
            
###----------------------------------------------------------------------------
def check_spectra(ifirst=0, ilast=3, Eref=100*u.GeV):
    
    ngrb = ilast - ifirst
    t_unit = u.d
    fig, axall = plt.subplots(nrows=ngrb,ncols=1,figsize=(15,3.5*ngrb),
                              sharex=True)
    for grb, ax in zip(grbpop[ifirst:ilast+1],axall):
        
        times = grb.tval[1:].to(t_unit)
        dts   = (grb.tval[1:]-grb.tval[:-1]).to(t_unit)
        
        idxref  = np.abs(grb.Eval-Eref).argmin() # Get flux close to Eref
        fluxref = grb.fluxval[:,idxref] # Obtain time spectrum for that energy       
        idtfmax, _ = flux_peak_time(grb,Eref=Eref)
        fluxmax   = grb.fluxval[idtfmax][idxref]
        
        with quantity_support():
            ax.plot(times,dts,marker=".",alpha=0.5)
            ax.grid(ls="--",lw=1)
            ax.set_ylabel("$\Delta t$ (" + ax.get_ylabel()+")")
            ax.set_xscale("log")
            ax.set_xlabel(None)        
    
        axx = ax.twinx()
        
        label = f"{grb.id:s} \n {grb.Eval[idxref].value:5.2f} {Eref.unit:s}"
        with quantity_support():
            axx.plot(times,fluxref[1:],color="red",
                     ls="-",label=label,marker=".",alpha=0.5)
        
        axx.axhline(y=fluxmax/10,color="tab:green",label="10%",alpha=0.5)
        axx.axhline(y=fluxmax/100,color="tab:orange",label="1%",alpha=0.5)
        axx.axvline(x=1*u.d,label="Day 1",color="black",ls=":")
        axx.axvline(x=2*u.d,label="Day 2",color="darkgrey",ls=":")
        axx.axvline(x=3*u.d,label="Day 3",color="lightgrey",ls=":")
        axx.set_yscale("log")
        single_legend(axx, bbox_to_anchor=[1.15,1])    
        plt.tight_layout()

###############################################################################
if __name__ == "__main__":


    import sys
    CODE = "../../"
    sys.path.append(CODE)
    
    # Generic output file name and folders
    infolder    = "D:/CTA/CTAA/SoHAPPy/input/"
    outfolder   = "D:/CTA/CTAA/SoHAPPy/output/pop_tests/"
    grb_folder  = Path(infolder,"lightcurves/LONG_FITS/")
    os.environ['GAMMAPY_DATA'] =infolder+r'gammapy-extra-master/datasets'

    outfilename = "grbpop.bin"
    
    
    # This flags allows creating and storing the GRB classes and/or
    # reading them
    store = True # Store and process immediately
    store = False # Read a previously stored population 
    
    ### Read data from exisiting population binary or create it beforehand
    popname = Path(outfolder, "grbpop.bin")
    if store:## Store population
        grbpop = get_population(folder=grb_folder)
        outfile  = open(popname,"wb")
        pickle.dump(grbpop,outfile)
        outfile.close()
        print(" Population of ",len(grbpop)," GRBs saved to : {}".format(popname))
    else:
        outfile  = open(popname,"rb")
        grbpop= pickle.load(outfile)
        outfile.close()
        print(" Population of ",len(grbpop)," GRBs read from : {}".format(popname))   
        
    # Check binning
    check_binning(grbpop)
    
    # Time at which the flux is maximal
    energies = [50*u.GeV, 100*u.GeV, 200*u.GeV, 500*u.GeV]
    flux_peak_time_all(energies)
    
    # Check time intervals and spectral shape for some grbs
    check_spectra(ifirst=0, ilast=8, Eref=100*u.GeV)
    
