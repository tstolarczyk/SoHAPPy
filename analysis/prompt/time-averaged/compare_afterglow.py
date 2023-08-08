# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 11:43:53 2022

@author: Stolar
"""
import os, sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import astropy.units as u

codefolder = "../../" 
sys.path.append(codefolder) 
os.environ['GAMMAPY_DATA'] =r'../../../input/gammapy-extra-master/datasets'

import warnings
#warnings.filterwarnings('error')
warnings.filterwarnings('ignore')

plt.style.use('seaborn-talk') # Make the labels readable
#plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!

from grb import GammaRayBurst
import grb_plot as gplt
from   utilities import  MyLabel

###----------------------------------------------------------------------------
def grblist_from_txt(filename, suffix="Event", debug=False):
    file = open(filename,"r")
    grblist = []
    for line in file.read().split('\n'):
        if line.find(suffix) !=-1:
            if debug: print(line, line[len(suffix):])
            grblist.append(int(line[len(suffix):]))
    print(len(grblist)," GRB found")
    return grblist
###----------------------------------------------------------------------------
def display(grb, save=False, res_dir= "./"):
    print(grb)
    gplt.time_spectra(grb)                
    gplt.energy_spectra(grb)   
    gplt.visibility_plot(grb, loc ="North")
    gplt.visibility_plot(grb, loc ="South")
    if (save): grb.write(res_dir) # Save GRB if requested    
    return
###----------------------------------------------------------------------------
def control_units(grb):
    print(32*"-")
    flx_unit = grb.spec_prompt[0].evaluate(100*u.GeV).unit
    print(" Prompt in : ",flx_unit,end="")
    if flx_unit != grb.spec_afterglow[0].evaluate(100*u.GeV).unit:
        print(" Afterglow : ",grb.spec_afterglow[0].evaluate(100*u.GeV).unit)
        sys.exit(" Change units")
    else:
        print(" same as Afterglow")
    print(" Energies  : ",grb.Eval.unit)
    print(32*"-")
    return flx_unit
###----------------------------------------------------------------------------
def integrate_dump_counts(grblist, 
                filename="dump_counts.txt", 
                grb_folder = "D:/CTAA/SoHAPPy/input/lightcurves/",
                visibility = None,
                prompt = None,
                e_unit = u.GeV,
                E_range = [30*u.GeV, 500*u.GeV],
                nmax=0, debug=False):
    
    import scipy.integrate as integrate
    from gammapy.irf import load_cta_irfs
    
    Emin = E_range[0].to(e_unit)
    Emax = E_range[1].to(e_unit)
    
    irf_dir   = "D:/CTAA/SoHAPPy/input/irf/Full/prod3-v2"
#     irf_dir   = "../../../input/irf/Full/prod3-v2"

    irf_fname = "FullArray/North/20deg/North_z20_average_30m/irf_file.fits.gz"
    irf_file  = Path(irf_dir,irf_fname)
    irf       = load_cta_irfs(irf_file)
    effarea   = irf["aeff"].data
    eff_unit  = effarea.evaluate(energy_true=100*u.GeV,offset=0*u.deg).unit
    
    out   = open(filename,"w")
    print("{:10} {:>10} {:>7} {:>10} {:>7} {:>10} {:>10}"
          .format("Name", "Prompt", "t90", "Afglw90", "dt90", "Afglw", 'dttot' ), file=out)
    N_prompt      = []
    N_afterglow   = []
    N_afterglow90 = []
    tobs          = []
    tobs90        = []
    
    first = True
    for j, item in enumerate(grblist[:nmax]):
        ### Read GRB
        filename = Path(grb_folder,"LONG_FITS","Event"+str(item)+".fits.gz")
        if debug: print(filename)
        else: print(j,end=" ")
        grb = GammaRayBurst.from_fits(filename, ebl = "dominguez",
                                  prompt = prompt, 
                                  vis = visibility,
                                  dbg = int(debug))
        # Show GRB characteristics
        if debug: display(grb)

        # Skip following if no prompt
        if grb.spec_prompt == None:
            print(filename, " Skipped (no prompt)")
            continue
        if first: 
            flx_unit = control_units(grb) 
            first=False

        ### Integrate, get prompt counts
        # E is a scalar, as well as the boundaries in integrate - ensure units are correct
        fval = lambda E : grb.spec_prompt[0].evaluate(E*e_unit).value \
                        * effarea.evaluate(energy_true=E*e_unit,offset=0*u.deg).value
        Nprompt = integrate.quad(fval,Emin.to_value(e_unit),Emax.to_value(e_unit))
        Nprompt = (Nprompt*grb.t90*flx_unit*eff_unit*e_unit).to("")
        N_prompt.append(Nprompt[0]) # Skip error

        ### Integrate afterglow
        Ntot  = 0
        N90   = 0
        dt90  = 0
        dttot = 0

        for i, flux in enumerate(grb.spec_afterglow[:-1]):

            if i!= 0:
                dt = (grb.tval[i]-grb.tval[i-1])
            else:
                dt = grb.tval[i]
                # print(" Obs point : ",grb.tval[i]," dt = ",dt)

            fval = lambda E : flux.evaluate(E*e_unit).value \
                            * effarea.evaluate(energy_true=E*e_unit,offset=0*u.deg).value
            F = integrate.quad(fval,Emin.to_value(e_unit),Emax.to_value(e_unit))
            F = (F*dt*flx_unit*eff_unit*e_unit).to("")    

            if i <= grb.id90: 
                N90 += F
                dt90 += dt

            Ntot += F
            dttot += dt

        if debug:

            print("COUNTS")
            print(" Prompt          : {:5.2e} +/- {:5.2e}, t90 = {:5.2f}, max tval={:5.2f}"
                  .format(Nprompt[0].value,Nprompt[1],grb.t90,grb.tval[grb.id90]))

            print(" Afterglow <t90  : {:5.2e} +/- {:5.2e}, tobs = {:5.2f}"
                  .format(N90[0].value,N90[1],dt90))
            print("           total : {:5.2e} +/- {:5.2e}, tobs = {:5.2f}"
                  .format(Ntot[0].value,Ntot[1],dttot))

        ### Write to file
        print("{:10s} {:10.2e} {:7.2f} {:10.2e} {:7.2f} {:10.2e} {:10.2f} "
          .format(grb.name, Nprompt[0], grb.t90.value, N90[0], dt90.value, Ntot[0], dttot.value ),file=out)

        ### Appends result to get the list
        N_afterglow.append(Ntot[0])
        N_afterglow90.append(N90[0])
        tobs.append(dttot)
        tobs90.append(dt90)
    print(" All done !")   
    out.close()
    
    # Go from list to numpy arrays (faster)
    N_prompt      = np.array(N_prompt)      
    N_afterglow   = np.array(N_afterglow)   
    N_afterglow90 = np.array(N_afterglow90) 
    tobs          = np.array([t.value for t in tobs])*tobs[0].unit
    tobs90        = np.array([t.value for t in tobs90])*tobs90[0].unit
    
    return N_prompt, N_afterglow, N_afterglow90, tobs, tobs90

###----------------------------------------------------------------------------
def read_counts(inputfile = "dump_counts.txt", t_unit = u.Unit("s"), debug=False):

    infile = open(inputfile,"r")
    lines = infile.readlines()

    N_prompt      = []
    N_afterglow   = []
    N_afterglow90 = []
    tobs          = []
    tobs90        = []

    for l in lines[1:]:
        [name, nP, t90, nA90, dt90, nA, dttot]  = [d for d in l.split()]
        N_prompt.append(float(nP)) # Skip error
        N_afterglow.append(float(nA))
        N_afterglow90.append(float(nA90))
        tobs.append(float(dttot))
        tobs90.append(float(dt90))
        if debug: print(name)

    N_prompt      = np.array(N_prompt)      
    N_afterglow   = np.array(N_afterglow)   
    N_afterglow90 = np.array(N_afterglow90) 
    tobs          = np.array(tobs)*t_unit
    tobs90        = np.array(tobs90)*t_unit

    print("All done !")
    return N_prompt, N_afterglow, N_afterglow90, tobs, tobs90 

###----------------------------------------------------------------------------
def plot(data, name="dummy", zoom=False,Nmin=0.1):
    from astropy.visualization import quantity_support
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
   
    N_prompt, N_afterglow, N_afterglow90, tobs, tobs90 = data
    
    fig,(ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3,figsize=(17,5))
    
    with quantity_support():
        cb = ax1.scatter(N_afterglow90, N_prompt, c=np.log10(tobs90.value) ,marker="o")
        plt.colorbar(cb, ax=ax1)
        if zoom:
            xmin = 1e-6
            ymin = Nmin
        else:
            xmin = np.min(N_afterglow90[N_afterglow90>0]) 
            ymin = np.min(N_prompt[N_prompt>0])
        xmax = np.max(N_afterglow90)
        ymax = np.max(N_prompt)
        print(xmin,xmax,ymin,ymax)
        ax1.plot([min(xmin,ymin), max(xmax,ymax)],
                [min(xmin,ymin), max(xmax,ymax)],ls="--",color="grey")
        if not zoom:
            ax1.add_patch(Rectangle((np.log10(xmin), np.log10(ymin)), 
                                    np.lo10(xmax)-np.log10(xmin), np.log10(ymax)-np.log1O(ymin),
                          edgecolor = 'red',
                          facecolor = None,
                          fill=False,
                          lw=1))

        ax1.set_xlim([xmin, xmax])
        ax1.set_ylim([ymin, ymax])
        ax1.set_xscale("log")    
        ax1.set_yscale("log")
        ax1.set_xlabel(" Afterglow counts, $t<t_{90}$ ") #"("+ax1.get_xlabel()+")")
        ax1.set_ylabel(" Prompt counts") # ("+ax1.get_ylabel()+")")
        ax1.grid(which="major",ls="-",alpha=0.5)
        ax1.grid(which="minor",ls=":",alpha=0.5)

        cb = ax2.scatter( N_afterglow90 , N_prompt/N_afterglow90,c=np.log10(tobs90.value) ,marker="o")
        plt.colorbar(cb, ax=ax2)
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        if zoom:
            ax2.set_xlim(1)
            ax2.set_ylim([1e-5,1])
        ax2.grid(which="both",alpha=0.5)
        ax2.set_xlabel(" Afterglow counts, $t<t_{90}$ ")
        ax2.set_ylabel(" Prompt / Afterglow count ratio")

        n, bins, _ = ax3.hist(np.log10(N_afterglow90[N_afterglow90>Nmin]),
                              alpha=0.5,bins=25,
                              label=MyLabel(N_afterglow90[N_afterglow90>Nmin]," Afterglow, $t<t_{90}$ ",stat="med"))
        ax3.hist(np.log10(N_prompt[N_prompt>Nmin]),alpha=0.5,bins=bins,
                 label=MyLabel(N_prompt[N_prompt>Nmin]," Prompt ",stat="med"))
        if zoom: ax3.set_xlim(Nmin)
        ax3.set_xlabel("log counts")
        ax3.axvline(x=1,label="N=10",ls=":",color="red")
        ax3.legend()
    fig.suptitle(name)
    plt.tight_layout()
    return


###############################################################################
if __name__ == "__main__":

    """
    This works with the 0-t90 prompt averaged flux.
    A standalone function to read prompt data and compare with afterglow.
    

    """