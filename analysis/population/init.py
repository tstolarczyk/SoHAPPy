# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:47:25 2020

@author: Stolar
"""
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from   astropy.table import Table
import pandas as pd

import setup
from utilities import MyLabel
from setup import col_3, col_5

__all__ = ["get_data", "create_csv"]

###-------------------------------------------------------------
def get_eff_lvl(grb, debug=False):
    """
    Compute absolute confidence level from the CL in percentage 
    and the number of events

    Parameters
    ----------
    grb : Pandas frame
        A GRB analysed population.

    Returns
    -------
    eff_lvl : Integer
        The absolute confidence level in counts.

    """
    
    import mcsim_config as mcf
    detlvl  = mcf.det_level # get this from the config file in case of changes
    niter   = max(set(grb.err)) 
    eff_lvl = detlvl * niter
    if debug:
        print(" Det. level = {}, niter = {} -> Eff. det. level is {}"
        .format(detlvl, niter, eff_lvl))

    return eff_lvl

###-------------------------------------------------------------
def get_data(filename,maxgrb=100000, debug=False):
    """
    
    Get data from a csv SoHAPPy output file
    Compute combinatory of site visibility.
    
    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    maxgrb : integer, optional
        Restrict to a maximal number of entries. The default is a lot.
    debug : Boolean, optional
        If True, lets' talk a bit. The default is False.

    Returns
    -------
    grb : Pandas dataframe
        The full GRB population (Entries North, South and Both).
    gn0 : Pandas dataframe
        DESCRIPTION.
    gs0 : Pandas dataframe
        DESCRIPTION.
    gn : Pandas dataframe
        DESCRIPTION.
    gs : Pandas dataframe
        DESCRIPTION.
    gb : Pandas dataframe
        DESCRIPTION.

    """

    grb    = pd.read_csv(filename)
        
    # Extract "not visible" flag and iteration from the data
    unvis     = min(set(grb.err))
    niter     = max(set(grb.err))
    niter_3s  = max(grb.d3s)
    niter_5s  = max(grb.d5s)

    # Site population
    gs = grb[grb.site=="South"][:maxgrb]
    gn = grb[grb.site=="North"][:maxgrb]
    gb = grb[grb.site=="Both" ][:maxgrb]

    # Get GRB names in the data
    if ( len(gs.name) != len(gn.name) or len(gs.name) != len(gb.name) ):
        print(" Inconstistency in GRB name list - CHECK ")
    else:
        names = gs[:maxgrb].name
    
    # Warn user in case the maximal GRB limit is below the data content
    if len(names) > maxgrb:
        print(" WARNING : limited statistics (",len(names),") >",maxgrb)

    # Add combinatory to data frame for the name list
    add_combinatory(names,grb,unvis=unvis,debug=False)
    suppinfo = ("N" in grb) and ("S" in grb) and ("B" in grb)

    # The iteration number of the simulation can be guessed
    # from the error code (except if all simulations failed!)
    g_ana = grb[ grb.err   == niter] # All iterations were simulated
    gn0   = g_ana[ (g_ana.site =="North") & (g_ana.N==1)] # North only
    gs0   = g_ana[ (g_ana.site =="South") & (g_ana.S==1)] # South only
    gn    = g_ana[  g_ana.site =="North"] # North and maybe elsewhere
    gs    = g_ana[  g_ana.site =="South"] # South and maybe elsewhere
    gb    = g_ana[ (g_ana.site =="Both")  & (g_ana.B==1)] # Seen both

    if (debug):
        print(" DATA READING from ",filename)
        if (suppinfo):
            print("Supplementary information is present")
            print(" grb.N==1 seen North only")
            print(" gbr.S==1 seen South only")
            print(" grb.B==1 seen on both")
            print()
        print("+-------------------------- Flags ---------------------------+")
        print(" Flags:")
        print("   No visible flag, unvis           = ",unvis)
        print("   Iteration # from error code, 3s and 5s counts : ",
              niter, niter_3s, niter_5s)

        print("+----------------------- Statistics -------------------------+")
        print(" {:^15s} {:^15s} {:^15s}"
          .format("Not visible","Fully analyzed","Aborted"))
        print(" {:^15d} {:^15d} {:^15d}"
          .format(len(grb[  grb.err == unvis]),
                  len(grb[  grb.err == niter]),
                  len(grb[ (grb.err != niter) & (grb.err!=unvis) ])))
        print()
        print(" Raw statistics - max per site =",maxgrb)
        print("  - total      : ",len(grb))
        print("  - analyzable : ",len(g_ana))
        print("  - North      : ",len(gn))
        print("  - South      : ",len(gs))
        print("  - Both sites : ",len(gb),
              "-> total = ",len(gn)+len(gs)+len(gb))

        print("  - North only : ",len(gn0))
        print("  - South only : ",len(gs0))

        print("+------------------------------------------------------------+")

    return (grb, gn0, gs0, gn, gs, gb)

###-------------------------------------------------------------------
def sanity_check(file, grb, gn0, gs0, gn, gs, gb,
                 maxgrb=2000, debug=False):
    print("+======================== Sanity checks =========================+")


    # Min altitude
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10,2))
    for axi, gpop, tag in zip(ax,[gn,gs,gb],["North","South","Both"]):
        axi.hist(gpop.altmx,bins=100,label=tag)
        axi.set_title("Altitude at max $\sigma$")
        axi.legend()
        if tag != "Both":
            print(" Estimated min altitude in ",tag," :",min(gpop.altmx))
    plt.show()

    # Slewing delay - check with values in the configuration file
    cf = get_config(file)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (10,2))

    for axi, gpop, loc in zip(ax,[gn,gs],["North","South"]):
        axi.hist(gpop[gpop.t3s>=0].t3s,bins=100,label=loc)
        delay = cf.dtslew[loc]
        axi.axvline(x = delay.value,
                    color="red",ls=":",
                    label=str(delay.value)+" "+str(delay.unit))
        delay = delay+cf.dtswift
        axi.axvline(x= delay.value,
                   color="green",ls=":",
                   label=str(delay.value)+" "+str(delay.unit))

        axi.set_title("Time delay to $3 \sigma$  - "+loc)
        axi.set_xlim(xmin=0,xmax=10+1.3*delay.value)
        #axi.set_yscale("log")
        axi.legend()
        print(" Estimated total delay in {:5s}".format(loc))
        print("    From visibility start : {:5.1f}"
              .format(min(gpop[gpop.t1>=0].t1)))
        print("    From 3s detection     : {:5.1f}"
              .format(min(gpop[gpop.t3s>=0].t3s)))

    print("+================================================================+")

    return

###-------------------------------------------------------------------
def rate(grb, nyears = 0, det_lvl=0.9, summary=False, header=True):
    """


    Parameters
    ----------
    grb : TYPE
        DESCRIPTION.
    nyears : TYPE, optional
        DESCRIPTION. The default is 1.
    det_lvl : TYPE, optional
        DESCRIPTION. The default is 0.9.

    Returns
    -------
    None.

    """
    niter  = max(set(grb.err))
    cl_min = det_lvl*niter
    unvis  = min(set(grb.err))

    # Check that N only an S only tags exist
    suppinfo = ("N" in grb) and ("S" in grb) and ("B" in grb)
    if (not suppinfo):
        print(" North/south combinatory does not exist")
        return

    # Header
    if header:
        print()
        print("",102*"-")
        if nyears:
            print(" Normalized to {} year".format(nyears),end="")
        if (nyears>1): print("s")
        else: print()
        
        print("",102*"-")
        print(" Rate : {:>15} {:>15}".format("N","S"),end="")
        print("{:>16} {:>15} {:>15} {:>15}".format("Nonly","Sonly","Both","Total"))

    #--------------------------------------------------------
    def separator():
         print(" ------ {:>15} {:>15}".format(14*"-",14*"-"),end="")
         print("{:>16} {:>15} {:>15} {:>15}".format(14*"-",14*"-",14*"-",14*"-"))       
         return
     
    # def stat_line(gn,gs,gb,gn0,gs0,tag="",ny=1):
    def stat_line(glist,tag="",ny=1):
        [gn,gs,gb,gn0,gs0] = glist
        glist2 = [gn,gs,gb,gn0,gs0, gn0+gs0+gb]
        
        def prt_line(ny=1,tag="dummy"):
            print(" {:5s}:".format(tag),end="")
            # for g in [gn,gs,gn0,gs0,gb,gn0+gs0+gb]:
            for g in glist2:
                print(" {:7.1f} +- {:4.1f}"
                      .format(len(g)/ny, np.sqrt(len(g))/ny),end="" )  
            return
        
        def prt_vis():
            print(" {:5s}:".format("@trig"),end="")
            for g in glist2:
                r = len(g[g.vis==1])/len(g) if len(g) else 0
                print(" {:7d}  {:5.1f}%"
                      .format(len(g[g.vis==1]),100*r),end="" ) 
            return

        if ny !=1: 
            separator()
            prt_line(ny=1,tag = tag)
            print()
            prt_line(ny=ny, tag=" ")
            print()
            prt_vis()
        else:
            separator()
            prt_line(ny=ny, tag=tag)
        
        print()
        return
    #--------------------------------------------------------
    # Population base - visible
    g_ana = grb[grb.err != unvis]
    gn =  g_ana[g_ana.site=="North"]
    gs =  g_ana[g_ana.site=="South"]
    gb =  g_ana[g_ana.site=="Both"]
    gn0 = g_ana[g_ana.N == 1]
    gs0 = g_ana[g_ana.S == 1]
    
    
    poplist = [gn,gs,gn0,gs0, gb]
    if not summary: stat_line(poplist,tag="Vis.",ny=nyears)
        
    # Analysed
    poplist = [g[g.err == niter] for g in poplist]
    if not summary: stat_line(poplist,tag="Ana.",ny=nyears)

    # 3 sigma detected
    poplist = [g[g.d3s >= cl_min] for g in poplist]
    if not summary: stat_line(poplist,tag="3s",ny=nyears)

    # 5 sigma detected
    poplist = [g[g.d5s >= cl_min] for g in poplist]
    stat_line(poplist,tag="5s",ny=nyears)


    print("",102*"-")

    return
###-------------------------------------------------------------------
def add_combinatory(names,grb,unvis=-999,debug=False):
    """
    Add combinatory to data frame
    """

    # If columns do not exist, create them
    if "N" not in grb: grb.insert(1,"N",0) # North only
    if "S" not in grb: grb.insert(1,"S",0) # South only
    if "B" not in grb: grb.insert(1,"B",0) # North and South

    if (debug):
        print("{:>10s} {:>3s} {:>3s} {:>3s} {:>3s} {:>3s}"
              .format("name","N","S","B","No","So"))

    for name in names:
        g = grb[grb.name==name]
        seen_n = (g[g.site=="North"].err!=unvis).bool()
        seen_s = (g[g.site=="South"].err!=unvis).bool()
        seen_b = (g[g.site=="Both"].err!=unvis).bool()
        seen_sonly = seen_s & ~seen_n
        seen_nonly = seen_n & ~seen_s
        if (debug):
            print("{:>10s} {:3d} {:3d} {:3d}{:3d} {:3d}"
                  .format(name,
                        int(seen_n),
                        int(seen_s),
                        int(seen_b),
                        int(seen_nonly),
                        int(seen_sonly)))
        #print(g.index)
        for idx in g.index:
            # Not in pandas 1.0.3
            # grb.set_value(idx,"N",int(seen_nonly))
            # grb.set_value(idx,"S",int(seen_sonly))
            # grb.set_value(idx,"B",int(seen_b))
            grb.at[idx,"N"] = int(seen_nonly)
            grb.at[idx,"S"] = int(seen_sonly)
            grb.at[idx,"B"] = int(seen_b)
    return

###-------------------------------------------------------------
def get_config(file, debug=False):
    """
    Get configuration file from the csv file name

    Parameters
    ----------
    file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    dirname = Path(file).parents[0].absolute() 
    conf_file = None
    
    for f  in dirname.iterdir():
        if f.suffix == ".yaml": 
            conf_file = f
            if debug: print(" Found configuration file :",conf_file)
            
    if conf_file == None:
        sys.exit(" No configuration file found in",dirname)
    else:
        from configuration import Configuration
        cf = Configuration([],conf_file=conf_file, vis_file=None)
        
    if debug:
        from utilities import Log
        cf.print(log=Log(name  = "tobedeleted.log", 
                         talk=not cf.silent) )
    
    return cf

###-------------------------------------------------------------
def create_csv(file="parameter.yaml", datafile="data.txt", debug=False):
    """
    From the current parameter file containg the folder to be analysed,
    create the csv file from default txt file.
    If the paramter file is not given, then the datafile is supposed to 
    contain the full path data file name.
    Get the simulation duration. 
    
    Parameters
    ----------
    file : String, optional
        Input parameter file name, can be None. The default is `paramter.yaml`.
    datafile : TYPE, optional
        DESCRIPTION. The default is "data.txt".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    nyears : FLOAT
    csvfilename : STRING

    None.

    """
    # Get input folder from the paramter file if defined - create filename
    if file != None:
        import yaml
        from yaml.loader import SafeLoader
        folder = (yaml.load(open(file), Loader=SafeLoader))["outfolder"]
        base   = (yaml.load(open(file), Loader=SafeLoader))["base_folder"]
        nyears = (yaml.load(open(file), Loader=SafeLoader))["duration"]
        txtfilename = Path(base+folder, datafile) 
    else:
        txtfilename = Path(datafile)
        
    if (debug): print("Full name :",txtfilename)
    
    # Build csv filename
    csvfilename = txtfilename.with_suffix('.csv')
    if (debug): print(" >>> ",csvfilename)

    # Check existence of csv file, try to create it from the txt file
    if not csvfilename.is_file():
        if txtfilename.is_file():
            print("Text file found, converting...")
            data = Table.read(txtfilename.as_posix(),format="ascii",guess=False)
            data.write(csvfilename.as_posix(), format="ascii.csv",overwrite=True)
            #data.write(filename+'.fits',format="fits",     overwrite=True)
            print(csvfilename," Created")
        else:
            sys.exit(" *** Requested data file is not available as .csv or .txt ***")
    else:
        if debug: print(csvfilename," exists")

    return nyears, csvfilename

###-------------------------------------------------------------------------
def computing_time(gpop,eff_lvl, nbin=25, ax=None, **kwargs):
    
    niter = max(gpop.err)

    ax = plt.gca() if ax is None else ax
    
    n, bins, _ = ax.hist(niter*gpop.mct,bins=nbin,
                          facecolor="none",
                          edgecolor="black",
                          label=MyLabel(niter*gpop.mct,stat="med"))

    ax.hist(niter*gpop[gpop.d3s>eff_lvl].mct,
             bins=bins,
             alpha=0.5,
             facecolor = col_3, 
             label=MyLabel(niter*gpop[gpop.d3s>eff_lvl].mct,
                           "$3\sigma$ det.",stat="med"),
             **kwargs)

    ax.hist(niter*gpop[gpop.d5s>eff_lvl].mct,
             bins=bins,
             alpha=0.5,
             facecolor = col_5, 
             label=MyLabel(niter*gpop[gpop.d5s>eff_lvl].mct,
                           "$5\sigma$ det.",stat="med"),
             **kwargs)
    
    ax.set_title(f"Simulation duration ({niter:3d} iter.)")
    ax.set_xlabel("Time (s)")
    ax.set_yscale("log")
    ax.legend(bbox_to_anchor=[1,1])

    return
###-------------------------------------------------------------
def detection_level(var,cl=0.9,nbin=25, ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax
    
    ax.hist(100*var, bins=nbin, range=[0,100], **kwargs)
    ax.set_yscale("log")
    ax.axvline(x=100*cl,
               color="red",ls=":",label="min. C.L. (%)")
    ax.legend()
    ax.set_xlabel("Confidence level (%)")
    ax.set_ylabel("Event counts")
    ax.grid("both")
    
    return

###-------------------------------------------------------------
if __name__ == "__main__":

    # file = create_csv(file="analysis/population/parameter.yaml",debug=True)
    file = create_csv(file="parameter.yaml",debug=True)
    (grb, gn0, gs0, gn, gs, gb) = get_data(file, debug=True)
    # sanity_check(file, grb, gn0, gs0, gn, gs, gb, debug=True) 
    # rate(grb)
    rate(grb,nyears=44)
    # computing_time(grb)
