# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:47:25 2020

@author: Stolar
"""
import os
import sys
import matplotlib.pyplot as plt
from   astropy.table import Table
import pandas as pd

__all__ = ["get_data", "create_csv"]
###-------------------------------------------------------------
def get_data(filename,maxgrb=2000, debug=False):
    """
    Get data from a csv SoHAPPy output file
    Compute combinatory of site visibility

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    grb : TYPE
        DESCRIPTION.
    g_ana : TYPE
        DESCRIPTION.
    gn : TYPE
        DESCRIPTION.
    gs : TYPE
        DESCRIPTION.
    gb : TYPE
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

    # Add cobinatory to data frame for the name list
    add_combinatory(names,grb,unvis=unvis,debug=False)
    suppinfo = ("N" in grb) and ("S" in grb) and ("B" in grb)

    # The iteration number of the simulation can be guessed
    # from the error code (except if all simulations failed!)
    g_ana = grb[ grb.err   == niter]
    gn0 = g_ana[ (g_ana.site=="North") & (g_ana.N==1)] # North only
    gs0 = g_ana[ (g_ana.site=="South") & (g_ana.S==1)] # South only
    gn  = g_ana[ g_ana.site=="North"] # North and maybe elsewhere
    gs  = g_ana[ g_ana.site=="South"] # South and maybe elsewhere
    gb  = g_ana[ (g_ana.site=="Both")  & (g_ana.B==1)] # Seen both

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
                  len(grb[  grb.err == niter_3s]),
                  len(grb[ (grb.err != niter_3s) & (grb.err!=unvis) ])))
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
def sanity_check(grb, gn0, gs0, gn, gs, gb,
                 filename,maxgrb=2000, debug=False):
    print("+----------------------- Sanity checks----------------------+")

    import sys
    # import config as cf # Load local config file


    # Min altitude
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10,2))
    for axi, gpop, tag in zip(ax,[gn,gs,gb],["North","South","Both"]):
        axi.hist(gpop.altmx,bins=100,label=tag)
        axi.set_title("Altitude at max $\sigma$")
        axi.legend()
        if tag != "Both":
            print(" Estimated min altitude in ",tag," :",min(gpop.altmx))
    plt.show()


    # Slewing delay
    outfolder = os.path.split(filename)[0]
    sys.path.append(outfolder)
    import config as cf

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (10,2))

    for axi, gpop, loc in zip(ax,[gn,gs],["North","South"]):
        axi.hist(gpop[gpop.t3s>=0].t3s,bins=100,label=loc)
        delay = cf.dtslew[loc]
        axi.axvline(x = delay.value,
                    color="red",ls=":",
                    label=str(delay.value)+" "+str(delay.unit))
        delay = cf.dtslew[loc]+cf.dtswift
        axi.axvline(x= delay.value,
                   color="green",ls=":",
                   label=str(delay.value)+" "+str(delay.unit))

        axi.set_title("Time delay to $3 \sigma$  - "+loc)
        axi.set_xlim(xmin=0,xmax=200)
        #axi.set_yscale("log")
        axi.legend()
        print(" Estimated total delay in ",loc," :",min(gpop[gpop.t3s>=0].t3s))
    return

###-------------------------------------------------------------------
def rate(grb, nyears = 1, det_lvl=0.9):
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
    print()
    print("",55*"-")
    print(" Normalized to {} year".format(nyears),end="")
    if (nyears>1): print("s")
    else: print()
    print("",55*"-")
    print(" Rate : {:>7} {:>7}".format("N","S"),end="")
    print("{:>8} {:>7} {:>7} {:>7}".format("Nonly","Sonly","Both","Total"))

    print(" ------ {:>7} {:>7}".format(7*"-",7*"-"),end="")
    print("{:>8} {:>7} {:>7} {:>7}".format(7*"-",7*"-",7*"-",7*"-"))

    #--------------------------------------------------------
    def stat_line(gn,gs,gb,gn0,gs0,tag="",ny=1):

        print(" {:5s}: {:7.1f} {:7.1f}"
              .format(tag,
                      len(gn)/ny,
                      len(gs)/ny),end="")

        print("{:>8.1f} {:>7.1f} {:>7.1f} {:>7.1f}"
                  .format(len(gn0)/ny,
                          len(gs0)/ny,
                          len(gb)/ny,
                          (len(gn0)+len(gs0)+len(gb))/ny ))
        return
    #--------------------------------------------------------
    # Population base - visible
    g_ana = grb[grb.err != unvis]
    gn =  g_ana[g_ana.site=="North"]
    gs =  g_ana[g_ana.site=="South"]
    gb =  g_ana[g_ana.site=="Both"]
    gn0 = g_ana[g_ana.N == 1]
    gs0 = g_ana[g_ana.S == 1]
    stat_line(gn,gs,gb,gn0,gs0,tag="Vis.",ny=nyears)

    # Analysed
    stat_line(gn[gn.err == niter],
              gs[gs.err == niter],
              gb[gb.err == niter],
              gn0[gn0.err == niter],
              gs0[gs0.err == niter],
              tag="Ana.",ny=nyears)

    # 3 sigma detected
    stat_line(gn[gn.d3s>cl_min],
              gs[gs.d3s>cl_min],
              gb[gb.d3s>cl_min],
              gn0[gn0.d3s>cl_min],
              gs0[gs0.d3s>cl_min],
              tag="3s",ny=nyears)

    # 5 sigma detected
    stat_line(gn[gn.d5s>cl_min],
              gs[gs.d5s>cl_min],
              gb[gb.d5s>cl_min],
              gn0[gn0.d5s>cl_min],
              gs0[gs0.d5s>cl_min],
              tag="5s",ny=nyears)
    print("",55*"-")

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
def create_csv(folder=None, datafile="data.txt", debug=False):
    """
    From the SoHAPPy outputfoler, create the csv file from default txt file
    if necessary, returns the filen ame

    Parameters
    ----------
    folder : TYPE, optional
        DESCRIPTION. The default is None.
    datafile : TYPE, optional
        DESCRIPTION. The default is "data.txt".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    csvfilename : STRING

    None.

    """
    prefix = os.path.splitext(datafile)[0] # Get prefix
    # print(os.listdir(outfolder))

    csvfilename = folder+prefix+".csv" # Build csv filename
    if (debug): print(" >>> ",csvfilename)

    # Check existence of file, try to create it from the txt file
    csv_exists = os.path.isfile(csvfilename)
    if not csv_exists:
        fullname =folder + datafile
        if (debug): print("Full name :",fullname)
        txt_exists =  os.path.isfile(fullname)
        if (txt_exists):
            print("Text file found, converting...")
            data = Table.read(fullname,format="ascii",guess=False)
            data.write(csvfilename, format="ascii.csv",overwrite="True")
            #data.write(filename+'.fits',format="fits",     overwrite="True")
            print(csvfilename," Created")
        else:
            sys.exit(" *** Requested data file is not available as .csv or .txt ***")
    else:
        print(csvfilename," exists")

    return csvfilename

###-------------------------------------------------------------
if __name__ == "__main__":
    print("mains")


if __name__ == "__main__":
    print("mains")

    import init as init

    outfolder= "../../../output/new17-107-24deg-newirf/"
    file = init.create_csv(outfolder)
    (grb, gn0, gs0, gn, gs, gb) = get_data(file)

    rate(grb)
    rate(grb,nyears=44)
