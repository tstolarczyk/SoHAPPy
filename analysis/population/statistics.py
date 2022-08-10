# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 13:50:29 2022

@author: Stolar
"""

import sys
import numpy as np

codefolder = "../../" 
sys.path.append(codefolder) 
from utilities import t_fmt

###-------------------------------------------------------------------
def detection_above_sigma(pop, s=5,cl=90):
    import astropy.units as u
    from astropy.time import Time
    from init import get_eff_lvl
    
    (grb, gn0, gs0, gn, gs, gb) = pop
    det_lvl = get_eff_lvl(grb)
    det_lvl = cl
            
    for g, tag in zip([gn0, gs0, gb], ["N only","S only", "Both"]):
        print("###------- ",tag,"-----")
        if s == 5:
            glist = g[g.d5s >= det_lvl ]
            print(" >= 5 sigma at {:2.0f}% CL :"
                  .format(100*det_lvl/max(set(grb.err))))
        elif s==3:
            glist = g[g.d3s >= det_lvl ]
            print(" >= 3 sigma at {:2.0f}% CL :"
                  .format(100*det_lvl/max(set(grb.err))))
        else:
            glist = g[ (g.d3s>=det_lvl) & (g.d5s<det_lvl)]
            print(" >= 3 and <5 sigma at {:2.0f}% CL :"
                  .format(100*det_lvl/max(set(grb.err))))        
            
        print("{:5s} {:3s}  {:5s} {:18s}  {:14s}    {:18s} {:23s} {:3s} {:3s}"
              .format(" ","vis","z","t5s","sigmax","tmax","date","%3s", "%5s"))

        for idx, g in glist.iterrows():
            t5s = t_fmt(g.t5s*u.s,digit=2)
            et5s = t_fmt( (g.et5s*u.s),digit=2).to(t5s.unit)
            
            tmx = t_fmt(g.tmx*u.s,digit=2)
            etmx = t_fmt( (g.etmx*u.s),digit=2).to(tmx.unit)
            
            ttrig = Time(g.ttrig,format="mjd").isot
            if tag=="Both":
                vis = str(int(gn[gn["name"]==g["name"]].vis))+"/"+str(int(gs[gs["name"]==g["name"]].vis))
            else:
                vis = g.vis
            print("{:5s} {:3s} {:5.1f} {:5.1f} +/- {:<5.1f} {:3s}  {:>5.1f} +/- {:>3.1f}   {:5.1f} +/- {:<5.1f} {:3s} {:23s} {:3.0f} {:3.0f}"
                  .format(g["name"][5:],str(vis),g.z,t5s.value,et5s.value,t5s.unit,g.sigmx,g.esigmx,tmx.value, etmx.value,tmx.unit,ttrig,g.d3s,g.d5s))
    return

###-------------------------------------------------------------------
def rate(grb, nyears = 0, det_lvl=0.9, summary=False, header=True):
    """


    Parameters
    ----------
    grb : Pandas data
        DESCRIPTION.
    nyears : TYPE, optional
        DESCRIPTION. The default is 1.
    det_lvl : TYPE, optional
        DESCRIPTION. The default is 0.9.

    Returns
    -------
    None.

    """
    
    #--------------------------------------------------------
    # Internal functions
    #--------------------------------------------------------
    def separator():
        """
        Just a separator line
        """
        print(" ----------- {:>15} {:>15}".format(14*"-",14*"-"),end="")
        print("{:>16} {:>15} {:>15} {:>15}".format(14*"-",14*"-",14*"-",14*"-"))       
        return
     
    #--------------------------------------------------------
    def stat_mean(glist,tag="",ny=1):
        """
        Print mean detection above a given significance

        Parameters
        ----------
        glist : population list
            list of populations to be analysed.
        tag : String, optional
            "3s" or "5s" for 2 and 5 sigma respectively. The default is "".
        ny : integer, optional
            Number of years, normalisation. The default is 1.

        Returns
        -------
        None.

        """
       
        [gn,gs,gb,gn0,gs0] = glist
        glist2 = [gn,gs,gb,gn0,gs0, [gn0, gs0, gb] ]

        #---------------------------------------------------    
        # Mean detection above a certain significance
        def prt_mean(sig, ny=1,tag=None):
            if tag!=None: 
                print(" {:9s} :".format(tag),end="") ### !!!
            else:    
                tag="yr-1"
                print(" {:>9s} :".format(tag),end="") ### !!!
            for i, g in enumerate(glist2):
                if sig == "3s": 
                    if type(g) == list:
                        var = sum([sum(x.d3s) for x in g])
                    else:
                        var = sum(g.d3s)
                elif sig == "5s": 
                    if type(g) == list:
                        var = sum([sum(x.d5s) for x in g])
                    else:
                        var= sum(g.d5s) 
                else:
                    sys.exit(" Should be '3s' or '5s'")
                
                nmean = var/niter
                # print(" ***",i,nmean)
                print(" {:7.1f} +- {:4.1f}"
                      .format(nmean/ny, np.sqrt(nmean)/ny),end="" )
            return
        
        if ny !=1: 
            separator()
            prt_mean(tag, ny=1,tag= tag)
            print()
            prt_mean(tag, ny=ny)
            print()
        else:
            separator()
            prt_mean(tag, ny=ny, tag=tag)
            print()
           
        return
    #--------------------------------------------------------
    def stat_line(glist,tag="",ny=1):
        """
        Print statistics in various conditions
        """
        
        [gn,gs,gb,gn0,gs0] = glist
        glist2 = [gn,gs,gb,gn0,gs0, gn0+gs0+gb]
        #--------------------------------------------------------
        # Stat line internal functions
        #--------------------------------------------------------        

    
        def prt_line(ny=1,tag="dummy"):
            if ny==1:
                print(" {:9s} :".format(tag),end="") ### !!!
            else:
                print(" {:>9s} :".format("yr-1"),end="") ### !!! 
            # for g in [gn,gs,gn0,gs0,gb,gn0+gs0+gb]:
            for g in glist2:
                print(" {:7.1f} +- {:4.1f}"
                      .format(len(g)/ny, np.sqrt(len(g))/ny),end="" )  
            return
        
        def prt_vis():
            print(" {:>9s} :".format("@trig"),end="")
            for g in glist2:
                r = len(g[g.vis==1])/len(g) if len(g) else 0
                print(" {:7d}  {:5.1f}%"
                      .format(len(g[g.vis==1]),100*r),end="" ) 
            return
        #--------------------------------------------------------        
   
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
            
        return # End of statline function
    #--------------------------------------------------------
    
    ### Back to the rate function
    
    # Get some parameters from the data
    niter  = max(set(grb.err))  # Number of MC iteration
    cl_min = det_lvl*niter      # Min. number of successful iteration
    unvis  = min(set(grb.err))  # Check if visible

    # Check that N only an S only tags exist
    suppinfo = ("N" in grb) and ("S" in grb) and ("B" in grb)
    if (not suppinfo):
        print(" North/south combinatory does not exist")
        return

    # Header
    if header:
        print()
        print("",107*"-")
        if nyears:
            print(" Normalized to {} year".format(nyears),end="")
        if (nyears>1): print("s")
        else: print()
        
        print("",107*"-")
        print(" Rate : {:>15} {:>15}".format("N","S"),end="")
        print("{:>16} {:>15} {:>15} {:>15}".format("Nonly","Sonly","Both","Total"))
            
    # Population base - visible
    g_ana = grb[grb.err != unvis]
    gn =  g_ana[g_ana.site=="North"]
    gs =  g_ana[g_ana.site=="South"]
    gb =  g_ana[g_ana.site=="Both"]
    gn0 = g_ana[g_ana.N == 1]
    gs0 = g_ana[g_ana.S == 1]
    
    poplist = [gn,gs,gn0,gs0, gb]
    if not summary: stat_line(poplist,tag="Vis.",ny=nyears)
        
    # # # Analysed
    # poplist = [g[g.err == niter] for g in poplist]
    # if not summary: stat_line(poplist,tag="Ana.",ny=nyears)

    # 3 sigma mean detection
    stat_mean(poplist,tag="3s",ny=nyears)
    
    # 5 sigma mean detection
    stat_mean(poplist,tag="5s",ny=nyears)

    # 3 sigma 90%CL detected
    poplist = [g[g.d3s >= cl_min] for g in poplist]
    if not summary: stat_line(poplist,tag="3s 90%CL",ny=nyears)

    # 5 sigma 90%CL detected
    poplist = [g[g.d5s >= cl_min] for g in poplist]
    stat_line(poplist,tag="5s 90%CL",ny=nyears)

    print("",107*"-")

    return

###-------------------------------------------------------------
if __name__ == "__main__":

    import sys
    codefolder = "../../" 
    sys.path.append(codefolder) 
    from init import create_csv, get_data
    
    nyears, file = create_csv(file="parameter.yaml",debug=False)
    (grb, gn0, gs0, gn, gs, gb) = get_data(file, debug=False)

    rate(grb,nyears=nyears)
    # detection_above_sigma((grb, gn0, gs0, gn, gs, gb),s=[3,5])