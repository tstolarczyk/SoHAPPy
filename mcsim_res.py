# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""
import numpy as np
import astropy.units as u
import mcsim_config as mcf
import ana_config as cf
from utilities import warning, failure, success, highlight,banner

###########################################################################
def welcome():
    banner("+================================================================+")
    banner("||                     LAUNCHING SIMULATION                     ||")
    banner("+================================================================+")    
    print("   Detection level    : ",mcf.det_level)
    print("   Alpha              : ",mcf.alpha)
    print("   FOV                : ",mcf.fov)
    print("   Bin size           : ",mcf.binsize)
    print("   Slewing time       : ",cf.dtslew)
    print("      Fixed           : ",cf.fixslew)
#    print("  Reduction factor    : ",cf.redfactor)    
    
###########################################################################
def mean_values(mc, idlist):
    
    """
    Extract mean values and RMS from the index list
    """
    t_det, e_t_det, alt_det, e_alt_det, az_det, e_az_det  = (-1,-1,-1,-1,-1,-1)
    if (len(idlist)):
        # t = [mc.tobs[i][0].value for i in idlist
        #      if mc.tobs[i][0].value>=0]
        # Bug corrected April 27th, 2020 - obs time is end time in Roma version
        t = [mc.tobs[i][1].value for i in idlist 
             if mc.tobs[i][1].value>=0]
        t_det   = np.mean(t)
        e_t_det = np.std(t)
        
       # Bug corrected April 27th, 2020 - obs time is end time in Roma version
       # alt =  [mc.altaz[i][0].alt.value for i in idlist
       #          if mc.altaz[i][0].alt.value>=0]       
       alt =  [mc.altaz[i][1].alt.value for i in idlist
                if mc.altaz[i][1].alt.value>=0]
        alt_det   = np.mean(alt)
        e_alt_det = np.std(alt)
        
        # Bug corrected April 27th, 2020 - obs time is end time in Roma version
        # az =  [mc.altaz[i][0].az.value for i in idlist
        #        if mc.altaz[i][0].az.value>=0]
        az =  [mc.altaz[i][1].az.value for i in idlist
               if mc.altaz[i][1].az.value>=0]
        az_det   = np.mean(az)
        e_az_det = np.std(az)

    return  t_det, e_t_det, alt_det, e_alt_det, az_det, e_az_det  

###########################################################################
def result(mc,popfile=None,write_header="False"):

    """
    Print out results.
    If only one GRB is simulated, store the present mc class on disks
    for further investigation.
    """

    smax_mean, e_smax, nex_smax, e_nex_smax, nb_smax, e_nb_smax= 6*(-1,)
    tmx,e_tmx,altmx,e_altmx,azmx, e_azmx = 6*(-1,)
    
    nex3s, e_nex3s, nb3s, e_nb3s = 4*(-1,)    
    t3s,e_t3s,alt3s,e_alt3s,az3s, e_az3s = 6*(-1,)
    
    nex5s, e_nex5s, nb5s, e_nb5s = 4*(-1,)    
    t5s,e_t5s,alt5s,e_alt5s,az5s, e_az5s = 6*(-1,)
    
    # General status message
    if  (mc.err < 0): 
        print("+--------------+-------------------------------------------------+")
        print("| Status       |  ===> ",end="")
        failure("Simulation not possible (GRB not visible)",end="")
        print(" |")
        print("+--------------+-------------------------------------------------+")         
    if (mc.err> 0 and mc.err != mc.niter): 
        print("+--------------+-------------------------------------------------+")
        print("| Status       |  ===> ",end="")
        failure("Simulation was aborted at trial {:>3d}".format(mc.err),end="")
        print("       |")
        print("+--------------+-------------------------------------------------+")         
    if (mc.err == mc.niter):
        if (mc.detect_5s/mc.niter >= mcf.det_level):
            message = "5 sigma detected"
        elif (mc.detect_3s/mc.niter >= mcf.det_level):
            message = "3 sigma detected"
        else: 
            message = "NOT detected    "            
        print("\n+--------------+-------------------------------------------------+")
        print("| Status       |  ===> ",end="")
        failure(message,end="")
        print("                          |")
        print("+--------------+-------------------------------------------------+")         
    
    # Compute results if everything went well
    if (mc.err == mc.niter):
        #print("+--------------+-------------------------------------------------+")
        highlight("  RESULTS      : {:>10s} - {:<10s}"
                  .format(mc.grb.name,mc.where))
        #print("+--------------+-------------------------------------------------+")         
        print("+--------------+")         
        print(" Observation from {:6.2f} to {:6.2f} ({:5d} slices)"
              .format(mc.tstart,mc.tstop,len(mc.tobs)))
        print("  - Alt. (zen) :  {:6.2f} - {:6.2f} ({:6.2f} - {:6.2f})"
              .format(mc.pos_start.alt,
                      mc.pos_stop.alt,
                      90*u.deg-mc.pos_start.alt,
                      90*u.deg-mc.pos_stop.alt))
        print("  - Azimuth    :  {:6.2f} to {:6.2f}"
              .format(mc.pos_start.az,
                      mc.pos_stop.az))
        print(" E range (IRF) : {:6.2f} {:6.2f}"
              .format(mc.Emin[0], mc.Emax[0]))

        # Max significance is reached
        smax_mean  = np.mean(mc.smax_list)
        e_smax     = np.std(mc.smax_list)
        nex_smax   = np.mean(mc.nex_smax_list)
        e_nex_smax = np.std(mc.nex_smax_list)        
        nb_smax    = np.mean(mc.nb_smax_list)
        e_nb_smax  = np.std(mc.nb_smax_list)
    
        tmx,e_tmx,altmx,e_altmx,azmx, e_azmx = mean_values(mc,mc.id_smax_list)    
           
        ### 3 sigma
        if (len(mc.id_3s_list)):
            nex = [n for n in mc.nex_3s_list if n >=0]
            nb  = [n for n in mc.nb_3s_list  if n >=0]
            nex3s    = np.mean(nex)
            e_nex3s  = np.std(nex)
            nb3s     = np.mean(nb)
            e_nb3s   = np.std(nb)
        
        t3s,e_t3s,alt3s,e_alt3s,az3s, e_az3s = mean_values(mc,mc.id_3s_list)    

        ### 5 sigma

        if (len(mc.id_5s_list)):
            nex = [n for n in mc.nex_5s_list if n >=0]
            nb  = [n for n in mc.nb_5s_list  if n >=0]
            nex5s    = np.mean(nex)
            e_nex5s  = np.std(nex)
            nb5s     = np.mean(nb)
            e_nb5s   = np.std(nb)

        t5s,e_t5s,alt5s,e_alt5s,az5s, e_az5s = mean_values(mc,mc.id_5s_list)    
            
        ### Simplify err_slices if too long (remove duplicated slice id)
        if (mc.niter > 10):
            mc.err_slice = list(set(mc.err_slice))
        t_unit = u.s
        alt_unit = u.deg   
        
        print("+----------------------------------------------------------------+") 
        print("  Max. sigma           = {:>8.1f}  +/- {:<8.1f}"
              .format(smax_mean,e_smax))
    
        print("                  time = {:>8.1f}  +/- {:<8.1f} {:5s}"
               .format(tmx,e_tmx,t_unit))
        print("                  alt. = {:>8.1f}  +/- {:<8.1f} {:5s}"
               .format(altmx,e_altmx,alt_unit))
        print("                   az. = {:>8.1f}  +/- {:<8.1f} {:5s}"
               .format(azmx,e_azmx,alt_unit))
        print("           S, B counts = {:>8.1f} {:0.1f}"
              .format(nex_smax,nb_smax))
        print("        on, off counts = {:>8.1f} {:0.1f}"
              .format(nex_smax+nb_smax,nb_smax/mcf.alpha))
        print()
    
        print("  3 sigma :       time = {:>8.1f}  +/- {:<8.1f} {:5s}"
              .format(t3s,e_t3s,t_unit))
        print("                  alt. = {:>8.1f}  +/- {:<8.1f} {:5s}"
              .format(alt3s,e_alt3s,alt_unit))
        print("                   az. = {:>8.1f}  +/- {:<8.1f} {:5s}"
              .format(az3s,e_az3s, alt_unit))
        print("            Det. level = {:>8.2f} %"
              .format(100*mc.detect_3s/mc.niter))
        print("            Det. level = {:>8.2f} %"
              .format(100*len(mc.id_3s_list)/mc.niter))
        print("           S, B counts = {:>8.1f}, {:<8.1f}"
              .format(nex3s,nb3s))
        print("        on, off counts = {:>8.1f}, {:<8.1f}"
              .format(nex3s+nb3s,nb3s/mcf.alpha))        
        print()
    
        print("  5 sigma :       time = {:>8.1f}  +/- {:<8.1f} {:5s}"
              .format(t5s,e_t5s,t_unit))
        print("                  alt. = {:>8.1f}  +/- {:<8.1f} {:5s}"
              .format(alt5s,e_alt5s,alt_unit))
        print("                   az. = {:>8.1f}  +/- {:<8.1f} {:5s}"
              .format(az5s,e_az5s, alt_unit))
        print("            Det. level = {:>8.2f} %"
              .format(100*mc.detect_5s/mc.niter))
        print("           S, B counts = {:>8.1f}, {:<8.1f}"
              .format(nex5s,nb5s))
        print("        on, off counts = {:>8.1f}, {:<8.1f}"
              .format(nex5s+nb5s,nb5s/mcf.alpha))
     
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++") 
        print("{:>}".format("SoHAPPy!"))
    # Dump result into population file - even if simulation was not
    # performed or aborted
    

    if ((write_header==True) and (popfile != None)): output_header(popfile)

    txt = "{:>10s} {:8.2e} {:5.2f} "  # GRB, Eiso, z
    txt+= "{:>5s} {:>6.2f} {:>6.2f} " # IRF : site, Ra, Dec
    txt+= "{:>15.8f} "                # Trigger time
    txt+= "{:>9.2f} {:>9.2f} "        # tstart, tstop, (ref:trigger date)
    txt+= "{:>5.2f} {:>5.2f} "        # alt-start, alt-stop
    txt+= "{:>6.2f} {:>6.2f} "        # az-start, az-stop, 
    txt+= "{:>3d} {:>3d} "            # nbins, ndof
    txt+= "{:>7.1f} {:>7.1f} {:>10.1f} {:>10.1f} {:>10.1f} {:>10.1f} " # Max
    txt+= "{:>9.2f} {:>9.2f} {:>5.2f} {:>6.2f} {:>6.2f} {:>6.2f} "     # Max
    txt+= "{:>8.1f} {:>8.1f} {:>8.1f} {:>8.1f} " #3s
    txt+= "{:>9.2f} {:>9.2f} {:>5.2f} {:>6.2f} {:>6.2f} {:>6.2f} " # 3s
    txt+= "{:>4d} " # 3s
    txt+= "{:>8.1f} {:>8.1f} {:>8.1f} {:>8.1f} " #5s
    txt+= "{:>9.2f} {:>9.2f} {:>5.2f} {:>6.2f} {:>6.2f} {:>6.2f} " # 5s
    txt+= "{:>4d} " # 5s
    txt+= "{:>5.2f} {:>4d} "
    txt+= "{}"

    if (popfile != None):
        print(txt.format(
              mc.grb.name, mc.grb.Eiso.value,     mc.grb.z,
              mc.where,    mc.grb.radec.ra.value, mc.grb.radec.dec.value,
              mc.grb.t_trig.mjd,
              mc.tstart.value, mc.tstop.value,
              mc.pos_start.alt.value, mc.pos_stop.alt.value,
              mc.pos_start.az.value, mc.pos_stop.az.value,
              mc.ntbins,mc.ndof,
              smax_mean, e_smax,  nex_smax, e_nex_smax, nb_smax, e_nb_smax, 
              tmx,       e_tmx,   altmx,    e_altmx,    azmx,    e_azmx, 
              nex3s,     e_nex3s, nb3s,     e_nb3s,
              t3s,       e_t3s,   alt3s,    e_alt3s,    az3s,    e_az3s,  
              mc.detect_3s,
              nex5s,     e_nex5s, nb5s,     e_nb5s,
              t5s,       e_t5s,   alt5s,    e_alt5s,    az5s,    e_az5s,
              mc.detect_5s,
              mc.mctime, mc.err,
              '"'+str(mc.err_slice)+'"'), # A trick to store it as one element
              file = popfile)

    return

###############################################################################      
def output_header(popfile):
    """
    Keywords 
    name  : source name (10s) 
    eiso  : Eiso (8.2e)
    z     : Redshift (5.2f)
    site  : South, North, Both (5s)
    ra    : Right ascension of the source (6.2f)
    dec   : Declination of the source (6.2f)
    ttrig : Trigger date (mjd) (15s)
       
    t1    : Observ. start since trigger (s) (9.2f) (incl. slewing time)
    t2    : Observ. stop since trigger (s) (9.2f)
    alt1  : Altitude at t1 (deg) (5.2f)
    alt2  : Altitude at t2 (deg) (5.2f)
    az1   : Azimuth at t1 (deg) (6.2f)
    az2   : Azimuth at t2 (deg) (6.2f)
    
    nt    : Number of time bins in analysis (3d)
    dof   : Number of degree of freedom in analysis (0=on-off) (3d)
       
    sigmx : Mean maximum significance reached (7.1f)
    nexmx : Mean number of excess counts at sigmx (10.1f)
    enexmx: Error on mean number of excess counts at sigmx (10.1f)
    nbmx  : Mean number of background counts at sigmx (10.1f)
    enbmx : Error on mean number of background counts at sigmx (10.1f)
    tmx   : Mean time at sigmx (discrete, related tot time binning) (9.2f)
    altmx : Mean altitude at sigmx (deg) (5.2f)
    azmx  : Mean azimuth at sigmx (deg) (6.2f)
    etmx  : Error on mean time at sigmx (discrete, related tot time binning) (9.2f)
    ealtmx: Error on mean altitude at sigmx (deg) (5.2f)
    eazmx : Error on mean azimuth at sigmx (deg) (6.2f)    
    
    nex3s : Mean number of excess counts at 3 sigma (8.1f)
    enex3s: Error on mean number of excess counts at 3 sigma (8.1f)
    nb3s  : Mean number of background counts at 3 sigma (8.1f)
    enb3s : Error on mean number of background counts at 3 sigma (8.1f)
    t3s   : Mean time at 3 sigma (discrete, related tot time binning) (9.2f)
    alt3s : Mean altitude at 3 sigma (deg) (5.2f)
    az3s  : Mean azimuth at 3 sigma  (deg) (6.2f)
    et3s  : Error on mean time at 3 sigma (related to time binning) (9.2f)
    ealt3s: Error on mean altitude at 3 sigma (deg) (5.2f)
    eaz3s : Error on mean azimuth at 3 sigma  (deg) (6.2f)    
    d3s   : Number of trials above 3 sigma (4d)
    
    nex5s : Mean number of excess counts at 5 sigma (8.1f)
    enex5s: Error on mean number of excess counts at 5 sigma (8.1f)
    nb5s  : Mean number of background counts at 5 sigma (8.1f)
    enb5s : Error on mean number of background counts at 5 sigma (8.1f)
    t5s   : Mean time at 5 sigma (discrete, related tot time binning) (9.2f)
    alt5s : Mean altitude at 5 sigma (deg) (5.2f)
    az5s  : Mean azimuth at 5 sigma  (deg) (6.2f)
    et5s  : Error on nean time at 5 sigma (related to time binning) (9.2f)
    ealt5s: Error on mean altitude at 5 sigma (deg) (5.2f)
    eaz5s : Error on mean azimuth at 5 sigma  (deg) (6.2f)    
    d5s   : Number of trials above 5 sigma (4d)
       
    mct   : computing time per trial (s) (5.2f)    
    err   : general error code (4d)
    err_slice : list of slices for which an error or a warning is reported

    Error codes :
    err : 
        If positive, give the time slice at which the simulation was stopped
        A complete simulation has err = number of trials requested
        If negative :
         -1n    : simulations completed but errors reported in n trials  
         -999   : simulation not done (GRB not visible)
    """
    
    txt = "{:>10} {:>8} {:>5} " # GRB, Eiso, z
    txt+= "{:>5} {:>6} {:>6} "  # site, Ra, Dec
    txt+= "{:>15} "             # Trigger date
    txt+= "{:>9} {:>9} "        # tstart, tstop,  (ref:trigger date)
    txt+= "{:>5} {:>5} "        # alt-start, alt-stop
    txt+= "{:>6} {:>6} "        # az-start, az-stop, 
    txt+= "{:>3} {:>3} "        # nbins, ndof
    txt+= "{:>7} {:>7} {:>10} {:>10} {:>10} {:>10} "   # sigmax nex_max nb_max 
    txt+= "{:>9} {:>9} {:>5} {:>6} {:>6} {:>6} " # tmax  alt_max az_max 
    txt+= "{:>8} {:>8} {:>8} {:>8} " # nex3s, nb3s
    txt+= "{:>9} {:>9} {:>5} {:>6} {:>6} {:>6} " # t3s, alt3s, az3s
    txt+= "{:>4} " # detect3s
    txt+= "{:>8} {:>8} {:>8} {:>8} " # nex_5s, nb_5s
    txt+= "{:>9} {:>9} {:>5} {:>6} {:>6} {:>6} " # t_5s, alt_5s, az_5s
    txt+= "{:>4} " #  detect_5s
    txt+= "{:>5} {:>4} " # computing time for this GRB
    txt+= "{}"
    print(txt.format("name" , "eiso"  , "z"     ,
                     "site" , "ra"    , "dec"   ,
                     "ttrig",
                     "t1"   , "t2"    , 
                     "alt1" , "alt2"  ,
                     "az1"  , "az2"   ,
                     "nt"   , "dof"   ,
                     "sigmx", "esigmx", "nexmx" , "enexmx", "nbmx"  , "enbmx",
                     "tmx"  , "etmx"  , "altmx" , "ealtmx", "azmx" , "eazmx",
                     "nex3s", "enex3s", "nb3s"  , "enb3s" ,
                     "t3s"  , "et3s"  , "alt3s" , "ealt3s", "az3s" , "eaz3s",
                     "d3s"  ,  
                     "nex5s", "enex5s", "nb5s"  , "enb5s" ,   
                     "t5s"  , "et5s"  , "alt5s" , "ealt5s", "az5s" , "eaz5s",
                     "d5s"  ,  
                     "mct"  ,   "err" ,
                     "err_slice"),
                      file = popfile)
    return