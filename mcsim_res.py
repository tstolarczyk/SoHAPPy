# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""
import numpy as np
import astropy.units as u
import mcsim_config as mcf
#from utilities import warning, failure, success, highlight,banner

###########################################################################
def welcome(subarray,log=None):
    log.banner("+================================================================+")
    log.banner("||                     LAUNCHING SIMULATION                     ||")
    log.banner("+================================================================+")
    log.prt("   Detection level    : {}".format(mcf.det_level))
    log.prt("   on/off regions     : {}".format(mcf.alpha))
    log.prt("   FOV                : {}".format(mcf.fov))
    log.prt("   On-region size     : N: {:5} -  S: {:5}".format(mcf.on_size[subarray["North"]],
                                                       mcf.on_size[subarray["South"]]))
    log.prt("   Offset from center : N: {:5} -  S: {:5}".format(mcf.offset[subarray["North"]],
                                                       mcf.offset[subarray["South"]]))
    log.prt("   Eff. area cont.    : {}".format(mcf.containment))
    log.prt("   Min on/off counts  : {}".format(mcf.nLiMamin))
    #log.prt("   Bin size           : {}".format(mcf.binsize))

    return
###########################################################################
def mean_values(mc, idlist):

    """
    Extract mean values and RMS from the index list
    """
    t_det, e_t_det = 2*(-1*u.s,)
    alt_det, e_alt_det, az_det, e_az_det  = 4*(-1*u.deg,)

    if (len(idlist)):
        t_unit = mc.slot.slices[0].tobs().unit

        t = [mc.slot.slices[i].tobs().value for i in idlist
             if mc.slot.slices[i].tobs().value>=0]
        t_det   = np.mean(t) *t_unit
        e_t_det = np.std(t)  *t_unit

        if (mc.slot.site != "Both"):
            altaz =  mc.slot.grb.altaz(dt=t*t_unit,loc=mc.slot.site)

            alt = altaz.alt
            alt_det   = np.mean(alt)
            e_alt_det = np.std(alt)

            az = altaz.az
            az_det   = np.mean(az)
            e_az_det = np.std(az)

    return  t_det, e_t_det, alt_det, e_alt_det, az_det, e_az_det

###########################################################################
def result(mc,grb, log=None, header=True,pop=None):

    """
    Print out results.
    If only one GRB is simulated, store the present mc class on disks
    for further investigation.

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

    ###----------------------------
    ### General status message
    ###----------------------------
    if  (mc.err < 0):
        log.prt("+--------------+-------------------------------------------------+")
        log.prt("| Status       |  ===> ",end="")
        log.failure("Simulation not possible (GRB not visible)",end="")
        log.prt(" |")
        log.prt("+--------------+-------------------------------------------------+")
    if (mc.err> 0 and mc.err != mc.niter):
        log.prt("\n+--------------+-------------------------------------------------+")
        log.prt("| Status       |  ===> ",end="")
        log.failure("Simulation was aborted at trial {:>3d}".format(mc.err),end="")
        log.prt("       |")
        log.prt("+--------------+-------------------------------------------------+")
    if (mc.err == mc.niter):
        if (mc.detect_5s/mc.niter >= mcf.det_level):
            message = "5 sigma detected"
        elif (mc.detect_3s/mc.niter >= mcf.det_level):
            message = "3 sigma detected"
        else:
            message = "NOT detected    "
        log.prt("\n+--------------+-------------------------------------------------+")
        log.prt("| Status       |  ===> ",end="")
        log.highlight(message,end="")
        log.prt("                          |")
        log.prt("+--------------+-------------------------------------------------+")

    ###-----------------------------------------
    ### Compute results if everything went well
    ###-----------------------------------------
    tstart    = -1*u.s
    tstop     = -1*u.s
    ntbins    = -1
    alt_start = -1*u.deg
    alt_stop  = -1*u.deg
    az_start  = -1*u.deg
    az_stop   = -1*u.deg

    smax_mean, e_smax, nex_smax, e_nex_smax, nb_smax, e_nb_smax= 6*(-1,)
    tmx, e_tmx = 2*(-1*u.s,)
    altmx, e_altmx, azmx, e_azmx = 4*(-1*u.deg,)

    nex3s, e_nex3s, nb3s, e_nb3s = 4*(-1,)
    t3s,e_t3s                    = 2*(-1*u.s,)
    alt3s,e_alt3s,az3s, e_az3s   = 4*(-1*u.deg,)

    nex5s, e_nex5s, nb5s, e_nb5s = 4*(-1,)
    t5s,e_t5s                    = 2*(-1*u.s,)
    alt5s,e_alt5s,az5s, e_az5s   = 4*(-1*u.deg,)

    if (mc.err == mc.niter):
        tstart = mc.slot.tstart
        tstop  = mc.slot.tstop
        ntbins = len(mc.slot.slices)

        loc = mc.slot.site
        if (loc == "North" or loc == "South"):
            altaz = [mc.slot.grb.altaz(dt=mc.slot.tstart,
                               loc=loc),
             mc.slot.grb.altaz(dt=mc.slot.tstop,
                               loc=loc)]
            alt_start = altaz[0].alt
            alt_stop = altaz[1].alt
            az_start = altaz[0].az
            az_stop = altaz[1].az

        #log.prt("+--------------+-------------------------------------------------+")
        log.highlight("  RESULTS      :              {:<10s}"
                  .format(mc.name))
        #log.prt("+--------------+-------------------------------------------------+")
        log.prt("+----------------------------------------------------------------+")
        log.prt(" Window : {:6.2f} - {:6.2f} * Delay: {:6.2f} * {:3d} slices"
              .format(tstart,tstop,mc.slot.delay,ntbins))

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

        ###----------------------------
        ### Dump results
        ###----------------------------
        log.prt("+----------------------------------------------------------------+")
        log.prt("  Max. sigma           = {:>8.1f}  +/- {:<8.1f}"
              .format(smax_mean,e_smax))

        log.prt("                  time = {:>8.1f}  +/- {:<8.1f} {:5s}"
               .format(tmx.value,e_tmx.value,tmx.unit))
        log.prt("                  alt. = {:>8.1f}  +/- {:<8.1f} {:5s}"
               .format(altmx.value,e_altmx.value,altmx.unit))
        log.prt("                   az. = {:>8.1f}  +/- {:<8.1f} {:5s}"
               .format(azmx.value,e_azmx.value,azmx.unit))
        log.prt("           S, B counts = {:>8.1f} {:0.1f}"
              .format(nex_smax,nb_smax))
        log.prt("        on, off counts = {:>8.1f} {:0.1f}"
              .format(nex_smax+nb_smax,nb_smax/mcf.alpha))

        log.prt("+----------------------------------------------------------------+")
        log.prt("                           3 sig                   5 sig")
        log.prt("+----------------------------------------------------------------+")
        log.prt(" Time ({:3s})     = {:>8.1f}  +/- {:<8.1f}  {:>8.1f}  +/- {:<8.1f}"
              .format(t3s.unit, t3s.value, e_t3s.value,
                                t5s.value , e_t5s.value))
        log.prt(" Alt.({:3s})      = {:>8.1f}  +/- {:<8.1f}  {:>8.1f}  +/- {:<8.1f}"
              .format(alt3s.unit,alt3s.value, e_alt3s.value,
                                 alt5s.value, e_alt5s.value))
        log.prt(" Az.({:3s})       = {:>8.1f}  +/- {:<8.1f}  {:>8.1f}  +/- {:<8.1f}"
              .format(az3s.unit , az3s.value, e_az3s.value,
                                  az5s.value, e_az5s.value))
        log.prt(" Det. level (%) =  {:>8.2f}                {:>8.2f}"
              .format(100*mc.detect_3s/mc.niter,100*mc.detect_5s/mc.niter))
        log.prt(" S, B           = {:>8.1f}, {:<8.1f}      {:>8.1f}, {:<8.1f}"
              .format(nex3s,nb3s,nex5s,nb5s))
        log.prt(" On, Off        = {:>8.1f}, {:<8.1f}      {:>8.1f}, {:<8.1f}"
              .format(nex3s+nb3s,nb3s/mcf.alpha,nex5s+nb5s,nb5s/mcf.alpha))
        log.prt("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        log.prt("{:>}".format("SoHAPPy!"))

    ###----------------------------
    ### Define dictionnary
    ###----------------------------
    # Note : grb = mc.slot.grb would avoid having grb passed in the function
    # but unvisble have no slot defined
    # Solution : do not print out unviible GRB ? But then output data file will
    #  have mmissing lines for some GRB
    records = {
            "name"     : {"v": grb.name,                      "f": ">10s"  },
            "Eiso"     : {"v": grb.Eiso.value,                "f": ">8.2e"  },
            "z"        : {"v": grb.z,                         "f": ">5.2f"  },
            "site"     : {"v": mc.name[mc.name.find("-")+1:], "f": ">5s"  },
            "ra"       : {"v": grb.radec.ra.value,  "f": ">6.2f" },
            "dec"      : {"v": grb.radec.dec.value, "f": ">6.2f" },
            "ttrig"    : {"v": grb.t_trig.mjd,      "f": ">15.8f"},
            "t1"       : {"v": tstart.value,      "f": ">9.2f" },
            "t2"       : {"v": tstop.value,       "f": ">9.2f" },
            "alt1"     : {"v": alt_start.value,   "f": ">5.2f" },
            "alt2"     : {"v": alt_stop.value,    "f": ">5.2f" },
            "az1"      : {"v": az_start.value,    "f": ">6.2f" },
            "az2"      : {"v": az_stop.value,     "f": ">6.2f" },
            "nt"       : {"v": ntbins,            "f": ">3d"   },
            "dof"      : {"v": mc.method,    "f": ">3d"   },
            "sigmx"    : {"v": smax_mean,  "f": ">7.1f" },
            "esigmx"   : {"v": e_smax,     "f": ">7.1f" },
            "nexmx"    : {"v": nex_smax,   "f": ">10.1f"},
            "enexmx"   : {"v": e_nex_smax, "f": ">10.1f"},
            "nbmx"     : {"v": nb_smax,    "f": ">10.1f"},
            "enbmx"    : {"v": e_nb_smax,  "f": ">10.1f"},
            "tmx"      : {"v": tmx.value,        "f": ">9.2f" },
            "etmx"     : {"v": e_tmx.value,      "f": ">9.2f" },
            "altmx"    : {"v": altmx.value,      "f": ">5.2f" },
            "ealtmx"   : {"v": e_altmx.value,    "f": ">6.2f" },
            "azmx"     : {"v": azmx.value,       "f": ">6.2f" },
            "eazmx"    : {"v": e_azmx.value,     "f": ">6.2f" },
            "nex3s"    : {"v": nex3s,         "f": ">8.1f" },
            "enex3s"   : {"v": e_nex3s,       "f": ">8.1f" },
            "nb3s"     : {"v": nb3s,          "f": ">8.1f" },
            "enb3s"    : {"v": e_nb3s,        "f": ">8.1f" },
            "t3s"      : {"v": t3s.value,           "f": ">9.2f" },
            "et3s"     : {"v": e_t3s.value,         "f": ">9.2f" },
            "alt3s"    : {"v": alt3s.value,         "f": ">5.2f" },
            "ealt3s"   : {"v": e_alt3s.value,       "f": ">6.2f" },
            "az3s"     : {"v": az3s.value,          "f": ">6.2f" },
            "eaz3s"    : {"v": e_az3s.value,        "f": ">6.2f" },
            "d3s"      : {"v": mc.detect_3s,  "f": ">4d"   },
            "nex5s"    : {"v": nex5s,         "f": ">8.1f" },
            "enex5s"   : {"v": e_nex5s,       "f": ">8.1f" },
            "nb5s"     : {"v": nb5s,          "f": ">8.1f" },
            "enb5s"    : {"v": e_nb5s,        "f": ">8.1f" },
            "t5s"      : {"v": t5s.value,           "f": ">9.2f" },
            "et5s"     : {"v": e_t5s.value,         "f": ">9.2f" },
            "alt5s"    : {"v": alt5s.value,         "f": ">5.2f" },
            "ealt5s"   : {"v": e_alt5s.value,       "f": ">6.2f" },
            "az5s"     : {"v": az5s.value,          "f": ">6.2f" },
            "eaz5s"    : {"v": e_az5s.value,        "f": ">6.2f" },
            "d5s"      : {"v": mc.detect_5s,              "f": ">4d"   },
            "mct"      : {"v": mc.mctime,                 "f": ">5.2f" },
            "err"      : {"v": mc.err,                    "f": ">4d"   },
            "err_slice": {"v": '"'+str(mc.err_slice)+'"', "f": ""      },
    }

    ###--------------------------------------
    ### Dump dictionnary into population file
    ###--------------------------------------
    # Even if simulation was not performed or aborted
    if (pop != None):

        if (header == True):
            # Print header the first time
            for r in records:
                x = records[r]
                # Get format from the dictionnary
                f = x["f"][:-1] # remove type
                pos = f.find(".") # Get dot position
                if (pos>= 0): f = f[:pos]
                f = f +"s"

                print("{:{}}".format(r,f),end=" ",file=pop)
            print(file=pop)


            header = False
        for r in records:
            x = records[r]
            print("{:{}}".format(x["v"],x["f"]),end=" ",file=pop)
        print(file=pop)

    return