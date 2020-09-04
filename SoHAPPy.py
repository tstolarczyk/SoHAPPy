"""
* Create a GRB list to be analysed from the config file or the command line (:func:`init`)


"""
# from . import __version__ # does not work
__version__ = "Sofia dev"

import os
import sys, getopt

import numpy as np
import time, datetime
import pickle

from   pathlib import Path
import astropy.units as u

import ana_config as cf # Steering parameters

# Transform warnings into errors - useful to find who is guilty !
# import warnings
# warnings.filterwarnings('error')

import gammapy # Just for the version number

from grb            import GammaRayBurst
from timeslot       import Slot
from   mcsim        import MonteCarlo
import mcsim_res    as mcres

from utilities import backup_file, Log, warning
import visibility as v

os.environ['GAMMAPY_EXTRA'] =r'../input/gammapy-extra-master'
os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
#print(os.getenv('GAMMAPY_EXTRA'))
#print(os.listdir(os.getenv('GAMMAPY_EXTRA')))

__all__ = ["main", "init", "summary", "get_grb_fromfile", "get_delay"]
###############################################################################
def get_grb_fromfile(i,prompt=False, afterglow=False, log=None):
    """
    Obtain data for the ith GRB file and create a GammaRayBurst instance.

    Parameters
    ----------
    i : integer
        GRB position in the list
    prompt : boolen, optional
        If True, reads the prompt GRB file, otherwise the afterglow (default).
    afterglow : boolean, optional
        If True, use the afterglow file generic information (e.g.) for the
        prompt data, oterwise use the default.
        The default is False.
    log : TextIO, optional
        Log file. The default is None.

    Returns
    -------
    grb : GammaRayBurst
        A GammaRayBurst instance

    """

    if (cf.test_prompt):
        # create a new object from the default (Visible in North)
        loc = Path('../input/lightcurves/prompt'
                   + "/events_"+str(i)+".fits")
        if (cf.use_afterglow):
            # use afterglow characteristics
            loc_glow = Path(cf.grb_dir + "/Event"+str(i)+".fits")
            glow = GammaRayBurst.read(loc_glow, ebl = cf.EBLmodel)
            grb = GammaRayBurst.read_prompt(loc, glow=glow, ebl = cf.EBLmodel)
        else:
            # use default visibility
            grb = GammaRayBurst.read_prompt(loc,
                                            glow=None,
                                            ebl = cf.EBLmodel,
                                            z=cf.redshift)
    else:
        if (cf.old_file):
            loc = Path(cf.grb_olddir + "LGRB"+str(i))
            grb = GammaRayBurst.read_old(loc,ebl = cf.EBLmodel)
        else:
            loc = Path(cf.grb_dir + "/Event"+str(i)+".fits")
            grb = GammaRayBurst.read(loc, ebl = cf.EBLmodel)

    return grb

###############################################################################
def init(argv):
    """
    Decode the command line arguments if any, treat some information from the
    configuration file, overwrite some configuration parameters with the
    command line parameters.
    Build the debugging flag for printout and plots.
    Create the result output folder.

    Parameters
    ----------
    argv : Command line arguments
        User command line arguments.

    Returns
    -------
    None.

    """

    # Read arguments from command line, supersede default ones
    try:
          opts, args = getopt.getopt(argv,
                                     "hn:f:d:N:o:",
                                     ["ngrb=",
                                      "first=",
                                      "niter=",
                                      "dbg=",
                                      "output=",])

          for opt, arg in opts:
             if opt == '-h':
                 print(" SoHAPPy.py "
                       + "-N <ngrb> "
                       + "-f <1st grb or list> "
                       + "-n <MC iterations> "
                       + "-o <Output folder> "
                       + "-d <debug> ")
                 sys.exit()
             elif opt in ("-N", "--ngrb"):
                 cf.ngrb =  int(arg)
             elif opt in ("-f", "--first"):
                 cf.ifirst = int(arg)
             elif opt in ("-n", "--niter"):
                 cf.niter = int(arg)
             elif opt in ("-o", "--output"):
                  cf.res_dir = arg
             elif opt in ("-d", "--debg"):
                 dbg = int(arg)
                 cf.dbg = abs(dbg)

    except getopt.GetoptError:
        print("No line arguments passed: Using default values")

    # Create the show debugging flag from the general debug flag
    if (cf.dbg < 0): cf.show = 0
    else: cf.show = abs(cf.dbg)

    if (cf.do_fluctuate == False): cf.niter = 1
    if (cf.niter == 1): cf.do_fluctuate=False
    if (cf.dbg>0): cf.silent = False


    # Check that the output folder exist, otherwise create it
    if (cf.res_dir[-1] != "/"): cf.res_dir = cf.res_dir+"/"
    if not os.path.isdir(cf.res_dir):
        warning("Creating {}".format(cf.res_dir))
        os.mkdir(cf.res_dir)

    return

###############################################################################
def summary(log=None):
    """
    Printout the main characteristics of the simulation.

    Parameters
    ----------
    log : TextIO, optional
        Log file. The default is None.

    Returns
    -------
    None.

    """

    log.prt("")
    log.prt("+----------------------------------------------------------------+")
    log.prt("|                                                                |")
    log.prt("|                    SoHAPPy with GammaPy {:4s}                   |"
          .format(gammapy.__version__))
    log.prt("|                            ({:4s})                              |"
          .format(__version__))
    log.prt("|  (Simulation of High-energy Astrophysics Processes in Python)  |")
    log.prt("|                                                                |")
    log.prt("+----------------------------------------------------------------+")
    log.prt(" Simulation:")
    if type(cf.ifirst)!=list:
        log.prt("     *Number of sources  : {:>5d}".format(cf.ngrb))
        log.prt("     *First source       : {:>5d}".format(cf.ifirst))
    else:
        log.prt("     * Source list       : {}".format(cf.ifirst))
    log.prt("     *Number of trials   : {:>5d}".format(cf.niter))
    log.prt(" EBL model               : {}".format(cf.EBLmodel))
    log.prt(" Input/output :")
    log.prt("     *Debug mode         : {:>5d}".format(cf.dbg))
    log.prt("     *Show plots         : {:>5d}".format(cf.show))
    log.prt("      Analysing files in : {}".format(cf.grb_dir))
    log.prt("      IRF files in       : {}".format(cf.irf_dir))
    log.prt("     *Result folder      : {}".format(cf.res_dir))
    if (cf.lkhd ==0):
        method = "On-Off"
    else:
        method = "{}D likelihood".format(cf.lkhd)
    log.prt(" Analysis (ndof)         : {}".format(method))
    if (cf.newvis):
        log.prt(" Vis. computed up to : {}".format(cf.depth))
    else:
        log.prt(" Visibility          : default")
    log.prt(" Skip up to night        : {}".format(cf.skip))

    log.prt("+----------------------------------------------------------------+")
    log.prt("|                 *: can be changed with command line (use -h)   |")
    log.prt("+----------------------------------------------------------------+")
    log.prt(" Developments:")
    if (cf.save_grb == True):
        log.highlight("Simulation saved to disk save_grb (save_grb = True)       ")
    if (cf.write_slices == True):
        log.highlight("Slice information saved to disk (write_slices=True)       ")
    if (cf.signal_to_zero == True):
        log.warning(  "Signal set to zero (signal_to_zero==True)                 ")
    if (cf.do_fluctuate == False):
        log.warning(  "No fluctuation in simulation (do_fluctuate==False)        ")
    if (cf.do_accelerate  == False):
        log.warning(  "No abortion if first 10% undetected (do_accelarate==False)")
    if (cf.fixed_zenith != False):
        log.warning(  "Zenith angle requested to be fixed at keyword '{:5s}'     "
               .format(cf.fixed_zenith))
    if (cf.niter == 0):
        log.failure(  " Cannot run simulation with ZERO trials")
        log.warning(  " Use other main specific scripts for tests")
        sys.exit( " At least one trial is requested")
    if (cf.test_prompt):
        log.warning(  "Test prompt simulation")

    log.prt("")

    return

###############################################################################
def get_delay():
    """
    Compute the overall delay to be applied to the start of detection
    (satelite and telescope slewing), according to the user parameters

    Returns
    -------
    dt : Quantity (time)
        Delay before the detection can start.

    """

    dt = 0*u.s
    if (cf.fixslew):  dt = cf.dtslew
    else:             dt = cf.dtslew*np.random.random()

    if (cf.fixswift): dt = dt + cf.dtswift # don't do += !!!
    else: sys.exit("Variable SWIFT delay not implemented)")

    return dt.to(u.s)

###############################################################################
def main(argv):
    """
    This is the main and it should be seriously documented

    Parameters
    ----------
    argv : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    init(argv)           # Steering parametets and welcome message
    if (cf.ngrb<=0):
        print(" NO ANALYSIS REQUIRED (ngrb<=0)")
        sys.exit(2)

    start_all = time.time() # Starts chronometer

    # GRB list to be analysed
    if type(cf.ifirst)!=list:
        grblist = list(range(cf.ifirst,cf.ifirst+cf.ngrb))
        first = str(cf.ifirst)
    else:
        grblist = cf.ifirst
        first = str(grblist[0])

    # Create population simulation and logfile name
    sim_filename    = cf.res_dir  + cf.datafile
    log_filename    = cf.res_dir  + cf.logfile
    log = Log(name  = log_filename, talk=not cf.silent)
    conf_filename   = "config.py"
    conf_filename   = backup_file(folder=cf.res_dir,dest=conf_filename)

    # Print Summary
    summary(log=log)

    start_pop = time.time() # Start chronometer

    with open(sim_filename, 'w') as pop:

        ##############################
        # Loop over GRB population   #
        ##############################

        mcres.welcome(log=log) # Remind simulation parameters

        first = True # Actions for first GRB only
        for i in grblist:

            grb = get_grb_fromfile(i,log=log) ### Get GRB

            # Create original slot (slices) and fix observation points
            origin = Slot(grb,
                          opt   = cf.obs_point,
                          name  = grb.name,
                          debug = bool(cf.dbg>1))

            # Recompute visbility windows
            if (cf.newvis):
                vis = v.Visibility(grb,loc="North") # Constructor
                vis.compute(altmin = cf.altmin,
                            depth  = cf.depth,
                            skip   = cf.skip)
                grb.update_visibility(vis)

                vis = v.Visibility(grb,loc="South") # Constructor
                vis.compute(altmin = cf.altmin,
                            depth  = cf.depth,
                            skip   = cf.skip)
                grb.update_visibility(vis)

            # Printout grb and visibility windows
            if cf.niter<=1 or cf.dbg>0 or cf.ngrb==1 :
                log.prt(grb)
                grb.print_visibility(loc="North",log=log)
                grb.print_visibility(loc="South",log=log)

            # Plot grb spectra and lightcurve and visibility windows
            if (cf.show >0):
                import grb_plot     as gplt
                gplt.spectra(grb,opt="Packed")
                gplt.visibility_plot(grb,loc="North")
                gplt.visibility_plot(grb,loc="South")

                gplt.pause()

                if (cf.dbg>2) :
                    gplt.animated_spectra(grb,savefig=True,outdir=cf.res_dir)

    #            plt.show(block=True)

            # Save GRB to file if requested
            if (cf.save_grb):
                grb_class_file = cf.res_dir + "/" \
                               + grb.name + ".bin"
                outfile  = open(grb_class_file,"wb")
                pickle.dump(grb,outfile)
                outfile.close()
                print(" Saving grb {} to file : {}"
                      .format(grb.name,grb_class_file))

            ###--------------------------------------------###
            #  Check individual sites - Loop over locations
            ###--------------------------------------------###
            for loc in grb.site_keys:

                name = grb.name + "-" + loc

                # Create a MC object
                log.banner(" SIMULATION  : {:<50s} ".format(name))
                mc = MonteCarlo(niter  = cf.niter,
                                ndof   = cf.lkhd,
                                debug  = cf.dbg,
                                name   = name)

                # If visible, run simulation
                if grb.vis_tonight[loc]:

                    slot = origin.copy(name="loc")
                    slot.apply_visibility(delay     = get_delay(),
                                          site      = loc)
                    mc.run(slot)

                # Get information and results even if not visible
                first = mcres.result(mc, grb, log=log, header=first, pop=pop)

                # If Simulation was not aborted, plot some results
                if (mc.err == mc.niter) and (cf.show > 0):
                    slot.plot()
                    import mcsim_plot as mplt
                    mplt.show(mc,loc=loc)

                # If requested save simulation to disk
                if (cf.save_simu):
                    mc.save_to_disk(cf.res_dir + "/" +name + ".bin")

            ###--------------------------------------------###
            #   Check GRB seen on both sites
            ###--------------------------------------------###
            name = grb.name + "-Both"

            # Create a MC object
            log.banner(" SIMULATION  : {:<50s} ".format(name))
            mc = MonteCarlo(niter  = cf.niter,
                            ndof   = cf.lkhd,
                            debug  = cf.dbg,
                            name   = name)

            # If visible, run simulation
            if grb.vis_tonight["North"] and grb.vis_tonight["South"]:

                slot = origin.both_sites(delay   = get_delay(),
                                           debug =(cf.dbg>1))
                mc.run(slot)

            # Get information and results even if not visible
            first= mcres.result(mc, grb, log=log, header=first, pop = pop)

            # If Simulation was not aborted, plot some results
            if (mc.err == mc.niter) and (cf.show > 0):
                if (cf.show>0): slot.plot()
                import mcsim_plot as mplt
                mplt.show(mc,loc="Both")

            # If requested save simulation to disk
            if (cf.save_simu):
                mc.save_to_disk(cf.res_dir + "/" +name + ".bin")

        # END of Loop over GRB

    # Stop chronometer
    end_pop = time.time()
    end_all = time.time()
    elapsed = end_pop-start_pop
    log.prt("\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    log.prt(" Duration   = {:8.2f} (s)".format(elapsed))
    log.prt("  per GRB   = {:8.2f} (s)".format( (elapsed)/cf.ngrb))
    log.prt("  per trial = {:8.3f} (s)".format( (elapsed)/cf.ngrb/cf.niter))
    log.prt("-*-*-*-*-*-*-*- End of full population simulation -*-*-*-*-*-*-*-*\n")
    log.prt(" ******* End of job - Total time = {:8.2f} min *****"
                 .format((end_all-start_all)/60))

    # Close log file
    log.close()

    # tar gzip outputs, delete originals if requested
    nw = datetime.datetime.now()
    from pathlib import Path
    outprefix = Path(cf.res_dir).parts[-1]
    filename  = outprefix + "_" \
                                + str(nw.hour) \
                                + str(nw.minute) \
                                + str(nw.second) \
                                +".tar.gz"
    import tarfile
    tar = tarfile.open(cf.res_dir + filename, "w:gz")
    tar.add(sim_filename,arcname=os.path.basename(sim_filename))
    tar.add(log_filename,arcname=os.path.basename(log_filename))
    tar.add(conf_filename,arcname=os.path.basename(conf_filename))
    if (cf.remove_tarred):
        os.remove(sim_filename)
        os.remove(conf_filename)
        os.remove(log_filename)

    tar.close()
    print("... completed")

###############################################################################
if __name__ == "__main__":
    main(sys.argv[1:])