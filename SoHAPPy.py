"""
Create a source list to be simulated and analysed from using parameters given
in a configuration file and the command line, and fiels on disks.

Notes for experts:
* Change `warnings.filterwarnings('ignore')` into
`warnings.filterwarnings('error')` to have the code stopped in case of warning,
and be able to identify its origin.
* The 'ignore' option is motivated by warnings issued by astropy for
deprecation or too distant dates with respect to the running date.
* IERS data are not refreshed as it can take long in case of bad Internet
or no connections at all. To refresh these data use:

..  code-block::

    # For refreshing
    print(" Refreshing IERS")
    from astroplan import download_IERS_A
    download_IERS_A
    print(" ->Done")

todo:
* Remove delay when analysing not the first night (skipping the first night)
* Why applying the delay only if trigger at night?
     if grb.vis[loc].vis_night: # Apply delays to original slot

"""

import sys, os

import time
from   datetime import datetime
from   pathlib  import Path

from configuration  import Configuration
from niceprint      import Log, failure, heading, warning

from grb import GammaRayBurst
from timeslot import Slot
from mcsim  import mc_welcome, MonteCarlo
from analyze import Analysis

# Do not refresh IERS data
from astropy.utils import iers
iers.conf.auto_download = False

import warnings
warnings.filterwarnings('ignore')
# warnings.filterwarnings('error')

###############################################################################
def welcome(log):
    """
    Good luck!

    Parameters
    ----------
    log : Log object
        See :class:`Log` for details.

    """
    import gammapy
    from __init__ import __version__

    log.prt(datetime.now())
    
    log.prt(f"+{78*'-':78s}+")
    log.prt(f"+{'':^78s}+")
    log.prt(f"+{'SoHAPPy with GammaPy '+gammapy.__version__:^78s}+")
    log.prt(f"+{'('+__version__+')':^78s}+")
    log.prt(f"+{'(Simulation of High-energy Astrophysics Processes in Python)':^78s}+")
    log.prt(f"+{'':^78s}+")
    log.prt(f"+{78*'-':78s}+")

###############################################################################
def main():
    """
    The SoHAPPy main function.

    1. Manage input/output
        - Source data identifier list
        - open output simulation and log files
        - load configuration parameters

    2. Loop over input identifier list
        - Get source data from the identifiers
        - Create original time slot from the source data
        - Update visibilities in N and S if requested
        - Save source if requested
        - Check individual sites (N and S)
            - Create a MC object
            - Modify the time slot fopr the visibility including the delays
            - If GRB is still vsisible:
                - Dress the GRB slot with physics (IRF and spectra)
                - Run the simulation
            - Display results even if not visible
        - Check the N+S case if GRB visible on both sites
            - Create a MC object
            - Modify the time slot fopr the visibility including the delays
            - Dress the GRB slot with physics (IRF and spectra)
            - Run the simulation
            - Display results

    3. Close files, terminate

    Returns
    -------
    None.

    """

    # Change gammapy logging to avoid warning messages
    import logging
    logging.basicConfig()
    log = logging.getLogger("gammapy.irf")
    log.setLevel(logging.ERROR)

    ### ------------------------------------------------
    ### Configuration and output files
    ### ------------------------------------------------
    # Build the Configuration, from the defaults, a configuration file and
    # the command line arguments (sys.argv) if any.
    cf = Configuration.command_line()

    data_path = Path(cf.infolder,cf.data_dir) # Input data folder
    
    # Create output folder
    # The subfolder name Follows the convention:
    # "population name"/"user keyword"/"visibility keyword and identifiers"
    res_dir   = cf.create_output_folder()  

    # This is required to have the EBL models read from gammapy
    os.environ['GAMMAPY_DATA'] = str(Path(cf.infolder,cf.extra_dir))

    # Backup the current configuration data for further use
    # This is potentially a modified version of the local file
    conf_filename = Path(res_dir, cf.filename.name)
    cf.write(out_name = conf_filename)

    # Output file names
    sim_filename = Path(res_dir, cf.datafile) # Population file (data.txt)
    log_filename = Path(res_dir, cf.logfile)  # Log file

    # Open log file - If Silent is True, only in file, otherwise on Screen too
    log = Log(name = log_filename, talk = not cf.silent)

    # Print welcome message and configuration summary
    welcome(log)
    cf.print(log)

    # Check if something can be analysed
    if cf.nsrc <= 0:
        sys.exit(" NO ANALYSIS REQUIRED (nsrc <= 0)")

    # Prepare expert output file for individual slices
    if cf.write_slices:
        dump_dir = res_dir
    else:
        dump_dir = None

    ### ------------------------------------------------
    ### Decode visibility info
    ### ------------------------------------------------
    # visinfo contains a string or a dictionnary to be treated later
    visinfo = cf.decode_visibility_keyword(res_dir)

    ### ------------------------------------------------
    ### Start processing
    ### ------------------------------------------------
    start_pop = time.time()   # Start chronometer

    with open(sim_filename, 'w') as pop:

        #################################
        # Loop over source population   #
        #################################
        mc_welcome(cf.arrays,log=log) # Say hello, remind simulation parameters

        first = True # Actions for first GRB only
        # for i, item in enumerate(srclist):
        for item in cf.srclist:

            # If silence required, keep at least the event number for crashes
            if cf.silent: 
                print("#",item)
                
            ### Get GRB
            if isinstance(item, int):
                
                fname = Path(data_path,cf.prefix+str(item)+cf.suffix)
                if not fname.is_file():
                    failure(f" SKIPPING - File not found {fname:}")
                    continue
                
                grb = GammaRayBurst.from_fits(fname,
                                              prompt  = cf.prompt_dir,
                                              ebl     = cf.EBLmodel,
                                              Emax    = cf.Emax,
                                              dt      = cf.tshift,
                                              magnify = cf.magnify)
                # if cfg.test_prompt:
                #     return get_time_resolved_prompt_fromfile()


            elif isinstance(item, str): # this is a GRB name string
                fname    = "lightcurves/historical/GRB_"+item+".yml"
                filename = Path(cf.infolder,fname)
                grb = GammaRayBurst.historical_from_yaml(filename,
                                                         ebl = cf.EBLmodel)

            # Assign visibilities
            for loc in ["North","South"]:
                grb.set_visibility(item, loc,info=visinfo)
                # grb.limit_time_range(cf.n_night, tmin, tmax)

            # Printout grb, visibility windows, display plots
            if (cf.niter<=1 and cf.do_fluctuate==True) \
                or cf.dbg>0 or cf.nsrc==1 :
                heading(grb.id)
                log.prt(grb)
                grb.vis["North"].print(log=log)
                grb.vis["South"].print(log=log)

            if cf.save_fig and cf.show > 0: # Init. pdf output
                from matplotlib.backends.backend_pdf import PdfPages
                pdf_out = PdfPages(Path(grb.id+"_booklet.pdf"))
            else: pdf_out = None

            if cf.show > 0 : grb.plot(pdf_out)
            if cf.save_grb : grb.write_to_bin(res_dir)

            ###--------------------------------------------###
            #  Loop over locations
            ###--------------------------------------------###
            delay = cf.get_delay()

            # Create original slot (slices) and fix observation points
            origin = Slot(grb,
                          opt   = cf.obs_point,
                          name  = grb.id,
                          debug = bool(cf.dbg>1))

            for loc in ["North","South","Both"]:

                name = grb.id + "-" + loc
                if not cf.silent: # For large production, be silent
                    log.banner(" SIMULATION  : {:<50s} ".format(name))

                # Create a MC object
                # It has dummy values that will be dumped to the output
                # even if the simulation is not possible (not visible)
                mc = MonteCarlo(niter     = cf.niter,
                                fluctuate = cf.do_fluctuate,
                                nosignal  = (cf.magnify==0),
                                seed      = cf.seed,
                                debug     = cf.dbg,
                                name      = name)

                still_vis = False # Assumed not visible

                ### ------------
                ### Both sites - create a slot
                ### ------------
                if loc=="Both":
                    if grb.vis["North"].vis_night \
                    and grb.vis["South"].vis_night:
                        slot = origin.both_sites(delay = delay,
                                                 debug = (cf.dbg>1))
                        if slot != None: still_vis = True

                ### ------------
                ### North or South, create a slot
                ### ------------
                else:
                    if grb.vis[loc].vis_night: # Apply delays to original slot

                        slot = origin.copy(name="loc")
                        still_vis = slot.apply_visibility(delay = delay[loc],
                                                          site  = loc)
                ### ------------
                ### Run simulation if still visible;, prepare analysis
                ### ------------

                if still_vis:
                    # Add IRF feature and run - Note that this can
                    # modify the number of slices (merging)
                    slot.dress(irf_dir = Path(cf.infolder,cf.irf_dir),
                               arrays  = cf.arrays,
                               zenith  = cf.fixed_zenith)

                    ana = Analysis(slot, nstat = mc.niter,
                                         alpha = cf.alpha, cl = cf.det_level)

                    if cf.dbg > 1: 
                        print(slot)
                        slot.plot() 
                            
                    mc.run(slot, ana,
                                 boost     = cf.do_accelerate,
                                 dump_dir  = dump_dir)
                else:
                    # Define a default analysis for dump_to_file
                    ana = Analysis(origin, nstat = mc.niter, loc = loc)

                # If requested save simulation to disk
                if cf.save_simu: mc.write(Path(res_dir,name + "_sim.bin"))

                # Display status - even if simulation failed (not visible)
                if cf.dbg : mc.status(log=log)

                ### ------------
                ### Analyze simulated data
                ### ------------
                if ana.err == mc.niter:  # Simulation is a success
                    ana.run()
                    if cf.dbg  : ana.print(log = log)
                    if cf.show : ana.show(pdf = pdf_out)

                # # Even if not detected nor visibile, dump to file
                first = ana.dump_to_file(grb, pop, header=first)

            if cf.save_fig and cf.show>0: pdf_out.close()

            # End of loop over sites
        # END of Loop over GRB
        if cf.silent:
            print("") # Line break
        
    # Stop chronometer
    end_pop = time.time()
    elapsed = end_pop-start_pop

    log.prt("\n""-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    log.prt(" Duration   = {:8.2f} (s)".format(elapsed))
    log.prt("  per GRB   = {:8.2f} (s)".format( (elapsed)/cf.nsrc))
    log.prt("  per trial = {:8.3f} (s)".format( (elapsed)/cf.nsrc/cf.niter))
    log.prt("-*-*-*-*-*-*-*-*- End of population simulation -*-*-*-*-*-*-*-*-*\n")
    log.prt(" ******* End of job - Total time = {:8.2f} min *****"
                 .format((end_pop-start_pop)/60))
    log.prt("")
    log.prt(datetime.now())

    # Close log file
    log.close()

    # tar gzip outputs, delete originals if requested
    nw = datetime.now()
    outprefix = Path(res_dir).parts[-1]
    filename  = outprefix + "_" + nw.strftime("%Y%m%d_%H%M%S") \
                                +".tar.gz"

    import tarfile
    tar = tarfile.open(Path(res_dir,filename), "w:gz")
    tar.add(sim_filename,  arcname=os.path.basename(sim_filename))
    tar.add(log_filename,  arcname=os.path.basename(log_filename))
    tar.add(conf_filename, arcname=os.path.basename(conf_filename))

    if cf.remove_tar:
        os.remove(sim_filename)

        os.remove(Path(res_dir,cf.filename.name)) # Remove the copy, not the original !

        # After CTRL-C in Spyder, or when the code crashes, the log file
        # cannot be removed (although it was possible to overwrite it)
        # It seems to be Windows specific
        if not log.log_file.closed: log.log_file.close()
        try:
            os.remove(log_filename)
        except IOError:
            print("{} removal failed: locked".format(log_filename))

    tar.close()
    print("... completed")

###############################################################################
if __name__ == "__main__":
    
    if len(sys.argv[1:]) <= 0:
        print("------------------> Execute examples")
        # sys.argv=["", "-c","myConfigs/config-LongFinalTest-omega.yaml"]
        # sys.argv=["", "-c","myConfigs/config_Long1000_strictmoonveto_1.yaml"]
        sys.argv=["", "-c","config.yaml"]
    
    main()
