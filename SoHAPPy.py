"""
Create a source list to be simulated and analysed using the parameters given
in a configuration file, the command line, and files on disks.

Notes for experts:

* Change :code:`warnings.filterwarnings('ignore')` into
  :code:`warnings.filterwarnings('error')` to have the code stopped in case of
  warning, and be able to identify its origin.

* The 'ignore' option is motivated by warnings issued by `astropy` for
  deprecation or too distant dates with respect to the running date.

* IERS data are not refreshed as it can take long in case of bad Internetor no
  connections at all. To refresh these data use:

    ..  code-block:: python

        # For refreshing
        print(" Refreshing IERS")
        from astroplan import download_IERS_A
        download_IERS_A
        print(" ->Done")

"""

__all__ = ["main"]

import sys
import os
import warnings
import tarfile
import logging

import time
from   datetime import datetime
from   pathlib  import Path

from matplotlib.backends.backend_pdf import PdfPages

from astropy.utils import iers

import gammapy
from __init__ import __version__

from configuration  import Configuration
from niceprint      import Log, failure, heading

from grb import GammaRayBurst
from timeslot import Slot
from mcsim  import MonteCarlo
from analyze import Analysis

# Do not refresh IERS data
iers.conf.auto_download = False

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
    The `SoHAPPy` main function.

    1. Manage input/output
        - Source data identifier list
        - Open output simulation and log files
        - load configuration parameters

    2. Loop over input object list
        - Get source data from the identifiers
        - Compute the visibility on all sites
        - Get the delays on all sites
        - Create original time slot from the source data
        - Loop over site configurations (N, S abd Both)
            - Create a MonteCarlo object
            - Modify the time slot fopr the visibility including the delays
            - If GRB is still vsisible:
                - Dress the GRB slot with physics (IRF and spectra).
                - Create and initialize and Analysis object.
                - Run the simulation.
                - Analyse the simulation.

    3. Close files, terminate

    Returns
    -------
    None.

    """

    # Change gammapy logging to avoid warning messages
    logging.basicConfig()
    log = logging.getLogger("gammapy.irf")
    log.setLevel(logging.ERROR)

    ### ------------------------------------------------
    ### Configuration and output files
    ### ------------------------------------------------

    # Retrieve the input and output base folder from environment variables
    if "HAPPY_IN"  in os.environ.keys():
        infolder = Path(os.environ["HAPPY_IN"])
    else:
        sys.exit("The HAPPY_IN environment variable should be defined")

    if "HAPPY_OUT" in os.environ.keys():
        resfolder = Path(os.environ["HAPPY_OUT"])
    else:
        sys.exit("The HAPPY_OUT environment variable should be defined")

    # Build the Configuration, from the defaults, a configuration file and
    # the command line arguments (sys.argv) if any.
    cf = Configuration.command_line()

    data_path = Path(infolder,cf.data_dir) # Input data folder

    if cf.prompt_dir is not None:
        cf.prompt_dir = Path(infolder, cf.prompt_dir)

    # Create output folder
    # The subfolder name Follows the convention:
    # "population name"/"user keyword"/"visibility keyword and identifiers"
    res_dir   = cf.create_output_folder(resfolder)

    # IRF folder
    irf_dir = Path(infolder,cf.irf_dir)

    # This is required to have the EBL models read from gammapy
    os.environ['GAMMAPY_DATA'] = str(Path(Path(__file__).absolute().parent,
                                          "data"))

    # Backup the current configuration data for further use
    # This is potentially a modified version of the local file
    conf_filename = Path(res_dir, cf.filename.name)
    cf.write(out_name = conf_filename)

    # Output file names
    sim_filename = Path(res_dir, cf.datafile) # Population file (data.txt)
    log_filename = Path(res_dir, cf.logfile)  # Log file

    # Open log file - If Silent is True, only in file, otherwise on Screen too
    log = Log(log_name = log_filename, talk = not cf.silent)

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
        MonteCarlo.welcome(cf.arrays, log = log) # Remind simulation parameters

        first = True # Actions for first GRB only

        for item in cf.srclist:

            # If silence required, keep at least the event number for crashes
            if cf.silent:
                if first is True:
                    print("Processing items :", end=" ")
                print(item,end=' ')

            ### Get GRB
            if isinstance(item, int): # from a number as an indentifier
                fname = Path(data_path,cf.prefix+str(item)+cf.suffix)

                if not cf.test_prompt: # Afterglow + time integrated prompt
                    if not fname.is_file():
                        failure(f" SKIPPING - File not found {fname:}")
                        continue
                    grb = GammaRayBurst.from_fits(fname,
                                              prompt  = cf.prompt_dir,
                                              ebl     = cf.ebl_model,
                                              emax    = cf.emax,
                                              dt      = cf.tshift,
                                              magnify = cf.magnify,
                                              tmax    = cf.tmax)
                else: # Prompt component alone

                    pname = Path(Path(cf.infolder, cf.prompt_dir,
                                      "events_"+str(item)+".fits"))
                    if not cf.use_afterglow:
                        fname = None

                    grb = GammaRayBurst.prompt(pname, fname,
                                               ebl     = cf.ebl_model,
                                               magnify = cf.magnify,
                                               emax    = cf.emax,
                                               tmax    = cf.tmax)

            elif isinstance(item, str): # this is a GRB name string
                if cf.visibility == "built-in":
                    sys.exit(" Error: yaml GRB file with `built-in` visibility")
                grb = GammaRayBurst.historical_from_yaml(item,
                                                         ebl = cf.ebl_model)

            # Assign visibilities
            for loc in ["North","South"]:
                grb.set_visibility(item, loc,
                                   info    = visinfo,
                                   n_night = cf.maxnight,
                                   n_skip  = cf.skip)

            # Printout grb, visibility windows, display plots

            if cf.save_fig and cf.show > 0:
                pdf_out = PdfPages(Path(grb.id+"_booklet.pdf"))
            else:
                pdf_out = None

            if (cf.niter<=1 and cf.do_fluctuate is True) \
                or cf.dbg>0 or cf.nsrc==1 :
                heading(grb.id)
                log.prt(grb)
                grb.vis["North"].print(log=log)
                grb.vis["South"].print(log=log)

            if cf.show > 0 :
                grb.plot(pdf = pdf_out)
            if cf.save_grb :
                grb.write_to_bin(res_dir)

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
                    log.banner(f" SIMULATION  : {name:<50s} ")

                # Create a MC object
                # It has dummy values that will be dumped to the output
                # even if the simulation is not possible (not visible)
                mc = MonteCarlo(niter     = cf.niter,
                                fluctuate = cf.do_fluctuate,
                                nosignal  = (cf.magnify==0),
                                seed      = cf.seed,
                                name      = name,
                                dbg       = cf.dbg)

                still_vis = False # Assumed not visible

                ### ------------
                ### Both sites - create a slot
                ### ------------
                if loc == "Both":
                    if  grb.vis["North"].vis_night \
                    and grb.vis["South"].vis_night:
                        slot = origin.both_sites(delay = delay,
                                                 debug = (cf.dbg>1))
                        if slot is not None:
                            still_vis = True

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
                    # Add IRF features and run - Note that this can
                    # modify the number of slices (merging)
                    slot.dress(irf_dir = irf_dir,
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
                if cf.save_simu:
                    mc.write(Path(res_dir,name + "_sim.bin"))
                    ana.write(Path(res_dir,name+"_ana.bin"))

                # Display status - even if simulation failed (not visible)
                if cf.dbg :
                    mc.status(log=log)

                ### ------------
                ### Analyze simulated data
                ### ------------
                if ana.err == mc.niter:  # Simulation is a success
                    ana.run()
                    if cf.dbg  :
                        ana.print(log = log)
                    if cf.show :
                        ana.show(pdf = pdf_out)

                # Even if not detected nor visibile, dump to file
                first = ana.dump_to_file(grb, pop, header=first)

            if pdf_out is not None:
                pdf_out.close()

            # End of loop over sites
        # END of Loop over GRB
        if cf.silent:
            print("") # Line break

    # Stop chronometer
    end_pop = time.time()
    elapsed = end_pop - start_pop

    log.prt("\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    log.prt(f" Duration   = {elapsed:8.2f} (s)")
    log.prt(f"  per GRB   = {elapsed/cf.nsrc:8.2f} (s)")
    log.prt(f"  per trial = {elapsed/cf.nsrc/cf.niter:8.3f} (s)")
    log.prt("-*-*-*-*-*-*-*-*- End of population simulation -*-*-*-*-*-*-*-*-*\n")
    log.prt(f" ******* End of job - Total time = {(end_pop-start_pop)/60:8.2f} min *****")

    log.prt("")
    log.prt(datetime.now())

    # Close log file
    log.close()

    # tar gzip outputs, delete originals if requested
    nw = datetime.now()
    outprefix = Path(res_dir).parts[-1]
    filename  = outprefix + "_" + nw.strftime("%Y%m%d_%H%M%S") \
                                +".tar.gz"

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
        if not log.log_file.closed:
            log.log_file.close()

        try:
            os.remove(log_filename)
        except IOError:
            print(f"{log_filename} removal failed: locked")

    tar.close()
    print("... completed")

###############################################################################
if __name__ == "__main__":

    os.environ["HAPPY_IN"] = "D:\\CTA\SoHAPPy\input"
    os.environ["HAPPY_OUT"] = "D:\\CTA\SoHAPPy\output"

    print(" argv : ",sys.argv," len = ",len(sys.argv))
    if len(sys.argv[1:]) <= 1:
        print("------------------> Execute examples")
        # sys.argv=["", "-c","myConfigs/config-LongFinalTest-omega.yaml"]
        # sys.argv=["", "-c","myConfigs/config_Long1000_strictmoonveto_1.yaml"]
        # sys.argv= ["", "-c","data/config_ref.yaml", ]
        sys.argv= ["", "-c","config_maximal_detection14h.yaml", ]
        # sys.argv= ["", "--first", "1", "--nsrc", "3",
        #            "--visibility", "strictmoonveto_9999_1_interactive_test",
        #            "-d", "0"]
                   # "--config", r"//dapdc5/Stolar/My Documents/CTA_Analysis/GRB paper/SoHAPPy/data/config_ref.yaml",
    main()
