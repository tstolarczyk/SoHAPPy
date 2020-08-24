"""

Input parameters:
-----------------

    The parameters are described in ana_congif.py

    Folders (by default) are :
        - the IRF are stored in :math:`irf`
        - the GRB files are in :math:`data`
        - the documentation is generated with sphynx in :math:`docs`
        - some plots are in :math:`out`
        - results are where you wish, typicall :math:`Result`


Description:
------------

- The parameters are read from the ana_congfi.py file which is renamed at each execution with the date (_hhmmss)

(-) Loop over the reconstruction chain (reco1, reco2) and the IRF observation duration

    - Opens an output file for that set

    (-) Loop over the site position, zenith and azimuth (this will be replaced by a GRB position later)
        - Compute the pointing (temporary)
        - Read the propoer IRF, stored in :math:`cta\_perf` (depending on the analysis chosen)
        - Read the absorption from the EBL model

        (-) Loop over the GRB list and the reduction factor
            - Create a :class:`grb.GammaRayBurst` object (includes the absorbed spectrum and the positionin the sky)
            - Prepare an output file for the given GRB
            - Create a :class:`montecarlo.GRBMonteCarlo` object from the GRB files
                - get the charcateristics (name, redshift...)

                (-) Loop over time slices
                    - Read differential flux as an array :math:`(E, df/dE)`

            (-) Loop over trials (:meth:`montecarlo.GRBMonteCarlo.loop_3D` or :meth:`montecarlo.GRBMonteCarlo.loop_OnOff`)

                - Initialize the list of results from each trials
                -  Depending on the analysis chain :
                    (-) loop over time slices :meth:`grb.GammaRayBurst.run_OnOffsimulation`
                        - Create a target :class:`grb.cta_grb_observation.GRBObservationParameters.GRBTarget`
                        - Create an observation :class:`cta_grb_observation.GRBObservationParameters`
                        - Simulate the current slice :meth:`cta_grb_observation.GRBObservationSimulation.simulate_OnOffobs`


                    (-) loop over time slices :meth:`grb.GammaRayBurst.run_3Dsimulation`
                        - Create a target :class:`grb.cta_grb_observation.GRBObservationParameters.GRBTarget`
                        - Create an observation :class:`cta_grb_observation.GRBObservationParameters`
                        - Simulate the current slice :meth:`cta_grb_observation.GRBObservationSimulation.simulate_3Dobs`

                -   Analyze results, get :math:`3 sigma` and :math:`5 sigma` responses

"""
###############################################################################
import os
import sys, getopt
import shutil
import time
import datetime
import pickle

from   pathlib import Path
from   astropy.table       import Table

# Gammapy
import gammapy
from   gammapy.spectrum.models import Absorption

from grb            import GammaRayBurst

from   mcsim        import MonteCarlo
import mcsim_res    as mcres

from utilities import warning, failure, success, highlight,banner

import ana_config as cf # Steering parameters
# from . import __version__ # does not work
__version__ = "Roma"
 
os.environ['GAMMAPY_EXTRA'] =r'../input/gammapy-extra-master'
os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
#print(os.getenv('GAMMAPY_EXTRA'))
#print(os.listdir(os.getenv('GAMMAPY_EXTRA')))

###############################################################################
def save_config(dbg=False):
    """
    Copy configuration file to result folder
    """

    if (not os.path.exists(cf.res_dir)):
            if (dbg): print(" *** Creating result folder ",cf.res_dir)
            os.makedirs(cf.res_dir)
    config_file = cf.res_dir+"/config.txt"

    if (os.path.exists(config_file)):
        now = datetime.datetime.now()
        newname = config_file+"_"+str(now.hour)+str(now.minute)+str(now.second)
        os.rename(config_file, newname)
        if (dbg): print("     --- config.txt exists, renamed to ",newname)

    shutil.copy('ana_config.py',cf.res_dir+'/config.txt')

    return

###############################################################################
def get_grb_file(i, ebl = None):
    """
    Obtain ith GRB file

    Parameters
    ----------
    i : TYPE
        DESCRIPTION.

    Returns
    -------
    name : TYPE
        DESCRIPTION.

    """
    if (cf.old_file):
        loc = Path(cf.grb_olddir + "LGRB"+str(i))
        grb = GammaRayBurst.read_old(loc,ebl = ebl)
    else:
        loc = Path(cf.grb_dir + "/Event"+str(i)+".fits")
        grb = GammaRayBurst.read(loc, ebl = ebl)
        
    if (cf.save_grb):
        grb_class_file = cf.res_dir + "/" \
                       + grb.name + ".bin"
        outfile  = open(grb_class_file,"wb")
        pickle.dump(grb,outfile)
        outfile.close()
        print(" Saving grb {} to file : {}"
              .format(grb.name,grb_class_file))

    if cf.niter<=1 or cf.dbg>0 or cf.ngrb==1 :
        print(grb)
        visible = grb.visibility()
        if (cf.show >0) and (visible):
            import grb_plot     as gplt
            gplt.spectra(grb,opt="Packed")
            gplt.visibility_chart(grb)
            gplt.visibility_alt_plot(grb)
            gplt.pause()

            if (cf.dbg>2) :
                gplt.animated_spectra(grb,savefig=True,outdir=cf.res_dir)
#            plt.show(block=True) 
            
    return grb

###############################################################################
def init(cf,argv):
    
    """
    Debugging is controlled as follows with the dbg variables
    dbg = 0 No debug, no plot
          1 First level of debug, summary plots
          2 Second level of debug, event plots
          3 Thors level of debug, at the trial level
    dbg < 0 No plots, absolute value as above      
    
    """
    
    # Create the show debugging flag from the general debug flag
    if (cf.dbg < 0):
        cf.show = 0
    else:
        cf.show = abs(cf.dbg)
            
    # Read arguments from command line, supersede default ones
    try:
          opts, args = getopt.getopt(argv,
                                     "hn:f:d:N:",
                                     ["ngrb=",
                                      "first=",
                                      "niter=",
                                      "dbg="])    
    except getopt.GetoptError:
        print("No line arguments passed: Using default values")
    
    for opt, arg in opts:
      if opt == '-h':
         print(" SoHAPPy.py -N <ngrb> -f <1st grb> -n <iterations> -d <debug>")
         sys.exit()
      elif opt in ("-N", "--ngrb"):
         cf.ngrb =  int(arg)
      elif opt in ("-f", "--first"):
         cf.ifirst = int(arg)
      elif opt in ("-n", "--niter"):
         cf.niter = int(arg)
      elif opt in ("-d", "--debg"):
         dbg = int(arg)
         cf.dbg = abs(dbg)
         if (dbg<0):
             cf.show = 0
         else:
            cf.show = abs(dbg)
        
    if (cf.do_fluctuate == False): cf.niter=1
        
    # Summarise the situation
    print("")
    print("+----------------------------------------------------------------+")
    print("|                                                                |")
    print("|                    SoHAPPy with GammaPy {:4s}                   |"
          .format(gammapy.__version__))
    print("|                            ({:4s})                              |"
          .format(__version__))
    print("|  (Simulation of High-energy Astrophysics Processes in Python)  |")
    print("|                                                                |")
    print("+----------------------------------------------------------------+")
    print(" Simulation:")
    print("     *Number of sources  : ",cf.ngrb)
    print("     *First source       : ",cf.ifirst)
    print("     *Number of trials   : ",cf.niter)
    print(" EBL model               : ",cf.EBLmodel)
    print(" Input/output :")
    print("     *Debug mode         : ",cf.dbg)
    print("     *Show plots         : ",cf.show)
    print("      Analysing files in : ",cf.grb_dir)
    print("      IRF files in       : ",cf.irf_dir)
    print("      Result folder      : ",cf.res_dir)
    if (cf.lkhd ==0):
        method = "On-Off"
    else:
        method = "{}D likelihood".format(cf.lkhd)
    print(" Analysis (ndof)         : ",method)
    print("+----------------------------------------------------------------+")
    print("|                 *: can be changed with command line (use -h)   |")
    print("+----------------------------------------------------------------+")
    print(" Developments:")
    if (cf.save_grb == True):
        highlight("Simulation saved to disk save_grb (save_grb = True)       ")        
    if (cf.write_slices == True):
        highlight("Slice information saved to disk (write_slices=True)       ")
    if (cf.signal_to_zero == True):
        warning(  "Signal set to zero (signal_to_zero==True)                 ")
    if (cf.do_fluctuate == False):
        warning(  "No fluctuation in simulation (do_fluctuate==False)        ")
    if (cf.do_accelerate  == False):
        warning(  "No abortion if first 10% undetected (do_accelarate==False)")
    if (cf.fixed_zenith != False):
        warning(  "Zenith angle requested to be fixed at keyword '{:5s}'     "
               .format(cf.fixed_zenith))
    if (cf.day_after == True):
        warning(  "Visibility windows the day after                          ")

    print()

    return

###############################################################################
def main(argv):

    init(cf,argv)           # Sterring parametets and welcome message
    if (cf.ngrb<=0):
        print(" NO ANALYSIS REQUIRED (ngrb<=0)")
        sys.exit(2)
        
    save_config()           # save the config file in the result folder
    start_all = time.time() # Starts chronometer
        
    # Create population summary file
    sumfile = cf.res_dir+"/Pop_" + str(cf.ifirst)+"-" \
                                    + str(cf.ifirst+cf.ngrb-1)+"GRB_" \
                                    + str(cf.lkhd)+"dof_"+ str(cf.niter)+"iter"
    print("======>  Population file :\:n",sumfile,"\n")

    start_pop = time.time() # Start chronometer

    with open(sumfile+".txt", 'w') as popfile:

        ##############################
        # Loop over GRB population   #
        ##############################
        
        eblabs = Absorption.read_builtin(cf.EBLmodel) # Initialise absorption
    
        first_simul = True
        for i in range(cf.ifirst,cf.ifirst+cf.ngrb):

            grb = get_grb_file(i,ebl = eblabs) ### Get GRB

            for location in ["North","South"]:

                # Create Monte Carlo default object
                mc = MonteCarlo(grb,niter  = cf.niter,
                                    where  = location,
                                    ndof   = cf.lkhd,
                                    debug  = cf.dbg)
                if first_simul : mcres.welcome()

                mc.run()
                mcres.result(mc,
                            popfile      = popfile,
                            write_header = first_simul)
                first_simul = False
                    
                # If requested save simulation to disk
                if (cf.save_simu):
                    mc_class_file = cf.res_dir + "/" \
                                  + grb.name + "-" \
                                  + location + "-" \
                                  + str(cf.niter)+".bin"
                    outfile  = open(mc_class_file,"wb")
                    pickle.dump(mc,outfile)
                    outfile.close()
                    print(" Saving simulation of {} to file : {}"
                          .format(mc.grb.name,mc_class_file))
                    
                # If GRB visible
                if (mc.err > -1) and (abs(cf.dbg)>1) : mc.compare_time_interval()
                
                # If Simulation was not aborted - time to look at results
                if (mc.err == mc.niter) and (cf.show > 0):
                    import mcsim_plot   as mplt
                    mplt.stat(mc,cf.show,
                              Path(cf.res_dir,grb.name))
                    mplt.story(mc,cf.show,
                               Path(cf.res_dir,grb.name),ref="VIS")
                    mplt.story(mc,cf.show,
                               Path(cf.res_dir,grb.name),ref="GRB")
                    mplt.pause() 
        # END of Loop over GRB

    end_pop = time.time()
    elapsed = end_pop-start_pop  # Stop chronometer

    # Convert population file into smarter format
    # Note : Guess=False is very useful to understand potential crashes
    data = Table.read(sumfile+".txt",format="ascii",guess=False)
    data.write(sumfile+'.csv', format="ascii.csv",overwrite="True")
    data.write(sumfile+'.fits',format="fits",     overwrite="True")

    print("\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    print(" Simulation duration   = {0:8.2f} (s)".format(elapsed))
    print("             per GRB   = {0:8.2f} (s)"
          .format( (elapsed)/cf.ngrb))
    print("             per trial = {0:8.3f} (s)"
          .format( (elapsed)/cf.ngrb/cf.niter))
    print("-*-*-*-*-*-*-*- End of full population simulation -*-*-*-*-*-*-*-*\n")

    end_all = time.time()
    print(" ******* End of job - Total execution time = ",
          round((end_all-start_all)/60,1),"  minutes *****")
                    
    # This prevents having losing the pictures at end of scripts


###############################################################################
if __name__ == "__main__":
    main(sys.argv[1:])