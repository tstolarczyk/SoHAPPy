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

Todo :
------
    - Add number of counts to the output;
    - Compute the fluence;
    - Interpolate time counts within the time slices;
    - in output file, separate constant information from variable information (use a header).

"""
###############################################################################
import os

import shutil
import time
import datetime
import matplotlib.pyplot   as plt
from   pathlib import Path

from   astropy.table       import Table

# Gammapy
import gammapy
from   gammapy.spectrum.models import Absorption

from grb            import GammaRayBurst
import grbplot      as gplt
from montecarlo     import MonteCarlo

# Steering parameters
import ana_config as cf

# Styles
#plt.rcdefaults()
#plt.rcParams['font.weight'] = 'bold'
#plt.rcParams['axes.labelweight'] = 'bold'
#plt.rcParams['axes.titleweight'] = 'bold'
#plt.style.use('seaborn-talk') # Make the labels readable
plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!
# print(plt.style.available)

###############################################################################
def save_config(dbg=False):
    """
    Copy configuration file to result folder
    """

    if (not os.path.exists(cf.res_folder)):
            if (dbg): print(" *** Creating result folder ",cf.res_folder)
            os.makedirs(cf.res_folder)
    config_file = cf.res_folder+"/config.txt"

    if (os.path.exists(config_file)):
        now = datetime.datetime.now()
        newname = config_file+"_"+str(now.hour)+str(now.minute)+str(now.second)
        os.rename(config_file, newname)
        if (dbg): print("     --- config.txt exists, renamed to ",newname)

    shutil.copy('ana_config.py',cf.res_folder+'/config.txt')

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
        loc = Path(cf.grb_oldfolder + "LGRB"+str(i))
        grb = GammaRayBurst.read_old(loc,ebl = ebl)
    else:
        loc = Path(cf.grb_folder + "/Event"+str(i)+".fits")
        grb = GammaRayBurst.read(loc, ebl = ebl)

    if cf.niter<=1 or cf.dbg_level>0 or cf.ngrb==1 :

        gplt.visibility(grb)
        gplt.spectra(grb,opt="Packed")
        print(grb)
    return grb

###############################################################################
def welcome(cf):
    print("")
    print("#################################################################################")
    print("################ ")
    print("################ SoHAPPy with GammaPy ",gammapy.__version__)
    print("################ (Simulation of High-energy Astrophysics Processes in Python)")
    print("################ ")
    print("#################################################################################")
    print("  Analysing files in     : ",cf.grb_folder)
    print("  Number of sources      : ",cf.ngrb)
    print("  First source index     : ",cf.ifirst)
    print("  IRF files in           : ",cf.irf_folder)
    print("  Reduction factor       : ",cf.redfactor)
    print("  Alpha                  : ",cf.alpha)
    print("  FOV                    : ",cf.fov[0])
    print("  Bin size               : ",cf.fov[1])
    print("  Slewing time           : ",cf.dtslew[0])
    print("                   Fixed : ",cf.dtslew[1])
    print("  Detection level        : ",cf.det_level)
    print("  EBL model              : ",cf.EBLmodel)
    print("")
    print("  Number of MC trials    : ",cf.niter)
    print("  Likelihood             : ",cf.lkhd)
    print("  Debug mode             : ",cf.dbg_level)
    print("  Show plots             : ",cf.showplots)
    print("  Result folder          : ",cf.res_folder)
    print("")
    return
###############################################################################
if __name__ == "__main__":

        SOLVE THE PROBLEM THAT IT CANNOT RUN IF NO INTERNET
        BECAUSE OF THE SITE DATABASE
    welcome(cf)  # Welcome message, with sterring parameters
    save_config() # save the config file in the result folder
    start_all = time.time() # Starts chronometer

    # Create population summary file
    sumfile = cf.res_folder+"/Pop_" + str(cf.ifirst)+"-" \
                                    + str(cf.ifirst+cf.ngrb-1)+"GRB_" \
                                    + str(cf.lkhd)+"dof_"+ str(cf.niter)+"iter"
    print("======>  Population file : ",sumfile,"\n")

    # Intialise timing and stats
    ndet_3s = 0
    ndet_5s = 0
    start_pop = time.time() # Start chronometer

    # Initialise absorption
    eblabs = Absorption.read_builtin(cf.EBLmodel)

    with open(sumfile+".txt", 'w') as popfile:

        ##############################
        # Loop over GRB population   #
        ##############################
        write_header = True
        for i in range(cf.ifirst,cf.ifirst+cf.ngrb):

            ### Get the source
            grb = get_grb_file(i,ebl = eblabs) ### Get GRB

            #for location in ["North","South"]:
            for location in ["North","South"]:

                print("\n*** SIMULATION - ", grb.name,
                      " - ",location," : ",end="\n")
                # Create Monte Carlo default object
                mc = MonteCarlo(grb,niter  = cf.niter,
                                    alpha  = cf.alpha,
                                    fov    = cf.fov,
                                    dtslew = cf.dtslew,
                                    where  = location,
                                    fixed  = cf.fixed,
                                    ndof   = cf.lkhd,
                                    clevel = cf.det_level)

                mc.run(cf.dbg_level)
                mc.result(debug       = cf.dbg_level,
                          popfile     = popfile,
                          write_header= write_header)
                write_header = False

                # If debugging level high or # evt low, be talkative
#                    if (cf.niter <2) or (cf.dbg_level > 1):
                    #mc.plot_trial()
                if (mc.abort==0) : # Simulation was not aborted
                    if (cf.showplots > 0) or (cf.dbg_level > 2):
                        # grb.quicklook(mc.simulations,plot=True ) # deprecated ?
                        mc.plot(cf.showplots,Path(cf.res_folder,grb.name))
                        mc.plot2(cf.showplots,Path(cf.res_folder,grb.name),ref="VIS")
                        mc.plot2(cf.showplots,Path(cf.res_folder,grb.name),ref="GRB")

        # END of Loop over GRB       #

    end_pop = time.time()
    elapsed = end_pop-start_pop  # Stop chronometer

    # Convert population file into smarter format
    # Note : Guess=False is very useful to understand potential crashes
    data = Table.read(sumfile+".txt",format="ascii",guess=False)
    data.write(sumfile+'.csv', format="ascii.csv",overwrite="True")
    data.write(sumfile+'.fits',format="fits",     overwrite="True")

    print("\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n")
    print(" Simulation duration   = {0:8.2f} (s)".format(elapsed))
    print("             per GRB   = {0:8.2f} (s)"
          .format( (elapsed)/cf.ngrb/len(cf.redfactor)))
    print("             per trial = {0:8.3f} (s)"
          .format( (elapsed)/cf.ngrb/len(cf.redfactor)/cf.niter))
    print("-*-*-*-*-*-*-*- End of full population simulation -*-*-*-*-*-*-*-\n")


    end_all = time.time()
    print(" ****************** End of job - Total execution time = ",
          round((end_all-start_all)/60,1),"  minutes *****")



    #     if (cf.niter<=1):

    #         if (cf.likelihood == 0):
    #             if (cf.dbg_level):
    #                 gplt.stats_detection(grb,
    #                                      mc.simulations,
    #                                      savefig=True,
    #                                      outdir='./out/')
    #             # PlotGRB.make_gif_from_models(grb, savefig=True, outdir='./out/')
    #         # plt.show(block=False)

