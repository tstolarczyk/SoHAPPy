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
import itertools
import datetime
import matplotlib.pyplot as plt

from astropy.table import Table
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

# Personnal utilities
import utilities as ut

# Gammapy
import gammapy
from gammapy.spectrum.models import Absorption
from gammapy.irf import load_cta_irfs # analysis3D - new

from cta_irf_onaxis import CTAPerf_onaxis # ...whereas this does work
from grb import GammaRayBurst
from grb_plot import PlotGRB
from montecarlo import MonteCarlo

# Steering parameters
import ana_config as cfg

# Styles
#plt.rcdefaults()
#plt.rcParams['font.weight'] = 'bold'
#plt.rcParams['axes.labelweight'] = 'bold'
#plt.rcParams['axes.titleweight'] = 'bold'
plt.style.use('seaborn-talk') # Make the labels readable
#plt.style.use('seaborn-poster') # Make the labels readable
# print(plt.style.available)

###############################################################################
def save_config(dbg=False):
    """
    Copy configuratio file to result folder
    """

    if (not os.path.exists(cfg.res_folder)):
            if (dbg): print(" *** Creating result folder ",cfg.res_folder)
            os.makedirs(cfg.res_folder)
    config_file = cfg.res_folder+"/config.txt"

    if (os.path.exists(config_file)):
        now = datetime.datetime.now()
        newname = config_file+"_"+str(now.hour)+str(now.minute)+str(now.second)
        os.rename(config_file, newname)
        if (dbg): print("     --- config.txt exists, renamed to ",newname)

    shutil.copy('ana_config.py',cfg.res_folder+'/config.txt')

    return

###############################################################################
def pop_out(out=None, grb=None, mc=None, hem="Unknown", theta="361", phi="361",
           header=False):
    """
    Write header and contents to the poupulation summary file
    """
    if (header):
        print(  "GRB "    + "fluence " + "z "
              + "Site "   + "Zenith "  + "Azimuth "
              + "ndof "
              + "sigmax " + "nex_max " + "nb_max "
              + "nex_3s " + "nb_3s "   + "t3s " + "det3s "
              + "nex_5s " + "nb_5s "   + "t5s " + "det5s ",
              file = out)
    else:
        print(grb.name, grb.fluence, grb.z,
              hem, theta, phi,
              mc.ndof,
              mc.sigmax, mc.nex_sigmax, mc.nb_sigmax,
              mc.nex_3sigma, mc.nb_3sigma, mc.t_3sigma,mc.detect_3sigma,
              mc.nex_5sigma, mc.nb_5sigma, mc.t_5sigma,mc.detect_5sigma,
              file = out)
    return
###############################################################################
def fig_out():
    return
###############################################################################
if __name__ == "__main__":

    print("")
    print("#################################################################################")
    print("################ ")
    print("################ SoHAPPy with GammaPy ",gammapy.__version__)
    print("################ (Simulation of High-energy Astrophysics Process in Python)")
    print("################ ")
    print("#################################################################################")
    print("  Analysing files in     : ",cfg.grb_folder)
    print("  IRF files in           : ",cfg.irf_folder)
    print("  Reduction factor       : ",cfg.redfactor)
    print("  Alpha                  : ",cfg.alpha)
    print("  FOV                    : ",cfg.fov)
    print("  Bin size               : ",cfg.binsize)
    print("  Detection level        : ",cfg.det_level)
    print("  EBL model              : ",cfg.EBLmodel)
    print("")
    print("  Number of MC trials    : ",cfg.niter)
    print("  Likelihood             : ",cfg.likelihood)
    print("  Debug mode             : ",cfg.dbg_level)
    print("  Show plots             : ",cfg.showplots)
    print("  Result folder          : ",cfg.res_folder)

    save_config() # save the config file in the result folder

    #------------------------------------------------------------------------
    ### Summarise the result of a population study in a fits file
    ### or series of Monte Carlo simulations with various IRF
    #------------------------------------------------------------------------

    start_all = time.time()

    for ch,dt in itertools.product(cfg.chain,cfg.duration):

        sumfile =  cfg.res_folder + "/PopSummary_" \
                                 + ch + "_" \
                                 + dt + "_" \
                                 + str(cfg.likelihood) + "dof" \
                                 + str(cfg.niter) + "iter"
        if (cfg.NSB): sumfile += "_NSB"
        with open(sumfile+".txt", 'w') as popfile:

            pop_out(popfile,header=True) # Write population file header

            #------------------------------------------------------------------------------
            ### Run a series of Monte Carlo simulations with various IRF
            #------------------------------------------------------------------------------
            for hem,theta,phi in itertools.product(cfg.hemisphere,cfg.zenith,cfg.azimuth):

                # Create poiting to the source
                print("                       Warning phi=",phi," -> phi=0")
                pointing =  SkyCoord(alt = 90-float(theta),
                                     az  = 0,unit = u.deg,
                                     obstime = Time([cfg.obstime], scale='utc'),
                                     frame = 'altaz',
                                     location = EarthLocation.of_site(cfg.site) )

                if (cfg.likelihood != 0):
                    # status     =
                    # irf_folder =
                    irf_file   = "irf_file-North_z20_average_100s.fits"
                    irf_file_full = cfg.irf_folder + "/" + irf_file
                    cta_perf = load_cta_irfs(irf_file_full) # Load IRF
                    reco_energy = cta_perf["bkg"].data.axis("energy").edges
                else:
                    status,folder,irf_file = ut.GetIrfFile(cfg.irf_folder,
                                                         ch, hem, theta, phi,
                                                         dt, cfg.NSB)
                    irf_file_full = folder+"/"+irf_file
                    cta_perf = CTAPerf_onaxis.read(irf_file_full)
                    # cta_perf.peek()
                    reco_energy = cta_perf.bkg.energy.edges

                Emin = min(reco_energy)
                Emax = max(reco_energy)
                absorption = Absorption.read_builtin(cfg.EBLmodel)

                print("========================================================")
                print("==== IRF : ",irf_file)
                print("========================================================")
                print("  Analysis chain : ",ch)
                if (cfg.likelihood !=0):
                    print("  Site           :  ",cfg.site)
                    print("  Date           :  ",cfg.obstime)
                print("  Hemisphere     : ",hem)
                print("  Zenith         : ",theta)
                print("  Azimuth        : ",phi)
                print("  Duration       : ",dt)
                print("  High NSB       : ",cfg.NSB)
                print("  Emin (IRF)     : {0:6.2f}".format(Emin.to(u.GeV)))
                print("  Emax (IRF)     : {0:6.2f}".format(Emax.to(u.TeV)))

                ndetected_3sigma = 0
                ndetected_5sigma = 0
                start_pop = time.time()

                for source,rd in itertools.product(cfg.grblist,cfg.redfactor):

                    grb_file_full   = cfg.grb_folder+"/"+source

                    ### Create the GRB time-interval spectra from the given data file
                    grb = GammaRayBurst.read(filepath   = grb_file_full,
                                             absorption = absorption,
                                             reduction_factor = rd,
                                             pointing = pointing,
                                             dbg = cfg.dbg_level)

                    #--- Observation investigation
                    if (cfg.dbg_level > 2) : grb.quicklook(plot=True )

                    # Modifiy name for pseudo-population
                    if (rd != 1): grb.name = grb.name+"-rd"+str(round(rd,1))

                    # Prepare outputfile for the present GRB and configuration
                    out_file_noext = ut.OutputPrefixName(grb.name,
                                                         ch, hem, theta, phi, dt,
                                                         cfg.niter, cfg.NSB)
                    result_file = open(cfg.res_folder+'/'+out_file_noext+".txt", 'w')

                    # Run simulations of the current grb
                    # Compute significance on the full observation,
                    # significance max on intervals,
                    # Time of siginificance max, time of 3 sigma reached

                    ### Run the Monte Carlo
                    print("\n*** SIMULATION - ", grb.name," : ",end="")
                    mc = MonteCarlo(cfg.niter,
                                    cfg.fov,
                                    cfg.binsize,
                                    cfg.alpha,
                                    Emin,
                                    Emax,
                                    cta_perf,
                                    grb,
                                    ndof = cfg.likelihood)

                    mc.run(cfg.likelihood,cfg.dbg_level) # Run the simulation

                    mc.result(out=result_file) # in file
                    mc.plot(cfg.showplots,cfg.res_folder+'/'+out_file_noext)

                    if (cfg.niter<=1):
                        print(grb)
                        if (cfg.likelihood == 0):
                            if (cfg.dbg_level):
                                PlotGRB.plot_stats_detection(grb,
                                                               mc.simulations,
                                                               savefig=True,
                                                               outdir='./out/')
                            # PlotGRB.make_gif_from_models(grb, savefig=True, outdir='./out/')
                        # plt.show(block=False)

                    result_file.close()

                    ### Collect some statistics
                    print()
                    if (mc.detect_3sigma >= cfg.det_level):
                        print(grb.name," => 3 sigma detected")
                        ndetected_3sigma += 1
                    if (mc.detect_5sigma >= cfg.det_level):
                        print(grb.name," => 5 sigma detected")
                        ndetected_5sigma += 1

                    pop_out(out=popfile, grb=grb, mc=mc)
                # For a source

            # For an observation position

        # --- end of population study
        end_pop = time.time()
        elapsed = end_pop-start_pop
        print(" Simulation duration   = {0:8.2f} (s)".format(elapsed))
        print("             per GRB   = {0:8.2f} (s)"
              .format( (elapsed)/len(cfg.grblist)/len(cfg.redfactor)))
        print("             per trial = {0:8.3f} (s)"
              .format( (elapsed)/len(cfg.grblist)/len(cfg.redfactor)/cfg.niter))
        print("-*-*-*-*-*-*-*- End of full population simulation -*-*-*-*-*-*-*-\n")

        # Convert population file into smarter format
        # Note that this is safer than accumulating data in a list and dumping them in
        # the end : the "with" loop creating thet txt file ensures that the data
        # are written on the fly and accummulated data are not lost if a crash occurs
        # Note : Guess=False is very useful to understand potential crashes
        data = Table.read(sumfile+".txt",format="ascii",guess=False)
        data.write(sumfile+'.csv', format="ascii.csv",overwrite="True")
        data.write(sumfile+'.fits',format="fits",     overwrite="True")


        # End : Open population file for the analysis chain and obs. time

    # End : for ch,dt  - analysis chain and observation time

    end_all = time.time()
    print(" ****************** End of job - Total execution time = ",
          round((end_all-start_all)/60,1),"  minutes *****")
