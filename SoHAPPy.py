###############################################################################
# Main script
#
# INPUT / OUPUT
#  - the parameters are passed through the ana_config.py script
#       - the IRF are stored in irf_repository
#       - the GRB files are stored in grb_repository
#       - the GRB to be analysed are in the grblist variable
#       - the results are stored in a folder res_folder, it contains :
#           - a text file for all GRB individually
#           - summary plots for each GRB on screen and in a jpg file if requested
#           - a population summary file in the following format :
#                - text
#                - fits
#                - csv
# RUNNING TIME
#    - the input configuration file is copied (name with the date _hhmmss)
#------------------------------------------------------------------------------
### Code Initialisations
#------------------------------------------------------------------------------
dbg_level = 0 # 0 : event counting, 1: event counting + some results, 2: + details for each event, 3: + plots 
showplots = 2 # 0 : not shown, not written, 1: shown on screen, 2 : shown on screen and written out

#------------------------------------------------------------------------------
### Simulation Initialisations
#------------------------------------------------------------------------------
niter = 100 # Number of iterations for the Monte Carlo
alpha = 0.01 # 5 zones (1/0.2) to define the off region for background estimation
#------------------------------------------------------------------------------
### Physics Initialisations
#------------------------------------------------------------------------------
det_level       = 0.9 # Will be detected if n=3,5 sigma reached in det_level of the population
#redfactor       = [1, 2, 3, 5, 10, 20, 50, 100] # Initial flux is divided by this value, for test
#redfactor       = [0.5,1, 2] # Initial flux is divided by this value, for test
redfactor = [1]

grb_repository  = 'data/long_grb_test/' # to be changed into a repository to be scanned
#grblist         = ["LGRB1","LGRB2","LGRB3","LGRB4","LGRB5","LGRB6","LGRB7","LGRB8","LGRB9","LGRB10"]
#grblist         = ["LGRB1","LGRB2"]
grblist         = ["LGRB1"]

EBLmodel        = "dominguez"

#------------------------------------------------------------------------------
### Instrument response functions
#------------------------------------------------------------------------------
irf_repository = "irf"
chain       = ["Reco1"] # or Reco2
hemisphere  = ["North"] # or South
#zenith      = ["20","40","60"] # or 40 deg...
zenith      = ["20"] # or 40 deg...
azimuth     = ["average"]  # or N, or S
#duration    = ["100s","30m"] # or 30m, 05h 50h
duration    = ["100s"] # or 30m, 05h 50h
NSB         = False # if True try to find a high NSB file

#
###############################################################################

import os
import shutil
import time
import itertools
import datetime
import matplotlib.pyplot as plt
#plt.rcdefaults()
#plt.rcParams['font.weight'] = 'bold'
#plt.rcParams['axes.labelweight'] = 'bold'
#plt.rcParams['axes.titleweight'] = 'bold'
plt.style.use('seaborn-talk') # Make the labels readable
#plt.style.use('seaborn-poster') # Make the labels readable
# print(plt.style.available)

from astropy.table import Table
# Windows is clever enough to replace / into the native \
# The r in front of the string menas that \ -used in windows- should be interpreted as such
os.environ['GAMMAPY_EXTRA']=r'../gammapy-extra-master'
os.environ['GAMMAPY_DATA'] =r'../gammapy-extra-master/datasets'
#print(os.getenv('GAMMAPY_EXTRA'))
#print(os.listdir(os.getenv('GAMMAPY_EXTRA')))

# Gammapy
import gammapy
from gammapy.spectrum.models import Absorption

# Note : this has been for a while in gammapy 0.6 and then was removed 
from cta_irf_onaxis import CTAPerf_onaxis # ...whereas this does work

# Utilities
from grb_utils import GammaRayBurst
from plot_utils import PlotGammaRayBurst
import utilities as ut
from montecarlo import GRBMonteCarlo

import ana_config as cfg

#
#Todo : 
#    add nb of counts to the output
#    compute the fluence
#    interpolate time counts
#    use full enclosure IRF
#    try likelihood fit
#    remode "deg" in zenith output - done
#    separate constant information from varaible information


#------------------------------------------------------------------------------
### Display 
#------------------------------------------------------------------------------
print("")
print("#################################################################################")
print("################ ")
print("################ SoHAPPy with GammaPy ",gammapy.__version__)
print("################ (Simulation of High-energy Astrophysics Process in Python)")
print("################ ")
print("#################################################################################")
print("  Analysing files in     : ",cfg.grb_repository)
print("  Reduction factor       : ",cfg.redfactor)
print("  Alpha                  : ",cfg.alpha)
print("  Detection level        : ",cfg.det_level)
print("  EBL model              : ",cfg.EBLmodel)
print("")
print("  Number of MC trials    : ",cfg.niter)
print("  Debug mode             : ",cfg.dbg_level)
print("  Show plots             : ",cfg.showplots)
print("  Result folder          : ",cfg.res_folder)

print("...copying config file to result folder")
if (not os.path.exists(cfg.res_folder)):
        print(" *** Creating result folder ",cfg.res_folder)
        os.makedirs(cfg.res_folder)
config_file = cfg.res_folder+"/config.txt"
if (os.path.exists(config_file)):
    now = datetime.datetime.now()
    newname = config_file+"_"+str(now.hour)+str(now.minute)+str(now.second)
    os.rename(config_file, newname)
    print("     --- config.txt exists, renamed to ",newname)
    
shutil.copy('ana_config.py',cfg.res_folder+'/config.txt')        
#------------------------------------------------------------------------------
### Summarise the result of a population study in a fits file
### or series of Monte Carlo simulations with various IRF
#------------------------------------------------------------------------------
start_all = time.time()
for ch,dt in itertools.product(cfg.chain,cfg.duration):
    summary_file_noext = cfg.res_folder+"/PopulationSummary_" \
                                            +ch+"_" \
                                            +dt+"_" \
                                            +str(cfg.niter)+"iter"
    if (cfg.NSB): summary_file_noext = summary_file_noext+"_NSB"
    with open(summary_file_noext+".txt", 'w') as popfile:
        print("GRB fluence z Site Zenith Azimuth sigmax nex_max nb_max nex_3s nb_3s t3s det3s nex_5s nb_5s t5s det5s",file=popfile)
    
        #------------------------------------------------------------------------------
        ### Run a series of Monte Carlo simulations with various IRF
        #------------------------------------------------------------------------------
        for hem,theta,phi in itertools.product(cfg.hemisphere,cfg.zenith,cfg.azimuth):
            #----------------------------------------------------------
            ### Run one Monte Carlo of niter trials
            #----------------------------------------------------------
            status, irf_folder, irf_file = ut.GetIrfFile(cfg.irf_repository, ch, hem, theta, phi, dt, cfg.NSB)
            irf_file_full = irf_folder+"/"+irf_file
            
            cta_perf = CTAPerf_onaxis.read(irf_file_full)
            # cta_perf.peek()
            
            absorption = Absorption.read_builtin(cfg.EBLmodel)
            
            print("=================================================================================")
            print("==== IRF : ",irf_file)
            print("=================================================================================")
            print("  Analysis chain         : ",ch)
            print("  Site                   : ",hem)
            print("  Zenith                 : ",theta)
            print("  Azimuth                : ",phi)
            print("  Duration               : ",dt)
            print("  High NSB               : ",cfg.NSB)
            
            ndetected_3sigma = 0
            ndetected_5sigma = 0
            start_pop = time.time()
            
    #        grb_detected_list = []
            for source,rd in itertools.product(cfg.grblist,cfg.redfactor):
                grb_file_full   = cfg.grb_repository+"/"+source                    
        
                ### Create the GRB time-interval spectra from the given data file
                grb = GammaRayBurst.from_file(filepath=grb_file_full, 
                                              absorption=absorption,
                                              reduction_factor=rd,
                                              dbg=cfg.dbg_level)
                
                if (cfg.dbg_level > 2) : grb.quicklook(plot=True ) #--- Observation investigation
                if (cfg.niter<=1):
                    print(grb)
                    PlotGammaRayBurst.plot_stats_detection(grb, savefig=True, outdir='./out/')
                    #PlotGammaRayBurst.make_gif_from_models(grb, savefig=True, outdir='./out/')
                    plt.show(block=False)
                if (rd != 1): # Modifiy name for pseudo-population
                    grb.name = grb.name+"-rd"+str(round(rd,1))
                    
                    
                # Prepare outputfile for the present GRB and configuration
                out_file_noext = ut.OutputPrefixName(grb.name, 
                                                     ch, hem, theta, phi, dt,
                                                     cfg.niter, cfg.NSB)
                result_file = open(cfg.res_folder+'/'+out_file_noext+".txt", 'w')
        
          
                # Run simulations of the current grb
                # Compute significance on the full observation, 
                # significance max on intervals,
                # Time of siginificance max, time of 3 sigma reached
                ### Create the GRB Monte Carlo object
        
                mc = GRBMonteCarlo(cfg.niter,cfg.alpha,cta_perf,grb)
                
                ### Run the Monte Carlo
                print("\n*** SIMULATION ", grb.name," : ",end="")
                mc.run(cfg.dbg_level)
                
                ### Display results for the current IRF
                # mc.result() # on screen
                mc.result(out=result_file) # in file
                if (cfg.showplots): mc.plot(cfg.showplots,cfg.res_folder+'/'+out_file_noext)
            
                result_file.close()
                
                
                ### Collect some statistics
                print()
                if (mc.detect_3sigma >= cfg.det_level): 
                    print(grb.name," => 3 sigma detected")
                    ndetected_3sigma += 1
                if (mc.detect_5sigma >= cfg.det_level):
                    print(grb.name," => 5 sigma detected")
                    ndetected_5sigma += 1
                
                print(grb.name, grb.fluence, grb.z, hem, theta, phi, 
                      mc.sigmax, mc.nex_sigmax, mc.nb_sigmax, 
                      mc.nex_3sigma, mc.nb_3sigma, mc.t_3sigma,mc.detect_3sigma,
                      mc.nex_5sigma, mc.nb_5sigma, mc.t_5sigma,mc.detect_5sigma,
                      file=popfile)
            # For a source
                
        # For an observation position    

    end_pop = time.time()
    

    
    print(" Simulation duration   = {0:8.2f} (s)".format(end_pop-start_pop))
    print("             per GRB   = {0:8.2f} (s)".format( (end_pop-start_pop)/len(cfg.grblist)/len(cfg.redfactor)))
    print("             per trial = {0:8.3f} (s)".format( (end_pop-start_pop)/len(cfg.grblist)/len(cfg.redfactor)/cfg.niter))
    print("-*-*-*-*-*-*-*- End of full population simulation -*-*-*-*-*-*-*-\n")
    

    # Convert population file into smarter format
    # Note that this is safer than accumulating data in a list and dumping them in 
    # the end : the "with" loop creating thet txt file ensures that the data
    # are written on the fly and accummulated data are not lost if a crash occurs
    data = Table.read(summary_file_noext+".txt",format="ascii")
    data.write(summary_file_noext+'.csv', format="ascii.csv",overwrite="True")
    data.write(summary_file_noext+'.fits',format="fits",overwrite="True")
            

    # Open population file for the analysis chain and obs. time        

# for ch,dt  - analysis chain and observation time                                   

end_all = time.time()
print(" ****************** End of job - Total execution time = ",round((end_all-start_all)/60,1),"  minutes *****")
#--- Flux
#import gammapy.spectrum as sp
#import astropy.units as u
#from gammapy.utils.energy import EnergyBounds 
#
## Also availablle ; LogParabola, PowerLaw, ExponentialCutoffPowerLaw
#model = sp.models.PowerLaw(
#    index=2.0,
#    amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
#    reference=0.1 * u.TeV,)
#
#fit = sp.SpectrumFit(obs_list=grb.stack_obs, 
#                     model=model,
#                    fit_range=[0.03*u.TeV,0.2*u.TeV])
#fit.run()
#fitresult = fit.result
#ax0, ax1 = fitresult[0].plot(figsize=(8, 8))
#ax0.set_ylim(0, 20)
#print(fitresult[0])

#ebounds = EnergyBounds.equal_log_spacing(0.03, 0.5, 4, unit=u.TeV)
#seg = sp.SpectrumEnergyGroupMaker(obs=grb.stack_obs)
#seg.compute_groups_fixed(ebounds=ebounds)
#fpe = sp.FluxPointEstimator(obs=grb.stack_obs, 
#                            groups=seg.groups, 
#                            model=model)
#fpe.compute_points()
#fpe.flux_points.table
#for idx, obs in enumerate(grb.stack_obs):
#   print("GRB ",grb.name," - Observation: ",idx)
#   print("e_reco=")
#    seg.compute_groups_fixed(ebounds=ebounds)
#    print(seg.groups_from_obs())
#    flux=sp.FluxPointEstimator(obs,seg,model)
#    flux.compute_points()
#SpectrumObservation()


###############
# EVTDisplay  # ---------------------------------------------------------------
###############
# North
#☺irf_file = 'irf/Reco1-EvtDisp/North/CTA-Performance-North-20deg-average-100s_20170627.fits.gz'
#irf_file = 'irf/Reco1-EvtDisp/North/CTA-Performance-North-20deg-average-30m_20170627.fits.gz'
#irf_file = 'irf/Reco1-EvtDisp/North/CTA-Performance-North-NSBx05-20deg-average-100s_20170627.fits.gz'

# South
#irf_file = 'irf/Reco1-EvtDisp/South/CTA-Performance-South-20deg-100s_20170627.fits.gz'
#irf_file = 'irf/Reco1-EvtDisp/South/CTA-Performance-South-20deg-30m_20170627.fits.gz'

###############
# MARS        # ---------------------------------------------------------------
###############
# North
#irf_file = 'irf/Reco2-MARS/North/CTA-Performance-North-20deg-average-onaxis-100s_20170627.fits.gz'
#irf_file = 'irf/Reco2-MARS/North/CTA-Performance-North-20deg-average-onaxis-30m_20170627.fits.gz'
#irf_file = 'irf/Reco2-MARS/North/CTA-Performance-North-20deg-average-onaxis-05h_20170627.fits.gz'

# South
#irf_file = 'irf/Reco2-MARS/South/CTA-Performance-South-20deg-100s_IRFreco2_20170627.fits.gz'
#irf_file = 'irf/Reco2-MARS/South/CTA-Performance-South-20deg-30m_IRFreco2_20170627.fits.gz'
#irf_file = 'irf/Reco2-MARS/South/CTA-Performance-South-20deg-05h_IRFreco2_20170627.fits.gz'

