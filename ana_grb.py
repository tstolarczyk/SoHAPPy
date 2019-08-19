import matplotlib.pyplot as plt
#plt.rcdefaults()
#plt.rcParams['font.weight'] = 'bold'
#plt.rcParams['axes.labelweight'] = 'bold'
#plt.rcParams['axes.titleweight'] = 'bold'
#plt.style.use('seaborn-talk') # Make the lables readable
plt.style.use('seaborn-poster') # Make the lables readable
# print(plt.style.available)


# SET GAMMAPY_EXTRA to a decent value before starting
# to be updated
import os
import sys

# Windows is clever enough to replace / into the native \
# The r in front of the string menas that \ -used in windows- should be interpreted as such
os.environ['GAMMAPY_EXTRA']=r'../gammapy-extra-master'
os.environ['GAMMAPY_DATA'] =r'../gammapy-extra-master/datasets'
#print(os.getenv('GAMMAPY_EXTRA'))
#print(os.listdir(os.getenv('GAMMAPY_EXTRA')))

import numpy as np

# Gammapy
import gammapy
from gammapy.spectrum.models import Absorption

# Note : this whas been for a while in gammapy 0.6 and then was removed 
from cta_irf_imported import CTAPerf_imported # ...whereas this does work

# Utilities
#from config_handler import ConfigHandler
from grb_utils import GammaRayBurst
from plot_utils import PlotGammaRayBurst
from utilities import GetNameNoExtension

#------------------------------------------------------------------------------
### Initialisations
#------------------------------------------------------------------------------
niter = 1000 # Number of iterations for the Monte Carlo
alpha = 0.2 # 5 zones (1/0.2) to define the off region for background estimation
dbg_level = 0 # 0 : event counting, 1: event counting + some results, 2: + details for each event, 3: + plots 
saveplots = True

EBLmodel = "dominguez"
filepath = 'data/long_grb_test/LGRB1/'
#filepath = 'data/long_grb_test/LGRB2/'
#filepath = 'data/long_grb_test/LGRB3/'
#filepath = 'data/long_grb_test/LGRB4/'
#filepath = 'data/long_grb_test/LGRB5/'
#filepath = 'data/long_grb_test/LGRB6/'
#filepath = 'data/long_grb_test/LGRB7/'
#filepath = 'data/long_grb_test/LGRB8/'
#filepath = 'data/long_grb_test/LGRB9/'
#filepath = 'data/long_grb_test/LGRB10/'

# IRF files
###############
# EVTDisplay  # ---------------------------------------------------------------
###############
# North
#irf_file = 'irf/Reco1-EvtDisp/North/CTA-Performance-North-20deg-average-100s_20170627.fits.gz'
#irf_file = 'irf/Reco1-EvtDisp/North/CTA-Performance-North-20deg-average-30m_20170627.fits.gz'
irf_file = 'irf/Reco1-EvtDisp/North/CTA-Performance-North-NSBx05-20deg-average-100s_20170627.fits.gz'



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



#------------------------------------------------------------------------------
### Reads IRF, reads EBL absorption Let's go
irf_name = GetNameNoExtension(irf_file) # Extract irf name for further use
cta_perf = CTAPerf_imported.read(irf_file)
#cta_perf.peek()
absorption = Absorption.read_builtin(EBLmodel)

#------------------------------------------------------------------------------
### Create the GRB time-interval spectra from the given data file
grb = GammaRayBurst.from_file(filepath=filepath, absorption=absorption,dbg=dbg_level)

#------------------------------------------------------------------------------
# Analysis summary
print("")
print("==========================  GammaPy ",gammapy.__version__)
print("  Analysing files in     : ",filepath)
print("    - GRB name           : ",grb.name)
print("    - redshift         z : ",grb.z)
print("  EBL model              : ",EBLmodel)
print("  Response function      : ",irf_name)
print("  alpha (on/off regions) : ",alpha)
print("  Number of MC trials    : ",niter)
print("  Debug mode             : ",dbg_level)
print("==========================")


#------------------------------------------------------------------------------
# Run simulations of the current grb
# Compute significance on the full observation, 
# significance max on intervals,
# Time of siginificance max, time of 3 sigma reached


sigmax_list = []
t_sigmax_list = []
t_3sigma_list = []
sigfull_list = []

sigma_sum = 0
sigma2_sum = 0 

iMC=1
prt_frequency = 10

while(iMC <= niter):
    if (niter <= 10) or (np.mod(iMC,prt_frequency) == 0): 
        print(" ----------------------------------------- Simulation ",iMC)
    grb.run_simulation(cta_perf,alpha=alpha)
    obs_stat = grb.get_cumulative_stats()

    # Best sensitivity
    sigmax     = max(obs_stat['sigma'])
    t_sigmax   = obs_stat['livetime'][ (np.where( obs_stat['sigma']==sigmax)[0][0]) ]    
    sigmax_list.append(sigmax)
    t_sigmax_list.append(t_sigmax)
    if (dbg_level) : print(" --- Max significance = ",sigmax," at tmax=",t_sigmax)
    
    # 3 sigma reached
    t_3sigma = obs_stat['livetime'][np.where(obs_stat['sigma']>=3)[0][0]]
    t_3sigma_list.append(t_3sigma)
    if (dbg_level) : print(" --- Alert time (3 sigma) =",t_3sigma)
    
    # Overall sensitivity
    sig_full = obs_stat['sigma'][len(obs_stat['sigma'])-1]
    sigfull_list.append(sig_full)
    if (dbg_level) : print(" --- Overal signif. = ",sig_full)
    
    # Compute mean significance at each time intervals
    sigma_sum += obs_stat['sigma']
    sigma2_sum += obs_stat['sigma']**2
    

    iMC+=1

#--- Observation investigation
if (dbg_level > 2) : grb.quicklook(plot=True )

if (niter<=1):
    print(grb)
    PlotGammaRayBurst.plot_stats_detection(grb, savefig=True, outdir='./out/')
    #PlotGammaRayBurst.make_gif_from_models(grb, savefig=True, outdir='./out/')
    plt.show(block=False)



##--- PLOTS (to be moved)
fig = plt.figure(figsize=(18,18))

# Mean significance dsitribution
sigma_mean = sigma_sum/niter
sigma2_mean = sigma2_sum/niter
sigma_rms = np.sqrt(sigma2_mean-sigma_mean**2)
#+print("Mean sigma =",sigma_mean,"rms=",sigma_rms)

a1 = plt.subplot(221)
a1.set_xlabel('Observation duration (s)')
a1.set_ylabel('Mean significance')
a1.set_title('Significance for '+str(niter)+' realisations')
a1.grid(which='both')
plt.errorbar(obs_stat['livetime'],sigma_mean,yerr=sigma_rms,fmt='.') # Use last simulation time interval - they are all the same

# Plot max significance
a2 = plt.subplot(222)
a2.set_xlabel('Standard deviation (Li&Ma)')
a2.set_ylabel('#Realisation')
a2.set_title('Max. significance- Mean = '
             + str( round(np.mean(sigmax_list),1))
             + ' +/- '
             + str( round(np.std(sigmax_list),1)))
a2.grid(which='both')
plt.hist(sigmax_list,color="grey")
#y,binEdges = np.histogram(sigmax_list)
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#plt.bar(bincenters, y, yerr=np.sqrt(y))

# Plot time to get 3 sigma
a3 = plt.subplot(223)
a3.set_xlabel('Observation duration (s)')
a3.set_ylabel('#Realisation')
a3.set_title('Time to claim 3s - Mean = '
             + str( round(np.mean(t_3sigma_list),1))
             + ' +/- '
             + str( round(np.std(t_3sigma_list),1)))
a3.grid(which='both')
plt.hist(t_3sigma_list, color="grey")


# Plot overall significance
a4 = plt.subplot(224)
a4.set_xlabel('Standard deviation (Li&Ma)')
a4.set_ylabel('#Realisation')
a4.set_title('Overall significance - mean = '
             + str( round(np.mean(sigfull_list),1))
             + ' +/- '
             + str( round(np.std(sigfull_list),1)))
a4.grid(which='both')
plt.hist(sigfull_list,color="grey")
plt.show()

if (saveplots):
    filename = grb.name+'_'+irf_name+'-'+EBLmodel+'-'+str(niter)+'iter'
    fig.savefig(filename+'.eps',transparent=True)
    fig.savefig(filename+'.jpg')


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