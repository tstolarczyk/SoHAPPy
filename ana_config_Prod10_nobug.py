# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import os
import astropy.units as u

#------------------------------------------------------------------------------
### Physics Initialisations
#------------------------------------------------------------------------------
ngrb            = 1000 # Number of GRB to be read
ifirst          = 1 # First GRB to be read
det_level       = 0.9 # Will be detected if n=3,5 sigma reached in det_level of the population
#redfactor       = [1, 2, 3, 5, 10, 20, 50, 100] # Initial flux is divided by this value, for test
#redfactor       = [0.5,1, 2] # Initial flux is divided by this value, for test
redfactor = [1]
EBLmodel        = "dominguez"
#------------------------------------------------------------------------------
### Simulation Initialisations
#------------------------------------------------------------------------------
niter  = 100 # Number of iterations for the Monte Carlo
alpha  = 0.2 # 5 zones (1/0.2) to define the off region for background estimation
fov    = [5.*u.deg, 0.125*u.deg] # Full Field-of-view and bin size
dtslew = [30*u.s, True] # Time to point the GRB,  If False, random < dt_slew
# for quick tests
#binsize  = 0.5*u.deg
#fov      = 2.5*u.deg # full width (i.e. if centered in zero, -fov/2, fov/2)
#------------------------------------------------------------------------------
### Code Initialisations
#------------------------------------------------------------------------------
dbg_level = 0  # 0 : evt counting, 1: evt counting + some results, 2: + details for each evt, 3: + plots
showplots = 0 # 0 : not shown, not written,
               # 1 : shown on screen
               # 2 : shown on screen and written out
               # 3 : written out, not shown on screen (for production)
lkhd      = 0 # If zero, do the on-off analysis, otherwise set the degree of freedom
#------------------------------------------------------------------------------
### Input
#------------------------------------------------------------------------------
res_folder = "../output/Result" # Folder for results
res_folder = "../output/Today" # Folder for results
res_folder = "../output/Prod10_nobug"

old_file   = False # Use old GRB file
grb_oldfolder = '../input/lightcurves/long_grb_test/' # to be changed into a repository to be scanned
grb_folder    = '../input/lightcurves/LONG_FITS'
os.environ['GAMMAPY_EXTRA'] =r'../input/gammapy-extra-master'
os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
#print(os.getenv('GAMMAPY_EXTRA'))
#print(os.listdir(os.getenv('GAMMAPY_EXTRA')))


#------------------------------------------------------------------------------
### Instrument response functions
#   Default values to force the source observation conditions
#------------------------------------------------------------------------------
fixed       = False # Unused - If True the source observation conditions are fixed
irf_folder  = "../input/irf"
#chain       = ["Reco1"] # or Reco2 (Reco1 = eventdisplay, Reco2 = MARS)
#hemisphere  = ["North"] # or South
#zenith      = ["20deg"] # ["20deg","40deg","60deg"]
#azimuth     = ["average"]  # or N, or S
#duration    = ["100s"] # or 30m, 05h 50h
#NSB         = False # if True try to find a high NSB file
