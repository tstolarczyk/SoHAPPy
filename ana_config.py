# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import os
import astropy.units as u
# Suppress RunTimeWarning in mc_3d background fitting
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


#------------------------------------------------------------------------------
### Code Initialisations
#------------------------------------------------------------------------------
dbg_level = 0 # 0 : evt counting, 1: evt counting + some results, 2: + details for each evt, 3: + plots
showplots = 2 # 0 : not shown, not written,
              # 1 : shown on screen
              # 2 : shown on screen and written out
              # 3 : written out, not shown on screen (for production)
likelihood = 2 # If zero, do the on-off analysis, otherwise set the degree fof freedom
#------------------------------------------------------------------------------
### Simulation Initialisations
#------------------------------------------------------------------------------
niter = 10 # Number of iterations for the Monte Carlo
alpha = 0.2 # 5 zones (1/0.2) to define the off region for background estimation
binsize  = 0.125*u.deg
fov      = 5.*u.deg # full width (i.e. if centered in zero, -fov/2, fov/2)
# for quick tests
#binsize  = 0.5*u.deg
#fov      = 2.5*u.deg # full width (i.e. if centered in zero, -fov/2, fov/2)
#------------------------------------------------------------------------------
### Input
#------------------------------------------------------------------------------
res_folder = "../output/Result" # Folder for results
res_folder = "../output/Today" # Folder for results
irf_folder = "../input/irf"
grb_folder = '../input/lightcurves/long_grb_test/' # to be changed into a repository to be scanned

os.environ['GAMMAPY_EXTRA'] =r'../input/gammapy-extra-master'
os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
#print(os.getenv('GAMMAPY_EXTRA'))
#print(os.listdir(os.getenv('GAMMAPY_EXTRA')))
#------------------------------------------------------------------------------
### Physics Initialisations
#------------------------------------------------------------------------------
det_level       = 0.9 # Will be detected if n=3,5 sigma reached in det_level of the population
#redfactor       = [1, 2, 3, 5, 10, 20, 50, 100] # Initial flux is divided by this value, for test
#redfactor       = [0.5,1, 2] # Initial flux is divided by this value, for test
redfactor = [1]
#grblist         = ["LGRB1","LGRB2","LGRB3","LGRB4","LGRB5","LGRB6","LGRB7","LGRB8","LGRB9","LGRB10"]
grblist         = ["LGRB1","LGRB2","LGRB3","LGRB4","LGRB5"]
grblist         = ["LGRB1"]
EBLmodel        = "dominguez"

#------------------------------------------------------------------------------
### Instrument response functions
#------------------------------------------------------------------------------
site        ='Roque de los Muchachos'
obstime     = '2000-01-01 22:00:04'

chain       = ["Reco1"] # or Reco2
hemisphere  = ["North"] # or South
#hemisphere  = ["South"] # or South
#zenith      = ["20","40","60"] # or 40 deg...
zenith      = ["20"] # or 40 deg...
azimuth     = ["average"]  # or N, or S
#duration    = ["100s","30m"] # or 30m, 05h 50h
duration    = ["100s"] # or 30m, 05h 50h
NSB         = False # if True try to find a high NSB file
