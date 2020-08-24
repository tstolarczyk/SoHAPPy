# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import astropy.units as u

#-----------------------------------------------------------------------------#
# These values can be superseded by the command line
niter  = 10 # Number of iterations for the Monte Carlo
#glist = [188,306,458,454,765,709]
ngrb = 1  # Number of GRB to be read
ifirst=1
#ifirst = 57 # Not detected N and S separately
#ifirst = 6 # Not detected
#ifirst = 10 # first GRB
ifirst = 188 # Mammuth, detected the second day also
#ifirst = 13 # Not detected
#ifirst = 815 # First GRB to be read
#ifirst = 398 # First GRB to be read - 29 slices
#ifirst = 100 # First GRB to be read
ifirst = 64
# Debugging : if negative or zero : no plot
# 0: evt counting, 1: + some results, 2: + event details
# 3: Plots for each trials - create an animated gif - heavy
dbg    = -2
#-----------------------------------------------------------------------------#
# Physics
EBLmodel        = "dominguez"
#-----------------------------------------------------------------------------#
# Detection
dtslew   = 30*u.s # Time delay to point the GRB,  
fixslew  = True   # If False, random < dt_slew
#-----------------------------------------------------------------------------#
# Analysis
lkhd      = 0   # If zero : on-off analysis; otherwise degrees of freedom
#-----------------------------------------------------------------------------#
# Development
save_simu      = False # (False) Simulation svaed to file for offline use
save_grb       = False  # (False) GRB class content is saved for offfline use
write_slices   = False # (False) Store detailed information on slices if True
signal_to_zero = False # (False) Keep only background, set signal to zero
do_fluctuate   = False  # (True) Statistical fluctuations are enabled 
do_accelerate  = False # (True) Simulation skipped if no detection in first 10%  
fixed_zenith   = False # If set to a value ("20deg") freeze zenith in IRF 
day_after      = 0     # Days to add to the simulation (default is 0)  

#-----------------------------------------------------------------------------#
# Input-Ouput
old_file   = False # Use old GRB file
grb_olddir = '../input/lightcurves/long_grb_test/'
grb_dir    = '../input/lightcurves/LONG_FITS'
irf_dir    = "../input/irf/OnAxis/"
# res_folder  = "../output/Result" # Folder for results
res_dir    = "../output/Today" # Folder for results
# res_folder  = "../output/Prod10_nobug"