# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import astropy.units as u
"""
This module contains most of the steering parameters of SoHAPPy.
"SoHAPPy -h" in a shell gives the list of parameters that can be superseded
on the command line 
"""

#-----------------------------------------------------------------------------#
# These values can be superseded by the command line
niter  = 100 # Number of iterations for the Monte Carlo
ngrb   = 1 # Number of GRB to be read - if 1 specific studies on one grb

# First GRB - can be a list then ngrb is useless
# Prompt
#ifirst = [2, 6, 12, 30, 139, 172, 191, 278, 490, 506] # Prompt test files 
#ifirst=12

# Afterglows
ifirst=457 # Remains at very low altitude in North - good to check altmin
ifirst = 6   # Seen in North, not seen in South
ifirst = 85
#ifirst=3
#ifirst = [54,416,457,926] # This ones have large slice number in N and S
#ifirst=901 # Nice overlap
#ifirst=204 
#ifirst = [11, 12]
#ifirst = [188,306,458,454,765,709]
#ifirst = 57  # Not detected N and S separately. A few slices only -> Check N+S
#ifirst = 10  # Not detected N and S separately. A few slices only
#ifirst = 188 # Mammuth, detected the second day also, large site overlap 
#ifirst = 13  # Not detected
#ifirst=64 # Detected in N and S with very few slices - good for debug !!!
#ifirst = 815 # By far the brightest one, not seen in South
#ifirst = 398 # 5 sigma detected in North, 29 slices. Not visible in South
#ifirst =399 # Not visible in North nor in South
#ifirst = 54
#ifirst = 913

# Debugging : if negative or zero : no plot
# 0: evt counting, 1: + some results, 2: + event details
# 3: Plots for each trials - create an animated gif - heavy
dbg    = 1


# prompt tests
# redshift = 4.0 # Prompt test
# res_dir    = "D:/000_Today/prompt-visible-nofluctuation-z"+str(redshift)

res_dir    = "D:/000_Today/testnew"
#-----------------------------------------------------------------------------#
# Physics
EBLmodel        = "dominguez"
#-----------------------------------------------------------------------------#
# Detection
dtslew    = 0*u.s # Time delay to point the GRB, usually 30*u.s 
fixslew   = True   # If False, random < dt_slew
dtswift   = 0*u.s # Mean SWIFT latency : mean value is 77*u.s
fixswift  = True   # If False, the latency is generated from real data file
swiftfile = "../input/swift/Swift_delay_times.txt" # Real data file
altmin    = 10*u.deg # Minimum altitude (original default is 10 degrees)
#-----------------------------------------------------------------------------#
# Analysis
lkhd      = 0     # If zero : on-off analysis; otherwise degrees of freedom
obs_point = "end" # Define the position of observation in the time slice
#-----------------------------------------------------------------------------#
# Development
silent         = False  # If True, nothing written on screen (output to log)
save_simu      = False # (False) Simulation svaed to file for offline use
save_grb       = False # (False) GRB saved to sik -> use grb.py main
write_slices   = False # (False) Store detailed information on slices if True
signal_to_zero = False # (False) Keep only background, set signal to zero
do_fluctuate   = False # (True) Statistical fluctuations are enabled 
do_accelerate  = True # (True) Simulation skipped if no detection in first 10%  
fixed_zenith   = False # If set to a value ("20deg") freeze zenith in IRF 
remove_tarred  = False # Remove tarred files, otherwise keep for faster access
test_prompt    = False # (False) If True test prompt alone (experimental)
use_afterglow  = False # Prompt characteritics from the afterglow with same id.
get_visibility = True  # (True) If True, visibility windows are recomputed
day_after      = 0   # Add up to that days (if negative only that day) 
#-----------------------------------------------------------------------------#
# Input-Ouput
old_file   = False # Use old GRB file
grb_olddir = '../input/lightcurves/long_grb_test/'
grb_dir    = '../input/lightcurves/LONG_FITS/'
irf_dir    = "../input/irf/OnAxis/"
datafile   = "data.txt"
logfile    = "analysis.log"

