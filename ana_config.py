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
ngrb   = 10 # Number of GRB to be read - if 1 specific studies on one grb
res_dir    = "D:/000_Today/testnew"
#-----------------------------------------------------------------------------#
# Physics
EBLmodel        = "dominguez"
#-----------------------------------------------------------------------------#
# Detection
dtslew    = 30*u.s # Time delay to point the GRB, usually 30*u.s
fixslew   = True   # If False, random < dt_slew
dtswift   = 0*u.s # Mean SWIFT latency : mean value is 77*u.s
fixswift  = True   # If False, the latency is generated from real data file
swiftfile = "../input/swift/Swift_delay_times.txt" # Real data file
altmin    = 10*u.deg # Minimum altitude (original default is 10 degrees)
newvis    = True  # (True) If True, visibility windows are recomputed
depth     = 3*u.day # Duration after trigger to compute the visibility.
skip      = 0 # If 1 skip first night, if zero consider it

ifirst=1
###-----------
### Afterglows
###-----------

### Interesting cases
#ifirst = 188  # Mammuth, 100s @10° and 24°(2nd night alone not det.)
#ifirst = 815 # By far the brightest one, not seen in South
#ifirst = 457 # Very low altitude in North - good to check altmin
#ifirst = 6   # Not vis. @24° -  @10° N: 3.5@30s, 2.9@107s, , S:not vis.
#ifirst = 85  # Has 4 windows over 2 days for altmin=, not vis. for 24°
#ifirst = 10  # Not detected N, S, B. A few slices only
#ifirst = [54,416,457,926] # These ones have large slice numbers in N and S

### North and South, not detected
#ifirst = 64  # Seen N&S, very few slcies, but not detected
#ifirst =3    # Seen N&S, few slices, not detected
#ifirst =901    #Seen N&S, Nice overlap, not detected

### North and South, detected
#ifirst = 204 # @10deg: N&S, B:12s, almost perfect overlap @24°:not seen S
#ifirst = 54 # N:197s S:21s B:198s

### Not detected N&S, almost detected B - would a likelihood help ?
#ifirst = 57  # Not detected N, S separately. Few slices  B:2.8S @10° 2.3S@24°

### Others
#ifirst = 398 # 5 sigma detected in North, 29 slices. Not visible in South

###-------
### Prompt
###-------
#ifirst = [2, 6, 12, 30, 139, 172, 191, 278, 490, 506] # Prompt test files
#ifirst=12

# Debugging : if negative or zero : no plot
# 0: evt counting, 1: + some results, 2: + event details
# 3: Plots for each trials - create an animated gif - heavy
dbg    = 2

# prompt tests
# redshift = 4.0 # Prompt test
# res_dir    = "D:/000_Today/prompt-visible-nofluctuation-z"+str(redshift)


#-----------------------------------------------------------------------------#
# Analysis
lkhd      = 0     # If zero : on-off analysis; otherwise degrees of freedom
obs_point = "end" # Define the position of observation in the time slice
#-----------------------------------------------------------------------------#
# Development
silent         = False # If True, nothing written on screen (output to log)
save_simu      = False # (False) Simulation svaed to file for offline use
save_grb       = False # (False) GRB saved to disk -> use grb.py main
write_slices   = False # (False) Store detailed information on slices if True
signal_to_zero = False # (False) Keep only background, set signal to zero
do_fluctuate   = False # (True) Statistical fluctuations are enabled
do_accelerate  = True  # (True) Simulation skipped if no detection in first 10%
fixed_zenith   = False # If set to a value ("20deg") freeze zenith in IRF
remove_tarred  = False # Remove tarred files, otherwise keep for faster access
test_prompt    = False # (False) If True test prompt alone (experimental)
use_afterglow  = False # Prompt characteritics from the afterglow with same id.
#-----------------------------------------------------------------------------#
# Input-Ouput
old_file   = False # Use old GRB file
grb_olddir = '../input/lightcurves/long_grb_test/'
grb_dir    = '../input/lightcurves/LONG_FITS/'
irf_dir    = "../input/irf/OnAxis/"
datafile   = "data.txt"
logfile    = "analysis.log"