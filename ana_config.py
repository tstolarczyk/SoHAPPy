# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import astropy.units as u
import gammapy
"""
This module contains most of the steering parameters of SoHAPPy.
"SoHAPPy -h" in a shell gives the list of parameters that can be superseded
on the command line
"""
#-----------------------------------------------------------------------------#
# Events to be processed
#-----------------------------------------------------------------------------#
# These values can be superseded by the command line

# ifirst : either the first event, then ngrb are processed, or a list of events
# See end of file for some interesting case
# ifirst=85 # 4 slices South, det. @3s, not at 24°
#ifirst=815 # 29 slices, N
#ifirst = 54
#ifirst=3
# ifirst=10 # 2 independent slices merged in B
# ifirst = 901
# ifirst=16
#ifirst = [3,85,815]
#ifirst = 1
# ifirst=398
# ifirst = 751 - sent for comparison with ctools
#ifirst  = 54
ifirst = 2
ngrb    = 1 # Number of GRB to be read - if 1 specific actions on one grb
seed    = 'random-seed' # use a fix number to get a fixed sequence sequence
#res_dir = "../output/testnew" # Output directory for results
res_dir = "D:/000_Today/testnew" # Output directory for results

#-----------------------------------------------------------------------------#
# Analysis
#-----------------------------------------------------------------------------#
niter     = 1  # Number of iterations for the Monte Carlo (1: no fluctuate)
method    = 0  # Zero : apert. photometry - the only analysis implemented
obs_point = "end" # Define the position of observation in the time slice

do_fluctuate   = False # (True) Statistical fluctuations are enabled,
do_accelerate  = True  # (True) Simulation skipped if no detection in first 10%
#-----------------------------------------------------------------------------#
# Debugging
#-----------------------------------------------------------------------------#
# if negative or zero : no plot
# 0: evt counting, 1: + some results, 2: + event details
# 3: Plots for each trials - create an animated gif - heavy
dbg    = 0
silent = False # If True, nothing written on screen (output to log)

#-----------------------------------------------------------------------------#
# Input
#-----------------------------------------------------------------------------#
grb_dir    = '../input/lightcurves/LONG_FITS/'
if (gammapy.__version__ == "0.12"):
    irf_dir = "../input/irf/OnAxis/prod3-v2"
if (gammapy.__version__ == "0.17"):
    # irf_dir    = "../input/irf/Full/prod3-v2"
    irf_dir = "D:\CTA\Analyse\SoHAPPY-IRF\prod3-v2"

# prompt tests
# redshift = 4.0 # Prompt test
# res_dir    = "D:/000_Today/prompt-visible-nofluctuation-z"+str(redshift)
#-----------------------------------------------------------------------------#
# Ouput
#-----------------------------------------------------------------------------#
datafile   = "data.txt" # Main output file (population study)
logfile    = "analysis.log"

# Special outout action
save_simu      = False  # (False) Simulation saved to file for offline use
save_grb       = False # (False) GRB saved to disk -> use grb.py main
save_dataset   = False  # Will be set to False if more than one MC trial
write_slices   = False # (False) Store detailed information on slices if True
remove_tarred  = False # Remove tarred files, otherwise keep for faster access

#-----------------------------------------------------------------------------#
# Visibility and Physics
# The absolute maximum repointing times for the CTA telescopes (to and from
# anywhere in the observable sky) will be 50 s for the LSTs and 90 s for the
# MSTs and SSTs, with the goal to reach shorter slewing times
# (Ref. CTA Science document).
#-----------------------------------------------------------------------------#
EBLmodel  = "dominguez" # Extragalactic Background Light model
#EBLmodel  = "built-in" # Default built-in EBL model

#arrays    = {"North":"FullArray", "South":"FullArray"} # "FullArray", "LST",...
arrays   = {"North":"FullArray", "South":"MST"} # "FullArray", "LST",...
#dtslew  = {"North":30*u.s, "South":30*u.s} # Maximum slewing time delay
dtslew  = {"North":30*u.s, "South":60*u.s} # Maximum slewing time delay
dtswift = 77*u.s # Fixed SWIFT latency (The mean value is 77*u.s)

fixslew   = True   # If False,generate a random delay < dtslew
fixswift  = True   # If False, the latency is generated from real data file
swiftfile = "../input/swift/Swift_delay_times.txt" # Swift latency file

altmin    = 24*u.deg # Minimum altitude (original default is 10 degrees)
altmoon   = 90*u.deg # Moon maximum altitude
moondist  = 0*u.deg  # Moon minimal distance
brightness= 1        # Moon maximum brigthness

newvis    = True   # (True) If True, visibility windows are recomputed
depth     = 3*u.day # Maximum duration after trigger to compute the visibility.
skip      = 0       # If 1 skip first obsercation night, if zero consider it

test_prompt    = False # (False) If True test prompt alone (experimental)
use_afterglow  = False # Prompt characteritics from the afterglow with same id.
signal_to_zero = False # (False) Keep only background, set signal to zero

fixed_zenith   = False # If set to a value ("20deg") freeze zenith in IRF
magnify        = 1 # Multiplicative factor of the flux, for tests (def is 1)
#-----------------------------------------------------------------------------#
# Some interesting events to be considered
#-----------------------------------------------------------------------------#
#-------------------------- Afterglows

### Interesting cases
#ifirst = 188  # Mammuth, 100s @10° and 24°(2nd night alone not det.)
#ifirst = 815 # By far the brightest one, not seen in South
#ifirst = 457 # Very low altitude in North - good to check altmin
#ifirst = 6   # Not vis. @24° -  @10° N: 3.5@30s, 2.9@107s, , S:not vis.
#ifirst = 85  # Has 4 windows over 2 days for altmin=, not vis. for 24°
#ifirst = 10  # Not detected N, S, B. A few slices only
#ifirst = [54,416,457,926] # These ones have large slice numbers in N and S

### North and South, not detected
#ifirst = 64  # Seen N&S, very few slices, but not detected
#ifirst =3    # Seen N&S, few slices, not detected
#ifirst =901  # Seen N&S, Nice overlap, not detected

### North and South, detected
#ifirst = 204 # @10deg: N&S, B:12s, almost perfect overlap @24°:not seen S
#ifirst = 54 # N:197s S:21s B:198s

### Not detected N&S, almost detected B - would a likelihood help ?
#ifirst = 57  # Not detected N, S separately. Few slices  B:2.8S @10° 2.3S@24°

### Others
#ifirst = 398 # 5 sigma detected in North, 29 slices. Not visible in South

### Bug reknown
#  In astroplan v0.7.dev0, risetime is not found because too close from trigger
# ifirst=[233, 644, 769]
#-------------------------- Prompt

#ifirst = [2, 6, 12, 30, 139, 172, 191, 278, 490, 506] # Prompt test files
#ifirst=12