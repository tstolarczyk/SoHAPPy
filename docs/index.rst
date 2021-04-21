#######
SoHAPPy
#######

** Simulation of High-energy Astrophysics Process in Python **

SoHAPPy computes the response of CTA to a population of GRB provided as fits
files. Each GRB is a set of time slices with their own respective spectral model.
The GRB position is assumed to be known and has a fixed position in the
field of view.

Repository structure
====================
The repository has 2 subfolders:

* analysis: contains notebooks and scripts to analyse the SoHAPPy outputs.
* docs: contains the Sphinx information to generate the documentation from the code

The irf and GRB files are expected to be obtained from specific folders.
The output consists in a few files packed together into a tar.gz archive.
They are written in a specific output folder

Launching the code
-----------------
The steering parametres are obtained from the ana_config.py module which is imported wherever it is necessary. Many of them are not intended to be relevant but for development and debugging.
The most important parametres can be passed on the command line (they supersede the ana_config.py parameters which are considered as the default).

Type:
``python SoHAPPy.py -h``
to have the list explicited (to be completed).

Concepts
--------
SoHAPPy takes fits file as input:

* GRB source file

GRB files, either afterglow or prompt, contains a series of flux measurement at certain energies for a series of time slices. This files also give the position, the redshift, some physical information and the default visibility over one night for a minimal altitude of 10 degrees (this visibility can be recomputed, see later)

Afterglow files have time slices from the physical trigger time to typically 2 days.

Prompt files... (to be described)

The GRB information is handled by the ``GammaRayBurst`` class in the ``grb.py`` module.
Each ``GammaRayburst`` object carries a list in time of ``AbsorbedSpectralModel`` models obtained from template models, which are interpolation of the initial flux points after applying a given EBL absorption.
Alternatively the initial file can contain already absorbed spectra and in that case the model as to be declared "built-in" in ana_congif.py.

* Time slots
A so-called time slot, the class::Slot:: in timeslot.py, is the key element of the SoHAPPy analysis.
A Slot contains a GRB object, a time delay, the start and end time of the slot, the site(s) of the observation and a collection of Slices.

The native slot is composed ot the GRB initial time slices. After application of delays or visibilities, it contains the effective time slices on which the analysis can be done. It can handle observations on one or several sites.

The Slot class has the tools to modifiy the initial slot into the slot attached to the GRB at a given site (function apply_visibility and both_sites). Slots can be merged with the merge function.

Initial naked slots that handles the time slices are dressed with physical information : this has for effect to associated the correct (or best) IRF and spectrum to each of the slices.

* Time slices
A slice carries a start an stop time, an observation time, one or several IRFs, a spectrum identifier from the GRB class spectra, a site tagging (A slice can belong to more than one site).
The observation time can be set at the start, stop or mean time of the slice.

The default is to set the observation at the "end" of the slice (can be changed in ana_config.py).
The associated IRFs are obtained according to the altitude (zenith) and slice duration at the beginning of the slice (see recommendation on how to choose the IRF).
The flux is the flux from the GRB which is the closest in time to the observation time

* Monte Carlo
The Monte Carlo class (MonteCarlo in mcsim.py) is the deepest class of SoHAPPy. it run simulations over a slot, perform the analysis, obtain the results.
It is associated to the mcsim_config.py and mcsim_res.py modules that contain simulation parameters and output tools respectively.
The MonteCarlo class carries a Slot instantiation and the following resulting variables:
- the maximum siginificance reached and the corresponding observation time and excess and background counts
- the time at which 3 and 5 sigma were reached and the associated excess and background count numbers.
- the evolution of the siginificance along time
- a list of slices for which warning were issued for further debugging (not intended to the normal user).

Running the simulation generate for a given number of iterations, on and off counts, fluctuate them, get the max siginificance and when 3 sigma and 5 sigma are reached. The user can request to not fluctuate the counts, which automatically set the number of iteration to one.


* Visibilities


What does SoHAPPy do ?
The various steps can in the main (SoHAPPY.py) are the following:

1. read the steering parmateres from ana_config.py and/or the command line

2. create the GRB list to be analysed (either a large number for a populations tudy or only a few for specific spectral analysis)

3. Open the log output file and make a copy of the ana_config file for further reference

4. Remind the simulation/analysis parameters, archive to the log

5. Loop over the GRB list
	5.1 Read the GRB information, create a GammaRayBurst object (list of interpolated abosrbed spectrum)
	5.2 Create a time slot with the original GRB time slices
	5.3 Update the default visibility, if requested, extending the observation window beyond one day for any given minimal altitude.
	5.4 Printout an plot a ceertain number of information managed so far. Save the GRB class to disk if requested

	5.5 For all involved sites (N, S)
		5.5.1 Create a Monte Carlo object


The steering parameters
-----------------------

Main parameters are defined in ana_config.py
The module ana_config.py  contains most of the steering parameters of SoHAPPy.
"SoHAPPy -h" in a shell gives the list of parameters that can be superseded
on the command line

+-----------+-------------------+---------------------------------------------+
| variable  | Default           | What is it ?                                |
+===========+===================+=============================================+
| Events to be processed
+-----------+-------------------+---------------------------------------------+
| ifirst    | 1                 | the first GRB file to be processed          |
+-----------+-------------------+---------------------------------------------+
| ngrb      | None              | Number of GRB to be processed               |
|           |                   | If 1, special actions are taken             |
+-----------+-------------------+---------------------------------------------+
| seed      | 'random-seed'     | Choose a fixed integer for a fixed sequence |
+-----------+-------------------+---------------------------------------------+
| res_dir   | ../output         | The output folder                           |
+-----------+-------------------+---------------------------------------------+
| variable  | Default           | What is it ?                                |
+-----------+-------------------+---------------------------------------------+
| variable  | Default           | What is it ?                                |
+-----------+-------------------+---------------------------------------------+

#-----------------------------------------------------------------------------#
# Analysis
#-----------------------------------------------------------------------------#
niter     = 1  # Number of iterations for the Monte Carlo (1: no fluctuate)
method    = 0  # Zero : apert. photometry - the only analysis implemented
obs_point = "end" # Define the position of observation in the time slice

#-----------------------------------------------------------------------------#
# Input-Ouput
#-----------------------------------------------------------------------------#

# Debugging : if negative or zero : no plot
# 0: evt counting, 1: + some results, 2: + event details
# 3: Plots for each trials - create an animated gif - heavy
dbg    = 0

# Data and result files
old_file   = False # Use old GRB file if True
grb_olddir = '../input/lightcurves/long_grb_test/' # Old GRB folder
grb_dir    = '../input/lightcurves/LONG_FITS/'
datafile   = "data.txt" # Main output file (population study)
logfile    = "analysis.log"
if (gammapy.__version__ == "0.12"):
    irf_dir = "../input/irf/OnAxis/prod3-v2"
if (gammapy.__version__ == "0.17"):
    # irf_dir    = "../input/irf/Full/prod3-v2"
    irf_dir = "D:\CTA\Analyse\SoHAPPY-IRF\prod3-v2"
array = {"North":"FullArray", "South":"FullArray"} # "FullArray", "LST",...

# prompt tests
# redshift = 4.0 # Prompt test
# res_dir    = "D:/000_Today/prompt-visible-nofluctuation-z"+str(redshift)

# Special outout action
silent         = False # If True, nothing written on screen (output to log)
save_simu      = False  # (False) Simulation saved to file for offline use
save_grb       = False # (False) GRB saved to disk -> use grb.py main
save_dataset   = False  # Will be set to False if more than one MC trial
write_slices   = False # (False) Store detailed information on slices if True
remove_tarred  = False # Remove tarred files, otherwise keep for faster access

#-----------------------------------------------------------------------------#
# Detection and Physics
#-----------------------------------------------------------------------------#
EBLmodel  = "dominguez" # Extragalactic Background Light model
#EBLmodel  = "built-in" # Default built-in EBL model

test_prompt    = False # (False) If True test prompt alone (experimental)
use_afterglow  = False # Prompt characteritics from the afterglow with same id.

The absolute maximumrepointing times for the CTA telescopes (to and from anywhere in the observable sky) will be 50 s for the LSTs and 90 s for the MSTs and SSTs, with the goal to reach shorter slewing times. Moreover the distance of the GRB to the presenthmy observed source can be small. The typical value of the LST slewing is set to 30 s, and 60 s for the MST.

dtslew    = {"North":30*u.s,"South":30*u.s} # Maximum time delay to point the GRB, usually 30*u.s
fixslew   = True   # If False,generate a random delay < dtslew
dtswift   = 77*u.s # Fixed SWIFT latency (The mean value is 77*u.s)
fixswift  = True   # If False, the latency is generated from real data file
swiftfile = "../input/swift/Swift_delay_times.txt" # Swift latency file
altmin    = 24*u.deg # Minimum altitude (original default is 10 degrees)
newvis    = True   # (True) If True, visibility windows are recomputed
depth     = 3*u.day # Maximum duration after trigger to compute the visibility.
skip      = 0       # If 1 skip first obsercation night, if zero consider it
signal_to_zero = False # (False) Keep only background, set signal to zero
do_fluctuate   = False # (True) Statistical fluctuations are enabled,
do_accelerate  = True  # (True) Simulation skipped if no detection in first 10%
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

* do_accelerate
When True, the simulation is stopped ifnone of the firs trials in the limit on
1 - det_level have rzeached hthe minimal siginificance. In a simulation with 100 trials requesting a 3 sigma detection in 90% of the iterations,the simulation will be stopped after 10 iteration if none of the trial have rzached 3 sigma.

Note that this bias he resulting popualtion since it articiially deplete the max significance population below the minimum required (e.g. 3 sigma).

Secondary parameters : mcsim_config.py

This parameters are in principle not to be changed.

The effect of the acceleration has been stidied in detail and has of course no effect on the population above 3 and 5 sigma.

* do_fluctuate
If False, the code runs the siumaltion once withut fluctuation. For bright sources with max significance much above 5 sigma this does not bais the result.
Tyhe subpopualtion near 3 and 5 sigma is biassed

DESCRIBE THE BIAS



Introduction
------------
On-off method
Spectral analysis is done outsid ethe code

.. toctree::
	:maxdepth: 2

    irf.rst

.. comment automodapi:: SoHAPPy

.. automodapi:: SoHAPPy
   :include-all-objects:
   :no-inheritance-diagram:

#.. automodapi:: ana_config
#   :include-all-objects:
#   :no-inheritance-diagram:

.. automodapi:: grb
       :no-inheritance-diagram:

.. automodapi:: grb_plot
       :no-inheritance-diagram:

#.. automodapi:: visibility
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: timeslot
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: obs_slice
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: mcsim
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: mcsim_config
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: mcsim_res
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: mcsim_plot
#   :include-all-objects:
#   :no-inheritance-diagram:



#.. automodapi:: irf_onaxis
#   :include-all-objects:
#   :no-inheritance-diagram:

#.. automodapi:: fit_onoff12
#   :include-all-objects:
#   :no-inheritance-diagram:
