Instrument response functions
#############################

The instrument response fucntion class and methods are in irf.py.
This module hosts the functionnalities to choose the best IRF from the available
files and the recipes on their validity domains in observation duration and
angles.


IRF data file organisation
==========================
prod5 v1
--------
This correspnd to teh alpha configuration, with 4 LST and 9 MST in NOrth, and 14 MST and 37 SST in South (no LST). The instrument response function (IRF) files are expected to be found in a
dedicated folder, with files organised as follows (the structure is defined for the code to find the best IRF).  
..  code-block::

+ FullArray 
	+-- North
		+-- Prod5-North-20deg-AverageAz-4LSTs09MSTs.1800s-v0.1.fits.gz
		+-- Prod5-North-20deg-AverageAz-4LSTs09MSTs.18000s-v0.1.fits.gz
    	+-- ...

    +-- South
    	+-- Prod5-South-20deg-AverageAz-14MSTs37SSTs.1800s-v0.1.fits.gz
    	+-- Prod5-South-20deg-AverageAz-14MSTs37SSTs.18000s-v0.1.fits.gz
    	+-- ...

+ MST

+ SST


prod3 v2
--------
The instrument response function (IRF) files are expected to be found in a
dedicated folder, with files organised as follows (the structure is defined for the code to find the best IRF). The IRF data are in gzip files named irf_file.fits.gz. Here is an excerpt of the organisation :

..  code-block::

+ FullArray 
	+-- North
		+-- 20deg
			+-- North_z20_average_5h
				+-- irf_file.fits.gz

			+-- North_z20_average_30m
			+-- North_z20_average_50h
			+-- North_z20_average_100s

    	+-- 40deg
    	+-- 60deg
    	
    +-- South
    	+-- 20deg
    	+-- 40deg
    	+-- 60deg


IRF validity ranges
===================

The recommendations for prod3v2 IRF are the following (Information from
G. Maier obtained by I. Sadeh on January 30th, 2020), and later private
discussion on March 31st, 2021.

Angular "interpolation"
-----------------------

*Zenith*

The zenith windows should be defined following a 1/cos(theta) law, which means that depending on the zenith angle the following IRF should be
used :
- zen < 33        : "20deg"
- 33 <= zen < 54  : "40deg"
- zen >= 54       : "60deg"

On top of this, it is foreseen that the instrument performance will not be
reliable above 66 degrees as implicitly mentionned in A-PERF-0720
(*Observable Sky: The system as a whole must be able to target any astrophysical
object in the sky which has an elevation >24 degrees*).

*Azimuth*

It is recommended to use the average azimuth (not a step function due to the North and South IRF).

Energy thresholds
-----------------

Although events should be generated from the minimum allowed value for a given
IRF, the analysis should restrict to a higher threshold, in order to avoid
cutoff effects.

Minimum *generated* energies are the following:
- 20deg —> E_gen >= 12.6 GeV
- 40deg —> E_gen >= 26.4 GeV
- 60deg —> E_gen >= 105.3 GeV

Minimum *reconstructed* energies are the following (suggested values):
- 20deg —> E_grb >= 30 GeV
- 40deg —> E_grb >= 40 GeV
- 60deg —> E_grb >= 110 GeV

This corresponds approximately to the generated energies to which one standard
deviation has been added (the expected energy resolution at low energies is
deltaE/E (68%) = 0.25.

Note that the CTA performance plots have a cutoff of E>20 GeV.
Requirements expect the IRF to be reliable above 30 GeV.

IRF available data
==================

Official public files
---------------------

The officla IRF, for the full array, can be found from here  :
    prod3b-v2 :
        https://www.cta-observatory.org/science/cta-performance/
        They have IRF for 20, 40 and 60 degrees, average azimuth

    prod3b-v1, computed in 2017 :
        https://www.cta-observatory.org/science/cta-performance/cta-prod3b-v1/
        They have only IRF for 20 deg and 40 deg, average azimuth

CTA consortium files
--------------------

Found and described from there :
https://forge.in2p3.fr/projects/cta_analysis-and-simulations/wiki/Prod3b_based_instrument_response_functions

Direct download (20181203): https://forge.in2p3.fr/attachments/download/62074

All IRF available for : 100s, 30m, 5h, 50h

..  code-block::

	- North:
		- FullArray, LST, MST, TS, TS_MST :
			- 20deg, 40 deg, 60 deg
			- N, S, average

	- South:
		- FullArray, LST, MST, MSTSST, SST, TS, TS_MST, TS_SST
			-  20deg, 40 deg, 60 deg
			- N, S, average
        

CTA consortium older files
--------------------------

IRF older than prodv3b2 can be used. They essentially differ in the angular
resolution at high energies. They include both full arrays, threshold arrays
and LST-only IRFs. They cover 20 and 40 degrees but not 60 degrees.
Some 20 degree IRF have 5x NSB (North and South) and for SST only 30x the nominal NSB.

See here :
https://forge.in2p3.fr/projects/cta_analysis-and-simulations/wiki/Prod3b_based_instrument_response_functions
And download there : https://forge.in2p3.fr/attachments/download/55554

ALL have 20deg and 40deg IRF, with azimuth N,S,average (except Reco2, average only)

- North:
    - FullArray, LST, MST : 20deg NSBx05, 20deg Reco2
    - TS, TS_MST : 20 deg reco2

- South:
    - FullArray : 20deg NSBX05 (average only), 20 deg Reco2, 50h only
    - LST, MST : 20deg NSBx05 (average only)
    - MSTSST, TS, TS_MST, TS_SST
    - SST : 20deg NSBx05 and NSBx30 (average only)




.. automodapi:: irf
   :include-all-objects:
   :no-inheritance-diagram: