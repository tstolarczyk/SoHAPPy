# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:34:07 2021

@author: Stolar
"""

This module hosts the functionnalities to choose the best IRF from the available
files and the recipes on their validity domains in observation duration and
angles.

The recommendations for prod3v2 IRF are the following (Information from
G. Maier obtained by I. Sadeh on January 30th, 2020), and later private
discussion on March 31st, 2021.

Angular "interpolation"
=======================

Zenith
------
The zenith windows should be defined following a 1/cos(theta) law, which means
(Iftach says) that depending on the zenith angle the following IRF should be
used :
    - zen < 33        : "20deg"
    - 33 <= zen < 54  : "40deg"
    - zen >= 54        : "60deg"

On top of this, it is foreseen that the instrument performance will not be
reliable above 66 degrees as implicitly mentionned in A-PERF-0720
(*Observable Sky: The system as a whole must be able to target any astrophysical
object in the sky which has an elevation >24 degrees*).

Azimuth
-------
It is recommended to use the average azimuth (not a step function due to the North and South IRF).

Energy thresholds
=================
Although events should be generated from the minimum allowed value for a given
IRF, the analysis should restrict to a higher threshold, in order to avoid
cutoff effects.

Minimum *generated* energies are the following:
    - 20deg —> E_gen >= 12.6 GeV
    - 40deg —> E_gen >= 26.4 GeV
    - 60deg —> E_gen >= 105.3 GeV

Minimum *reconstructed* energies are the following (Iftach suggests these
values but does not say it this is his conclusion or more generally accepted
values):
    - 20deg —> E_grb >= 30 GeV
    - 40deg —> E_grb >= 40 GeV
    - 60deg —> E_grb >= 110 GeV
This corresponds approximately to the generated energies to which one standard
deviation has been added (the expected energy resolution at low energies is
deltaE/E (68%) = 0.25.

Note that the CTA performance plots have a cutoff of E>20 GeV.
Requirements expect the IRF to be reliable above 30 GeV.

IRF data
========
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
Dircet download (20181203): https://forge.in2p3.fr/attachments/download/62074

All IRF available for : 100s, 30m, 5h, 50h
North:
    FullArray, LST, MST, TS, TS_MST :
        20deg, 40 deg, 60 deg
        N, S, average
South:
    FullArray, LST, MST, MSTSST, SST, TS, TS_MST, TS_SST
    20deg, 40 deg, 60 deg
    N, S, average

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
North:
    FullArray, LST, MST : 20deg NSBx05, 20deg Reco2
    TS, TS_MST : 20 deg reco2

South:
    FullArray : 20deg NSBX05 (average only), 20 deg Reco2, 50h only
    LST, MST : 20deg NSBx05 (average only)
    MSTSST, TS, TS_MST, TS_SST
    SST : 20deg NSBx05 and NSBx30 (average only)


On-axis files
-------------
SoHAPPy has been originally build to use the so called On-Axis IRF. This files
are extracted from the more general "Full enclosure off-axis" files (see below).
They consist in the IRF selected on-axis (center of the field of view) and
an agular extension corresponding to the energy dependednt 68% PSF containment.
They have been produced for both CTA-Mars and Evnt-Display analysis. In this
latter it corresponds to the official CTA performace curves (that are obtained
with Event-Display that indeed use a 68% containment whereas the MARS analysis
has an optimised containment cut).
At the date to the SoHAPPy development the following configuration were
available:

(Descritpion of arrays requested)
FullaArray:
    CTA Southern Site (covered area ~4 km2):
        4 Large-Sized Telescopes
        25 Medium-Sized Telescopes
        70 Small-Sized Telescopes )
    CTA Northern Site: (covered area ~0.6 km2))
        4 Large-Sized Telescopes
        15 Medium-Sized Telescopes
LST : LST only
MST : MST only
TS : Threshold
TS MST
MST-SST


All IRF available for : 100s, 30m, 5h, 50h

Reco1:
    North:
        ALL ARE NEW VERSIONS : (20181203)

        FullArray LST, MST, TS, TS-MST
            20deg : average, N, S
                    old versions (20170627) exist, same options

            20deg_NSBx05 : average, N, S, old versions (20170627)

            40deg: average, N, S
                    old versions (20170627) exist, same options
            60deg: average, N, S

    South:
        ALL ARE OLD VERSIONS : (20170627)

        FullArray
            20deg: average (not mentionned), N, S
            40deg: average (not mentionned), N, S
            Not available : 60deg
        LST, MST-SST, SST, TS, TS-MST, TS-MSTSST
            20deg: average only (not mentionned)
            40deg: average only (not mentionned)
            60deg : not available

Reco2 :
    ALL ARE OLD VERSIONS : (20170627)
    FullArray, LST, MST, TS, TS-MST

    North:
            20deg: average, N, S
            40deg not available
            60deg not available
    South:
            20deg; average (not mentionned)
            40deg: average (not mentionned)

#---------------------------

.. automodapi:: irf
   :include-all-objects:
   :no-inheritance-diagram: