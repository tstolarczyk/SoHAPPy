SoHAPPy with GammaPy,
Simulation of High-energy Astrophysics Processes in Python
roma version (gammapy 0.12)

Installation
============
SoHAPPy will easily run under an anaconda environement, requiring essentially
the gammapy package.
The Roma and Sofia versions are compatible with gammapy 0.12.

To create and activate the proper environment do :
	conda create -c conda-forge -n gpy0.12 gammapy=0.12
	conda activate gpy0.12
	conda install matplotlib

A created environment can be removed as follows :
conda env list
conda env remove -n gpy0.12

dependecy
pandas 1.03
...

Run the code
============
Type python SoHAPPy.py -h to get help on the command line
Any omitted argument is taken from the ana_config.py file.
The recommanded defautf values are listed below.

For a large number of events the debugging should set to zero.
In that case no plots are produced and the matplotlib is not required

Default recommanded parameters
==============================
(subject to change without notice)

import astropy.units as u
#-----------------------------------------------------------------------------#
# These values can be superseded by the command line
niter  = 100 # Number of iterations for the Monte Carlo
ngrb = 1 # Number of GRB to be read - if 1 specific studies on one grb

ifirst=30 # First GRB - can be a list then ngrb is useless
#ifirst = [11, 12]

# Debugging : if negative or zero : no plot
# 0: evt counting, 1: + some results, 2: + event details
# 3: Plots for each trials - create an animated gif - heavy
dbg    = 2
res_dir    = "D:/000_Today" # Folder for results

#-----------------------------------------------------------------------------#
# Physics
EBLmodel        = "dominguez"
#-----------------------------------------------------------------------------#
# Detection
dtslew    = 30*u.s # Time delay to point the GRB, usually 30*u.s
fixslew   = True   # If False, random < dt_slew
dtswift   = 77*u.s # Mean SWIFT latency : mean value is 77*u.s
fixswift  = True   # If False, the latency is generated from real data file
swiftfile = "../input/swift/Swift_delay_times.txt" # Real data file
#-----------------------------------------------------------------------------#
# Analysis
lkhd      = 0     # If zero : on-off analysis; otherwise degrees of freedom
obs_point = "end" # Define the position of observation in the time slice
#-----------------------------------------------------------------------------#
# Development
save_simu      = False # (False) Simulation svaed to file for offline use
save_grb       = False # (False) GRB saved to sik -> use grb.py main
write_slices   = False # (False) Store detailed information on slices if True
signal_to_zero = False # (False) Keep only background, set signal to zero
do_fluctuate   = True  # (True) Statistical fluctuations are enabled
do_accelerate  = True  # (True) Simulation skipped if no detection in first 10%
fixed_zenith   = False # If set to a value ("20deg") freeze zenith in IRF
remove_tarred  = True  # Remove tarred files, otherwise keep for faster access
test_prompt    = False # (False) If True test prompt alone (experimental)
use_afterglow  = True  # Prompt characteritics from the afterglow with same id.
day_after      = 0     # Add up to that days (if negative only that day)
#-----------------------------------------------------------------------------#
# Input-Ouput
old_file   = False # Use old GRB file
grb_olddir = '../input/lightcurves/long_grb_test/'
grb_dir    = '../input/lightcurves/LONG_FITS'
irf_dir    = "../input/irf/OnAxis/"



History
=======
Since Feb. 2020

SOFIA version
-------------
The So-called SOFIA version, named after the cancelled CTA meeing in Sofia (May 2020)
intends to add to the Roma version the functionnalities requested to stack observations
from several site or several days. It has also to take in charge the Gammapy 0.16
version and the use of a likelihood analysis.
Developement started on March 24, 2020.

2020 09 07 : The IRF valifity windows in zenith have been changed, following
             the recommandations. In particular no IRF are made avaialble above
             66° (24° of altitude). The time validity of the IRF have also
             be changed to reflect the logarithmic binning of the IRF sampling
             in observation time.

2020 09 04 : The multi-window visibility computation is completed. The
             day_after (skip the first day and compute detection adding 1 day
             to all dates) has become useless.
             It is some how replaced bythe varaible depth that fix the number
             of days after the trigger the visibility windows are searched for,
             and the variable skip that defines the number of first nights to
             be skipped (skip=0 to consider first night, skip=1, to consider
             from first night on).

2020 08 27 : Following the discrepancies found on visibility, MGB has rerun
             his code for the GRB that have a very short period seen at the end
             of the day. The repository has been updated on :
             https://drive.google.com/drive/folders/1sFYlG0ZFa2GpDOh-MuM_a5c2_Sn-6Jes
             and the agreement is obtained except for 197 N for whihc a 4
             second long time window is not seen by MGB.

2020 08 26 : The GammaRayBurst t_start and t_stop dictionnaries now become a
             list in order to handle lore than one visibility period.

2020 08 25 : Create a "docs" file to generate a SPHINX documentation
             This folder contains :
             - conf.py : general configuration
             - index.rst: Is the seed for the forthcoming index.html entry
             page. It lists the modules to be added to the documentation
             - make.bat, the main engine for the documentation creation
             - Makefile
             - README.rst

             The documentation is generated as follows from the docs
             folder:
             cmd.exe /c .\make.bat html

             In order to handle the numpydoc style of function documentation,
             numpydc has to be added to the extension list in conf.py

2020 08 25 : BUG in Astroplan
             An issue has been opened on the Astroplan github to report a
             misbehaviour in the search of the rise-time.
             For three events, 233 S, 644 N and 769 S, the correct rise
             period is missed because it s very close to th etrigger time
             and falls into a first time slice in Astroplan whihch is
             systematiccally and wrongly discarded.
             See : https://github.com/astropy/astroplan/issues/464

2020 08 25 : Further studies investigate the discrepancies between the default
             visibilities provided by Maria Grazia. The situation is explained
             in the text below, adapted from a from a mail sent to MGB.

             *** Short windows missed by MGB

             Using the same site coordinate as LMGB (see earlier comments), and
             a search limited to 24h (end of visibility)), some missed short
             windows of are found that make some GRB becomint Observable
             tonight:

             North (14): 129, 174, 197, 214, 309, 380, 381, 479, 488, 502, 609,
                         673, 933, 965
             South (14): 24, 30, 182, 186, 319, 338, 444, 475, 734, 772, 826,
                         879, 917, 986

             *** The search of one day after the trigger should be applied to
             the start date

             The event 522 N was not considered visible in SoHAPPy because the
             depth (duration afte rthe trigger) to accept a window was on the
             end. It has been changed to the beginning (i.e. a visibility
             window can extend much after the "depth", typically one day)

            *** Some events have 2 visibility windows in one night
            There are 19 events in which the GRB is above the horizon at dusk,
            goes below during the night but rise again before dawn. They were
            not expected by MGB that looked for the first window only.
            The list is the following:
            North : 75, 137, 139, 151, 153, 239, 270, 353, 507, 727, 859, 906,
                    923, 948
            South : 85*, 115, 414, 754, 755
            Event 85 has even 3 visibility windows !

-
2020 07 30 : Add scale="utc" to all Time definitions.
2020 07 28 : Using the dedicate function check_visibility in grb.py, one check
        the visibilities of the GRB for the delault altitude of 10°.
        The defauly values were obtained by maria Grazia Bernardini (MGB).
        A few discrepancies were found an understood (MGB, July 2020).

        The observatory coordinates
        ---------------------------
        In SoHAPPy they are obtained from EarthLocation.of_site() at the sites
        'Roque de los Muchachos' and  'Paranal Observatory'. In order to
        avoid accessing the external database, the coordinates have been
        copied and harcoded as follows.

        SoHAPPy uses, (x,y,z) from geocentric :
        North : ( 5327448.9957829, -1718665.73869569, 3051566.90295403) m
        South : ( 1946404.34103884,-5467644.29079852,-2642728.20144425) m

        MGB uses (lon,lat,elevation) from geodetic :
        North : ('342.1184', '28.7606', 2326. * u.meter)
        South : (289.5972',  '-24.6253', 2635. * u.meter)
        or after applying EarthLocation_from_geodetic :
        North :  (5327285.09211954, -1718777.11250295, 3051786.7327476) m
        South : (1946635.7979987, -5467633.94561753, -2642498.5212285) m

        The distance discrepancy between the positions is 296 and 326 m.

        This creates variation in the rise/set time of the target which is
        usually at the level of a few seconds. Note that the discrepancy can be
        bigger of the numbr of grid step in Astroplan is not at default (150).
        When  using MGB values, I still find discrepancies at the
        level of second or less for a very few GRBs.


2020 07 02 : - Astroplan 0.7 dev version is used because of a bug in the
               released astroplan v0.6 (see astroplan issue #460 on Astroplan
               github).
               Installation was done from the github repository, see here :
               https://astroplan.readthedocs.io/en/latest/installation.html#id2

2020 07 01 : - New function visibility that compute the start and stop of the
               visibility taking into account any horizonand any number of days
               after the trigger. In particular it complies wuth the native
               visibility values given by Maria Grazia Bernardini that are stored
               in the fits file.
               New variables t_start and t_stop are added to the GammaRayBusrt
               class and store the start and end point of the vsisbility in
               seconds (time relative to the GRB trigger)
2020 06 18 : - the function creating the time slices are moved from SoHAPPy.py to
               timeslot.py and become slot methods.
2020 06 12 : - new flag "silent" to avoid any print on screen (output
               redirected ot the log file)
2020 06 11 : - The output data, configuration and log files have now a fixed
               name, and the results are tarred in a single file identified
               from the output folder name.
2020 06 09 : - MAJOR bug found. In the core of the program, when the signal
               count is computed, the ith spectrum is used. The slice concerned
               is the ith slice of the selected time window. The bug was that i
               does not refer to this but to the ith spectrum in the list,
               i.e. whathever the time slice organisation is, the first slice
               is always the first grb slice, the one with the highest flux.
               This had obviously dramatic consequences on the counting rate.
2020 05 26 : - Add a function that restruit the visibilty to the GRB avalaible
               data points when masking the slices.
2020 05 25 : - Prompt file reading and corresponding flag are added for tests.
               The GRB information is obtained from the afterglow files.
2020 05 21 : - EBL model in ana_config is passed to the GRB class (instead of
               the Absorption object)
2020 05 18 : - The copy of the config file is moved from SoHAPPy.py to
               utilities.py and is more general.
             - A log file is now created that replace the screen output to
               a text file. A class, Log, is created in utilities.py to
               dump messaged on both the screen and the log file
             - A log file is now created, sensible message are sent to it.
2020 05 04 : - Add back the possibility to simulate the ith day after the first
               one.
2020 04 27 : - After re-implementaion of the chain to have the possibility to
               merge sites and subsequent days, a slight difference in results
               was found(initially on GRB 188) between the present code and
               the Roma version. This has been believed to be due to the fact
               that the slewing delay was strictly apply to the GRB trigger
               time (i.e. GRB observation starts at t_trig + 30 s in case of
               a fixed slew duration), whereas in the Roma version the delay
               was applied from he first measurement point.
               (===> But this was not the reason ! See the bug corrected in the
               Roma version on April 27th in the corresponding Readme file).
               In the cases studied, the first observation slice was
               31.02-37.20 s for the Roma version and 30.00-37.20 for the Sofia
               version. The difference in start time is exacltly the time of
               the first measurement point : (1.02e+00, folowed by 1.27e+00,
               1.60e+00).
               This difference leads to several comments :

               1. The flux at zero is not known. The first measurement is here
               a bit more than 1 s.
               This was commented by Lara Nava in a mail on December 18th 2019:

               "In GRB afterglow modeling, t=0 is the time a ‘’virtual'' photon
                emitted at R=0 reaches the observer. So no photons can arrive
                at t=0. Arrival times can be 10^ minus something. However, to
                start simulations from very small times, such as 0.0001 or so
                increases a lot the computation time. Since there is a delay
                due to the notice arrival + repointing time, I guess that to
                start simulations at 0.01, 0.1 or 1-2 seconds is irrelevant for
                 the simulations of detectability with CTA. Am I correct?
                In any case, the starting time of simulations is calculated by
                the simulation for each GRB, accounting for its specific
                parameters, by requiring that simulations start well before the
                 deceleration time, when the flux is small (and then
                 negligible)."

               As a consequence, in the case of a null slewing time, the
               minimal delay cannot be less than the time at which the first
               measurement point is know.

               2. Even in the case of a GRB immediately in the field-of-view,
               which can be only very rare, there is a minimal time between
               the GRB aws detected by SWIFT and the alert was received on
               ground. The time distribution of these delays was received from
               Maria Grazia Bernardini on March 2nd, 2020. The mean Swift
               delay is 77 seconds.

               The Sofia version takes into account these delays and
               limitation.

2020 04 09 : - Now that stacked site observations are possible, the simulation
               is no more related to a unique altiute and azimuth start and
               stop. These should be recomputed at the analysis time from
               the tstart and tstop variables, although this can be problematic
               when several days are stacked, in particular when it exceeds the
               maximulm GRB simulation time (tstop would then correspond to
               the lighcurve end).

2020 04 08 : - Because iof dates more than 30 years in the future, astropy
               generate Erfa warnings (unable to ensure precisionon leap
               seconds on time beyond 30 years). For that reason, warning are
               remporarily disabled where needed :
                       with warnings.catch_warnings():
                            warnings.filterwarnings("ignore")
                            'code generating the warning'
2020 04 03 : - Modifies grb.py and grbplot.py so that is odes not refer anymore
               to measurement points, that are an analysis decision; The GRB
               GammaRayBurst class now contains monly information obtained from
               the input files.

2020 03 27 : - change the position of the Astropy cache to avoid warning when
               the cache data are written in a protected directory.
               The error is reported below :

[...] astropy\config\configuration.py:557: ConfigurationMissingWarning:
Configuration defaults will be used due to FileExistsError:17 on None  warn(ConfigurationMissingWarning(msg))
WARNING: CacheMissingWarning: Remote data cache could not be accessed due to FileExistsError: [WinError 183]
WARNING: CacheMissingWarning: Not clearing data cache - cache inaccessible due to FileExistsError: [WinError 183]
WARNING: CacheMissingWarning: ('Cache directory cannot be read or created,
providing data in temporary file instead.',
'C:\\Users\\stolar\\AppData\\Local\\Temp\\astropy-download-6212-out__lab') [astropy.utils.data]

2020 03 24 : - Start to implement the posibility to havethe code running with
               bioth gammapy 0.12 and 0.16 using a test on gammapy.__version__
ROMA version
------------

This version has for basis the code used to produce the results at the Roma face-to-face meeting.
It got a few modifications as below. This is essentially a code missing the essential features:
- stacking North and South simulations, second and later days contributions,
and possibly a second spectrum due to the prompt.
- Using Gammapy 0.16 and in particular Full enclosure off-axis IRF, requiring
a major change in the fitting process.

History of changes between Feb. 2020 and 2020 03 23

2020 03 18 : - The 'analysis' folder, that contain a series of notebook to
               analyse the data, is removed from the SoHAPPy folder.
               It is moved one level up, next to the 'input' and 'output'
               folders.  This requires a few changes in the path of certain notebooks.
               It is intended that the analysis folder is not part of the
               SoHAPPy distribution anymore (in particular its development
               speed has nothing to do with the SoHAPPy ones))
             - Minor changes and clean-up made to remove the dependency to matplotlib
			   when the run does not require it (debug flag less or equal to zero). This
			   is required to run SoHAPPy with non graphical environements like for
			   GRID productions (where matplotlib is not installed in principle).
2020 03 14 : - when running SoHappy from a shell script and requesting plot to be
			   shown, the script does not pause anymore at each plot
			   (use of plt.show(block=False).
2020 04 22 : - Class ObservationParameter removed. Most of the information is
               in the MonteCarlo class. The alpha parameter becomes an
               internal parameter of fit_onoff and cannot be changed
             - Protection agains GERB for whihc all sigma would be Nan
             - On-Off siginificance are Nan if non or noff below limit (e.g. 10)
2020 02 04 : - Changed cta_irf_onaxis in irf_onaxis
2020 02 03 : - MC specific parameters now in mcsim_config.
             - Changes made on variables
             - Parameters removed from the MonteCarlo class as they are not
             MC dependent
             - Add a global boolean variable to fit_onoff to force S=0
             - Add config files to mcsim and grb to handle global variables and
             options
             - Improve grb display with fewer selected spectra and LC displayed


