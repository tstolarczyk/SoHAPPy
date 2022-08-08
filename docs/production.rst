Productions
===========
Recommendations for running large simulations and analysis.

Introduction
------------

A prodcution starts from a set of files with a collection of energy spectra for a series of time slices, plus some physical information from the sources (redshift, trigger time, sky position). The data re usually stored as fits file with a specific format but other formats are possible (like yaml file to reproduce simplistic lighcurve with two time and energy indices.

In some cases, the GRB trigger (explosion) time should be changed, either for specific studies or because the 
population should be generated between two given dates. In that case, a yaml file with trigger dates should be produced using the script `trigger_dates.py`.

At this stage the vsibilities will have changed. If they are not computed on the fly, i.e. read from the input grb files or from files read of disk, then they should be recomputed taking into account the new trigger times. This can be done with the visibility_write script main function.
Computing the visibility of 1000 GRB takes 2 hours (7.5 second for each source, 3.2 seconds per site).

A full production with 100 iteration per all trials takes 2 hours for 1000 GRBs. 

Configuration
-------------

At this stage, the configuration file can be edited. It  will pass all simulation and analysis parameters to SoHAPPy.
An example configration file, config.yaml, is available in the code repository and can be copied and modified for that purpose.


* INPUT/OUPUT
  - choose the folders:
		- input folder, where the data (data_dir) and irf (irf_dir) subfolder should be found.
			- the irf folder depends on the detector simulation used (e.g. prod3 or prod5)
  - set the first file to be read and the number of files; Files are supposed to be named Eventxyz.fits, with xyz from 0 or 1 depending on the data set
  - set the output folder (res_dir)
  - the trigger (or explosion) dates ca be changed either by adding a constant delay in days to the existing dates or by providing a text file with one number od days per source to supersed the existing date. The text file can be produced from the script explosion.py (see later for the usage). Note that the visibility depends on the explostion dates and they should be computed beforehand, or on the fly accordingly (in particular the dfault visibility cannot be used)
  
* SIMULATION PARAMETERS

Results can be obtained quickly with niter=1, but fluctuations will be considered only above 1. Expereice shows that 100 iterations give a decent convergence to the results
The seen should be left untouched in order to have results reproduced otherwise use "random-seed"
When the number f iteration is >1, do_fluctuate should be True
do_accelerate help saving time but prevent frmm comparing two simlations as the number of trials varyis (the simualtion stops as soon as 10% of the total number of iteratiosn ahve not given a 5 sigma detection)

* SIMULATION PHYSICS

