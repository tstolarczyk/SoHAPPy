SoHAPPy with GammaPy, 
Simulation of High-energy Astrophysics Processes in Python
roma version (gammapy 0.12)

Installation:
SoHAPPy will easy run under an anaconda environement, requiring essentially the gammapy package.
The roma version is compatible with gammapy 0.12.

To create and activate the proper environment do :
	conda create -c conda-forge -n gpy0.12 gammapy=0.12
	conda activate gpy0.12
	conda install matplotlib


A created environment can be removed as follows :
conda env list
conda env remove -n gpy0.12

-----
SoHAPPy History
Since Feb. 2020
2020 04 27 : - BUG CORRECTED 
               The observation time and the alt-az were wrongly considered to 
               be the beginning of the time slice, which is contrary to the 
               initial assumption.
               This was corrected in mcsim_res.py in mean_velues
               The bug was introduced when the function mean_values was 
               created, i.e. after the Roma face to face meeting (i.e.
               the master version in GitHub is correct and has been tagged 
               v3.1 for further reference).
               Note that by default,the best performance is found with the
               time slice duration at the alt-az at the beginning of the
               slice which is somehow incoherent with the above asumption.

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
             
             
    
To do :
    - Extrapolate Fermi GRBs to CTA, compute response
	- Inject Fermi data and IRF, compute response
	- min zenith should be a chosen parameter
    - reject slices for which zenith are below the min slices
    - Add number of counts to the output;
    - Compute the fluence;
    - Interpolate time counts within the time slices;
    - in output file, separate constant information from variable information (use a header).
    - DONE : ensure "show" for plot showing depends on debug like on the command line, 
    remove 'show", from ana_config'
    - DONE : not necessary as SoHappy depends on gammapy "only" : create an environment to be able to reinstall the code to another 
    platform (conda env export)
    - use logging for warning and debug messages
    - create test units and run them (e.g. a dummy simplistic GRB)
"""


