Configuration file
==================
A simulation and analysis parameter file can be passed to the main script with the `-c` options (some of these options can be superseded on the command line, use  `-h` to get the list of parameters).
If not file is given, the parameters are read from the *config.yaml* in the `SoHAPPy` repository.
The excecution of `SoHAPPy` automatically create a copy (backup) of the chosen confugration file to the output folder. This copy can be used to reproduce the results later. It is also necessary to run the provided Jupyter notebooks (`analysis sub folder).
Modifications of SoHAPPy may lead to new parameters that have to be added "by hand" in the former backup configura`tion, either to run a simulation/analysis or use the Jupyter Notebooks.
  
The following parameters are defined and mandatory in the configuration file.

Input and output
----------------

.. tabularcolumns:: |l|c|p{5cm}|

+------------------+-------------------+--------------------------------------------------+
| variable         | Default/Suggested | What is it ?                                     |
+==================+===================+==================================================+
| ifirst           | 1                 | the first GRB file to be processed               |
+------------------+-------------------+--------------------------------------------------+
| ngrb             | None              | | Number of GRB to be processed                  |
|                  |                   | | If 1, special actions are taken                |
+------------------+-------------------+--------------------------------------------------+
| trigger          | 0                 | | Time shift in days applied to the original     |
|                  |                   | | trigger dates.                                 |
|                  |                   | | If a float, is added to the existing dates.    |
|                  |                   | | Can be any duration and/or a number of Julian  |
|                  |                   | | days to have absolute dates from relative ones.|
|                  |                   | | If a file, read and overwrite the dates from   |
|                  |                   | | the file, starting from line 3 (First and      |
|                  |                   | | second lines are reserved). The number of dates|
|                  |                   | | should match the number of sources (ngrb).     |
+------------------+-------------------+--------------------------------------------------+
| in_folder        | ../input          | | Input main folder                              |
|                  |                   | | Contains the `visibility`, `irf`, Swift        |
|                  |                   | | latency data                                   |
+------------------+-------------------+--------------------------------------------------+
| res_dir          | ../output         | Output folder                                    |
+------------------+-------------------+--------------------------------------------------+
| data_dir         | None              | Data (GRB files) subfolder                       |
+------------------+-------------------+--------------------------------------------------+
| irf_dir          | None              | IRF subfolder (specific organisation)            |
+------------------+-------------------+--------------------------------------------------+

(1): e.g. `Time('1800-01-01T00:00:00', format='isot', scale='utc').jd = 2378496.5` 
    
Simulation parameters
---------------------

.. tabularcolumns:: |l|c|p{5cm}|


+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| niter                 | 1                      | Number of Monte Carlo trials                |
+-----------------------+------------------------+---------------------------------------------+
| seed                  | 2021                   | Choose 'random-seed' to randomize           |
+-----------------------+------------------------+---------------------------------------------+
| do_fluctuate          | False                  | Statistical fluctuations are enabled        |
+-----------------------+------------------------+---------------------------------------------+
| do_accelerate         | True                   | | When True, the simulation is stopped if   |
|                       |                        | | none of the first trials in the limit     |
|                       |                        | | of 1 - det_level have reached the minimal |
|                       |                        | | significance. In a simulation with 100    |
|                       |                        | | trials requesting a 3 sigma detection in  |
|                       |                        | | 90% of the iterations, the simulation will|
|                       |                        | | be stopped after 10 iterations*           |
+-----------------------+------------------------+---------------------------------------------+


Simulation physics
------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| EBLmodel              | "dominguez"            | | Extragalactic Background Light model as   |
|                       |                        | | defined in gammapy or 'built-in'          |
+-----------------------+------------------------+---------------------------------------------+
| prompt                | False                  | | If True read prompt model from disk       |
|                       |                        | | (see doc. for accepted data)              |
+-----------------------+------------------------+---------------------------------------------+
| array_North           | "FullArray"            | IRF North array                             |
+-----------------------+------------------------+---------------------------------------------+
| array_South           | "FullArray"            | IRF South array                             |
+-----------------------+------------------------+---------------------------------------------+
| dtslew_North          | "30 s"                 | Maximum slewing time delay in North         |
+-----------------------+------------------------+---------------------------------------------+
| dtslew_South          | "30 s"                 | Maximum slewing time delay in South         |
+-----------------------+------------------------+---------------------------------------------+
| fixslew               | True                   | If False generate a random delay < dtslew   |
+-----------------------+------------------------+---------------------------------------------+
| dtswift               | 77*u.s                 | | Fixed SWIFT latency**:                    |
|                       |                        | | - minimum : 12 s;                         |
|                       |                        | | - mean value : 77 s;                      |
|                       |                        | | - median : 34 s;                          |
|                       |                        | | - 65% of GRBs have delays < 52 s;         |
|                       |                        | | - 90% of GRBs have delays < 130 s.        |
+-----------------------+------------------------+---------------------------------------------+
| fixswift              | True                   | | If False, the latency is generated from   |
|                       |                        | | real data file (see below)                |
+-----------------------+------------------------+---------------------------------------------+
| swiftfile             | | swift/               | Swift latency file                          |
|                       | | Swift_delay_times.txt|                                             |
+-----------------------+------------------------+---------------------------------------------+
| visibility            | strictmoonveto         | | Can be "Null" (read from the data files if|
|                       |                        | | it exists), a subfolder where to find the |
|                       |                        | | pre-computed visibility files, a keyword  |
|                       |                        | | corresponding to a dictionnary entry in   |
|                       |                        | | `visibility.yaml` to compute the          | 
|                       |                        | | visibility on the fly.                    |
+-----------------------+------------------------+---------------------------------------------+

Debugging and bookkeeping
-------------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| dbg                   | 0                      | From 0 to 3, increasingly verbosy           |
+-----------------------+------------------------+---------------------------------------------+
| save_simu             | False                  | Simulation saved to file for offline use    |
+-----------------------+------------------------+---------------------------------------------+
| save_grb              | False                  | GRB class saved to disk -> use grb.py main  |
+-----------------------+------------------------+---------------------------------------------+
| datafile              | "data.txt"             | Population study main output file           |
+-----------------------+------------------------+---------------------------------------------+
| logfile               | "analysis.log"         | Text file with results, status and warning  |
+-----------------------+------------------------+---------------------------------------------+
| remove_tar            | False                  | | Remove tarred files, otherwise keep for   |
|                       |                        | | faster access                             |
+-----------------------+------------------------+---------------------------------------------+


Experts and developpers only
----------------------------

.. tabularcolumns:: |l|c|p{5cm}|


+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| method                | 0                      | Not used (detection method)                 |
+-----------------------+------------------------+---------------------------------------------+
| obs_point             | "end"                  | Observation position in the time slice      |
+-----------------------+------------------------+---------------------------------------------+
| test_prompt           | False                  | If True test prompt alone (experimental)    |
+-----------------------+------------------------+---------------------------------------------+
| use_afterglow         | False                  | | Prompt characteritics from the afterglow  |
|                       |                        | | with same id.                             |
+-----------------------+------------------------+---------------------------------------------+
| signal_to_zero        | False                  | Keep only background, set signal to zero    |
+-----------------------+------------------------+---------------------------------------------+
| fixed_zenith          | None                   | If a value ("20*u.deg") freeze zenith in IRF|
+-----------------------+------------------------+---------------------------------------------+
| magnify               | 1                      | Multiplicative factor of the flux, for tests|
+-----------------------+------------------------+---------------------------------------------+
| silent                | False                  | If True, nothing on screen (output to log)  |
+-----------------------+------------------------+---------------------------------------------+
| write_slices          | False                  | Store detailed information on slices if True|
+-----------------------+------------------------+---------------------------------------------+
| save_dataset          | False                  | Not implemented (save datasets)             |
+-----------------------+------------------------+---------------------------------------------+
| forced_visible        | False                  | | If True, the GRB is always visible        |
|                       |                        | | (infinite nights)                         |
+-----------------------+------------------------+---------------------------------------------+
| n_night               | Null                   |  Limit data to a maximal number of nights   |
+-----------------------+------------------------+---------------------------------------------+
| Emax                  | Null                   |  Limit data energy bins to Emax             |
+-----------------------+------------------------+---------------------------------------------+


(*) Note that this bias the resulting population since it articiially deplete the max significance population below the minimum required (e.g. 3 sigma).

(**) M. Grazia Bernardini, private communication, February 28th, 2020