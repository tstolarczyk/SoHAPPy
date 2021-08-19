

Steering parameters
===================

Here are the steering parameters.



Events to be processed
----------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| ifirst                | 1                      | the first GRB file to be processed          |
+-----------------------+------------------------+---------------------------------------------+
| ngrb                  | None                   | | Number of GRB to be processed             |
|                       |                        | | If 1, special actions are taken           |
+-----------------------+------------------------+---------------------------------------------+

Simulation and analysis
-----------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| niter                 | 1                      | Number of Monte Carlo trials                |
+-----------------------+------------------------+---------------------------------------------+
| method                | 0                      | Not used (detection method)                 |
+-----------------------+------------------------+---------------------------------------------+
| obs_point             | "end"                  | Observation position in the time slice      |
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
| EBLmodel              | "dominguez"            | | Extragalactic Background Light model as   |
|                       |                        | | defined in gammapy or 'built-in'          |
+-----------------------+------------------------+---------------------------------------------+
| arrays                | {"North":"FullArray",  | Indicate the IRF subarrays in each site     |
|                       |  "South":"FullArray"}  |                                             |
+-----------------------+------------------------+---------------------------------------------+

visibility
----------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| depth                 | 3                      | Max. number of nights                       |
+-----------------------+------------------------+---------------------------------------------+
| skip                  | 0                      | Number of first nights to be skipped        |
+-----------------------+------------------------+---------------------------------------------+
| altmin                | 24*u.deg               | Minimum altitude                            |
+-----------------------+------------------------+---------------------------------------------+
| altmoon               | -0.25*u.deg (horizon)  | Moon maximum altitude                       |
+-----------------------+------------------------+---------------------------------------------+
| moondist              | 30*u.deg               | Moon minimal distance                       |
+-----------------------+------------------------+---------------------------------------------+
| moonlight             | 0.6                    | Moon maximum brigthness                     |
+-----------------------+------------------------+---------------------------------------------+
| dtslew                | {"North":30*u.s,       | Maximum slewing time delay                  |
|                       |  "South":30*u.s}       |                                             |
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
| swiftfile             | | ../input/swift/      | Swift latency file                          |
|                       | | Swift_delay_times.txt|                                             |
+-----------------------+------------------------+---------------------------------------------+
| vis_cmp               | True                   | If True, visibility is computed on the fly  |
+-----------------------+------------------------+---------------------------------------------+
| vis_dir               | None                   | If defined, visibility are read from folder |
+-----------------------+------------------------+---------------------------------------------+

Special simulation actions
--------------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
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

Debugging
---------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| dbg                   | 0                      | From 0 to 3, increasingly verbosy           |
+-----------------------+------------------------+---------------------------------------------+
| silent                | False                  | If True, nothing on screen (output to log)  |
+-----------------------+------------------------+---------------------------------------------+
| **INPUT**                                                                                    |
+-----------------------+------------------------+---------------------------------------------+
| grb_dir               | None                   | Folder hosting the GRB files                |
+-----------------------+------------------------+---------------------------------------------+
| irf_dir               | None                   | IRF folder (specific subfolder organisation)|
+-----------------------+------------------------+---------------------------------------------+

Output
------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| res_dir               | ../output              | The output folder                           |
+-----------------------+------------------------+---------------------------------------------+
| datafile              | "data.txt"             | Population study main output file           |
+-----------------------+------------------------+---------------------------------------------+
| logfile               | "analysis.log"         | Text file with results, status and warning  |
+-----------------------+------------------------+---------------------------------------------+

Special output actions
----------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| save_simu             | False                  | Simulation saved to file for offline use    |
+-----------------------+------------------------+---------------------------------------------+
| save_dataset          | False                  | Not implemented (save datasets)             |
+-----------------------+------------------------+---------------------------------------------+
| save_grb              | False                  | GRB class saved to disk -> use grb.py main  |
+-----------------------+------------------------+---------------------------------------------+
| write_slices          | False                  | Store detailed information on slices if True|
+-----------------------+------------------------+---------------------------------------------+
| remove_tarred         | False                  | | Remove tarred files, otherwise keep for   |
|                       |                        | | faster access                             |
+-----------------------+------------------------+---------------------------------------------+

(*) Note that this bias he resulting popualtion since it articiially deplete the max significance population below the minimum required (e.g. 3 sigma).

(**) M. Grazia Bernardini, private communication, February 28th, 2020