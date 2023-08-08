Configuration file
==================
A simulation and analysis parameter file can be passed to the main script with
the `-c` options (some of these parameters can be superseded on the command
line, use  `-h` to get the list). If no file name is given, the parameters
are read from the *config.yaml* in the `SoHAPPy` repository.
The execution of `SoHAPPy` automatically creates a copy (backup) of the chosen
configuration file to the output folder, with the initial parameters
overwritten by the ones on the command line, so that this backuped
configuration can be used to reprocess the data in the exact same conditions.

The configuration data are stored and processed in the
:class:`configuration.Configuration` in :obj:`configuration.py`.

The following parameters are defined and mandatory in the configuration file.
Extra parameters can be added and are notified as ignored at run time.

Physics parameters
------------------

.. tabularcolumns:: |l|c|p{5cm}|

+------------------+-------------------+--------------------------------------------------+
| variable         | Default/Suggested | What is it ?                                     |
+==================+===================+==================================================+
| ``ifirst``       | 1                 | | The first GRB identifier to be processed       |
|                  |                   | | or a list of identifiers. The identifier can   |
|                  |                   | | include a string referring to a `yaml` file    |
|                  |                   | | with parameters to generate spectra from an    |
|                  |                   | | analytical function.                           |
|                  |                   | | (See the `data/historical` folder)             |
+------------------+-------------------+--------------------------------------------------+
| ``nsrc``         | None              | | Number of GRB to be processed if ``ifirst`` is |
|                  |                   | | an integer. If 1, special actions are taken.   |
|                  |                   | | Not used if ``ifirst`` is a list.              |
+------------------+-------------------+--------------------------------------------------+
| ``visibility``   | None              | | Can be:                                        |
|                  |                   | | * `built-in` (read from the data file if it    |
|                  |                   | | exists),                                       |
|                  |                   | | * a keyword pointing to a folder containing    |
|                  |                   | | `json` files with a colection of visibilities; |
|                  |                   | | * a keyword corresponding to a dictionnay entry|
|                  |                   | | in `visibility.yaml` to compute the visibility |
|                  |                   | | on the fly                                     |
|                  |                   | | * `forced` to force a single infinite night    |
|                  |                   | | * `permanent` to force a permanent visibility  |
+------------------+-------------------+--------------------------------------------------+
| ``EBLmodel``     | "dominguez"       | | Extragalactic Background Light model as        |
|                  |                   | | defined in gammapy or `built-in` if provided   |
|                  |                   | | or `None` if ignored (no absortpion)           |
+------------------+-------------------+--------------------------------------------------+

Input / Output
--------------

.. tabularcolumns:: |l|c|p{5cm}|

+------------------+-------------------+--------------------------------------------------+
| variable         | Default           | What is it ?                                     |
+==================+===================+==================================================+
| ``infolder``     | ../input          | Base input folder                                |
+------------------+-------------------+--------------------------------------------------+
| ``resfolder``    | ../output         | Base output main folder                          |
+------------------+-------------------+--------------------------------------------------+
| ``out_dir``      | None              | | Specific output subfolder, with a name         |
|                  |                   | | related to tests or IRF                        |
+------------------+-------------------+--------------------------------------------------+
| ``data_dir``     | None              | Data subfolder                                   |
+------------------+-------------------+--------------------------------------------------+
| ``prompt_dir``   | None              | | if not Null, read time-integrated              |
|                  |                   | | prompt from this subfolder                     |
+------------------+-------------------+--------------------------------------------------+
| ``irf_dir``      | None              | | IRF subfolder - should contain a tag           |
|                  |                   | | for the production used                        |
+------------------+-------------------+--------------------------------------------------+
| ``prefix``       | None              | | The prefix before the source identifier        |
|                  |                   | | to create the full file name                   |
+------------------+-------------------+--------------------------------------------------+
| ``suffix``       | None              | | The suffix after the source identifier to      |
|                  |                   | | create the file name (including the extension) |
+------------------+-------------------+--------------------------------------------------+

Simulation parameters
---------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+----------------+-------------------------------------------------+
| variable              | Default        | What is it ?                                    |
+=======================+================+=================================================+
| ``niter``             | 1              | Number of Monte Carlo trials                    |
+-----------------------+----------------+-------------------------------------------------+
| ``seed``              | 2021           | Choose 'random-seed' to randomize               |
+-----------------------+----------------+-------------------------------------------------+
| ``do_fluctuate``      | False          | Statistical fluctuations are enabled            |
+-----------------------+----------------+-------------------------------------------------+
| ``do_accelerate``     | True           | | When True, the simulation is stopped if       |
|                       |                | | none of the first trials in the limit         |
|                       |                | | of ``1 - det_level`` have reached the minimal |
|                       |                | | significance. In a simulation with 100        |
|                       |                | | trials requesting a 3 sigma detection in      |
|                       |                | | 90% of the iterations, the simulation will    |
|                       |                | | be stopped after 10 iterations*               |
+-----------------------+----------------+-------------------------------------------------+

(*) Note that this bias the resulting population since it artificially depletes the maximal significance population below the minimum required (e.g. 3 sigma).

Detection parameters
--------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| ``array_North``       | "FullArray"            | IRF North array                             |
+-----------------------+------------------------+---------------------------------------------+
| ``array_South``       | "FullArray"            | IRF South array                             |
+-----------------------+------------------------+---------------------------------------------+
| ``dtslew_North``      | "30 s"                 | Maximum slewing time delay in North         |
+-----------------------+------------------------+---------------------------------------------+
| ``dtslew_South``      | "30 s"                 | Maximum slewing time delay in South         |
+-----------------------+------------------------+---------------------------------------------+
| ``fixslew``           | True                   | If False generate a random delay less       |
|                       |                        |   than ``dtslew``. If True use ``dtswift``  |
+-----------------------+------------------------+---------------------------------------------+
| ``dtswift``           | 77*u.s                 | | Alert latency Fixed SWIFT latency         |
|                       |                        | | (e.g. SWIFT latency, with                 |
|                       |                        | | a mean value of 77 s)                     |
+-----------------------+------------------------+---------------------------------------------+
| ``fixswift``          | True                   | | If False, the latency is generated from   |
|                       |                        | | real data file (see below)                |
+-----------------------+------------------------+---------------------------------------------+
| ``swiftfile``         | | data/swift/          | Swift latency file                          |
|                       | | Swift_delay_times.txt|                                             |
+-----------------------+------------------------+---------------------------------------------+

Debugging and bookkeeping
-------------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| ``dbg``               | 0                      | From 0 to 3, increasingly verbosy           |
+-----------------------+------------------------+---------------------------------------------+
| ``silent``            | False                  | If True, nothing on screen (output to log)  |
+-----------------------+------------------------+---------------------------------------------+
| ``save_simu``         | False                  | Simulation saved to file for offline use    |
+-----------------------+------------------------+---------------------------------------------+
| ``save_grb``          | False                  | GRB class saved to disk -> use grb.py main  |
+-----------------------+------------------------+---------------------------------------------+
| ``datafile``          | "data.txt"             | Population study main output file           |
+-----------------------+------------------------+---------------------------------------------+
| ``logfile``           | "analysis.log"         | Text file with results, status and warning  |
+-----------------------+------------------------+---------------------------------------------+
| ``remove_tar``        | False                  | | Remove tarred files, otherwise keep for   |
|                       |                        | | faster access                             |
+-----------------------+------------------------+---------------------------------------------+


Experts and developpers only
----------------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| ``det_level``         | 0.9                    | Confidence level for detection              |
+-----------------------+------------------------+---------------------------------------------+
| ``alpha``             | 0.2                    | One over the number of on/off regions       |
+-----------------------+------------------------+---------------------------------------------+
| ``obs_point``         | "end"                  | Observation position in the time slice      |
+-----------------------+------------------------+---------------------------------------------+
| ``test_prompt``       | False                  | If True test prompt alone (experimental)    |
+-----------------------+------------------------+---------------------------------------------+
| ``use_afterglow``     | False                  | | Prompt characteritics from the afterglow  |
|                       |                        | | with same id.                             |
+-----------------------+------------------------+---------------------------------------------+
| ``tshift``            | 0                      | Shift in days applied to all trigger dates  |
+-----------------------+------------------------+---------------------------------------------+
| ``signal_to_zero``    | False                  | Keep only background, set signal to zero    |
+-----------------------+------------------------+---------------------------------------------+
| ``fixed_zenith``      | None                   | If a value (`20*u.deg`) freeze zenith in IRF|
+-----------------------+------------------------+---------------------------------------------+
| ``magnify``           | 1                      | Multiplicative factor of the flux, for tests|
+-----------------------+------------------------+---------------------------------------------+
| ``write_slices``      | False                  | Store detailed information on slices if True|
+-----------------------+------------------------+---------------------------------------------+
| ``save_dataset``      | False                  | Not implemented (save datasets)             |
+-----------------------+------------------------+---------------------------------------------+
| ``n_night``           | Null                   |  Limit data to a maximal number of nights   |
+-----------------------+------------------------+---------------------------------------------+
| ``skip``              | Null                   |  Skip the firts nights                      |
+-----------------------+------------------------+---------------------------------------------+
| ``Emax``              | Null                   |  Limit data energy bins to Emax             |
+-----------------------+------------------------+---------------------------------------------+

.. automodapi:: configuration
   :include-all-objects:
   :no-inheritance-diagram:

