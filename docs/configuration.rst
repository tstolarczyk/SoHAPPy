Configuration file
==================
A simulation and analysis parameter file can be passed to the main script with
the `-c` options (some of these parameters can be superseded on the command
line, use  `-h` to get the list). If no file name is given, the parameters
are read from the *data/config_ref.yaml* default file provided with the code
sources.

The base input and output folders are defined through the `HAPPY_IN` and
`HAPPY_OUT` environment variables respectively.

The execution of `SoHAPPy` automatically creates a copy (backup) of the chosen
configuration file to the output folder, with the initial parameters
overwritten by the ones on the command line, so that this backuped
configuration can be used to reprocess the data in the exact same conditions.

The configuration data are stored and processed in the
:class:`configuration.Configuration` class in :file:`configuration.py`.

The following parameters are defined and mandatory in the configuration file.
Extra parameters can be added and are notified as ignored at run time. The
defauilt in the following tables are the one in the *data/config_ref.yaml*
file.

Physics parameters
------------------

.. tabularcolumns:: |l|c|p{5cm}|

+----------------+------------------+--------------------------------------------------+
| variable       | Default          | What is it ?                                     |
+================+==================+==================================================+
| ``ifirst``     | 1                | | The first GRB identifier to be processed       |
|                |                  | | or a list of identifiers or a json file with   |
|                |                  | | paths to files to be processed.                |
|                |                  | | The identifiers can include a string referring |
|                |                  | | to a `yaml` file with parameters to generate   | 
|                |                  | | spectra from an analytical function.           |
|                |                  | | (See the `data/historical` folder)             |
+----------------+------------------+--------------------------------------------------+
| ``nsrc``       | 1                | | Number of GRB to be processed if ``ifirst`` is |
|                |                  | | an integer. If 1, special actions are taken.   |
|                |                  | | Not used if ``ifirst`` is a list.              |
+----------------+------------------+--------------------------------------------------+
| ``visibility`` | "strictmoonveto" | | Can be:                                        |
|                |                  | | * `built-in` (read from the data file if it    |
|                |                  | | exists),                                       |
|                |                  | | * a keyword pointing to a folder containing    |
|                |                  | | `json` files with a collection of visibilities;|
|                |                  | | * a keyword corresponding to a dictionnay entry|
|                |                  | | in `visibility.yaml` to compute the visibility |
|                |                  | | on the fly (example given)                     |
|                |                  | | * `forced` to force a single infinite night    |
|                |                  | | * `permanent` to force a permanent visibility  |
+----------------+------------------+--------------------------------------------------+
| ``ebl_model``  | "dominguez"      | | Extragalactic Background Light model as        |
|                |                  | | defined in ``gammapy`` or `built-in` if        |
|                |                  | | provided or `None` if ignored (no absortpion)  |
+----------------+------------------+--------------------------------------------------+
| ``maxnight``   | Null             |  Limit data to a maximal number of nights        |
+----------------+------------------+--------------------------------------------------+
| ``skip``       | Null             |  Skip the first nights                           |
+----------------+------------------+--------------------------------------------------+


Input / Output file names and directories
-----------------------------------------

.. tabularcolumns:: |l|c|p{5cm}|

See file for examples.

+------------------+-------------------+--------------------------------------------------+
| variable         | Default           | What is it ?                                     |
+==================+===================+==================================================+
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
| ``prefix``       | "event_"          | | The prefix before the source identifier        |
|                  |                   | | to create the full file name                   |
+------------------+-------------------+--------------------------------------------------+
| ``suffix``       | "long_ISM.fits"   | | The suffix after the source identifier to      |
|                  |                   | | create the file name (including the extension) |
+------------------+-------------------+--------------------------------------------------+
| ``dgt``          | 5                 | | Number of digit on which the integer refereing |
|                  |                   | | to anidentifier is encoded (e.g. 00001 for 1)  |
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

(*) Note that this bias the resulting population since it artificially depletes
the maximal significance population below the minimum required (e.g. 3 sigma).

Detection parameters
--------------------

.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default                | What is it ?                                |
+=======================+========================+=============================================+
| ``observatory``       | "CTAO"                 | Observatory name (for position on Earth)    |
+-----------------------+------------------------+---------------------------------------------+
| ``array_North``       | "4LSTs09MSTs"          | IRF North array, subfolder name             |
+-----------------------+------------------------+---------------------------------------------+
| ``array_South``       | "14MSTs37SSTs"         | IRF South array, subfolder name             |
+-----------------------+------------------------+---------------------------------------------+
| ``dtslew_North``      | "30 s"                 | Maximum slewing time delay in North         |
+-----------------------+------------------------+---------------------------------------------+
| ``dtslew_South``      | "90 s"                 | Maximum slewing time delay in South         |
+-----------------------+------------------------+---------------------------------------------+
| ``fixslew``           | True                   | If False generate a random delay less       |
|                       |                        |   than ``dtslew``. If True use ``dtswift``  |
+-----------------------+------------------------+---------------------------------------------+
| ``dtswift``           | "77 s"                 | | Alert latency Fixed SWIFT latency         |
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
| ``dbg``               | 0                      | From 0 to 3, increasingly talkative         |
+-----------------------+------------------------+---------------------------------------------+
| ``silent``            | False                  | If True, nothing on screen (output to log)  |
+-----------------------+------------------------+---------------------------------------------+
| ``save_simu``         | False                  | Simulation saved to file for offline use    |
+-----------------------+------------------------+---------------------------------------------+
| ``save_grb``          | False                  | GRB class saved to disk -> use grb.py main  |
+-----------------------+------------------------+---------------------------------------------+
| ``save_fig``          | False                  | Figures saved to disk                       |
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
| ``elimit``            | Null                   | Limit energy bins for CPU saving (Quantity) |
+-----------------------+------------------------+---------------------------------------------+
| ``tlimit``            | Null                   | Limit time bins for CPU saving (Quantity)   |
+-----------------------+------------------------+---------------------------------------------+
| ``emin``              | Null                   | Analysis Energy lower limit (Quantity)      |
+-----------------------+------------------------+---------------------------------------------+
| ``emax``              | Null                   | Analysis Energy upper limit (Quantity)      |
+-----------------------+------------------------+---------------------------------------------+
| ``e_dense``           | Null                   | Use a denser E-binning for spectral analysis|
+-----------------------+------------------------+---------------------------------------------+
| ``tmin``              | Null                   | Analysis time lower limit (Quantity)        |
+-----------------------+------------------------+---------------------------------------------+
| ``tmax``              | Null                   | Analysis time  upper limit (Quantity)       |
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
| ``fixed_zenith``      | None                   | If a value (`20*u.deg`) freeze zenith in IRF|
+-----------------------+------------------------+---------------------------------------------+
| ``magnify``           | 1                      | Multiplicative factor of the flux, for tests|
+-----------------------+------------------------+---------------------------------------------+
| ``write_slices``      | False                  | Store detailed information on slices if True|
+-----------------------+------------------------+---------------------------------------------+
| ``save_dataset``      | False                  | Not implemented (save datasets)             |
+-----------------------+------------------------+---------------------------------------------+


.. automodapi:: configuration
   :include-all-objects:
   :no-inheritance-diagram:

