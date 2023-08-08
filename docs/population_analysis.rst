Population analysis
===================

Introduction
------------

The population is analysed from the devoted `SoHAPPy` output file
(named `data.txt` by default).

The available tools consist in a python script files and a few jupyter
notebooks. The user is intended to indicate the `SoHAPPy` folder containing
the population file in a file named `parameter.yaml`, with the variable
`outfoler`, like shown here:

..  code-block:: python

	# outfolder : "../../../output/pop_vis24_moonlight/"
	outfolder : "../../../output/pop_vis24_fullmoonveto/"

* **pop_io.py**: Contains functions to read and convert the input data file and
  the associated configuration parameter file, either from disk or extracted
  from atarred compressed archive.

* **population.py** :
  Contains the definition of the class Population. The contructor function
  (__init__) opens the input file and store the read data in various
  pandas table that are reused later. The main function reads the data,
  perform a sanity check comparing the results to the requirement in the
  original configuration file, check for negative significances, the fraction
  of GRB for which the prompt component can be detected etc.

* **rate.py**:
  Compute detection rates for various site-related population; normalised to
  the data taking duration to get mean yearly rates above 3 and 5 standard
  deviations of detection significance.

* **pop_plot.py**:
  Show distributions and population coverage as functions of Eiso, Epeak and z.

* **sigmas.py**:
  Detection rates as a function of the mean maximal significance thresholds,
  evolution of rates with this threshold, distribution of mean maximal significances.

* **times.py**:
  Characterizes the time needed to get a detection or the mean maximal
  significance. Plots in particular for...

* **skymaps.py**:
  Show skymap coverage in ra-dec and in altitude-azimuth

* **universe.py**:
  Show universe coverage.

* **setup.py**:
  Gathers some defintions for plot colors


Function documentation
----------------------

.. automodapi:: analysis.population.pop_io
   :include-all-objects:
   :no-inheritance-diagram:


.. automodapi:: analysis.population.population
   :include-all-objects:
   :no-inheritance-diagram:


.. automodapi:: analysis.population.rate
   :include-all-objects:
   :no-inheritance-diagram:


.. automodapi:: analysis.population.pop_plot
   :include-all-objects:
   :no-inheritance-diagram:


.. automodapi:: analysis.population.sigmas
   :include-all-objects:
   :no-inheritance-diagram:


.. automodapi:: analysis.population.times
   :include-all-objects:
   :no-inheritance-diagram:

.. automodapi:: analysis.population.skymaps
   :include-all-objects:
   :no-inheritance-diagram:

.. automodapi:: analysis.population.universe
   :include-all-objects:
   :no-inheritance-diagram:

.. automodapi:: analysis.population.setup
   :include-all-objects:
   :no-inheritance-diagram:
