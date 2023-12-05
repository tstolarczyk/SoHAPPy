===================
Population analysis
===================

Introduction
============

The population is analysed from the devoted `SoHAPPy` output files
(named `data.txt` by default).

The available tools consist in python script files with the most common
analysis algorithms and plot facilities. The user is intended to indicate
in the file `parameter.yaml` the `SoHAPPy` subfolder where the data can be
found as well as the duration in years associated to the simulation. This
duration is used to  normalize the detection rates.

..  code-block:: python

    duration:
        - 44
    outfolders:
        - "long/alpha_prod5/pop_vis24_strictmoonveto-100iter"
    tags:
        - "Example"

The subfolder can contain a data file, or a `tar.gz` archive containing it. The
subfolder can also contain folders resulting from the batch processing of a
series of subpopulation as described in `<production.rst>_`. In that case the
data will be found in all folders and merged. This does mean that they have
to come from the same unique population simulation.

The initial data file are readable text files. They are converted into a `.csv`
file opened with `pandas`.

The input and output process is handled by the functions in ``pop_io.py``
(see documentation below).

The following scripts are available for the analysis:


.. automodapi:: analysis.population.population
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.rate
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.pop_plot
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.sigmas
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.times
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.skymaps
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.universe
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.compare_sets
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.compare_ebl_models
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.explore_simulation
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.population.setup
   :no-inheritance-diagram:
   :include-all-objects:



