Monte Carlo
===========

The Monte Carlo class (MonteCarlo in mcsim.py) is the deepest class of SoHAPPy. It runs simulations over a slot, perform the analysis, obtain the results. It is associated to a site or several (two) sites at once.
It is associated to the `mcsim_config.py` and `mcsim_res.py` modules that contain simulation parameters and output tools respectively.
It carries a Slot instantiation and the following variables resulting from the chosen number of trials:
- the list of maximum significances reached and the corresponding observation times and excess and background counts;
- the list of times at which 3 and 5 sigma were reached and the associated excess and background count numbers.
- a list of slices for which warning were issued for further debugging (not intended to the normal user).
- a list of datasets (i.e. a python list of gammapy `Dataset` objects, not a `Datasets` gammapy object).

The dataset list is used to store the `Slot` time slice information, compute the expected and generated signal and background counts, and hence the on and off counts and the corresponding significances.
When only one realisation is requested the MonteCarlo instantiation can be stored on disk and re-used within a dedicated analysis script (see later).

There are two main ways to use the :meth:`run` method:
	- Compute only once the number of signal and background for each time slice. This is in pricniple done without fluctuation and allow getting rapidly results on a large population. In particular the 5 sigma detected sub-population is fairly unbiased compared to the one obtained with many iterations since it usually concerns luminous objects with large number of signal events (however the number of background event can be small enough to marginally change the result. This can be checked more quantitavely). This kind of simulation can be saved to disk for further use;
	- The normal use is to generated many trials -of the order of 100- fluctuating the signal and background at each time slice of the `Slot`, in order to address, in particular, objects for which the 5 sigma detection is not granted depending on the count fluctuations. 

The associated `mcsim_config.py` script handles parameter that have been tuned in the case of the official CTA prod3b
production and are not intended to be changed without a thorough independent study.

.. automodapi:: mcsim
   :include-all-objects:
   :no-inheritance-diagram:
   
.. automodapi:: mcsim_config
   :include-all-objects:
   :no-inheritance-diagram:

.. automodapi:: mcsim_res
   :include-all-objects:
   :no-inheritance-diagram:

.. automodapi:: mcsim_plot
   :include-all-objects:
   :no-inheritance-diagram: