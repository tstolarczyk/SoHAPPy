How it works ?
==============

Here is an overview of the processes handled by the main function
(in `SoHAPPy.py` ):

* Read the steering parameters from a `yaml` file (`data/config_ref.yaml` by
  default) and/or the command line (see `configuration <configuration.html>`_).
* Create the specific output folder.
* Copy an updated version of the configuration *yaml*.
* Create the object list to be processed.
* Open the output files (log and population files).
* Create the object list to be analysed, either a large number for a
  populations study or only a few for specific spectral analysis.
* Decode the `visibility <visibility.html>`_ information returning either a
  python dictionnay or a keyword to be treated later.
* Loop over the object list:

	* Read the object information, create a :class:`grb.GammaRayBurst` object
	  (a list of interpolated EBL absorbed spectra).
	* Set and/or compute the visibility for that object for the all sites
	  (see :class:`visibility.Visibility`).
	* If requested save the :class:`grb.GammaRayBurst` instance on disk.
	* Get the delays on all sites from the slewing (CTA) and latency
	  (external alert system) delays.
	* Create an original time Slot with the original object time slices
	  (:class:`timeslot.Slot` and :class:`obs_slice.Slice`, read more on
	  time `slots and slices <time_slice.html>`_ ).
	* Loop over the site configuration (*North*, *South*, *Both*), see below.
	* Close the population and log file.
	* `tar` the output files in the output folder.

For each site and the combination of the two, simulate the GRB as many times
as required and analyse the results. Here are the steps
(see more on `simulation and analysis <simulation.html>`_):

* Create a :class:`mcsim.MonteCarlo` class instance.
* Copy the original :class:`timeslot.Slot` and modify the content to account
  for the computed visibility (:func:`timeslot.Slot.apply_visibility` ).
* If the GRB is still visible on the current site:

	* Associate a physical flux (:func:`timeslot.Slot.dress`) to each time
	  slices of the :class:`timeslot.Slot`.
	* Create and initialize an analysis object (:class:`analyze.Analysis`)
	* Run the simulation on the resulting dressed :class:`timeslot.Slot` and
	  prepare the data to be analysed.
	* If requested store the simulation class instance on disk.
	* If the simulation is successful, analyse the resulting data.
	* Dump the result of the analysis into the population file.

.. automodapi:: SoHAPPy
   :include-all-objects:
   :no-inheritance-diagram:
