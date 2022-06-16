How it works ?
==============
Here is an overview of the processes handled by the main function (in `SoHAPPy.py`):

* Read the steering parmaters from a *yaml* file (`config.yaml`) and/or the command line;

* Copy the configuration *yaml* file to the output folder and open the output files (log and population files) and make a copy of the configuraton file for further reference.

* Create the object list to be analysed, either a large number for a populations study or only a few for specific spectral analysis;

* Loop over the object list:
	* Read the object information, create a :class:`GammaRayBurst` object (see :class:`grb.GammaRayBurst`, a list of interpolated EBL absorbed spectra) including the associated visibility (see :class:`visibility.Visibility`) either from the default given in the input file, computed on the fly from the start date, or read from a previously computed visibilty stored on disk..

	* Create an original time :class:`Slot` with the original object time slices (:class:`timeslot.Slot` and :class:`obs_slice.Slice`)

	* If requested save the :class:`GammaRayBurst` instance on disk.

	* Get the delays on all sites from the slewing (CTA) and latency (external alert system) delays.

	* Loop over the sites (*North*, *South*):

		- Create a :class:`MonteCarlo` class instance (see :class:`mcsim.MonteCarlo`)
		- Copy the original :class:`Slot` and modify the content to account for the computed visibility (:func:`timeslot.Slot.apply_visibility`;
		
		- If the GRB is still visible on the current site:
		
			- Associate a physical flux (:func:`timeslot.Slot.dress`) to each time slices of the :class:`Slot`.
			- Run the simulation on the resulting dressed :class:`Slot`
			
		- Display and store the results for the current object into the population file and/or the spectral analysis file (Valid even if the GRB was not visible and the simulation was not run).
		- If requested store the simulation class instance on disk.

	* For objects visible on both sites, run a combined simulation:

		- Create a :class:`MonteCarlo` class instance (see :class:`mcsim.MonteCarlo`)
		- If the source is visible on all sites:
		
			- Copy the original :class:`Slot` and modify the content to account for the computed visibilities;
			- Associate a physical flux (:func:`timeslot.Slot.dress`) to each time slices of the :class:`Slot`.
			- Run the simulationon the resulting dressed :class:`Slot`
			
		- Display and store the results for the current object into the population file and/or the spectral analysis file (even if the simulation was not run).		
		- If requested store the simulation class instance on disk.

	* Close the population and log file.
	* `tar` the output files in the output folder.

.. comment automodapi:: SoHAPPy

.. automodapi:: SoHAPPy
   :include-all-objects:
   :no-inheritance-diagram:
   