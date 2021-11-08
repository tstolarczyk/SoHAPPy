How it works ?
==============
Here is an overview of the processes handled by the main function (in `SoHAPPy.py`):

* Read the steering parmaters from a *yaml* file (`config.yaml`) and/or the command line;

* Create the object list to be analysed, either a large number for a populations study or only a few for specific spectral analysis;

* Copy the configuration *yaml* file to the output folder and open the output files (log and population files) and make a copy of the configuraton file for further reference.

* Loop over the object list:
	* Read the object information, create a :class:`GammaRayBurst` object (see :class:`grb.GammaRayBurst`, a list of interpolated EBL absorbed spectra)

	* Create an original time :class:`Slot` with the original object time slices (:class:`timeslot.Slot` and :class:`obs_slice.Slice`)

	* For each site (*North*, *South*) assign a visibility to the current object (see :class:`visibility.Visibility`), either:

		- from the default given in the input file;
		- computed on the fly from the start date;
		- read from a previously computed visibilty on disk.

	* If requested save the :class:`GammaRayBurst` instance on disk.

	* Loop over the sites (*North*, *South*):

		- Create a :class:`MonteCarlo` class instance (see :class:`mcsim.MonteCarlo`)
		- Copy the original :class:`Slot` and modify the content to account for the computed visibility (:func:`timeslot.Slot.apply_visibility`;
		- Associate a physical flux (:func:`timeslot.Slot.dress`) to each time slices of the :class:`Slot`.
		- Run the simulationon the resulting dressed :class:`Slot`
		- Display and store the results for the current object into the population file and/or the spectral analysis file.
		- If requested store the simulation class instance on disk.

	* For objects visible on both sites, run a combined simulation:

		- Create a :class:`MonteCarlo` class instance (see :class:`mcsim.MonteCarlo`)
		- Copy the original :class:`Slot` and modify the content to account for the computed visibilities;
		- Associate a physical flux (:func:`timeslot.Slot.dress`) to each time slices of the :class:`Slot`.
		- Run the simulationon the resulting dressed :class:`Slot`
		- Display and store the results for the current object into the population file and/or the spectral analysis file.		
		- If requested store the simulation class instance on disk.

	* Close the population and log file.

.. comment automodapi:: SoHAPPy

.. automodapi:: SoHAPPy
   :include-all-objects:
   :no-inheritance-diagram:
   