#######
SoHAPPy
#######

*Simulation of High-energy Astrophysics Process in Python*

`SoHAPPy` computes the response of `CTA <https://www.cta-observatory.org>`_ to a population of astrophysical objects 
-initially GRBs- described in *fits* or *yaml* files in the simplest cases. 
Each object is a set of time slices with their respective spectral model.
The object position is assumed to be known with an uncertainty below the instrument precision. It has a fixed position in the field of view, i.e. the telescope array constantly follows he source during the observation window.
`SoHAPPy` handles the visibility from a start (trigger) date.
Its main output are either a summary *csv* file with the characterisitics of the detection (significance and position in the sky) and individual output files allowing a time and spectral analysis. Both output files canbe analysed thanks to dedicated python notebooks.

.. toctree::
	:maxdepth: 1

	steering.rst
	irf.rst

.. toctree::
	:maxdepth: 1
	:caption: Concepts and classes
	
	grb.rst
	time_slice.rst
	visibility.rst
	montecarlo.rst
	
.. toctree::
	:maxdepth: 1
	:caption: Analysis
	
	population_analysis.rst
	source_analysis.rst

Dependencies
============

The code uses the following external components and code:
    * gammapy: current supported version is 0.18.2
    * astroplan version 0.8
    * Intrument response function (IRF) files organised in a particular subfolder structure (see below)
    * input object files : at the moment they follow a format as described in [ref].

Repository structure
====================
The SoHAPPy repository has 2 subfolders:
    - analysis: contains notebooks and scripts to analyse the SoHAPPy outputs.
    - docs: contains the Sphinx information to generate the present documentation.

The IRF and source files are expected to be obtained from specific folders.
The IRF file are organised in a particular way (see IRF section).

The output consists in a few files packed together into a tar.gz archive.
They are written in a specific output folder.


Launching the code
==================
The steering parametres are obtained from the config.yaml file. 
The most important parametres can be passed on the command line (they supersede the ana_config.py parameters which are considered as the default).
Type:
``python SoHAPPy.py -h``
to have the list explicited (to be completed).

"SoHAPPy -h" in a shell gives the list of parameters that can be superseded
on the command line

The processing of 1000 GRB files with one iteration per trial takes roughly 2 
hours on a I7 laptop, including the visibility computation.
With 100 iterations the computing time goes up to 

How it works ?
==============
Here is an overview of the processes handled by the main function (in `SoHAPPy.py`):

* Read the steering parmateres from a *yaml* file (`config.yaml`) and/or the command line;

* Create the object list to be analysed, either a large number for a populations study or only a few for specific spectral analysis;

* Copy the configuration *yaml* file to the output folder and open the output files (log and population files) and make a copy of the configuraton file for further reference.

* Loop over the object list:
	* Read the object information, create a :class:`GammaRayBurst` object (see :class:`grb.GammaRayBurst`, a list of interpolated EBL abosorbed spectrum)

	* Create an original time :class:`Slot` with the original object time slices (:class:`timeslot.Slot` and :class:`obs_slice.Slice`)

	* For each site (*North*, *South*) assign a visibility to the current object (see :class:`visibility.Visibility`), either:

		- from the default given in the input file;
		- computed on the fly from the start date;
		- read from a previously computed visibilty on disk.

	* If requested save the :class:`GammaRayBurst` instance on disk.
	
	* Loop over the sites (*North*, *South*):
	
		- Create a :class:`MonteCarlo` class instance (see :class:`mcsim.MonteCarlo`)
		- Copy the original :class:`Slot` and modify the content to account for the computed visibility;
		- Associate a physical flux (:class:`timeslot.Slot.dress`) to each time slices of the :class:`Slot`.
		- Run the simulationon the resulting dressed :class:`Slot`
		- Display and store the results for the current object into the population file and/or the spectral analysis file.
		- If requested store the simulation class instance on disk.
	
	* For objects visible on both sites, run a combined simulation:
		
		- Create a :class:`MonteCarlo` class instance (see :class:`mcsim.MonteCarlo`)
		- Copy the original :class:`Slot` and modify the content to account for the computed visibilities;
		- Associate a physical flux (:class:`timeslot.Slot.dress`) to each time slices of the :class:`Slot`.
		- Run the simulationon the resulting dressed :class:`Slot`
		- Display and store the results for the current object into the population file and/or the spectral analysis file.		
		- If requested store the simulation class instance on disk.

	* Close the population and log file.
	
.. comment automodapi:: SoHAPPy

.. automodapi:: SoHAPPy
   :include-all-objects:
   :no-inheritance-diagram:



Producing the documentation
===========================
- See a `RST syntax tutorial <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
- From the docs SOHAPPy folder, execute:
    - windows : .\make.bat html
    - Linux : make html 