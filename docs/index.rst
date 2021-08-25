#######
SoHAPPy
#######

*Simulation of High-energy Astrophysics Process detection in Python*

`SoHAPPy` computes the response of `CTA <https://www.cta-observatory.org>`_ to a population of astrophysical objects 
-initially GRBs- described in *fits* or *yaml* files in the simplest cases. 
Each object is a set of time slices with their respective spectral model.
The object position is assumed to be known with an uncertainty below the instrument precision.

`SoHAPPy` handles the visibility from a start (trigger) date.
Its main output are either a summary *csv* file with the characterisitics of the detection (significance and position in the sky) and individual output files allowing a time and spectral analysis. Both output files can be analysed thanks to dedicated python notebooks.

The SoHAPPy project started mid 2019. 
Check the `gihub repository <https://github.com/tstolarczyk/SoHAPPy>`_ for the latest improvements.

More on SoHAPPy
===============

.. toctree::
	:maxdepth: 1

	introduction.rst
	sohappy.rst
	configuration.rst
	irf.rst
	file_format.rst

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
	
.. toctree::
	:maxdepth: 1
	:caption: Tools

	tools.rst

.. toctree::
	:maxdepth: 1
	:caption: The developer corner

	developer.rst

.. rubric:: References
	[tot] toto



Producing the documentation
===========================

* See a `RST syntax tutorial <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
* From the docs SOHAPPy folder, execute:
    * windows : .\\make.bat html
    * Linux : make html 

