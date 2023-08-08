#######
SoHAPPy
#######

*Simulation of High-energy Astrophysics Process detection in Python*

`SoHAPPy` computes the response of `CTA <https://www.cta-observatory.org>`_ to a population of astrophysical objects
-initially GRBs- described in *fits* or *yaml* files in the simplest cases.
Each object is a set of time slices with their respective spectral model.
The object position is assumed to be known with an uncertainty below the instrument precision.

`SoHAPPy` handles the visibility from a start (trigger) date.
Its main output are either a summary text file with the characteristics of the detection (significance and position in the sky) or individual files allowing a time and spectral analysis. Both output files can be analysed thanks to dedicated python notebooks.

The `SoHAPPy` project started mid 2019.
Check the `gihub repository <https://github.com/tstolarczyk/SoHAPPy>`_ for the latest improvements.

More on SoHAPPy
===============

.. toctree::
   :maxdepth: 1
   :caption: Introduction

   introduction.rst
   sohappy.rst
   configuration.rst
   irf.rst
   file_format.rst

..
   test.rst

.. toctree::
   :maxdepth: 1
   :caption: Concepts and classes

   grb.rst
   time_slice.rst
   visibility.rst
   simulation.rst

.. toctree::
   :maxdepth: 1
   :caption: Population analysis

   population_analysis.rst

.. toctree::
   :maxdepth: 1
   :caption: Spectral analysis

   spectral_analysis.rst

.. toctree::
   :maxdepth: 1
   :caption: Tools

   tools.rst

.. toctree::
   :maxdepth: 1
   :caption: Running large simulations and analysis

   production.rst

.. toctree::
   :maxdepth: 1
   :caption: The developer corner

   developer.rst

.. rubric:: References
   [ref Put here a reference



Producing the documentation
===========================

* See a `RST syntax tutorial <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
* From the docs SOHAPPy folder, execute:
    * windows : .\\make.bat html
    * Linux : make html

