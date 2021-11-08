Introduction
############

Dependencies
============

The code uses the following external components and data:
    * `gammapy <https://gammapy.org/>`_: current supported version is 0.18.2
    * `astroplan <https://pypi.org/project/astroplan/>`_ : version 0.8
    * `Official CTA Intrument response function <https://www.cta-observatory.org/cta-performance-prod3b-v2/>`_ (IRF) files, organised in a particular subfolder structure (see `IRF <irf.rst>`_)
    * input object files with a `dedicated format <file_format.rst>`_.

Repository structure
====================
The SoHAPPy repository has 2 subfolders:
    - *analysis*: contains notebooks and scripts to analyse the SoHAPPy outputs.
    - *docs*: contains the Sphinx information to generate the present documentation.

| The IRF and source files are expected to be obtained from specific folders.
| The output consists in a few files packed together into a tar.gz archive. They are written in a specific output folder.

Launching the code
==================
The steering parametres are obtained from the yaml configuration file (`conf.yaml` is used by default if it is present). 
Two example are given by default:
* config_population.yaml, for a large population study, with minimal debugging output
* config_singlesource.yaml, with the all processing messages and the infiviadual simulation saved on disk for further use. 

The most important parametres can be passed on the command line, in particular the configuration file name (they supersede parameters of the configuration fle).
Type:
``python SoHAPPy.py -h``
to have the list explicited.

The processing of 1000 GRB files with one iteration per trial takes roughly 2 
hours on a I5 -16GB laptop, including the visibility computation.
With 100 iterations the computing time is only slightly increased.

Producing the visibility alone and writing it to disk takes 2hr.
Running SoHAPPy with pre-computed visibilities takes 35 minutes to one hour (when no Moon veto)with one trial and 50 minutes with 100 trials.