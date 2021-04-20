SoHAPPy with GammaPy,
Simulation of High-energy Astrophysics Processes in Python
roma version (gammapy 0.12)

Installation
============
SoHAPPy will easily run under an anaconda environement, requiring essentially
the gammapy package.
The Roma and Sofia versions are compatible with gammapy 0.12.

To create and activate the proper environment do :
	conda create -c conda-forge -n gpy0.12 gammapy=0.12
	conda activate gpy0.12
	conda install matplotlib

A created environment can be removed as follows :
conda env list
conda env remove -n gpy0.12

dependecy
pandas 1.03
...

Run the code
============
Type python SoHAPPy.py -h to get help on the command line
Any omitted argument is taken from the ana_config.py file.
The recommanded defautf values are listed below.

For a large number of events the debugging should set to zero.
In that case no plots are produced and the matplotlib is not required

Creating the documentation
==========================
The documentation is created in the docs folder using the Sphynx package.
The docs/index.rst file contains the documentation header and the list
allowing the modules to be part of the documentation.

Help on rst syntax can be found here:
https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

Install sphinx and some extensions sphinx-automodapi:
conda install -c conda-forge sphinx

conda install -c astropy sphinx-automodapi
pip install nbsphinx
pip install recommonmark
pip install sphinx_rtd_theme

Then execute:


