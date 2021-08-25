Input file format
=================

Fits file
---------
The fits files used in SoHAPPy follow the following format:

.. image:: figures/01-fits_content.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

File contents

.. image:: figures/02-fits_header.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

Header

.. image:: figures/03-fits_energy_bins.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

Energy bins

.. image:: figures/04-fits_time_spectra.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

Time spectra

Yaml file
---------

An object parametrised as a time and energy powerlaws can be described in a yaml file. here is an example for one of the brightest Fermi/LAT GRB.


.. code-block::

	# Fermi GRB catalog : https://fermi.gsfc.nasa.gov/ssc/observations/types/grbs/lat_grbs/
	---
	name  : 080916C
	E0    : 1 TeV
	t0    : 35 s 
	gamma : 1.85
	beta  : 1.4
	K     : 0.993e-7 1/(TeV s cm2)
	z     : 4.3
	t_trig: 2008-09-16T00:12:45.0
	ra    : 119.84717 deg
	dec   : -56.638333 deg
	Eiso  : 8.80E+54 erg
	Epeak : 0 keV
	tmin  : 1 s
	tmax  : 2 d
	Emin  : 10 GeV
	Emax  : 100 TeV
	ref   : "Abdo et al., Science 323, 2009"
	gcn   : "https://gcn.gsfc.nasa.gov/other/080916C.gcn3"

