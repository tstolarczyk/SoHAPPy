Input file format
=================

This should be revised as the column name have changed, as well as the format.

Fits file
---------
The fits files used in SoHAPPy follow the following format:

.. figure:: figures/01-fits_content.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

   File contents

.. figure:: figures/02-fits_header.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

   Header

.. figure:: figures/03-fits_energy_bins.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

   Energy bins

.. figure:: figures/04-fits_time_spectra.jpg
   :alt: alternate text
   :scale: 80 %
   :align: center

   Time spectra

Yaml file
---------

If the source filename is a string (raher than an integer identifying a source
in a list), then `SoHAPPy` try to process a `yaml` file containing information
to parameterise the flux as a time and energy powerlaws. These files are
expected to be in the `data/historical` `SoHAPPy` subfolder. More files can be
added. They are usually built from real data found in article analysing past
GRBs.

The flux is described as:

.. code-block:: python

    flux = K*(E/E0)**-gamma*(t/t0)**-beta

Here is an example for one of the brightest Fermi/LAT GRB.

.. code-block:: python

	# Fermi GRB catalog : https://fermi.gsfc.nasa.gov/ssc/observations/types/grbs/lat_grbs/
	---
	name  : 080916C
	E0    : 1 TeV # reference energy
	t0    : 35 s # reference time
	gamma : 1.85 # Energy index
	beta  : 1.4  # Time index
	K     : 0.993e-7 1/(TeV s cm2) # Normalisation
	z     : 4.3 # Redshift
	t_trig: 2008-09-16T00:12:45.0 # Explosion time
	ra    : 119.84717 deg # Right ascencion, galactic coordinate
	dec   : -56.638333 deg # Declination, galactic ccordinate
	Eiso  : 8.80E+54 erg # Emitted energy assuming and isotropic emission
	Epeak : 0 keV # Prompt peak energy
	tmin  : 1 s # Start of simulation
	tmax  : 2 d # End of simulation
	Emin  : 10 GeV # Minimum simulated energy
	Emax  : 100 TeV # Maximal simulated energy
	ref   : "Abdo et al., Science 323, 2009" # reference if any
	gcn   : "https://gcn.gsfc.nasa.gov/other/080916C.gcn3" # Link to GCN alert

