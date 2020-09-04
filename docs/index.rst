SoHAPPy
=======
Simulation of High-energy Astrophysics Process in Python

SoHAPPy compute the response of cta to a population of GRB provided as fits files.
Each GRB is a set of time slices with their own spectral model.
The GRB position is assumed to be known and lies in the center of the field of view.

Introduction
------------
Two analysis methods are implemented :
    - a On-Off method using the on-axis IRF obtained for a fixed angular cut of 68% of containment. This IRF are used to produce the official CTA performance curve. They have been extracted by Julien Lefaucheur from the oriinal root files. The Off counts is thus driectly obtained from the IRF background distribution, that depends only on the energy E. The signal is obtained form the convolution of the spectrum model with the effective area. The siginificance is computed from the signal and the background counts using the Li & Ma formula with a fixed signal to background surface ratio, :math:`alpha`.
    - a 3D likelihood meeting as implemented in Gammapy, based on the full enclosure off axis IRF (the one that are available if fits for the whole consortium). The siginificance is computed as the difference inthe log-likelihood of a fit with signal and background and a fit with background only.

.. toctree::
	:maxdepth: 3

.. comment automodapi:: SoHAPPy

.. automodapi:: SoHAPPy
   :include-all-objects:
   :no-inheritance-diagram:

.. automodapi:: grb
       :no-inheritance-diagram:

.. automodapi:: visibility
   :include-all-objects:
   :no-inheritance-diagram:



