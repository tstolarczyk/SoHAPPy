Astrophysical objects (GRB)
===========================

The input files contain a series of flux measurement at certain energies for a series of time slices. This files also give the position, the redshift, some physical information and the default visibility over one night for a minimal altitude of 10 degrees (this visibility can be recomputed, see later).Alternatively the object can described as  energy and time powerlaws in extranl yaml file that will be used to generate energy spectra and time slices.

The source information is handled by the `~grb.GammaRayBurst` class in the ``grb.py`` module.
Each ``GammaRayburst`` object carries a list in time of *WRONG* ``AbsorbedSpectralModel`` models obtained from template models, which are interpolation of the initial flux points after applying a given EBL absorption.
Alternatively the initial file can contain already absorbed spectra (in that case the model as to be declared "built-in" in ana_congif.py).

To each object are attached the visibilities at each site (typically North and South for CTA). They are used to compute the time slot used for the simulmation and analysis. 

The grb_plot.py modules contains the fucntions displaying the grb charcateristics, including the visibility windows.

.. automodapi:: grb
   :no-inheritance-diagram:
  
.. automodapi:: grb_plot
    :no-inheritance-diagram: