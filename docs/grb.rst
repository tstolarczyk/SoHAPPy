===========================
Astrophysical objects (GRB)
===========================

The input data files contain a series of flux measurements at certain energies
for a series of time slices. The files also give the position, the redshift,
some physical information. In some cases it also contain a default visibility
over one night for a minimal altitude of 10 degrees (this visibility can be
recomputed, see later). Alternatively the source spectrum can be  described
as  energy and time powerlaws in an external `yaml` file that will be used to
generate energy spectra and time slices.

The source information is handled by the :class:`grb.GammaRayBurst` class.
Each instance carries a list in time of models with a default Spatial model
(point-like source at the position given, :obj:`PointSpatialModel`) and an
energy spectrum stored as template models (:obj:`TemplateSpectralModel`),
which are interpolations of the initial flux points. If required an EBL
absorption is applied at the redshift value given (usual case). Alternatively
the initial file can contain already absorbed spectra (in that case the model
as to be declared `built-in` in the configuration file).

To each object are attached the visibilities at each site, thanks to the :func:`grb.GammaRayBurst.set_visibility` function. They are used to compute
the time slot used for the simulation and analysis.


The `grb` module
================

This module contains only the :class:`grb.GammaRayBurst` class and a
main function with examples.

.. autoclass:: grb.GammaRayBurst
    :special-members: __init__
    :members:

The `ebl` module
================

.. automodapi:: ebl
   :no-inheritance-diagram:
   :include-all-objects: