
============================================
Spectral analysis of an astrophysical source
============================================

Some tools exist to analyse the `SoHAPPy` outputs for a single source
(or a few sources individually) stored as binary file(s).
They can be found in the `analysis/spectra` subfolder.
The `spectral.py` script shows an example on doing a full analysis of
this kind of files.

The binary file is simply a binary image of the :class:`MonteCarlo` instances
for that source. They are too large to be stored in the code repository, and
could be replaced by the storage of the `Gammapy dataset` object with a lighter
format in future versions.

Introduction
============
It is possible to run `SoHAPPy` on a limited number of sources and
obtain for each of them additionnal output files to be used for a spectral
analysis. The simulation can be saved using the :code:`save_simu` keyword in
the `configuration file <configuration.html>`_ ).

The following command line can be used (One source, starting at `343`) :

``python SoHAPPy.py -f 343 -N 1 --save``

and a list of sources can be given :

``python SoHAPPy.py -f 343 980 465 676 785  -N 1 --save``

Note that the number of Monte Carlo trials is automatically forced to one.

The output file is a dump of the `mcsim` class for each of the considered site
(usually `North`, `South` and `Both`) and is therefore strictly related to the
`SoHAPPy` version that created it.
It is stored in an output folder as defined in the
`introduction <introduction.html>`_ chapter.

Analysing a source output file
==============================
The `spectral.py` script uses most of the tools available for analysing source
files.
It requires the input file to be defined through the :code:`INFILE` environment
variable otherwise it will try to load the file used for the development (which
is not distributed).
The tuning was done with this development file, and each source requires
specific analysis cuts and options. It is therefore recommended to develop the
steps suited for each source in a `jupyter notebook`, using the functions
demonstrated in the script, and the other ones in the provided modules.

The philosophy of `SoHAPPy` is to analyse the data acquired along
the observation time. For this reason, stacking the time slices can be useful.
This can be done with the function:

:code:`dataset_tools.stacking`

It returns a `datasets` with an equal number of `dataset`, each of them being
the sum of the initial preceding ones.
Note that when the `dataset`s in a `datasets` are stacked, the energy masking
is automatically applied, and the models applied to the stacked sets are all
the model of the first set. This choice is not suitable to analyse slices with
evoluting energy spectra as it is the case for GRBs. The stacking in `Gammapy`
is primarily proposed to stack consecutive observations of a source with
a unique assumption on its underlying physical model.

On top of these considerations, the data are made of slices for each of which
some information can be extracted (count numbers, mean fluxes, energy spectra).
Some of slices can be unrealistically very short (like a few seconds).
For that reason, it is possible to stack them into a minimal duration
(e.g. 10 minutes). This is done with tthe fucntion:

:code:`dataset_tools.compactify`

See more on Stacking in `Stacking`_.


Avalaible modules and functions
===============================

Some of the tools used to process the `datasets` can be found in the
`dataset_tools.py` module at the level of the `SoHAPPy` main sources (see
`<simulation.rst>`_).
The `spectral.py` script calls most of the available functions stored in the
following modules:

.. automodapi:: analysis.spectra.dataset_counts
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.spectra.dataset_flux
   :no-inheritance-diagram:
   :include-all-objects:

Technical modules
=================
These are base modules to handle the various physical objects.

.. automodapi:: analysis.spectra.dataset_plot
   :no-inheritance-diagram:
   :include-all-objects:

.. automodapi:: analysis.spectra.dataset_compare
   :no-inheritance-diagram:
   :include-all-objects:

.. _Stacking:

Consideration on stacking, masking and associated models
========================================================

A dataset, as defined in `gammapy` can carry models, in particular a spectral
model, and a `mask_safe` that defines the energy region on which the dataset
data are considered valid. It is in particular used during a fit and when
several datasets are stacked together.

The model attached to a dataset can be used to simulate the count content. It
can also be used for fitting. A model can be an analytical expression like a
powerlaw or a `TemplateSpectraModel` obtained from point measurements. In that
case the model is an interpolation through the points.

Not that fitting the data with a model different from the one used for
simulation is possible. However the initial model should be saved beforehand,
as shown:

.. code-block:: python

    model_init = ds.models
    model_fit = SkyModel(spectral_model=PowerLawSpectralModel(...))
    ds.models = model_fit
    Fit(ds)
    ds.models = model_init # Go back to initial model

The predcited counts in a dataset are computed using the model (and the IRF),
which is a continous function, at the center of the bins of the reconstrcuted
axis (although the model is defined upon true energies).
When the mask is applied, the predicted counts in the masked bins are set to
zero. The model is unchanged.

Stacking is usually intended to group several observations of a same
phenomemon. Each observation can have a different masking (e.g. varying
energy thresholds), but all observations are intended to be interpreted from
the same fit model.

When simulatiing a source with a varying energy spectrum in time, this is the
case with GRB simulations, each time slice carries its own spectral model (in
the most general case a template spectral model).

In Gammapy 0.18.2 the model in a stacked model `dsstack` obtained from two
datasets `ds1` and `ds2` is the model attached to `ds1` (as it is supposed to
be a fit model for a set of observations).

However it is possible to compute an equivalent mean model as the sum of the
original models weighted by the dataset livetime using a template spectral
model.

.. code-block:: python

    flux1 = ds1.models[0].spectral_model(e_axis.center)
    dt1   = ds1.livetime
    flux2 = ds2.models[0].spectral_model(e_axis.center)
    dt2   = ds2.livetime

    flux_new  = (dt1*flux1 + dt2*flux2)/(dt1+dt2)
    spec      = TemplateSpectralModel(energy = e_axis.center,
                             values = flux_org*mask_org,
                             norm   = 1.,
                             interp_kwargs ={"values_scale": "log"})
    dsstack.models = SkyModel(spectral_model = spec, name  = "Stack"+"-"+ds1.name)

This works perfectly as long as the two datasets are defined on the same
energy range, i.e. have the same `mask_safe` attribute. With this, the
predicted counts, that are computed from the model, can still be compared to
the generated counts in a simulation

If the two datasets have NOT the same `mask_safe` attribute, one could set to
zero the flux point in the definition of the template spectral model as
follows. Note that this requires to define the model on the reconstructed
energy axis (axis on whih the `mask _safe` is defined):

.. code-block:: python

    e_axis = ds1.background.geom.axes[0] # E reco sampling from the IRF
    mask1 = ds1.mask_safe.data.flatten()
    flux1 = ds1.models[0].spectral_model(e_axis.center)*mask1 # masked bins at zero
    dt1   = ds1.livetime

However the resulting interpolated spectrum can/will have edge effects that
will not guarantee that the flux is stricly at zero. As a consequence it is
not strictly possible to propagate a mean effective model when stackig two
datasets.

The solution to check that the simulated count numbers follow the original
modelling is to cumulate by hand the predicted counts from the original
individual datasets, this is how ot os implemented in `SoHAPPy`.

#.. automodapi:: singlesource_analysis
#   :include-all-objects:
#   :no-inheritance-diagram: