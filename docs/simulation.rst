========================
Simulation and Analysis
========================

    - `Simulation`_ of the instrumentation response
        - The `mcsim`_ module
        - `mcsim_config`_ module

    - `Analysis`_ of simulated data
        - The `analyze`_ module

    -`Dataset`_ tools

The simulation and analysis is based on a `on-off` method, where a background
below a source is estimated from regions without signal with the
same acceptance as the signal (and underneath background) region.
In this method, all regions are placed a the same distance (radius) of the
center of the field-of-view, assuming that the acceptance is indeed
axi-symmetrical (see figure).
This method works only for situation in which there is only
one point-like source in the field-of-view (and not other sources to be
avoided).

.. figure:: figures/OnOff_method.png
   :alt: alternate text
   :scale: 50 %
   :align: center

   Illustration of the on-off method

`SoHAPPy` simulates the number of counts in the `on-off` regions within the
field of view, and cumulate these counts along time, i.e. the slices of a time
slot (respectivley handled by the :class:`timeslot.Slot` and :class:`obs_slice.Slice` classes).
This is done using some `gammapy` objects (:obj:`Dataset`). It is repeated a
certain number of times (the Monte Carlo trials). Each time the counts are
fluctuated.

The counts are passed to a :class:`analyze.Analysis` instance to compute the running significances along slices, using the Li&Ma formula. From these running
significance the mean time for a 3 sigma and 5 sigma detection is obtained, as
well as the mean time at which the significance is maximal, and
the mean maximal significance.

The output is stored inside a text file for each source studied. This file can
then be analysed through dedicated `SoHAPPy` tools found in the `analysis`
subfolder.

.. _Simulation:

Simulation of the instrument response
=====================================

The Monte Carlo class (:class:`mcsim.MonteCarlo` in :obj:`mcsim.py`) is the
deepest class of `SoHAPPy`. It runs simulations over a time slot and is
associated to a site or several sites at once. The simulation uses
parameters that have been determined independently and are store in the
:obj:`mcsim_config.py` module.

The counts are obtained from a python list of `Gammapy` :obj:`dataset` (not a
`Gammapy` :obj:`datasets` object) obtained from each time slice, using the
attached IRF, and that are fluctuated along the Monte Carlo trials.

There are two main ways to use the :func:`mcsim.MonteCarlo.run` methods:

* Compute only once the number of signal and background for each time slice.
  This is in principle done without fluctuation and allows getting rapidly
  results on a large population. In particular the 5 sigma detected
  sub-population is fairly unbiased compared to the one obtained with many
  iterations since it usually concerns luminous objects with large number of
  signal events (however the number of background event can be small enough to
  marginally change the result. This can be checked more quantitavely).
  This kind of simulation can be saved to disk for further use and in
  particular spectral analysis (see `here <spectral_analysis.html>`_);

* The normal use is to generated many trials -of the order of 100- fluctuating
  the signal and background at each time slice of the :class:`timeslot.Slot`,
  in order to address, in particular, objects for which the 5 sigma detection is
  not granted depending on the count fluctuations. This use is for studying the
  response of the instrument to populations.

..
    .. automodapi:: mcsim
       :include-all-objects:
       :no-inheritance-diagram:

.. _mcsim:

The `mcsim` module
------------------
This module contains only the :class:`mcsim.montecarlo` class.


.. autoclass:: mcsim.MonteCarlo
    :special-members: __init__
    :members:

.. _mcsim_config:

.. automodapi:: mcsim_config
   :include-all-objects:
   :no-inheritance-diagram:

.. _Analysis:

Analysis of the simulated data
==============================

The functions of the :class:`analyze.Analysis` class handle the list over time
of the on and off-counts provided by the :class:`mcsim.MonteCarlo` class. It
gets from the Monte Carlo trials:

    * the list of maximum significances reached and the corresponding
      observation times and excess and background counts;
    * the list of times at which 3 and 5 sigma were reached and the associated
      excess and background count numbers.
    * a list of slices for which warning were issued for further debugging
      (not intended to the normal user).

and from this it obtains the mean times at which 3 and 5 sigma detection were
obatined as well as the mean time at which the significance was found maximal
together with the maximal significance mean value.

.. _analyze:

The `analyze` module
--------------------
This module contains only the :class:`analyze.Analysis` class.

.. autoclass:: analyze.Analysis
    :special-members: __init__
    :members:


.. _Dataset:

Dataset tools
=============

.. automodapi:: dataset_tools
   :no-inheritance-diagram:
   :include-all-objects: