Time slices and slots
=====================

Time slots are defined through the :class:`timeslot.Slot` class in
:file:`timeslot.py` and contain a collection of time slices handled by
the :class:`obs_slice.Slice` class in :file:`obs_slice.py`. The time slots
reflect the GRB time slices after taking into account the various delays and
visibility periods.

Time Slot
---------
A so-called time slot is handled by the :class:`timeslot.Slot` class in
:file:`timeslot.py`. It is the key element of the `SoHAPPy` analysis.
A time slot refers to a GRB object, a time delay, the start and end time of
the slot, the site(s) of the observation and a collection of slices as defined
in the class :class:`obs_slice.Slice`.

The native slot is composed ot the GRB initial time slices. After application
of delays or visibilities, it contains the effective time slices on which the
analysis can be done. It can handle observations on one or several sites.

The :class:`timeslot.Slot` class holds the tools to modifiy the initial slot
into the slot attached to the GRB at a given site (function
:func:`timeslot.Slot.apply_visibility`for a given site and
:func:`timeslot.Slot.both_sites` to have both site taken into account).
Slots can be merged with the :func:`timeslot.Slot.merge` function.

Initial naked slots that handles the time slices are dressed with physical
information with :func:`timeslot.Slot.dress`: this has for effect to
associate the correct (or best) IRF and spectrum to each of the slices.

At that stage the time slot can be passed to the :class:`mcsim.MonteCarlo`
object.

.. automodapi:: timeslot
   :include-all-objects:
   :no-inheritance-diagram:

Time slices
-----------
A time slice, handled by the :class:`obs_slice.Slice` in :file:`obs_slice.py`,
carries a start an stop time, an observation time, one or several IRFs, a
spectrum identifier from the GRB class spectra, a site tagging (A slice can
belong to more than one site). The observation time can be set at the start,
stop or mean time of the slice (more options could be implemented).

The default is to set the observation at the "end" of the slice (can be changed
by the user in the `configuration file <configuration.rst>`). The associated
IRFs are obtained according to the altitude (zenith) and slice duration at the beginning of the slice where the flux is expected to be higher. The flux is
the flux from the GRB which is the closest in time to the observation time.

.. automodapi:: obs_slice
   :include-all-objects:
   :no-inheritance-diagram: