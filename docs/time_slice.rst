Time slices and slots
=====================

Time Slot
---------
A so-called time slot, the class::Slot:: in timeslot.py, is the key element of the SoHAPPy analysis.
A Slot contains a GRB object, a time delay, the start and end time of the slot, the site(s) of the observation and a collection of Slices.


The native slot is composed ot the GRB initial time slices. After application of delays or visibilities, it contains the effective time slices on which the analysis can be done. It can handle observations on one or several sites.

The Slot class hold the tools to modifiy the initial slot into the slot attached to the GRB at a given site (function apply_visibility and both_sites). Slots can be merged with the merge function.

Initial naked slots that handles the time slices are dressed with physical information : this has for effect to associated the correct (or best) IRF and spectrum to each of the slices.

.. automodapi:: timeslot
   :include-all-objects:
   :no-inheritance-diagram:

Time slices
-----------
A time slice, object Slice in obs_sliece.py, carries a start an stop time, an observation time, one or several IRFs, a spectrum identifier from the GRB class spectra, a site tagging (A slice can belong to more than one site).
The observation time can be set at the start, stop or mean time of the slice (more options could be implemented).

The default is to set the observation at the "end" of the slice (can be changed by user).
The associated IRFs are obtained according to the altitude (zenith) and slice duration at the beginning of the slice where the flux is expected to be higher (see recommendation on how to choose the IRF).
The flux is the flux from the GRB which is the closest in time to the observation time.

.. automodapi:: obs_slice
   :include-all-objects:
   :no-inheritance-diagram: