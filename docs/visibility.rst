Visibilities
============

The time slices and slots are derived from the original *fits* file slices after applying the visibility. 
The visibility is defined by the :class:`visibility.Visibility` class in `visibility.py`. It is either read from original provided data (:func:`visibility.Visibility.from_fits`) or pre-computed values (:func:`visibility.Visibility.read`), or computed on the fly (:func:`visibility.Visibility.compute`).
Once computed it can be save on disk (:func:`visibility.Visibility.write`) to be further re-used.

The visibility consists in time periods that will be applied to the original object time slice to create new time slices on which the observation will be simulated. *It does not account for the delay due to the alert or the telescope slewing* so that it can be re-used to studies the effect of these delays.

A function (:func:`visibility.Visibility.check`) allows comparing two visibilities obtained in different ways and give the compatibility of the two within a time difference.

For each site location, the visibility it depends on:
	- the abscence of Sun (astronomical twilight, Sun below -18°);
	- the altitude of the object above the horizon (typicall values are 10°, or 24° following the CTA requirements);
	- the presence of Moon light as defined by :

		- the Moon alitude in the sky;
		- the moon phase (or illumination, between 0 and 1);
		- the angular distance of the object to the Moon.
    
All this is obtained from astroplan functions.

.. automodapi:: visibility
   :include-all-objects:
   :no-inheritance-diagram: