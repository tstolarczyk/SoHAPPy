Visibilities
============

The time slice and slot are derived from the original fit file slices after applying the visibilty. 
The visibility is defined by the Visibilty class in visibility.py. It is either read from original provided data (from_fits method) or pre-computed values (read method), or computed on the fly (compute method).
Once computed it can be save on disk (write method) to be further read.
The visibility consists in time periods that will be applied to the original object time slice to create new time slices on which the observation will be simulated.
A method (check) allows comparing two visibilities obtained in different ways and give the compatibility of the two within a time difference.
A specific visibility_write.py module can be used to compute and save on disks the visibilities for a large population of objects.

For each site location, it depends on:
- the abscence of Sun (astronomical twilight, Sun below -18°)
- the alitude of the object above the horizon (typicall values are 10°, or 24° following the CTA requirements).
- the presence of Moon light as defined by :

    - the Moon alitude in the sky;
    - the moon phase (or illumination, between 0 and 1);
    - the distance of the object to the Moon.
    
All this is obtained from astroplan functions.


.. automodapi:: visibility
   :include-all-objects:
   :no-inheritance-diagram: