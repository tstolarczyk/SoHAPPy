Visibilities
============
The time slices and slots are derived from the original *fits* file slices after applying the visibility. 

For each site location, the visibility depends on:
	- the abscence of Sun (astronomical twilight, Sun below -18°);
	- the altitude of the object above the horizon (typicall values are 10°, or 24° following the CTA requirements);
	- the presence of Moon light as defined by :

		- the Moon alitude in the sky;
		- the moon phase (or illumination, between 0 and 1);
		- the angular distance of the object to the Moon.

All this is obtained from astroplan functions.

The visibility is defined by the :class:`visibility.Visibility` class in `visibility.py`. It is either read from original provided data (:func:`visibility.Visibility.from_fits`) or pre-computed values (:func:`visibility.Visibility.read`), or computed on the fly (:func:`visibility.Visibility.compute`). In that case the visibility parameters are read from a `visibility.yaml`file in the SoHAPPy local folder (see below)
Once computed it can be save on disk (:func:`visibility.Visibility.write`) to be further re-used.

The visibility consists in time periods that will be applied to the original object time slice to create new time slices on which the observation will be simulated. *It does not account for the delay due to the alert or the telescope slewing* so that it can be re-used to studies the effect of these delays.

A function (:func:`visibility.Visibility.check`) allows comparing two visibilities obtained in different ways and give the compatibility of the two within a time difference.

Configuration file
------------------

The visibility parameters are stored in dedicated `visibility.yaml`. They are used when the SoHAPPy configuration file (the default is `config.yaml`) requires that the visibility are recomputed on the fly. In that case the visibility entry takes on of these parameters:

* strictmoonveto : the source can be considered visible only if the Moon is below the horizon (alitude less than 0.25°).
* nomoonveto : the source can be considered visible even if the Moon is above the horizon whatever its distance to the soource or it phase.
* moonlight : the source canbe considered visible f the moon phase is small enough andthe distance to the source large enough.
* alwaysvisible : the source is visible as soon as it is abovethe horizon (no altitude lower limit) and the Moon is ignored (Unrealistic, used to explore the parameter space coverage).


The visibility parameters are the following:


.. tabularcolumns:: |l|c|p{5cm}|

+-----------------------+------------------------+---------------------------------------------+
| variable              | Default for moonlight  | What is it ?                                |
+=======================+========================+=============================================+
| where                 | "CTA"                  | Observatory position keyword                |
+-----------------------+------------------------+---------------------------------------------+
| depth                 | 3                      | Max. number of nights                       |
+-----------------------+------------------------+---------------------------------------------+
| skip                  | 0                      | Number of first nights to be skipped        |
+-----------------------+------------------------+---------------------------------------------+
| altmin                | 24 deg                 | Minimum altitude                            |
+-----------------------+------------------------+---------------------------------------------+
| altmoon               | -0.25 deg   (horizon)  | Moon maximum altitude                       |
+-----------------------+------------------------+---------------------------------------------+
| moondist              | 30 deg                 | Moon minimal distance                       |
+-----------------------+------------------------+---------------------------------------------+
| moonlight             | 0.6                    | Moon maximum brigthness                     |
+-----------------------+------------------------+---------------------------------------------+



.. automodapi:: visibility
   :include-all-objects:
   :no-inheritance-diagram: