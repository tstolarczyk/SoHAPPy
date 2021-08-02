# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import astropy.units as u
"""
This module hold the visibility computation parameters.

"""

depth     = 3        # (3) Maximum number of nights to compute the visibility.
skip      = 0        # (0) Number of first nights to be skipped

# Suggested best default
altmin    =  24*u.deg   # Minimum altitude (original default is 10 degrees)
altmoon   = -0.25*u.deg # Moon maximum altitude (-0.25 u.deg for horizon)
moondist  =  30*u.deg   # Moon minimal distance
moonlight =  0.6        # Moon maximum brigthness

# In tests, to maximise the visibility use the following values:
# altmin    =  0*u.deg  # ensure that the source is always above horizon
# altmoon   =  90*u.deg # Ensure that he moon never vetoes the visibility
# moondist  =  0*u.deg  # The Moon distance do not veto the visibility
# moonlight =  1.0      # The Moon brightness is not a limitation


def print(log=None):
    log.prt(" Vis. computed up to     : {} night(s)".format(depth))
    log.prt(" Skip up to              : {} night(s)".format(skip))
    log.prt(" Minimum altitude        : {}".format(altmin))
    log.prt(" Moon max. altitude      : {}".format(altmoon))
    log.prt(" Moon min. distance      : {}".format(moondist))
    log.prt(" Moon max. brightness    : {}".format(moonlight))

    return

if __name__ == "__main__":
    from utilities import Log
    log = Log("dummy.log")
    print(log)
