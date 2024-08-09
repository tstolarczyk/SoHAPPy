# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 10:26:34 2022

This module contains information on various observatory locations.

@author: Stolar
"""
import astropy.units as u
from   astropy.coordinates   import Angle, EarthLocation

__all__ = ["xyz"]

### Site positions
xyz = {
   "CTAO":
         { "North": EarthLocation.from_geocentric( float(5327448.9957829),
                                                   float(-1718665.73869569),
                                                   float(3051566.90295403),
                                                   unit="m"),
           "South": EarthLocation.from_geocentric( float(1946404.34103884),
                                                   float(-5467644.29079852),
                                                   float(-2642728.20144425),
                                                   unit="m")
           # Requires to have an internet connection
           # Get all possible sites with : EarthLocation.get_site_names()
           # "North": EarthLocation.of_site('Roque de los Muchachos'),
           # "South": EarthLocation.of_site('Paranal Observatory')
           # Values taken from Maria Grazia Bernardini (July 28, 2020)
           # "North": EarthLocation.from_geocentric( 5327285.09211954,
           #                                        -1718777.11250295,
           #                                         3051786.7327476,
           #                                         unit="m"),
           # "South": EarthLocation.from_geocentric(1946635.7979987,
           #                                       -5467633.94561753,
           #                                       -2642498.5212285,
           #                                        unit="m")
           # Alternative :
           # "North": EarthLocation.from_geodetic('342.1184',
           #                                       '28.7606',2326.* u.meter)
           # "South": EarthLocation.from_geodetic('289.5972',
           #                                      '-24.6253',2635.* u.meter)
         },

  # from default values in
  # https://www.mpi-hd.mpg.de/hfm/HESS/pages/home/visibility/
  # Warning : use CTA-North for virtual HESS North
  "HESS":
        { "North": EarthLocation.from_geocentric( float(5327448.9957829),
                                                  float(-1718665.73869569),
                                                  float(3051566.90295403),
                                                  unit="m"),
          "South": EarthLocation.from_geodetic(lat=Angle('-23d16m18.0s'),
                                      lon=Angle('16d30m0.0s'),
                                      height=1800*u.m)
        }
 }
