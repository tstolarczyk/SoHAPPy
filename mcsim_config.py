# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:58:16 2020

@author: Stolar
"""
import astropy.units as u

<<<<<<< HEAD
# These are the parameters specific to the Monte Carlo simulation
=======
# This are the parameters specific to the Monte Carlo simulation
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76

det_level  = 0.9 # Fraction declaring a detection above the 3,5 sig. threshold
alpha      = 0.2 # 5 zones (1/0.2) to define the off region for B estimate
n_lima_min = 10 # Count number below which Li&Ma cannot be trusted anymore

fov        = 5.*u.deg    # Full Field-of-view
binsize    = 0.125*u.deg #  bin size
# fov      = 2.5*u.deg # Quicker dirty tests
<<<<<<< HEAD
# binsize  = 0.5*u.deg # Quicker dirty tests
=======
# binsize  = 0.5*u.deg # Quicker dirty tests

from   astropy.coordinates   import EarthLocation
# The following information can be obtained from EarthLocation.of_site
# but this requires to have an internet connection
# Get all available site names :
# sitelist = EarthLocation.get_site_names()

# Get coordinattes of Paranal and LaPalma
# Note : with the Gammapy-0.12 envirnment,this requires to 
#        downgrade numpy to version 1.16.2 (instead on >1.18)
#        See discussion here :
#        https://github.com/astropy/astropy/issues/9306

#xyz_north = EarthLocation.of_site('Roque de los Muchachos')
xyz_north = EarthLocation.from_geocentric( 5327448.9957829,
                                          -1718665.73869569,
                                           3051566.90295403,
                                           unit="m")

# xyz_south = EarthLocation.of_site('Paranal Observatory')
xyz_south = EarthLocation.from_geocentric( 1946404.34103884,
                                          -5467644.29079852,
                                          -2642728.20144425,
                                          unit="m")

pos_site ={"North":xyz_north, "South":xyz_south} 
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
