# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:58:16 2020

@author: Stolar
"""
import astropy.units as u

# These are the parameters specific to the Monte Carlo simulation

det_level  = 0.9 # Fraction declaring a detection above the 3,5 sig. threshold
alpha      = 0.2 # 5 zones (1/0.2) to define the off region for B estimate
n_lima_min = 10 # Count number below which Li&Ma cannot be trusted anymore

fov        = 5.*u.deg    # Full Field-of-view
binsize    = 0.125*u.deg #  bin size
# fov      = 2.5*u.deg # Quicker dirty tests
# binsize  = 0.5*u.deg # Quicker dirty tests