# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:58:16 2020

@author: Stolar

This modules contain the parameters of the anmysis and simulation handled
by mcsim.py;It is also used by mcsim_plot and mcsim_res.
- det_level (0.9): Fraction for declaring a detection above the 3 or 5 sigma thresholds
- alpha (0.2) : related to the number of off N zones in the on-off analysis, alpha = 1/N
- on_size : the on (and off-) region radius
- offset : the on and off region distance to the center of the field of view (Should be greater than the on-size radius)
- containment : the fraction of the background area corresponding to on/off regions
- n_lima_min (10) : Count number below which Li & Ma formula cannot be trusted anymore
- fov (5.*u.deg) : Full Field-of-view for the sky map
- binsize (0.125*u.deg) : square pixel length for the sky map

"""
import astropy.units as u

det_level   = 0.9 # Fraction declaring a detection above the 3,5 sig. threshold
nLiMamin    = 10  # Count number below which Li&Ma cannot be trusted anymore
alpha       = 0.2 # 5 zones (1/0.2) to define the off region for B estimate

containment = 0.68

on_size     = { "FullArray": 0.4*u.deg,
                "LST"      : 0.4*u.deg,
                "MST"      : 0.2*u.deg } # On region size

offset      = { "FullArray": 0.75*u.deg,
                "LST"      : 0.75*u.deg,
                "MST"      : 0.4*u.deg} # Offset, so that regions do not overlap
# on_size     = 1.3*u.deg # On region size
# offset      = 1.3*u.deg # Offset

fov        = 5.*u.deg    # Full Field-of-view
binsize    = 0.125*u.deg #  bin size if building a map
# fov      = 2.5*u.deg # Quicker dirty tests
# binsize  = 0.5*u.deg # Quicker dirty tests