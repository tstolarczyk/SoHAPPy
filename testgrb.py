# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 17:38:42 2019

@author: Stolar
"""

from grb            import GammaRayBurst
import grbplot      as gplt
from   gammapy.spectrum.models import Absorption
import matplotlib.pyplot   as plt
from   pathlib import Path

import ana_config as cf

i=1
eblabs = Absorption.read_builtin(cf.EBLmodel)
loc = Path(cf.grb_folder + "/Event"+str(i)+".fits")
grb = GammaRayBurst.read_new(loc, ebl = eblabs)
gplt.spectra(grb,"Packed")
# print(grb.site)
# print(grb.pointing)

# print("Eval =",grb.Eval)
# print("tval =",grb.tval)

# print(" tval length ",len(grb.tval))

# print("interva length = ",len(grb.time_interval))

# print("flux length = ",(grb.fluxval.shape))