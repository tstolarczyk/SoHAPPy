# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:24:09 2021

@author: Stolar

A standalone function to compute and dumpp visibilities to disk
"""
import os
import astropy.units as u
from   utilities import Log

os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'

from SoHAPPy import get_grb_fromfile
import grb_plot as gplt

#import ana_config as cf # Steering parameters
import time


# Overseed GRB parameters
altmin     = 24*u.deg
altmoon    = 90*u.deg
moondist   = 0*u.deg  # Moon minimal distance
moonlight  = 1        # Moon maximum brigthness
depth      = 3*u.day # Maximum duration after trigger to compute the visibility.
skip       = 0
vis_folder = "."

dbg  = 0
show = 0

ngrb    = 2 # 250
ifirst  = 1

readvis = False # If True, recovery testing

log = Log(name  = "visibility_write.log",talk=True)

##--- Lets' start
start = time.time() # Starts chronometer


print(">>>> Output folder : ",vis_folder)
os.makedirs(vis_folder, exist_ok=True)

# GRB list to be analysed
if type(ifirst)!=list:
    grblist = list(range(ifirst,ifirst+ngrb))
else:
    grblist = ifirst

# Loop over GRB list
for i in grblist:

    print(" Processing ",i)
    grb = get_grb_fromfile(i,log=log)

    for loc in ["North","South"]:
        grb.vis[loc].compute(altmin  = altmin,
                                  altmoon = altmoon,
                                  depth   = depth,
                                  skip    = skip,
                                  debug=False)
    if dbg > 0:
        print(grb)
        grb.vis["North"].print(log=log)
        grb.vis["South"].print(log=log)

    save_vis = True
    if save_vis:
        for loc in ["North","South"]:
            grb.vis[loc].write(folder=vis_folder,debug=False)

    if (show > 0):
        gplt.visibility_plot(grb, loc="North")
        gplt.visibility_plot(grb, loc="South")


stop = time.time() # Starts chronometer
log.prt("Completed in {:8.2f} s ({:4.2f} s per source)".format(stop-start,(stop-start)/ngrb))

read_vis = True
if (read_vis):
    from visibility import Visibility
    print("\n .................... RECOVERING ...............\n")
    for loc in ["North","South"]:

        grb.vis[loc]= Visibility.read(grb.name+"_"+loc+"_vis.bin",debug=True)
        grb.vis[loc].print(log)
        if (show >0): gplt.visibility_plot(grb, loc=loc)

