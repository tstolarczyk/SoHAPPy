# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:24:09 2021

@author: Stolar

A standalone function to compute and dump visibilities to disk
Suggested parameters

 
* To maximise the visibility use the following values:
  - altmin    =  0  deg   # Ensure that the source is always above horizon
  - altmoon   =  90 deg   # Ensure that the moon never vetoes the visibility
  - moondist  =  0  deg   # The Moon distance do not veto the visibility
  - moonlight =  1.0      # The Moon brightness is not a limitation

* To minimize the visibility (i.e. veto as soon as the moon is above horizon)
  - altmin    =  24 deg   # Ensure that the source is always above horizon
  - altmoon   = -0.25 deg # Moon above horizon
  - moondist  = 180 deg   # The Moon vetoes even if far away
  - moonlight =  0.0      # The Moon vetoes even if new Moon

"""
import os
import time
import astropy.units as u

from SoHAPPy import get_grb_fromfile
import grb_plot as gplt
import visibility as vis
from   utilities import Log

# GRB to be processed
ifirst  = 1
ngrb    = 1000 # 250

conditions = "maximum"
save_vis   = True # Save to disk in vis_folder
read_vis   = False # Read from disk in vis_folder

if conditions == "normal":
    vis_folder = "../output/vis_24_moonlight/"
    altmin    =    24*u.deg # CTA requirement horizon
    altmoon   = -0.25*u.deg # Moon above horizon
    moondist  =    30*u.deg # Acceptable Moon distance
    moonlight =   0.6       # Acceptable Moon brightness    
elif conditions == "minimum":
    vis_folder = "../output/vis_24_nomoonveto/"
    altmin    =   24*u.deg  # CTA requirement horizon
    altmoon   =  90*u. deg  # Ensure that the moon never vetoes the visibility
    moondist  =   0*u.deg   # The Moon distance do not veto the visibility
    moonlight = 1.0         # The Moon brightness is not a limitation    
elif conditions == "maximum":
    vis_folder = "../output/vis_24_strictmoonveto/"
    altmin    =    24*u.deg # CTA requirement horizon
    altmoon   = -0.25*u.deg # Moon above horizon
    moondist  =   180*u.deg # The Moon vetoes even if far away
    moonlight =   0.0       # The Moon vetoes even if new Moon
    
depth = 3        # (3) Maximum number of nights to compute the visibility.
skip  = 0        # (0) Number of first nights to be skipped
dbg   = 0
show  = 0

# GRB list to be analysed
if type(ifirst)!=list:
    grblist = list(range(ifirst,ifirst+ngrb))
else:
    grblist = ifirst

# Compute and save visibility
if save_vis:
    start = time.time() # Starts chronometer
    os.makedirs(vis_folder, exist_ok=True) # Create output folder
    print("Writing visibility to Output folder : ",vis_folder)
    log        = Log(name  = vis_folder+"visibility_write.log",talk=True)

    # Loop over GRB list
    for i in grblist:
        print(" Processing ",i)
        grb = get_grb_fromfile(i, grb_folder = "../input/lightcurves/" , log=log)
    
        for loc in ["North","South"]:
            
             grb.vis[loc] = vis.Visibility.compute(grb,
                                                   loc,
                                                   altmin    = altmin,
                                                   altmoon   = altmoon,
                                                   moondist  = moondist,
                                                   moonlight = moonlight,
                                                   depth     = depth,
                                                   skip      = skip,
                                                   debug     = bool(dbg>2))         
             grb.vis[loc].write(folder=vis_folder,debug=False)


        if dbg > 0:
            print(grb)
            grb.vis["North"].print(log=log)
            grb.vis["South"].print(log=log)

    
        if (show > 0):
            gplt.visibility_plot(grb, loc="North")
            gplt.visibility_plot(grb, loc="South")

    stop = time.time() # Starts chronometer
    log.prt("Completed in {:8.2f} s ({:4.2f} s per source)".format(stop-start,(stop-start)/ngrb))

if (read_vis):
    print("recovering visibility from input folder : ",vis_folder)
    for i in grblist:

        for loc in ["North","South"]:
            grb = get_grb_fromfile(i, grb_folder = "../input/lightcurves/" , log=log)

            grb.vis[loc]= vis.Visibility.read(vis_folder+grb.name+"_"+loc+"_vis.bin",debug=True)
            grb.vis[loc].print(log)
            if (show >0): gplt.visibility_plot(grb, loc=loc)
