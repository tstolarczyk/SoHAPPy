# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:24:09 2021

@author: Stolar

A standalone function to compute and dump visibilities to disk.
See visibility.yaml fro possible parameters.
"""
import os, sys
import time
from pathlib import Path
import yaml
from yaml.loader import SafeLoader
import grb_plot as gplt

from grb import GammaRayBurst
import visibility as vis
from   utilities import Log

# Transform warnings into errors - useful to find who is guilty !
import warnings
#warnings.filterwarnings('error')
warnings.filterwarnings('ignore')

###---------------------------
### change your conditons here
###---------------------------
def_vis     = "visibility.yaml" # Visibility parameter fils
conditions  = "strictmoonveto" # One condition in the default visibility file
vis_folder  = conditions+"_DC-2028_01_01_000000-2034_12_31_235959"
# vis_archive = "../input/visibility/short_vis_24_strictmoonveto"
# vis_archive = conditions

###---------------------------
### Your actions
###---------------------------
save_vis   = True # Save to disk in vis_folder
read_vis   = False # Read from disk in vis_folder
debug      = False

###---------------------------
### Input GRB data
###---------------------------
ifirst  = 1
ngrb    = 2 # 250
# grb_folder = "../input/lightcurves/SHORT_FITS/"
grb_folder = "D:/CTAA/SoHAPPy/input/lightcurves/LONG_FITS/"
triggers = "D:/CTAA/SoHAPPy/input/visibility/long/Trigger_1000-2028_01_01_000000-2034_12_31_235959.yaml"

if type(ifirst)!=list: grblist = list(range(ifirst,ifirst+ngrb))
else: grblist = ifirst

###---------------------------
### Get new trigger dates if requested
###---------------------------
from trigger_dates import get_from_yaml
data = get_from_yaml(triggers)

###---------------------------
### Let's go
###---------------------------
start = time.time() # Starts chronometer

# Prepare writing on disk
if save_vis:
    os.makedirs(vis_folder, exist_ok=True) # Create output folder
    log = Log(name  = vis_folder+"/visibility_write.log",talk=True)
    log.prt("Writing visibility to Output folder : {} ".format(vis_folder))
    
    # Read the visibility parameters 
    log.prt(">>> Read Visibility configuration from {}".format(def_vis))    
    with open(def_vis) as f:
        visdict  = yaml.load(f, Loader=SafeLoader)
        if conditions in visdict.keys(): visibility = visdict[conditions]
        else: sys.exit("Conditions are not found in {}".fomart(def_vis))    
else:
    visibility = None
    
# Loop on files
for i, item in enumerate(grblist):
    fname = "Event"+str(item)+".fits.gz"
    if save_vis:
        grb = GammaRayBurst.from_fits(Path(grb_folder,fname),
                                  ebl     = None,
                                  prompt  = False, 
                                  dt      = data[fname],
                                  dt_abs  = True,
                                  vis     = visibility)

        for loc in ["North","South"]:
    
            if save_vis:
                log.prt(" Saving visibility for GRB {} in {}".format(item,loc))
                grb.vis[loc].write(folder=vis_folder,debug=False)
                if debug:
                     grb.vis[loc].print()
                     gplt.visibility_plot(grb, loc=loc)               
    elif read_vis: 
        log=Log()
        for loc in ["North","South"]:
            visname = Path(vis_archive,filename.stem+"_"+loc+"_vis.bin")
            print(" Reading visibility from {}".format(visname))
            myvis = vis.Visibility.read(visname,debug=True)
            if debug:
                 myvis.print()
                 # Requires GRB
                 # gplt.visibility_plot(grb, loc=loc)    

# Goodbye !
stop = time.time() # Starts chronometer
log.prt("Completed in {:8.2f} s ({:4.2f} s per source)"
        .format(stop-start,(stop-start)/ngrb))