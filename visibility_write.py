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
conditions = "moonlight"
def_vis    = "visibility.yaml"
vis_folder = conditions

###---------------------------
### Your actions
###---------------------------
save_vis   = True # Save to disk in vis_folder
read_vis   = False # Read from disk in vis_folder
debug      = False
log = Log(name  = vis_folder+"visibility_write.log",talk=True)

###---------------------------
### Input GRB data
###---------------------------
ifirst  = 1
ngrb    = 2 # 250
grb_folder = "../input/lightcurves/"

if type(ifirst)!=list: grblist = list(range(ifirst,ifirst+ngrb))
else: grblist = ifirst

###---------------------------
### Let's go
###---------------------------
start = time.time() # Starts chronometer

# Prepare writing on disk
if save_vis:
    os.makedirs(vis_folder, exist_ok=True) # Create output folder
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
for item in grblist:
    filename = Path(grb_folder,"LONG_FITS","Event"+str(item)+".fits")

    grb = GammaRayBurst.from_fits(filename,
                                  ebl     = None,
                                  prompt  = False, 
                                  vis     = visibility)

    for loc in ["North","South"]:

        if save_vis:
            log.prt(" Saving visibility for GRB {} in {}".format(item,loc))
            grb.vis[loc].write(folder=vis_folder,debug=False)
            if debug:
                 grb.vis[loc].print()
                 gplt.visibility_plot(grb, loc=loc)               
        
        if read_vis: 
            filename = Path(vis_folder,grb.name+"_"+loc+"_vis.bin")
            log.prt(" Reading visibility from {}".format(filename))
            grb.vis[loc]= vis.Visibility.read(filename,debug=True)
            if debug:
                 grb.vis[loc].print()
                 gplt.visibility_plot(grb, loc=loc)    

# Goodbye !
stop = time.time() # Starts chronometer
log.prt("Completed in {:8.2f} s ({:4.2f} s per source)"
        .format(stop-start,(stop-start)/ngrb))