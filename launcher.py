# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 15:02:47 2022

@author: Stolar
"""
import os
import subprocess
datasets = [[1,100],[101,200]]

# datasets = [[1,100],[101,200],[201,300]]
# from pathlib import Path
# direxe = Path(__file__).parent.resolve() # For the directory of the script being run
direxe = "./"
# print(direxe)
# pathlib.Path().resolve() # For the current working directory:

for dset in datasets:
    first      = dset[0]
    Nsrc       = 100 #  dset[1]-dset[0]+1
    niter      = 100
    debug      = 0
    config     = "config_prod.yaml"
    outdir     = "test_"+str(dset[0])+"_"+str(Nsrc+dset[0])
    visibility = "visibility/long/vis_"+str(dset[0])+"_"+str(dset[1])+"_strictmoonveto.json"
    # trigger    = "visibility/long/Trigger_1000-2028_01_01_000000-2034_12_31_235959.yaml"
    trigger    = 0

    command  = "python"
    command += " SoHAPPy.py"
    command += " -d " + str(debug)
    command += " -c " + config
    command += " -o " + outdir
    command += " -n " + str(niter)
    command += " -N " + str(Nsrc)
    command += " -f " + str(first)
    command += " -V " + visibility
    command += " -D " + str(trigger)
    print("cmd = ",command)
    #subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    os.system(command)

