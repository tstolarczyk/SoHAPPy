# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 17:19:19 2024

@author: Stolar
"""
import os
import pickle
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from gammapy.utils.random import get_random_state
from niceplot import MyLabel
from dataset_tools import onoff_from_simulation, sigmax

sys.path.append("../")
sys.path.append("../../../SoHAPPy")


# ##---------------------------------------------------------------------------
def datasets_from_binary(filepath=None):

    if filepath is None:
        if "INFILE" in os.environ.keys():
            filepath = Path(os.environ["INFILE"])
        else:
            sys.exit("No file given and no INFILE environement variable")

    if not filepath.exists():
        sys.exit(f"File {filepath:} not found")

    # Get the MC class instance
    infile = open(filepath, "rb")
    mc = pickle.load(infile)
    infile.close()

    # Re-create the on-off dataset list (datasets) from the simulation
    # Note that this is another MC realisation"
    dsets = onoff_from_simulation(mc, random_state=get_random_state(2021),
                                  debug=False)

    print(f"Number of sets in datasets : {len(dsets)}")

    return dsets, mc.slot.grb


# ##---------------------------------------------------------------------------
def check_max_sigma(mc, sigmx):

    siglist = []
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    for i in range(100):
        print(i, " ", end="")
        dset_test = onoff_from_simulation(mc, debug=False,
                                          random_state=get_random_state(2021))
        siglist.append(sigmax(dset_test))
    print()
    ax.hist(siglist, label=MyLabel(siglist, label="MC Siginificances"))
    ax.axvline(sigmx, ls="--", color="red")
    ax.legend()
