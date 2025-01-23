# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:50:34 2023

This script shows how to analyse the time and energy spectra of a simulated
GRB.
Note that in some cases the available data are unsufficient to get a
spectrum fitted, in particular when observation starts early and the time
slices are very short. In that case it is recommened to stack the slices within
a minimal duration.

@author: Stolar
"""

import sys
import os
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import astropy.units as u

from gammapy.utils.random import get_random_state

from dataset_plot import windows, panels
from dataset_counts import excess_counts, excess_versus_time, lightcurve, \
    residuals
from dataset_tools import onoff_from_simulation, check_datasets, sigmax

from niceprint import heading
from niceplot import stamp

sys.path.append("../")
sys.path.append("../../../SoHAPPy")

###############################################################################
# Get data from a GRB simulated file
# Use the INFILE environment variable if it exists

if "INFILE" in os.environ.keys():
    file = Path(os.environ["INFILE"])
else:
    base = Path(r"D:\CTAO\SoHAPPy\HISTORICAL\180720B\Night_1")

    GRB_id = "180720B"
    site = "South"
    file = Path(base, GRB_id+"-"+site+"_sim.bin")

# This is required to have the EBL models read from gammapy
os.environ['GAMMAPY_DATA'] =\
 str(Path(Path(__file__).absolute().parent.parent.parent, "data"))

if not file.exists():
    sys.exit(f"File {file:} not found")

# Get the MC class instance
infile = open(file, "rb")
mc = pickle.load(infile)
infile.close()

# GRB characteristics
heading("GRB")
print(mc.slot.grb)
# mc.slot.grb.plot()

# Re-create the on-off dataset list (datasets) from the simulation
# Note that this is another MC realisation"
deeper = True  # Display Energy bins in each set
heading("DATASETS")
dset_init = onoff_from_simulation(mc,
                                  random_state=get_random_state(2021),
                                  debug=False)

dsets = dset_init.copy()
print(f"Number of sets in datasets : {len(dsets)}")
check_datasets(dsets, masked=False, deeper=deeper, header=True)

# Check significance
heading("Max significance stability")

sigmx = sigmax(dset_init)
print(f" Significance from this realization = {sigmx:5.1f}")

# Mean significance from many realisations
# fromutils import check_max_sigma
# check_max_sigma(mc, sigmx)

# Display the observation start times and windows
heading("Observation windows")
windows(dsets)
stamp(file.stem)

# ## ------------------
# ##   Excess counts
# ## ------------------

# Excess count numbers
# Min and max number of counts can be changed, as well as Emin and Emax, see
# excess_count function documentation
heading(" Excess count numbers")
merging = False
fig = panels(dsets,
             excess_counts,
             xsize=6, yscale="linear",
             nmaxcol=3,
             stacked=merging,
             tref=mc.slot.grb.t_trig,
             fixrange=True)
stamp(file.stem, axis=fig)
plt.tight_layout(h_pad=0, w_pad=0)

# Excess counts residuals - uses the Gammapy function - to betuned
heading(" Excess count residuals")
fig = panels(dsets, residuals, ysize=5, nmaxcol=4, yscale="linear")
stamp(file.stem, axis=fig)
plt.tight_layout(h_pad=0, w_pad=0)

# Excess versus time
heading(" Excess versus time")
xscale = "linear"  # log, linear
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15, 8))
excess_versus_time(dsets, ax=ax[0], debug=True, rate=False,
                   xscale=xscale, yscale="log", unmasked=False)
excess_versus_time(dsets, ax=ax[1], debug=True, rate=True,
                   xscale=xscale, yscale="log")
stamp(file.stem, axis=fig)

# Light curve simulation, compared to prediction/theory
heading(" Light curve")

dtbin = 225*u.s
xscale = "linear"
yscale = "linear"

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 6))
ax = lightcurve(dsets,
                date_ref=mc.slot.grb.t_trig, tag="Excess", style="line",
                binwidth=dtbin,
                ax=ax, xscale=xscale, yscale=yscale, marker=".")
lightcurve(dsets,
           date_ref=mc.slot.grb.t_trig,
           tag="Prediction", style="bar",
           color="tab:blue", alpha=0.2,
           binwidth=dtbin,
           ax=ax, xscale=xscale, yscale=yscale)
stamp(file.stem)
plt.tight_layout()
