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

from gammapy.utils.random import get_random_state


from niceprint import heading
from niceplot import stamp, MyLabel

import matplotlib.pyplot as plt
import astropy.units as u

from dataset_tools import createonoff_from_simulation, check_datasets, sigmax
from dataset_plot import windows, panels
from dataset_counts import excess_counts, excess_versus_time, lightcurve, residuals
from dataset_flux import extract_spectrum, flux_versus_time

sys.path.append("../")
sys.path.append("../../../SoHAPPy")

###############################################################################

# Get data from a GRB simulated file
# Use the INFILE environment variable if it exists
if "INFILE"  in os.environ.keys():
    file = Path(os.environ["INFILE"])
else:
    base = Path("D:/CTA/SoHAPPy/output/long_1_1000/test_omega/strictmoonveto/strictmoonveto_343")
    GRB_id  = 343
    site    = "South"
    file    = Path(base,"Event"+str(GRB_id)+"-"+site+"_sim.bin")

if not file.exists():
    sys.exit(f"File {file:} not found")

# Get the MC class instance
infile  = open(file,"rb")
mc      = pickle.load(infile)
infile.close()

# GRB characteristics
heading("GRB")
print(mc.slot.grb)
mc.slot.grb.plot()

# Re-create the on-off dataset list (datasets) from the simulation
# Note that this is another MC realisation"
deeper = True # Display Energy bins in each set
heading("DATASETS")
dset_init = createonoff_from_simulation(mc,
                                        random_state=get_random_state(2021),
                                        debug=False)

dsets = dset_init.copy()
print(f"Number of sets in datasets : {len(dsets)}")
check_datasets(dsets, masked=False, deeper=deeper, header=True)

# Check siginificance
heading("Max significance stability")

sigmx = sigmax(dset_init)
print(f" Significance from this realization = {sigmx:5.1f}")

# Mean significance from many realisations
siglist = []
fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(8,8))
for i in range(100):
    print(i," ", end="")
    dset_test = createonoff_from_simulation(mc,debug=False,
                                    random_state=get_random_state(2021))
    siglist.append(sigmax(dset_test))
    ax.hist(siglist, label=MyLabel(siglist,label="MC Siginificances"))
    ax.axvline(sigmx,ls="--",color="red")
    ax.legend()

# Display the observation start times and windows
heading("Observation windows")
windows(dsets)
stamp(file.stem)

# Excess count numbers
# Mina nd max number of counts can be changed, as well as Emin and Emax, see
# excess_count function documentation
heading(" Excess count numbers")
merging = False
fig = panels(dsets,excess_counts,
                    xsize   = 6, yscale="log",
                    nmaxcol = 3,
                    stacked = merging,
                    max_margin = 1.1,
                    tref       = mc.slot.grb.t_trig,
                    fixrange   = True)
stamp(file.stem,axis=fig)
plt.tight_layout(h_pad=0, w_pad=0)

# Excess counts residuals - uses the Gammapy fucntion - to betuned
heading(" Excess count residuals")
fig = panels(dsets, residuals, ysize=5, nmaxcol=4, yscale="linear")
stamp(file.stem,axis=fig)
plt.tight_layout(h_pad=0, w_pad=0)

# Excess versus time
heading(" Excess versus time")
xscale = "log" # linear
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15,8))
excess_versus_time(dsets, ax = ax[0], debug=False, rate=False, xscale=xscale)
excess_versus_time(dsets, ax = ax[1], debug=False, rate=True,  xscale=xscale)
stamp(file.stem, axis=fig)

# Light curve simulation, compared to prediction/theory
heading(" Light curve")

tmin   = 0*u.s
tmax   = 1500*u.s
dtbin  = 10*u.s
xscale = "linear"
yscale = "linear"

fig, ax= plt.subplots(nrows=1, ncols=1, figsize=(15,6))
ax = lightcurve(dsets,tag="excess",style="line",
                          binwidth=dtbin,tmax=tmax,
                          ax=ax, xscale=xscale, yscale=yscale, marker=".")
lightcurve(dsets,tag="prediction",style="bar",color="tab:blue",alpha=0.2,
            binwidth=dtbin,tmax=tmax,
            ax=ax, xscale=xscale, yscale=yscale )
stamp(file.stem)


# Use a subset when useful, e.g. dsets[6:10]
heading(" Energy spectra")
fig = panels(dsets, func = extract_spectrum, stacked=False,
                    xsize     = 6, yscale="log",
                    nmaxcol   = 4,
                    e_ref     = 1*u.TeV,
                    tref      = mc.slot.grb.t_trig,
                    e_unit    = "TeV",
                    flux_unit = "TeV-1 cm-2 s-1",
                    flux_min  = 1.e-14,
                    flux_max  = 1.e-6,
                    tag = "",
                    debug     = True)
stamp(file.stem)

# Plot the flux time evolution for a given range (limited to the reco axis range)
heading("Flux versus time")
fig,ax = plt.subplots(figsize=(20,8))
flux_versus_time(dsets,
                 emin=100*u.GeV, emax=5*u.TeV,
                 tmin=300*u.s, tmax=2*u.d,
                 stacked=False, fit=True, debug=True)
stamp(file.stem)
