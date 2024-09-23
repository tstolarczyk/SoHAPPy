# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:27:45 2024

@author: Stolar
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

from pathlib import Path
from niceplot import MyLabel
from niceprint import heading

from population import Pop
from pop_io import get_data

from configuration import Configuration


# -----------------------------------------------------------------------------
def det_proba(tmax_slew=30, dtmin=30, ntrials=1000, plot=True, debug=False):
    """
    Get detection probability

    Parameters
    ----------
    tmax_slew : float, optional
        maximal slewing time. The default is 30.
    dtmin : float, optional
        Minimal required observation time for physics. The default is 30.
    ntrials : integer, optional
        Number of trials. The default is 100.
    debug : boolean, optional
        If True, talkative. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # Loop over t90s, generate slewing delays and Swift delays
    success = []
    for t90 in t90s:

        # Generate Swift delays
        t_alert = np.random.choice(tproba, ntrials, p=proba)

        # Generate slewing delay
        t_slew = np.random.random(ntrials)*tmax_slew

        # Total delay
        t_delay = t_alert + t_slew

        # Find indices of t90 > t_delay and check there is enough time
        # for physics
        idxs = np.where(t_delay + dtmin <= t90)

        # Count the number of time delay left dtmin for observation before the
        # prompt was over
        passed = 0 if len(idxs) == 0 else len(idxs[0])
        success.append(passed/ntrials)

        if debug:
            print(f" t90 = {t90:7.2f}"
                  f" <tdelay> = {np.mean(t_delay):7.2f}"
                  f" passed = {passed:6d}"
                  f" fraction = {100*passed/ntrials:3.1f} %")

    if plot:
        label = f"Slewing max= {tmax_slew:4.1f} - Det. min={dtmin:4.1f}"
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(t90s, success, label=label)
        ax.set_xlabel("t90")
        ax.set_ylabel("Detectable")
        ax.legend()

    print(f"Prompt mean detection probability: {np.mean(success):5.2f} "
          f"[dtmin={dtmin:4.1f}, tmax_slew={tmax_slew:4.1f}]")

    return np.mean(success)


# #############################################################################

# Bigger texts and labels
sns.set_context("notebook")  # poster, talk, notebook, paper

# ##--------------
# ## SWIFT  latency data
# ##--------------
heading("Swift latency data")

cfg = Configuration()
print("File : ", cfg.swiftfile)

dt_swift = []
with open(Path("../../", cfg.swiftfile)) as input_file:
    for line in input_file:
        # print(line)
        dt_swift.append(float(line.strip().split(" ")[-1]))
dt_swift = np.array(dt_swift)

# Compute occurence probability of Swift delay from an histogram
nbin = 50
proba, t_edges = np.histogram(dt_swift, bins=nbin, density=True, weights=None)
proba = proba/sum(proba)  # Renormalise to 1
tproba = (t_edges[1:] + t_edges[:-1])/2


# ##--------------
# ## GRB t90s
# ##--------------

# Read any kind ot population with the original GRB T90 information
nyears, files, tag = get_data(debug=False)
pop = Pop(files, tag=tag, nyrs=nyears)
t90s = pop.ref.t90
Ndata = len(t90s)


heading("Statistics")
print(" File : ", files[0])
print(" - Ndata = ", Ndata)
print(" - Minimal t90 : ", min(t90s))
print(" - Maximal t90 : ", max(t90s))
print()

# Critical values
txtlist = ["Short GRB limit",
           "LST max. slew",
           "Swift mean latency",
           "MST max. slew",
           "Max. delay (LST)",
           "Max. delay (MST)"]
delays = [2, 30, 70, 90, 107, 137]
ntxtmax = max([len(txt) for txt in txtlist])

# Statistics
print("Statistics for t90 above")
for txt, val in zip(txtlist, delays):
    n = len(t90s[t90s >= val])
    print(f" - {txt:{str(ntxtmax + 1)}s} [{val:3.0f} s]:"
          f" {100*n/Ndata:8.2f} %")

# ##--------------
# ## PLOTS
# ##--------------

# #-------------------------------------
color = cm.rainbow(np.linspace(0, 1, len(delays)))


def plot_ref_time(axx):
    i = 0
    for txt, val in zip(txtlist, delays):
        label = f"{txt:s} [{val:3.0f} s]"
        ax.axvline(np.log10(val), label=label, ls=":", color=color[i])
        i += 1
# #-------------------------------------


# Display t90
fig, ax = plt.subplots(figsize=(12, 5))
n, bins, _ = ax.hist(np.log10(pop.ref.t90), bins=25, alpha=0.5,
                     label=MyLabel(np.log10(pop.ref.t90), "$t_{90}$"))

plot_ref_time(ax)
ax.set_xlabel("$t_{90}$")
ax.legend()

# Display t90 and swift delays
fig, ax = plt.subplots(figsize=(12, 5))
n, bins, _ = ax.hist(np.log10(pop.ref.t90), bins=100, alpha=0.5,
                     label=MyLabel(np.log10(pop.ref.t90), "$t_{90}$"))
ax.hist(np.log10(dt_swift), bins=bins,
        label=MyLabel(np.log10(dt_swift), " Swift delays"))
plot_ref_time(ax)

ax.legend()

# ##--------------
# Compute probability to have enough observation time
# ##--------------
heading("Detection probability")

det_times = [10, 30, 60, 90, 120, 180, 360, 600, 3000]
det_times = [10, 30, 60, 90, 120]

# --------------------------
# LST max slewing time
# --------------------------


# for the plot versus t90
det_proba(tmax_slew=30)

# Scan detection time
fLST = []
fLST20 = []
for dtmin in det_times:
    frct = det_proba(tmax_slew=30, dtmin=dtmin, plot=False)
    frct20 = det_proba(tmax_slew=20, dtmin=dtmin, plot=False)
    fLST.append(frct)
    fLST20.append(frct20)

# --------------------------
# SST max slewing time
# --------------------------
# for the plot versus t90
det_proba(tmax_slew=60)

# Scan detection time
fSST = []
for dtmin in det_times:
    frct = det_proba(tmax_slew=60, dtmin=dtmin, plot=False)
    fSST.append(frct)

# --------------------------
# MST max slewing time
# --------------------------
# for the plot versus t90
det_proba(tmax_slew=90)

# Scan detection time
fMST = []
for dtmin in det_times:
    frct = det_proba(tmax_slew=90, dtmin=dtmin, plot=False)
    fMST.append(frct)

# --------------------------
# Summary plot
# --------------------------
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(det_times, fLST20, label=" Max. slew time:\n20s (LST)",
        marker="o", color="tab:blue", ls=":")
ax.plot(det_times, fLST, label=" 30s (LST)",
        marker="o", color="tab:blue")
ax.plot(det_times, fMST, label=" 90s (MST)",
        marker="o", color="tab:green")
ax.plot(det_times, fSST, label="  60s (SST)",
        marker="o", color="tab:orange")
ax.set_xlabel("Time for physics (s)")
ax.set_ylabel("Fraction of characterizable prompts")
ax.grid(which="both")
ax.set_title("Prompt signal accessible to detection")
ax.legend()
