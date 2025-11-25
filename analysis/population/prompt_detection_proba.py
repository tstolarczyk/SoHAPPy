# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:27:45 2024.

@author: Stolar
"""
import sys
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


# Critical values in seconds
critical_times = {"Short GRB limit": 2,
                  "LST max. slew": 30,
                  "MST max. slew": 90,
                  "SST max. slew": 60,
                  "Max. ref. delay (LST)": 107,
                  "Max. ref. delay (MST)": 167
                  }


# -----------------------------------------------------------------------------
def det_proba(tproba, proba, tmax_slew=30, dtmin=30, ntrials=1000,
              plot=True, bins=None, debug=False):
    """
    Get detection probability.

    For each GRB, generate ntrials Swift and slewing delays.
    The Swift alert delays is obatined from the probability distribution.
    The slewing delay is assumed to be flat until tmax_slew.
    Compute the fraction of trails where the remaining time to t90, dtmin, is
    considered enough for doing an analysls.

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
        if bins is None:
            bins = np.logspace(-1, 5.2, 50)

        label = f"Slewing max= {tmax_slew:4.1f} - Det. min={dtmin:4.1f}"
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(t90s, success, marker=".", label=label)
        ax.set_xlabel("t90")
        ax.set_ylabel("Detectable")
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.grid(axis="both")
        ax.legend()

    print(f"Prompt mean detection probability: {np.mean(success):5.2f} "
          f"[dtmin={dtmin:4.1f}, tmax_slew={tmax_slew:4.1f}]")

    return np.mean(success)  # Return the main success


# ##---------------------------------------------------------------------------
def Swift_delay(filename, debug=True):
    """Get Swift delays and probability distribution."""
    dt_swift = []
    with open(Path("../../", filename)) as input_file:
        for line in input_file:
            # print(line)
            dt_swift.append(float(line.strip().split(" ")[-1]))
    dt_swift = np.array(dt_swift)

    # Compute occurence probability of Swift delay from an histogram
    nbin = 50

    proba, t_edges = np.histogram(dt_swift, bins=nbin,
                                  density=True, weights=None)
    proba = proba/sum(proba)  # Renormalise to 1
    tproba = (t_edges[1:] + t_edges[:-1])/2

    if debug:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot(tproba, proba, color="orange")
        ax.set_ylabel("Probability")

    return dt_swift, tproba, proba


# ##---------------------------------------------------------------------------
def Swift_delay_plot(times, bins=None, ax=None, **kwargs):
    """Plot Swift delays."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    if bins is None:
        bins = 25

    ax.hist(times, bins=bins, label=MyLabel(times, "Swift delays"), **kwargs)
    ax.set_xlabel(" Swift/BAT alert delay (s)")
    ax.set_ylabel("Counts")
    ax.grid(axis="both")
    ax.legend()


# ##---------------------------------------------------------------------------
def GRB_t90s(parfilename, fpeakmin=1, mask=None):
    """Get GRB t90's from any output file."""
    nyears, files, tag = get_data(parpath=parfilename, debug=False)

    pop = Pop(files, tag=tag, nyrs=nyears, fpeakmin=fpeakmin)

    if mask is None:
        t90s = pop.ref.t90
    elif mask == "d5s":
        t90s = pop.g_tot[pop.g_tot.d5s >= pop.eff_lvl].t90
    else:
        sys.exit("This mask isnot implemented")

    print(" File : ", files[0])
    print(" - Ndata = ", len(t90s))
    print(" - Minimal t90 : ", min(t90s))
    print(" - Maximal t90 : ", max(t90s))
    print()

    return t90s


# ##---------------------------------------------------------------------------
def GRB_t90s_plot(times, bins=None, ax=None, **kwargs):
    """Plot t90s distributions."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    if bins is None:
        bins = np.logspace(-1, 5.2, 50)

    ax.hist(times, bins, label=MyLabel(t90s, r"$t_{90}$"), **kwargs)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("Counts")
    ax.set_xlabel("Analysed GRB t90s")
    ax.grid(axis="x", ls=":")
    ax.legend()

    # axx = ax.twinx()

    # hist, edges = np.histogram(times, bins=bins)
    # xcenter = 0.5*(edges[1:] + edges[:-1])
    # vals = np.cumsum(hist)/np.sum(hist)
    # axx.plot(xcenter, vals, marker="", color="tab:orange")
    # axx.set_ylabel("Cumulated fraction")
    # axx.grid(axis="both", ls=":")


# ##---------------------------------------------------------------------------
def specific_times_plot(ax):
    """Plot some critical time values."""
    if ax is None:
        print("Specific times - cannot plot in the abscence of an axis")
        return

    for key, val in critical_times.items():
        print(key, val)
        # label = f"{key:s} [{val:3.0f} s]"
        ax.axvline(val, ls=":", color="grey")
        ax.text(val*1.05, ax.get_ylim()[1]*0.92, s=key,
                va="top", rotation=90, size=7, color="black")


# #############################################################################

# Bigger texts and labels
sns.set_context("notebook")  # poster, talk, notebook, paper

# ##--------------
# ## SWIFT  latency data
# ##--------------
heading("Swift latency data")

cfg = Configuration()
print("File : ", cfg.swiftfile)
dtswift, tproba, proba = Swift_delay(cfg.swiftfile)

# # ##--------------
# # ## GRB t90s
# # ##--------------

# Read any kind ot population with the original GRB T90 information
heading("GRB t90's")

# nyears, files, tag = get_data(parpath=None, debug=True)
parpath = "parameter_100k_ISM.yaml"
# parpath = "parameter_100k_ISM-max.yaml"

t90s = GRB_t90s(parpath, fpeakmin=1, mask="d5s")
GRB_t90s_plot(t90s)

# # ##--------------
# # ## GRB t90s and Swift delays plot
# # ##--------------

fig, ax = plt.subplots(figsize=(8, 5))
bins = np.logspace(-1, 4.2, 100)

GRB_t90s_plot(t90s, ax=ax, bins=bins, alpha=0.5)
Swift_delay_plot(dtswift, bins=bins, ax=ax, color="tab:orange", alpha=0.5)
specific_times_plot(ax)
ax.set_xlabel("Time (s)")  # Supersede existing label

# Statistics
print()
print("Statistics for t90 above")
for key, val in critical_times.items():
    n = len(t90s[t90s >= val])
    print(f" - {key:22s} [{val:3.0f} s]:"
          f" {100*n/len(t90s):8.2f} %")


# ##--------------
# Compute probability to have enough observation time
# ##--------------
heading("Detection probability")

det_times = [10, 30, 60, 90, 120, 180, 360, 600, 3000]
det_times = [10, 30, 60, 90, 120]

# --------------------------
# LST max slewing time
# --------------------------
# Display the fraction of detectable GRB from their t90 for a chosen minimal
# detection duration.
det_proba(tproba, proba, dtmin=30, tmax_slew=30)

# Do the same for various detection times - Store the results
fLST = []
fLST20 = []
for dtmin in det_times:
    frct = det_proba(tproba, proba, tmax_slew=30, dtmin=dtmin, plot=False)
    frct20 = det_proba(tproba, proba, tmax_slew=20, dtmin=dtmin, plot=False)
    fLST.append(frct)
    fLST20.append(frct20)

# # --------------------------
# # SST max slewing time
# # --------------------------
# # for the plot versus t90
# det_proba(tmax_slew=60)

# Scan detection time
fSST = []
for dtmin in det_times:
    frct = det_proba(tproba, proba, tmax_slew=60, dtmin=dtmin, plot=False)
    fSST.append(frct)

# # --------------------------
# # MST max slewing time
# # --------------------------
# # for the plot versus t90
# det_proba(tmax_slew=90)

# # Scan detection time
fMST = []
for dtmin in det_times:
    frct = det_proba(tproba, proba, tmax_slew=90, dtmin=dtmin, plot=False)
    fMST.append(frct)

# # --------------------------
# # Summary plot
# # --------------------------
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
