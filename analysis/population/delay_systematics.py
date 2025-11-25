# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:38:43 2025.

@author: Stolar
"""

import sys

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from niceprint import heading

from configuration import Configuration

import rate
from pop_io import get_data_from_folder
from population import Pop

sys.path.append("../../")


###############################################################################
if __name__ == "__main__":

    # Bigger texts and labels
    import seaborn as sns
    sns.set_context("talk")  # poster, talk, notebook, paper

    # Build a default Configuration to use the sourc_ids function
    sys.argv = ["", "-c", "../../data/config_ref.yaml", "-f", "1"]
    cf = Configuration.command_line()

    base = Path(r"D:\CTAO\SoHAPPy\output\100k_long_ISM")
    delays = [0, 15, 30, 60, None, 200, 500, 1000, 2000]

    # Choose whihc stat to display
    d5s90 = False  # 5 sigma at 90%CL, otherwise mean rate at 5 sigma

    ndet5s_tot = np.zeros(len(delays))
    ndet5s_n = np.zeros(len(delays))
    ndet5s_s = np.zeros(len(delays))

    for item, dt in enumerate(delays):

        # Build data folder name
        if dt is not None:
            dttext = "_dt" + str(dt).zfill(3)
        else:
            dttext = ""
            id_def = item

        folder = Path(base, "prod5_100000_std" + dttext + "/strictmoonveto")
        print(folder, end=" ")
        if folder.exists():
            print("found")
        else:
            print(" NOT FOUND - Skipped")
            continue

        # Create population from current files in folder
        heading(" " + str(item) + ": dt=" + str(dt))
        tag = str(dt)
        nyears = 210.6
        filelist = get_data_from_folder(folder)

        pop = Pop(filelist, tag=tag,  nyrs=nyears, fpeakmin=1, debug=False)
        pop.print()
        # pop.sanity_check()
        rate.compute(pop)

        if d5s90:
            ndet5s_tot[item] = len(pop.g_tot[pop.g_tot.d5s >= pop.eff_lvl])
            ndet5s_n[item] = len(pop.g_n[pop.g_n.d5s >= pop.eff_lvl])
            ndet5s_s[item] = len(pop.g_s[pop.g_s.d5s >= pop.eff_lvl])
        else:
            ndet5s_tot[item] = np.sum(pop.g_tot.d5s)/pop.niter
            ndet5s_n[item] = np.sum(pop.g_n.d5s)/pop.niter
            ndet5s_s[item] = np.sum(pop.g_s.d5s)/pop.niter

    # ## Plot variations
    dt_north = (cf.dtslew_north + cf.dtswift).value
    dt_south = (cf.dtslew_south + cf.dtswift).value
    print(" Standard delays : N=", dt_north, " South=", dt_south)

    fig, ax = plt.subplots(figsize=(10, 10))

    # ## Total
    # Raplace zero into 1 and remove None form the data"
    dts = [dt if dt != 0 else 1 for dt in delays]  # Remove zero
    none_idx = dts.index(None)
    dts.pop(none_idx)
    ndet = np.delete(ndet5s_tot, none_idx)
    ax.errorbar(dts,
                ndet/pop.nyears,
                yerr=np.sqrt(ndet)/pop.nyears,
                marker="o", label="Combined", color="tab:green")

    # South
    dts = [dt if dt is not None else dt_south for dt in delays]
    dts = [dt if dt != 0 else 1 for dt in dts]
    ax.errorbar(dts, ndet5s_s/pop.nyears,
                yerr=np.sqrt(ndet5s_n)/pop.nyears,
                marker="o", color=pop.col_s, label="South")
    ax.axvline(dt_south, ls="--", color=pop.col_s,
               label="Standard delay South")

    # North
    dts = [dt if dt is not None else dt_north for dt in delays]
    dts = [dt if dt != 0 else 1 for dt in dts]
    ax.errorbar(dts, ndet5s_n/pop.nyears,
                yerr=np.sqrt(ndet5s_s)/pop.nyears,
                marker="o", color=pop.col_n, label="North")
    ax.axvline(dt_north, ls="--", color=pop.col_n,
               label="Standard delay North")
    ax.minorticks_on()
    ax.tick_params(which='both', width=2)
    ax.grid(which="major", axis="both", color="grey", ls="-")
    ax.grid(which="minor", axis="both", color="grey", alpha=0.5,  ls=":")
    if d5s90:
        ax.set_ylabel(r"$ Detection \ rate \ (5 \sigma \ 90 \%) \ yr^{-1}$")
    else:
        ax.set_ylabel(r"$ Detection \ rate \ (5 \sigma) \ yr^{-1}$")

    ax.set_xlabel(r"$Overall \ time \ delay \ (s)$")
    ax.set_xscale("log")
    ax.legend(bbox_to_anchor=(1.0, 1.00))

    print("...completed")
