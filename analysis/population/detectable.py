# -*- coding: utf-8 -*-

"""
Investigate the GRB subpopulation which is never detected.

This can be used in particular on a population which was simulated and analysed
with a configuration of maximal dectection (infinite nights, no delay,
high altitude). Once the undetected population is know, it can be avoided
to simulated and analysed in standard conditions.

Created on Wed Jul 10 11:45:41 2024

@author: Stolar
"""
import os
import sys
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

from population import Pop
from pop_io import get_data
from pop_plot import plot_generated

from niceplot import single_legend, draw_contours, lower_limit
from niceprint import heading

# Bigger texts and labels
sns.set_context("notebook")  # poster, talk, notebook, paper

codefolder = "../../"
sys.path.append(codefolder)


# -----------------------------------------------------------------------------
def dump_detectable(namelist, subfolder=None, outname=None, chunk=2000, dgt=5):
    """
    Obtain and store the full paths of a a series of GRBs.

    GRBs identified from their name and the location of the folder where they
    shall be found. The found path are stored in json files in specific folders
    created in chunk of a fixed size in order to ease further use on batch
    nodes.
    Note that the HAPPY_IN variable shall be set in oreder to get a list of
    files with a relative path that can be used in another environement.

    Parameters
    ----------
    namelist : numpy array of strings
        A list of source names (can be obtained from a popualtion file).
    folder : Path, optional
        THe Path where to find the original files. The default is None.
    outname : Path or string, optional
        The output file name. The default is None.
        Any given extension is replaced by the json extension.
    dgt: integer
        Number of digits on which integers are shown in the filename

    Returns
    -------
    None.

    """
    if "HAPPY_IN" not in os.environ:
        sys.exit("HAPPY_IN environement variable shall be defined")
    else:
        base = Path(os.environ["HAPPY_IN"])  # This base can be different

    if outname is None:
        sys.exit(" dump_detectable: provide a valid output file name")
    else:
        outname = Path(Path(outname).stem + ".json")

    if subfolder is None:
        sys.exit(" dump_detectable: provide a valid subfolder name")

    # ## -----------------
    # ## Original folder
    # ## -----------------

    # Get full file names from the original folder
    # look in the given folder and below, try to find fits files, and if
    # necessary compressed fits files.
    folder = Path(base, subfolder)
    filelist = list(folder.glob(r"**/*.fits"))
    if len(filelist) == 0:
        filelist = list(folder.glob(r"**/*.fits.gz"))
        if len(filelist) == 0:
            sys.exit(f"Dump_detectable: no fits nor fits.gz files found "
                     f"below {folder:}")
        else:
            ext = ".fits.gz"
    else:
        ext = ".fits"

    # ## -----------------
    # ## Find matches
    # ## -----------------
    # Loop over the two lists and find matching
    # Add the extension to the names to avoid mutiple matches

    istart = 1
    nfound = 0
    fullnames = []

    for ifile, file in enumerate(filelist, start=1):

        for name in namelist:

            # If found, store the file name in its subfolder
            # Remove the base (HAPPY_IN) but keep the subfolder
            if str(file.name).find(name+ext) != -1:
                idx = file.as_posix().find(base.as_posix())\
                     + len(base.as_posix())
                fullnames.append(file.as_posix()[idx+1:])
                nfound += 1
                # print('++', nfound, "  ", str(file.name), " * ", name+ext)
                break

        if np.mod(ifile, chunk) == 0:
            filename = outname.stem + "_"\
                    + str(istart).zfill(dgt) + "_" + str(ifile).zfill(dgt)
            print(f" Create file {filename:}"
                  f" with {len(fullnames):} entries")

            # Dump to json file
            with open(filename+".json", 'w') as fj:
                json.dump(fullnames, fj)
            istart = ifile + 1
            fullnames = []

    if nfound != len(namelist):
        print(f" Sources in name list = {len(namelist):}, "
              f"but {nfound:} found in {folder:}")
        sys.exit(" dump_detectatble: matching problem")
    else:
        print(f" Sources in name list = {len(namelist):}, "
              f"and {nfound:} found in {folder:} - all is OK")


# -----------------------------------------------------------------------------
def statistics(pops, tags, cut="d3s", cut_min=0):
    """
    Dump detecttion statitics of a series of population.

    Parameters
    ----------
    pops : List of pandas table
        5sub)population to be analysed.
    tags : List of String
        tags describing each (sub)population.
    cut : String, optional
        Variable name on which a cut is applied. Default is "d3s".
    cut_min: float, optional
        Cut value (strict), default is 0.

    Returns
    -------
    None.

    """
    def sep():
        """Display a simple separator."""
        print(f"{10*'-':>10s}+ {8*'-':>8s} {8*'-':>8s} {8*'-':>8s}"
              f" {8*'-':>8s}")

    sep()
    print(f"{'Sub pop.':<10s}: {'Missed':>8s} {'Fraction':>8s} "
          f"{'Det.':>8s} {'Total':>8s}")
    sep()
    for gsub, tag in zip(pops, tags):
        n_tot = len(gsub)
        n_miss = len(gsub[gsub[cut] <= cut_min])
        n_det = len(gsub[gsub[cut] > cut_min])
        if n_tot != n_miss + n_det:
            sys.exit(" Problems in counting detected and not detected")

        print(f"{tag:10s}: {n_miss:8d} {100*n_miss/n_tot:7.2f}% "
              f"{n_det:8d} {n_tot:8d}")

        sep()


# -----------------------------------------------------------------------------
def coverage_and_fit(pops, tags,
                     varx, vary, logx=False, logy=True, nbins=50,
                     cut="d3s", cut_min=0,
                     alt="z", alt_max=5):
    """
    Plot the (x,y) coverage of a population and find lower limit.

    Parameters
    ----------
    pops : List of pandas table
        5sub)population to be analysed.
    tags : List of String
        tags describing each (sub)population.
    varx : String
        Abscissa variable name.
    vary : String
        Ordinate variable name.
    logx : Boolean, optional
        If True, x axis in log scale. The default is False.
    logy : Boolean, optional
        If True, y axis in log scale. The default is True.
    nbins : Integer, optional
        Number of bins of the 2D histrogram axis. The default is 50.
    cut : String, optional
        Variable name on which a cut is applied. Default is "d3s".
    cut_min: float, optional
        Cut value (strict), default is 0.
    alt : String. optional.
        The name of the variable for which the values will be limited.
        The default is "z".
    alt_max: float. Optional.
        The value not to be eexceeeded by the variable above.

    Returns
    -------
    None.

    """
    print("   x    : ", varx)
    print("   y    : ", vary)

    # Lower limit finding
    fitdeg = 4  # Fit order for lower limit
    thrs = 0.1  # Threshold for the limit
    norm = "sum"  # alternative is "max"

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6), sharey=True)

    first = True

    for gsub, tag, ax in zip(pops, tags, axes):

        # All values
        x_all = gsub[gsub[alt] <= alt_max][varx]
        y_all = gsub[gsub[alt] <= alt_max][vary]

        # Superimpose contours (they are transparent, unlike hsitrograms)
        levels = [50, 100, 300, 1000]
        HH, xe, ye = np.histogram2d(gsub.z, np.log10(gsub.Eiso), bins=25)
        grid = HH.transpose()
        midpoints = (xe[1:] + xe[:-1])/2, (ye[1:] + ye[:-1])/2
        cntr = ax.contour(*midpoints, grid,
                          colors="grey", linestyles="dashed",
                          levels=levels)
        clbl = cntr.clabel()  # Labels
        if logy is False:
            ax.set_yscale("log")
        ax.grid("both")
        ax.set_title(tag)

        # ax.scatter(x_all, y_all, alpha=0.1, color="grey", marker=".")
        ax.hist2d(x_all, y_all, alpha=0.5, cmap="binary", bins=nbins)
        ax.set_xlabel(("$log_{10}$" if logx else "") + varx)

        mask = (gsub[cut] > cut_min) & (gsub[alt] <= alt_max)
        x_det = gsub[mask][varx]
        y_det = gsub[mask][vary]

        if logx:
            x_det = np.log10(x_det)
            x_all = np.log10(x_all)
        if logy:
            y_det = np.log10(y_det)
            y_all = np.log10(y_all)

        if first:
            pfit0 = lower_limit(x_det, y_det,
                                norm=norm, fitdeg=fitdeg, threshold=thrs,
                                ax=ax, nbins=nbins,
                                color="red", cbar=False, cmap="PuBu")

            xlim = ax.get_xlim()
            xsample = np.linspace(xlim[0], xlim[1], 20)
            ax.set_ylabel(("$log_{10}$" if logy else "") + vary)
            first = False

        else:
            _ = lower_limit(x_det, y_det,
                            norm=norm, fitdeg=fitdeg, threshold=thrs,
                            ax=ax, nbins=nbins,
                            color="green", cmap="PuBu",
                            cbar=False, plot_lim=False)

            ax.plot(xsample, np.polyval(pfit0, xsample), ls="--", alpha=0.5,
                    color="red", label="Reference")

        # ax.set_xlim(xmin=np.min(x_all))
        ax.set_ylim(ymin=50, ymax=55)

    single_legend(fig)
    plt.tight_layout()

    return pfit0


# #############################################################################
if __name__ == "__main__":

    # Bigger texts and labels
    sns.set_context("talk")  # poster, talk, notebook, paper

    # ---------------------------------------------
    # Max detection default dataset (1000 GRBs)
    # ---------------------------------------------
    # parpath = "../../data/samples/max_detection_parameter.yaml"
    # os.environ["HAPPY_IN"] = r"D:\\CTAO\SoHAPPy\input"
    # data_dir = Path(r"lightcurves\long_1_1000")

    # ---------------------------------------------
    # Standard stricmoonveto demo dataset (1000 GRBs)
    # ---------------------------------------------
    # parpath = None
    # os.environ["HAPPY_IN"] = r"D:\\CTAO\SoHAPPy\input"
    # data_dir = Path(r"lightcurves\long_1_1000")

    # ---------------------------------------------
    # Max detection new sample(>10000 GRBs)
    # ---------------------------------------------
    parpath = "parameter_100k_ISM-max.yaml"
    data_dir = Path(r"lightcurves\100k_long_ISM")

    nyears, popfiles, tag = get_data(parpath=parpath, debug=False)

    # ## -------------------
    # ## Definitions
    # ## -------------------

    heading(" Fraction of never detected GRB")
    cut = "d3s"
    cut_min = 0.
    alt = "z"
    alt_max = 5.0

    print(" Analysing population with :")
    print("   cut           : ", cut)
    print("   min. (strict) : ", cut_min)
    print("   alt. cut      : ", alt)
    print("   alt. cut max  : ", alt_max)
    print()

    # Read population
    heading(" Population reading")

    pop = Pop(popfiles, tag=tag, nyrs=nyears, debug=False)

    # Stat on total (combined) population
    print()
    print(f"--- File:  {popfiles[0].parent:}...")
    print("  Detection level at 3 sigma:")
    print(f"   Range: [{np.min(pop.g_tot.d3s):4.1f}, "
          f"{np.max(pop.g_tot.d3s):4.1f}]")
    print(f"   Mean : {np.mean(pop.g_tot.d3s):4.1f} "
          f"+- {np.std(pop.g_tot.d3s):4.1f}")

    # Define population list for multiple plots
    poplist = [pop.g_n, pop.g_s, pop.g_tot]
    taglist = ["North", "South", "Combined"]

    # ## -------------------
    # ## Dump combined with full path
    # ## -------------------
    # heading("Combined subpopulation file")

    # mask = pop.g_tot[cut] > cut_min
    # namelist = [x["name"] for _, x in pop.g_tot[mask].iterrows()]

    # dump_detectable(namelist, subfolder=data_dir,
    #                     outname="combined_detected")

    # # ## -------------------
    # # ## Display generated
    # # ## -------------------
    plot_generated(pop)

    # ## -------------------
    # ## (Non)detectable subpopulation statistics
    # ## -------------------
    heading("Detected Subpopulations")
    statistics(poplist, taglist, cut=cut, cut_min=cut_min)

    # # z coverage - that's done already in the pop_plot.distri function?
    # varx = "z"
    # print("   x    : ", varx)
    # logx = False
    # nbins = 50

    # first = True
    # fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 4), sharey=True)

    # for gsub, tag, ax in zip(poplist, taglist, axes):

    #     # All values
    #     x_all = gsub[gsub[varx] <= alt_max][varx]
    #     mask = (gsub[cut] > cut_min) & (gsub[varx] <= alt_max)
    #     x_det = gsub[mask][varx]

    #     _, bins, _ = ax.hist(x_all, bins=nbins,
    #                          facecolor="none", edgecolor="black")
    #     ax.hist(x_det, bins=bins, color="tab:blue")
    #     ax.set_xlabel("Redshift z)")
    #     ax.set_yscale("log")
    #     if first:
    #         ax.set_ylabel("Counts")
    #         first = False
    #     ax.grid(which="both", ls=":")
    # plt.tight_layout()

    # ## -------------------
    # ## Lower limit
    # ## -------------------
    # Fit a function representing the lower limit of the populated plane
    # defined by the follwoing variables. For a series of data, superimpose
    # the first fit to the subsequent one.
    # ## Z, and Eiso coverages
    

# WRONG ????
#     heading("Lower  limits")

#     first = True
#     fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 4), sharey=True)
#     for gsub, tag, ax in zip(poplist, taglist, axes):

#         # All values
#         y_all = gsub[gsub[alt] <= alt_max][vary]

#         # Detected values
#         mask = (gsub[cut] > cut_min) & (gsub[alt] <= alt_max)
#         y_det = gsub[mask][vary]

#         _, bins, _ = ax.hist(np.log10(y_all), bins=nbins,
#                              facecolor="none", edgecolor="black")
#         ax.hist(np.log10(y_det), bins=bins, color="tab:blue")
#         ax.set_xlabel(r"$log_{10} (Eiso) (erg)$")
#         ax.set_yscale("log")
#         if first:
#             first = False
#             ax.set_ylabel("Counts")
#         ax.grid(which="major", ls=":")
#         ax.set_title(tag)
#     plt.tight_layout()

    heading("Lower detection limits")
    fit_lowerlimit = coverage_and_fit(poplist, taglist, "z", "Eiso")

    # # ## -------------------
    # # ## Subpopulation contours
    # # ## -------------------
    # logx = False
    # logy = True

    # fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5), sharey=True)

    # first = True
    # for gsub, tag, ax in zip(poplist, taglist, axes):

    #     mask = (gsub[cut] <= cut_min) & (gsub[alt] <= alt_max)
    #     x_miss = gsub[mask][varx]
    #     y_miss = gsub[mask][vary]
    #     mask = (gsub[cut] > cut_min) & (gsub[alt] <= alt_max)
    #     x_det = gsub[mask][varx]
    #     y_det = gsub[mask][vary]

    #     xval = x_miss
    #     yval = y_miss
    #     if logx:
    #         xval = np.log10(xval)
    #     if logy:
    #         yval = np.log10(yval)

    #     draw_contours(xval, yval, ax=ax, nbins=10,
    #                   levels=15)
    #     ax.set_title(tag)
    #     # zorder=100, norm="linear",
    #     #                   colors="black", alpha=0.2, ls=":")

    #     ax.set_xlabel(("$log_{10}$" if logx else "") + varx)
    #     if first:
    #         ax.set_ylabel(("$log_{10}$" if logy else "") + vary)
    #         first = False

    #     if logy is False:
    #         ax.set_yscale("log")
    #     ax.grid("both")

    #     # ax.set_xlim(xmin=np.min(np.log10(pop.g_tot[varx])))
    #     # ax.set_ylim(ymin=np.min(gsub[vary]))
    plt.tight_layout(h_pad=0, w_pad=0)
