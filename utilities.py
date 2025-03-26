# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019.

A bucnh of utility functions used in `SoHAPPy`.

@author: Stolar
"""

import tarfile
from pathlib import Path

import sys
import math
import os
import shutil
import datetime
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

__all__ = ["log_and_flat_intervals", "subset_ids",
           "get_filename",
           "file_from_tar", "backup_file", "Df", "Dp"]


# ##---------------------------------------------------------------------------
def log_and_flat_intervals(xmin, xmax, nbin=10, distmax=0., debug=True):
    """
    Generate logarithmically-spaced interavlls until a limit, then flat.

    Parameters
    ----------
    xmin : float
        Minimal value.
    xmax : float
        maximal value.
    nbin : integer, optional
        Initial number of bins. The default is 10.
    distmax : float, optional
        Limit in time interval length. The default is 0..
    debug : boolean, optional
        If True illustrate the process. The default is True.

    Returns
    -------
    Numpy array
        Final list.

    """
    # First generated nbin log-spaced intervals
    xsample = np.logspace(xmin, xmax, nbin)

    if distmax == 0:
        return xsample

    if xmin <= 0:
        print(" Starting value shall be >0")
        return False

    # Original log-spaced intervals
    xsamp = np.logspace(np.log10(xmin), np.log10(xmax), nbin)
    dxsamp = np.diff(xsamp)

    # Find intervals exceeding the maximal distance
    idx = np.where(dxsamp > distmax)[0][0]
    # From that value,generate flat spacing
    xflat = np.arange(xsamp[idx], xmax, distmax)
    # Concatenate the two arrays
    xfinal = np.concatenate((xsamp[:idx], xflat))

    if debug:
        print(" Setting intervals in [", xmin, xmax, "]")
        print("   - number of bins        = ", nbin)
        print("   - max. allowed distance = ", distmax)
        print(" Distance is exceeded at index ", idx, " val= ", xsamp[idx+1])

        fig, ax = plt.subplots(figsize=(8, 3))
        xorig = xsamp[xsamp < max(xfinal)]  # Limit original sample
        xorig = xsamp
        print(" sample = ", xsamp)
        print(" limited xorig = ", xorig)
        print(" final = ", xfinal)
        ax.plot(xfinal, np.ones(len(xfinal)), marker=".",
                label="Final")
        ax.plot(xsamp[:idx], 2*np.ones(len(xsamp[:idx])), marker=".",
                label="log seq.")
        ax.plot(xflat, 1.5*np.ones(len(xflat)), marker=".",
                label="lin seq.")
        ax.set_yticklabels([])
        ax.axvline(xsamp[idx], ls=":")
        ax.legend()
        axx = ax.twinx()
        axx.plot(xfinal[1:], np.diff(xfinal), ls=":", color="red",
                 label="Final intervals", marker=".")
        axx.plot(xorig[1:], np.diff(xorig), ls=":", color="grey",
                 label="Initial intervals", marker=".")
        axx.axhline(distmax, ls=":", color="tab:orange",
                    label="Max. allowed")
        axx.set_yscale("log")
        axx.set_ylabel("Interval spacing")
        axx.legend()

        fig.legend(bbox_to_anchor=[1.2, 0.9])
        for axis in fig.axes:
            lgd = axis.get_legend()
            if lgd is not None:
                lgd.remove()

    return xfinal


# ##---------------------------------------------------------------------------
def subset_ids(istart, nmax, nsets, debug=False):
    """
    This is a simple code defining `nsets` interval for a list of integer
    starting at 1 up to `nmax`. If `nmax` is not divisible by `nsets`, the rest
    gives and extra set.

    Parameters
    ----------
    nmax : integer
        Maximal nubmer.
    nsets : integer
        Number of sets.
    debug : Boolean, optional
        If True, verbosy. The default is False.

    Returns
    -------
    dsets : list of List
        The list of concsecutiveintervals.

    """
    if nsets >= nmax:
        return [1, nmax]

    delta = int(nmax/nsets)  # regular set size
    rest = nmax % nsets  # Last set size

    # Define low and high values of the intervals
    ids_low = list(range(1, nmax, delta))
    ids_high = list(range(delta, min(nmax+delta, nmax), delta))
    if len(ids_high) < len(ids_low):
        ids_high += [nmax]
    elif ids_high[-1] < nmax:
        ids_high[-1] = nmax

    if debug:
        print("steps = ", delta, " uncovered = ", rest)
        print(ids_low, " -> ", len(ids_low), " items")
        print(ids_high, " -> ", len(ids_high), " items")

    # Create interval list
    dsets = [[idl+istart-1, idh+istart-1]
             for (idl, idh) in zip(ids_low, ids_high)]

    # Check counts
    icount = 0
    for ds in dsets:
        icount += len(range(ds[0], ds[1]+1))
    if icount != nmax:
        print(dsets)
        sys.exit(f"Total entries = {icount:}")

    return dsets


# ##---------------------------------------------------------------------------
def get_filename(filename):
    """
    Check if file was compressed and/or if it exists.
    Returns exiting file.

    Parameters
    ----------
    filename : Path
        The current assumed file name.

    Returns
    -------
    filename: Path
    The name of the existing file

    """

    # Strangely path.is_file() does not give False but an error
    if not filename.is_file():  # Current file does not exist
        if filename.suffix == ".gz":  # If a gz file, try the non gz file
            filename = filename.with_suffix("")
        else:  # If this not a gz file, try the gz file
            filename = filename.with_suffix(filename.suffix+".gz")
    if not filename.is_file():  # Check that the new file exist
        sys.exit("utilities.get_filename : {filename} not found")

    return filename


# ##---------------------------------------------------------------------------
def file_from_tar(folder=None, tarname=None, target=None):
    """
    Extract a given file from a tar file.

    Parameters
    ----------
    folder : string, optional
        The folder to scan. The default is "None".
    tarname : string, optional
        The archive name in case more than one in the folder.
        The default is None.
    target : string, optional
        The file to be found in the archive, or a file name with the same
        extension. The default is None.

    Returns
    -------
    datafile: pointer to file
        A pointer to a file in memory

    """

    tag = "file_from_tar: "

    # ## ------------------------
    # ## Get archive file name
    # ## ------------------------
    if not Path(folder).is_dir():
        sys.exit(f"{tag} Folder not found")

    if target is None:
        sys.exit("f{tag} Specify a target in the archive")

    if tarname is None:
        # Find the archive in the folder
        p = Path(folder).glob('*.tar.gz')
        files = [x for x in p]  # if x.is_file()] <- Crash when not connected
        if len(files) > 1:
            sys.exit(f"{tag} More than one .tar.gz file, specify a name")
        elif len(files) == 0:
            return None
        else:
            tarname = files[0]
    print(f"{tag} found {tarname}")

    # ## ------------------------
    # ## Open tar file, get members, check data.txt exists
    # ## ------------------------
    tar = tarfile.open(tarname, "r:gz")

    # If the target is not found explicitely, tires a file with same extension
    if target not in [member.name for member in tar.getmembers()]:

        print(f"{tag} No {target} in archive. Tries with extension")

        # find members with correct extension
        extmatch = \
            [(m.name if Path(m.name).suffix == Path(target).suffix else None)
             for m in tar.getmembers()]
        extmatch = list(filter(None, extmatch))  # Remove None

        if len(extmatch) > 1:
            sys.exit(f"{tag} More than one file matching in archive "
                     "(try largest?)")
        elif len(extmatch) == 0:
            sys.exit(f"{tag} No file matching extension in the archive")
        else:
            print(f"{tag} {extmatch[0]} matches {target}")
            target = extmatch[0]

    # At this stage, ether the target was found or deduced from the extension
    print(f"{tag} A file matching {target} was found")
    datafile = tar.extractfile(target)

    return datafile


# ##---------------------------------------------------------------------------
def backup_file(filename, folder=None, dbg=False):
    """
    Copy a file to a result folder.
    If it already exists, make a backup with the date.

    Parameters
    ----------
    file : string
        File name.
    folder : String, optional
        Output folder. The default is None.
    dbg : Boolean, optional
        If True, display messages. The default is False.

    Returns
    -------
    None.

    """

    # Create result folder if not exisitng
    if not os.path.exists(folder):
        if dbg:
            print(" *** Creating output folder ", folder)
        os.makedirs(folder)

    output_file = folder+"/"+filename

    if os.path.exists(output_file):
        nw = datetime.datetime.now()
        output_file = output_file + "_" \
                                  + str(nw.hour) \
                                  + str(nw.minute) \
                                  + str(nw.second)

    shutil.copy(filename, output_file)
    if dbg:
        print("   ----", filename, " copied to ", output_file())


# ##---------------------------------------------------------------------------
def Df(x):
    """
    Returns a time in Julian day format if the argument is finite
    """

    if math.isfinite(x):
        return Time(x, format="mjd")
    return x


# ##---------------------------------------------------------------------------
def Dp(x):
    """
    Returns a time in ISO format if the argument is an astropy
    :class:`Time`, and if not returns the argument (unchanged).
    """
    if isinstance(x, Time):
        return x.iso
    return x
