# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019

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

from   astropy.time import Time

__all__ = ["subset_ids","get_filename","file_from_tar","backup_file",
           "Df", "Dp"]

###----------------------------------------------------------------------------
def subset_ids(nmax, nsets, debug=False):
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
        return [1,nmax]

    delta = int(nmax/nsets) # regular set size
    rest  = nmax%nsets # Last set size

    # Define low and high values of the intervals
    ids_low  = list(range(1,nmax,delta))
    ids_high = list(range(delta,min(nmax+delta,nmax),delta))
    if len(ids_high) < len(ids_low):
        ids_high += [nmax]
    elif ids_high[-1]<nmax:
        ids_high[-1] = nmax

    if debug:
        print("steps = ",delta, " uncovered = ",rest)
        print(ids_low, " -> ",len(ids_low)," items")
        print(ids_high, " -> ",len(ids_high)," items")

    # Create interval list
    dsets = [[idl, idh] for (idl, idh) in zip(ids_low, ids_high)]

    # Check counts
    icount=0
    for ds in dsets:
        icount += len(range(ds[0],ds[1]+1))
    if icount != nmax:
        print(dsets)
        sys.exit(f"Total entries = {icount:}")

    return dsets

###----------------------------------------------------------------------------
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
    The name of the existing file

    """

    # Strangely path.is_file() does not give False but an error
    if not filename.is_file(): # Current file does not exist
        if filename.suffix == ".gz": # If a gz file, try the non gz file
            filename = filename.with_suffix("")
        else: # If this not a gz file, try the gz file
            filename = filename.with_suffix(filename.suffix+".gz")
    if not filename.is_file(): # Check that the new file exist
        sys.exit("utilities.get_filename : {filename} not found")

    return filename

###----------------------------------------------------------------------------
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
    A pointer to a file in memory

    """

    tag = "file_from_tar: "

    ### ------------------------
    ### Get archive file name
    ### ------------------------
    if not Path(folder).is_dir():
        sys.exit(f"{tag} Folder not found")

    if target is None:
        sys.exit("f{tag} Specify a target in the archive")

    if tarname is None:
        # Find the archive in the folder
        p = Path(folder).glob('*.tar.gz')
        files = [x for x in p if x.is_file()]
        if len(files) > 1:
            sys.exit(f"{tag} More than one .tar.gz file, specify a name")
        elif len(files) == 0:
            return None
        else:
            tarname = files[0]
    print(f"{tag} found {tarname}")

    ### ------------------------
    ### Open tar file, get members, check data.txt exists
    ### ------------------------
    tar = tarfile.open(tarname, "r:gz")

    # If the target is not found explicitely, tires a file with same extension
    if not target in [member.name for member in tar.getmembers()]:

        print(f"{tag} No {target} in archive. Tries with extension")

        # find members with correct extension
        extmatch = \
            [(m.name if Path(m.name).suffix==Path(target).suffix else None) \
             for m in tar.getmembers()]
        extmatch = list(filter(None,extmatch)) # Remove None

        if len(extmatch) > 1:
            sys.exit(f"{tag} More than one file matching in archive (try largest?)")
        elif len(extmatch) == 0:
            sys.exit(f"{tag} No file matching extension in the archive")
        else:
            print(f"{tag} {extmatch[0]} matches {target}")
            target = extmatch[0]

    # At this stage, ether the target was found or deduced from the extension
    print(f"{tag} A file matching {target} was found")
    datafile = tar.extractfile(target)

    return datafile

###----------------------------------------------------------------------------
def backup_file(filename,folder=None, dbg=False):

    """
    Copy a file to a result folder.
    If it already exists, make a backup with the date.

    Parameters
    ----------
    file : TYPE
        DESCRIPTION.
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
            print(" *** Creating output folder ",folder)
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
        print("   ----",filename," copied to ",output_file())

###------------------------------------------------------------------------
def Df(x):
    """
    Returns a time in Julian day format if the argument is finite
    """

    if math.isfinite(x):
        return Time(x,format="mjd")
    return x

###------------------------------------------------------------------------
def Dp(x):
    """
    Returns a time in ISO format if the argument is an astropy
    :class:`Time`, and if not returns the argument (unchanged).
    """
    if isinstance(x,Time):
        return x.iso
    return x
