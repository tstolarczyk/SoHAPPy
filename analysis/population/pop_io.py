# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:47:25 2020.

Contain the functions to obtain population data from disk, either directly or
from a text file to be converted into a `.csv` file. The text file is in some
cases extracted from a `.tar.gz` archive.

The file or folder names, the simulation duration and the tags associated to
each files are expected to be found in a `.yaml` file. Here is an example.

@author: Stolar
"""
import sys
import os
from pathlib import Path

from itertools import chain

import yaml
from yaml.loader import SafeLoader

from astropy.table import Table
from utilities import file_from_tar
from niceprint import heading

sys.path.append("../../")

__all__ = ["get_config_data", "get_data"]


# ##---------------------------------------------------------------------------
def get_config_data(dirname, debug=False):
    """
    Get data from the configuration file.

    Data are obtained either from opening the file from
    disk or by getting the file from the archive on disk.

    Parameters
    ----------
    filename : string
        The data full file name, used to find back the corresponding
        configuration file.

    Returns
    -------
    cfile_dict : dictionnary
        The content of the configuration file in the form of a python
        dictionnary

    """
    print("Extracting configuration parameters from ", dirname)

    # Get configuration file from the current folder
    for fname in dirname.iterdir():
        if fname.suffix == ".yaml":
            if debug:
                print(" Found configuration file :", fname)
            cfile_dict = yaml.load(open(fname, "r"), Loader=SafeLoader)
            return cfile_dict  # Manages missing file later

    # If it failed, try to get it from the archive
    print(" No configuration file found in", dirname, ". Try archive")
    cfile = file_from_tar(folder=dirname,
                          tarname=None,
                          target="config.yaml")

    cfile_dict = yaml.load(cfile.read(), Loader=SafeLoader)

    return cfile_dict


# ##---------------------------------------------------------------------------
def check_and_convert(folder, target, debug=False):
    """
    Check whether the current folder contains a data file.

    Th efile can be a text file (.txt) or a converted file (.csv). If not, try
    to extract from a tar file in that folder.
    Convert the text file if needed.
    Returns the converted file (.csv) or None if no file found.


    Parameters
    ----------
    folder : Path
        Current folder.
    target : string or Path
        File to be found
    debug: boolean, optional
        If `True`, gives details. The default is `False`

    Returns
    -------
    Path or None
        Converted file if found or created.

    """
    datafile = None

    # Build the text and cvs file Path
    txtfile = Path(folder, target)
    csvfile = Path(folder, target).with_suffix(".csv")

    # Check if converted target exists
    # if csvfile.exists(): # Problematic when not connected on the server
    # Since it convert the variable to an absolute path
    if csvfile.name in os.listdir(csvfile.parent):
        if debug:
            print(" Found ", csvfile, " already converted")
        return csvfile

    # If not csv file, try to convert txt file
    # if txtfile.exists(): # See above
    if txtfile.name in os.listdir(txtfile.parent):

        if debug:
            print(txtfile, " exists, convert it")
        data = Table.read(txtfile.as_posix(), format="ascii", guess=False)
        data.write(csvfile.as_posix(), format="ascii.csv", overwrite=True)
        if debug:
            print(csvfile, " Created")
        return csvfile

    # No .csv nor .txt file : last chance in the tar file
    if debug:
        print(txtfile, " not found, try to extract from tar.gz")
    datafile = file_from_tar(folder=folder, target=target)  # ExFileObject

    # Nothing worked!
    if datafile is None:
        if debug:
            print("Failed")
        return None

    if debug:
        print("data from tar ", datafile.name)
    data = Table.read(datafile, format="ascii", guess=False)
    data.write(csvfile, format="ascii.csv", overwrite=True)

    return csvfile


# ##---------------------------------------------------------------------------
def get_data_from_folder(folder, dataname="data.txt", debug=False):
    """
    Get the .csv file from the folders, either existing or to be created.

    Parameters
    ----------
    folders : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    filelist = []

    filepath = check_and_convert(folder, dataname)

    # No file found. Look into the subfolders
    if filepath is None:

        if debug:
            print(" Not found ", dataname, "in ", folder)
            print(" Scanning the folder")

        dirnames = [f for f in folder.iterdir() if f.is_dir()]
        for newdir in dirnames:
            filepath = check_and_convert(newdir, dataname)
            if filepath is not None:
                if debug:
                    print(" Found ", dataname, " in ", newdir)
                # Found one, add to the list
                filelist.append(filepath)
            else:
                sys.exit(f"No file found in {str(folder):s}")

    # File found, add to the list
    else:
        filelist.append(filepath)

    return filelist


# ##---------------------------------------------------------------------------
def get_data(parpath=None,
             dataname="data.txt",
             debug=False):
    """
    Get data from a file path.

    If `file` does not point to an existing parameter file, then the function
    runs in demo mode, using the samples in the `SoHAPPy/data/samples`
    subfolder.
    From the parameter file containing the folder to be analysed,
    create the csv file from the initial text file.
    If the parameter file is not given, then the `file` is supposed to
    contain the full path data file name.
    Get the simulation duration, and tgas associated to the datafile to be
    displayed on the plots (Otherwise build some from the file names).

    Parameters
    ----------
    file : Path, optional
        Input parameter path, can be None.
    dataname : string, optional
        Data file name. The default is "data.txt".
    debug : Boolean, optional
        If True, display processing information. The default is False.

    Returns
    -------
    nyears : Float
        Number of years of simulation
    csvfilepaths : list of Path
        List of path for the files found
    tags : String
        Filename for later identifications

    """
    # Find the parameter file
    local = False
    if parpath is None:
        heading(" RUNNING in DEMO mode")
        print(" No `.yaml` file was given")
        parpath = Path("../../data/samples", "pop_parameter.yaml")
        local = True
    else:
        parpath = Path(parpath)  # If the user forgot

    if not parpath.exists():
        sys.exit(f" The parameter file {parpath:} was not found")

    # Get the folder names and duration from the parameter file.
    # Add the output foder path if defined
    xdict = yaml.load(open(parpath.as_posix()), Loader=SafeLoader)

    if "HAPPY_OUT" in os.environ and local is False:
        base = os.environ["HAPPY_OUT"]
    else:
        print(" Data supposed in local folder")
        base = Path("../../data/samples")
    folders = [Path(base, folder) for folder in xdict["outfolders"]]
    nyears = xdict["duration"][0]

    # Check if the folders exist
    for curdir in folders:
        if not curdir.exists():
            sys.exit(f"{curdir} does not exist\nPlease correct {parpath:} "
                     "or consider defining the SoHAPPy environment variables.")

    csvfilepaths = []
    for curdir in folders:
        filepath = get_data_from_folder(curdir, dataname=dataname)
        csvfilepaths.append(filepath)

    # Flatten the list
    csvfilepaths = list(chain.from_iterable(csvfilepaths))

    # Final folder list with data
    if debug:
        print("--------------")
        print("Data found in ", len(csvfilepaths), "folders")
        for curfiles in csvfilepaths:
            print(curfiles)

    # Create defaults tags or get them from the file
    if "tags" not in xdict.keys():
        tags = csvfilepaths[0].parent.parent.name  # Might be long!
    else:
        tags = xdict["tags"]

    return nyears, csvfilepaths, tags


###############################################################################
if __name__ == "__main__":

    # parpath = "parameter.yaml"  # Specific series
    # parpath = None
    parpath = "parameter_100k_ISM_alpha.yaml"

    get_data(parpath=parpath, debug=True)  # Local file
    print("... completed")
