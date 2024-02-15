# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:47:25 2020

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

import yaml
from yaml.loader import SafeLoader

from astropy.table import Table
from utilities import file_from_tar
from niceprint import heading

sys.path.append("../../")

__all__ = ["get_config_data", "get_data"]

###----------------------------------------------------------------------------
def get_config_data(dirname, debug=False):
    """
    Get data from the configuration file, either from opening the file from
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

    print("Extracting configuration parameters from ",dirname)

    # Get configuration file from the current folder
    for fname  in dirname.iterdir():
        if fname.suffix == ".yaml":
            if debug:
                print(" Found configuration file :",fname)
            cfile_dict = yaml.load(open(fname,"r"), Loader=SafeLoader)
            return cfile_dict

    # If it failed, try to get it from the archive
    print(" No configuration file found in",dirname,". Try archive")
    cfile = file_from_tar(folder  = dirname,
                          tarname = None,
                          target  = "config.yaml")

    cfile_dict = yaml.load(cfile.read(), Loader=SafeLoader)

    return cfile_dict

###----------------------------------------------------------------------------
def check_and_convert(folder, target, debug=False):
    """
    Check whether the current folder contains a data file (.txt) or a
    converted file (.csv). If not, try to extract from a tar file in that
    folder.
    Convert the text file if needed
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
    txtfile = Path(folder,target)
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
        data = Table.read(txtfile.as_posix(),format="ascii",guess=False)
        data.write(csvfile.as_posix(), format="ascii.csv",overwrite=True)
        if debug:
            print(csvfile," Created")
        return csvfile

    # No .csv nor .txt file : last chance in the tar file
    if debug:
        print(txtfile," not found, try to extract from tar.gz")
    datafile = file_from_tar(folder=folder, target=target) # ExFileObject

    # Nothing worked!
    if datafile is None:
        if debug:
            print("Failed")
        return None

    if debug:
        print("data from tar ",datafile.name)
    data = Table.read(datafile,format="ascii",guess=False)
    data.write(csvfile, format="ascii.csv",overwrite=True)

    return csvfile

###----------------------------------------------------------------------------
def get_data(parpath  = None,
             dataname = "data.txt",
             debug    = False):

    """
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
    datafile : string, optional
        Data file name. The default is "data.txt".
    debug : Boolean, optional
        If True, display processing information. The default is False.

    Returns
    -------
    nyears : FLOAT
        Number of years of simulation
    csvfilepaths : list of Path
        List of path for the files found
    None.

    """

    if parpath is None:
        heading(" RUNNING in DEMO mode")
        print(" No `.yaml` file was given")
        # parpath = Path("../../data/samples", "pop_parameter.yaml").resolve()
        # base    = Path("../../data/samples").resolve()
        parpath = Path("../../data/samples", "pop_parameter.yaml")
        base    = Path("../../data/samples")
    else:
        # Check if parameter file exists
        parpath = Path(parpath).resolve() # If teh user forgot
        if not parpath.exists():
            sys.exit(f" The parameter file {parpath:} was not found")

        # Check if HAPPY_OUT was defined, otherwise use the local folder
        if "HAPPY_OUT" not in os.environ.keys():
            base = Path("../../data/samples").resolve()
            heading(" DEMO mode")
        else:
            base = Path(os.environ["HAPPY_OUT"])

    csvfilepaths = []

    # Get the folder names from the parameter file.
    xdict   = yaml.load(open(parpath.as_posix()), Loader=SafeLoader)
    folders = [Path(base,dir) for dir in xdict["outfolders"] ]
    nyears  = xdict["duration"][0]

    # Check if the folders exist
    for curdir in folders:
        if not curdir.exists():
            sys.exit(f"{curdir} does not exist, please correct {parpath:}")

    # Get the .csv file from the folders, either existing or to be created
    for curdir in folders:

        filepath = check_and_convert(curdir, dataname)

        # No file found. Look into the subfolders
        if filepath is None:

            if debug:
                print(" Not found ", dataname, "in ", curdir)
                print(" Scanning the folder")

            dirnames = [f for f in curdir.iterdir() if f.is_dir()]
            for newdir in dirnames:
                filepath = check_and_convert(newdir, dataname)
                if filepath is not None:
                    if debug:
                        print(" Found ", dataname, " in ",newdir)
                    # Found one, add to the list
                    csvfilepaths.append(filepath)

        # File found, add to the list
        else:
            csvfilepaths.append(filepath)

    # Final folder list with data
    if debug:
        print("--------------")
        print("Data found")
        for curfile in csvfilepaths:
            print(curfile)

    # Create defaults tags or get them from the file
    if "tags" not in xdict.keys():
        tags = csvfilepaths[0].parent.parent.name # Migght be long!
    else:
        tags = xdict["tags"]

    return nyears, csvfilepaths, tags

###############################################################################
if __name__ == "__main__":

    get_data(parpath=None, debug=True) # Demo mode
    # get_data(parpath="parameter.yaml", debug=True) # Local file
    print("... completed")
