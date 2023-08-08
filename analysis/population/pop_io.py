# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:47:25 2020

@author: Stolar
"""

import sys
sys.path.append("../../")

from pathlib import Path
from   astropy.table import Table

import yaml
from yaml.loader import SafeLoader

from utilities import file_from_tar

__all__ = ["get_config_data", "create_csv"]
###-------------------------------------------------------------
def get_config_data(filename, debug=False):
    """
    Get data from the configuration file, either from opening the file from
    disk or by getting the file from the archive on disk

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
    dirname = Path(filename).parents[0].absolute()

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

###-------------------------------------------------------------
def create_csv(file     = "parameter.yaml",
               dataname = "data.txt",
               debug=False):

    """
    From the current parameter file containing the folder to be analysed,
    create the csv file from default txt file.
    If the parameter file is not given, then the datafile is supposed to
    contain the full path data file name.
    Get the simulation duration.

    Parameters
    ----------
    file : String, optional
        Input parameter file name, can be None. The default is `paramter.yaml`.
    datafile : TYPE, optional
        DESCRIPTION. The default is "data.txt".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    nyears : FLOAT
    csvfilename : STRING

    None.

    """

    nyears = 1

    if file is not None:
        xdict  = yaml.load(open(file), Loader=SafeLoader)

        folders = xdict["outfolders"]
        base    = xdict["base_folder"]
        nyears  = xdict["duration"]

        # Tags are optionall
        if "tags" not in xdict.keys():
            tags = len(folders)*[""]
        else:
            tags    = xdict["tags"]

        txtfilenames = [Path(base,Path(folder,dataname)) for folder in folders]
    else: # If the population name is given, it is unique, not a list
        txtfilenames = [Path(dataname)]
        base=""
        folders = [txtfilenames.parent]

    csvfilenames = []
    for txtfilename, folder in zip(txtfilenames, folders):
        csvfilename = txtfilename.with_suffix('.csv')
        if debug:
            print(" Full name :",txtfilename)
            print(" >>> ",csvfilename)
        ### Check existence of csv file
        if csvfilename.is_file():
            csvfilenames.append(csvfilename)
            if debug:
                print(csvfilename," exists")
        else: ### No file found.
            # Try to create from existing text file otherwise check the archive
            print("csvfile not found - should be created")
            if txtfilename.is_file():
                print("Text file found")
                datafile = txtfilename.as_posix()
            else: ### Extract data from the archive
                print(" Text file not found, try to extract from tar.gz")
                datafile = file_from_tar(folder=Path(base,folder), target="data.txt")
            # At this stage, data are in hands and can be converted
            print("Converting data...")
            data = Table.read(datafile,format="ascii",guess=False)
            data.write(csvfilename.as_posix(), format="ascii.csv",overwrite=True)
            print(csvfilename," Created")
            csvfilenames.append(csvfilename)

    if len(csvfilenames) > 1:
        return nyears, csvfilenames, tags
    else:
        return nyears[0], csvfilenames[0], tags[0]