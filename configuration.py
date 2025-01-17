# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:17:27 2020

This module is organised around the :class:`Configuration` class.

@author: Stolar

"""
import sys
import os
import ast
import json

import argparse
from pathlib import Path
import numpy as np

import yaml
from yaml.loader import SafeLoader

import astropy.units as u

from niceprint import warning, failure, Log
from visibility import Visibility

__all__ = ['Configuration']


###############################################################################
class Configuration():
    """
    SoHAPPy steering parameters and input/output.

    This class handles all the external parameters required to run a simulation
    and the corresponding analysis. It also handles the tools to build the
    input and output folder names, and create the necessary output folders
    if they do not exist.
    """

    ignore = ["filename", "dtslew", "arrays", "dbg"]
    """Keyword to be ignored when dumping the Configuration class instance to
    a file on disk"""

    quantities = ["dtslew_south", "dtslew_north", "dtswift"]
    """ Class members that should be handled as quantities and not strings"""

    def_conf = Path("data", "config_ref.yaml")
    """Default configuration file name"""

    # ##----------------------------------------------------------------------
    def __init__(self):
        """
        Parameter definitions and default values.

        Returns
        -------
        None.

        """
        self.filename = None

        # ## -----------------
        # ## PHYSICS PARAMETERS
        # ## -----------------
        self.ifirst = 1  # id to start with, or a list (then nsrc is not used)
        self.nsrc = 1  # integer, number of sources to be read

        # visibility can be :
        # - "built-in" (read from the data file if it exists);
        # - a subfolder where to find the pre-computed visibility files;
        # - a json file with a colection of visibilities;
        # - a keyword corresponding to a dictionnay in visibility.yaml
        #   to compute the visivility on the fly.
        self.visibility = "strictmoonveto"

        # Possible EBL models are from gammapy or Gilmore data on disk.
        self.ebl_model = "dominguez"
        self.maxnight = None  # Limit data to a maximal number of nights
        self.skip = None  # Number of first nights to be skipped

        # ## -----------------
        # ## INPUT/OUPUT
        # ## -----------------
        self.data_dir = "lightcurves/long_1_1000/"  # Afterglow subfolder
        self.out_dir = "test"  # Result subfolder
        self.prompt_dir = None  # Prompt subfolder (if None, not considered)
        self.irf_dir = "irf/Full/prod3-v2"  # IRF subfolder
        self.datafile = "data.txt"  # Population main output file name
        self.prefix = "Event"  # Prefix file name, followed by a number
        self.suffix = ""  # Suffix file name, following a number
        self.dgt = 0  # Add leading zero to have a coding over "dgt" characters

        # ## -----------------
        # ## SIMULATION PARAMETERS
        # ## -----------------
        self.niter = 1  # Number of Monte Carlo trials
        self.seed = 2021  # Choose ‘random-seed’ to randomize

        # If True Statistical fluctuations are enabled,
        # If False niter forced to one
        self.do_fluctuate = False

        # When True, the simulation is stopped if none of the first
        # trials in the limit of 1 - det_level have reached the minimal
        # significance (3 sigma).
        self.do_accelerate = True

        # Observation position in the time slice
        self.obs_point = "end"

        # ## -----------------
        # ## DETECTION PARAMETERS
        # ## -----------------

        # Fraction of trials above significance threshold to declare a
        # detection at that threshold
        self.det_level = 0.9

        # Number of on over off regions for detection
        self.alpha = 0.2

        self.observatory = "CTAO"
        self.array_north = "FullArray"  # IRF subarrays in North
        self.array_south = "FullArray"  # IRF subarrays in South
        self.dtslew_north = 30*u.s  # Maximum slewing time delay in North
        self.dtslew_south = 30*u.s  # Maximum slewing time delay in South

        self.dtslew = {"North": self.dtslew_north,
                       "South": self.dtslew_south}
        self.arrays = {"North": self.array_north,
                       "South": self.array_south}

        # If True use max. slewing value, otherwise a random delay < dtslew
        self.fixslew = True

        # Alert latency (e.g. SWIFT latency, with amean value of 77 s)
        self.dtswift = 77*u.s

        # If True, use above value. If False, latency generated from Swift data
        self.fixswift = True

        # SWIFT Latency distribution data file in input base folder
        self.swiftfile = "data/swift/Swift_delay_times.txt"

        # ## -----------------
        # ## DEBUGGING / BOOKKEEPING
        # ## -----------------

        # From 0 to 3, increasingly talkative
        # if negative or zero : no plot (show=0)
        # 0: evt counting, 1: + some results, 2: + event details
        # 3: details for each trials
        self.dbg = 1
        self.show = 1
        self.logfile = "analysis.log"  # log file

        self.save_simu = False  # If True, Simulation class saved to file
        self.save_grb = False  # If True, GRB class saved to disk
        self.save_vis = False  # If True, save computed visibility
        self.save_fig = False  # If True, plots saved to pdf file
        self.remove_tar = False  # If True, remove tarred output files
        self.silent = True   # If True, nothing on screen, to log (if dbg=0))

        self.cmd_line = ""     # Command line arguments

        # ## -----------------
        # ## EXPERTS/DEVELOPPERS ONLY
        # ## -----------------
        self.test_prompt = False  # If True test prompt alone (experimental)

        # Prompt characteristics from the afterglow with same id.
        self.use_afterglow = False

        # Shitf in time - Useful for tests
        self.tshift = 0.    # float (days)

        # Force fixed zenith
        self.fixed_zenith = None  # If a value ("20 deg"), freezes zenith

        # Multiplicative factor of the input flux
        self.magnify = 1

        # Store detailed information on slices if True
        self.write_slices = False

        # Limit data energy bins to Elimit
        self.elimit = None

        # Limit lightcurve time bins to tlimit
        self.tlimit = None

        # Analysis minimal and maximal energies - supersede IRF threshold
        self.emin = None
        self.emax = None

        # Observation minimal and maximal times - supersede GRB data
        self.tmin = None
        self.tmax = None

        # If True, use a denser binning for spectral analysis
        self.edense = False

    # ##-----------------------------------------------------------------------
    @classmethod
    def command_line(cls):
        """
        Update an instance from the command line parameters.

        Parameters
        ----------
        debug : Boolean, optional
            If True, talk a bit. The default is False.

        Returns
        -------
        Configuration object
            Current instance.

        """
        inst = cls()  # Initialize default

        # Define the command line arguments - default values will come from
        # the configuration file and are not known at this stage
        # When using getopt, it was possible to treat the reading of
        # the configuraion file separately.
        # As a consequence, the default values from the __init__ constructor
        # are not used here.
        parser = argparse.ArgumentParser(description="SoHAPPy", epilog="---")

        # ifirst can be an integer, a string (either a GRB identifier or a
        # filename) or a list of string or integer. using ast.literal_evl is
        # tricky. So ifirst is considered to be a string and is analysed
        # later when assigned to the ifirst class variable.
        # parser.add_argument('-f', '--first', nargs='+',
        #                     help="First source id, or a list, or a file",
        #                     type=ast.literal_eval,  # Any type
        #                     default="None")

        parser.add_argument('-f', '--first', nargs='+',
                            help="First source id, or a list, or a file",
                            default="None")

        parser.add_argument('-N', '--nsrc',
                            help="Number of source files",
                            default=1,
                            type=int)

        parser.add_argument('-n', '--niter',
                            help="Number of Monte Carlo iteration",
                            type=int)

        parser.add_argument('-c', '--config',
                            default=inst.def_conf,
                            help="Configuration file name")

        parser.add_argument('-m', '--maxnight',
                            help="Maximal number of nights")

        parser.add_argument('-s', '--skip',
                            help="Number of nights to skip")

        parser.add_argument('-V', '--visibility',
                            help="Visibility keyword")

        parser.add_argument('-d', '--debug',
                            help="Debugging level",
                            default=None,
                            type=int)

        parser.set_defaults(save=None)
        parser.add_argument('--save',
                            dest='save',
                            action='store_true',
                            help="Save simulation")
        parser.add_argument('--nosave',
                            dest='save',
                            action='store_false',
                            help="Do not save simulation")

        args, _ = parser.parse_known_args()

        # Find a configuration file, load the data
        # This replaces the constructor defaults
        inst.find_and_read(args.config, debug=args.debug)

        # Supersede parameters if given
        if args.first is not None:
            # If a list is given but it has only one item, get the item
            if len(args.first) != 1:
                inst.ifirst = args.first
            else:
                warning("Unique item list converted to the item")
                inst.ifirst = args.first[0]
                if inst.ifirst.isdigit():
                    inst.ifirst = int(inst.ifirst)
                elif "[" in inst.ifirst and "]" in inst.ifirst:
                    inst.ifirst = ast.literal_eval(inst.ifirst)
        if args.nsrc is not None:
            inst.nsrc = args.nsrc
        if args.niter is not None:
            inst.niter = args.niter
        if args.maxnight is not None:
            inst.maxnight = int(args.maxnight)
        if args.skip is not None:
            inst.skip = int(args.skip)
        if args.visibility is not None:
            inst.visibility = args.visibility
        if args.debug is not None:
            inst.dbg = args.debug
        if args.save is not None:
            inst.save_simu = args.save

        # If prompt_dir is not None, change to Path
        if inst.prompt_dir is not None:
            inst.prompt_dir = Path(inst.prompt_dir)

        # Show plots or not
        inst.show = 0 if inst.dbg < 0 else abs(inst.dbg)

        # If debugging is requested, cannot be silent
        if inst.dbg > 0:
            inst.silent = False

        # If the simulation is saved, it is not fluctuated
        if inst.save_simu:
            inst.do_fluctuate = False

        # If no fluctuation, one simulation is enough !
        if inst.do_fluctuate is False:
            inst.niter = 1

        # Bulld the source list to be processed
        # inst.filelist = inst.source_ids()

        # Extract quantities
        if inst.elimit is not None:
            inst.elimit = u.Quantity(inst.elimit)
        if inst.tlimit is not None:
            inst.tlimit = u.Quantity(inst.tlimit)
        if inst.emin is not None:
            inst.emin = u.Quantity(inst.emin)
        if inst.emax is not None:
            inst.emax = u.Quantity(inst.emax)
        if inst.tmin is not None:
            inst.tmin = u.Quantity(inst.tmin)
        if inst.tmax is not None:
            inst.tmax = u.Quantity(inst.tmax)

        # Redefine array and dtslew variaables
        inst.dtslew = {"North": inst.dtslew_north,
                       "South": inst.dtslew_south}
        inst.arrays = {"North": inst.array_north,
                       "South": inst.array_south}

        # Check contradictions in visibility request
        if inst.visibility in ["permanent", "forced"] \
           and inst.maxnight is not None:
            failure(" Number of nights cannot be limited"
                    " if 'permanent' or 'forced'")
            sys.exit(" Please correct inconsistency in configuration input")

        # Generate command line
        vals = parser.parse_args()
        inst.cmd_line = "python SoHAPPy.py "
        for (k, v) in vars(vals).items():

            # if k in ["debug"]:
            #     inst.cmd_line += "--no"+k+" "  if v is False else "--"+k+" "

            # For reasons that I have not understood, the first source id. is
            # taken as a list even if there is one item
            if isinstance(v, list) and len(v) == 1:
                v = v[0]
            if k in ["config"]:
                inst.cmd_line += '--' + k + " " + '"' + str(v) + '"' + " "
            else:
                if v is not None:
                    inst.cmd_line += "--" + k + " " + str(v) + " "

        return inst

    # ##-----------------------------------------------------------------------
    def find_and_read(self, name=None, debug=True):
        """
        Read steering parameters from a file.

        Check existence of the given configuration file name, and read data if
        it exists.

        Parameters
        ----------
        name : None, String or Path
            A tentative configuration file name.
        debug : boolean, optional
            If True, talk a bit. The default is True.

        Returns
        -------
        None.

        """
        self.filename = None  # Assume not given

        if name is None:  # No argument was given on the command line
            basefolder = Path(__file__).parent  # Where the code is
            self.filename = Path(basefolder, Configuration.def_conf)
            if debug:
                print(" Using default config file")

        else:  # A name was given
            if Path(name).is_file():  # The argument is an existing file
                self.filename = Path(name)
            else:
                sys.exit(f"{__name__:}.py/build: "
                         f"{name:} configuration file does not exist")

        # Read configuration file from file name
        self.read_from_yaml(debug=debug)

    # ##-----------------------------------------------------------------------
    def read_from_yaml(self, filename=None, debug=False):
        """
        Read the configuration data from a `yaml` file on disk.

        Parameters
        ----------
        filename : string, optional
            Configuration file name. The default is None.

        debug : boolean, optional
            If True, talk a bit. The default is True.

        Returns
        -------
        None.

        """
        if filename is not None:
            self.filename = Path(filename)

        # ---------------------------------------------------
        def obj_dic(record):

            # top = type('new', (object,), d)
            seqs = tuple, list, set, frozenset
            for key, val in record.items():

                # Do not read wrong/deprecated keywords
                if key not in self.__dict__:
                    print(f"{__name__:}.py: warning: {key:} ignored")
                    continue

                if isinstance(val, dict):
                    setattr(self, key, obj_dic(val))
                elif isinstance(val, seqs):
                    setattr(self, key,
                            type(val)(obj_dic(sv) if isinstance(sv, dict)
                                      else sv for sv in val))
                else:
                    if key in self.quantities:
                        setattr(self, key, u.Quantity(val))
                    else:
                        setattr(self, key, val)
        # ---------------------------------------------------

        if debug:
            print(">>> Read configuration from ", self.filename)

        try:
            file = open(self.filename, "r", encoding="utf-8")
        except IOError:
            sys.exit(f"{__name__:}.py : {self.filename:} does not exist.")

        data = yaml.load(file, Loader=SafeLoader)
        obj_dic(data)

        # Bulld the source list to be processed
        # self.filelist = self.source_ids()

    # ##-----------------------------------------------------------------------
    def write(self, out_name=None):
        """
        Write current configuration to a `yaml` file for further use.

        Note that this does not copy the possible comments of the original
        file.

        Parameters
        ----------
        out_name : string, optional
            Configuration file name. The default is None.

        Returns
        -------
        None.

        """
        if out_name is None:
            out_name = Path("config_backup.yaml")

        # Create a dict with unnecessary keywords dropped and no quantities
        newdict = {}
        for key, val in self.__dict__.items():

            if key in self.ignore:
                continue
            if isinstance(val, (u.Quantity, Path)):
                newdict[key] = str(val)
            else:
                newdict[key] = val

        with open(out_name, 'w', encoding="utf-8") as file:
            yaml.dump(newdict, file, sort_keys=False)

    # ##-----------------------------------------------------------------------
    def print(self, out=None):
        """
        Print configuration class instance to screen and 'out' file.

        Parameters
        ----------
        out : Log class instance
            Logbook output.

        Returns
        -------
        None.

        """
        if out is None:
            out = Log()

        # ----------------------------------------------------
        def title(txt):
            out.prt("")
            out.prt(f"=== {txt:<27s} {40*'=':40s}")
        # ----------------------------------------------------

        out.prt(72*"=")
        out.prt(f" Configuration file*     : {self.filename:} ")

        # ## -----------------
        # ## PHYSICS PARAMETERS
        # ## -----------------
        title("Physics parameters")
        if not isinstance(self.ifirst, list) \
           and not isinstance(self.ifirst, str):
            out.prt(f" Number of sources*        : {self.nsrc:>5d}")
            out.prt(f" First source*             : {self.ifirst:>5d}")
        else:
            out.prt(f" Source list*              : {self.ifirst:}")
        out.prt(f" Visibility*               : {self.visibility:}")
        out.prt(f" EBL model                 : {self.ebl_model:}")
        nnight = self.maxnight if self.maxnight is not None else "Inf"
        out.prt(f" Number of nights limit*   : {nnight}")
        nskip = self.skip if self.skip is not None else "0"
        out.prt(f" Night(s) skipped*         : {nskip}")
        out.prt(f" Energy max. limit         : {self.elimit:}", end=" ")
        if self.elimit is not None:
            out.warning(" !")
        else:
            out.prt("")

        out.prt(f" Time max. limit           : {self.tlimit:}", end=" ")
        if self.tlimit is not None:
            out.warning(" ! ")
        else:
            out.prt("")

        # ## -----------------
        # ## INPUT/OUPUT
        # ## -----------------
        title("Input/output")
        out.prt(f" Data subfolder            : {self.data_dir:}")
        out.prt(f" IRF subfolder             : {self.irf_dir}")

        if self.prompt_dir is not None:
            out.prt(f" Prompt data subfolder     : {self.prompt_dir:}")
        else:
            out.prt(" Prompt data               : not considered")
        out.prt(f" Output population         : {self.datafile:}")
        out.prt(f" Source file prefix        : {self.prefix:}")
        out.prt(f"             suffix        : {self.suffix:}")

        # ## -----------------
        # ## SIMULATION PARAMETERS
        # ## -----------------
        title("Simulation")
        if not self.niter:
            sys.exit(f"{__name__:}.py: At least 1 trial is requested")
        else:
            out.prt(f" Number of trials*          : {self.niter:>5d}")
        out.prt(f" Seed                       : {self.seed:}")

        out.prt(f" Counts randomised          : {self.do_fluctuate:} ", end="")
        if not self.do_fluctuate:
            out.warning("-> Only 1 trial")
        else:
            out.prt("")

        out.prt(f" Stop if cannot be detected : {self.do_accelerate:}")

        # ## -----------------
        # ## DETECTION PARAMETERS
        # ## -----------------
        title("Detection")
        out.prt(f" Detection level            : {self.det_level}")
        out.prt(f" Alpha (1/n)                : {self.alpha}")
        out.prt(f" Observatory                : {self.observatory}")
        out.prt(f" Site sub-arrays            : N:{self.array_north:} "
                f"S:{self.array_south:}")
        out.prt(f" Slewing time               : N:{self.dtslew_north:}"
                f" S:{self.dtslew_south:}")
        out.prt(f"    Fixed                   : {self.fixslew}")
        out.prt(f" SWIFT delay fixed          : {self.fixswift}")
        if self.fixswift:
            out.prt(f"    Value                   : {self.dtswift}")
        else:
            out.prt(f"            Read from : {self.swiftfile}")

        # ## -----------------
        # ## DEBUGGING / BOOKKEEPING
        # ## -----------------
        title("Debugging / bookkeeping")
        out.prt(f" Debug mode*                : {self.dbg:>5d}")
        out.prt(f" out file                   : {self.logfile}")
        out.prt(f" Show plots                 : {self.show:>5d}")
        out.prt(f" Save simulation            : {self.save_simu} ", end="")
        if self.save_simu:
            out.warning("-> Only 1 trial")
        else:
            out.prt("")
        out.prt(f" Save computed visibility   : {self.save_vis}")
        out.prt(f" Save source class          : {self.save_grb}")
        out.prt(f" Save figures to pdf        : {self.save_fig}")
        out.prt(f" Remove tarred file         : {self.remove_tar}")

        out.prt("+"+66*"-"+"+")
        out.prt(" *: can be changed with command line (use -h)")
        out.prt("+"+66*"-"+"+")
        out.prt(" Developments:")

        # ## -----------------
        # ## EXPERTS/DEVELOPPERS ONLY
        # ## -----------------
        if self.tlimit is not None:
            out.warning(f"Light curve time samples limited to {self.tlimit:}")
        if self.elimit is not None:
            out.warning(f"Energy samples limited to {self.elimit:}")
        if self.emin is not None:
            out.warning(f"Mimimal energy for analysis is {self.emin:}")
        if self.emax is not None:
            out.warning(f"Maximal energy for analysis is {self.emax:}")
        if self.edense is True:
            out.warning("Dense energy sampling for spectral analysis")
        if self.tmin is not None:
            out.warning(f"Mimimal observation time is {self.tmin:}")
        if self.tmax is not None:
            out.warning(f"Maximl observation time is {self.tmax:}")
        if self.test_prompt:
            out.warning(f"{'Test prompt simulation':60s}")
            if self.use_afterglow:
                out.warning(f"{'-> use afterglow information':60s}")
        if self.fixed_zenith is not None:
            out.warning(f"Zenith angle fixed at value '{self.fixed_zenith}'")
        if self.magnify != 1:
            out.warning(f"GRB flux values are multiplied by {self.magnify}")
        if self.tshift != 0:
            out.warning(f"Times shifted by {self.tshift} days")
        if not self.silent:
            out.warning("Maximise screen output")

        if self.write_slices is True:
            out.highlight("Slice information disk-saved (write_slices=True)")

        out.prt("    --end\n")

    # ##-----------------------------------------------------------------------
    def create_output_folder(self, resfolder):
        """
        Create the code outptut folder.

        The convention is that the output folder refers to the folder name
        containing the lightcurves, a subfolder refering to the visibility
        which will be used.

        Parameters
        ----------
        resfolder : String
            Output file name

        Returns
        -------
        resdir : pathlib Path
            output folder.

        """
        # Classical case of a series of integer from `ifirst` and `nsrc`
        if isinstance(self.ifirst, int):
            ext = "_" + str(self.ifirst).zfill(self.dgt)\
                + "_"+str(self.ifirst+self.nsrc-1).zfill(self.dgt)
        # ifirt is a string, and possibly a filenmae
        elif isinstance(self.ifirst, str):
            # If this is an existing file, its name (without path nor
            # extension) is added to the lowest subfolder
            if Path(self.ifirst).exists():
                ext = "_" + str(Path(self.ifirst).stem)
            else:
                ext = "_" + str(self.ifirst)
        else:
            ext = ""

        resdir = Path(resfolder,
                      Path(self.data_dir).name,
                      self.out_dir,
                      self.visibility,
                      self.visibility+ext).resolve()

        # Check that the output folder exists, otherwise try to create it
        if not resdir.exists():
            warning(f"Creating {resdir}")
            try:
                Path.mkdir(resdir, parents=True)
            except AssertionError:
                sys.exit(f"{__name__:}.py: Could not create {resdir:}")
        else:
            warning(f" Already exists :{resdir:}")

        return resdir

    # ##-----------------------------------------------------------------------
    def source_ids(self, infolder, bunchsize=5000, debug=False):
        """
        Obtain the data file list from the input parameters.

        In the configuration file, 'first' can be:
            - an integer, giving the file identifier (number) in a data set.
              The list is build up to 'nsrc' which is also given in the
              configuration file.
            - a list of items, including integers and strings. The integers are
              the file identifiers in a data set, whereas the string refers
              to a GRB in the local 'data/samples/historical' folder.
            - a json file absolute path containing the relative path file list
              with respect to the HAPPY_IN variable (The json files can be
              obtained for instance from the
              'analysis/population/detectable.py' module).
        When the given information is an integer, the file name is build from
        the 'prefix', 'suffix' and 'data_dir' information in the configuration
        file.
        The file is searched for in the 'data_dir' subfolder, first in
        subfolders corresponding to a larget set chunked in 'bunchsize'
        subpopulation and then directly in that 'data_dir'.

        Returns
        -------
        A list of filenames to be simulated and analysed.

        """
        srclist = None

        if not isinstance(self.ifirst, list):
            # Start with the possibility to have a file, as a filename is also
            # a string.
            if Path(str(self.ifirst)).exists():

                with open(str(self.ifirst), encoding="utf-8") as fdata:
                    srclist = json.load(fdata)  # Complete file name
                    self.nsrc = len(srclist)
            # It is a string but not a filename (historical file)
            elif isinstance(self.ifirst, str):
                srclist = [self.ifirst]
                self.nsrc = 1

            # It is an integer, and lead to a list of integers
            elif isinstance(self.ifirst, int):
                srclist = list(range(self.ifirst, self.ifirst+self.nsrc))
        else:
            # It is a list
            srclist = self.ifirst
            self.nsrc = len(srclist)

        # Build the file paths in-situ from the list of identifiers
        filelist = []
        fname = None

        for src in srclist:

            # from a number as an identifier
            if isinstance(src, int):

                name = self.prefix + str(src).zfill(self.dgt) + self.suffix

                # First form the name with large population subfolders
                # Large population have priority

                fmin = str(int((src-1)/bunchsize)*bunchsize + 1)
                fmax = str((int((src-1)/bunchsize) + 1)*bunchsize)
                folder = fmin.zfill(self.dgt) + "_" + fmax.zfill(self.dgt)
                fname = Path(infolder, self.data_dir, folder, name)
                if debug:
                    print("Try to find within bunched of ",
                          bunchsize, ": ", fname)
                    print(fname)

                if not fname.exists():
                    # Same name without subfolder
                    fname = Path(infolder, self.data_dir, name)

                    # We tried, but we failed
                    if not fname.exists():
                        sys.exit(" configuration/source_ids : "
                                 "file not found even in bunched subfolders")

            # From a string giving a relative path wrt HAPPY_IN
            elif Path(infolder, src).exists():
                fname = Path(infolder, src)

            # From a string pointing to an historical GRB yaml file
            elif isinstance(src, str):  # Historical GRB
                fname = Path(Path(__file__).absolute().parent,
                             "data/historical/GRB_"+src+".yaml")

            # not understood, something wrong
            else:
                sys.exit(" configuration/source_ids : input encoding failed")

            filelist.append(fname.as_posix())

        print(f" Number of files to be processed : {len(filelist):}")

        return np.array(filelist)

    # ##-----------------------------------------------------------------------
    def decode_visibility_keyword(self, folder=None):
        """
        Decode the `visibility` keyword value.

        Return an information
        to be handled at running time. The visibility keyword can be:

        * a visibility folder that should contain a suitable visibility `json`
          file where all visibility instances of a data subset are stored.
          By definition the folder structure is the following:

            * output
            * population name
            * a subfolder name as described in skygen,
              i.e `keyword_year1_Nyears_version`

        * `forced` to force a single infinite night
        * `permanent` to force a permanent visibility
        * `built-in` : will use the default visibility stored in the source
          files, if any (for backward compatibility).
        * a keyword giving access to information in a `yaml` file, as a
          python dictionnary.

        Parameters
        ----------
        debug : boolean, optional
            If True, talk a bit. The default is False.

        Returns
        -------
        a visibility information
            The information necessary to get the visibility at running time.
            Various type of information in a single variable. It can be a
            string or a dictionnary.
        """
        # Check predefined possibilities
        if self.visibility in ["built-in", "forced", "permanent"]:
            print(f" Visibilities from keyword: '{self.visibility}' ")
            return self.visibility

        # Check if referenced in visibility.yaml
        vispar = Visibility.params_from_key(self.visibility)
        if vispar is not None:
            return vispar

        # Check if the keyword corresponds to a visibility subfolder
        if folder is None:
            sys.exit("configuration/secode_visibility_keyword:"
                     " Please give a folder name")

        test = Path(folder.parent)
        print(f"Accessing visibility in : {test:}")

        if not test.is_dir():
            sys.exit(f"{__name__}.py : visibility Keyword supposed to be "
                     "a folder and is not")

        # A folder is given
        # Find file name pattern from source identifiers
        idmin = self.ifirst
        idmax = self.ifirst + self.nsrc - 1

        # Complete version that find at least a valid file for the source
        # id range

        candidates = list(test.glob("*.json"))
        visfile = None

        for fname in candidates:
            # Get ids in the filename
            bnds = [int(s) for s in fname.stem.split("_") if s.isdigit()]

            if len(bnds) == 4:  # Range
                bndmin = bnds[-2]
                bndmax = bnds[-1]
            elif len(bnds) == 3:   # One source only
                bndmin = bndmax = bnds[-1]
            # print(fname.name,": ",bnds,bndmin, bndmax )

            if bndmin <= idmin and idmax <= bndmax:
                visfile = fname
                break

        # Load data from the visibility file
        if visfile is not None:
            print(" Visibility data :", visfile)
            with open(visfile, "r", encoding="utf-8") as f:
                return json.load(f)
        else:
            sys.exit(f"{__name__}.py : Visibility: Recomputation "
                     "seems required")

    # ##-----------------------------------------------------------------------
    def get_delay(self):
        """
        Compute the overall delay to be applied to the start of detection.

        Includes the satellite and telescope slewing, according to the user
        parameters.

        Returns
        -------
        delay : Quantity (time)
                Delay before the detection can start.

        """
        delay = {"North": 0*u.s, "South": 0*u.s}

        for loc in ["North", "South"]:

            delta = 0*u.s
            if self.fixslew:
                delta = self.dtslew[loc]
            else:
                delta = self.dtslew[loc]*np.random.random()

            if self.fixswift:
                delta = delta + self.dtswift  # don't do += !!!
            else:
                sys.exit(f"{__name__}.py: Variable SWIFT delay "
                         "not implemented)")
            delay[loc] = delta.to(u.s)

        return delay


###############################################################################
if __name__ == "__main__":

    # A standalone function to read a configuration from a configuration file
    # with wome parameters possibly ovsersed by the command line.
    # * Directly from the file on disk, e.g?

    # os.environ["HAPPY_IN"] = "D:\\CTAO\SoHAPPy\input"
    # os.environ["HAPPY_OUT"] = "D:\\CTAO\SoHAPPy\output"

    log = Log()

    testfile = Path("data/config_ref.yaml")

    # sys.argv = ["", "-f", "343", "-N", "10","-n", "100",
    #             "-V", "nomoonveto"]
    cf1 = Configuration()  # Default
    cf1.read_from_yaml(filename=testfile)
    res_dir = cf1.create_output_folder(Path(os.environ["HAPPY_OUT"]))
    cf1.decode_visibility_keyword(res_dir)
    cf1.print(log)

    # cf1.read_from_yaml(filename="config.yaml")
    # # cf1.print(log)
    # cf1.write(out_name = "config_bck.yaml")

    # cf2 = Configuration()
    # cf2.read_from_yaml(filename="config_bck.yaml")
    # cf2.print(log)

    # cf4 = Configuration.build(sys.argv[1:])
    # cf4.print(log)

    # test automatic output file naming

    # data_sub = Path(cf1.data_dir).name
    # if cf1.trigger ==0:
    #     trigpos_sub = "TP_default"
    # else:
    #     trigpos_sub = "TP_tobedefined"

    # print(" output subfolder for this run :")
    # out_sub = Path(data_sub,trigpos_sub)
    # print(Path(out_sub))
    # print(Path(cf1.resfolder,out_sub))
