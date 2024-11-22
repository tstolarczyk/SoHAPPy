# -*- coding: utf-8 -*-
"""
This module is used to generate visibilities of sources appearing at a
certain rate during a given period. It can also generate random positions
in the sky. Alternatively, positions and dates can be read from exiting files.
It contains a :class:`Skies` class and a main function with the following
usage:

The parameters are initialised with the `Skies` constructor which is called
from the generate_sky or sky_from_source function if randomly generated or
read from exisiting data files respectively.

Created on Tue Feb 21 13:16:50 2023

@author: Stolar
"""
import sys
import os
from pathlib import Path

import numpy as np
import json
import yaml
from yaml.loader import SafeLoader
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib import cm

import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord

from visibility import Visibility
import observatory as obs
from configuration import Configuration

from niceprint import heading, warning, highlight, failure
from niceplot import MyLabel, draw_sphere

__all__ = ["Skies"]


###############################################################################
class Skies():
    """
    In what follows, the terminology "sky" refers to positions and explosion
    dates while visibility gives the source visibility from these parameters
    under defined circumstances (e.g: minimal altitude, Moon ligth etc).

    This class handles the parameters and function to generate visibilities
    from a number of sources within a certain time period defined by 2 years
    (from beginning of first year to the end of last year).
    The output files are written in a `folder/visibility` subfolder where
    folder as a fixed naming convention:

        ``keyword_year1_Nyears_version``

    * ``keyword`` is one of the dictionnary keyword of the `visibility.yaml`
      file or any other keyword (line `forced` or `permanent`).
    * ``year1`` is the first year;
    * ``Nyears`` the number of years;
    * ``version`` is a version tag to differentiate visibilities generated
      from different seeds.

    In this folder two kind of files can be created:

    * `yaml` files containing the generated dates and positions. These files
      starts with `DP` (for Dates and Position). They can be used to
      generate other visibilities from the same dates and position.
    * `json` file containing the visibility class content of all sources
      in the identifier range.

    Both files reproduce the name of the folder above with the source range
    added (`id1` and `id2` is the source identifier range):

    * Date and position files: ``DP_keyword_year1_Nyears_version_id1_id2.yaml``
    * Visibility class dump: ``keyword_year1_Nyears_version_id1_id2.json``

    To produce the visibility files for data with positon and dates defined,
    use ``sky_from_source`` (with the possibility to change the dates).


    """
    prfx = "ev"  # Prefix before the id number

    # -------------------------------------------------------------------------
    def __init__(self,
                 year1=9999, nyears=1,
                 first=1, Nsrc=1,
                 version="1",
                 duration=3.0,
                 visibility="strictmoonveto",
                 cfg_path=None,
                 output=Path("skygen_vis"),
                 seed=2022,
                 newpos=False,
                 newdate=False,
                 debug=False):
        """
        Create de default object from external parameters.

        Parameters
        ----------
        year1 : integer
            First year. The default if 9999.
        nyears : integer
            Number of years. The default is 1.
        first : integer
            First source identifier. The default is 1.
        Nsrc : integer
            Number of sources. The default is 1.
        version : string, optional
            Visibility version tag, chosen by the user. The default is "1".
        duration : float, optional
            Duration in days on which the visibility is computed. The default
            is 3.0.
        visibility : string, optional
            A dictionnary entry to be found in `visibility.yaml`.
            The default is `strictmoonveto`.
        cfg_path : pathlib Path, optional
            Configuration file name full path. The default is 'congif.yaml'.
        output : pathlib Path, optional
            Output base folder. The default is `Path("./skygen_vis")`.
        seed : integer, optional
            Seed for random generator. The default is 2022.
        newpos: boolean
            If True, if the positions are read from source file, the positions
            are re-generated isotropically on the sky.
        newdate: boolean
            If True, if the data are read from source file, the dates are
            re-generated from the given source range.
        debug : boolean, optional
            If True, print out various information. The default is False.

        Returns
        -------
        None.

        """

        # Initialise default parameters
        self.dbg = debug
        self.version = version

        self.seed = seed
        np.random.seed(self.seed)

        self.id1 = first
        self.id2 = first + Nsrc - 1
        self.Nsrc = Nsrc

        self.year1 = year1
        self.year2 = year1 + nyears - 1

        self.newdate = newdate
        self.newpos = newpos
        self.viskey = visibility
        self.duration = duration

        # Input parameters (backward compatibilty)
        self.cfg = Configuration()
        if cfg_path is None:
            cfg_path = self.cfg.def_conf

        self.config = cfg_path
        self.cfg.read_from_yaml(filename=self.config)

        # Output
        self.basedir = output if output is not None else '.'
        self.out_folder = None

        # Final command line
        self.cmd_line = ""

        # List of computed visibilities
        self.vis_list = None

    # -------------------------------------------------------------------------
    @classmethod
    def command_line(cls):
        """
        Supersede default argument with the command line arguments and generate
        the command line to be further used for Batch submission

        Returns
        -------
        Skies instance
            Instance filled with the command line argument values.

        """

        import argparse
        inst = cls()  # Initialize default

        parser = argparse.ArgumentParser(description="Generate visibilities",
                                         epilog="---")

        parser.add_argument('-y', '--year1',
                            help="First year",
                            default=inst.year1,
                            type=int)

        parser.add_argument('-n', '--nyears',
                            help="Number of years",
                            default=inst.year2 - inst.year1 + 1,
                            type=int)

        parser.add_argument('-f', '--first',
                            help="First source identifier",
                            default=inst.id1,
                            type=int)

        parser.add_argument('-N', '--Nsrc',
                            help="Number of sources",
                            default=inst.Nsrc,
                            type=int)

        parser.add_argument('-v', '--version',
                            help="version number",
                            default=inst.version)

        parser.add_argument('-D', '--days',
                            help="Visibility range in days",
                            default=inst.duration)

        parser.add_argument('-V', '--visibility',
                            help="Visibility keywords",
                            default=inst.viskey)

        parser.add_argument('-c', '--config',
                            help="Configuration file name",
                            default=Path(inst.config)
                            if (inst.config is not None)
                            else Path("data/config_ref.yaml"))

        parser.add_argument('-o', '--output',
                            help="Output base folder (path)",
                            default=Path(inst.basedir)
                            if (inst.basedir is not None)
                            else Path("."))

        parser.add_argument('-s', '--seed',
                            help="Seed from random generator",
                            default=inst.seed,
                            type=int)

        parser.add_argument('--debug',
                            dest='debug',
                            action='store_true',
                            help="Display processing details")

        parser.add_argument('--nodebug',
                            dest='debug',
                            action='store_false',
                            help="Does not display processing details")

        parser.set_defaults(trigger=inst.newdate)
        parser.add_argument('--trigger',
                            dest='trigger',
                            action='store_true',
                            help="(re)generate dates")
        parser.add_argument('--notrigger',
                            dest='trigger',
                            action='store_false',
                            help="Do not generate dates if already existing")

        parser.set_defaults(position=inst.newpos)

        parser.add_argument('--position',
                            dest='position',
                            action='store_true',
                            help="Generate new positions in the sky")
        parser.add_argument('--noposition',
                            dest='position',
                            action='store_false',
                            help="No new sky positions if existing")

        # Decode command line
        vals = parser.parse_args()

        # Fill the class instance
        inst.dbg = vals.debug
        inst.version = vals.version

        inst.seed = vals.seed
        np.random.seed(inst.seed)

        inst.id1 = vals.first
        inst.id2 = vals.first + vals.Nsrc - 1
        inst.Nsrc = vals.Nsrc

        inst.year1 = vals.year1
        inst.year2 = vals.year1 + vals.nyears - 1

        inst.newdate = vals.trigger
        inst.newpos = vals.position
        inst.viskey = str(vals.visibility)
        inst.duration = vals.days

        # Input parameters
        if vals.config is not None:
            inst.config = str(Path(vals.config).resolve())
            vals.config = inst.config

        # Output
        # Note : The difference between resolve and absolute is that absolute()
        # does not replace the symbolically linked (symlink) parts of the path,
        # and it never raises FileNotFoundError . It does not modify the case
        # either.
        if vals.output is not None:
            inst.basedir = str(Path(vals.output).absolute())
            vals.output = inst.basedir

        # Generate command line
        inst.cmd_line = 'python "' + str(Path(__file__)) + '" '

        # It would be more logical to loop over the class content
        for (k, v) in vars(vals).items():

            if k in ["trigger", "position", "debug", "batch"]:
                inst.cmd_line += "--no"+k+" " if v is False else "--"+k+" "
            elif k in ["config", "output"]:
                inst.cmd_line += '--' + k + " " + '"' + str(v) + '"' + " "
            else:
                if v is not None:
                    inst.cmd_line += "--" + k + " " + str(v) + " "

        # Note that dates are float in MJD and ra, dec float in degrees
        inst.ra = np.zeros(inst.Nsrc)
        inst.dec = np.zeros(inst.Nsrc)
        inst.dates = np.zeros(inst.Nsrc)

        return inst

    # -------------------------------------------------------------------------
    def sky_from_source(self, debug=False):
        """
        Get the sky positions and dates from original source files when
        this information exists.

        Returns
        -------
        None.

        """

        heading("Dates and positon from source files")

        # retrieve data input folder
        if "HAPPY_IN" in os.environ.keys():
            infolder = Path(os.environ["HAPPY_IN"])
        else:
            sys.exit("The HAPPY_IN environment variable should be defined")

        found_position = False
        found_trigger = False

        # Get information from data files
        for i, item in enumerate(range(self.id1, self.id2+1)):

            if (self.Nsrc <= 10) or (np.mod(i, 10) == 0):
                print("#", item, " ", end="")

            fname = Path(infolder,
                         self.cfg.data_dir,
                         self.cfg.prefix + str(item) + self.cfg.suffix)
            try:

                if debug:
                    print(f"Accessing {fname:}")

                hdul = fits.open(fname)
                hdr = hdul[0].header
                keys_0 = list(hdul[0].header.keys())

                if "RA" in keys_0 and "DEC" in keys_0:
                    self.ra[i] = hdr['RA']
                    self.dec[i] = hdr['DEC']
                    found_position = True

                if "GRBJD" in keys_0:  # Not in SHORTFITS
                    date = Time(hdr['GRBJD']*u.day,
                                format="jd", scale="utc").mjd
                    found_trigger = True
                elif "GRBTIME" in keys_0:
                    date = Time(hdr['GRBTIME'],
                                format="jd", scale="utc").mjd
                    found_trigger = True
                else:
                    date = 0  # MJD
            except Exception:
                failure(f" SKIPPING - File not found {fname:}\n")
                date = 0

            self.dates[i] = date

            if self.dbg:
                print("Found: ", i, self.ra[i], self.dec[i], self.dates[i])

        year1 = Time(np.min(self.dates),
                     format="mjd", scale="utc").datetime.year
        year2 = Time(np.max(self.dates),
                     format="mjd", scale="utc").datetime.year
        print("\n Year range in source files : ", year1, " -",
              year2, "(used for the output folder if dates kept)")

        # If requested, supersede the dates
        if self.newdate:
            if found_trigger:
                highlight(" Regenerating dates")
            else:
                highlight(" Generating dates (absent from input files")
            self.generate_dates()

        # If requested, supersede the positions
        if self.newpos:
            if found_position:
                highlight(" Regenerating positions")
            else:
                highlight(" Generating positions (absent from input files")
            self.generate_positions()

    # -----------------------------------------------------------------------------
    @classmethod
    def sky_from_yaml(cls, filename, version=None):
        """
        Read dates and postions from an existing "DP" yaml file.
        The dates are stored as mjd in the file and Time object in the
        instance.

        Parameters
        ----------
        filename : pathlib Path
            Input file name.
        version: integer
            New version tag f needed

        Returns
        -------
        Skies instance
            New instance.

        """

        heading(f" Dates and position load from {filename.name:s}")

        infile = open(filename, "r")
        data = yaml.load(infile, Loader=SafeLoader)

        print("File created: ", data["created"])

        year1 = data["start"]
        year2 = data["stop"]

        if version is None:
            version = data["version"]

        inst = cls(year1=year1, Nyears=year2-year1+1,
                   first=data["id1"], nsrc=data["nsrc"],
                   version=version,
                   duration=data["duration"],
                   visibility=data["key"],
                   output=data["basedir"],
                   seed=data["seed"],
                   debug=False)

        for i, item in enumerate(range(inst.id1, inst.id2+1)):
            key = inst.prfx+str(item)
            [inst.dates[i], inst.ra[i], inst.dec[i]] = data[key].split()

        return inst

    # -------------------------------------------------------------------------
    def generate_dates(self):
        """
        Generate dates in MJD from given year intervall, from January 1st of
        the first year at midnight to December 31st of last year at 23:59:59

        Returns
        -------
        None.

        """

        tstart = Time(datetime(self.year1, 1, 1, 0, 0, 0)).mjd
        tstop = Time(datetime(self.year2, 12, 31, 23, 59, 59)).mjd

        days = np.random.random(self.Nsrc)*(tstop-tstart)
        self.dates = tstart + days  # MJD

    # -------------------------------------------------------------------------
    def generate_positions(self):
        """
        Generate random uniform (ra,dec) positons in the sky in degrees.

        Returns
        -------
        None.

        """

        self.ra = 360*np.random.random(self.Nsrc)
        self.dec = np.arcsin(2*np.random.random(self.Nsrc) - 1)
        self.dec = self.dec*180/np.pi

    # -------------------------------------------------------------------------
    def generate_sky(self):
        """
        Generate dates and position from the seed.

        Returns
        -------
        None.

        """
        heading("Dates and positon - random")

        self.generate_positions()  # RA, DEC in degrees
        self.generate_dates()      # Dates

    # -------------------------------------------------------------------------
    def create_output_folder(self):
        """
        Create the folder containing the corresponding files

        Returns
        -------
        None.

        """
        heading("Create output folder")

        prfx = self.prefix()
        self.out_folder = Path(self.basedir,
                               Path(self.cfg.data_dir).stem,
                               self.cfg.out_dir,
                               prfx)

        self.basename = prfx + "_" + str(self.id1)

        if self.id2 > self.id1:
            self.basename += "_" + str(self.id2)

        # Check if folder exists, otherwise create it
        if not self.out_folder.exists():
            warning("Creating {}".format(self.out_folder))
            try:
                Path.mkdir(self.out_folder, parents=True)
                print(f"Created: {self.out_folder:}")
            except Exception:
                sys.exit(f"{__name__:}.py: Could not create {self.folder:}")
        else:
            warning(f"{self.out_folder:} Already exists")

    # -------------------------------------------------------------------------
    def create_vis(self, paramfile="visibility.yaml", observatory="CTA"):
        """
        Compute visibilities, store in an array.

        Returns
        -------
        None

        """
        heading("Creating visibilities")

        # Check that dates and position have been generated
        if self.ra.all == 0:
            sys.exit("Dates and position are missing")

        # Get visibility paramters from default file
        param = (True, Visibility.params_from_key(self.viskey,
                                                  parfile=paramfile))[1]

        vislist = []

        # Loop over items
        for i, item in enumerate(range(self.id1, self.id2+1)):

            if (self.Nsrc <= 10) or (np.mod(i, 10) == 0):
                print("#", item, " ", end="")

            # print(item, self.ra[i], self.dec[i], self.dates[i])
            radec = SkyCoord(self.ra[i]*u.deg, self.dec[i]*u.deg, frame='icrs')
            tvis1 = Time(self.dates[i], format="mjd", scale="utc")
            tvis2 = tvis1 + self.duration

            for loc in ["North", "South"]:

                vis = Visibility(pos=radec,
                                 site=obs.xyz[observatory][loc],
                                 window=[tvis1, tvis2],
                                 name=str(item)+"_"+loc,
                                 status="")
                vis.compute(param=param)

                if self.dbg:
                    vis.print()

                vislist.append(vis)
        print(" - Done")
        self.vis_list = np.array(vislist)

    # -------------------------------------------------------------------------
    def sky_to_yaml(self, version=None):
        """
        Dump dates and sky position into a yaml file for further use.
        posiitons are stored as they are (float) whereas dates are stored as
        modified Julian Days from the original Time object.

        Returns
        -------
        filename: pathlib Path
            output filename

        """

        self.create_output_folder()

        filename = Path(self.out_folder, "DP_" + self.basename+".yaml")

        with open(filename, "w") as out:

            heading(" Dumping generated dates and posiiton")
            print(f" Output: {filename}")

            print(f"created: {datetime.now()}", file=out)
            print(f"id1: {self.id1}", file=out)
            print(f"nsrc: {self.Nsrc}", file=out)
            print(f"seed: {self.seed:d}", file=out)
            print(f"start: {self.year1}", file=out)
            print(f"stop: {self.year2}", file=out)

            if self.config is not None:
                print(f"config: {str(Path(self.config).parent.parent)}",
                      file=out)
            else:
                print(f"config: {self.config}", file=out)

            print(f"basedir: {str(self.out_folder.parent.parent)}", file=out)
            print(f"key: {self.viskey}", file=out)
            print(f"duration: {self.duration}", file=out)
            print(f"version: {self.version:s}", file=out)

            for i, item in enumerate(range(self.id1, self.id2+1)):
                date = self.dates[i]
                dstr = Time(self.dates[i], format="mjd", scale="utc").isot
                ra = self.ra[i]
                dec = self.dec[i]
                print(f"ev{item:d}: {date:20.10f} {ra:20f} {dec:20f} # {dstr}",
                      file=out)

            out.close()
            print("Done!")

        return filename

    # -------------------------------------------------------------------------
    def prefix(self, version=None):
        """
        This create the prefix string of both yaml (DP) and json (visibility)
        files.

        Parameters
        ----------
        version : string, optional
            A version tag from the user. The default is None.

        Returns
        -------
        string
            The name (no extension) of the yaml and json files.

        """

        if version is None:
            version = self.version
        if version.isnumeric():
            version = "v"+version

        return self.viskey + "_"\
            + str(self.year1) + "_"\
            + str(self.year2 - self.year1+1) + "_"\
            + str(version)

    # -------------------------------------------------------------------------
    def vis2json(self):
        """
        Dump the list of visibility instances to a json file

        Returns
        -------
        None.

        """

        if not self.out_folder.is_dir():
            self.create_output_folder()

        filename = Path(self.out_folder, "VS_"+self.basename+".json")
        heading(" Dumping generated visibilities ")
        print("Output:", filename)

        with open(filename, "w") as f_all:
            json.dump({v.name: v for v in self.vis_list}, f_all,
                      default=Visibility.object_to_serializable, indent=None)
            f_all.close()

    # -------------------------------------------------------------------------
    def plot_sky(self, nbin=25):
        """
        Plot result of the generation of dates and positions

        Parameters
        ----------
        nbin : integer, optional
            Number of bins in histograms. The default is 25.

        Returns
        -------
        None.

        """

        # Generated dates
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 5))
        dt = Time(self.dates, format="mjd", scale="utc")

        ax.hist(dt.datetime, bins=self.nsrc, alpha=0.8,
                label=MyLabel([t.mjd for t in dt]))
        ax.set_xlabel("Date")
        ax.grid(which="both")
        ax.legend()

        # Generated posiitons - ra and dec
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(20, 5))
        ax1.hist(self.ra, bins=nbin)
        ax1.set_xlabel("RA (°)")

        # Generated posiitons - dec versus ra
        ax2.hist(self.dec, bins=nbin)
        ax2.set_xlabel("DEC (°)")

        fig = plt.figure(figsize=(15, 6))
        # ax = fig.add_subplot(111, projection='aitoff')
        ax = fig.add_subplot(111)
        ra = [Angle(x*u.deg).wrap_at(180*u.deg).value for x in self.ra]
        dec = [Angle(x*u.deg).wrap_at(180*u.deg).value for x in self.dec]
        ax.scatter(ra, dec, s=5)
        ax.grid("both")
        ax.set_xlabel("ra (°)")
        ax.set_ylabel("dec (°)")

        # Transform to cartesian coordinates
        radius = 1
        x = radius*np.cos(self.ra*np.pi/180)*np.cos(self.dec*np.pi/180)
        y = radius*np.sin(self.ra*np.pi/180)*np.cos(self.dec*np.pi/180)
        z = radius*np.sin(self.dec*np.pi/180)

        # Check 2D projections
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3,
                                            figsize=(15, 5), sharey=True)
        ax1.plot(x, y, ls="", marker=".", markersize=10, alpha=0.5)
        c1 = plt.Circle((0, 0), radius, fill=False, color="red")
        ax1.add_patch(c1)

        ax2.plot(x, z, ls="", marker=".", markersize=10, alpha=0.5)
        c2 = plt.Circle((0, 0), radius, fill=False, color="red")
        ax2.add_patch(c2)

        ax3.plot(y, z, ls="", marker=".", markersize=10, alpha=0.5)
        c3 = plt.Circle((0, 0), radius, fill=False, color="red")
        ax3.add_patch(c3)

        plt.tight_layout()

        # Generated posiitons - on the sphere
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zlabel("z")
        draw_sphere(radius=0.99*radius, colormap=plt.cm.viridis,
                    ax=ax, alpha=0.1)

        color = cm.cool((z-min(z))/(max(z) - min(z)))

        ax.scatter(x, y, z, color=color, alpha=0.2, s=50)
        ax.view_init(40, 140)

    # -------------------------------------------------------------------------
    def print(self):
        """
        Print class contents with comments

        Returns
        -------
        None.

        """
        print(50*"-")
        print(" Generation number (version): ", self.version)
        print(" Source identifiers from ", self.id1, " to ", self.id2,
              " (", self.Nsrc, ")")
        print(" Generated dates ")
        print("  - from               : ", self.year1)
        print("  - to                 : ", self.year2)
        print("  - Duration (yr)      : ", (self.year2 - self.year1)+1)
        if self.config is not None:
            print(" Reading dates and positions from source files")
            print("  - Configuration file  : ", self.config)
            print("  - New dates           : ", self.newdate)
            print("  - New positions       : ", self.newpos)
        print(" Visibility:")
        print("  - Visibility keyword : ", self.viskey)
        print("  -            range   : ", self.duration)
        print("  - Output folder      : ", self.basedir)
        print(" Debugging : ", self.dbg)
        print()
        print("Command line:")
        print(self.cmd_line)
        print(50*"-")


###############################################################################
if __name__ == "__main__":

    # A standalone function to generate and visibility information from trigger
    # times and sky position either generated randomly or obtained from data
    # in source files.

    import time

    # If no command line, use examples - useful for debugging
    if len(sys.argv[1:]) == 0:
        heading("Running examples")
        # Define command line arguments
        sys.argv = ["skygen.py", "-h"]
        sys.argv = ["skygen.py", "-y", "2004", "-n", "10", "-f", "8",
                    "-N", "5", "-V", "nomoonveto"]
        sys.argv = [r"J:\My Documents\CTA_Analysis\GRB paper\SoHAPPy\skygen.py",
                    "-f", "1",
                    "-N", "3",
                    "-v", "default",
                    "-c", r".\data\config_ref.yaml"]

        sys.argv = ["skygen.py",
                    "-f", "1",
                    "-N", "5",
                    "-v", "tests",
                    "-c", r"data/config_ref.yaml",
                    "--notrigger", "--debug"]
        # ## If no configuration file is given, generates ex-nihilo
        # sys.argv = ["skygen.py",
        #         "-f", "1",
        #         "-N", "5",
        #         "-y","2018",
        #         "-n","10",
        #         "-v", "default",
        #         "--trigger",
        #         "--debug"]
        # sys.argv = ["skygen.py",
        #             "--year1", "2023", "--nyears", "10",
        #             "--first", "1", "--Nsrc", "46",
        #             "--version", "test",
        #             "--days", "4",
        #             "--visibility", "strictmoonveto",
        #             "--config","config-LongFinalTest-omega.yaml",
        #             "--output","skygen_vis",
        #             "--seed","2022",
        #             "--nodebug",
        #             "--trigger", "--noposition"]

    # Extract command line arguments
    gvis = Skies.command_line()
    gvis.print()

    # Start the process
    start_pop = time.time()   # Start chronometer

    # If not configutation file is given, the dates and positons are
    # computed ex-nihilo
    if gvis.config is None:  # Computed ex-nihilio
        gvis.generate_sky()  # Generate positions and dates

    else:  # This is for backward compatibility only
        gvis.sky_from_source(debug=True)

    gvis.sky_to_yaml()
    gvis.create_vis()
    gvis.vis2json()

    # Stop chronometer
    end_pop = time.time()
    elapsed = end_pop-start_pop

    print("\n-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*")
    print(" Duration      = {:8.2f} (s)".format(elapsed))
    print("  per source   = {:8.2f} (s)".format((elapsed)/gvis.Nsrc))
    print(" ******* End of job - Total time = {:8.2f} min *****"
          .format((end_pop-start_pop)/60))
    print("")
    print(datetime.now())
