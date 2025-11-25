# -*- coding: utf-8 -*-

"""
Created on Thu Dec 12 09:01:41 2019.

@author: Stolar

"""
import gammapy
import sys
import itertools
from pathlib import Path

import numpy as np

import seaborn as sns

import mcsim_config as mcf
from niceplot import single_legend
from niceprint import t_str

import astropy.units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt

from gammapy.maps import MapAxis
if gammapy.__version__ < "1.2":
    from gammapy.irf import load_cta_irfs
else:
    from gammapy.irf import load_irf_dict_from_file


__all__ = ['IRF']


###############################################################################
class IRF():
    """Handle the Instrument Response Function information and utilities."""

    zenith_list = {"20deg": 20*u.deg,
                   "40deg": 40*u.deg,
                   "60deg": 60*u.deg
                   }
    """ List of available zenith angles in the IRF file names and associated
    values."""

    dt_list = {"prod3":
               {"100s": (100*u.s),
                "30m": (30*u.min).to(u.s),
                "5h": (5*u.h).to(u.s),
                "50h": (50*u.h).to(u.s)
                },
               "prod5":
               {"1800s": (1800*u.s),   # 30 min
                "18000s": (18000*u.s),  # 5h
                "180000s": (180000*u.s)
                }  # 50h
               }

    """ List of available time durations in the IRF file names and
    associated values."""

    dtl = {"prod3":
           {"100s": np.log10((100*u.s).value),
            "30m": np.log10((30*u.min).to(u.s).value),
            "5h": np.log10((5*u.h).to(u.s).value),
            "50h": np.log10((50*u.h).to(u.s).value)
            },
           "prod5":
           {"1800s": np.log10((1800*u.s).value),    # 30 min
            "18000s": np.log10((18000*u.s).value),   # 5h
            "180000s": np.log10((180000*u.s).value)  # 50h
            }
           }
    """ List of time intervals in logarithm corresponding to the durations in
    the IRF file names."""

    zenith_valid = {"20deg": [0*u.deg, 33*u.deg],
                    "40deg": [33*u.deg, 54*u.deg],
                    "60deg": [54*u.deg, 90*u.deg]  # Will be limited by altmin
                    }
    """Validity range in zenith of the IRFs (see documentation of this module).
    The `60deg` IRF is allowed to be used down to alitude zero for tests
    Its use is foreseen to be limited by the `altmin` variable."""

    dt_log_valid = \
        {"prod3":
         {
          "100s":
          [0,
           10**(0.5*(dtl["prod3"]["100s"] + dtl["prod3"]["30m"]))
           ],
          "30m":
          [10**(0.5*(dtl["prod3"]["100s"] + dtl["prod3"]["30m"])),
           10**(0.5*(dtl["prod3"]["30m"] + dtl["prod3"]["5h"]))
           ],
          "5h":
          [10**(0.5*(dtl["prod3"]["30m"] + dtl["prod3"]["5h"])),
           10**(0.5*(dtl["prod3"]["5h"] + dtl["prod3"]["50h"]))],
          "50h":
          [10**(0.5*(dtl["prod3"]["5h"] + dtl["prod3"]["50h"])),
           np.inf]
         },
         "prod5":
         {"1800s":
          [0,
           10**(0.5*(dtl["prod5"]["1800s"] + dtl["prod5"]["18000s"]))],
          "18000s":
          [10**(0.5*(dtl["prod5"]["1800s"] + dtl["prod5"]["18000s"])),
           10**(0.5*(dtl["prod5"]["18000s"] + dtl["prod5"]["180000s"]))],
          "180000s":
          [10**(0.5*(dtl["prod5"]["18000s"] + dtl["prod5"]["180000s"])),
           np.inf]
          }
         }
    """ Validity range of IRF in time, taking into account that the validity
    intervals are somehow in logarithmic scale.
    The edge values are the following :
    0, 424.26 s, 5692.01 s (94.9'), 56921.0 s (15.8h) """

    egen_min = {"20deg": 12.6*u.GeV,
                "40deg": 26.4*u.GeV,
                "60deg": 105.3*u.GeV}
    """ Minimal acceptable True (i.e. generated) energies depending on the
    IRF zenith angle.
    Note: The masking later will remove the bin containing a given E value.
    If the E value is an edge, the subsequent bin is lost.
    The minimal and maximal energies need therefore to be slighlty
    before, resp. after the first, resp last edge of interest."""

    etrue_max = {"20deg": 17*u.TeV,
                 "40deg": 17*u.TeV,
                 "60deg": 17*u.TeV}
    """ Minimal acceptable reconstructed energies depending on the IRF zenith
    angle.
    Note: The masking later will remove the bin containing the E value.
    If the E value is an edge, the subsequent bin is lost.
    The minimal and maximal energies need therefore to be slighlty
    before, resp. after the first, resp. last edge of interest."""

    nbin_per_decade = 4
    """Bins per decade for the true energy axis"""

    # ##-----------------------------------------------------------------------
    def __init__(self,
                 filename=None,
                 irf=None,
                 subarray=None,
                 etrue=None,
                 kzen=None,
                 kaz=None,
                 kdt=None):
        """
        Initialise the object. Not in use.

        Parameters
        ----------
        name : String, optional
            The IRF file name on disk. The default is "none".
        irf : TYPE, optional
            DESCRIPTION. The default is None.
        array : String
            The array configuration. default is None
        etrue : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.


        """
        self.filename = filename
        self.irf = irf
        self.subarray = subarray
        self.etrue = etrue
        self.kzen = kzen
        self.kaz = kaz
        self.kdt = kdt

    # ##-----------------------------------------------------------------------
    @classmethod
    def from_observation(cls,
                         zenith=0*u.deg,
                         azimuth=None,
                         obstime=0*u.h,
                         subarray=None,
                         loc=None,
                         nsb=None,
                         irf_dir=None,
                         closest=False):
        """
        Get the IRF data from the characteristics of an observation.

        Notes
        -----
        * `MapAxis` accepts ony a normalised list of axis type as described
        `here
          <https://docs.gammapy.org/dev/irf/index.html#irf-axis-naming>`_.
        * The true Energy Axis can be defined on-the-fly, but it is not
          optimal for masking except if the number of bins is very large:

        ..  code-block:: python

            erec_axis  = MapAxis.from_energy_bounds(min(erec_min.values()),
                                                    max(erec_max.values()),
                                                    nbin = nbin_per_decade,
                                                    per_decade=True,
                                                    name="Rec. energy")

        Parameters
        ----------
        cls : IRF class instance
            Current object.
        zenith : Quantity angle, optional
            Observation zenith angle. The default is 0*u.deg.
        azimuth : Quantity angle, optional
            Observation azimuth angle. The default is 0*u.deg.
        obstime : Quantity time, optional
            Observation time. The default is 0*u.h.
        subarray : String, optional
            Subarray within the simulation. The default is None.
        loc : String, optional
            Site location ("North", "South"). The default is None.
        nsb : TYPE, optional
            Handle special NSB cases. The default is None.
        irf_dir : string, optional
            IRF folder. The default is None.
        closest : Boolean, optional
            Choose IRF closes to IRF range bounds. The default is False.

        Returns
        -------
        IRF ojbect
            Initialise an IRF object with the obtained values.

        """
        inst = cls()

        if loc is None:
            sys.exit("from_observation : location should be defined")

        # Find best keys (IRF interpolation)
        # Assumes irf_dir is a Path ending with the correct keyword
        prod = irf_dir.name[:5]

        kzen, kaz, kdt = inst.find_best_keys(zenith, azimuth, obstime,
                                             closest=closest,
                                             prod=prod)

        # Find base folder from the keys
        if prod == "prod3":
            folder = Path(irf_dir, subarray, loc, kzen)
            if subarray != "FullArray":
                kaz = kaz+"_"+subarray
            subfolder = loc+"_z"+kzen[:2]+"_"+kaz+"_"+kdt
            inst.filename = Path(folder, subfolder, "irf_file.fits.gz")

        elif prod == "prod5":
            folder = Path(irf_dir, subarray)
            if kaz == "average":
                filename = "Prod5-"+loc+"-"+kzen+"-AverageAz-"\
                         + subarray+"."+kdt+"-v0.1.fits.gz"
            else:
                sys.exit(f"Azimuth is {kaz:}, but only average is implemented")

            inst.filename = Path(folder, filename)

        else:
            sys.exit("Not implemented")

        if inst.filename.exists() is False:
            sys.exit(f" This file does not exist : {inst.filename:}")

        if gammapy.__version__ < "1.2":
            inst.irf = load_cta_irfs(inst.filename)
        else:
            inst.irf = load_irf_dict_from_file(inst.filename)

        if gammapy.__version__ < "1":
            eirf_min = min(inst.irf["aeff"].data.axes["energy_true"].edges)
        else:
            eirf_min = min(inst.irf["aeff"].axes["energy_true"].edges)

        inst.etrue = MapAxis.from_energy_bounds(eirf_min,
                                                inst.etrue_max[kzen],
                                                nbin=inst.nbin_per_decade,
                                                per_decade=True,
                                                name="energy_true")

        inst.subarray = subarray
        inst.kzen = kzen
        inst.kaz = kaz
        inst.kdt = kdt

        return inst

    # ##-----------------------------------------------------------------------
    def find_best_keys(self, zenith, azimuth, obstime,
                       closest=False, prod=None):
        """
        Find the best keys to later identify the IRF data file.

        Parameters
        ----------
        zenith : Quantity, angle
            zenith angle.
        azimuth : Quantity, angle
            Azimuth angle.
        obstime : Quantity, time
            Observation duration.
        closest : Boolean, optional
            If True, obtain zenith and observation time from the closest
            available sampling, instead of using the predefined validity value.
            The default is False.
        prod: String
            Monte Carlo production identifier. Default is None.

        Returns
        -------
        String, String, String
            The three strings defining the IRF file and/or folder.

        """
        # ##--------------------------
        def find_closest_key(mydict, x):
            dict_values = np.asarray([x.value for k, x in mydict.items()])
            closest_val = dict_values[(np.abs(dict_values - x.value)).argmin()]
            closest_key = [k for k, x in mydict.items()
                           if x.value == closest_val]
            return closest_key[0]
        # ##--------------------------

        # Zenith

        if closest:
            kzen = find_closest_key(self.zenith_list, zenith)
        else:
            found = False
            for key, val in self.zenith_valid.items():
                if val[0] <= zenith < val[1]:
                    kzen = key
                    found = True
                    continue
            if not found:
                sys.exit(f"irf.find_best_keys: zenith= {zenith} "
                         "=> range not found")

        # Azimuth - implement here N, S choice
        kaz = "average" if azimuth is None else azimuth

        # Observation time
        if closest:
            kdt = find_closest_key(self.dt_list[prod], obstime.to(u.s))
        else:
            found = False
            for key, val in self.dt_log_valid[prod].items():
                if obstime.to(u.s).value >= val[0] \
                   and obstime.to(u.s).value < val[1]:
                    kdt = key
                    found = True
                    continue
            if not found:
                sys.exit("irf.find_best_keys: obstime range not found")

        return kzen, kaz, kdt

    # ##-----------------------------------------------------------------------
    def print(self):
        """
        Print out some IRF class contents.

        Returns
        -------
        None.

        """
        print("IRF             : ", self.filename)
        print(" Etrue          : ", self.etrue)
        print("          edges : ", self.etrue.edges)
        print(" Sub-array      : ", self.subarray)
        print(" Zenith         : ", self.kzen)
        print(" Azimuth        : ", self.kaz)
        print(" Duration       : ", self.kdt)

# #############################################################################
# ## Utilities and check plots
# #############################################################################


# ##------------------------------------------------------------------------
def containment_plot(irf,
                     eunit="GeV", erec_min=10*u.GeV, erec_max=100*u.TeV,
                     subarray=None, tag=None, ax=None):
    """
    Display containment plots.

    Parameters
    ----------
    irf : IRF instance
        The current IRF instance.
    eunit : astropy.unit, optional
        Energy unit. The default is "GeV".
    erec_min : astropy.Quantity, optional
        Minimal energy. The default is 10*u.GeV.
    erec_max : astropy.Qauntity, optional
        Maximal energy. The default is 100*u.TeV.
    subarray : String, optional
        Subarray. The default is None.
    tag : String, optional
        Plot label. The default is None.
    ax : matplotlib.axes, optional
        Current axis. The default is None.

    Returns
    -------
    None.

    """
    if ax is None:
        ax = plt.subplots()[1]
    # irfname =  self.filename.parts[-2]

    # Add lower bin to exisiting array
    unit = "GeV"
    e_edges = np.append(np.array([10., 20.]),
                        mcf.erec_edges[subarray].to(unit).value)
    e2plot = MapAxis.from_edges(e_edges, unit=unit,
                                name="energy", interp="log")
    if gammapy.__version__ == "0.18.2":
        radii = irf.irf['psf'].containment_radius(energy=e2plot.edges,
                                                  theta=mcf.offset[subarray],
                                                  fraction=mcf.containment)[0]
    elif gammapy.__version__ == "1.2":
        radii = irf.irf['psf'].containment_radius(energy_true=e2plot.edges,
                                                  offset=mcf.offset[subarray],
                                                  fraction=mcf.containment)
    ax.plot(e2plot.edges.to(eunit).value, radii.value,
            marker="o", alpha=0.5, label=tag)
    ax.axvline(erec_min.to(eunit).value, ls=":")
    ax.axvline(erec_max.to(eunit).value, ls=":")
    ax.set_xlabel("Energy ("+eunit+")")
    ax.set_ylabel("Containment radius (°)")
    ax.set_xscale("log")
    ax.legend()


# ##-----------------------------------------------------------------------
def onoff_sketch_plot(irf, emin=30*u.GeV, subarray=None, tag=None,
                      nmaxcol=4,
                      alpha=0.2, fov=2.5*u.deg,
                      debug=False):
    """
    On-off geometry sketches versus reconstructed energy edges.

    Parameters
    ----------
    irf : IRF instance
        Current IRF.
    emin : astropy.Quantity, optional
        Minimal energy. The default is 30*u.GeV.
    subarray : String, optional
        Subarray. The default is None.
    tag : String, optional
        Legend text. The default is None.
    nmaxcol : integer, optional
        Maximum number of columns in the plot. The default is 4.
    alpha : float
        One over the number of on/off regions. The default is 0.2.
    fov: Astropy Quantity
        Field-of-view. The default is 2.5 degrees.
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    None.

    """
    # Retrieve radii at EDGES
    # Assumes that Erec = Etrue as the PSF is given versus Etrue
    if gammapy.__version__ == "1.2":
        radii = irf.irf['psf']\
             .containment_radius(energy_true=mcf.erec_edges[subarray],
                                 offset=mcf.offset[subarray],
                                 fraction=mcf.containment)
    else:
        radii = irf.irf['psf']\
            .containment_radius(energy=mcf.erec_edges[subarray],
                                theta=mcf.offset[subarray],
                                fraction=mcf.containment)[0]
    if debug:
        print(72*"=")
        print(radii.value)

    # ##-----------------------------
    def onoffsketch(radius, energy, ax=None):
        """
        Plot the on-off sketch n the field of view.

        Use the containment radius at a given energy (used for labelling only)
        """
        if ax is None:
            ax = plt.subplots()

        # On region
        x_on = mcf.offset[subarray].value
        y_on = 0

        # Default on region (useless)
        on_reg = plt.Circle((x_on, y_on), mcf.on_size[subarray].value,
                            fill=True, color="tab:blue", alpha=0.1)
        # 68% Aeff contained region
        on68 = plt.Circle((x_on, y_on), radius.value,
                          fill=True, color="green", alpha=0.5)
        # ax.text(x_on,y_on,s="on")

        # Build label
        txt = str(round(energy.value, 1)) + " " + str(emin.unit)
        # txt+= " - "+ str(round(100*mcf.containment,0)) + "%"
        # ax.set_title('Field of view -'+txt )
        ax.add_artist(on_reg)
        ax.add_artist(on68)
        ax.legend([on68], [txt])

        # Equal acceptance circle
        accept = plt.Circle((0, 0), mcf.offset[subarray].value,
                            fill=False, color="black", ls="--", alpha=1)
        ax.add_artist(accept)

        # Off regions
        for i in range(1, noff_reg):
            theta = i*dtheta
            x0 = x_on*np.cos(theta)
            y0 = x_on*np.sin(theta)
            off_reg = plt.Circle((x0, y0), radius.value,
                                 fill=True, color="red", alpha=0.5)
            ax.add_artist(off_reg)

        # Mark center of FOV
        ax.axvline(x=0, ls=":")
        ax.axhline(y=0, ls=":")

        # Set field of view, and aspect ratio
        view = fov.value/2
        ax.set_xlim([-view, view])
        ax.set_ylim([-view, view])

        ax.set_aspect(1)  # Nice but conflicts with grid spacing

    # ##-----------------------------

    noff_reg = int(1/alpha) + 1
    dtheta = 360*u.degree/noff_reg

    xsize = 3
    ysize = 3
    erec = mcf.erec_edges[subarray]
    mask = mcf.erec_edges[subarray] >= emin
    nplots = len(erec[mask]) - 1
    ifirst = np.where(mask)[0][0]
    ncols = min(nmaxcol, nplots)  # If nplots < nmaxcol, take nplots
    nrows = int(nplots / ncols) + 1 * (nplots % ncols != 0)

    f_adj = 0.92  # adjusts the size when plot scale is forced to square
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows,
                           figsize=(f_adj*xsize*ncols, ysize*nrows),
                           sharex=True, sharey=True)

    iplot = 0

    for jrow, icol in itertools.product(range(nrows), range(ncols)):

        ax0 = ax[jrow][icol] if nrows > 1 else ax[icol]

        if iplot < nplots:
            radius = radii[iplot+ifirst]
            energy = mcf.erec_edges[subarray][iplot+ifirst].to(emin.unit)
            onoffsketch(radius, energy, ax=ax0)
        else:
            ax0.axis('off')
            continue  # Next plot

        # Compactify
        if jrow+1 != nrows:
            ax0.set_xlabel(None)
        if icol != 0:
            ax0.set_ylabel(None)
        ax0.tick_params(which='major', length=10, width=2, direction='in')
        ax0.tick_params(which='minor', length=5, width=2, direction='in')

        iplot += 1

    # Figure title
    fig.suptitle(irf.filename.parts[-2]+tag, fontsize=18, y=1.02)

    # known feature/bug : does not take supttile into account !
    # fig.suptitle("Title centered above all subplots", fontsize=14)
    fig.tight_layout(h_pad=0, w_pad=0)

    # Adjust AFTER tight_layout
    plt.subplots_adjust(hspace=0, wspace=0, top=0.95)


# ##-----------------------------------------------------------------------
def aeff_plot(irf, axs=None,
              ethreshold=100*u.GeV,
              min_fraction=0.05,
              unit="GeV",
              tag="",
              yscale="log"):
    """
    Display effective areas.

    Parameters
    ----------
    irf : IRF instance
        Current IRF.
    axs: Matplotlib axes
        Set of 3 axes on a line
    ethreshold : astropy.Quantity, optional
        Energy threshold to be displayed. The default is 100*u.GeV.
    min_fraction : float, optional
        A fraction of the max. effective area to be displayed.
        The default is 0.05.
    unit : astropy.unit, optional
        The plot energy unit. The default is "GeV".
    tag : String, optional
        A text for the label. The default is "".

    Returns
    -------
    None.

    """
    if axs is None:
        axs = plt.subplots(nrows=1, ncols=2)

    # Effective area
    effarea = irf.irf["aeff"].data

    if gammapy.__version__ < "1.2":
        e_edges = effarea.axes[0].center
    else:
        e_edges = irf.irf["aeff"].axes[0].center

    # Loop over offsets in the field of view
    for joff, off in enumerate([0*u.deg, 0.5*u.deg, 1*u.deg]):

        ax = axs[joff]

        if gammapy.__version__ < "1.2":
            effoff = irf.irf["aeff"].to_effective_area_table(off)
            effmax = effoff.max_area
            emin = effoff.find_energy(effmax*min_fraction)
            print(f" {emin[0].to(unit).value:5.1f} ", end="")

        else:
            effoff = irf.irf["aeff"].evaluate(offset=off)
            effmax = max(effoff)[0]
            emin = e_edges[np.where(effoff >= effmax*min_fraction)[0][0]]
            print(f" {emin.to(unit).value:5.1f} ", end="")

        with quantity_support():

            label = tag + "\n"+str(off.value)+"° "

            # effarea.evaluate(energy_true = e_edges,offset =off))
            if gammapy.__version__ < "1.2":
                p = ax.plot(e_edges,
                            effoff.data.evaluate(energy_true=e_edges)/1e6,
                            label=label)
                ax.axhline(y=min_fraction*effmax/1e6, ls="--",
                           color=p[0].get_color(),
                           label=str(100*min_fraction)+"%")
            else:
                p = ax.plot(e_edges,
                            irf.irf["aeff"].evaluate(offset=0.5*u.deg,
                                                     energy_true=e_edges),
                            label=label)
                ax.axhline(y=min_fraction*effmax, ls="--",
                           color=p[0].get_color(),
                           label=str(100*min_fraction)+"%")

            ax.axvline(ethreshold, color="black", ls="--", label="Threshold")
            ax.axvline(emin,
                       ls="--",
                       color=p[0].get_color())
            ax.set_xscale("log")
            ax.set_yscale(yscale)
            ax.legend()

            if joff > 0:
                ax.set_ylabel(None)

            ax.grid("both", which="major", alpha=0.5)
            ax.grid("both", which="minor", alpha=0.3)
            single_legend(fig)

    print()
    plt.tight_layout(h_pad=0, w_pad=0)


###############################################################################
if __name__ == "__main__":

    # Code example to use the IRF class

    # Bigger texts and labels
    sns.set_context("notebook")  # poster, talk, notebook, paper

    print(" Running with gammapy ", gammapy.__version__)

    # Suppress invalid unit warning when reading IRF files
    import logging
    logging.basicConfig()
    log = logging.getLogger("gammapy.irf")
    log.setLevel(logging.ERROR)

    infolder = "D:/CTAO/SoHAPPy/input/"
    irf_dir = "irf/Full/prod3-v2"
    irf_dir = "irf/Full/prod5-v0.1"

    irf_dir = Path(infolder, irf_dir)

    array = {"North": "FullArray", "South": "FullArray"}
    array = {"North": "4LSTs09MSTs", "South": "14MSTs37SSTs"}
    # array = {"North": "4LSTs09MSTs", "South": "4LSTs14MSTs40SSTs"}

    prod = irf_dir.name[:5]
    print(" Processing ", prod, " files")

    # Choose plot to show in the following
    show = ["containment", "onoff", "effarea", "generic", 'edisp']
    # show = ['containment']

    print(IRF.dt_log_valid)

    # ##-------------------------
    # ## E dispersion
    # ##-------------------------
    if "edisp" in show:

        npoints = 25
        ratio = np.linspace(0.25, 1.75, npoints)

        # This list allows investigating the E-dispersion matrix
        etrue = np.array([30, 40, 60, 110, 200, 350])*u.GeV

        # These lists are used to sample the IRF file parameters
        zenlist = [21*u.deg, 41*u.deg, 61*u.deg]
        dtlist = [100*u.s, 1*u.h, 10*u.h, 100*u.h]

        for loc in ["North", "South"]:

            offset = mcf.offset[array[loc]]
            print(loc, " -> Offset = ", offset)

            for dt in dtlist:
                fig, ax = plt.subplots(nrows=1, ncols=3,
                                       figsize=(17, 6), sharey=True)
                iplot = 0

                # Display the 3 zenith plots
                for ax0, zenith in zip(ax, zenlist):

                    irf = IRF.from_observation(loc=loc,
                                               subarray=array[loc],
                                               zenith=zenith,
                                               azimuth=None,
                                               obstime=dt,
                                               irf_dir=irf_dir)
                    print(" ", irf.filename)

                    # Display E dispersion for true energies
                    irf.irf["edisp"].plot_migration(ax0,
                                                    offset=offset,
                                                    #  migra= ratio,
                                                    energy_true=etrue,
                                                    alpha=0.5,
                                                    marker=".")

                    # Superimpose the dispersion for the threshold energy
                    ethreshold = [mcf.erec_min[array[loc]][irf.kzen]
                                  + mcf.safe_margin]

                    irf.irf["edisp"].plot_migration(ax0,
                                                    offset=offset,
                                                    # migra=ratio,
                                                    energy_true=ethreshold,
                                                    alpha=1.0,
                                                    marker="o")

                    # Change color of threshold and simplify labels
                    lines, labels = ax0.get_legend_handles_labels()
                    iline = 0
                    elist = np.append(etrue, ethreshold)
                    for line, E in zip(lines, elist):
                        if iline < len(lines) - 1:
                            txt = str(E.value) + " " + str(E.unit)
                        else:
                            txt = "Threshold"
                            line.set_color("black")

                        line.set_label(txt)
                        # print("***",iline," ->", txt)
                        iline += 1

                    # Display acceptable +/- 25% area
                    ax0.axvspan(xmin=0.75, xmax=1.25, alpha=0.2, color="grey",
                                label=r"$\pm 25\%$")
                    ax0.axvline(x=1.0, alpha=0.5, color="grey")

                    # Titles
                    ax0.set_title(loc + " - " + irf.kzen+" -  " + irf.kdt
                                  + " - " + str(offset.value) + "°")

                    # Decoration
                    # ax0.set_yscale("log")
                    if iplot < 2:
                        ax0.get_legend().remove()
                    if iplot != 0:
                        ax0.set_ylabel(None)
                    iplot += 1

                ax[2].legend(bbox_to_anchor=[1, 0.5], fontsize=12)

                # Rearrange the plots
                plt.tight_layout(w_pad=0)

    # ##-------------------------
    # ## Generic
    # ##-------------------------
    if "generic" in show:
        for loc in ["North", "South"]:
            irf = IRF.from_observation(loc=loc,
                                       subarray=array[loc],
                                       zenith=20*u.deg,
                                       azimuth=None,
                                       obstime=100*u.s,
                                       irf_dir=irf_dir)
            irf.print()
            irf.irf["psf"].peek()
            irf.irf["edisp"].peek()

        print(irf.dt_log_valid)

    # ##-------------------------
    # ## Show containment radii
    # ##-------------------------
    if "containment" in show:

        for loc in ["North", "South"]:
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 5),
                                   sharex=True, sharey=True)
            iplot = 0
            for ztag, axi in zip(IRF.zenith_list, ax):
                for dt in IRF.dt_list[prod].values():
                    irf = IRF.from_observation(loc=loc,
                                               subarray=array[loc],
                                               zenith=IRF.zenith_list[ztag],
                                               azimuth=None,
                                               obstime=dt,
                                               irf_dir=irf_dir)
                    print("---->", irf.filename, dt)
                    containment_plot(irf, subarray=array[loc],
                                     erec_min=mcf.erec_min[array[loc]][ztag],
                                     erec_max=mcf.erec_max[ztag],
                                     ax=axi,
                                     tag=ztag[:2]+"° - "+t_str(dt))

                    axi.grid("both", which="major", alpha=0.8)
                    axi.grid("both", which="minor", alpha=0.5)
                    axi.set_ylim(ymax=0.5)

                    if iplot > 0:
                        axi.set_ylabel(None)
                    iplot += 1

            fig.suptitle(array[loc] + " " + loc, fontsize=18, y=0.95)
            plt.tight_layout(w_pad=0)
            # plt.subplots_adjust(top=0.95)

    # ##-------------------------
    # ## Show on-off sketch / containment radii
    # ##-------------------------
    nmaxcol = 6
    if "onoff" in show:
        for loc in ["North", "South"]:
            for ztag in IRF.zenith_list.keys():
                for dt in IRF.dt_list[prod].values():

                    irf = IRF.from_observation(loc=loc,
                                               subarray=array[loc],
                                               zenith=IRF.zenith_list[ztag],
                                               azimuth=None,
                                               obstime=dt,
                                               irf_dir=irf_dir)
                    print(f"---->{irf.filename}")
                    emin = mcf.erec_min[array[loc]][ztag]
                    onoff_sketch_plot(irf, emin=emin, nmaxcol=nmaxcol,
                                      subarray=array[loc],
                                      tag=" - "+ztag[:2]+"° - "+t_str(dt))

    # ##-------------------------
    # ## Effectiva area plots
    # ##-------------------------
    if "effarea" in show:

        unit = "GeV"
        min_fraction = 0.05  # 0.05, 0.1

        print(f" *** Threshold for {100*min_fraction:2.0f}% "
              f"eff.max ({unit:3s})  *** ")

        for loc in ["North", "South"]:

            fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 12),
                                   sharex=True, sharey=True)
            fig.suptitle(array[loc] + "-" + loc, fontsize=24)

            print(f"{loc:5s} {array[loc]:15s} {'0°':7s} {'0.5°':7s} {'1°':7s}")

            for z, ax3 in zip([20, 40, 57], ax):

                for dt in [100*u.s, 0.5*u.h, 5*u.h, 50*u.h]:

                    irf = IRF.from_observation(loc=loc,
                                               subarray=array[loc],
                                               zenith=z*u.deg,
                                               azimuth=None,
                                               obstime=dt,
                                               irf_dir=irf_dir)

                    print(f"{str(z)+'° '+str(dt.value)+' '+str(dt.unit):20s}:",
                          end="")

                    aeff_plot(irf, axs=ax3, unit=unit,
                              min_fraction=min_fraction,
                              ethreshold=mcf.erec_min[irf.subarray][irf.kzen],
                              tag=irf.kzen + "\n" + irf.kdt,
                              yscale="log")

            plt.subplots_adjust(top=0.95)
