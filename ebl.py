# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 15:24:21 2021

This module contains the functions to decode and plot the EBL absorption model
beyond waht exist in `Gammapy`, as well as a main function illustrating the
use of these functions.

@author: Stolar
"""

import sys
import os

from pathlib import Path
from itertools import cycle

import numpy as np
import astropy.units as u

import matplotlib.pyplot as plt

from scipy.interpolate import interp2d
from matplotlib import cm
from astropy.visualization import quantity_support

from gammapy.modeling.models import EBLAbsorptionNormSpectralModel

from niceplot import single_legend

# Bigger texts and labels
import seaborn as sns

__all__ = ["EBL_from_file", "EBL_plot", "plot_dominguez_2011",
           "compare_EBL_models"]


# -----------------------------------------------------------------------------
def EBL_from_file(filename, debug=False):
    """
    Read a specific EBL model from file.

    Note that the formats are not
    normalised. They are different for each file, depending on the provider.

    Parameters
    ----------
    file : String
        Input file name.
    debug : Boolean, optional
        If True display information. The default is False.

    Returns
    -------
    attenuation : interpd2d object
        Interpolated attenuation.

    """
    if filename.find("gilmore") != -1:
        data = np.loadtxt(filename)
        z_dat = data[0, 1:]
        E_dat = (data[1:, 0]/1000)*u.GeV
        EBLabs = np.array(data[1:, 1:]).T
        EBLabs = np.exp(-EBLabs)
        if debug:
            print(" Gilmore from Lara")

    elif filename.find("dominguez") != -1:
        data = np.loadtxt(filename)
        E_dat = data[0, 1:]*u.GeV
        z_dat = data[1:, 0]
        EBLabs = data[1:, 1:]
        if debug:
            print(" Dominguez from Renaud")
    else:
        sys.exit("ebl.py: cannot process this EBL data fle")

    return interp2d(E_dat, z_dat, EBLabs, kind="cubic")


# -----------------------------------------------------------------------------
def EBL_plot(eblmodels,
             redshifts=0.1, energies=100*u.GeV, tags="",
             axa=None, xsize=8, ysize=8,
             color=None, ls=None,
             **kwargs):
    """
    Plot various EBL model absortpion curve.

    Parameters
    ----------
    eblmodels : list of interp2d objects
        List of EBL models.
    redshifts : list of float or float, optional
        Redshift or list of redshifts. The default is 0.1.
    energies : List of Quantity or Quantity, optional
        Energy or list of energies with unit. The default is 100*u.GeV.
    tags : List of strings or string, optional
        Model tags. The default is "".
    ax : matplotlib.axes, optional
        Current axis. The default is None.
    xsize : float, optional
        Figure horizontal size. The default is 8.
    ysize : float, optional
        Figre vertical size. The default is 8.
    color : Boolean, optional
        If True, colors from a colormap. The default is None.
    ls : string, optional
        Matplotlib line style. The default is None.
    **kwargs : pointer list
        Extra arguments.

    Returns
    -------
    ax : Matplolib.axes
        Current plot axis.

    """
    # If not a list but a single element
    if not isinstance(eblmodels, list):
        eblmodels = [eblmodels]
    if not isinstance(tags, list):
        tags = [tags]
    if not isinstance(redshifts, np.ndarray):
        redshifts = [redshifts]
    if not isinstance(energies, np.ndarray):
        energies = [energies]

    # A line style per attenuation model
    linecycler = cycle(["-", "--", "-.", ":"])

    set_color = False
    set_line = False

    if axa is None:
        _, axa = plt.subplots(figsize=(xsize, ysize))

    if color is None:
        set_color = True

    if ls is None:
        set_line = True

    with quantity_support():

        for att, tag in zip(eblmodels, tags):
            if set_line:
                ls = next(linecycler)

            for i, z in enumerate(redshifts):

                if set_color:
                    color = cm.cool(i/len(redshifts))  # One color per z
                else:
                    color = "grey"

                if isinstance(att, EBLAbsorptionNormSpectralModel):
                    attenuation = att.evaluate(energies, redshift=z,
                                               alpha_norm=1)
                else:
                    attenuation = att(energies, z)
                axa.plot(energies, attenuation,
                         label=tag+" z="+str(round(z, 2)),
                         ls=ls, color=color, **kwargs)

        axa.set_xscale("log")
        axa.set_ylim(ymin=1e-2, ymax=2)
        single_legend(axa)

    return axa


# -----------------------------------------------------------------------------
def plot_dominguez_2011(redshifts, energies, alpha=0.8, yscale="linear"):
    """
    Plot Dominqguez absorption curve. Probably redundant with EBL_plot.

    To be optimised

    Parameters
    ----------
    redshifts : list of float
        List of redshifts.
    energies : list of Quantity or Qauntity
        List of energies.
    alpha : float, optional
        Matplotlib color transparency. The default is 0.8.

    Returns
    -------
    None.

    """
    _, axa = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

    for i, z in enumerate(redshifts):
        do_z = dominguez.evaluate(energies, redshift=z, alpha_norm=1)
        do_z = np.array([(x == max_att if x < max_att else x) for x in do_z])
        color = cm.rainbow(i/len(zlist))

        axa.plot(Elist, do_z, color=color, alpha=alpha, label=str(round(z, 1)))

    axa.set_yscale(yscale)
    axa.set_xscale("log")
    axa.set_ylabel("Absorption")
    axa.set_xlabel("Energy (GeV)")
    axa.set_title("Dominguez 2011", size=18)
    axa.legend()


# -----------------------------------------------------------------------------
def compare_EBL_models(redshifts, energies, ratio=False, alpha=0.5):
    """
    Compare all availbale EBL models.

    Parameters
    ----------
    redshifts : list of float
        List of redshifts.
    energies : list of Quantity or Qauntity
        List of energies.
    ratio : Boolean, optional
        If True show ratii. The default is False.
    alpha : float, optional
        Matplotlib color transparency. The default is 0.5.

    Returns
    -------
    None.

    """
    if ratio:
        _, ((ax11, ax12, ax13),
            (ax21, ax22, ax23)) = plt.subplots(nrows=2, ncols=3,
                                               figsize=(15, 6),
                                               sharey=False, sharex=True,
                                               gridspec_kw={'height_ratios':
                                                            [3, 2]})
    else:
        _, (ax11, ax12, ax13) = plt.subplots(nrows=1, ncols=3,
                                             figsize=(15, 5),
                                             sharey=False, sharex=True)

    for i, z in enumerate(redshifts):
        # color = cm.cool(i/len(zlist))
        # color = cm.viridis(0.7-0.7*i/len(zlist))
        color = cm.rainbow(i/len(zlist))

        # Evaluate model at given z
        do_z = dominguez.evaluate(energies, redshift=z, alpha_norm=1)
        do_z = np.array([(x == max_att if x < max_att else x) for x in do_z])

        fr_z = franceschini.evaluate(Elist, redshift=z, alpha_norm=1)
        fr_z = np.array([x == max_att if x < max_att else x for x in fr_z])

        fi_z = finke.evaluate(Elist, redshift=z, alpha_norm=1)
        fi_z = np.array([x == max_att if x < max_att else x for x in fi_z])

        gi_z = gilmore_ln(Elist, z)
        gi_z = np.array([x == max_att if x < max_att else x for x in gi_z])

        fr_ratio = (fr_z - do_z) / do_z
        fi_ratio = (fi_z - do_z) / do_z
        gi_ratio = (gi_z - do_z) / do_z

        # First row, plot attentaion
        ax11.plot(Elist, do_z, color=color, alpha=alpha, label="Dominguez")
        ax11.plot(Elist, fr_z, color=color, ls="--", label=str(round(z, 1)))
        if not ratio:
            ax11.set_xlabel("Energy (GeV)")

        ax12.plot(Elist, do_z, color=color, alpha=alpha, label="Dominguez")
        ax12.plot(Elist, fr_z, color=color, ls="--", label=str(round(z, 1)))
        if not ratio:
            ax12.set_xlabel("Energy (GeV)")

        ax13.plot(Elist, do_z, color=color, alpha=alpha, label="Dominguez")
        ax13.plot(Elist, gi_z, color=color, ls="--", label=str(round(z, 1)))
        if not ratio:
            ax13.set_xlabel("Energy (GeV)")

        if ratio:
            ax21.plot(Elist, fr_ratio, color=color, alpha=0.7)
            ax22.plot(Elist, fi_ratio, color=color, alpha=0.7)
            ax23.plot(Elist, gi_ratio, color=color, alpha=0.7)

    # first row decoration
    ymin = 1e-2
    ymax = 1.05

    yscale = "linear"
    ax11.set_ylabel("Absorption")
    ax11.set_title("Franceschini 2008", size=18)

    ax12.set_ylabel(None)
    ax12.set_title("Finke 2010", size=18)
    ax12.grid("both", which="minor", ls=":", alpha=0.5)

    ax13.set_ylabel(None)
    ax13.set_title("Gilmore 2012", size=18)
    ax13.grid("both", which="minor", ls=":", alpha=0.5)
    single_legend(ax13, bbox_to_anchor=[1, 1.0])

    ax12.tick_params(left=False, labelleft=False, bottom=False)
    ax13.tick_params(left=False, labelleft=False, bottom=False)

    for ax in [ax11, ax12, ax13]:
        ax.set_ylim(ymin=ymin)
        ax.set_ylim(ymax=ymax)
        ax.set_yscale(yscale)
        ax.set_xscale("log")
        ax.grid("both", which="minor", ls=":", alpha=0.5)
        ax.grid("both", which="major", ls="-", alpha=0.5)

    if ratio:
        # Second row decoration
        yscale = "linear"
        ax21.set_ylabel(r"$Relative \ \Delta$")
        ax22.set_ylabel(None)
        ax23.set_ylabel(None)

        ax22.tick_params(left=False, labelleft=False)
        ax23.tick_params(left=False, labelleft=False)
        for ax in [ax21, ax22, ax23]:
            ax.set_xlabel("Energy (GeV)")
            ax.set_ylim(ymin=-1.5)
            ax.set_ylim(ymax=1.5)
            ax.set_yscale(yscale)
            ax.set_xscale("log")
            ax.grid("both", which="minor", ls=":", alpha=0.5)
            ax.grid("both", which="major", ls="-", alpha=0.5)
            ax.set_yticks(np.arange(-1.5, 1.5, 0.5))

    plt.tight_layout(h_pad=0, w_pad=0)


# -----------------------------------------------------------------------------
if __name__ == "__main__":

    sns.set_context("poster")  # talk, notebook, paper
    # Compare absorption models versus energy for various redshifts.
    # This is quite under development and can certainly be optimised

    # Read SoHAPPY absorptions - These are interp2d function!
    # Do not pass Quantity as argument
    gilmore_ln = EBL_from_file("data/ebl/others/"+"ebl_gilmore12.dat",
                               debug=True)
    file_rb = "data/ebl/others/EBL_abs_RBelmont-dominguez-20170425.dat"
    dominguez_rb = EBL_from_file(file_rb, debug=True)

    # Read Gammapy absorptions
    os.environ['GAMMAPY_DATA'] = str(Path(Path(__file__).absolute().parent,
                                          "data"))
    dominguez = EBLAbsorptionNormSpectralModel.read_builtin("dominguez")
    franceschini = EBLAbsorptionNormSpectralModel.read_builtin("franceschini")
    finke = EBLAbsorptionNormSpectralModel.read_builtin("finke")

    # Define redshift, energy space
    zmin, zmax, nzbin = 1, 5, 4
    emin, emax, nebin = 10*u.GeV, 10*u.TeV, 100

    zlist = np.append([0], np.logspace(np.log10(zmin), np.log10(zmax), nzbin))
    Elist = np.logspace(np.log10(emin.value),
                        np.log10(emax.to(emin.unit).value),
                        nebin)*emin.unit
    max_att = 1e-3  # Limit attenutaion value to avoid insane small numbers

    # Plot Dominguez, the reference model in CTA
    plot_dominguez_2011(zlist, Elist)

    # Display models
    compare_EBL_models(zlist, Elist, ratio=False)

    # Plot Model ratii compared to Dominquez
    compare_EBL_models(zlist, Elist, ratio=True)

    # Compare several models
    taglist = ["gilmore", "dominguez"]
    ax = EBL_plot([gilmore_ln, dominguez],
                  redshifts=zlist, energies=Elist, tags=taglist)

    ax.set_yscale("log")

    print("All done")
