# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:12:17 2025.

@author: Stolar
"""
import os
import sys
import numpy as np
import pickle
from pathlib import Path

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from astropy.table import Table, QTable
from astropy.visualization import quantity_support

from niceplot import MyLabel
from configuration import Configuration


class Fluxes():
    """Blabla."""

    # ##-----------------------------------------------------------------------
    def __init__(self):
        """Initialize."""
        self.Eval = None
        self.tval = None
        self.fluxes = None
        self.Fpeaks = None

    # ##-----------------------------------------------------------------------
    @classmethod
    def from_fits(cls, filelist):
        """
        Get the times.

        Parameters
        ----------
        filelist : TYPE
            DESCRIPTION.

        Returns
        -------
        tval : TYPE
            DESCRIPTION.
        Eval : TYPE
            DESCRIPTION.
        fluxval : TYPE
            DESCRIPTION.
        """
        cls = Fluxes()

        ngrb = len(filelist)
        first = True

        fluxes = []
        fpeaks = []

        for item, fname in enumerate(filelist):

            # Read file, header and keys
            hdul = fits.open(fname)
            hdr = hdul[0].header
            keys_0 = list(hdul[0].header.keys())

            # It can be long , give news from time to time
            if first is True:
                print("Processing items :", end=" ")
            if np.mod(item, int(ngrb/50)) == 0 or ngrb < 10:
                print(item, end=' ')

            # Get the time and energy bins from the first file only
            # It assumes that all data files have the same binning
            if first:

                tval = QTable.read(hdul["TIMES (AFTERGLOW)"])
                cls.tval = np.array(tval[tval.colnames[0]].value)\
                    * tval[0][0].unit

                tab_key = "Energies (afterglow)"
                col_key = Table.read(hdul[tab_key]).colnames[0]
                Eval = Table.read(hdul[tab_key])[col_key].quantity
                cls.Eval = np.array(Eval)*Eval[0].unit

                first = False

            # Peak flux - used for seclecting the data
            if "PHFLUX" in keys_0:
                Fpeak = hdr['PHFLUX']*u.Unit("cm-2.s-1")
            elif "PHFX" in keys_0:
                Fpeak = hdr['PHFX']*u.Unit("cm-2.s-1")
            fpeaks.append(Fpeak.value)

            # Afterglow flux
            flux = QTable.read(hdul["SPECTRA (AFTERGLOW)"])
            flux_unit = u.Unit(flux.meta["UNITS"])
            if str(flux_unit).find("ph") > -1:
                flux_unit = flux_unit/u.Unit("ph")  # Removes ph

            # Store the flux. Note the transposition
            itmax = len(tval) - 1
            jEmax = len(Eval) - 1
            fluxval = np.zeros((itmax+1, jEmax+1))

            for i in range(0, itmax+1):
                for j in range(0, jEmax+1):
                    fluxval[i][j] = flux[j][i]
            fluxes.append(fluxval)

        cls.fluxes = np.array(fluxes)*flux_unit
        cls.Fpeaks = np.array(fpeaks)*Fpeak.unit

        return cls

    # ##-----------------------------------------------------------------------
    def FluxPeakTime(self, Eref, debug=True):
        """
        Get the maximal flux in time at the reference energy.

        Parameters
        ----------
        Eref : Astropy Quantity, optional
            Reference energy. The default is 100*u.GeV.

        Returns
        -------
        idmax : integer
            Time index of the maximal flux.
        idxref : integer
            Energy index of teh reference energy.

        """
        # Get energy closest to Eref
        idxref = np.abs(self.Eval.to(Eref.unit).value - Eref.value).argmin()
        if debug:
            print(32*"-")
            print(f" Closest energy ({Eref:}) : idx={idxref:}, "
                  f"E={self.Eval[idxref]:5.2f}")
        # Loop over the stored flux

        flx_max = []
        t_flx_max = []
        for item, flx in enumerate(self.fluxes):
            # Obtain time spectrum for that energy
            fluxref = flx[:, idxref]

            # Get the index of maximal flux
            idmax = fluxref.argmax()  # Get time index of max value

            if debug:
                print(f" Maximal flux {item:5d}: {fluxref[idmax]:5.2e} "
                      f"@index {idmax:2d} tmax={self.tval[idmax]:5.2f}")

            flx_max.append(fluxref[idmax].value)
            t_flx_max.append(self.tval[idmax].value)

        self.flux_max = np.array(flx_max)*fluxref.unit
        self.t_flux_max = np.array(t_flx_max)*self.tval.unit

    # ##-----------------------------------------------------------------------
    def dump_to_file(self, filename):
        """Dimp class instance to disk."""
        with open(filename, "wb") as outfile:
            pickle.dump(self, outfile)
            print(" Class instance was dumped to ", filename)

    # ##-----------------------------------------------------------------------
    @classmethod
    def read_from_disk(cls, filename):
        """Read class instance from disk."""
        with open(Path(filename), "rb") as f:
            cls = pickle.load(f)
        return cls


###############################################################################
if __name__ == "__main__":

    sys.path.append("../../")

    outfilename = "explore_simulation_from_fits.bin"
    save = False

    # Bigger texts and labels
    import seaborn as sns
    sns.set_context("poster")  # poster, talk, notebook, paper

    if save is True:
        # Build a default Configuration to use the sourc_ids function
        sys.argv = ["", "-c", "../../data/config_ref.yaml", "-f", "1"]
        cf = Configuration.command_line()

        cf.ifirst = 1
        cf.nsrc = 10

        filelist = cf.source_ids(os.environ["HAPPY_IN"])
        flxs = Fluxes.from_fits(filelist)

        print(" Constant number of bins: ")
        print("   Energies : ", len(flxs.Eval))
        print("   Times    : ", len(flxs.tval))
        flxs.dump_to_file(outfilename)

    else:
        flxs = Fluxes.read_from_disk(outfilename)

    # Plot Fpeak distribution
    with quantity_support():
        fig, ax = plt.subplots()
        ax.hist(np.log10(flxs.Fpeaks.value))
        ax.set_xlabel("log 10 Fpeak (" + str(flxs.Fpeaks[0].unit) + ")")

    # Get Maximal flux time and plot
    bins = np.logspace(0, 5, 50)
    energies = [50*u.GeV, 100*u.GeV, 500*u.GeV]
    fig, axes = plt.subplots(nrows=1, ncols=len(energies),
                             figsize=(len(energies)*10, 10), sharey=True)

    axes[0].set_ylabel("Counts")
    for Eref, ax in zip(energies, axes):
        flxs.FluxPeakTime(Eref=Eref)
        with quantity_support():
            _, bins, _ = ax.hist(flxs.t_flux_max, bins=bins,
                                 label=MyLabel(flxs.t_flux_max,
                                               label=str(Eref)))

            mask = np.where(flxs.Fpeaks.value >= 1)
            ax.hist(flxs.t_flux_max[mask], bins=bins, alpha=1.0,
                    label=MyLabel(flxs.t_flux_max[mask],
                                  label=r"$F_{peak} > 1$"))

            ax.set_xlabel(r"$T_{ \phi_{max} }$ ("
                          + ax.get_xlabel()+")")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend()
            ax.grid(axis="both",  ls=":")
    plt.tight_layout()

    print("...completed")
