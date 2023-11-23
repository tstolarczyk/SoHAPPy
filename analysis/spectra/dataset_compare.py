# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 12:10:37 2021

This module has functions to compare two dataset (It can be used to check the
dataset stacking fucntions).

This module has to be re-arranged and is not functionnal as it is!!!

@author: Stolar
"""
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import astropy.units as u

from gammapy.utils.random import get_random_state
import yaml
from yaml.loader import SafeLoader

from dataset_counts import excess_counts
from dataset_flux import extract_spectrum
from dataset_tools import createonoff_from_simulation, check_datasets, sigmax , compactify

from niceplot import stamp

__all__ = ['comparison_plot']
###----------------------------------------------------------------------------
def comparison_plot(GRB_id,site="South",
                    folder = "./",
                    param_file="parameter.yaml",
                    merging=True, counts=True, flux=True,
                    fit_tag = None,
                    e_min = None, e_max = None, e_fit_max = None,
                    arrays = ("omega","alpha"), debug=False):
    """
    Plot counts and/or flux for several configurations (arrays).

    Parameters
    ----------
    GRB_id : String
        GRB identifier (xyz for Eventxyz.fits).
    site : String, optional
        "South", "North". The default is "South".
    folder: String
        Source file folder; default is local folder.
    param_file: String
        Name of the file where the plot and fit parameters are read.
    merging : Boolean, optional
        Indicated if the dataset has been stacked. The default is True.
    counts : Boolean, optional
        If True, compare counts. The default is True.
    flux : Boolean, optional
        If True, compare flux. The default is True.
    fit_tag: String
        If None, will be read from the paramter file.
    e_min: astropy.Quantity
        Minimum energy in the count and flux plots.
    e_max: astropy.Quantity
        maximum energy in the count plots.
    e_fit_max: astropy.Quantity
        maximum energy in the flux plots.
    arrays : List of String, optional
        The arrays or configuration names to be compared.
        The default is ["omega","alpha"].
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    None.

    """

    #plt.style.use('seaborn-talk') # Make the labels readable
    plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!

    ###----------------------
    ### Get data
    ###----------------------
    base = Path(folder)

    # Read parameters
    par = (yaml.load(open(param_file), Loader=SafeLoader))
    GRB = GRB_id+site[0] # The key to find back the parameters

    if fit_tag is None:
        fit_tag = par[GRB]["fit_tag"]
    if e_min is None:
        e_min  = 20*u.GeV
    if e_max is None:
        e_max  = 10*u.TeV
    if e_fit_max is None:
        e_fit_max = u.Quantity(par[GRB]["e_fitmax"])
    # Open Omega file
    # filetag_1 = "GRB "+GRB_id+" - South ($\\"+array +"$)"
    array   = arrays[0]
    file    = Path(base,array,"Event"+GRB_id+"-"+site+"_sim.bin")
    infile  = open(file,"rb")
    mc_1    = pickle.load(infile)
    infile.close()
    dset_init_1 = createonoff_from_simulation(mc_1,
                                              random_state=get_random_state(2021))

    # Open Alpha file
    # filetag_2 = "GRB "+GRB_id+" - South ($\\"+array +"$)"
    array=arrays[1]
    file    = Path(base,array,"Event"+GRB_id+"-"+site+"_sim.bin")
    infile  = open(file,"rb")
    mc_2    = pickle.load(infile)
    infile.close()
    dset_init_2 = createonoff_from_simulation(mc_2,
                                              random_state=get_random_state(2021))

    ###----------------------
    ### Merge slices if requested
    ###----------------------
    if merging: # Stacking is also possible, Cf. dataset_tools
        dset_list_1 = compactify(dset_init_1,
                                 dtmin=u.Quantity(par[GRB]["dtmin"]))
        dset_list_2 = compactify(dset_init_2,
                                 dtmin=u.Quantity(par[GRB]["dtmin"]))
    else:
        dset_list_1 = dset_init_1.copy()
        dset_list_2 = dset_init_2.copy()

    ## Check datasets
    if debug:
        check_datasets(dset_init_1,deeper=True,masked=False)
        check_datasets(dset_list_1,deeper=True,masked=False)
        check_datasets(dset_init_2,deeper=True,masked=False)
        check_datasets(dset_list_2,deeper=True,masked=False)
    print(" Max. Significance before merging")
    print(f"   - {arrays[0]:6s} = {sigmax(dset_init_1,stacked=False):5.1f}")
    print(f"   - {arrays[1]:6s} = {sigmax(dset_init_2,stacked=False):5.1f}")
    print(" Max. Significance after (maybe) merging")
    print(f"   - {arrays[0]:6s} = {sigmax(dset_list_1, stacked=merging):5.1f}")
    print("   - {arrays[1]:6s} = {sigmax(dset_list_2, stacked=merging):5.1f}")

    nslots = len(dset_list_1)
    if nslots != len(dset_list_2):
        print(" Problem : the two sets are different")

    ###----------------------
    ### Plot excess counts
    ###----------------------
    if counts:
        fig, ax = plt.subplots(nrows=1, ncols=nslots, figsize=(nslots*5.5,6),sharey=True)

        # If model not shown, use steps for the darta points
        ds = "steps-mid" if merging else "default"
        iplot = 0
        lw = 2 if merging else 2

        for ds_1, ds_2, ax  in zip(dset_list_1, dset_list_2, ax):

            elapsed = ds_1.gti.time_start[0] - mc_1.slot.grb.t_trig
            excess_counts(ds_1, elapsed = elapsed,
                          ax=ax, stacked = merging, model_bar = True,
                          n_min = 1., n_max  = float(par[GRB]["n_max"]),
                          emin  = e_min, emax = e_max,
                          color ="tab:blue", alpha = 0.5, fillstyle="full",
                          model_label=None,
                          drawstyle=ds, lw=lw,
                          tag = arrays[0],
                          debug = False)
            elapsed = ds_2.gti.time_start[0] - mc_2.slot.grb.t_trig
            excess_counts(ds_2, elapsed = elapsed,
                          ax=ax, stacked = merging, model_bar = True,
                          n_min = 1., n_max = float(par[GRB]["n_max"]),
                          emin = e_min, emax = e_max,
                          color = "tab:orange", alpha = 0.5, fillstyle="full",
                          model_label=None,
                          drawstyle=ds, lw = lw,
                          tag=arrays[1],
                          debug=False)
            if iplot:
                ax.set_ylabel(None)
            iplot+=1
        fig.suptitle(GRB, fontsize=16)
        plt.tight_layout(w_pad=0)

    ###----------------------
    ### Plot flux and fit
    ###----------------------
    if flux:
        tag1 = "\\"+arrays[0]
        tag2 = "\\"+arrays[1]
        ysize=6 if fit_tag is None else 6
        xsize = nslots*6.5
        lw=3 if merging else 2

        fig, ax = plt.subplots(nrows=1, ncols=nslots,
                               figsize=(xsize,ysize),sharey=True)

        iplot=0
        for ds_1, ds_2, ax  in zip(dset_list_1, dset_list_2, ax):
            elapsed = ds_1.gti.time_start[0] - mc_1.slot.grb.t_trig
            extract_spectrum(ds_1,
                             elapsed = elapsed, stacked=merging,
                             ax = ax, model_style ="bar",
                             color = "tab:blue",  lw=lw,
                             yscale    = "log",
                             alpha_model=0.1, alpha_fit = 0.3, fit_color="purple",
                             e_ref     = u.Quantity(par[GRB]["e_ref"]),
                             e_unit    = "TeV",
                             flux_unit = "TeV-1 cm-2 s-1",
                             flux_min  = 1.e-13,
                             flux_max  = float(par[GRB]["flux_max"]),
                             fit_tag   = fit_tag, tag=tag1,
                             e_min     = 25*u.GeV,
                             e_max    = u.Quantity(par[GRB]["e_fitmax"]),
                             debug = False)
            elapsed = ds_2.gti.time_start[0] - mc_2.slot.grb.t_trig
            extract_spectrum(ds_2,
                             elapsed = elapsed, stacked = merging,
                             ax = ax, model_style ="bar",
                             color = "tab:orange", lw=lw,
                             yscale    = "log",
                             alpha_model = 0.1, alpha_fit = 0.3, fit_color = "red",
                             e_ref     = u.Quantity(par[GRB]["e_ref"]),
                             e_unit    = "TeV",
                             flux_unit = "TeV-1 cm-2 s-1",
                             flux_min  = 1.e-13,
                             flux_max  = float(par[GRB]["flux_max"]),
                             fit_tag   = fit_tag,   tag=tag2,
                             e_min     = e_min, e_max = e_fit_max,
                             debug = False)
            # if fit_tag != None:
            #     ax.legend(bbox_to_anchor=(0.5,-0.8,0.5,0.5),ncol=2,fontsize=12)
            # else:
            #     ax.legend(bbox_to_anchor=(0.5,1.02),loc="upper center",ncol=2,fontsize=12)

            if iplot:
                ax.set_ylabel(None)
            iplot+=1
        fig.suptitle(GRB)
        plt.tight_layout(w_pad=0)

        stamp(GRB, axis=fig)

###############################################################################
if __name__ == "__main__":

    folder = "D:/000_Today/SoHAPPy_tests/single_sources_denser/"
    gidlist = [343, 980, 465, 676, 785]
    for gid in gidlist:
        print(" **************** ",gid," *****************")
        comparison_plot(str(gid), site="South", merging=True, folder=folder)
        plt.show()
