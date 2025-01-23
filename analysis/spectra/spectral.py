# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:50:34 2023

This script shows how to analyse the time and energy spectra of a simulated
GRB.
Note that in some cases the available data are unsufficient to get a
spectrum fitted, in particular when observation starts early and the time
slices are very short. In that case it is recommened to stack the slices within
a minimal duration.

@author: Stolar
"""

import sys
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time
from astropy.visualization import quantity_support

from gammapy.modeling.models import EBLAbsorptionNormSpectralModel
from gammapy.modeling.models import SkyModel

from gammapy.estimators import LightCurveEstimator

from dataset_plot import windows, panels
from dataset_tools import stacking, check_datasets

from dataset_flux import (get_fluxes, fit_and_plot_spectra,
                          fit_model, plot_fitted_spectrum,
                          print_fitted_spectrum)

from dataset_flux import (flux_versus_time)

from niceprint import heading
from niceplot import stamp

from utils import datasets_from_binary

sys.path.append("../")
sys.path.append("../../../SoHAPPy")


# #############################################################################
class FitModels():
    
    # -----------------------------------------------------------------------------
    def __init__(self, datas, e_model, abs_model):
        
        self.e_model = e_model
        self.abs_model = abs_model
        self.dsets = datas
        
    # -----------------------------------------------------------------------------
    def cumulated_energy_spectra(self, debug=False):
        
        ds_stacked = stacking(self.dsets)
        
        # Get flux points
        model = SkyModel(spectral_model=self.abs_model, name=GRB_id)
        flx_obs, flx_intrinsic = get_fluxes(ds_stacked[-1],
                                            model=model,
                                            intrinsic=self.e_model)
        if debug:
            print(flx_obs.reference_model)
            print(flx_intrinsic.reference_model)
        
        # Get flux parameters - note that indices are the same with or without
        # deabsoprtion
        pdict = flx_obs.reference_spectral_model.parameters.to_dict()
        
        index = pdict["name=" == "index"]
        
        # Plot results
        
        label = (rf"$\gamma = {index['value']:5.2f}"
                 rf" \pm {index['error']:5.2f}$")
        
        e_range = [flx_obs.energy_min[0], flx_obs.energy_max[-1]]
        ax = plot_fitted_spectrum(flx_obs, e_range=e_range,
                                  sed_type="dnde", label=label)
        plot_fitted_spectrum(flx_intrinsic, e_range=e_range,
                             color="red", ax=ax, sed_type="dnde")
        ax.legend()
        
        
        # Print flux values at sampled energy
        if debug:
            print_fitted_spectrum(flx_obs)
            print_fitted_spectrum(flx_intrinsic)
        
    # -----------------------------------------------------------------------------
    def individual_energy_spectra(self, debug=False):
        
        # Fit a powerLaw with expoential cutoff
        fig = panels(self.dsets, fit_and_plot_spectra,
                     tref=grb.t_trig, nmaxcol=4,
                     model=fit_model("cutoff"), intrinsic=None)
        fig.suptitle("Individual energy spectra")
        plt.tight_layout()
        
        # Fit an absorbed powerla
        fig = panels(dsets, fit_and_plot_spectra,
                     tref=grb.t_trig, nmaxcol=4,
                     model=model, intrinsic=spec_model)
        
        fig.suptitle("Individual energy spectra")
        plt.tight_layout()
        
    # -----------------------------------------------------------------------------
    def energy_index(self, t0 = 0, debug=False):
        
        # Fit a powerLaw with expoential cutoff
        indices = []
        index_errors = []
        times = []
        errtimes = []
        first = True
        t0 = 0
        
        for ds in dsets:
            if first:
                t0 = ds.gti.time_start[0]
                t_unit = u.s
                first = False
                
            # Time intevals  
            err_tds = ds.gti.time_sum/2
            tds = (ds.gti.time_start[0]-t0).sec*u.s + err_tds
            
            # Fit flux
            model = SkyModel(spectral_model=self.abs_model, name=ds.name)
            
            flx_obs, flx_intrinsic = get_fluxes(ds,
                                                model=model,
                                                intrinsic=self.e_model)
                  
            # intrinsic indices
            pdict = flx_intrinsic.reference_spectral_model.parameters.to_dict()
        
            indices.append(pdict[0]["value"]) 
            index_errors.append(pdict[0]["error"])
            
            times.append(tds.to(t_unit).value)
            errtimes.append(err_tds.to(t_unit).value)    
            
        indices = np.array(indices)
        index_errors = np.array(index_errors)
        times = np.array(times)
        errtimes = np.array(errtimes)
        
        idx_mean = [np.mean(indices), np.std(indices)]
        
        fig, ax = plt.subplots(figsize=(12, 5))
        
        with quantity_support():  # errorbar do not suport units

            ax.errorbar(x=times, y=indices,
                        xerr=errtimes, yerr=index_errors,
                        color="tab:blue", ls="", label="Intrinsic")
            
            ax.set_xlabel("Elapsed time ("+ str(t_unit)+")")
            ax.set_ylabel("Energy index")
            ax.grid("both")
            ax.axhline(idx_mean[0], ls="--")
            ax.axhspan(idx_mean[0] - idx_mean[1],
                       idx_mean[0] + idx_mean[1],
                       alpha=0.2, color="tab:blue")
            ax.legend()
            
        print(f" Mean index  = {idx_mean[0]:5.2f} +/- {idx_mean[1]:5.2f}")
        print(f" Mean errors = {np.mean(index_errors):5.2f}")
            
        return ax
        
###############################################################################
# Get data from a GRB simulated file
# Use the INFILE environment variable if it exists

# This is required to have the EBL models read from gammapy
os.environ['GAMMAPY_DATA'] =\
  str(Path(Path(__file__).absolute().parent.parent.parent, "data"))

# File to be analysed
base = Path(r"D:\CTAO\SoHAPPy\HISTORICAL\180720B\Night_1")
GRB_id = "180720B"
site = "South"
filepath = Path(base, GRB_id+"-"+site+"_sim.bin")

# ## -----------------------------------
# ## Get back the datasets and the GRB
# ## -----------------------------------
heading("DATASETS & GRB")
dsets, grb = datasets_from_binary(filepath=filepath)

# dsets = dset_init.copy()  # Why do I copy. Can't remember...
# Display count number of each energy bin
# check_datasets(dsets, masked=False, deeper=True, header=True)

print(grb)
# grb.plot()

# ## -----------------------------------
# ##   Energy spectra
# ## -----------------------------------

#  Define the emodel
# To get the fit with no EBLabsoprtion, create a compound model
spec_model = fit_model("powerlaw")
model = spec_model * \
            EBLAbsorptionNormSpectralModel.read_builtin("dominguez",
                                                        redshift=grb.z)
# print(model)
fit = FitModels(dsets, spec_model, model)

# Cumulated Energy spectra
heading("Cumulated Energy spectra")
fit. cumulated_energy_spectra()

##  Individual Energy spectra
heading("Individual Energy spectra")
# fit.individual_energy_spectra()

# Enegy indices
heading("Energy indices with time")
# fit.energy_index()

# ## ------------------
# ##  Light curve in flux
# ## ------------------

# Plot the flux time evolution for a given range
# (limited to the reco axis range)
heading("Flux versus time")
fig, ax = plt.subplots(figsize=(20, 8))
flux_versus_time(dsets,
                 emin=100*u.GeV, emax=5*u.TeV,
                 tmin=300*u.s, tmax=2*u.d,
                 stacked=False, fit=True, debug=True,
                 model=model)
stamp(filepath.stem)
# debug = True


# e_unit = "GeV"
# b_min = np.where(dsets[0].mask_safe.data.flatten())[0][0]
# b_max = np.where(dsets[0].mask_safe.data.flatten())[0][-1] + 2
# e_edges = dsets[0].counts.geom.axes["energy"].edges.to(e_unit)[b_min:b_max]
# e_range = [e_edges[0], e_edges[-1]]
# print(" Energy range = ", e_range)

# # For unknown reasons,
# dts = []
# for ids, ds in enumerate(dsets):
#     tstart = ds.gti.time_start
#     tstop = ds.gti.time_stop
#     # if ids != 0:
#     #     tstart = tstart + 1e-6*u.s
#     dts.append(Time([tstart.mjd[0], tstop.mjd[0]], format="mjd", scale="utc"))
#     if debug:
#         print(" - ", ids)
#         print("   start ", tstart)
#         print("   stop  ", tstop)

# lc_maker_1d = LightCurveEstimator(
#                                    energy_edges=e_range,
#                                    time_intervals = dts,
#                                    reoptimize=False,
#                                    selection_optional="all",
#                                    )
# lc_1d = lc_maker_1d.run(dsets)
# print(lc_1d.to_table(sed_type="flux", format="lightcurve"))
# fig, ax = plt.subplots()
# lc_1d.plot(ax=ax, marker="o", label="1D", sed_type="e2dnde")










