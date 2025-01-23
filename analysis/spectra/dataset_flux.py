# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

@author: Stolar
"""
import sys
import numpy as np


import gammapy

from copy import deepcopy
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.visualization import quantity_support

from scipy.optimize import curve_fit

from niceplot import single_legend

from gammapy.modeling import Fit
from gammapy.modeling.models import SkyModel
from gammapy.estimators import FluxPointsEstimator
from gammapy.modeling.models import ExpCutoffPowerLawSpectralModel
from gammapy.modeling.models import LogParabolaSpectralModel
from gammapy.modeling.models import PowerLawSpectralModel

__all__ = ['extract_spectrum', 'flux_versus_time', 'fluxes_versus_time']


# -----------------------------------------------------------------------------
def fit_model(tag,
              amplitude=1e-12*u.Unit("cm-2 s-1 TeV-1"),
              e_ref=1*u.TeV):
    """
    Define the spectral model for the fit.

    Parameters
    ----------
    tag : String
        A tag for choosing the model.
    amplitude : astropy.Quantity
        The amplitude of the model flux.
    e_ref : astropy.Quantity
        The reference energy of the model flux.
    eblmodel: string
        EBL model name
    z: float
        Source redshift

    Returns
    -------
    gammapy Model
        The model to be fitted.

    """

    model = None

    if tag == "cutoff":
        lmbd_ = 1 * u.Unit("TeV-1")
        model = ExpCutoffPowerLawSpectralModel(index=2.0,
                                               amplitude=amplitude
                                               * np.exp(lmbd_*e_ref),
                                               lambda_=lmbd_,
                                               reference=e_ref,
                                               name=tag)
    elif tag == "powerlaw":
        model = PowerLawSpectralModel(index=2.0,
                                      amplitude=amplitude,
                                      reference=e_ref,
                                      name=tag)
    elif tag == "logparabola":
        model = LogParabolaSpectralModel(alpha=2.0,
                                         beta=0.5,
                                         amplitude=amplitude,
                                         reference=e_ref,
                                         name=tag)
    else:
        sys.exit(" Fit ="+str(tag)+" is not implemented")

    return model


# ## --------------------------------------------------------------------------
# ## ENERGY DSITRIBUTIONS
# ## --------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def get_fluxes(ds, model=None, intrinsic=None, e_unit="GeV", debug=False):
    """
    Obtain flux points from the data and the chosen models

    Parameters
    ----------
    ds : Dataset object
        One dataset.
    model : Gammapy SkyModel, optional
        Sky model. The default is None.
    e_unit : String, optional
        Energy unit. The default is "GeV".

    Returns
    -------
    flx_pts_obs : FluxPoints instance
        Observed flux points
    flx_pts_deabs : FluxPoints instance
        DE EBL-absorbed fux points.

    """

    if model is None:
        sys.exit("A Sky Model is required")

    # ##------------------------------------------
    # ## Get the valid energy range
    # ##------------------------------------------
    # There are n+1 edges for n bins. The safe_mask is for bins.
    # To get the minimal valid energy, I take the edges of the first
    # valid bin
    b_min = np.where(ds.mask_safe.data.flatten())[0][0]
    b_max = np.where(ds.mask_safe.data.flatten())[0][-1] + 2
    e_edges = ds.counts.geom.axes["energy"].edges.to(e_unit)[b_min:b_max]
    e_range = [e_edges[0], e_edges[-1]]
    if debug:
        print(" Energy range = ", e_range)

    # Observed flux
    flx_pts_obs = fit_spectrum(ds, model, e_edges=e_edges)

    # If intrinsic flux is passed (i.e. compound model)
    flx_pts_deabs = None
    if intrinsic is not None:
        flx_pts_deabs = deepcopy(flx_pts_obs)
        flx_pts_deabs._reference_model = \
            SkyModel(spectral_model=intrinsic.copy())

    return flx_pts_obs, flx_pts_deabs


# ----------------------------------------------------------------------------
def fit_and_plot_spectra(ds, model=None, intrinsic=None, e_unit="GeV",
                         ax=None, tag=None, plot=True, sed_type="e2dnde"):
    """


    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    model : TYPE, optional
        DESCRIPTION. The default is None.
    intrinsic : TYPE, optional
        DESCRIPTION. The default is None.
    e_unit : TYPE, optional
        DESCRIPTION. The default is "GeV".
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    tag : TYPE, optional
        DESCRIPTION. The default is None.
    sed_type : TYPE, optional
        DESCRIPTION. The default is "e2dnde".

    Returns
    -------
    ax : TYPE
        DESCRIPTION.

    """

    if ax is None:
        ax = plt.gca()

    flx_obs, flx_intrinsic = get_fluxes(ds,
                                        model=SkyModel(spectral_model=model,
                                                       name=ds.name),
                                        intrinsic=intrinsic)

    pdict = flx_obs.reference_spectral_model.parameters.to_dict()

    index = pdict["name=" == "index"]
    label = (rf"$\gamma = {index['value']:5.2f}"
             rf" \pm {index['error']:5.2f}$")

    e_range = [flx_obs.energy_min[0], flx_obs.energy_max[-1]]
    ax = plot_fitted_spectrum(flx_obs, e_range=e_range,
                              sed_type=sed_type, label=label, ax=ax)
    if flx_intrinsic is not None:
        plot_fitted_spectrum(flx_intrinsic, e_range=e_range,
                             color="red", ax=ax, sed_type=sed_type)
    ax.legend()

    return ax, pdict


# ----------------------------------------------------------------------------
def fit_spectrum(ds, model, e_edges=None,  e_unit="GeV", debug=False):
    """


    Parameters
    ----------
    ds : TYPE
        DESCRIPTION.
    model : TYPE
        DESCRIPTION.
    e_edges : TYPE, optional
        DESCRIPTION. The default is None.
    e_unit : TYPE, optional
        DESCRIPTION. The default is "GeV".
    debug : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    flux : TYPE
        DESCRIPTION.

    """

    ds.models = [model]  # SkyModel

    fit = Fit()
    result = fit.run(datasets=ds)

    # Get the flux at the given energy values
    fpe = FluxPointsEstimator(energy_edges=e_edges,
                              selection_optional="all")
    flux = fpe.run(datasets=ds)

    if debug:
        print(result.models.to_parameters_table())
        print(flux.reference_model)
        print(flux.to_table(sed_type="dnde", formatted=True))

    return flux


# ----------------------------------------------------------------------------
def print_fitted_spectrum(flux_points):
    """
    Print the content of the fitted flux points. This a a way to show
    how to access the data. All elements are ontained with
    print(flux_points)

    Parameters
    ----------
    flux_points : FluxPoints object
        FluxPoints instance.

    Returns
    -------
    None.

    """

    nbins = len(flux_points["npred"].data)
    if gammapy.__version__ < "1.2":
        for i in range(nbins-1):
            print("{:7.2f} - {:7.2f} : {:5.2e} -{:5.2e} +{:5.2e}"
                  .format(flux_points.table["e_min"].quantity[i],
                          flux_points.table["e_max"].quantity[i],
                          flux_points.table["dnde"].quantity[i].value,
                          flux_points.table["dnde_errn"].quantity[i].value,
                          flux_points.table["dnde_errp"].quantity[i]))
    else:
        for i in range(nbins-1):
            print("{:7.2f} - {:7.2f} : {:5.2e} -{:5.2e} +{:5.2e}"
                  .format(flux_points.energy_min[i],
                          flux_points.energy_max[i],
                          flux_points.dnde.data[i][0][0],
                          flux_points.norm.data[i][0][0],
                          flux_points.norm_err.data[i][0][0]))


#  ----------------------------------------------------------------------------
def plot_fitted_spectrum(flux_points, e_range=None, sed_type="e2dnde",
                         color="tab:orange", ax=None, **kwargs):
    """
    Plot flux point and corresponding fitted model and associated errors

    Parameters
    ----------
    flux_points : FluxPoints
        FluxPoints instance.
    e_range : Astropy Qauntity list, optional
        Energy range for the plot. The default is None.
    color : String, optional
        Plot color. The default is "tab:orange".
    ax : Matplotlib Axes, optional
        Current axis. The default is None.

    Returns
    -------
    ax : Matplotlib axes
        Current axis.

    """

    ax = flux_points.plot(ax=ax, sed_type=sed_type, color=color)
    flux_points.reference_spectral_model.plot(ax=ax,
                                              energy_bounds=e_range,
                                              sed_type=sed_type, color=color,
                                              **kwargs)
    flux_points.reference_spectral_model.plot_error(ax=ax,
                                                    energy_bounds=e_range,
                                                    sed_type=sed_type,
                                                    facecolor=color)

    # Some fluxes can be suprinsigly negative (because of excess counts < 0)
    flux_max = np.nanmax(flux_points.e2dnde.data.flatten())
    flux_min = max(0, np.nanmin(flux_points.e2dnde.data.flatten()))
    print(" Flux between ", flux_min, "and", flux_max)
    ax.set_ylim(bottom=max(5.e-14, flux_min), top=min(flux_max, 5.e-10))

    return ax

# ## --------------------------------------------------------------------------
# ## TIME DSITRIBUTIONS
# ## --------------------------------------------------------------------------


# -----------------------------------------------------------------------------
def flux_versus_time(dsets, ax=None,
                     emin=None, emax=None,
                     tmin=None, tmax=None,
                     stacked=False,
                     style="bar", fit=False,
                     e_unit="GeV",
                     flux_min=None,
                     xscale="log",
                     color="tab:blue",
                     model=None,
                     debug=False):
    """
    Plot extracted flux versus time within the given energy range.

    Parameters
    ----------
    dsets : Dataset list
        list of datasets.
    ax : matplotlib.axes, optional
        Current matplotlib axis. The default is None.
    emin : astropy.Quantity, optional
        Minimal energy for extracting the flux. The default is None.
    emax : astropy.Quantity, optional
        Maximal energy for flxu extraction. The default is None.
    tmin : astropy.Quantity, optional
        minimal observation time. The default is None.
    tmax : astropy.Qauntity, optional
        maximal observation time. The default is None.
    stacked : Boolean, optional
        If True, this is a stacked dataset (flux information may be lost).
        The default is False.
    style : String, optional
        The style of th eplot, "bar" or "line". The default is "bar".
    fit : Boolean, optional
        If True perform the fit. The default is False.
    e_unit : astropy.unit, optional
        Energy unit. The default is "GeV".
    flux_min : float, optional
        Minimal flux to be plotted (in the units of the flux extracted).
        The default is 1.e-23.
    xscale : String, optional
        x-scale, "linear" or "log". The default is "log".
    color : String, optional
        The data color. The default is "tab:blue".
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    bool
        True if everything went well.

    """

    t0 = dsets[0].gti.time_start[0]
    ftheory, time, errtime, dnde, errn, errp, erange \
        = [], [], [], [], [], [], []

    # ## -----------------------------------------------
    # ## Extract the flux for each slice
    # ## -----------------------------------------------

    # Loop over time slices, get the flux for the given energy range
    for i, ds in enumerate(dsets):

        # Time and duration of the current slice
        time.append((ds.gti.time_start[0]-t0).sec*u.s + ds.gti.time_sum/2)
        errtime.append(0.5*ds.gti.time_sum)

        # Energy range
        # If Energy boundaries at not given, use the range from the safe mask
        # If energy range is outside the dataset energy range, returns False
        e0 = ds.energy_range[0].data.flatten()[0]*ds.energy_range[0].unit
        e1 = ds.energy_range[1].data.flatten()[0]*ds.energy_range[1].unit
        e_min = e0 if emin is None else max(emin, e0)
        e_max = e1 if emax is None else min(emax, e1)

        if e_min < e0 or e_max > e1 or e_max <= e_min:
            if debug:
                print(" E boundaries out of range")
            return False
        if debug:
            print(f"--{ds.name:}-- Eff. E range: {e_min:5.2f} {e_max:5.2f}")

        # Extracted flux
        flx_obs, _ = get_fluxes(ds,
                                model=SkyModel(spectral_model=model,
                                               name=ds.name),
                                intrinsic=None)

        # Fill the resulting arrays
        if gammapy.__version__ <= "0.18.2":
            e_min = flx_obs.table["e_min"].quantity[0].to(e_unit)
            e_max = flx_obs.table["e_max"].quantity[0].to(e_unit)
            e_ref = flx_obs.table["e_ref"].quantity[0].to(e_unit)
            flx = flx_obs.table["dnde"].quantity[0]
            ep = flx_obs.table["dnde_errp"].quantity[0]
            en = flx_obs.table["dnde_errn"].quantity[0]
        else:
            e_min = flx_obs.energy_min[0].to(e_unit)
            e_max = flx_obs.energy_max[0].to(e_unit)
            e_ref = flx_obs.energy_ref.to(e_unit)
            flx = flx_obs.dnde.data.flatten()[0]*flx_obs.dnde.unit
            ep = 0*flx_obs.dnde.unit  # flx_obs.dnde_errp.flatten.quantity[0]
            en = 0*flx_obs.dnde.unit  # flx_obs.dnde_errn.flatten().quantity[0]
        # Mean theoretical flux at reference energy in the bin
        if stacked:
            fth = flx
        else:
            fth = ds.models[0].spectral_model(e_ref)[0]

        if debug:
            print(f"#{i:2} {e_min.value:6.1f} - {e_max:6.1f}"
                  f" : F= {flx.value:5.1e} +"
                  f"{ep.value:5.1e} - {en.value:5.1e} {flx.unit:s}")
            print(f"   {'':6s} - {'':6s}   T= {fth:5.2e} "
                  f"T/F= {fth.value/flx.to(fth.unit).value:}")

        # Store flux of each slice
        dnde.append(flx)
        errp.append(ep)
        errn.append(en)
        ftheory.append(fth)
        erange.append([e_min, e_max])

    # Move everything to numpy arrays with units
    time = np.asarray([x.value for x in time])*time[0].unit
    errtime = np.asarray([x.value for x in errtime])*errtime[0].unit
    dnde = np.asarray([x.item().value for x in dnde])*dnde[0].unit
    errn = np.asarray([x.item().value for x in errn])*errn[0].unit
    errp = np.asarray([x.item().value for x in errp])*errp[0].unit
    ftheory = np.asarray([f.value for f in ftheory])*ftheory[0].unit
    erange = np.asarray([[e[0].value, e[1].to(e[0].unit).value]*e[0].unit
                         for e in erange])

    # ## -----------------------------------------------
    # ## Plot light curves
    # ## -----------------------------------------------
    ax = plt.gca() if ax is None else ax

    label = f"{e_min.value:5.1f} - {e_max.to(e_min.unit).value:5.1f} "\
            f"{str(e_min.unit):s}"

    with quantity_support():
        ax.errorbar(x=time, y=dnde, xerr=errtime, yerr=[errn, errp],
                    color=color, ls="", label=label)

    # # ## -----------------------------------------------
    # # ## Energy limits
    # # ## -----------------------------------------------
    # axx = ax.twinx()

    # with quantity_support():
    #     eb = axx.errorbar(x=time, y=[e[:][0] for e in erange],
    #                       xerr=errtime, yerr=0,
    #                       color="grey", ls="", alpha=0.5,
    #                       label="$E_{min}, E_{max}$")
    #     eb[-1][0].set_linestyle('dashdot')
    #     eb = axx.errorbar(x=time, y=[e[:][1] for e in erange],
    #                       xerr=errtime, yerr=0,
    #                       color="grey", ls="", alpha=0.5)
    #     eb[-1][0].set_linestyle('dashdot')

    # axx.set_ylabel("Range ("+e_unit+")")
    # axx.set_yscale("log")
    # axx.legend()

    # # ## -----------------------------------------------
    # # ## Fit a t**-beta dependence if required, and plot
    # # ## -----------------------------------------------
    # if fit is True:
    #     # Check existence of data an results
    #     # Remove undefined dnde values
    #     t_min = tmin if (tmin is not None) else min(time)
    #     t_max = tmax if (tmax is not None) else max(time)
    #     mask = (time > t_min) & (time < t_max) & np.isfinite(dnde)

    #     # Check number of slices remaining, at least 3 are necessary
    #     if (len(dsets) - len(mask[mask is False])) <= 3:
    #         print(" Too few slices have flux a estimate for "
    #               "the fit to be performed")
    #     else:
    #         _, _ = fit_flux_versus_time(time[mask], dnde[mask], ax=ax)
    #         ax.axvline(x=t_min, ls="--", color="brown", label="Fit limits")
    #         ax.axvline(x=t_max, ls="--", color="brown")

    # ## -----------------------------------------------
    # ## Display theory if slices were not stacked
    # ## -----------------------------------------------
    if not stacked:
        if style == "bar":
            ax.bar(time, ftheory, width=2*errtime,
                   alpha=0.2, color=color, label="Model")
        elif style == "line":
            ax.plot(time, ftheory,
                    alpha=0.5, color=color, ls=":", marker="o", label="Model")

    # ## Decoration
    with quantity_support():
        if flux_min is not None:
            ax.set_ylim(ymin=flux_min*dnde[0].unit)
        else:
            ax.set_ylim(ymin=0.5*min(dnde))
        ax.set_ylim(ymax=2.*max(dnde))

        if time[-1] > 1*u.d:
            ax.axvline(x=1*u.d, ls=":", color="grey", label="One day")

    # ax.set_xlabel("Elapsed time ("+ax.xaxis.get_label_text()+")")
    # ax.set_yscale("log")
    # ax.set_xscale(xscale)
    # ax.legend(ncol=2)

    return True


# -----------------------------------------------------------------------------
def fit_flux_versus_time(times, fluxes, ax=None):
    """
    Fit a power law on data and plot

    Parameters
    ----------
    times : Astropy Qauntity array
        Abscissae.
    fluxes : Astropy Quantity array
        Values to be fitted.
    ax : matplotlib Axes, optional
        Current axis. The default is None.

    Returns
    -------
    float, float
        Result of the fit.

    """

    # --- Function to be fitted
    def func(x, a, b):
        return a*x**-b   # +c

    # This is valid for Scipy 1.4 (1.9 and later have more output values)
    popt, pcov = curve_fit(func, times, fluxes)
    a = popt[0]
    da = np.sqrt(pcov[0][0])
    b = popt[1]
    db = np.sqrt(pcov[1][1])
#   c = popt[2]
#   dc = np.sqrt(pcov[2][2])

    print(" a=", a, "+/-", da)
    print(" b=", b, "+/-", db)
#   print(" c=",c,"+:-",dc)

    if ax is not None:
        label = r"$\beta$= " + str(round(b, 2)) + r"$\pm$"+str(round(db, 2))
        with quantity_support():
            ax.plot(times, func(times.value, a, b),
                    color="tab:orange", alpha=0.9, label=label)

    return a, b


# -----------------------------------------------------------------------------
def fluxes_versus_time(dsets,
                       xscale="linear",
                       xsize=14, ysize=3,
                       stacked=False):
    """
    Plot the extracted flux on a given energy range, for all available time
    slices.

    Parameters
    ----------
    dsets : List of Dataset
        Current Dataset list.
    tmin : astropy.Quantity, optional
        Minimal plot time. The default is None.
    tmax : astropy.Quantity, optional
        Maximal plot time. The default is None.
    xscale : String, optional
        "log" or "linear". The default is "linear".
    xsize : float, optional
        Figure width. The default is 14.
    ysize : float, optional
        Figure height. The default is 3.
    stacked : Boolean, optional
        If True, the datasets have been stacked. The default is False.

    Returns
    -------
    None.

    """
    print("Consider using tmin and tmax")

    e_edges = dsets[0].excess.geom.axes[0].edges

    for i, _ in enumerate(e_edges[:-1]):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(xsize, ysize))
        status = flux_versus_time(dsets,
                                  emin=e_edges[i], emax=e_edges[i+1],
                                  style="line",
                                  xscale=xscale,
                                  ax=ax,  # color=color,
                                  flux_min=None,
                                  stacked=stacked,
                                  debug=False)
        if not status:
            fig.clear()
        else:
            ax_last = ax
            plt.tight_layout()
            fig.tight_layout(h_pad=0, w_pad=0)

        ax.set_xlabel(None)
        ax.set(xticklabels=[])
        ax.tick_params(bottom=False)

    ax_last.set_xlabel("Elapsed time (" + ax_last.xaxis.get_label_text() + ")")
    ax_last.tick_params(bottom=False)


# -----------------------------------------------------------------------------
def fit_result(fit, result):
    """
    Display the fit results

    Parameters
    ----------
    fit : gammapy fit
        The fit.
    result : gammapt fit result
        The result of the fit.

    Returns
    -------
    text : String
        The results in a readable format.

    """

    mtag = fit.datasets.models[0].spectral_model.tag

    if mtag[0] == "ExpCutoffPowerLawSpectralModel":
        res_index = fit.confidence("index")  # get index error
        res_lambda_ = fit.confidence("lambda_")

        text = (rf"{mtag[1]:}" "\n" r"$\gamma$ = "
                f"${result.parameters['index'].value:3.1f}  $"
                f"$^{{+{res_index['errp']:3.2f}}} $"
                f"$_{{-{res_index['errn']:3.2f}}}$")
        text = text + (rf"\n $\lambda = "
                       f"{result.parameters['lambda_'].value:3.1f}  "
                       f"^{{+{res_lambda_['errp']:3.2f}}} "
                       f"_{{-{res_lambda_['errn']:3.2f}}}$")

    elif mtag[0] == "PowerLawSpectralModel":
        res_index = fit.confidence("index")  # get index error

        text = (rf"{mtag[1]}:" "\n" r"$\gamma$ = "
                f"{result.parameters['index'].value:3.1f}  "
                f"^{{+{res_index['errp']:3.2f}}} "
                f"_{{-{res_index['errn']:3.2f}}}$")

    elif mtag[0] == "LogParabolaSpectralModel":
        res_alpha = fit.confidence("alpha")  # get index error
        res_beta = fit.confidence("beta")

        text = (f"{mtag[1]}:"
                "\n" + r"$\alpha = "
                f"{result.parameters['alpha'].value:3.1f}  "
                f"^{{+{res_alpha['errp']:3.2f}}} "
                f"_{{-{res_alpha['errn']:3.2f}}}$")
        text = text +\
            ("\n" + rf"$\beta = {result.parameters['beta'].value:3.1f} "
             f"^{{+{res_beta['errp']:3.2f}}} "
             f"_{{-{res_beta['errn']:3.2f}}}$")
    else:
        text = ""

    return text


# -----------------------------------------------------------------------------
def fit_spectrum2(ds, fex, e_ref=1*u.TeV, fit_tag=None):
    """
    Fit the spectrum from the extracted flux points

    Parameters
    ----------
    ds : Dataset
        The current dataset.
    fex : FluxPointsEstimator
        Extracted flux point.
    e_ref : astropy.qauntity, optional
        Refernce energy for the model. The default is 1*u.TeV.
    fit_tag : String, optional
        A keyword to select the model to be fitted. The default is None.

    Returns
    -------
    fit : Gammapy Fit
        The Fit instance.
    result : Gammapy Fit result
        The result of the fit.

    """

    # Get a reference point : flux value(xE2) at reference energy
    e_center = ds.background.geom.axes[0].center
    idx = (np.abs(e_center.to(e_ref.unit).value - e_ref.value)).argmin()
    amplitude = (fex.table["ref_dnde"].quantity)[idx]\
        * fex.table["norm"].quantity[idx]
    # if (debug):
    #     print(" Best energy : ",e_center[idx])
    #     print(" Amplitude   : ",amplitude)
    #     ax.scatter(e_center[idx],amplitude*e_center[idx]**2,marker="x",color="red")
    model_fit = fit_model(fit_tag, amplitude, e_ref)

    # Replace default model by the fit model
    ds.models = SkyModel(spectral_model=model_fit, name="Fit to data")

    # Fit the flux points
    fit = Fit(ds)

    minuit_opts = {}
    # minuit_opts = {"tol": 10000, "strategy": 1,"print_level": 2}
    result = fit.run(optimize_opts=minuit_opts)
    # result = fit.optimize(optimize_opts=minuit_opts)
    print(result)

    return fit, result


# -----------------------------------------------------------------------------
def extract_spectrum(ds,
                     elapsed=0*u.s,
                     index=2, flux_unit="TeV-1 cm-2 s-1", e_unit="TeV",
                     ax=None,
                     model_style="bar", color="tab:blue", lw=2,
                     alpha_model=0.2, alpha_fit=0.5, fit_color="blue",
                     xscale="log", yscale="log",
                     theory=True,     stacked=False,
                     fit_tag="cutoff", tag=None,
                     flux_min=None,     flux_max=None,
                     e_min=10*u.GeV, e_max=20*u.TeV,
                     e_ref=1*u.TeV,
                     count_min=10,
                     debug=False):
    """
    Extract the flux from the data in a Dataset

    Parameters
    ----------
    ds : Dataset
        Current dataset.
    elapsed : astropy.Quantity, optional
        Elapsed observation time. The default is 0*u.s.
    index : float, optional
        Index for plotting E**-index*Flux. The default is 2.
    flux_unit : astropy.Quantity, optional
        The flux unit. The default is "TeV-1 cm-2 s-1".
    e_unit : astropy.unit, optional
        Energy unit. The default is "TeV".
    ax : matplotlib.axes, optional
        Current matplotlib axis. The default is None.
    model_style : String, optional
        The model plot style, "bar" or "scatter". The default is "bar".
    color : String, optional
        Flux point color. The default is "tab:blue".
    lw : float, optional
        Line width of the flux plot. The default is 2.
    alpha_model : float, optional
        model transparency. The default is 0.2.
    alpha_fit : float, optional
        Fit plot transparency. The default is 0.5.
    fit_color : String, optional
        Fit plot color. The default is "blue".
    xscale : String, optional
        "log" or "linear". The default is "log".
    yscale : String, optional
        '"log" or "linear". The default is "log".
    theory : Boolean, optional
        If True, plot the theory if possible. The default is True.
    stacked : Boolean, optional
        If True, this is a stacked dataset (theory may be wrong).
        The default is False.
    fit_tag : String, optional
        A tag giving the model to be fitted. The default is "cutoff".
    tag : String, optional
        A text tagging the label. The default is None.
    flux_min : float, optional
        Minimal flux in unit of retrieved flux. The default is None.
    flux_max : float, optional
        Maximal flux in unit of retrieved flux. The default is None.
    e_min : astropy.Quantity, optional
        Minimum energy. The default is 10*u.GeV.
    e_max : astropy.Quantity, optional
        Maximal energy. The default is 20*u.TeV.
    e_ref : astropy.Qauntity, optional
        Reference energy for the fit model. The default is 1*u.TeV.
    count_min : float, optional
        Count number below whihc the fit is not performed. The default is 10.
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    (float, float)
        The minimal and maximal flux values (used by dataset_plot.panels).

    """

    ax = plt.gca() if ax is None else ax

    # It can be long, let's speak a bit
    print(" ------------------------ ", ds.name,
          " -----------------------------")

    # If not enough count, do not attempt to extract spectrum
    # count_max=max(ds.excess.data.flatten())
    # if count_max <= count_min:
    #     ax.text(0.05,0.95,
    #             "Counts too low ("+str(round(count_max,1))+")",
    #             color=color,
    #             transform=ax.transAxes)
    #     print(" Counts too low")
    #     return (-np.nan, np.nan)

    # Define label
    # if tag is None or len(tag) == 0 :
    #     label="$t$=" + t_str(elapsed) \
    #       + r"\n$\Delta$t="+t_str(ds.gti.time_sum)
    # else:
    #     label="$t_{"+tag+"}$=" + t_str(elapsed) +"\n" \
    #       + r"$\Delta t_{"+tag+"}$="+t_str(ds.gti.time_sum)

    # If the datasets were not stacked, display the original spectrum as bars
    # otherwise link measurement points with a dashed line
    ls = ":" if stacked else ""

    # Reconstrcuted energy axis
    # axis_reco=ds.background.geom.axes[0]
    # e_center =axis_reco.center.to(e_unit)
    # e_edges  =axis_reco.edges.to(e_unit)

    # ##---------------------
    # ## Plot extracted flux
    # ##---------------------
    # To extract the flux in the dataset E range use
    # e_edges[ (e_edges>=ds.energy_range[0]) & (e_edges<=ds.energy_range[1])]
    fex = extract_flux_points(ds, debug=debug)
    fex.plot(energy_power=index,
             flux_unit=flux_unit,
             energy_unit=e_unit,
             ax=ax, ls=ls, label=label, lw=lw, color=color)

    # ##---------------------
    # ## Plot theory (at reconstructed energies !)
    # ##---------------------
    if theory and not stacked:  # Stacked dataset have no proper spectrum
        model = ds.models[0].spectral_model
        ftheory = model(e_center).to(flux_unit)*e_center.to(e_unit)**index
        width_r = e_edges[1:] - e_center
        width_l = e_center - e_edges[:-1]

        with quantity_support():
            if model_style == "bar":
                ax.bar(e_center.to(e_unit), ftheory,
                       width=-width_l, align="edge",
                       alpha=alpha_model, color=color,
                       label=ds.name+" model")
                ax.bar(e_center.to(e_unit), ftheory,
                       width=width_r, align="edge",
                       alpha=alpha_model, color=color)
            if model_style == "scatter":
                ax.scatter(e_center.to(e_unit), ftheory, marker="o",
                           facecolor="white", edgecolors="red",
                           label=ds.name+" model")

    # ##---------------------
    # ## Fit the data and plot
    # ##---------------------
    # MINUIT reference:
    # https://iminuit.readthedocs.io/en/latest/reference.html

    if fit_tag is not None:
        # Copy the initial model to reset the dataset in the end
        model_init = ds.models
        fit, result = fit_spectrum(ds, fex, e_ref=e_ref, fit_tag=fit_tag)

    # ##---------------------
    # ## Plot the fit and the result text
    # ##---------------------
        # Result of fit text

        # If fit failed, just mention it on the plot
        if not result.success or result.parameters["amplitude"].value < 0:
            ax.text(0.05, 0.95, "FIT FAILED ",
                    color=color,
                    transform=ax.transAxes)
            print("Fit failed!")
        else:
            fit_text = fit_result(fit, result)
            print(fit_text)
            if debug:
                print(" FITTING : ", fit_tag)
                print("*** ", result)
                print("*** ", result.parameters)
                print(fit_text)
                print(ds.models[0].spectral_model)

        # model_fit.plot_error(ax=ax,
        ds.models[0].spectral_model.plot_error(ax=ax,
                                               energy_range=ds.energy_range,
                                               energy_power=index,
                                               energy_unit=e_unit,
                                               alpha=0.2,
                                               color=fit_color)

        # Note that ds.models[0].spectral_model.plot also works
        # model_fit.plot(ax=ax,energy_range    =ds.energy_range,
        if tag is not None and len(tag) != 0:
            mlabel = "$" + tag + "$" + " " + fit_tag
        else:
            mlabel = fit_tag
        ds.models[0].spectral_model.plot(ax=ax,
                                         energy_range=ds.energy_range,
                                         energy_power=index,
                                         energy_unit=e_unit,
                                         flux_unit=flux_unit,
                                         label=mlabel,
                                         alpha=alpha_fit,
                                         color=fit_color)
        ds.models = model_init

    # Decoration
    with quantity_support():
        if not stacked:  # Plot energy range
            ax.axvline(ds.energy_range[0], ls=":", color="grey")
            ax.axvline(ds.energy_range[1], ls=":", color="grey")

    # This does not work with units within quantity_support
    ax.set_xlim(xmin=e_min.to(e_unit).value)
    ax.set_xlim(xmax=e_max.to(e_unit).value)
    if flux_min is not None:
        flux_min = flux_min*u.Unit(flux_unit)*u.Unit(e_unit)**index
        ax.set_ylim(ymin=flux_min.value)
    else:
        flux_min = -np.Inf

    if flux_max is not None:
        flux_max = flux_max*u.Unit(flux_unit)*u.Unit(e_unit)**index
        ax.set_ylim(ymax=flux_max.value)
    else:
        flux_max = np.Inf

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid("both")
    single_legend(ax, fontsize=12)

    return (flux_min, flux_max)  # expect ymin, ymax
