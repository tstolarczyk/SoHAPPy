# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

@author: Stolar
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

import astropy # for type
from astropy.time import Time
from astropy.visualization import quantity_support
from gammapy.estimators import FluxPointsEstimator

import gammapy
from utilities import t_str


__all__ = ['extract_spectrum','excess_counts','residuals',
           'excess_versus_time','stacked_versus_time','flux_versus_time',
           'fluxes_versus_time','excess_versus_E_and_time','lightcurve',
           'panels', 'windows']

###############################################################################
### Models
###############################################################################
def models_plot(dsets):

    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(5,4))
    for i,ds in enumerate(dsets):
        spectrum = ds.models[0].spectral_model
        spectrum.plot(energy_range=ds.energy_range,
                      flux_unit='cm-2 s-1 erg-1',
                      energy_power=2,
                      energy_unit='TeV',
                      n_points=10,
                      ls=":",marker="",
                      ax=ax,
                      label=str(i))
        #print(spectrum)
    ax.legend()
###############################################################################
### Time windows
###############################################################################
def windows(dsets, nbin=25, ysize=5, unit="s", plot=True):
    """
    In a stacked dataset, GTI are merged (only) if they are contiguous
    """
    for i, ds in enumerate(dsets):
        dt = [t_str(t) for t in ds.gti.time_delta]
        print(" Set {:3d}: sum:{:>10} | start={}, stop={}, dt={}".format(i,
                                                                   t_str(ds.gti.time_sum),
                                                                   ds.gti.time_start,
                                                                   ds.gti.time_stop,
                                                                   dt))

    if plot:
        import matplotlib.cm as cm

        with quantity_support():

            fig, (ax,ax1) = plt.subplots(nrows=2,ncols=1,figsize=(10,2*ysize))
            
            ax.hist([(ds.gti.time_delta.to(unit).value) for ds in dsets],
                    bins=nbin)
            dt_short = np.array([(ds.gti.time_delta.to(unit).value) for ds in dsets])
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Occurences")
            ax1.set_xlabel("Observation periods")
            tmin = np.concatenate(np.array([(ds.gti.time_start).jd for ds in dsets]))
            tmax = np.concatenate(np.array([(ds.gti.time_stop).jd for ds in dsets]))
            it = 0
            for t1,t2 in zip(tmin,tmax):
#                     print(t1,t2)
                color = cm.cool(it/len(tmin))
                ax1.axvspan(xmin=t1, xmax=t2, ymin=0, ymax=1-it/len(tmin),color=color,alpha=0.5,label=str(it))
                it+=1
            if len(tmin)<6:
                ax1.legend()

            dt_short = np.concatenate(dt_short)
            dt_short = dt_short[dt_short<500]
            if (len(dt_short)):
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                axx = inset_axes(ax, width="50%", height=1.2,loc="upper right")
                axx.hist(dt_short,bins=nbin)
                axx.set_xlim(xmax=200)
                axx.set_xlabel("Time (s)")
        plt.tight_layout()
    return fig

###############################################################################
### Counts
###############################################################################
def excess_counts(ds,ax   = None,
                  elapsed = 0*u.s,
                  stacked = False, bar = True, debug = False,
                  alpha   = 0.2, min_counts = 0.1, 
                  emin    = 10*u.GeV, emax=20*u.TeV,
                  **kwargs):

    """
    This is a modified version of dataset plot_counts method.
    It does not explicitely mask the energy bins
    """
    
    ax = plt.gca() if ax is None else ax
    label = "$t$ = " + t_str(elapsed) \
          + "\n$\Delta$t = "+t_str(ds.gti.time_sum)

    ls=":" if stacked else ""

    ds.excess.plot(ax=ax, label=label,ls=ls)
    max_counts = np.nanmax(ds.counts.data).item()
    
    # Plot theory if possible
    if (not stacked): 
        npred      = ds.npred_signal().data.flatten()
        max_counts = np.max([max_counts, np.nanmax(npred).item()])

        if (not bar): # Use default function
            p   = ds.npred_signal().plot(ax=ax, label="mu_src")
            clr = p.get_lines()[1].get_color()

        else: # Use better function
            # Add bar for better reading
            # in log scale center is not center, widths are assymetric
            clr     = "tab:orange"

            Eedges  = ds.background.geom.axes[0].edges.flatten()
            Ecenter = ds.background.geom.axes[0].center.flatten()
            width_r = Eedges[1:]-Ecenter
            width_l = Ecenter - Eedges[:-1]

            ax.bar(Ecenter, npred, width = -width_l,
                   alpha = alpha, align="edge", color=clr, label="$\mu_{src}$")
            ax.bar(Ecenter, npred, width = width_r,
                   alpha = alpha, align="edge", color=clr)
            
        ds._plot_energy_range(ax=ax) # Show the dataset energy masking

    ax.set_xlim(emin.to(Eedges[0].unit).value,emax.to(Eedges[0].unit).value)
    ax.grid("both",ls="--",alpha=0.5)
    ax.legend(numpoints=1)
    ax.set_ylabel("Excess count number")
    ax.set_ylim(ymin=min_counts)
    ax.set_title("")

    return (min_counts, max_counts)
#------------------------------------------------------------------------------
def residuals(ds,ax=None,**kwargs):
    ds.plot_residuals(ax=ax)
    return

#------------------------------------------------------------------------------
def excess_versus_time(dsets, rate=False, unmasked=False, ax=None,
                       xscale="log", yscale="log",
                       color1 = "tab:blue",color2 = "tab:orange",
                       debug  = False):
    """
    Plot both masked and unmasked count evolutions

    Parameters
    ----------
    dsets : TYPE
        DESCRIPTION.
    rate : TYPE, optional
        DESCRIPTION. The default is True.
    unmasked : TYPE, optional
        DESCRIPTION. The default is False.
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    xscale : TYPE, optional
        DESCRIPTION. The default is "linear".
    yscale : TYPE, optional
        DESCRIPTION. The default is "linear".
    color1 : TYPE, optional
        DESCRIPTION. The default is "tab:blue".
    color2 : TYPE, optional
        DESCRIPTION. The default is "tab:orange".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    ax = plt.gca() if ax is None else ax

    with quantity_support():

        t0 = dsets[0].gti.time_start[0]
        tmax = 0*u.d

        for i,ds in enumerate(dsets):

            time      = (ds.gti.time_start[0]-t0).sec*u.s + ds.gti.time_sum/2
            errtime   = 0.5*ds.gti.time_sum
            if gammapy.__version__ == "0.17":
                npred     = ds.npred_sig().data.sum()
                npred_msk = ds.npred_sig().data[ds.mask_safe].sum()
            else: # 0.18.2
                npred     = ds.npred_signal().data.sum()
                npred_msk = ds.npred_signal().data[ds.mask_safe].sum()
            tmax      = max(tmax,time)
            xs        = ds.excess.data.sum()
            xs_msk    = ds.excess.data[ds.mask_safe].sum()

            norm      = 1/2/errtime if rate else 1*u.dimensionless_unscaled

            if (debug):
                print("{:2} : t = [{:10.2f} {:10.2f}] - y={:8.2f} ym={:8.2f} / th={:8.2f} thm={:8.2f}"
                      .format(i,(time-errtime).value, time+errtime,
                              xs*norm,xs_msk*norm,npred*norm,npred_msk*norm))
                ax.text(time,
                        max(1.2*abs(xs_msk)*norm,10*norm),
                        str(i))

            # Handling of quantities by Errorbar is strange
            # - requires numpy arrays - quantities are not scalars !
            ax.errorbar(x = [time.value]*time.unit,
                        y = [xs_msk*norm.value]*norm.unit,
                        xerr = [errtime.value]*errtime.unit,
                        yerr = [np.sqrt(xs_msk)*norm.value]*norm.unit,
                        label="Excess counts",color=color1)
            ax.bar( time,xs_msk*norm,width=2*errtime,alpha=0.2,color=color1)

            # Before masking
            if (unmasked):
                ax.errorbar(x     = [time.value]*time.unit,
                            y     = [xs*norm.value]*norm.unit,
                            xerr  = [errtime.value]*errtime.unit,
                            yerr  = [np.sqrt(xs)*norm.value]*norm.unit,
                            label ="Total excess counts",color=color2)
                ax.bar(time,npred*norm,width=2*errtime,alpha=0.2,color=color2)

        ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
        if rate:
            ax.set_ylabel("Excess Count rate (" + ax.yaxis.get_label_text()+ ")")
        else:
            ax.set_ylabel("Excess counts")
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if tmax>1*u.d: ax.axvline(x=1*u.d,ls=":",
                                      color="grey",label="One day")

        import collections
        handles, labels = ax.get_legend_handles_labels()
        by_label = collections.OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    return
#------------------------------------------------------------------------------
def excess_versus_E_and_time(dsets):

    """


    Parameters
    ----------
    dsets : Datasets object
        A collection of Gammapy datasets

    Returns
    -------
    None.

    """
    from dataset_tools import get_axis
    dt, E = get_axis(dsets, tunit=u.h, Eunit=u.GeV)

    if gammapy.__version__ == "0.17":
        data_pred   = np.asarray([ds.npred_sig().data.flatten() for ds in dsets])
    else: #0.18.2
        data_pred   = np.asarray([ds.npred_signal().data.flatten() for ds in dsets])

    data_excess = np.asarray([ds.excess.data.flatten() for ds in dsets])

    fig = plt.figure(figsize=(20,8))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    plot3D(dt, E, data_pred,  tag="Prediction",ax3d = ax1)
    plot3D(dt, E, data_excess, tag="Excess counts", ax3d = ax2)

    plt.tight_layout()
    plt.show()

    return

###############################################################################
### Fluxes
###############################################################################
from gammapy.modeling.models import PowerLawSpectralModel
from gammapy.modeling.models import SkyModel

def extract_flux_points(ds, e_edges, index=2, sigm_ul=2, debug=False):
    """
    The flux is extracted from the execess number assimg a certin model
    The difference in results between two models is larger when the energy binning is larger.
    """

# Replace default model by the fit model - Change for an E-2 model
# It was checked that it changes slighlty the extracted flux values
# The effect is larger for larger E bins
    model_init = ds.models # Will be put back in place later
    model_fit = PowerLawSpectralModel(
        index     =   2,
        amplitude = 1e-13 * u.Unit("cm-2 s-1 GeV-1"),
        reference =     1000 * u.GeV,name="pl")
    ds.models = SkyModel(spectral_model=model_fit, name="Fit to data")

    # Extact the flux assuming the generic model
    if gammapy.__version__== "0.17":
        fpe = FluxPointsEstimator(e_edges=e_edges, norm_min=0.2, norm_max=5, norm_n_values=11,
                                                             norm_values=None,
                                                             sigma=1, sigma_ul=2, reoptimize=False )
    else: # 0.18.2
        fpe = FluxPointsEstimator(energy_edges=e_edges, norm_min=0.2, norm_max=5, norm_n_values=11,
                                                             norm_values=None,
                                                             n_sigma=1, n_sigma_ul=2, reoptimize=False )
    fex = fpe.run(datasets = ds)
    fex.table["is_ul"] = fex.table["ts"] < sigm_ul**2

    # Put the original model back in place
    ds.models = model_init

    if (debug):
        print(" --- ",ds.name," --- assuming a powerlaw with index = 2")
        for i in range(len(e_edges)-1):
            print("{:7.2f} - {:7.2f} : {:5.2e} -{:5.2e} +{:5.2e}"
                  .format(fex.table["e_min"].quantity[i],
                          fex.table["e_max"].quantity[i],
                          fex.table["dnde"].quantity[i].value,
                          fex.table["dnde_errn"].quantity[i].value,
                          fex.table["dnde_errp"].quantity[i]))

    return fex
#------------------------------------------------------------------------------
def extract_spectrum(ds,
                     elapsed = 0*u.s,
                     index = 2, flux_unit = "GeV-1 cm-2 s-1", e_unit = "GeV",
                     ax = None, style ="bar", color = "tab:blue",
                     xscale   = "linear", yscale   = "log",
                     theory   = True,     stacked  = False, fit_tag="cutoff",
                     flux_min = None,     flux_max = None,
                     e_min    = 10*u.GeV, e_max    = 20*u.TeV, 
                     e_ref    = 1000*u.GeV,
                     count_min=10,
                     debug = False):
    """
    Extract the flux of one dataset
    Note that the flux is extracted on the masked dataset.
    """

    ax = plt.gca() if ax is None else ax

    # It can be long, let's speak a bit
    print(" ------------------------ ",ds.name," -----------------------------")

    # If not enough count, do not attempt to extract spectrum
    count_max = max(ds.excess.data.flatten())
    if count_max <= count_min:
        ax.text(0.5,0.5,"Counts too low ("+str(round(count_max,1))+")",transform=ax.transAxes)
        print(" Counts too low")
        return (0, 0)

    # Plot label : Elpased time since a reference, livetime
    label = "$t$ = " + t_str(elapsed) \
              + "\n$\Delta$t = "+ t_str(ds.gti.time_sum)
    # If the datasets were not stacked, display the original spectrum as bars
    # otherwise link measurement points with a dashed line
    ls = ":" if stacked else ""

    # Reconstrcuted energy axis
    axis_reco = ds.background.geom.axes[0]
    e_center  = axis_reco.center.to(e_unit)
    e_edges   = axis_reco.edges.to(e_unit)

    ###---------------------
    # Plot extracted flux
    ###---------------------
#      To extract the flux in the dataset E range use
#      e_edges[ (e_edges>=ds.energy_range[0]) & (e_edges<=ds.energy_range[1])]
    fex = extract_flux_points(ds, e_edges,debug=debug)
    fex.plot(energy_power = index,
             flux_unit    = flux_unit,
             energy_unit  = e_unit,
             ax = ax, ls = ls, label = label)

    ###---------------------
    # Plot theory (at reconstructed energies !)
    ###---------------------
    if theory and not stacked: # Stacked dataset have no proper spectrum
        model   = ds.models[0].spectral_model
        ftheory = model(e_center).to(flux_unit)*e_center.to(e_unit)**index
        width_r = e_edges[1:] -e_center
        width_l = e_center - e_edges[:-1]

        with quantity_support():
            if style== "bar":
                ax.bar(e_center.to(e_unit), ftheory, width = -width_l, align="edge",
                       alpha=0.2, color = color, label=ds.name+" model")
                ax.bar(e_center.to(e_unit), ftheory, width = width_r, align="edge",
                       alpha=0.2, color = color)
            if style == "scatter":
                ax.scatter(e_center.to(e_unit), ftheory, marker="o",
                           facecolor="white", edgecolors="red",
                           label=ds.name+" model")

    ###---------------------
    # Fit the data and plot
    ###---------------------
    # MINUIT reference : https://iminuit.readthedocs.io/en/latest/reference.html

    # Get a reference point : flux value(xE2) at reference energy
    idx = (np.abs(e_center.to(e_ref.unit).value - e_ref.value)).argmin()
    amplitude = amplitude = (fex.table["ref_dnde"].quantity)[idx] * fex.table["norm"].quantity[idx]
    if (debug):
        print(" Best energy : ",e_center[idx])
        print(" Amplitude   : ",amplitude)
        ax.scatter(e_center[idx],amplitude*e_center[idx]**2,marker="x",color="red")

    if (fit_tag != None):
        # Copy the initial model to reset the dataset in the end
        model_init = ds.models

        # Create the appropriate model
        if fit_tag=="cutoff":
            from gammapy.modeling.models import ExpCutoffPowerLawSpectralModel
            lambda_ = 1 * u.Unit("TeV-1")
            model_fit = ExpCutoffPowerLawSpectralModel( index     =   2.0,
                                                        amplitude = amplitude*np.exp(lambda_*e_ref), # Correct phi0 for cutoff
                                                        lambda_   =   lambda_,
                                                        reference =   e_ref)
        elif fit_tag=="powerlaw":
            from gammapy.modeling.models import PowerLawSpectralModel
            model_fit = PowerLawSpectralModel( index     =   2.0,
                                                amplitude = amplitude,
                                                reference = e_ref,name="pl")
        elif fit_tag=="logparabola":
            from gammapy.modeling.models import LogParabolaSpectralModel
            model_fit = LogParabolaSpectralModel(alpha     =   2.0,
                                                 beta      = 0.5,
                                                 amplitude = amplitude,
                                                 reference = e_ref,name="lpl")
        else:
            import sys
            sys.exit(" Fit = "+str(fit_tag)+" is not implemented")

        # Plot the model starting point before fitting
#         if debug:
#             model_fit.plot(ax=ax,energy_range=ds.energy_range,
#                      energy_power =index,
#                      energy_unit = e_unit,
#                      flux_unit   = flux_unit,
#                       label="Start",
#                       ls=":")

        # Replace default model by the fit model
        from gammapy.modeling.models import SkyModel
        ds.models = SkyModel(spectral_model=model_fit, name="Fit to data")

        # Fit the flux points
        from gammapy.modeling import Fit
        fit = Fit(ds)

        minuit_opts = {}
        # minuit_opts = {"tol": 10000, "strategy": 1,"print_level": 2}
        result = fit.run(optimize_opts=minuit_opts)
        # result = fit.optimize(optimize_opts=minuit_opts)
        print(result)

        # If fit failed, just mention it on the plot
        if not result.success or result.parameters["amplitude"].value<0:
            ax.text(0.5,0.5,"FIT FAILED ",color="red",transform=ax.transAxes)
            print("Fit failed!")
        else:
            if fit_tag=="cutoff":
                res_index   = fit.confidence("index") # get index error
                res_lambda_ = fit.confidence("lambda_")
                label = "{}:\n$\gamma = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                    .format(fit_tag,result.parameters["index"].value, \
                    res_index["errp"],res_index["errn"])
                label = label+"\n $\lambda = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                        .format(result.parameters["lambda_"].value, \
                            res_lambda_["errp"],res_lambda_["errn"])

            elif fit_tag == "powerlaw":
                res_index   = fit.confidence("index") # get index error
                label = "{}:\n$\gamma = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                    .format(fit_tag,result.parameters["index"].value, \
                    res_index["errp"],res_index["errn"])

            elif fit_tag == "logparabola":
                res_alpha   = fit.confidence("alpha") # get index error
                res_beta    = fit.confidence("beta")
                label = "{}:".format(fit_tag) +"\n"+r"$\alpha = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                    .format(result.parameters["alpha"].value, \
                    res_alpha["errp"],res_alpha["errn"])
                label = label+"\n"+r"$\beta = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                        .format(result.parameters["beta"].value, \
                            res_beta["errp"],res_beta["errn"])
            if (debug):
                print(" FITTING : ",fit_tag)
                print("*** ",result)
                print("*** ",result.parameters)
                print(label)
                print(ds.models[0].spectral_model)

            model_fit.plot_error(ax=ax,
                                 energy_range = ds.energy_range,
                                 energy_power = index,
                                 energy_unit  = e_unit)

            # Note that ds.models[0].spectral_model.plot also works
            model_fit.plot(ax=ax,energy_range     = ds.energy_range,
                                     energy_power = index,
                                     energy_unit  = e_unit,
                                     flux_unit    = flux_unit,
                                     label        = label)

        ds.models = model_init # reset models
   
    with quantity_support():
        if (not stacked):  # Plot energy range
            ax.axvline(ds.energy_range[0],ls=":",color="grey")
            ax.axvline(ds.energy_range[1],ls=":",color="grey")
            
        ax.set_xlim(xmin = e_min)
        ax.set_xlim(xmax = e_max)
        
        if flux_min != None : 
            flux_min = flux_min*u.Unit(flux_unit)*u.Unit(e_unit)**index
            ax.set_ylim(ymin = flux_min.value)
        if flux_max != None : 
            flux_max = flux_max*u.Unit(flux_unit)*u.Unit(e_unit)**index
            ax.set_ylim(ymax = flux_max.value)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid("both")
    ax.legend()

    return (flux_min,flux_max) # expect ymin, ymax
#------------------------------------------------------------------------------
def flux_versus_time(dsets, ax=None,
                          emin=None, emax=None,
                          tmin=None, tmax=None,
                          stacked = False,
                          style="bar", fit=False,
                          e_unit = "GeV",
                          flux_min = 1.e-23,
                          xscale="log",
                          color="tab:blue",debug=False):

    t0 = dsets[0].gti.time_start[0]
    ftheory, time, errtime, dnde, errn, errp, erange = [], [], [], [], [], [], []

    # Loop over time slices, get the flux for the given energy range
    for i,ds in enumerate(dsets):

        time.append((ds.gti.time_start[0]-t0).sec*u.s + ds.gti.time_sum/2)
        errtime.append(0.5*ds.gti.time_sum)


        # If Energy boundaries at not given, use the range from the safe mask
        # If energy range is outside the dataset energy range, returns
        e_min = max(emin,ds.energy_range[0]) if emin!=None else ds.energy_range[0]
        e_max = min(emax,ds.energy_range[1]) if emax!=None else ds.energy_range[1]

        if e_min < ds.energy_range[0] or e_max>ds.energy_range[1] or e_max<=e_min:
            if debug: print(" E boundaries out of range")
            return False
#         if (debug):
#             print("--- ",ds.name," ---")
#             print("Eff. E range = {:5.2f} {:5.2f}".format(e_min,e_max))

        # Extracted flux - it will be limited to existing ereco edges
        # Change emin and emax to reflect this
        fex  = extract_flux_points(ds,[e_min,e_max],debug=debug)

        e_min = fex.table["e_min"].quantity[0].to(e_unit)
        e_max = fex.table["e_max"].quantity[0].to(e_unit)
        e_ref = fex.table["e_ref"].quantity[0].to(e_unit)
        flx  = fex.table["dnde"].quantity[0]
        ep   = fex.table["dnde_errp"].quantity[0]
        en   = fex.table["dnde_errn"].quantity[0]

        # Mean theoretical flux at reference energy in the bin
        if stacked :
            fth = flx
        else:
            fth  = ds.models[0].spectral_model(e_ref)

        if (debug):
            print("#{:2} {:6.1f} - {:6.1f} : F= {:5.1e} +{:5.1e} -{:5.1e} {:s}"
                  .format(i,e_min.value,e_max, flx.value, ep.value, en.value, flx.unit))
            print("   {:6s} - {:6s}   T= {:5.2e} T/F= {}"
                  .format("","",fth, fth.value/flx.to(fth.unit).value))

        dnde.append(flx)
        errp.append(ep)
        errn.append(en)
        ftheory.append(fth)
        erange.append([e_min,e_max])

    time    = np.asarray([x.value for x in time])*time[0].unit
    errtime = np.asarray([x.value for x in errtime])*errtime[0].unit
    dnde    = np.asarray([x.item().value for x in dnde])*dnde[0].unit
    errn    = np.asarray([x.item().value for x in errn])*errn[0].unit
    errp    = np.asarray([x.item().value for x in errp])*errp[0].unit
    ftheory = np.asarray([f.value for f in ftheory])*ftheory[0].unit
    erange  = np.asarray( [[e[:][0].value,e[:][1].value] for e in erange]*u.Unit(e_unit))

    # Fit a t**-beta dependence if required
    if fit:
        t_min = tmin if (tmin!=None) else min(time)
        t_max = tmax if (tmax!=None) else max(time)
        mask = (time>t_min) & (time<t_max) & np.isfinite(dnde) # Remove undefined dnde values

        from scipy.optimize import curve_fit
        def func(x, a, b):
            return a*x**-b   #+c
        popt, pcov = curve_fit(func, time[mask], dnde[mask])
        a  = popt[0]
        da = np.sqrt(pcov[0][0])
        b  = popt[1]
        db = np.sqrt(pcov[1][1])
#         c = popt[2]
#         dc = np.sqrt(pcov[2][2])
        print(" a=",a,"+/-",da)
        print(" b=",b,"+/-",db)
#         print(" c=",c,"+:-",dc)

    # The plots
    ax = plt.gca() if ax is None else ax

    with quantity_support():
        label= "{:5.1f}-{:5.1f} {}".format(e_min.value,e_max.to(e_min.unit).value,str(e_min.unit))
        ax.errorbar(x = time, y = dnde, xerr = errtime, yerr = [errn,errp],
                    color = color, ls="", label=label)

        axx = ax.twinx()
        eb = axx.errorbar(x = time ,y = [e[:][0] for e in erange], xerr = errtime, yerr=0,
                     color="grey",ls="",alpha=0.5,label="$E_{min}, E_{max}$")
        eb[-1][0].set_linestyle('dashdot')
        eb = axx.errorbar(x = time ,y = [e[:][1] for e in erange], xerr = errtime, yerr=0,
                      color="grey",ls="",alpha=0.5)
        eb[-1][0].set_linestyle('dashdot')
        axx.set_ylabel("Range ("+e_unit+")")
        axx.set_yscale("log")
        axx.legend()

        if fit:
            label = r"$\beta$= "+str(round(b,2))+"$\pm$"+str(round(db,2))
            ax.plot(time[mask], func(time[mask].value,a,b),color="tab:orange",alpha=0.9,label=label)
            ax.axvline(x=t_min,ls="--",color="brown",label="Fit limits")
            ax.axvline(x=t_max,ls="--",color="brown")

        if not stacked:
            if style=="bar":
                ax.bar(time,ftheory,width=2*errtime,
                       alpha=0.2,color=color,label="Model")
            elif style=="line":
                ax.plot(time,ftheory, alpha=0.5,color=color,ls=":",marker="o",label="Model")

        if (flux_min != None): ax.set_ylim(ymin=flux_min*dnde[0].unit)
        if time[-1] > 1*u.d:  ax.axvline(x=1*u.d,ls=":",
                                          color="grey",label="One day")

    ax.set_xlabel("Elapsed time ("+ax.xaxis.get_label_text()+")")
    ax.set_yscale("log")
    ax.set_xscale(xscale)
    ax.legend()

    return

#------------------------------------------------------------------------------
def fluxes_versus_time(dsets,
                       tmin=None, tmax=None,
                       xscale="linear",yscale="log",
                       xsize=14,ysize=3,
                       stacked = False,
                       debug=False):

    e_edges = dsets[0].excess.geom.axes[0].edges

    for i,e in enumerate(e_edges[:-1]):
        fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(xsize,ysize))
        status = flux_versus_time(dsets,emin = e_edges[i], emax = e_edges[i+1],
                         style="line",
                         xscale = xscale,
                         ax=ax, #color=color,
                         flux_min=None,
                         stacked = stacked,
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

    return

###############################################################################
### A posteriori simulation
###############################################################################
def lightcurve(dsets,
               t_unit = u.s, tmin=0*u.s, tmax=1000*u.s, binwidth=10*u.s,
               tag="excess",
               ax=None, xscale = "linear", yscale="linear",
               style="bar", marker="o",
               alpha=1, color="black", fillcolor=None,
               debug=False, **kwargs):

    from scipy.interpolate import interp1d

    tmin     = tmin.to(t_unit)
    tmax     = tmax.to(t_unit)
    binwidth = binwidth.to(t_unit)

    # Create the time interpolation of the counting rate
    tsamp, dndt, tslice = [], [], []
    t0 = dsets[0].gti.time_start[0]
    duration = 0*t_unit

    for i,ds in enumerate(dsets):

        livetime = ds.gti.time_sum

        tsamp.append( ((ds.gti.time_start[0]-t0).sec*u.s
                           + livetime/2).to(t_unit) )
        tslice.append(duration)
        duration = duration + livetime.to(t_unit)

        if (tag == "excess"):
            dndt.append(ds.excess.data[ds.mask_safe].sum()/livetime.to(t_unit))
        elif (tag=="prediction"):
            if gammapy.__version__ == "0.17":
                dndt.append(ds.npred_sig().data[ds.mask_safe].sum()/livetime.to(t_unit))
            else: #0.18.3
                dndt.append(ds.npred_signal().data[ds.mask_safe].sum()/livetime.to(t_unit))
        else:
            print("light_curve: tag not implemented")

    # From now on, all variables have values corresponding to t_unit
    tsamp  = np.asarray([x.value for x in tsamp])
    dndt   = np.asarray([x.value for x in dndt])
    tslice = np.asarray([x.value for x in tslice])

    tmin = tmin.value
    tmax = tmax.value
    binwidth = binwidth.value

    # Warning : f is not a quantity - Therefore all time units should be
    # the same before the interpolation - f should have unit 1/t_unit like dndt
    f = interp1d(tsamp,dndt,fill_value="extrapolate")
    #print(tsamp,dndt,f(tsamp))

    # Generate a number of counts for each bin
    nbin    = int((tmax-tmin)/binwidth)
    t_edges = np.linspace(tmin,tmax,nbin+1)
    t_bins  = t_edges[:-1] + 0.5*(t_edges[1:]-t_edges[:-1])
    # print(t_bins)

    # Compute count number in each bin = rate * bin width"
    n_random = []
    for t in t_bins:
        counts = f(t)*binwidth
        # print(t,"rate=",f(t),"counts=",counts)
        if (tag != "prediction"):
            # If excess is negative, do not fluxtuate
            if (counts>=0):
                n_random.append(np.random.poisson(int(counts)))
            else:
                print(" Negative counts : ",counts," not fluctuated")
                n_random.append(int(counts))
        else:
            n_random.append(counts)
    #print(n_random)

    ax = plt.gca() if ax is None else ax
    
    if (style=="bar"):
        ax.bar(t_bins, n_random, width=binwidth, edgecolor=color, color=fillcolor, alpha=alpha,label=tag, **kwargs)
    if (style=="line"):
        ax.plot(t_bins, n_random, color=color,marker=marker,lw=1, alpha=alpha,label=tag, **kwargs)
    if (debug):
        ax.plot(t_bins,f(t_bins)*binwidth, ls=":",color="red")
        ax.plot(tsamp[tsamp<tmax],dndt*binwidth,marker="o",ls="",color="green",label="Initial data")

    ax.set_xlabel("Observation duration ("+str(t_unit)+")")
    ax.set_ylabel("Counts per bin (" + str(binwidth) + str(t_unit)+")")
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)

    for xl in tslice[tslice<tmax]:
        ax.axvline(x=xl,ls="--",lw=1,alpha=0.5, color="grey")
    ax.axhline(y=0)

    ax.legend()

    return ax
#------------------------------------------------------------------------------
def stacked_versus_time(dstacked, dsets=None,
                              rate=True, ax=None, debug=False,
                              xscale="linear", yscale="linear",
                              color="tab:blue"):

    """
    Display stacked statictics and model statistics if original datasets is 
    given (since model stacking is irrelevant)

    """
    print("to be finalised !")
    return

    if ax == None : fig,ax = plt.subplots(figsize=(7,7))

    duration = 0*u.s
    npred = 0
    with quantity_support():
        t0 = dstacked[0].gti.time_start[0]
        i = 0
        for dst, ds in zip(dstacked,dsets):

            livetime= dst.gti.time_delta[i]
            time       = (dst.gti.time_start[i]-t0).sec*u.s + livetime/2

            print(i," >>> ",livetime)
            errtime    = 0.5*(livetime-duration)
            i+=1

            xs     = dst.excess.data.sum()
            if gammapy == "0.17":
                npred += ds.npred_sig().data[ds.mask_safe].sum()
            else: # 0.18.2
                npred += ds.npred_signal().data[ds.mask_safe].sum()

            norm   = 1*u.dimensionless_unscaled if rate is False else 1/2/errtime

            if (debug):
                print(" At {:8.2f} +/-{:8.2f} : nxs = {:8.2f} +/-{:8.2f} Ratio={:5.2f}"
                      .format(time.value,errtime,xs*norm.value,np.sqrt(xs)*norm,npred/xs))

            # Handling of quantities by Errorbar is strange - requires numpy arrays - quantities are not scalars !

            ax.errorbar(x = [time.value]*time.unit,
                        y = [xs*norm.value]*norm.unit,
                        xerr = [errtime.value]*errtime.unit,
                        yerr = [np.sqrt(xs)*norm.value]*norm.unit,
                        label="Excess counts",color=color)
            ax.bar(time,npred*norm,width=2*errtime,alpha=0.2,color=color,label="Theory")
            duration += livetime

        ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
        ax.set_ylabel("Counts (" + ax.yaxis.get_label_text()+ ")")
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        import collections
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = collections.OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    return

##############################################################################
###  UTILITIES
##############################################################################
def panels(dsets,func = None,
           nmaxcol  = 4, xsize = 5, ysize = 5,
           xscale   = "log", yscale ="log",
           max_margin = 1.1,
           fixrange = True,
           tref     = Time('2000-01-01T00:00:00',format='isot', scale='utc'),
            **kwargs):

    """
    Display the plots obtained from a function applied to each individual
    dataset in an optimised panel.

    Parameters
    ----------
    dsets : Datasets object
        A dataset collection.
    func : Function, optional
        The pointer ot a function with first argument a dataset object, and
        a matplotib axis. The rest of the arguments are passed through the
        keyword arguments (**kwargs). The default is None.
    nmaxcol : Integer, optional
        The maximum number of columns in the panel. The default is 5.
    xsize : Float, optional
        The wifth of the columns in the panel. The default is 5.
    ysize : Float, optional
        The height of the rows in the panel. The default is 7.
    fixrange : Boolean, optional
        If True, all plot y axis are identical. The default is False.
    **kwargs : Keywor arguments
        Extra arguments for the function (func)

    Returns
    -------
    None.

    """

    # Panel geometry
    nplots = len(dsets)
    ncols  = min(nmaxcol,nplots) # If nplots < nmaxcol, take nplots
    nrows  = int(nplots/ncols)+ 1*(nplots%ncols != 0)
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows,
                           figsize=(xsize*ncols,ysize*nrows),
                           sharex=True, sharey=fixrange)
    # Plot common y limit - if fixrange os True
    ymax = -np.Inf
    ymin =  np.Inf

    with quantity_support():
        iplot = 0
        import itertools

        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            if nplots !=1:
                ax0 = ax[jrow][icol] if (nrows>1) else ax[icol]
            else:
                ax0=ax

            # If nothing to plot, blank space instead of empty plot
            if (iplot >= nplots):
                ax0.axis('off')
                continue # Next plot

            # Plot function here
            t_since_trigger = dsets[iplot].gti.time_start[0] - tref
            (y1, y2) = func(dsets[iplot], ax=ax0, elapsed=t_since_trigger, **kwargs)
            if y1 != y2:
                if (y2 > ymax): ymax = y2
                if (y1 < ymin): ymin = y1

            # Compactify
            if (jrow+1 != nrows): ax0.set_xlabel(None)
            if (icol !=0): ax0.set_ylabel(None)
            ax0.tick_params(which='major', length=10, width=2, direction='in')
            ax0.tick_params(which='minor', length=5, width=2, direction='in')

            iplot+=1

            ax0.set_xscale(xscale)
            ax0.set_yscale(yscale)
            
        if (fixrange):
            plt.subplots_adjust(right=0.8,hspace=0.05,wspace=0.05)
            for ax in fig.get_axes():
                ymax = max_margin*ymax # If scale is log, it is applied to the log value
                if type(ymax) == astropy.units.quantity.Quantity:
                    ymax = ymax.value
                if type(ymin) == astropy.units.quantity.Quantity:
                     ymin = ymin.value                   
                ax.set_ylim(bottom=ymin, top=ymax)
        else:
            plt.subplots_adjust(right=0.8,hspace=0.05)
    
        # use first label handles as the last plot can be empty
    #     h0, lbl0 = fig.axes[0].get_legend_handles_labels()
    #     fig.legend(h0, lbl0, loc="center right", borderaxespad=0.1)
        if (fixrange): fig.tight_layout(h_pad=0)
        else: fig.tight_layout(h_pad=0, w_pad=0)

    return fig

#------------------------------------------------------------------------------
def counts_from_plot(ax, dset, plot_line=False):
    """
    THis a utility to get the counts form the plot and check the count_str
    displays the content of the dataset

    Parameters
    ----------
    ax : Matplotlib axis
        The Gammapy dataset count plot to be analysed
    dest: Dataset object
        The original Dataset object to be compared with
    plot_line : Boolean, optional
        If True, plot vertical and horizontal lines to the extracted points.
        The default is False.

    Returns
    -------
    None.

    """

    Emu = ax.lines[0].get_xdata()
    ymu = ax.lines[0].get_ydata()
    #Edata = ax.lines[1].get_xdata()
    ydata = ax.lines[1].get_ydata()
    i=0
    print("--- Counts from plot ", 35*"-")
    print("{:3s} - {:>8s}     {:>8s} {:>8s}".format("bin","E"," Model", "Excess"))
    for E,y,yd in zip(Emu,ymu,ydata):
        print("{:3d} - {:8.2f} {:8.2f} {:8.2f}".format(i, E,y,yd))
        if (plot_line):
            ax.axhline(y,ls=":",color="grey")
            ax.axhline(yd,ls=":",color="red")
        i+=1

    print(" Total  excess counts = ",ydata.sum())
    print(" Masked excess counts = ",ydata[3:-1].sum(),
          "from dset ",dset.excess.data[dset.mask_safe].sum())
    print(" Total counts         = ",ymu.sum())
    print(" Masked counts        = ",ymu[3:-1].sum())

    return

#------------------------------------------------------------------------------
def plot3D(x, y, data,
                 tag    = "Unknown",
                 ax3d   = None,
                 zscale = "linear",
                 cmap   = "plasma",
                 debug=False):
    """
    Create a 3D plot with x is energy, y is time, z is data

    Parameters
    ----------
    x : 1D Numpy array or list
        x values
    y : 1D Numpy array or list
        y values
    data : 2D Numpy array
        The data
    tag : String, optional
        Used as a title. The default is "Unknown".
    ax3d : Matplotlib 3D axis, optional
        Axis used for plotting. The default is None.
    zscale : String, optional
        z axis scale, either "liner" or "log". The default is "linear".
    cmap : String, optional
        A colormap name. The default is "plasma".
    debug : Boolean, optional
        If True print some information. The default is False.

    Returns
    -------
    None.

    """

    if (ax3d == None):
        fig = plt.figure()
        ax3d = fig.add_subplot(111, projection='3d')

    yE, xt = np.meshgrid(y,x)

    zlbl = "Counts"
    if (zscale=="log"):
        data = np.log10(data)
        zlbl = zlbl + "(log)"

    ax3d.plot_surface(xt, yE, data,cmap=cmap, alpha=1.0,edgecolor='none')

    ax3d.set_title(tag)
    ax3d.set_ylabel('log(E) ('+str(y.unit)+')')
    ax3d.yaxis.labelpad=20
    ax3d.set_xlabel('Elapsed time ('+str(x.unit)+')')
    ax3d.xaxis.labelpad=20
    ax3d.set_zlabel(zlbl)
    ax3d.zaxis.labelpad=20
    ax3d.view_init(30, 45)

    return
