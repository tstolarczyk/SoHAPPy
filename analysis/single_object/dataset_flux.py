# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

@author: Stolar
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

from utilities import t_str, single_legend

from gammapy.modeling.models import PowerLawSpectralModel
from gammapy.modeling.models import SkyModel
from gammapy.estimators import FluxPointsEstimator

__all__ = ['extract_spectrum','flux_versus_time', 'fluxes_versus_time']

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
    """
    Plot flux versus time

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
        DESCRIPTION.

    """

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
        if stacked : fth = flx
        else: fth  = ds.models[0].spectral_model(e_ref)

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
    ax.legend(ncol=2)

    return

#------------------------------------------------------------------------------
def fluxes_versus_time(dsets,
                       tmin=None, tmax=None,
                       xscale="linear",yscale="log",
                       xsize=14,ysize=3,
                       stacked = False,
                       debug=False):
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
    yscale : String, optional
        "log" or "linear". The default is "log".
    xsize : float, optional
        Figure width. The default is 14.
    ysize : float, optional
        Figure height. The default is 3.
    stacked : Boolean, optional
        If True, the datasets have been stacked. The default is False.
    debug : Boolean, optional
        IF True, let's talk a bit. The default is False.

    Returns
    -------
    None.

    """

    e_edges = dsets[0].excess.geom.axes[0].edges

    for i,e in enumerate(e_edges[:-1]):
        fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(xsize,ysize))
        status = flux_versus_time(dsets,
                                  emin = e_edges[i], emax = e_edges[i+1],
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

#------------------------------------------------------------------------------
def models_plot(dsets):
    """
    Plot the models of the individual datasets in the list

    Parameters
    ----------
    dsets : List of Dataset
        Current dataset list.

    Returns
    -------
    None.

    """

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
    
    return
#------------------------------------------------------------------------------
def extract_flux_points(ds, index=2, sigm_ul=2, debug=False):
    """
    Extract flux points.
    The flux is extracted from the execess number assuming a certain model
    The difference in results between two models is larger when the energy 
    binning is larger.
    
    Parameters
    ----------
    ds : dataset list.
        DESCRIPTION.
    index : float, optional
        Index of the power law spectrum used in the fit. The default is 2.
    sigm_ul : float, optional
        Significance below which an upper limit is computed. The default is 2.
    debug : Boolean, optional
        If True let's talk a bit. The default is False.

    Returns
    -------
    fex : FluxPointsEstimator
        Extracted flux points.

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
    fpe = FluxPointsEstimator(
            energy_edges = ds.background.geom.axes[0].edges, 
            norm_min = 0.2, norm_max  =5, 
            norm_n_values = 11, norm_values = None,
            n_sigma=1, n_sigma_ul=2, 
            reoptimize=False )
    fex = fpe.run(datasets = ds)
    
    # Replace with a limit of significance lower than 2
    fex.table["is_ul"] = fex.table["ts"] < sigm_ul**2

    # Put the original model back in place
    ds.models = model_init

    if (debug):
        print(" --- ",ds.name," --- assuming a powerlaw with index = 2")
        for i in range(len(ds.background.axis[0].edges)-1):
            print("{:7.2f} - {:7.2f} : {:5.2e} -{:5.2e} +{:5.2e}"
                  .format(fex.table["e_min"].quantity[i],
                          fex.table["e_max"].quantity[i],
                          fex.table["dnde"].quantity[i].value,
                          fex.table["dnde_errn"].quantity[i].value,
                          fex.table["dnde_errp"].quantity[i]))

    return fex

#------------------------------------------------------------------------------
def fit_model(fit_tag, amplitude, e_ref):
    """
    Define the model for the fit.

    Parameters
    ----------
    fit_tag : String
        A tag for choosing the model.
    amplitude : astropy.Quantity
        The amplitude of the model flux.
    e_ref : astropy.Quantity
        The reference energy of the model flux.

    Returns
    -------
    model_fit : gammapy SkyModel
        The model to be fitted.

    """
    
    model_fit = None
    
    # Create the appropriate model
    if fit_tag=="cutoff":
        from gammapy.modeling.models import ExpCutoffPowerLawSpectralModel
        lambda_ = 1 * u.Unit("TeV-1")
        model_fit = ExpCutoffPowerLawSpectralModel( 
                    index     = 2.0,
                    amplitude = amplitude*np.exp(lambda_*e_ref), # Correct phi0 for cutoff
                    lambda_   = lambda_,
                    reference = e_ref)
    
    elif fit_tag=="powerlaw":
        from gammapy.modeling.models import PowerLawSpectralModel
        model_fit = PowerLawSpectralModel(index     =   2.0,
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

    return model_fit

#------------------------------------------------------------------------------
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
        res_index   = fit.confidence("index") # get index error
        res_lambda_ = fit.confidence("lambda_")
        
        text = "{}:\n$\gamma = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
            .format(mtag[1],result.parameters["index"].value, \
            res_index["errp"],res_index["errn"])
        text = text+"\n $\lambda = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                .format(result.parameters["lambda_"].value, \
                    res_lambda_["errp"],res_lambda_["errn"])

    elif mtag[0] == "PowerLawSpectralModel":
        res_index   = fit.confidence("index") # get index error
        
        text = "{}:\n$\gamma = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
            .format(mtag[1],result.parameters["index"].value, \
            res_index["errp"],res_index["errn"])

    elif mtag[0] == "LogParabolaSpectralModel":
        res_alpha   = fit.confidence("alpha") # get index error
        res_beta    = fit.confidence("beta")
        
        text = "{}:".format(mtag[1]) +"\n"+r"$\alpha = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
            .format(result.parameters["alpha"].value, \
            res_alpha["errp"],res_alpha["errn"])
        text = text+"\n"+r"$\beta = {:3.1f}  ^{{+{:3.2f}}} _{{-{:3.2f}}}$" \
                .format(result.parameters["beta"].value, \
                    res_beta["errp"],res_beta["errn"]) 
    else:
        text=""
                    
    return text
#------------------------------------------------------------------------------
def fit_spectrum(ds, fex, e_ref=1*u.TeV, fit_tag=None, debug=False):
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
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

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
    amplitude = (fex.table["ref_dnde"].quantity)[idx] * fex.table["norm"].quantity[idx]
    # if (debug):
    #     print(" Best energy : ",e_center[idx])
    #     print(" Amplitude   : ",amplitude)
    #     ax.scatter(e_center[idx],amplitude*e_center[idx]**2,marker="x",color="red")
    model_fit = fit_model(fit_tag , amplitude,  e_ref)

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
    
    return fit, result
#------------------------------------------------------------------------------
def extract_spectrum(ds,
                     elapsed = 0*u.s,
                     index = 2, flux_unit = "TeV-1 cm-2 s-1", e_unit = "TeV",
                     ax = None, 
                     model_style ="bar", color = "tab:blue", lw=2,
                     alpha_model = 0.2, alpha_fit = 0.5, fit_color = "blue",
                     xscale   = "log", yscale   = "log",
                     theory   = True,     stacked  = False, 
                     fit_tag  = "cutoff", tag      = None,
                     flux_min = None,     flux_max = None,
                     e_min    = 10*u.GeV, e_max    = 20*u.TeV, 
                     e_ref    = 1*u.TeV,
                     count_min=10,
                     debug = False):
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
    print(" ------------------------ ",ds.name," -----------------------------")

    # If not enough count, do not attempt to extract spectrum
    count_max = max(ds.excess.data.flatten())
    if count_max <= count_min:
        ax.text(0.05,0.95,
                "Counts too low ("+str(round(count_max,1))+")",
                color = color,
                transform=ax.transAxes)
        print(" Counts too low")
        return (0, 0)

    # Define label
    if tag == None or len(tag)==0 :
        label = "$t$ = " + t_str(elapsed) \
          + "\n$\Delta$t = "+t_str(ds.gti.time_sum)
    else: 
        label = "$t_{"+tag+"}$ = " + t_str(elapsed) +"\n" \
          + "$\Delta t_{"+tag+"}$ = "+t_str(ds.gti.time_sum)
    
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
    fex = extract_flux_points(ds,debug=debug)
    fex.plot(energy_power = index,
             flux_unit    = flux_unit,
             energy_unit  = e_unit,
             ax = ax, ls = ls, label = label, lw=lw, color=color)

    ###---------------------
    # Plot theory (at reconstructed energies !)
    ###---------------------
    if theory and not stacked: # Stacked dataset have no proper spectrum
        model   = ds.models[0].spectral_model
        ftheory = model(e_center).to(flux_unit)*e_center.to(e_unit)**index
        width_r = e_edges[1:] -e_center
        width_l = e_center - e_edges[:-1]

        with quantity_support():
            if model_style == "bar":
                ax.bar(e_center.to(e_unit), ftheory, 
                       width = -width_l, align="edge",
                       alpha=alpha_model, color = color, 
                       label=ds.name+" model")
                ax.bar(e_center.to(e_unit), ftheory, 
                       width = width_r, align="edge",
                       alpha=alpha_model, color = color)
            if model_style == "scatter":
                ax.scatter(e_center.to(e_unit), ftheory, marker="o",
                           facecolor="white", edgecolors="red",
                           label=ds.name+" model")

    ###---------------------
    # Fit the data and plot
    ###---------------------
    # MINUIT reference : https://iminuit.readthedocs.io/en/latest/reference.html

    if (fit_tag != None):
        # Copy the initial model to reset the dataset in the end
        model_init = ds.models       
        fit, result = fit_spectrum(ds, fex, e_ref=e_ref, fit_tag=fit_tag)                

    ###---------------------
    # Plot the fit and the result text
    ###---------------------    
        # Result of fit text
        
        # If fit failed, just mention it on the plot
        if not result.success or result.parameters["amplitude"].value<0:
            ax.text(0.05,0.95,"FIT FAILED ",
                    color=color,
                    transform=ax.transAxes)
            print("Fit failed!")
        else:
            fit_text = fit_result(fit, result)  
            print(fit_text)
            if (debug):
                print(" FITTING : ",fit_tag)
                print("*** ",result)
                print("*** ",result.parameters)
                print(fit_text)
                print(ds.models[0].spectral_model)

        # model_fit.plot_error(ax=ax,
        ds.models[0].spectral_model.plot_error(ax=ax,
                             energy_range = ds.energy_range,
                             energy_power = index,
                             energy_unit  = e_unit,
                             alpha = 0.2, 
                             color = fit_color)
    
        # Note that ds.models[0].spectral_model.plot also works
        # model_fit.plot(ax=ax,energy_range     = ds.energy_range,
        if tag != None and len(tag) != 0:
            mlabel = "$"+tag+"$"+" "+fit_tag
        else:
            mlabel = fit_tag
        ds.models[0].spectral_model.plot(ax=ax,
                                         energy_range = ds.energy_range,
                                         energy_power = index,
                                         energy_unit  = e_unit,
                                         flux_unit    = flux_unit,
                                         label        = mlabel,
                                         alpha = alpha_fit, 
                                         color = fit_color)   
        ds.models = model_init
        
    # Decoration
    with quantity_support():
        if (not stacked):  # Plot energy range
            ax.axvline(ds.energy_range[0],ls=":",color="grey")
            ax.axvline(ds.energy_range[1],ls=":",color="grey")

    # This does not work with units within quantity_support
    ax.set_xlim(xmin = e_min.to(e_unit).value)
    ax.set_xlim(xmax = e_max.to(e_unit).value)
    if flux_min != None : 
        flux_min = flux_min*u.Unit(flux_unit)*u.Unit(e_unit)**index
        ax.set_ylim(ymin = flux_min.value)
    if flux_max != None : 
        flux_max = flux_max*u.Unit(flux_unit)*u.Unit(e_unit)**index
        ax.set_ylim(ymax = flux_max.value)
    
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid("both")
    single_legend(ax, fontsize=12)

    return (flux_min,flux_max) # expect ymin, ymax