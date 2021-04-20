# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 16:49:31 2020

@author: Stolar
"""
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.visualization import quantity_support
from gammapy.estimators import FluxPointsEstimator

__all__ = ['extracted_spectrum','excess_counts','residuals',
           'excess_versus_time','stacked_versus_time','flux_versus_time',
           'fluxes_versus_time','excess_versus_E_and_time','lightcurve',
           'panels']
#------------------------------------------------------------------------------
def extracted_spectrum(ds, index=0,
                       flux_unit = "GeV-1 cm-2 s-1", e_unit = "GeV",
                       style ="bar", color ="tab:blue", ax=None,
                       xscale = "linear", yscale = "log",
                       minflux = 10**-18*u.Unit("GeV-1 cm-2 s-1"),
                       theory  =False, debug = False):
    """
    Note that flux is extracted on the masked datsets
    """
    ax = plt.gca() if ax is None else ax

    axis_reco = ds.excess.geom.axes[0]

    # Plot extracted flux from built-in function - Includes the errorbar
    fex = extract_flux_points(ds, axis_reco.edges)

    fex.plot(energy_power = index,
             flux_unit    = flux_unit,
             energy_unit  = e_unit,
             ax = ax,
             label = "$\Delta$t="+tobs_plot(ds)+"\nExtracted")

    # Plot theory
    if (theory == False): # Stacked dataset have no proper spectrum
        e_reco  = axis_reco.center
        de_reco = axis_reco.bin_width
        model = ds.models[0].spectral_model
        ftheory = model(e_reco).to(flux_unit)*e_reco.to(e_unit)**index

        if style== "bar":
            ax.bar(e_reco.to(e_unit), ftheory, width = de_reco.to(e_unit),
                   alpha=0.2, color = color, label=ds.name+" model")

        if style == "scatter":
            ax.scatter(e_reco.to(e_unit), ftheory, marker="o",
                       facecolor="white", edgecolors="red",
                       label=ds.name+" model")

    if (index != 0): # Rename ambiguous label
        ax.set_ylabel("$E^"+str(index)+"*$"+ax.yaxis.get_label_text())

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_ylim(ymin=minflux.value)
    ax.legend()

    return

#------------------------------------------------------------------------------
def extract_flux_points(ds,e_edges, min_siginificance=4,debug=False):

    fpe = FluxPointsEstimator(e_edges=e_edges)
    flux_points = fpe.run(datasets=ds)
    flux_points.table["is_ul"] = flux_points.table["ts"] < min_siginificance

    if (debug):
        print(flux_points.table["dnde","dnde_errn","dnde_errp"])

    return flux_points

#------------------------------------------------------------------------------
def excess_counts(ds,ax=None, stacked=False, mincounts = 0.1, **kwargs):

    """
    This is a modified version of dataset plot_counts method
    """
    ax = plt.gca() if ax is None else ax

    ds.excess.plot(ax=ax, label="$\Delta$t="+tobs_plot(ds)+"\nExcess")
    if (not stacked): # Weeir predcietd counts when stacked
        ds.npred_sig().plot(ax=ax, label="mu_src")
        # print(">>>>> ",ds.npred_sig().data)
    ds._plot_energy_range(ax=ax)

    # If stacked, force plotting beyongthe energy range
    if (stacked):
        eaxmin, eaxmax = ax.get_xlim() # Present plot limits
        emin = 0.95*ds.excess.geom.axes[0].edges[0]
        emax = 1.05*ds.excess.geom.axes[0].edges[-1]

        ax.set_xlim(min(emin.value,eaxmin),max(emax.value,eaxmax))
    ax.legend(numpoints=1)
    ax.set_ylim(ymin=mincounts)
    ax.set_title("")
    return
#------------------------------------------------------------------------------
def residuals(ds,ax=None,**kwargs):
    ds.plot_residuals(ax=ax)
    return

#------------------------------------------------------------------------------
def excess_versus_time(dsets, rate=True, unmasked=False, ax=None,
                       xscale="linear", yscale="linear",
                       color1 = "tab:blue",color2 = "tab:orange",
                       debug  = False):

    ax = plt.gca() if ax is None else ax

    with quantity_support():
        t0 = dsets[0].gti.time_start[0]
        for i,ds in enumerate(dsets):
            time       = (ds.gti.time_start[0]-t0).sec*u.s + ds.livetime/2
            errtime    = 0.5*ds.livetime

            xs        = ds.excess.data.sum()
            xs_msk    = ds.excess.data[ds.mask_safe].sum()
            npred     = ds.npred_sig().data.sum()
            npred_msk = ds.npred_sig().data[ds.mask_safe].sum()

            if (rate==True):
                norm = 1/2/errtime
            else:
                norm=1*u.dimensionless_unscaled
            if (debug):
                print(time, errtime)
                print(xs*norm,np.sqrt(xs)*norm)

            # Handling of quantities by Errorbar is strange - requires numpy arrays - quantities are not scalars !

            ax.errorbar(x = [time.value]*time.unit,
                        y = [xs_msk*norm.value]*norm.unit,
                        xerr = [errtime.value]*errtime.unit,
                        yerr = [np.sqrt(xs_msk)*norm.value]*norm.unit,
                        label="Excess counts",color=color1)
            ax.bar(time,npred_msk*norm,width=2*errtime,alpha=0.2,color=color1)

            # Before masking
            if (unmasked):
                ax.errorbar(x    = [time.value]*time.unit,
                            y    = [xs*norm.value]*norm.unit,
                            xerr = [errtime.value]*errtime.unit,
                            yerr = [np.sqrt(xs)*norm.value]*norm.unit,
                            label="Total excess counts",color=color2)
                ax.bar(time,npred*norm,width=2*errtime,alpha=0.2,color=color2)


        ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
        ax.set_ylabel("Counts (" + ax.yaxis.get_label_text()+ ")")

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        import collections
        handles, labels = ax.get_legend_handles_labels()
        by_label = collections.OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    return

#------------------------------------------------------------------------------
def stacked_versus_time(dstacked, dsets=None,
                              rate=True, ax=None, debug=False,
                              xscale="linear", yscale="linear",
                              color="tab:blue"):

    """
    Display stacked statictics and model statistics if original datasets is given (sonce model stacking is irrelevant)

    """
    if (ax==None): fig,ax = plt.subplots(figsize=(7,7))

    duration = 0*u.s
    npred = 0
    with quantity_support():
        t0 = dstacked[0].gti.time_start[0]
        for dst, ds in zip(dstacked,dsets):
            time       = (dst.gti.time_start[0]-t0).sec*u.s + duration + (dst.livetime - duration)/2
            errtime    = 0.5*(dst.livetime-duration)

            xs     = dst.excess.data.sum()
            npred += ds.npred_sig().data[ds.mask_safe].sum()
            norm   = 1*u.dimensionless_unscaled if rate is False else 1/2/errtime

            if (debug):
                print(time, errtime)
                print(xs*norm,np.sqrt(xs)*norm)
                print("npred/xs",npred/xs)

            # Handling of quantities by Errorbar is strange - requires numpy arrays - quantities are not scalars !

            ax.errorbar(x = [time.value]*time.unit,
                        y = [xs*norm.value]*norm.unit,
                        xerr = [errtime.value]*errtime.unit,
                        yerr = [np.sqrt(xs)*norm.value]*norm.unit,
                        label="Excess counts",color=color)
            ax.bar(time,npred*norm,width=2*errtime,alpha=0.2,color=color,label="Theory")
            duration = dst.livetime

        ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
        ax.set_ylabel("Counts (" + ax.yaxis.get_label_text()+ ")")
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        import collections
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = collections.OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def flux_versus_time(dsets, ax=None, emin=None, emax=None,
                     style="bar",
                     color="tab:blue",debug=False):

    t0 = dsets[0].gti.time_start[0]
    if (emin==None): emin = dsets[0].excess.geom.axes[0].edges[0]
    if (emax==None): emax = dsets[0].excess.geom.axes[0].edges[-1]
    emin = max(emin,dsets[0].excess.geom.axes[0].edges[0])
    emax = min(emax,dsets[0].excess.geom.axes[0].edges[-1])

    ftheory, time, errtime, dnde, errn, errp  = [], [], [], [], [], []

    # Loop over time slices, get the flux for the give energy range
    for i,ds in enumerate(dsets):

        # Date
        time.append((ds.gti.time_start[0]-t0).sec*u.s + ds.livetime/2)
        errtime.append(0.5*ds.livetime)

        # Extracted flux
        fex  = extract_flux_points(ds,[emin,emax])
        dnde.append(fex.table["dnde"].quantity)
        errp.append(fex.table["dnde_errp"].quantity)
        errn.append(fex.table["dnde_errn"].quantity)

        # Mean theoretical flux
        ftheory.append(ds.models[0].spectral_model.integral(emin,emax)/(emax-emin))

    time    = np.asarray([x.value for x in time])*time[0].unit
    errtime = np.asarray([x.value for x in errtime])*errtime[0].unit
    dnde    = np.asarray([x.item().value for x in dnde])*dnde[0].unit
    errn    = np.asarray([x.item().value for x in errn])*errn[0].unit
    errp    = np.asarray([x.item().value for x in errp])*errp[0].unit
    ftheory = np.asarray([f.value for f in ftheory])*ftheory[0].unit

    if (debug):
        for i,ds in enumerate(dsets):
            print(" Dataset ",i,"--------------" )
            print(" Flux   : {:5.2e} +{:5.2e} -{:5.2e} {:s}"
                      .format(dnde[i].value,
                              errp[i].value,
                              errn[i].value,dnde[i].item().unit))
            print(" Theory : {:5.2e} Ratio theory/data: {}"
                  .format(ftheory[i], ftheory[i].value/dnde[i].value))

    # The plot
    ax = plt.gca() if ax is None else ax
    label= "{:5.0f}-{:5.10} {}" \
           .format(emin.value,emax.to(emin.unit).value,str(emin.unit))

    with quantity_support():
            ax.errorbar(x = time, y = dnde,
                        xerr = errtime,
                        yerr = [errn,errp],
                        color=color, ls="", label=label)

            if (style=="bar"):
                ax.bar(time,ftheory,width=2*errtime,
                       alpha=0.2,color=color,label="Model")
            if (style=="line"):
                ax.plot(time,ftheory,
                        alpha=0.5,color=color,ls=":",marker="o",label="Model")

    import collections
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = collections.OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),bbox_to_anchor=(1,1), loc="upper left")

    return ax

#------------------------------------------------------------------------------
def fluxes_versus_time(dsets, nstep=1,
                       xscale="linear",yscale="log",cmap="brg",debug=False):

    fig,ax = plt.subplots(figsize=(12,8))
    cm = plt.get_cmap(cmap)

    e_edges = dsets[0].excess.geom.axes[0].edges
    nmax = len(e_edges)-1

    for i in range(0,nmax,nstep):
        color = cm(1.*i/nmax)
        if (debug): print("i/nmax",i,"/",nmax,"E:",e_edges[i],e_edges[i+1],"c=",color)

        flux_versus_time(dsets,emin = e_edges[i], emax = e_edges[i+1],
                         style="line",
                         ax=ax, color=color,
                         debug=debug)

    ax.set_xlabel("Elapsed time (" + ax.xaxis.get_label_text() + ")")
    ax.set_ylabel("Flux (" + ax.yaxis.get_label_text()+ ")")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

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

    data_pred   = np.asarray([ds.npred_sig().data.flatten() for ds in dsets])
    data_excess = np.asarray([ds.excess.data.flatten() for ds in dsets])

    fig = plt.figure(figsize=(20,8))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    count_plot3D(dt, E, data_pred,  tag="Prediction",ax3d = ax1)
    count_plot3D(dt, E, data_excess, tag="Excess", ax3d = ax2)

    plt.tight_layout()
    plt.show()

    return

#------------------------------------------------------------------------------
def lightcurve(dsets,
               t_unit = u.s, tmin=0*u.s, tmax=1000*u.s, binwidth=10*u.s,
               tag="excess",style="bar", alpha=1, color="black", fillcolor=None, ax=None, debug=False):

    from scipy.interpolate import interp1d

    tmin = tmin.to(t_unit)
    tmax = tmax.to(t_unit)
    binwidth = binwidth.to(t_unit)

    # Create the time interpolation of the counting rate
    tsamp, dndt, tslice= [], [], []
    t0 = dsets[0].gti.time_start[0]
    duration = 0*t_unit

    for i,ds in enumerate(dsets):
        tsamp.append( ((ds.gti.time_start[0]-t0).sec*u.s + ds.livetime/2).to(t_unit) )
        tslice.append(duration)
        duration = duration + ds.livetime.to(t_unit)

        if (tag == "excess"):
            dndt.append(ds.excess.data[ds.mask_safe].sum()/ds.livetime.to(t_unit))
        elif (tag=="prediction"):
            dndt.append(ds.npred_sig().data[ds.mask_safe].sum()/ds.livetime.to(t_unit))
        else:
            dndt.append(ds.counts.data[ds.mask_safe].sum()/ds.livetime.to(t_unit))

    # From now on, all variables have values corresponding to t_unit
    tsamp = np.asarray([x.value for x in tsamp])
    dndt  = np.asarray([x.value for x in dndt])
    tslice= np.asarray([x.value for x in tslice])

    tmin = tmin.value
    tmax = tmax.value
    binwidth = binwidth.value

    # Warning : f is not a quantity - Therefore all time units should be the same before the interpolation
    # f should have unit 1/t_unit like dndt
    f = interp1d(tsamp,dndt,fill_value="extrapolate")
    #print(tsamp,dndt,f(tsamp))

    # Generate a number of counts for each bin
    nbin = int((tmax-tmin)/binwidth)
    t_edges = np.linspace(tmin,tmax,nbin+1)
    t_bins = t_edges[:-1] + 0.5*(t_edges[1:]-t_edges[:-1])
    #print(t_bins)

    # Compute count number in each bin = rate * bin width"
    n_random = []
    for t in t_bins:
        counts = f(t)*binwidth
        #print(t,"rate=",f(t),"counts=",counts)
        if (tag != "prediction"):
            n_random.append(np.random.poisson(int(counts)))
        else:
            n_random.append(counts)
    #print(n_random)

    ax = plt.gca() if ax is None else ax
    if (style=="bar"):
        ax.bar(t_bins, n_random, width=binwidth, edgecolor=color, color=fillcolor, alpha=alpha,label=tag)
    if (style=="line"):
        ax.plot(t_bins, n_random, color=color,marker="o",lw=1, alpha=alpha,label=tag)
    if (debug):
        ax.plot(t_bins,f(t_bins)*binwidth, ls=":",color="red")
        ax.plot(tsamp[tsamp<tmax],dndt*binwidth,marker="o",ls="",color="green",label="Initial data")

    ax.set_xlabel("Observation duration ("+str(t_unit)+")")
    ax.set_ylabel("Counts per bin (" + str(binwidth) + str(t_unit)+")")
    for xl in tslice[tslice<tmax]:
        ax.axvline(x=xl,ls=":")
    ax.axhline(y=0)
    ax.legend()

    return ax
#------------------------------------------------------------------------------
####  UTILITIES
#------------------------------------------------------------------------------
def tobs_plot(dset):
    """
    A utility to have reasonable duration format  displayed

    Parameters
    ----------
    dset : Dataset
        A current Gammapy Dataset

    Returns
    -------
    tobs : Quantity
        A time with an adapted unit and rounding for plotting

    """
    # Get and format livetimes
    t = dset.livetime.to(u.s)
    if (t.value > 1800): t = t.to(u.h)
    if (t.value > 3600*24): t = t.to(u.day)
    tobs = str( round(t.value,2)) +" "+ str(t.unit)

    return tobs

#------------------------------------------------------------------------------
def count_plot3D(x, y, data,
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

#------------------------------------------------------------------------------
def counts_from_plot(ax, dset, plot_line=False):
    """
    THis a utility to get the counts form the plot and check the count_plot
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

#------------------------------------------------------------------------------
def panels(dsets,func=None,
               nmaxcol=5, xsize=5,ysize=7,
               xscale ="log", yscale="log",
               fixrange=False, **kwargs):
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
                           figsize=(min(20,xsize*ncols),ysize*nrows),
                           sharex=True, sharey=fixrange)

    with quantity_support():
        iplot = 0

        import itertools
        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            ax0 = ax[jrow][icol] if (nrows>1) else ax[icol]

            # If nothing to plot, blank space instead of empty plot
            if (iplot >= nplots):
                ax0.axis('off')
                continue # Next plot

            # Plot function here
            func(dsets[iplot],ax=ax0,**kwargs)

            # Compactify
            if (jrow+1 != nrows): ax0.set_xlabel(None)
            if (icol !=0): ax0.set_ylabel(None)
            ax0.tick_params(which='major', length=10, width=2, direction='in')
            ax0.tick_params(which='minor', length=5, width=2, direction='in')

            ax0.set_xscale(xscale)
            ax0.set_yscale(yscale)

            iplot+=1

    if (fixrange):
        # Adjust the scaling factor to fit your legend te
        plt.subplots_adjust(right=0.8,hspace=0.05,wspace=0.05)
    else:
        plt.subplots_adjust(right=0.8,hspace=0.05)

    # use first label handles as the last plot can be empty
#     h0, lbl0 = fig.axes[0].get_legend_handles_labels()
#     fig.legend(h0, lbl0, loc="center right", borderaxespad=0.1)
    if (fixrange): fig.tight_layout(h_pad=0)
    else: fig.tight_layout(h_pad=0, w_pad=0)
    plt.show()

    return