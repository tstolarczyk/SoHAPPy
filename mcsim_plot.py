# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes
# from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import astropy.units as u
from   astropy.visualization import quantity_support
from   astropy.coordinates   import AltAz

# from utilities import MyLabel, stamp

# Defaults colors for the plots
col3  = "tab:orange" # 3 sigma
col5  = "tab:green"  # 5 sigma
colmx = "black"      # Sigma max

# Bigger texts and labels
#plt.style.use('seaborn-talk')
plt.style.use('seaborn-poster') # Bug with normal x marker !!!

# If False, avoid scirpt to be paused when a plot is popped-up (plt.show)
block = False

###############################################################################
def show(mc,ax=None, loc="nowhere",pdf=None):
    
    fig_sig = sigma_vs_time(mc)
    # if (mc.niter > 1): fig_n = non_vs_noff(mc,ax)
    # else: fig_n = None # Otherwise plot one single point
        
    # fig_vis = story(mc,loc=loc,ref="VIS")
    # fig_grb = story(mc,loc=loc,ref="GRB")
    
    if pdf != None:
        pdf.savefig(fig_sig)
        # if fig_n != None: pdf.savefig(fig_n)
        # pdf.savefig(fig_vis)
        # pdf.savefig(fig_grb)
        
    # used to pause plot display in interactive mode on a shell script. In the
    # abscence of a call to that function figures will stack on the screen during
    # the run and all disappear at the end of the run.
    # Using this, figures will be stacked on screen and displayed for each event
    # until they are closed by the user.
    plt.show(block=True)

    return

###----------------------------------------------------------------------------
def sigma_vs_time(mc):
        
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15,5),
                                       gridspec_kw={'width_ratios': [2, 1]})
   
    # Measurement points
    t_s = np.asarray([s.tobs().value for s in mc.slot.slices])
    
    # Max siginificance and error and time
    tmax     = [mc.slot.slices[i].tobs().value for i in mc.id_smax_list]
    smax     = np.mean(mc.smax_list)
    err_smax = np.std(mc.smax_list)
    
    # 3 sigma and 5 sigma times
    t3s = [mc.slot.slices[i].tobs().value for i in mc.id_3s_list]
    t5s = [mc.slot.slices[i].tobs().value for i in mc.id_5s_list]    
    
    with quantity_support():

        ### Mean significance at each slice
        ax1.errorbar(t_s, mc.sigma_mean, yerr = mc.sigma_std, fmt  ='o')
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()

        tlist    = [tmax, t3s, t5s]
        siglist  = [smax, 3, 5]
        errlist  = [err_smax, 0, 0]
        collist  = [colmx, col3, col5]
        taglist  = ["$\sigma_{max}$","$3\sigma$","$5\sigma$"]
        
        for t, sig, errs, c, tag  in zip(tlist,
                                         siglist,
                                         errlist,
                                         collist,
                                         taglist):
            ax1.errorbar(x    = np.mean(t), y    = sig,
                         xerr = np.std(t),  yerr = errs,
                         fmt="o",color=c,
                         label = tag)
            ax1.vlines(np.mean(t), ymin = ymin, ymax = smax,
                       alpha=0.5,ls=":",color=c)
            ax1.hlines(sig, xmin=xmin,ls=":", xmax=np.mean(t),
                       alpha=0.5,color=c)

    
        ax1.set_xlabel('Observation duration (s)')
        ax1.set_ylabel("Significance $\sigma$")
        # ax1.legend(bbox_to_anchor=[1,1])
        ax1.legend()

        import matplotlib
        if matplotlib.__version__ >= "3.5.0":
            ax1.set_xscale("log", nonpositive='clip') # 3.5.0
        else:
            ax1.set_xscale("log", nonposx='clip') # 3.1.1
        ax1.grid(which='both',alpha=0.2)    
        
        # Sigma max distribution (if niter > 1)
        # if (mc.niter >1):
        # axx = ax.inset_axes([0.3,0.15,0.2,0.3]) # lower corner x,y, w, l
        n, bins,_ = ax2.hist(mc.smax_list,
                    bins  = max(int(mc.niter/2),1), # Cannot be < 1
                    range = [smax-3*err_smax,smax+3*err_smax],
                    alpha=0.5,
                    color = "grey",
                    label = " {:5.1} $\pm$ {:5.1}".format(smax,err_smax))
        ax2.set_xlabel("$\sigma_{max}$")
        ax2.set_ylabel('Trials')

        
        ax1.set_title(mc.name +' ('+str(mc.niter)+' iter.)',loc="right")
        
    plt.tight_layout()
    
    return fig

###----------------------------------------------------------------------------
def non_vs_noff(mc,ax,alpha=1/5,logscale=True):
    """
    Plot Non versus Noff from the background and excess counts.
    Draw the error bars from the variances
    """
    if ax == None:
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    else:
        fig = ax.get_figure()

    nexlist = [mc.nex_3s_list, mc.nex_5s_list, mc.nex_smax_list ]
    nblist  = [mc.nb_3s_list,  mc.nb_5s_list,  mc.nb_smax_list]
    collist = [col3, col5, colmx]
    taglist = ["$3\sigma$","$5\sigma$","$\sigma_{max}$"]
        
    for nex, nb, col, tag in zip(nexlist, nblist, collist, taglist):

        non       = np.add(nex,nb)
        noff      = np.true_divide(nb,alpha)
        non_mean  = np.mean(non)
        non_std   = np.std(non)
        noff_mean = np.mean(noff)
        noff_std  = np.std(noff) 
        xlabel    = "$N_{off}$"
        ylabel    = "$N_{on}$"
        if logscale:
            non       = np.log10(non)
            non_mean  = np.log10(non_mean)
            non_std   = np.log10(non_std)
            noff      = np.log10(noff)
            noff_mean = np.log10(noff_mean)
            noff_std  = np.log10(noff_std)  
            xlabel    = "$log_{10}$" + xlabel
            ylabel    = "$log_{10}$" + ylabel
        
        ax.plot(noff, non, 
                alpha=0.3, marker= ".", ls ="", color=col,
                label=tag )
        ax.errorbar(x    = noff_mean,    y    = non_mean,
                xerr = np.std(noff), yerr = np.std(non),
                color = col)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()       
        
    plt.tight_layout()

    return fig

###############################################################################
def story(mc, loc="nowhere", ref="VIS"):
    """
    Show altitude versus time and points used for computation

    """
    if (loc != "North" and loc!= "South"): return None

    # Plot sigma versus time on the full
    # GRB lifetime together with visibility
    grb = mc.slot.grb

    coord = grb.vis[loc].site
    dtplot = 500*u.s # Add some space for the plot


    # This are the visibilities in the sense of the Dark time
    # Allows seing trigger within the night
    # Create points from the stat of teh firts night
    # to the end of the last night
    tvis1 = grb.vis[loc].t_twilight[0][0]
    tvis2 = grb.vis[loc].t_twilight[-1][1]
    dtvis = (tvis2-tvis1).to(u.s)
    tvis  = np.linspace(0,dtvis.value,100)*dtvis.unit

    # Set the reference time
    event_name = grb.name+" - " + loc
    if (ref=="VIS"):
        tref = tvis1
        #text_ref =event_name + " (Ref: vis. start)"
        text_ref =event_name + " (Ref: Dark time start)"
    else:
        tref = grb.t_trig
        text_ref = event_name+ " (Ref: GRB trigger)"

    t_trig_rel = grb.t_trig - tref

    # Recompute time sampling in the two reference systems
    tvis_abs = tvis1 + tvis
    tvis_rel = tvis_abs - tref

    # Compute alt-az points for visibility time sampling
    altazvis = grb.radec.transform_to( AltAz(obstime= tvis_abs,
                                                  location=coord))

    # Compute GRB measurement points in the two refrence systems
    tgrb_abs   = grb.t_trig + grb.tval # Absolute time
    tgrb_rel   = tgrb_abs - tref
    altazgrb   = grb.radec.transform_to( AltAz(obstime= tgrb_abs,
                                                    location=coord))

    # Compute GRB observation points
    t_s = [s.tobs() for s in mc.slot.slices]*mc.slot.slices[0].tobs().unit
    tgrb_obs_abs = grb.t_trig + t_s # Absolute time
    tgrb_obs_rel = tgrb_obs_abs - tref
    altazobs     = grb.radec.transform_to( AltAz(obstime= tgrb_obs_abs,
                                                    location=coord))

    # Get the max significancae and the detection time if they exist
    from mcsim_res import mean_values
    tmx,e_tmx,altmx,e_altmx,azmx, e_azmx = mean_values(mc,mc.id_smax_list)
    t3s,e_t3s,alt3s,e_alt3s,az3s, e_az3s = mean_values(mc,mc.id_3s_list)
    t5s,e_t5s,alt5s,e_alt5s,az5s, e_az5s = mean_values(mc,mc.id_5s_list)

    with quantity_support():
        fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(10,6))

        # Plot limits
        xmin = -dtplot
        xmax = (tvis2 -tref + dtplot ).sec
        ymin = 0*u.deg

        # Relative to reference
        ax1.plot(tvis_rel.sec, altazvis.alt,
                 linestyle="--",
                 lw=1,
                 color="tab:blue",
                 label="Altitude")

        ax1.plot(tgrb_rel.sec,altazgrb.alt,
                 linestyle="",
                 marker="*",
                 markerfacecolor="b",
                 label="End of intervals")

        ax1.plot(tgrb_obs_rel.sec,altazobs.alt,
                 linestyle="",
                 marker="o",
                 markersize=10,
                 markerfacecolor="none",
                 markeredgewidth = 1.5,
                 markeredgecolor ="r",
                 label="Measurements")

        # Trigger
        ax1.axvline(x=t_trig_rel.sec,
                    ls=(0,(1,1)), # densely dotted
                    color="grey",
                    label="Trigger")

        # Dark time
        first = True

        for elt in grb.vis[loc].t_twilight:
            if (first):
                label="Night"
                first= False
            else:
                label = None
            ax1.axvspan((elt[0] - tref).sec,
                       (elt[1] - tref).sec,
                       alpha=0.2,color="black",label=label)


        # minimal altitude in simulation
        ax1.axhline(y=mc.slot.grb.vis[loc].altmin,
                    ls="dashdot",
                    color="grey")

        # Sig max
        ax1.errorbar(tmx.to(u.s).value + t_trig_rel.sec,
             altmx.value,
             xerr=e_tmx.to(u.s).value,
             yerr=e_altmx.value,
             fmt='o',color=colmx,
                     label="$\sigma_{max}$")
        ax1.vlines(tmx.to(u.s).value + t_trig_rel.sec,
                   ymin=0,
                   ymax=altmx.value,
                   alpha=0.5,lw=1,color=colmx)

        ax1.hlines(altmx.value,
                   xmin=xmin,
                   xmax=tmx.to(u.s).value + t_trig_rel.sec,
                   alpha=0.5,lw=1,color=colmx)

        # 3 sigma - if reached
        if (t3s.to(u.s).value > 0):
            ax1.errorbar(t3s.to(u.s).value + t_trig_rel.sec,
                 alt3s.value,
                 xerr=e_t3s.to(u.s).value,
                 yerr=e_alt3s.value,
                 fmt='o',color=col3,
                 label="$3\sigma$")

            ax1.vlines(t3s.to(u.s).value + t_trig_rel.sec,
                       ymin=0,
                       ymax=alt3s.value,lw=1,color=col3)

            ax1.hlines(alt3s.value,
                       xmin=xmin,
                       xmax=t3s.to(u.s).value + t_trig_rel.sec,
                       alpha=0.5,lw=1,color=col3)

        # 5 sigma - if reached
        if (t5s.to(u.s).value > 0):
            ax1.errorbar(t5s.to(u.s).value + t_trig_rel.sec,
                         alt5s.value,
                         xerr=e_t5s.to(u.s).value,
                         yerr=e_alt5s.value,
                         fmt='o',color=col5,
                         label="$5\sigma$")

            ax1.vlines(t5s.to(u.s).value + t_trig_rel.sec,
                       ymin=0,
                       ymax=alt5s.value,lw=1,color=col5)

            ax1.hlines(alt5s.value,
                       xmin=xmin,
                       xmax=t5s.to(u.s).value + t_trig_rel.sec,
                       alpha=0.5,lw=1,color=col5)

        # Title and limits
        ax1.set_title(text_ref,fontsize=12,loc="right")
        ax1.set_xlim(xmin=xmin,xmax= xmax)
        ax1.set_ylim(ymin=ymin)

        if (ref=="VIS"):
            ax1.set_xlabel("Elapsed time since visible (s)")
        else:
            ax1.set_xlabel("Elapsed time since Trigger (s)")
        ax1.set_ylabel("Altitude")

        #ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.legend(fontsize=10)

        # Display second axis
        ax = ax1.twiny()
        ax.plot(tvis_abs.datetime,altazvis.alt,
                 linestyle="",
                 color="b")
        ax.tick_params(axis='x', rotation=70)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

        if (ref=="VIS"): ax.set_xlabel("Date since twilight (hh:mm)")
        else:            ax.set_xlabel("Date since Trigger (hh:mm)")

    plt.tight_layout()
    
    return fig

###############################################################################
#
# Plot one trial
#
###############################################################################
def onetrial(dsets, slot):
    """
    

    Parameters
    ----------
    dsets : TYPE
        DESCRIPTION.
    slot : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from gammapy.estimators import FluxPointsEstimator

    print(" One trial plots")
    nplots = len(dsets)
    ncols = 5
    nrows = int((nplots)/ncols)+1
    fig, ax = plt.subplots(ncols=ncols,nrows=nrows,figsize=(20,5*nrows))

    iplot = 0
    for jrow in range(nrows):
        for icol in range(ncols):
        #print("col=",i," row=",j)
            if (iplot < nplots):
                ax0 = ax[jrow][icol]
                dsets[iplot].plot_counts(ax=ax0)
                ax0.legend()
                if (jrow<nrows-1):
                    ax0.axes.get_xaxis().set_visible(False)
                iplot+=1
    fig.tight_layout(h_pad=0)
    plt.show()

    fig, ax = plt.subplots(ncols=ncols,nrows=nrows,figsize=(15,5*nrows))
    iplot = 0
    for jrow in range(nrows):
        for icol in range(ncols):
        #print("col=",i," row=",j)
            if (iplot < nplots):
                dset = dsets[iplot]
                slice = slot.slices[iplot]
                ax0 = ax[jrow][icol]

                fpe = FluxPointsEstimator(
                    e_edges=slice.irf()[0].ereco.edges)
                flux_points = fpe.run(datasets=dset)
                flux_points.table["is_ul"] = flux_points.table["ts"] < 4
                flux_points.plot(energy_power=2,
                          flux_unit="erg-1 cm-2 s-1",
                          color="tab:blue",ax=ax0)
                spectrum = slot.grb.spectra[slice.fid()]
                t = slice.tobs().to(u.s)
                if (t.value > 3600):
                    t = t.to(u.h)
                if (t.value > 3600*24):
                    t = t.to(u.day)
                tobs = str( round(t.value,2)) + str(t.unit)
                spectrum.plot(energy_range=dset.energy_range,
                                  flux_unit='cm-2 s-1 erg-1',
                                  energy_power=2,
                                  energy_unit='TeV',
                                  n_points=10,
                                  ax=ax0,ls=":",color="red",marker="",
                                  label=tobs)
                ax0.legend()
                if (jrow<nrows-1):
                    ax0.axes.get_xaxis().set_visible(False)
                iplot+=1
                if (icol !=0): ax0.set_ylabel(None)
    fig.tight_layout(h_pad=0)
    plt.show()
    # for dset, slice in zip(dsets, slot.slices):
    #     fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols= 3, figsize=(16,5))

    #     # Extract spectrum
    #     # flux_points.to_sed_type("e2dnde").plot_ts_profiles(ax=ax3,
    #     #                                                    cmap="Blues",
    #     #                                                    alpha=0.5)
    #     spectrum = slot.grb.spectra[slice.fid()]
    #     spectrum.plot(energy_range=dset.energy_range,
    #                               flux_unit='cm-2 s-1 erg-1',
    #                               energy_power=2,
    #                               energy_unit='TeV',
    #                               n_points=10,
    #                               ax=ax3,ls=":",color="red",marker="",
    #                               label="Theory")
    #     ax3.legend()

    return



#    @staticmethod
#    def plot(simu, target):
#        """Plot some simulation results"""
#
#        import matplotlib.pyplot as plt
#        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
#                                       figsize=(10, 5))
#
#        # Spectrum plot
#        energy_range = [0.01 * u.TeV, 100 * u.TeV]
#        target.model.plot(ax=ax1, energy_range=energy_range,
#                          label='Model')
#        plt.text(0.55, 0.65, target.__str__(),
#                 style='italic', transform=ax1.transAxes, fontsize=7,
#                 bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
#        ax1.set_xlim([energy_range[0].value, energy_range[1].value])
#        ax1.set_ylim(1.e-17, 1.e-5)
#        ax1.grid(which='both')
#        ax1.legend(loc=0)
#
#        # Counts plot
#        on_off = simu.on_vector.data.data.value
#        off = 1. / simu.off_vector.backscal * simu.off_vector.data.data.value
#        excess = on_off - off
#        bins = simu.on_vector.energy.lo.value
#        x = simu.on_vector.energy.nodes.value
#        ax2.hist(x, bins=bins, weights=on_off,
#                 facecolor='blue', alpha=1, label='ON')
#        ax2.hist(x, bins=bins, weights=off,
#                 facecolor='green', alpha=1, label='OFF')
#        ax2.hist(x, bins=bins, weights=excess,
#                 facecolor='red', alpha=1, label='EXCESS')
#        ax2.legend(loc='best')
#        ax2.set_xscale('log')
#        ax2.set_xlabel('Energy [TeV]')
#        ax2.set_ylabel('Expected counts')
#        ax2.set_xlim([energy_range[0].value, energy_range[1].value])
#        ax2.set_ylim([0.0001, on_off.max() * (1 + 0.05)])
#        ax2.vlines(simu.lo_threshold.value, 0, 1.1 * on_off.max(),
#                   linestyles='dashed')
#        ax2.grid(which='both')
#        ax2.text(0.55, 0.05, simu.__str__(),
#                 style='italic', transform=ax2.transAxes, fontsize=7,
#                 bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
#        plt.tight_layout()
#
#------------------------------------------------------------------------------
if __name__ == "__main__":

    import pickle
    import mcsim_res as mcres

    mcdatafile = "../output/Today/Event398-North-100.bin"
#    mcdatafile = "../output/Today/Event501-South-1000.bin"
    print(" Simulation read from ",mcdatafile)
    infile  = open(mcdatafile,"rb")
    mc      = pickle.load(infile)
    infile.close()

#    story(mc, ref="VIS",saveplots="False",outfile="")
#    story(mc, ref="notVIS",saveplots="False",outfile="")
    stat(mc,saveplots="False",outfile="")
    mcres.result(mc)
