# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from   pathlib import Path

import astropy.units as u
from   astropy.visualization import quantity_support
from   astropy.coordinates   import AltAz

from gammapy.stats import background_error

import mcsim_config as mcf
import ana_config as cf

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
def pause():
    """
    Used to pause plot display in interactive mode on a shell script. In the
    abscence of a call to that function figures will stack on the screen during
    the run and all disappear at the end of the run.
    Usingthis, figures will be stacked on screen and displayed for each event
    until they are closed by the user.

    Returns
    -------
    None.

    """
    plt.show(block=True)
    return

###############################################################################
def show(mc,show=2,loc="nowhere"):
    filename = Path(cf.res_dir,mc.name)
    stat(mc,show,filename)
    story(mc,show,filename,loc=loc,ref="VIS")
    story(mc,show,filename,loc=loc,ref="GRB")
    pause()
    
    return
###############################################################################
def plot_on_off(ax,nex,nb,
                color = None,
                desc  = None,
                marker=".",
                log   = False):
    """
    Plot Non versus Noff from the background and excess counts.
    Draw the error bars from the variances
    """
    non       = np.add(nex,nb)
    noff      = np.true_divide(nb,mcf.alpha)
    
    non_mean  = np.mean(non)
    non_std   = np.std(non)
    noff_mean = np.mean(noff)
    noff_std  = np.std(noff)
    
    if (log == True):
        non       = np.log10(non)
        non_mean  = np.log10(non_mean)
        non_std   = np.log10(non_std)
        noff      = np.log10(noff)
        noff_mean = np.log10(noff_mean)
        noff_std  = np.log10(noff_std)  
        
        
    ax.plot(noff, non,
            alpha=0.5,
            marker=marker,
            markersize=20,
            markerfacecolor="none",
            markeredgewidth = 1.5,
            markeredgecolor =color,
            ls="",
            label=desc)
    ax.errorbar(x    = noff_mean,    y    = non_mean,
                xerr = np.std(noff), yerr = np.std(non),
                color = color)        

    return
    
###############################################################################
def stat(mc, saveplots, outfile):
    """
    Plot result of the simulations of the current grb:
        - evolution of sigma with observation time
        - distibutions of Non versus Noff
    """

    if (saveplots == 0): return
    if (mc.ndof != 0):
        sig = '$\sigma_{\sqrt{TS}}$'
    else:
        sig = '$\sigma_{Li&Ma}$'

    with quantity_support():
        
        # Compute relevant quantities
        eventid = mc.name
        t_s = [s.tobs() for s in mc.slot.slices]

        tmax = [mc.slot.slices[i].tobs().value
                 for i in mc.id_smax_list
                 if mc.slot.slices[i].tobs().value>=0]
        smax = np.max(mc.smax_list)
        err_smax = np.std(mc.smax_list)
        t3s = [mc.slot.slices[i].tobs().value
                 for i in mc.id_3s_list
                 if mc.slot.slices[i].tobs().value>=0]
        t5s = [mc.slot.slices[i].tobs().value
                 for i in mc.id_5s_list
                 if mc.slot.slices[i].tobs().value>=0]

        fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,figsize=(10,10))
        
        ### Mean significance at each slice
        ax1.errorbar(t_s,
                     mc.sigma_mean,
                     yerr=mc.sigma_std,
                     fmt='o')
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        
        ### Max significance
        ax1.errorbar(x    = np.mean(tmax), y    = smax,
                     xerr = np.std(tmax),  yerr = err_smax,
                     fmt="+",color=colmx,
                     label = "$\sigma_{max}$")
        ax1.vlines(np.mean(tmax), ymin = ymin, ymax = smax,
                   alpha=0.5,ls=":",color=colmx) 
        ax1.hlines(smax, xmin=xmin,ls=":", xmax=tmax,
                   alpha=0.5,color=colmx)
        
        ### 3 sigma
        ax1.errorbar(x = np.mean(t3s),y = 3.,xerr = np.std(t3s),
                     fmt="o",color=col3,
                     label = "$3\sigma$")
        ax1.vlines(np.mean(t3s), ymin = ymin, ymax = 3,
                   alpha=0.5,ls=":",color=col3) 
        ax1.hlines(3, xmin=xmin,ls=":", xmax=np.mean(t3s),
                   alpha=0.5,color=col3)        
        
        ### 5 sigma
        ax1.errorbar(x = np.mean(t5s),
                     y = 5.,
                     xerr = np.std(t5s),
                     fmt="o",color=col5, 
                     label = "$5\sigma$")
        ax1.vlines(np.mean(t5s), ymin = ymin, ymax = 5,
                   alpha=0.5,ls=":",color=col5) 
        ax1.hlines(5, xmin=xmin,ls=":", xmax=np.mean(t5s),
                   alpha=0.5,color=col5) 
               
        ax1.set_xlabel('Observation duration (s)')
        ax1.set_ylabel(sig)
        ax1.set_title(eventid +' ('+str(mc.niter)+' iter.)')
        ax1.set_xscale("log", nonposx='clip')
        ax1.grid(which='both',alpha=0.2)
        
        if (mc.niter >1):
            ax11 = ax1.inset_axes([0.1,0.5,0.3,0.5]) # lower corner x,y, w, l
            ax11.hist(mc.smax_list,
                      bins  = max(int(mc.niter/2),1), # Cannot be < 1
                      range = [smax-3*err_smax,smax+3*err_smax],
                      alpha=0.5,
                      color = "grey",
                      label = " {:5.1} $\pm$ {:5.1}".format(smax,err_smax))

            ax11.set_xlabel(sig,fontsize=12)
            ax11.set_ylabel('Trials',fontsize=12)
        
        ax1.legend(loc="lower right")
    
        ### Non, Noff - Log
        logscale = True
        plot_on_off(ax2,mc.nex_3s_list,mc.nb_3s_list,
                    color=col3,
                    desc = "$3\sigma$",
                    log=logscale)
        
        plot_on_off(ax2,mc.nex_5s_list,mc.nb_5s_list,
                    color=col5,
                    desc = "$5\sigma$",
                    log=logscale)
        
        plot_on_off(ax2,mc.nex_smax_list,mc.nb_smax_list,
                    color=colmx,
                    desc = "$\sigma_{max}$",
                    marker="+",
                    log=logscale)
        
        if (logscale):
            ax2.set_xlabel("$log<N_{off}>$")
            ax2.set_ylabel("$log<N_{on}>$") 
        else:
            ax2.set_xlabel("$<N_{off}>$")
            ax2.set_ylabel("$<N_{on}>$") 
        ax2.legend()
        
        plt.tight_layout()
        plt.show(block=block)

        if (saveplots==2 or saveplots==3):
            fig.savefig(Path(outfile).with_suffix('.pdf'))
            # unsupported fig.savefig(outfile.with_suffix('.jpg'))
    return

###############################################################################
def story(mc, saveplots, outfile,loc="nowhere",ref="VIS"):
    """
    Show altitude versus time
    
    """
    if (loc != "North" and loc!= "South"): return
    
    # Plot sigma versus time on the full
    # GRB lifetime together with visibility
    grb = mc.slot.grb
    
    coord = grb.pos_site[loc]
    dtplot = 500*u.s # Add some space for the plot

    # Visibility altitude curve
    # This is the "visibility" related to the appairion of the GRB
    # i.e. in dark time it cannot start before the trigger
    #tvis1 = grb.t_true[loc][0] # End of visibility
    #tvis2 = grb.t_true[loc][1] # End of visibility
    
    # This is the visibility in the sense of the Dark time
    # Allows seing trigger within the night
    tvis1 = grb.t_twilight[loc][0][0]
    tvis2 = grb.t_twilight[loc][0][1]
    dtvis = (tvis2-tvis1).to(u.s)
    tvis = np.linspace(0,dtvis.value,100)*dtvis.unit

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
    tgrb_obs_abs   = grb.t_trig + t_s # Absolute time
    tgrb_obs_rel   = tgrb_obs_abs - tref
    altazobs   = grb.radec.transform_to( AltAz(obstime= tgrb_obs_abs,
                                                    location=coord))
    
    # Get the max significancae and the detection time if they exist
    from mcsim_res import mean_values
    tmx,e_tmx,altmx,e_altmx,azmx, e_azmx = mean_values(mc,mc.id_smax_list)   
    t3s,e_t3s,alt3s,e_alt3s,az3s, e_az3s = mean_values(mc,mc.id_3s_list)
    t5s,e_t5s,alt5s,e_alt5s,az5s, e_az5s = mean_values(mc,mc.id_5s_list)

    with quantity_support():
        fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(10,8))
        
        # Plot limits
        xmin = -dtplot
        xmax = (tvis2 -tref + dtplot ).sec
        ymin = 0*u.deg

        # Relative to reference
        ax1.plot(tvis_rel.sec, altazvis.alt,
                 linestyle="--",
                 color="b",
                 label="Altitude")

        ax1.plot(tgrb_rel.sec,altazgrb.alt,
                 linestyle="",
                 marker="*",
                 markerfacecolor="b",
                 label="GRB end of intervals")

        ax1.plot(tgrb_obs_rel.sec,altazobs.alt,
                 linestyle="",
                 marker="o",
                 markersize=10,
                 markerfacecolor="none",
                 markeredgewidth = 1.5,
                 markeredgecolor ="r",
                 label="GRB measurements")

        # Trigger
        ax1.axvline(x=t_trig_rel.sec,
                    ls=(0,(1,1)), # densely dotted
                    color="grey",
                    label="Trigger")

        # Dark time
        ax1.axvline(x=(tvis1 - tref).sec,
                    ls=":",
                    color="green",
                    label="Start dark time")
        ax1.axvline(x=(tvis2 - tref).sec,
                    ls=":",
                    color="red",
                    label="Stop dark time")
        
        # minimal altitude in simulation
        ax1.axhline(y=mc.slot.grb.altmin,
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
        
        # 3 sigma
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

        # 5 sigma
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
        
        
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))        
        
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
    plt.show(block=block)
    return

###############################################################################
#
# Plot one trial
#
###############################################################################
def onetrial(mc):
    
    """
    Plot some detection statistics for the current trials 
    From J. Lefaucheur, a long time ago
     
    The plots can be rearranged giving the various axis values
    """

    from fit_onoff import cumulative_OnOffstats   
    from gammapy.stats import excess_error
    stats = cumulative_OnOffstats(mc.simulations)
     
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(18, 18))
    
    obstime = stats['livetime']
    
    # 1- Cumulated Non
    ax1.errorbar(obstime,
                 stats['n_on'], 
                 yerr = stats['n_on']/np.sqrt(mc.niter),
                 color='black', 
                 fmt='o')
    ax1.set_xlabel('Obs. time [s]', fontweight='bold')
    ax1.set_ylabel('$N_{on}$', fontweight='bold')
    ax1.set_title('On region count evolution', fontweight='bold')
    ax1.axhline(y=10,color="red",ls="--")
    ax1.grid(which='both')
    ax1.set_xscale('log')  
    
    # 2 - Cumulated excess
    yerr = excess_error(stats['n_on'],stats['n_off'],stats['alpha'])
    ax2.errorbar(obstime, stats['excess'],
                       yerr=yerr,
                       color='black', fmt='o')
    ax2.set_xlabel('Livetime [s]', fontweight='bold')
    ax2.set_ylabel('Excess counts', fontweight='bold')
    ax2.set_title('Cumulated excess', fontweight='bold')
    ax2.set_ylim(ymax=(stats['excess'] + yerr).max() * (1.1))
    ax2.grid(which='both')
    ax2.set_xscale('log')

    # 3- Cumulated Noff
    ax3.errorbar(obstime,
                 stats['n_off'], 
                 yerr = stats['n_off']/np.sqrt(mc.niter),
                 color='black', 
                 fmt='o')
    ax3.set_xlabel('Obs. time [s]', fontweight='bold')
    ax3.set_ylabel('$N_{off}$', fontweight='bold')
    ax3.set_title('Off region count evolution', fontweight='bold')
    ax3.axhline(y=10,color="red",ls="--")
    ax3.grid(which='both')
    ax3.set_xscale('log')  
    
    # 4- Cumulated background
    yerr = background_error(stats['n_off'],stats['alpha'])
    ax4.errorbar(obstime, stats['bkg'],
                     yerr=background_error(
                         stats['n_off'],
                         stats['alpha']
                     ),
                     color='black', fmt='o')

    ax4.set_xlabel('Obs. time [s]', fontweight='bold')
    ax4.set_ylabel('Background count', fontweight='bold')
    ax4.set_title('Cumulated background', fontweight='bold')
    ax4.set_ylim(stats['bkg'].min() * 0.9, (stats['bkg'] + yerr).max() * (1.1))
    ax4.grid(which='both')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    
    plt.show(block=block)
    
    # Significance evolution
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8, 6))
    ax.errorbar(obstime, stats['sigma'],color='black', fmt='o')
    ax.set_xlabel('Obs. time [s]', fontweight='bold')
    ax.set_ylabel('Significance', fontweight='bold')
    ax.set_title('Significance (Li & Ma)', fontweight='bold')
    #♠ax.set_ylim(0., stats['sigma'].max() * (1.1))
    ax.set_ylim(ymax=np.nanmax(stats['sigma']) * (1.1))
    ax.grid(which='both')
    ax.set_xscale('log') # CHANGED
    
    plt.show(block=block)

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


