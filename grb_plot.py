# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates

import numpy as np
import warnings

import astropy.units as u
from   astropy.coordinates   import AltAz
from   astropy.coordinates import get_moon

from   astropy.visualization import quantity_support

from   astroplan import moon_illumination, moon_phase_angle

from utilities import single_legend, t_str

# Bigger texts and labels
#plt.style.use('seaborn-talk')
plt.style.use('seaborn-poster') # Bug with normal x marker !!!

# If False, avoid scirpt to be paused when a plot is popped-up (plt.show)
block = False
#-----------------------------------------------------------------------------#
coln = "blue" # color for North site plots
cols = "red"  # color for South site plots
#-----------------------------------------------------------------------------#

__all__ = ['spectra',
           'energy_over_timeslices',
           'time_over_energyband',
           'energy_and_time_2d',
           'animated_spectra',
           'visibility_plot'
           ]

###-------------------------
def F(x):
    return x.datetime
###-------------------------

###############################################################################
def pause():
    """
    Used to pause plot display in interactive mode on a shell script. In the
    abscence of a call to that function figures wil stack on the screen during
    the run and all disappear at the end of the run.
    Using this, figures will be stacked on screen and displayed for each event
    until they are closed by the user.

    Returns
    -------
    None.

    """
    plt.show(block=True)
    return
###############################################################################
#
# Plot Energy and time spectra
#
###############################################################################

###---------------------------------------------------------------------------
def energy_spectra(grb, n_t_2disp = 5, ymin = 1e-16,  e_min = 1*u.GeV, ax=None):
    """
    Energy spectra for various measurement points in time
    Note: the plot methods handle the labels, should not be overwritten   

    Parameters
    ----------
    grb : GammaRayBurst instance
        A GRB.
    n_E_2disp : TYPE, optional
        Number of curves to show in energy spectra. The default is 5.
    ax : matplotlib.axes, optional
        Current matplotlib axis. The default is None.

    Returns
    -------
    None.

    """

    if ax == None:
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,6))

    # Compute number of spectra to be shown
    nspectra = len(grb.tval)-1
    if (nspectra > n_t_2disp):
        dnt = int(round(nspectra/n_t_2disp))
        tidx = list(range(0,nspectra,dnt))
        # tidx = list(range(0,grb.id90))
        # tidx = [12, 13, 14, 15] + tidx # Typical Prompt times
    else:
        dnt=1
        tidx = [1]
    
 
    # Plot in color some of the energy spectra and the initial data points
    with quantity_support():
        e_unit = grb.Eval[0].unit
        # In the present version, there is only one prompt spectrum up to t90
        if grb.prompt: # if no prompt id90 =-1
            ax.plot(grb.E_prompt.to(e_unit), grb.flux_prompt, 
                    alpha=0.5, ls="--", color = "black", 
                    label="$Prompt \ 0-t_{90}$")
            
        for i in tidx:
            t    = grb.tval[i]
            flux = grb.fluxval[i,:] # Afterglow
            
            # Plot the GRB interpolated spectra
            grb.spectra[i].plot([min(grb.Eval),max(grb.Eval)], 
                                ax, energy_unit=e_unit,
                                label="t= {:>s}".format(t_str(t)))
            
            # Plot the initial data points (should fit)
            c = ax.lines[-1].get_color() # Last color
            ax.plot(grb.Eval,
                    flux.to(1/u.cm**2/u.s/u.TeV) ,
                      ls='--', lw=1., marker=".", alpha=0.5, color=c)   
        
        # # Plot the rest faded out
        for i in range(0,nspectra):
            t = grb.tval[i]
            grb.spectra[i].plot([min(grb.Eval),max(grb.Eval)],ax,
                                energy_unit = e_unit,
                                alpha=0.2, color="grey",lw=1)
            
        ax.axvline(10*u.GeV, color="grey",alpha=0.2)
        ax.text(x = 10*u.GeV, y=ymin*50, s="10 GeV",
                  rotation=270, fontsize=14,va="bottom")
        title = "{}: {:>2d} Flux points".format(grb.name,len(grb.tval))
            
    if (n_t_2disp<=15): ax.legend(fontsize=12) # Too many t slices
    ax.set_title(title,fontsize=12)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(ymin=ymin)
    ax.set_xlim(xmin=e_min.to(e_unit).value)
    
    return
###---------------------------------------------------------------------------
def time_spectra(grb, n_E_2disp = 6, 
                      e_min=10*u.GeV, e_max=10*u.TeV, 
                      ymin=1.e-20, ax=None, debug=False):

    if ax == None: fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,6))

    E = np.logspace(np.log10(e_min.value),
                    np.log10(e_max.to(e_min.unit).value),
                    n_E_2disp)*e_min.unit   
    # print(E)


    ### -----------------------------------
    ### Compute the partial and total fluxes
    ### -----------------------------------
    
    # Afterglow is always here
    unit     = grb.spec_afterglow[0](100*u.GeV).unit
    flx_glow =   np.array([f(E).value for f in grb.spec_afterglow])
    flx_glow = flx_glow*unit

    if grb.prompt:
        # Prompt
        unit = grb.spec_prompt[0](100*u.GeV).unit
        flx_prompt = np.array([f(E) for f in grb.spec_prompt])
        flx_prompt = flx_prompt*unit

        # Prompt + afterglow
        flx_tot = flx_prompt.value + flx_glow[:(grb.id90+1),:].value
        flx_tot = np.concatenate( (flx_tot,flx_glow[grb.id90+1:,:].value) )*unit
    else:        
        flx_tot = flx_glow
        
    max_flx_tot = np.max(flx_tot)

    # Attenuated flux
    unit = grb.spectra[0](100*u.GeV).unit
    flx_tot_att = np.array([f(E).value for f in grb.spectra])
    flx_tot_att = flx_tot_att*unit
    
    ### -----------------------------------
    ### Plot the various fluxes for the E series
    ### -----------------------------------
    with quantity_support():

        for i, energy in enumerate(E):
            # color = cm.rainbow(i/len(E)) 
            # color = cm.cool(i/len(E)) 
            color = cm.Dark2(i/len(E)) 

            ax.plot(grb.tval,flx_tot[:,i],color=color,ls="-", alpha=0.8, lw=1,
                    label="Total")  
            
            if grb.prompt:
                ax.plot(grb.tval,flx_glow[:,i],color= color,alpha=0.7,ls=":",
                        label="afterglow")

                ax.plot(grb.tval[:(grb.id90+1)],flx_prompt[:,i],color=color,ls="--",alpha=0.5,
                            label="prompt")



            ax.plot(grb.tval,flx_tot_att[:,i],color=color,ls="-",marker=".",
                    label=str(round(E[i].value)*E[i].unit))           

        if grb.prompt:
            ax.axvline(grb.t90,color="grey",ls=":",alpha=0.2)
            ax.text(x=grb.t90, y= 1.2*ymin,
                     s= "$t_{90}$="+str(round(grb.t90.value,1)*grb.t90.unit),
                     rotation=0, fontsize=14)
        else:
            ax.axvline(30*u.s,color="black", ls="--", alpha=0.8)
            ax.text(x=30*u.s, y=1.2*ymin, s="30 s",rotation=0, fontsize=14)
            
    ax.set_xlabel("Time (s)")
    title = "{}: {:>2d} E points".format(grb.name,len(grb.Eval))
    ax.set_title(title,fontsize=12)    
    ax.set_ylim(ymin = ymin, ymax = 2*max_flx_tot)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(which="both",ls="-",color="lightgrey",alpha=0.5)
    single_legend(ax,loc='upper left',fontsize=11,bbox_to_anchor=[1.02,0.5])
    
    return

###---------------------------------------------------------------------------
def spectra(grb,opt="None"):
    """
    Plot time evolution (lightcurve) and energy spectra

    """
    if (opt == "None"):
        print("Choose : '2D', 'Packed', 'Time' or 'Energy'")
    elif (opt == "2D" ):
        energy_and_time_2d(grb)
        return
    elif (opt == "Time"):
        time_over_energyband(grb)
    elif (opt == "Energy"):
        energy_over_timeslices(grb)

    return

###---------------------------------------------------------------------------
def energy_over_timeslices(grb):
    """
    Plot E spectra along the time bins.
    Time bin number is variable from 39 to 55
    """

    # Compute grid size
    idxmax = len(grb.tval)
    ncols=6
    nrows = int(idxmax/ncols)+1
    if (idxmax%ncols == 0): nrows=nrows-1

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15,15),
                           sharex=True,
                           sharey=True)

    for icol in range(0,ncols):
        for irow in range(0,nrows):

            # print(irow,icol,icol+ncols*irow)
            idx = icol+ncols*irow
            a = ax[irow][icol]
            if (idx<idxmax):
                a.plot(np.log10(grb.Eval.value),
                       np.log10(grb.fluxval[idx].value),
                       marker='.',
                       markerfacecolor='r',
                       markeredgecolor='r',
                       label="t= {:06.2f}".format(grb.tval[idx]))
                a.grid("both")
                a.legend()
            a.set_xlim(-0.5,4.5)
            a.set_ylim(-22,-4)
            # This has been thought ot be necessary - strange it works without
            #if (irow != nrows-1): a.set_xticklabels([])
            # if (icol != 0):       a.set_yticklabels([])

    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none',
                    top=False, bottom=False,
                    left=False, right=False)
    plt.xlabel("$log(E (GeV))$",size=20)
    plt.ylabel("$log(Flux (Gev^{-1} cm^{-2} s^{-1}) )$",size=20)

    fig.suptitle(grb.name,size=20,y=0.9)
    plt.subplots_adjust(hspace = .001, wspace=0.001)
    plt.show(block=block)

    return
###---------------------------------------------------------------------------
def time_over_energyband(grb):
    """
    Plot t spectra along the E bins - E bin number is fixed :-)
    """

    tb_n1 = grb.t_start["North"]
    tb_n2 = grb.t_stop["North"]
    tb_s1 = grb.t_start["South"]
    tb_s2 = grb.t_stop["South"]

    dt_n = (tb_n2 - tb_n1)/3600
    dt_s = (tb_s2 - tb_s1)/3600
    # print(tbound_n, tbound_s)

    fig, ax = plt.subplots(nrows=8, ncols=5, figsize=(20,20))
    for irow in range(0,8):
        for icol in range(0,5):

            #print(irow,icol,icol+8*irow)
            idx = icol+5*irow
            a = ax[irow][icol]
            label = "E= {:6.2f}".format(grb.Eval[idx])
            #print(idx,grb.Eval[idx])
            # Artificially adds 1 s to avoid log(0)
            a.plot(np.log10(grb.tval.value+1),np.log10(grb.fluxval[:,idx].value),
                   marker='.',markerfacecolor='grey',markeredgecolor='grey',
                   label=label)
            # Display visibility boundaties - Add articifially 1 second  to
            # avoid log(0)
            if (dt_n):
                a.axvline(x=np.log10(tb_n1+1),linestyle =":",color=coln)
                a.axvline(x=np.log10(tb_n2+1),linestyle =":",color=coln)
            if (dt_s):
                a.axvline(x=np.log10(tb_s1+1),linestyle =":",color=cols)
                a.axvline(x=np.log10(tb_s2+1),linestyle =":",color=cols)

            a.set_xlim(-0.5,5.5)
            a.set_ylim(-22,-4)
            if (irow != 7): a.set_xticklabels([])
            if (icol != 0): a.set_yticklabels([])
            a.grid("both")
            a.legend()
            #a.set_xscale("log")
            #a.set_yscale("log")

    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none',
                    top=False, bottom=False,
                    left=False, right=False)
    plt.xlabel("$log(t + 1$ $(s))$",size=20)
    plt.ylabel("$log(Flux$ $(Gev^{-1} cm^{-2} s^{-1}) )$",size=20)

    title = grb.name \
    + " Obs. : North = {:6.2f}h - South = {:6.2f}h".format(dt_n,dt_s)

    fig.suptitle(title,size=16,y=0.9)
    plt.subplots_adjust(hspace = .001, wspace=0.001)
    plt.show(block=block)

    return
###---------------------------------------------------------------------------
def energy_and_time_2d(grb):

    with quantity_support():
        fig, (a1) = plt.subplots(nrows=1,ncols=1,figsize=(12,10))

        h = a1.contourf(grb.tval,grb.Eval,np.log10(grb.fluxval.value.T) )
        a1.set_xscale('log')
        a1.set_yscale('log')
        a1.set_xlabel("Time (s)")
        a1.set_ylabel("E (GeV)")
        cbar = fig.colorbar(h, ax=a1)
        cbar.set_label(str((grb.fluxval[0][0]).unit), rotation=270)

        return
###---------------------------------------------------------------------------
def animated_spectra(grb, emin=0.02 * u.TeV,
                         emax=10 * u.TeV,
                         savefig=False,
                         outdir='./out/'):

    """
    Create a gif animation of time slices,
    """
    #print(" grb_plot.animated_spectra seems corrupted ")
    return

    from matplotlib.animation import FuncAnimation

    def get_model(i):
        return grb.spectra[i]

    def get_time(i):
        return grb.time_interval[i]

    # Initialise plot
    fig_kw = dict(num=grb.name + ' models')
    fig, ax = plt.subplots(**fig_kw)
    model_init = get_model(1)
    fmin, fmax = model_init(energy=emax), model_init(energy=emin)
    x = np.linspace(emin.value, emax.value, 100) * u.TeV
    y = model_init(energy=x)
    ax.set(xlim=(emin.value,emax.value), ylim=(fmin.value, fmax.value))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy [TeV]')
    ax.set_ylabel('Flux [1 / (cm2 s TeV)]')
    ax.set_ylim([1.e-20, 1.e-6])
    ax.grid(which='both')
    line = ax.plot(x,y, color='k', lw=2)[0]

    def animate(i):
        model = get_model(i)
        y = model(energy=x)
        line.set_ydata(y)
        ax.set_title(grb.name + '; z={:.2f}; dt={:6.2f}--{:6.2f} s'
                     .format(grb.z,get_time(i)[0].value,get_time(i)[1].value))
    anim = FuncAnimation(fig, animate, interval=500,
                         frames=len(grb.spectra))
    if savefig == True:
        if not os.path.exists(outdir): os.makedirs(outdir)
        anim.save(outdir + '/' + grb.name + '_animate.gif',
                  writer='imagemagick')
    else:
        plt.draw()

    return
###############################################################################
#
# Plot visibility
#
###############################################################################

###---------------------------------------------------------------------------
def period_plot(twindow,
                ax=None,
                alpha = 0.2, color="black", color2="black",tag="?",
                tshift=0*u.s,
                **kwargs):

    ax = plt.gca() if ax is None else ax

    if (len(twindow[0]) == 0): return ax
    for elt in twindow:
        if elt[0].value !=-9 and elt[1].value != -9:  # Check if -9, i.e. undefined
            t1  = (elt[0] - tshift).datetime
            t2  = (elt[1] - tshift).datetime
            if isinstance(tag, list):
                ax.axvline(t1,color=color,label=tag[0])
                ax.axvline(t2,color=color2,label=tag[1])
            else:
                ax.axvspan(t1,t2,
                           alpha = alpha, color=color,label=tag,**kwargs)

    return ax

###---------------------------------------------------------------------------
def source_alt_and_flux(grb, vis, times, site="None", tshift=0*u.s,
                     ax=None):

    ax = plt.gca() if ax is None else ax

    ### Altitude sampling for plots - absolute time
    altaz   = grb.radec.transform_to(AltAz(obstime  = times,
                                           location = site))

    ### Change time reference if requested
    tref = grb.t_trig - tshift
    times = times - tshift # Change reference

    ### GRB altitude and minimum
    ax.plot(F(times),altaz.alt.value, color="darkblue",alpha=0.5, marker="+",label="Altitude")
    ax.axhline(y =vis.altmin.value,
               ls=":",color="tab:green",label="Min. Alt.")

    ### FLux points
    axx  = ax.twinx()
    Eref = 100 *u.GeV
    iref = np.argmin(np.abs(Eref-grb.Eval))

    axx.plot((grb.t_trig + grb.tval -tshift).datetime,
              #grb.fluxval[:,iref],
              [ g(Eref).value for g in grb.spectra]*grb.spectra[0](Eref).unit,
              marker=".",
              ls = "--",
              lw = 1,
              color="tab:purple",
              label = r"$E \sim {}$".format(Eref))
    axx.set_yscale("log")
    axx.legend(loc='center left', bbox_to_anchor=(1.1, 0.9),fontsize=12)


    ### Trigger (dot and vertical line)
    alttrig = grb.radec.transform_to(AltAz(obstime  = grb.t_trig,
                                           location = site)).alt.value
    ax.plot(F(tref),alttrig,label="Trigger",marker="o",markersize=10,
            color="tab:orange")
    ax.axvline(F(tref),ls=":",color="tab:orange")

    ### Trigger + 1 day (dot and vertical line)
    altend = grb.radec.transform_to(AltAz(obstime  = grb.t_trig + 1*u.day,
                                          location = site)).alt.value
    ax.plot(F(tref+1*u.day),altend,
            label="Trigger + 1 day",marker="o",markersize=10,color="black")
    ax.axvline(F(tref+1*u.day),ls=":",color="grey")

    ax.set_ylim(ymin=0,ymax=1.2*max(altaz.alt.value))
    ax.set_ylabel("altitude (°)")

    return ax

###---------------------------------------------------------------------------
def moon_alt_plot(vis, times, alt, ax=None, alpha=1, color="darkblue"):

    if (ax==None): fig, ax = plt.subplots(figsize=(21,5))

    with quantity_support():
        ax.plot(times.datetime,alt,
                color=color,alpha=alpha,label="Moon altitude")
        ax.axhline(y=vis.moon_maxalt , 
                   color= color, alpha=alpha, ls=":")
    ax.set_ylabel("Alt.(°)")
    ax.legend()

    return ax
###---------------------------------------------------------------------------
def moon_dist_plot(vis, grb,times, radec, site="None", 
                   ax=None, alpha=1, color="purple"):

    if (ax==None): fig, ax = plt.subplots(figsize=(21,5))

    dist = radec.separation(grb.radec)

    with quantity_support():
        ax.plot(times.datetime,dist,
                color=color,alpha=alpha, label="Moon distance")
        ax.axhline(y=vis.moon_mindist , 
                   color= color, ls=":", alpha=alpha, label="Min. distance")
    ax.set_ylabel("Dist.")
    ax.legend()

    return ax

###---------------------------------------------------------------------------
def moonlight_plot(vis, times, ax=None, color="tab:orange"):

    if (ax==None): fig, ax = plt.subplots(figsize=(21,5))

    with quantity_support():
        ax.plot(times.datetime, moon_illumination(times),
                color=color,label="Illumination")
        ax.axhline(y=vis.moon_maxlight , color= color, ls=":",
                   label="Max. illumination")

    ax.set_ylabel("Illumination")
    ax.set_ylim(ymin=0,ymax=1.)
    ax.legend(loc="lower right")

    return ax
###---------------------------------------------------------------------------
def moonphase_plot(times, ax=None, color="red"):

    if (ax==None): fig, ax = plt.subplots(figsize=(21,5))

    with quantity_support():
        ax.plot(times.datetime, moon_phase_angle(times),
                color=color,label="Brightness")
    ax.set_ylabel("Brightness")
    ax.set_ylim(ymin=0,ymax=1.)
    ax.legend(loc="upper right")

    return ax
###---------------------------------------------------------------------------
def visibility_plot(grb,
                    loc  = None,
                    moon_alt  = True,
                    moon_dist = True,
                    ax   = None,
                    dt_before = 0.25*u.day, # Before trigger
                    dt_after  = 0.25*u.day, # After visibility window
                    nalt = 250):
    """
    Plot the night, above-the-horizon and visibility periods, the altitude
    evolution with time as well as a lightcurve for a fixed reference
    energy.

    Parameters
    ----------
    ax : An axis instance
        Plots on this axis
    dt_before : astropy Quantity time, optional
        Time period for plotting before the trigger.
        The default is 0.25*u.day.
    dt_after : astropy Quantity time, optional
        Time period for plotting after the trigger.
        The default is 0.5*u.day.
    nalt : int, optional
        Number of points to sample the altitude evolution with time.
        The default is 25.
    inset : bool, optional
        IfTrue, the figure will be plotted in an inset of the main plot.
        They are simplified and the time is referred to the trigger time.

    """
    # if not vis.vis:
    #     print(" Not visible - plot not shown")
    #     return

    vis = grb.vis[loc]

    if moon_alt:
        if moon_dist:
            ratio = {'height_ratios': [7,2,2]}
            nrows = 3
            ysize = 10
        else:
            ratio = {'height_ratios': [7,2]}
            nrows = 2
            ysize = 9
    else:
        ratio = {'height_ratios': [1]}
        nrows = 1
        ysize = 7        
    
    fig, ax = plt.subplots(nrows=nrows, ncols=1, sharex=True,
                           gridspec_kw=ratio,
                           figsize=(17, ysize))

    with warnings.catch_warnings() and quantity_support() :
        warnings.filterwarnings("ignore")

        duration = (vis.tstop - vis.tstart + dt_after + dt_before).to(u.d)
        dt   = np.linspace(0,duration.value, nalt)
        tobs = vis.tstart  - dt_before + dt*duration.unit

        ### GRB (main plot)
        source_alt_and_flux(grb, vis, tobs, site=vis.site, ax=ax[0])
        ax[0]. set_title(vis.name)
        ax[0].grid("both",ls="--",alpha=0.5)

        ### GRB above horizon
        period_plot(vis.t_event,
                    ax = ax[0],color="tab:blue",alpha=0.2, tag="Above horizon",
                    ymin=0.,  ymax= 0.5,)

        ### Moon if requested
        if moon_alt or moon_dist:
            radec = get_moon(tobs, vis.site) # Moon position along time

            if moon_alt:
                ### Moon veto periods
                period_plot(vis.t_moon_up,
                            ax =ax[1],color="tab:orange",alpha=0.2, tag="Moon veto")
                ax[1].grid("both",ls="--",alpha=0.5)

                ### Moon altitude
                alt = radec.transform_to(AltAz(location=vis.site, obstime=tobs)).alt
                moon_alt_plot(vis, tobs, alt, ax = ax[1], alpha=0.5)

                if moon_dist:
                    ### Moon distance and brigthness
                    moon_dist_plot(vis, grb,tobs, radec,site=vis.site,
                                   ax = ax[2],alpha=0.5) # Distance
                    period_plot(vis.t_moon_up,
                    ax =ax[2],color="tab:orange",alpha=0.2, tag="Moon veto")
                    ax[2].grid("both",ls="--",alpha=0.5)

                    axx = ax[2].twinx()
                    moonlight_plot(vis, tobs,ax=axx) # Brightness
                    #moonphase_plot(tobs,ax=axx) # Phase (correlated to Brightness)

        ### Nights on all plots
        for axis in ax:
            period_plot(vis.t_twilight,
                        ax =axis,color="black",alpha=0.1, tag="Night")

        ### Visibility windows - if the GRB is visible
        if vis.vis:
            for axis in ax:
                period_plot(vis.t_true,
                            ax=axis,color="tab:green",color2="red",
                            tag=["Start","Stop"])

        # Reduce legends on all plots
        import collections
        for axis in ax:
            handles, labels = axis.get_legend_handles_labels()
            by_label = collections.OrderedDict(zip(labels, handles))
            axis.legend(by_label.values(), by_label.keys(),
                      loc='center left', bbox_to_anchor=(1.1, 0.5),fontsize=12)

        axis = ax[nrows-1]

        axis.xaxis.set_major_formatter(mdates.DateFormatter('%d-%H:%M'))
        axis.set_xlabel("Date (DD-HH:MM) UTC")
        axis.tick_params(axis='x', rotation=45)
        axis.set_xlabel("Date")
        axis.grid("both",ls="--",alpha=0.5)
        axis.set_xlim([F(min(tobs)),F(max(tobs))])

    fig.tight_layout(h_pad=0, w_pad=0)
    return
#------------------------------------------------------------------------------
if __name__ == "__main__":

    """
    A standalone function to read a GRB and plot the spectra
    """

    import matplotlib.pyplot as plt
    #plt.style.use('seaborn-talk') # Make the labels readable
    plt.style.use('seaborn-poster') # Make the labels readable - bug with normal x marker !!!
    from SoHAPPy import get_grb_fromfile
    import grb_plot as gplt
    
    dbg        = 0
    ngrb       = 1 # 250
    ifirst     = [980] # ["190829A"]
    grb_folder = "D:/CTAA/SoHAPPy/input/lightcurves/LONG_FITS"
    save_grb   = False # (False) GRB saved to disk -> use grb.py main
    res_dir    = "."
    prompt     = True
    
    # GRB list to be analysed
    if type(ifirst)!=list: grblist = list(range(ifirst, ifirst + ngrb))
    else: grblist = ifirst

    # Loop over GRB list
    for item in grblist:

        grb = get_grb_fromfile(item, grb_folder = grb_folder, 
                                     prompt=prompt,log=None)
        print(grb)
        #gplt.spectra(grb,opt="Packed")
        gplt.energy_spectra(grb)        
        gplt.time_spectra(grb,e_max=1*u.TeV)
        
        # # Read/compute a visibility and plot
        # import sys
        # from pathlib import Path
        # from configuration import Configuration
        # from visibility import Visibility
        # cf = Configuration(sys.argv[1:])
        # for loc in ["North","South"]:
        #     name = Path(cf.vis_dir,grb.name+"_"+loc+"_vis.bin")
        #     grb.vis[loc] = Visibility.read(name)
        #     gplt.visibility_plot(grb, loc =loc)
