# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import os
import matplotlib.pyplot as plt
import numpy as np
<<<<<<< HEAD
import warnings

import astropy.units as u
from   astropy.coordinates   import AltAz
from astropy.time import Time
from   astropy.visualization import quantity_support
import matplotlib.dates as mdates

=======
import astropy.units as u
from   astropy.coordinates   import AltAz, EarthLocation
from   astropy.visualization import quantity_support
import matplotlib.dates as mdates

import ana_config as cf

>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
# Bigger texts and labels
# plt.style.use('seaborn-talk') 
plt.style.use('seaborn-poster') # Bug with normal x marker !!!

# If False, avoid scirpt to be paused when a plot is popped-up (plt.show)
block = False
#-----------------------------------------------------------------------------#
n_t_2disp = 8 # Number of curves to show in lighcurve plots
n_E_2disp = 5 # Number of curves to show in energy spectra
coln = "blue" # color for North site plots
cols = "red"  # color for South site plots
#-----------------------------------------------------------------------------#

__all__ = ['energy_and_time_packed',
           'spectra',
           'energy_over_timeslices',
           'time_over_energyband',
           'energy_and_time_2d',
           'animated_spectra',
<<<<<<< HEAD
           'visibility_plot'
=======
           'visibility_chart',
           'visibility_alt_plot'
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
           ]

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
def energy_and_time_packed(grb):

    with quantity_support():
        
        ### Energy spectra for various measurement points in time
        # Note: the plot methods handle the labels, should not be overwritten        
<<<<<<< HEAD
        nspectra = len(grb.tval)-1
        if (nspectra > n_E_2disp): 
            dnt = int(round(nspectra/n_E_2disp))
            tlist = list(range(0,nspectra,dnt))
            #tlist = [12, 13, 14, 15] + tlist # Typical Prompt times
=======
        nspectra = len(grb.t_s)
        if (nspectra > n_E_2disp): 
            dnt = int(round(nspectra/n_E_2disp))
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
        else:
            dnt=1
        
        ymin = 1e-16
        
        fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,figsize=(10,12))

<<<<<<< HEAD
        for i in tlist:
            #print(" energy_and_time_packed",i)
            t = grb.tval[i]
=======
        for i in range(0,nspectra,dnt):
            t = grb.t_s[i]
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
            grb.spectra[i].plot([min(grb.Eval),max(grb.Eval)],
                                ax1,
                                label="t={:>8.1f}".format(t))
            c = ax1.lines[-1].get_color() # Last color
            ax1.plot(grb.Eval.to(u.TeV),
                     grb.fluxval[i,:].to(1/u.cm**2/u.s/u.TeV),
                     ls='--',
                     lw=1.,
                     marker=".",
                     alpha=0.5,
                     color=c)
            
        for i in range(0,nspectra):
<<<<<<< HEAD
            t = grb.tval[i]
=======
            t = grb.t_s[i]
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
            grb.spectra[i].plot([min(grb.Eval),max(grb.Eval)],
                                ax1,
                                alpha=0.2,
                                color="grey",
                                lw=1)
            
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylim(ymin=ymin)
        
        ax1.axvline(10*u.GeV, color="grey",alpha=0.2)
        ax1.text(x=10*u.GeV,
                y=ymin*50,
                s="10 GeV",
                rotation=270,
                fontsize=14,va="bottom")
<<<<<<< HEAD
        title = "{}: {:>2d} Flux points".format(grb.name,len(grb.tval))
=======
        title = "{}: {:>2d} meas. points".format(grb.name,len(grb.t_s))
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
        ax1.set_title(title,fontsize=12)
        if (n_E_2disp<=15): ax1.legend(fontsize=12) # Too many t slices

        
        # Light curve for some energy bins
        nlightcurves = len(grb.Eval)
        if (nlightcurves> n_t_2disp):
            dnE = int(round(nlightcurves/ n_t_2disp))
        else:
            dnE=1
            
        
        for i in range(0,nlightcurves,dnE):
<<<<<<< HEAD
            flux = [f[i].value for f in grb.fluxval ]* grb.fluxval[0].unit
            ax2.plot(grb.tval,
=======
            flux = [f[i].value for f in grb.fluxval][:-1] * grb.fluxval[0].unit
            ax2.plot(grb.t_s,
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
                     flux,
                     marker=".",
                     label="E= {:>8.2f}".format(grb.Eval[i]))
        
        ax2.axvline(30*u.s,color="grey",alpha=0.2)
        ax2.text(x=30*u.s,
                 y=(flux[len(flux)-1]), 
                 s="30 s",
                 rotation=0,
                 fontsize=14)           
        
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlabel("Time (s)")
        title = "{}: {:>2d} E points".format(grb.name,len(grb.Eval))
        ax2.set_title(title,fontsize=12)
            
        if (n_t_2disp<=8): 
            ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)
<<<<<<< HEAD
=======

        # Set visibility limits
        t1 = (grb.t_true["North"][0]-grb.t_trig).sec*u.s
        t2 = (grb.t_true["North"][1]-grb.t_trig).sec*u.s
        ax2.axvspan(t1,t2, color=coln,alpha=0.2,)

        t1 = (grb.t_true["South"][0]-grb.t_trig).sec*u.s
        t2 = (grb.t_true["South"][1]-grb.t_trig).sec*u.s
        ax2.axvspan(t1,t2, color=cols,alpha=0.2)
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
        
        plt.tight_layout()
        plt.show(block=block)

    return
###############################################################################
def spectra(grb,opt="None"):
    """
    Plot time evolution (lightcurve) and energy spectra

    """
    if (opt == "None"):
        print("Choose : '2D', 'Packed', 'Time' or 'Energy'")
    elif (opt == "2D" ):
        energy_and_time_2d(grb)
    elif (opt == "Packed" ):
        energy_and_time_packed(grb)
        return
    elif (opt == "Time"):
        time_over_energyband(grb)
    elif (opt == "Energy"):
        energy_over_timeslices(grb)

    return

###############################################################################
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
###############################################################################
def time_over_energyband(grb):
    """
    Plot t spectra along the E bins - E bin number is fixed :-)
    """

<<<<<<< HEAD
    tb_n1 = grb.t_start["North"]
    tb_n2 = grb.t_stop["North"]
    tb_s1 = grb.t_start["South"]
    tb_s2 = grb.t_stop["South"]

    dt_n = (tb_n2 - tb_n1)/3600
    dt_s = (tb_s2 - tb_s1)/3600
=======
    tb_n = [ (grb.t_true["North"][0] - grb.t_trig).datetime.total_seconds() ,
             (grb.t_true["North"][1] - grb.t_trig).datetime.total_seconds() ]
    tb_s = [ (grb.t_true["South"][0] - grb.t_trig).datetime.total_seconds(),
             (grb.t_true["South"][1] - grb.t_trig).datetime.total_seconds() ]
    dt_n = (tb_n[1] - tb_n[0])/3600
    dt_s = (tb_s[1] - tb_s[0])/3600
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
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
<<<<<<< HEAD
                a.axvline(x=np.log10(tb_n1+1),linestyle =":",color=coln)
                a.axvline(x=np.log10(tb_n2+1),linestyle =":",color=coln)
            if (dt_s):
                a.axvline(x=np.log10(tb_s1+1),linestyle =":",color=cols)
                a.axvline(x=np.log10(tb_s2+1),linestyle =":",color=cols)
=======
                a.axvline(x=np.log10(tb_n[0]+1),linestyle =":",color=coln)
                a.axvline(x=np.log10(tb_n[1]+1),linestyle =":",color=coln)
            if (dt_s):
                a.axvline(x=np.log10(tb_s[0]+1),linestyle =":",color=cols)
                a.axvline(x=np.log10(tb_s[1]+1),linestyle =":",color=cols)
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76

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
###############################################################################
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
############################################################################
def animated_spectra(grb, emin=0.02 * u.TeV,
                         emax=10 * u.TeV,
                         savefig=False,
                         outdir='./out/'):

    """
    Create a gif animation of time slices,
<<<<<<< HEAD
    """
    print(" grb_plot.animated_spectra seems corrupted ")
    return
=======
    
    """
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76

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
<<<<<<< HEAD
def visibility_plot(grb, 
                    ax=None, loc=None,
                    depth = 1,
                    dt_before   = 0.25*u.day, 
                    dt_after    = 1.25*u.day, 
                    nalt = 25, 
                    inset= False):
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
    if (loc == None):
        print(" visibility_plot : A location should be defined")
        return
    
    if (ax == None):
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,6))
        
    ###-------------------------
    def F(x): 
        if (inset): return x.sec
        return x.datetime
    ###-------------------------
    
    if (inset): # To be corrected
        tshift = grb.t_trig
        tref  = grb.t_true[0][0]
    else: 
        tshift = 0*u.s
        tref = grb.t_trig - tshift
    
    ### Plot limits
    tmin    = grb.t_trig - dt_before -tshift
    tmax    = grb.t_trig + dt_after - tshift        
    ax.set_xlim([F(tmin),F(tmax)])

    ### Altitude
    # Altaz sampling for plots - absolute time
    dt = np.linspace(0,(dt_before+dt_after).value, nalt)
    t = grb.t_trig  - dt_before + dt*dt_before.unit
    where = grb.pos_site[loc]
    altaz   = grb.radec.transform_to(AltAz(obstime  = t, 
                                                location = where))
    t = t - tshift # Change reference
    
    ax.plot(F(t),altaz.alt.value,
            ls=":", color="darkblue",alpha=0.5, marker="+",label="Altitude")
    ax.axhline(y =grb.altmin.value,
              ls=":",color="tab:green",label="Min. Alt.")

    # Trigger
    alttrig = grb.radec.transform_to(AltAz(obstime  = grb.t_trig, 
                                                location = where)).alt.value
    ax.plot(F(tref),alttrig,label="Trigger",marker="o",markersize=10,color="tab:orange")
    ax.axvline(F(tref),ls=":",color="tab:orange")
    
    altend = grb.radec.transform_to(AltAz(obstime  = grb.t_trig
                                               + depth*u.day, 
                                               location = where)).alt.value
    ax.plot(F(tref+depth*u.day),altend,
            label="End of search",marker="o",markersize=10,color="black")
    ax.axvline(F(tref+depth*u.day),ls=":",color="black")
    
    ### FLux points
    if (not inset):
        axx = ax.twinx()
        Eref = 100 *u.GeV
        iref = np.argmin(np.abs(Eref-grb.Eval))
        
        axx.plot((grb.t_trig + grb.tval -tshift).datetime,
                  grb.fluxval[:,iref],
                  marker=".",
                  ls = "--",
                  lw = 1,
                  color="tab:purple",
                  label = r"$E \sim {}$".format(Eref))
        axx.set_yscale("log")
        axx.legend(loc='center left', bbox_to_anchor=(1.1, 0.9),fontsize=12)
    
    ### Nights
    first=True
    for elt in grb.t_twilight[loc]:
        if isinstance(elt[0],Time): 
            t_dusk  = elt[0] - tshift
        else:
            t_dusk = tmin
        if isinstance(elt[1],Time):
            t_dawn  = elt[1] - tshift
        else:
            t_dawn = tmax
        if (first):
            label="Night"
            first = False
        else: label= None
        ax.axvspan(F(t_dusk),F(t_dawn), alpha=0.2,color="black",label=label)  
    
    ### Above horizon periods
    first=True
    for elt in grb.t_event[loc]:
        if isinstance(elt[0],Time): 
            t_rise = elt[0] - tshift
        else: 
            t_rise = tmin
        if isinstance(elt[1],Time):
            t_set   = elt[1] -tshift
        else:
            t_set = tmax
            
        if (first):
            label="Above horizon"
            first = False
        else: label= None
        ax.axvspan(xmin=F(t_rise), xmax=F(t_set),
                   ymin=0.,     ymax= 0.5,
                   alpha=0.2,color="tab:blue",label=label)  
    
    ### Visibility windows
    if (grb.vis_tonight[loc]):
        first = True
        for elt in grb.t_true[loc]:
            t_start = elt[0] -tshift
            t_stop  = elt[1] -tshift
            if (first):
                label1="Start"
                label2= "Stop"
                first = False
            else: label1= label2 = None
            ax.axvline(F(t_start),label=label1,color="tab:green")
            ax.axvline(F(t_stop),label=label2,color="tab:red")  
    
    if (inset==False):
        ax. set_title(grb.name + " -" + loc)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%H:%M'))
        ax.set_xlabel("Date (DD-HH:MM) UTC")
        ax.tick_params(axis='x', rotation=70)
        ax.set_ylabel("altitude (°)")
        ax.set_ylim(ymin=0,ymax=1.2*max(altaz.alt.value))
        ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5),fontsize=12)

    else:
        ax.set_ylim(ymax=1.05*alttrig)
        ax.set_xlabel("Time (s) wrt trigger")

    return     
=======
def visibility_chart(grb):
    """
    A schematic display of the various visibility intervals
    """
    
    fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1,figsize=(8,7))

    site  = ["North","South"]
    for loc, ax  in zip(site,[ax1,ax2]):     
#        if (grb.vis_tonight[loc]): # Seen within 24hr after the trigger
        if (grb.vis_tonight[loc]): # Seen within 24hr after the trigger
            t_event1 = grb.t_event[loc][0]
            t_event2 = grb.t_event[loc][1]
            t_twil1  = grb.t_twilight[loc][0]
            t_twil2  = grb.t_twilight[loc][1]
            t_true1  = grb.t_true[loc][0]
            t_true2  = grb.t_true[loc][1]
            
            tmin = min(t_event1,t_twil1,t_true1)
            tmax = max(t_event2+ 1*u.d*cf.day_after,
                       t_twil2+ 1*u.d*cf.day_after,
                       t_true2+ 1*u.d*cf.day_after)

            
           # Event : above horizon
            ax.plot([t_event1.datetime,t_event2.datetime],
                    ["Event","Event"],
                    color="tab:blue",
                    label=loc+"\nAbove horizon")
            ax.axvline(t_event1.datetime,ls=":",lw=1,color="tab:blue")
            ax.axvline(t_event2.datetime,ls=":",lw=1,color="tab:blue")
            
            # Twilight : dark time             
            ax.plot([t_twil1.datetime,t_twil2.datetime],
                    ["Twilight","Twilight"],
                    color="black",
                    label="Dark time")
            ax.axvline(t_twil1.datetime,ls=":",lw=1,color="black")
            ax.axvline(t_twil2.datetime,ls=":",lw=1,color="black")
            
            # True : observable (above horizon, dark time, triggered )
            ax.plot([t_true1.datetime,t_true2.datetime],
                    ["True","True"],
                    color="tab:green",
                    label="Observable")
            ax.axvline(t_true1.datetime,ls=":",lw=1,color="tab:green")
            ax.axvline(t_true2.datetime,ls=":",lw=1,color="tab:green")

            # Trigger
            ax.axvline(x=grb.t_trig.datetime,
                       ls=":",color="tab:green",label="Trigger")                 

            if (cf.day_after):
               
                # Event : above horizon
                ax.plot([(t_event1+1*u.d).datetime,(t_event2+1*u.d).datetime],
                        ["Event","Event"],
                        color="tab:blue")   
                ax.axvline((t_event1 + 1*u.d).datetime,
                           ls=":",lw=1,color="tab:blue")
                ax.axvline((t_event2 + 1*u.d).datetime,
                           ls=":",lw=1,color="tab:blue")
                
                # Twilight : dark time             
                ax.plot([(t_twil1+1*u.d).datetime,(t_twil2+1*u.d).datetime],
                        ["Twilight","Twilight"],
                        color="black")    
                ax.axvline((t_twil1 + 1*u.d).datetime,
                           ls=":",lw=1,color="black")
                ax.axvline((t_twil2 + 1*u.d).datetime,
                           ls=":",lw=1,color="black")

                # True : observable (above horizon, dark time, triggered )
                tstart = max(t_twil1+1*u.d,t_event1+1*u.d)
                ax.plot([tstart.datetime,(t_true2+1*u.d).datetime],
                        ["True","True"],
                        color="tab:green")
                ax.axvline(tstart.datetime,
                           ls=":",lw=1,color="tab:green")
                ax.axvline((t_true2 + 1*u.d).datetime,
                           ls=":",lw=1,color="tab:green")                
            
            
            ax.tick_params(axis='x', rotation=70)
            ax.tick_params(axis='y', rotation=45)
        
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            ax.set_xlim(xmin=tmin.datetime, xmax=tmax.datetime)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)

    plt.tight_layout()
    plt.show(block=block)
    return

###########################################################################
def visibility_alt_plot(grb):
    
    """
    Plot visibility from the souce position and times around which
    the GRB is indeed visible.
    The time scale extend from the trigger time to 1 hr after the last
    visibility point.
    Optionnaly, display the same for the next day in an approximate way 
    (adding one day to dates)

    """
    site  = ["North","South"]
    color = {"North":coln,"South":cols}
    
    if (grb.vis_tonight[site[0]]==False 
        and grb.vis_tonight[site[1]]==False):
        return
   
    tgrb_max = grb.t_trig + grb.tval[-1] # last point in GRB lightcurve

    # Compute time limits from the GRB visibility - add one day if requested
    if grb.vis_tonight["North"]:
        tmin = grb.t_event["North"][0]
        tmax = grb.t_event["North"][1] + 1*u.d*cf.day_after
        if grb.vis_tonight["South"]:
            tmin = min(tmin,grb.t_event["South"][0])
            tmax = max(tmax,grb.t_event["South"][1]+ 1*u.d*cf.day_after) 
            # tmax = max(tmax,grb.t_event["South"][1]+ 1*u.d*cf.day_after,tgrb_max) 
    else:
        tmin = grb.t_event["South"][0]
        tmax = grb.t_event["South"][1] + 1*u.d*cf.day_after    
            
    dt_visible = tmax - tmin
    dt = np.linspace(0,dt_visible.value,100)
    t =  (tmin + dt)

    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,5))
    with quantity_support():
        
        ax.tick_params(axis='x', rotation=70)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%D:%H:%M'))
        ax.set_ylim(ymin=0*u.deg,ymax=90*u.deg)
                    
        for loc in site:
            where = EarthLocation.of_site(grb.site[loc])

            if (grb.vis_tonight[loc]):
                
                # GRB trajectory in the sky (altitude)
                altaz =  grb.radec.transform_to( AltAz(obstime=t,
                                                  location=where ))
                ax.plot(t.datetime,altaz.alt,ls=":",color=color[loc])

                # Show the night(s)
                t1 = grb.t_twilight[loc][0].datetime
                t2 = grb.t_twilight[loc][1].datetime
                ax.axvspan(t1,t2, color=color[loc],alpha=0.2,)   
                if (cf.day_after):
                    t1 = (grb.t_twilight[loc][0]+ 1*u.d).datetime
                    t2 = (grb.t_twilight[loc][1]+ 1*u.d).datetime
                    ax.axvspan(t1,t2, color=color[loc],alpha=0.2,)   
                        
                tstart = grb.t_true[loc][0]
                tstop  = grb.t_true[loc][1]
                duration = tstop - tstart
                t10 = tstart +  np.linspace(0,duration.value,100)
                altaz =  grb.radec.transform_to( AltAz(obstime=t10,
                                                    location=where ))
                ax.plot(t10.datetime,altaz.alt, color=color[loc],label=loc)
                
                if (cf.day_after):
                    # Trigger occured, recompute the min time for that day
                    tstart = max(grb.t_twilight[loc][0] + 1*u.d,
                                 grb.t_event[loc][0] + 1*u.d)
                    tstop  = grb.t_true[loc][1] + 1*u.d
                    duration = tstop - tstart
                    t10 = tstart +  np.linspace(0,duration.value,100)
                    altaz =  grb.radec.transform_to( AltAz(obstime=t10,
                                                        location=where ))
                    idlow = np.where(altaz.alt < grb.altmin)[0]
                    if (idlow.size ==0):idlow = t10.size - 1
                    else : idlow = idlow[0]
                    ax.plot(t10.datetime[:idlow],altaz.alt[:idlow], 
                            color=color[loc])              

            else:
                print("Not visible in :",loc," within 24h after trigger")

 
        ax.axhline(10*u.deg,ls ="--",color="grey",label="Min alt.")
        ax.axvline(x=grb.t_trig.datetime,color="grey",ls=":")
        
        axx = ax.twinx()
        # Plot the GRB spectrum as an illustration
        Eref = 100 *u.GeV
        dE = np.inf
        for i,E in enumerate(grb.Eval):
            if np.abs(E - Eref) < dE:
                dE = np.abs(E - Eref)
                iref = i
                
        axx.plot((grb.t_trig + grb.tval).datetime,
                 grb.fluxval[:,iref],
                 marker=".",
                 ls = "--",
                 lw = 1,
                 color="tab:green",
                 label = "$E \sim {}$".format(Eref)) 
        axx.set_yscale("log")
        axx.legend()
        
        ax.set_title(grb.name)
        ax.set_xlabel("Date (DD-HH-MM) UTC")
        ax.set_xlim(xmax=(tmax + 4*u.h).datetime)
        ax.legend()
        
        plt.show(block=block)
       
       
    

    return
>>>>>>> 735c550f8a1ddaf3003080a0cc90781d95a86a76
