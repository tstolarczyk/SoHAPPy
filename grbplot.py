# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from   astropy.coordinates   import SkyCoord, AltAz, EarthLocation
from   astropy.visualization import quantity_support


from gammapy.stats.poisson import (excess_error, background_error)

__all__ = ['spectra','visibility','stats_detection',]

###############################################################################
#
# Plot Energy and time spectra
#
###############################################################################
def spectra(grb,opt="None"):
    """

    Returns
    -------
    None.

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
    # print(idxmax," ->",nrows,"*", ncols," = ",nrows*ncols)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20,20),sharex=True, sharey=True)

    for icol in range(0,ncols):
        for irow in range(0,nrows):

            # print(irow,icol,icol+ncols*irow)
            idx = icol+ncols*irow
            a = ax[irow][icol]
            if (idx<idxmax):
                a.plot(np.log10(grb.Eval.value),np.log10(grb.fluxval[idx].value),
                       marker='.',markerfacecolor='r',markeredgecolor='r',
                       label="t= {:06.2f}".format(grb.tval[idx]))
                a.grid("both")
                a.legend()
            a.set_xlim(--0.5,4.5)
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
    plt.show()

    return
###############################################################################
def time_over_energyband(grb):
    """
    Plot t spectra along the E bins - E bin number is fixed :-)
    """

    tb_n = [ (grb.t_true["North"][0] - grb.t_trig).datetime.total_seconds() ,
             (grb.t_true["North"][1] - grb.t_trig).datetime.total_seconds() ]
    tb_s = [ (grb.t_true["South"][0] - grb.t_trig).datetime.total_seconds(),
             (grb.t_true["South"][1] - grb.t_trig).datetime.total_seconds() ]
    dt_n = (tb_n[1] - tb_n[0])/3600
    dt_s = (tb_s[1] - tb_s[0])/3600
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
                a.axvline(x=np.log10(tb_n[0]+1),linestyle =":",color="blue")
                a.axvline(x=np.log10(tb_n[1]+1),linestyle =":",color="blue")
            if (dt_s):
                a.axvline(x=np.log10(tb_s[0]+1),linestyle =":",color="red")
                a.axvline(x=np.log10(tb_s[1]+1),linestyle =":",color="red")

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

    title = grb.name
    + " Obs. : North = {:6.2f}h  -  South = {:6.2f}h".format(dt_n,dt_s)

    fig.suptitle(title,size=16,y=0.9)
    plt.subplots_adjust(hspace = .001, wspace=0.001)
    plt.show()

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

###############################################################################
def energy_and_time_packed(grb):

    with quantity_support():
        fig, (a2,a3) = plt.subplots(nrows=2,ncols=1,figsize=(12,15))

        for i,t in enumerate(grb.t_s):
            # Unabsorbed  spectra, not very exciting
            # a2.plot(grb.Eval,grb.fluxval[i,:],label="t="+str(t),ls='--')
            grb.spectra[i].plot([min(grb.Eval),max(grb.Eval)],a2,
                                       label="t="+str(t))
            a2.set_xscale("log")
            a2.set_yscale("log")
            a2.set_ylim(ymin=1e-16)
            a2.set_xlabel("Energy "+ str(grb.Eval[0].unit))
            a2.set_title(grb.name + " "+str(len(grb.tval)) + " t slices")

            #a2.set_ylabel("Flux (/GeV/cm2/s)")
            if (len(grb.tval)<15): a2.legend() # Too many t slices

        for i in range(0,len(grb.Eval)):
            flux = [f[i].value for f in grb.fluxval][:-1] * grb.fluxval[0].unit
            a3.plot(grb.t_s,flux,
                    label="E="+str(round(grb.Eval[i].value,2)))
            a3.set_xscale("log")
            a3.set_yscale("log")
            a3.set_xlabel("Time (s)")
            a3.set_title(grb.name
                         + " "+str(len(grb.Eval)) + " E bands")
        t1 = (grb.t_true["North"][0]-grb.t_trig).sec*u.s
        t2 = (grb.t_true["North"][1]-grb.t_trig).sec*u.s
        a3.axvline(x=t1,ls=":",color="blue")
        a3.axvline(x=t2,ls=":",color="blue")

        t1 = (grb.t_true["South"][0]-grb.t_trig).sec*u.s
        t2 = (grb.t_true["South"][1]-grb.t_trig).sec*u.s
        a3.axvline(x=t1,ls=":",color="red")
        a3.axvline(x=t2,ls=":",color="red")

    plt.show()

    return
###############################################################################
#
# Plot visibility
#
###############################################################################
def visibility(grb,opt="None"):

    site  = ["North","South"]

    print("=============================================================")
    print(" VISIBILITY ")
    print("Minimum altitude for detection = ",grb.altmin)

    for loc in site:
        print(loc," -------------------------------------------------------------")
        print('  Visible    : {} - tonight, prompt : {}, {}'
              .format(grb.vis[loc],
                      grb.vis_tonight[loc],
                      grb.vis_prompt[loc]))
        if (grb.vis_tonight[loc]): # Seen within 24hr after the trigger
            print(" Event {}    :  {} * {}".format(loc,
                  grb.t_event[loc][0].datetime,
                  grb.t_event[loc][1].datetime))

            print(" TWilight {} :  {} * {}".format(loc,
                  grb.t_twilight[loc][0].datetime,
                  grb.t_twilight[loc][1].datetime))

            print(" True {}     :  {} * {}".format(loc,
                  grb.t_true[loc][0].datetime,
                  grb.t_true[loc][1].datetime))

    if (opt == "plot"):
        visibility_plot(grb)
    else:
        print(" Use opt='plot' for plots")

    print("===========================================================================\n")

    return
###########################################################################
def visibility_plot(grb):
     """
     Plot visibility from the souce position and times around which
     the GRB is indeed visible.
     The time scale extend from the trigger time to 1 hr after the last
     visibility point

     """

     site  = ["North","South"]
     color = {"North":"blue","South":"red"}
     tmin = grb.t_trig # GRB trigger time
     tmax = max(grb.t_event["South"][1],grb.t_event["North"][1])
     duration = tmax - tmin
     dt = np.linspace(0,duration.value,100)
     t =  tmin + dt

     with quantity_support():

         plt.xticks(rotation=70)

         for loc in site:
             where = EarthLocation.of_site(grb.site[loc])

             if (grb.vis_tonight[loc]):
                 altaz =  grb.radec.transform_to( AltAz(obstime=t,
                                                  location=where ))
                 plt.plot(t.datetime,altaz.alt,
                          linestyle=":",
                          color=color[loc])

                 # For the given site, above 10°
                 tstart = grb.t_true[loc][0]
                 tstop  = grb.t_true[loc][1]
                 duration = tstop - tstart
                 t10 = tstart +  np.linspace(0,duration.value,100)
                 altaz =  grb.radec.transform_to( AltAz(obstime=t10,
                                                     location=where ))
                 plt.plot(t10.datetime,altaz.alt, color=color[loc],label=loc)

             else:
                 print("Not visible in :",loc," within 24h after trigger")

         if (grb.vis_tonight[site[0]] or grb.vis_tonight[site[1]]):
             plt.hlines(y=10*u.deg,xmin=tmin.datetime,xmax=tmax.datetime,
                        linestyles='--',colors="grey",label="Min alt.")
             plt.title(grb.name)
             plt.xlabel("Date (MM-DD-HH) UTC")
             plt.ylim(ymin=0*u.deg,ymax=90*u.deg)
             plt.legend()
             plt.show()


     return


###############################################################################
#
# Statistics (obsolete)
#
###############################################################################
def stats_detection(grb, simulations, ax_excess=None, ax_relative=None, ax_bkg=None, ax_sigma=None, savefig=False, outdir='./out/'):
     """Plot some detection statistics"""

     stats = grb.get_cumulative_stats(simulations)

     import matplotlib.pyplot as plt

     plt.figure(num=grb.name + ' stats', figsize=(18, 6)) # Figure frame

     if ax_excess is None:
         ax_excess = plt.subplot2grid((1, 4), (0, 0))
     yerr = excess_error(
         stats['n_on'],
         stats['n_off'],
         stats['alpha']
     )
     ax_excess.errorbar(stats['livetime'], stats['excess'],
                        yerr=yerr,
                        color='black', fmt='o')
     ax_excess.set_xlabel('Livetime [s]', fontweight='bold')
     ax_excess.set_ylabel('#Evts', fontweight='bold')
     ax_excess.set_title('Cumulated excess', fontweight='bold')
     ax_excess.set_ylim(0., (stats['excess'] + yerr).max() * (1.1))
     ax_excess.grid(which='both')
     ax_excess.set_xscale('log')

     if ax_relative is None:
         ax_relative = plt.subplot2grid((1, 4), (0, 1))
     yerr = excess_error(
         stats['n_on'],
         stats['n_off'],
         stats['alpha']
     )
     ax_relative.errorbar(stats['livetime'], yerr/stats['excess'],
                        yerr=0,
                        color='black', fmt='o')
     ax_relative.set_xlabel('Livetime [s]', fontweight='bold')
     ax_relative.set_ylabel('Excess relative error', fontweight='bold')
     ax_relative.set_title('Relative error evolution', fontweight='bold')
     ax_relative.grid(which='both')
     ax_relative.set_xscale('log')

     if ax_bkg is None:
         ax_bkg = plt.subplot2grid((1, 4), (0, 2))
     yerr = background_error(
         stats['n_off'],
         stats['alpha']
     )
     ax_bkg.errorbar(stats['livetime'], stats['bkg'],
                     yerr=background_error(
                         stats['n_off'],
                         stats['alpha']
                     ),
                     color='black', fmt='o')

     ax_bkg.set_xlabel('Livetime [s]', fontweight='bold')
     ax_bkg.set_ylabel('#Evts', fontweight='bold')
     ax_bkg.set_title('Cumulated background', fontweight='bold')
     ax_bkg.set_ylim(stats['bkg'].min() * 0.9, (stats['bkg'] + yerr).max() * (1.1))
     ax_bkg.grid(which='both')
     ax_bkg.set_xscale('log')
     ax_bkg.set_yscale('log')

     if ax_sigma is None:
         ax_sigma = plt.subplot2grid((1, 4), (0, 3))
     ax_sigma.errorbar(stats['livetime'], stats['sigma'],
                       color='black', fmt='o')
     ax_sigma.set_xlabel('Livetime [s]', fontweight='bold')
     ax_sigma.set_ylabel('Significance', fontweight='bold')
     ax_sigma.set_title('Significance (Li & Ma)', fontweight='bold')
     ax_sigma.set_ylim(0., stats['sigma'].max() * (1.1))
     ax_sigma.grid(which='both')
     ax_sigma.set_xscale('log') # CHANGED

     plt.tight_layout()
     if savefig == True:
         if not os.path.exists(outdir):
             os.makedirs(outdir)
         plt.savefig(outdir + '/' + grb.name + '.png')
     return ax_excess, ax_relative, ax_bkg, ax_sigma


############################################################################
def make_gif_from_models(grb, emin=0.02 * u.TeV,
                         emax=10 * u.TeV,
                         savefig=False,
                         outdir='./out/'):

    """
    Create a gif animation of time slices
    """

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
        ax.set_title(grb.name + '; z={:.2f}; dt={}--{} s'.format(grb.z,
                                                            get_time(i)[0].value,
                                                            get_time(i)[1].value))
    anim = FuncAnimation(fig, animate, interval=500,
                         frames=len(grb.spectra))
    plt.draw()
    if savefig == True:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        anim.save(outdir + '/' + grb.name + '_animate.gif', writer='imagemagick')