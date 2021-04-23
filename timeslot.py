# -*- coding: utf-8 -*-
"""
This code handle the observation windows in CTA North and South sites and
combine observation of both site and/or on several days

Created on Mon Mar 30 10:36:09 2020

@author: Stolar
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from copy import deepcopy
import warnings

from obs_slice import Slice

vis_keys = ["North","South","Both"]

###----------------------------------------------------------------------------
class Slot():
    """
    A set of time slices, called "slot" is related to a GRB
    lightcurve given as a set of measurement points, and a visibility,
    namely an observation site.
    Slots are essentially a list of Slices (see the Slice class).
    The slice number starts at zero, and a a slice "i" covers the GRB
    from measurement i to measurement i+1.

    Time slots can be the results of merging initial slots, allowing for
    managing overlapping or consecutive observations.

    Slots get a name for clarity, and can be associated to a site, or
    a combination of sites explcitely. They follow the site naming convention
    of slices, which is mandatory for slot merging.

    The normal process for building time slots is the following:
        1. Create a Slot with 'naked' slices,
        with only the time information from the GRB. This requires
        deciding how the observation point is fixed.

        2. Apply a visibility constraint from a given site,
        change the slice collection, tag the slices with the site.

        3. If requested, create more slots for more visibilities
        (e.g. the following days for a given site)

        4. Merge all above cases in a single Slot.

        5. In the end, "dress" the slot with the physical information
        (flux, energy and sky position boundaries).

    The slot is then ready to be sent to the simulation Class.

    """
    #------------------------------------------------------------
    def __init__(self,grb, name="naked",
                 site="?",delay=-1*u.s,opt="end",debug=False):
        """
        Create a naked observation slot consisting of a GRB, a visibility
        window, naked slices, and a flag mentionning if the slot was dressed
        with physics information.
        The time reference is the GRB trigger.

        Parameters
        ----------
        grb :  GRB class
            a GRB class instantiation
        site : string, optional
            The site to which this slot is associated
            Can be North, South or Both
        name : String, optional
            A name for the slot. The default is "naked".
        opt : string, optional
            Allows choosing how the observation points are chosen.
            The default is "end".
        Returns
        -------
        None.
        """

        self.grb   = grb
        self.name  = grb.name
        self.site  = site
        self.delay = delay # Has so far no sense if merged sites
        self.opt   = opt
        self.tstart= np.NINF*u.h # Negative infinite
        self.tstop = np.Inf*u.h  # Positive infinite
        self.phys  = False       # Initially not associated to flux values
        self.dbg   = debug
        slices = []
        for i in range(1,len(grb.tval)): # Starts at 1, but uses i-1
            t1 = grb.tval[i-1]
            t2 = grb.tval[i]
            s  = Slice(i-1,t1,t2)
            slices.append(s)

        self.slices   = np.asarray(slices) # Works also it the list is already an array
        return

    #------------------------------------------------------------
    def dress(self, name="dressed", irf_dir="./", arrays=None, debug=False):
        """
        Dress slices with physics information

        Returns
        -------
        None.

        """

        self.phys = True
        for s in self.slices: s.dress(self.grb,
                                      irf_dir = irf_dir,
                                      arrays  = arrays,
                                      opt=self.opt)

        # After dressing, two consecutives slices can be associated to the
        # same spectrum id and should be merged
        to_be_compacted = True
        while (to_be_compacted):

            fidlist = [s.fid() for s in self.slices]
            if (len(fidlist) != len(set(fidlist))):
                #print(" Go to compacting\n", self)
                self.compact()
                #print("compacted!!!")
            else:
                to_be_compacted = False

        return

    #------------------------------------------------------------
    def compact(self):

        prev_fid = -999
        newslices = []
        idt = 0
        for i, s in enumerate(self.slices):
            fid = s.fid()
            if (fid == prev_fid):
                # Change end-time of prev-spec to this spec
                # print(" Slice ",i," merged with slice ",i-1)
                self.slices[i-1].merge(s)
            else:
                s.set_id(idt)
                newslices.append(s)
                idt+=1
            prev_fid = fid

        self.slices = newslices

        return

    #------------------------------------------------------------
    def copy(self,name="a_copy"):

        # Could change name and site here
        self.name = name
        slot_copy = deepcopy(self)

        return slot_copy

    #------------------------------------------------------------
    def apply_visibility(self, site="?",delay = 0*u.s,debug=False):
        """

        Apply visibility window and delay to the slot, i.e. change time
        boundaries and re-dress when needed.

        Parameters
        ----------
        site : String, optional
            The site identifier. The default is "?".
        delay : Quantity, optional
            The delay before detection can start. The default is 0*u.s.
        debug : Boolean, optional
            Verbose mode if True. The default is False.

        Returns
        -------
        None.

        """

        self.site = site
        unit = self.grb.tval[0].unit

        #------------------------------------------
        ### Check if t is within a visibility window
        def visible(t):
            ok = False
            for t1, t2 in zip(tstart,tstop):
                if t >= t1 and t <= t2:
                    ok = True
                    continue
            return ok
        #------------------------------------------

        # Get the list of all slices edges
        ticks = []

        for s in self.slices: # Original slices
            ticks.append(s.ts1().value)
            ticks.append(s.ts2().value)

        # Create visibility start and stop time and store
        shift = True
        tstart =[]
        tstop = []
        for tvis in self.grb.vis[self.site].t_true:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                t0 = (tvis[0] - self.grb.t_trig).sec + delay.to(unit).value*shift
                t1 = (tvis[1] - self.grb.t_trig).sec
            if (t1 < t0):
                print(" Delay made this slice vanished, duration =",t1-t0)
            else:
                tstart.append(t0)
                tstop.append(t1)
                ticks.append(t0)
                ticks.append(t1)
                if shift : shift=False

        if (len(tstart) == 0):
            # If delays are such that all slices disappear, then
            # the GRB become not visible
            return False

        ticks=list(set(ticks)) # FIRST remove duplicated entries
        ticks.sort() # THEN Sort (because 'set' change the order)

        # Now loop over all intervals and add to slot
        newslices = []
        idx = 0
        for i in range(len(ticks)-1):
            t1 = ticks[i]
            t2 = ticks[i+1]
            if visible(0.5*(t1+t2)):
                s  = Slice(idx,t1*unit,t2*unit,site=self.site)
                idx+=1
                newslices.append(s)

        # Create a slot from existing and replace content
        self.slices = newslices
        self.name   = "Visible"
        self.tstart = self.slices[0].ts1()
        self.tstop  = self.slices[-1].ts2()
        self.delay  = delay # Loose the delay information after merging
        #self.dress(name = self.name)
        return True
    #------------------------------------------------------------
    def merge(self,slot,name="Merged",debug=False):

        """
        Merge current slot (self) with another slot.

        Parameters
        ----------
        slot : Slot instance
            The slot to be merged with the current slot
        name : TYPE, optional
            DESCRIPTION. The default is "Merged".
        debug : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """

       # Check that the slots are attributed to an authorised site keyword
        if not (self.site in self.grb.site_keys and
                slot.site in self.grb.site_keys) :
            print("merge : site required -> ",self.grb.site_keys)
            print("Slot unchanged")
            return

        # Get the list of all slices edges
        ticks = []
        unit = self.slices[0].ts1().unit

        for obj in [self,slot]:
            for s in obj.slices:
                ticks.append(s.ts1().value)
                ticks.append(s.ts2().value)
        ticks=list(set(ticks)) # FIRST remove duplicated entries
        ticks.sort() # THEN Sort (because 'set' change the order)

        newslices = []
        icount=0
        for i in range(len(ticks)-1):
            t1 = ticks[i]
            t2 = ticks[i+1]
            tmid = 0.5*(t1+t2)
            idxa = self.find(tmid*unit)
            idxb = slot.find(tmid*unit)
            if (idxa >=0):
                if (idxb >=0): loc = "Both"
                else:          loc = self.site
            else:
                if (idxb >=0): loc = slot.site
                else:          loc = None
            # Append only if associated to one or two sites
            if (loc != None):
                s  = Slice(icount,t1*unit,t2*unit,site=loc)
                s.set_id(icount)
                icount+=1
                newslices.append(s)

        # Create a slot from existing and replace content
        self.slices = newslices
        #self.dress(name="Merged")
        self.tstart = min(ticks)*unit
        self.tstop  = max(ticks)*unit
        self.site   = "Both"
        self.delay  = -1*u.s # Loose the delay information after merging
        return

    #------------------------------------------------------------
    def find(self,time):
        """
        Find the slice containing time

        Parameters
        ----------
        time : Quantity (time)
            An observation time

        Returns
        -------
        Index of the slice or -1, -2 if not within the slices

        """
        time = time.to(self.slices[0].ts1().unit)
        for s in self.slices:
            if (time >= s.ts1() and time <= s.ts2()):
                return s.idt()
        return -1

    #------------------------------------------------------------
    def plot(self,ax=None,Eref= 100*u.GeV, **kwargs):
        """
        Shows the slices of the time slot, displaying the measurement
        points (if the slot was dressed).

        Parameters
        ----------
        ax : matplotlib axes
            Axes on which the plots are shown. The default is None.
        Eref : Quantity (energy), optional
            The reference energy used to display the lightcurve
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.
        """

        if (ax==None): ax=plt.subplots()[1]

        ### GRB lightcurve near Eref ---
        iref = np.argmin(np.abs(Eref-self.grb.Eval))
        ax.plot(self.grb.tval,self.grb.fluxval[:,iref],
                marker="o", ls=":", lw=0,color="grey",
                markeredgecolor = "black",
                markeredgewidth = 1.,
                markerfacecolor = "white",
                alpha=0.5,
                label="grb",**kwargs)

        # Compute display interval - takes delay into accont
        # Because a log scale does not accept zero
        tmin = max(1.0,self.slices[0].ts1().value - self.delay.to(u.s).value)
        tmax = self.slices[-1].ts2().value
        ax.set_xlim(0.9*tmin,1.10*tmax)

        ### Show delay window
        if(self.delay.value >= 0):
            label = "Delay {}".format(self.delay)
            ax.axvspan( tmin, tmin + self.delay.to(u.s).value,
                            alpha=0.5,color="grey",label=label)

        ### Slices and flux points
        if self.phys: # If dressed
            color = {"North":"blue", "South":"red", "Both":"purple"}
            for s in self.slices:
                tobs = s.tobs().value
                flux = self.grb.fluxval[s.fid(),iref].value # Not absorbed
                ts1  = s.ts1().value
                ts2  = s.ts2().value
                site = s.site()
                ax.plot(tobs,flux,
                        marker="o",
                        markerfacecolor = color[site],
                        markeredgecolor = color[site],
                        markersize=8,
                        ls="",
                        label=site, **kwargs,)

                ax.hlines(xmin=ts1,xmax=ts2,y=flux,ls="-",color="green",**kwargs)

        ax.set_xlabel("Time since trigger ("+ str(self.grb.tval[0].unit)+")")
        ax.set_ylabel("Flux at " + str(Eref) + " - "
                      + str(self.grb.fluxval[0][0].unit) )
        ax.set_xscale("log")
        ax.set_yscale("log")

        ### Visibility windows
        tref = self.grb.t_trig

        if (self.site == "North" or self.site == "Both"):
            for elt in self.grb.vis["North"].t_true:
                ax.axvspan( (elt[0]- tref).sec,
                            (elt[1]- tref).sec,
                            alpha=0.2,color="tab:blue",label="vis. North")

        if (self.site == "South" or self.site=="Both"):
            for elt in self.grb.vis["South"].t_true:
                ax.axvspan( (elt[0]- tref).sec,
                            (elt[1]- tref).sec,
                            alpha=0.2,color="tab:red",label="vis. South")


        ax.grid("both",linestyle=':')
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        ax.set_title(self.grb.name + " - " + self.site + " - " + self.name)

        return

    #------------------------------------------------------------
    def __str__(self):
        """
        Print our the content of a slice set.
        """

        print(" Observation set : ",self.name," site=",self.site)
        print("   Visibility : ",self.tstart," to ",self.tstop)
        print("    including : ",self.delay, "of delay")
        show_header = True
        for s in self.slices:
            s.print(header=show_header)
            if (show_header):
                show_header = False
        return " "

    #------------------------------------------------------------
    def both_sites(self,delay = {"North":0*u.s, "South":0*u.s},debug = False):
        """
        Find the vibilities in N and S be given the delays for observation.
        In particular in many case the delay will be applied to the first
        in time site (except if the original visibilities differ by a time
        difference less than the delay).

        Parameters
        ----------
        delay : Time Quantity, optional
            Delay for observation. The default is 0*u.s.
        debug : Boolean, optional
            Verbose mode if True. The default is False.

        Returns
        -------
        None.

        """

        # Compute the delays to be applied to the sites
        delta = (self.grb.vis["North"].t_true[0][0]
                 -  self.grb.vis["South"].t_true[0][0] ).sec*u.s

        # Add delay to the first one plus the difference to the second
        if (delta.value < 0):
            if debug: print("North is before South by ",-delta)
            delay_n = delay["North"]
            delay_s = max(0*u.s,delay["South"]-abs(delta))
        elif (delta.value > 0):
            if debug: print("South is before North by ",delta)
            delay_n = max(0*u.s,delay["North"]-delta)
            delay_s = delay["South"]
        elif (delta.value == 0):
            if debug: print("North and South are simultaneaous ",delta)
            delay_n = delay["North"]
            delay_s = delay["South"]

        # Get slots for both sites and merge
        slot_n = self.copy()
        slot_n.apply_visibility(site="North",delay = delay_n)
        slot_s = self.copy()
        slot_s.apply_visibility(site="South",delay = delay_s)
        slot_n.merge(slot_s)

        # Should not we check that one of the two is still visible after
        # the delays are applied ? In fact invisbility is very rare: this
        # would require the two sites to be quasi-simultaneous, with a short
        # initial visibility window

        # Dress with physics
        if (debug): print(slot_n)

        return slot_n

###---------------------------------------------------------------------------
### TESTS
###----------------------------------------------------------------------------
def test_visibility(slotet):
    """
    Test visibility
    """
    tvis =[[18*u.h,40*u.h],
            [3*u.h,60*u.h],
            [3*u.h,44*u.h],
            [3*u.h,39*u.h],
            [16*u.h,60*u.h] ,
            [7*u.h,13*u.h]]
    title = ["Inside",
              "Outside",
              "Before1",
              "Before2",
              "After",
              "Within"]

    for t, label in zip(tvis,title):
        print(" Apply visibility {} - {}".format(t[0],t[1]))
        fig, ax = plt.subplots(nrows=1,ncols=1)
        myset = slotet.copy(name=label)
        myset.mask(t[0],t[1])
        print(myset)
        myset.plot(ax,name=label)
        ax.axvline(t[0].value,ls=":",color="tab:green")
        ax.axvline(t[1].value,ls=":",color="tab:red")
        plt.show()

    return

###----------------------------------------------------------------------------
def test_flux(slotet):
    t = 23*u.h
    flux = slotet.get_flux(t)
    print(" flux at {:8.2f} is {:8.3f}".format(t,flux))

    return

###----------------------------------------------------------------------------
### MAIN
###----------------------------------------------------------------------------
if __name__ == "__main__":

    """
    Create time slices collection for a series of GRB.
    Print them, plot them, test visibility application, merging etc.
    """

    import os
    os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'

    from SoHAPPy import get_grb_fromfile

    import ana_config as cf # Steering parameters
    from utilities import Log

    logfilename = "timeslot.log"  # Log file
    log = Log(name  = logfilename,talk=True)


    cf.dbg_level = 2
    cf.newvis    = True
    cf.altmin    = 10*u.deg # Minimum altitude (original default is 10 degrees)

    delay = 0*u.s
    ifirst= [85] # 1, 85, 204
    # ifirst = 100
    ngrb = 1

    # GRB list to be analysed
    if type(ifirst)!=list:
        grblist = list(range(ifirst,ifirst+ngrb))
    else:
        grblist = ifirst

    # Loop over GRB list
    for i in grblist:

        grb = get_grb_fromfile(i) # Get GRBs

        # Recompute visbility windows
        if (cf.newvis):
            grb.vis["North"].compute(altmin  = cf.altmin,
                                     altmoon = cf.altmoon,
                                     depth   = cf.depth,
                                     skip    = cf.skip,
                                     debug   = False)
            grb.vis["South"].compute(altmin  = cf.altmin,
                                     altmoon = cf.altmoon,
                                     depth   = cf.depth,
                                     skip    = cf.skip,
                                     debug   = False)

        # Printout grb and visibility windows
        if cf.niter<=1 or cf.dbg>0 or cf.ngrb==1 :
            log.prt(grb)
            grb.vis["North"].print(log=log)
            grb.vis["South"].print(log=log)

        # Create the GRB  slot
        origin = Slot(grb,opt="end")
        #print(origin)
        origin.plot(ax=plt.subplots()[1])


        for loc in ["North","South"]:
            if grb.vis[loc].vis_tonight:
                slot = origin.copy(name="loc")
                still_vis = slot.apply_visibility(delay = delay,
                                                  site  = loc)
                #print(slot)
                if (still_vis):
                    slot.dress(irf_dir = cf.irf_dir,arrays=cf.arrays)
                else:
                    print(" Slot vanished after aplying delays")
                print(slot)

        if grb.vis["North"].vis_tonight and grb.vis["South"].vis_tonight:
            slot = origin.both_sites(delay   = delay,
                                       debug =(cf.dbg>1))
            slot.dress(irf_dir = cf.irf_dir,arrays=cf.arrays)

            print(slot)

    log.close()
    print("c'est fini")







