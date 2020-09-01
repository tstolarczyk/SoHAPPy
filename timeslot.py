# -*- coding: utf-8 -*-
"""
This script is a development and test code for manupulating the observation
windows in CTA North and South sites and in particular to to combine
observation of both site and/or on several days

Created on Mon Mar 30 10:36:09 2020

@author: Stolar
"""
import sys

import warnings

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from copy import deepcopy

from obs_slice import Slice
from utilities import warning

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
        self.name  = name
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
    def dress(self,name="dressed"):
        """
        Dress a slice with physiscs information

        Returns
        -------
        None.

        """
        
        self.phys = True
        for s in self.slices: s.dress(self.grb,opt=self.opt)

        return
        
    #------------------------------------------------------------
    def copy(self,name="a_copy"):

        # Could change name an site here
        self.name = name
        slot_copy = deepcopy(self)
 
        return slot_copy
    
    #------------------------------------------------------------
    def apply_visibility(self,t1,t2,site=None):
        """
        Restrict the time slice to the visibility windows [t1,t2].
        The time reference is the GRB trigger time.

        Parameters
        ----------
        t1 : TYPE
            DESCRIPTION.
        t2 : TYPE
            DESCRIPTION.
        site : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        
        if not (site in self.grb.site_keys):
            print("apply_visibility : site required -> ",self.grb.site_keys)
            print("Slot unchanged")
            return
        self.__mask(t1,t2,site=site)
        
        return
    #------------------------------------------------------------
    def __mask(self,tstart,tstop,site=None):
        """
        Recompute the time slices of a slot to account for the 
        visibility.
        The time refernce is the GRB trigger time.
        
        Parameters
        ----------
        tstart : Quantity (Date)
            Start of visibility period
        tstop : Quantity (Date)
            Stop of visibility period
        debug : Boolean, optional

        Returns
        -------
        None.

        """
        
        self.phys = False # Phys information deprecated 

        if tstop < tstart:
            if self.dbg : print("Reverse start and end time of visibility")
            tmp= tstop
            tstop = tstart
            tstart = tmp
            
        # Retrieve start and stop times of the slices
        t_unit = self.grb.tval[0].unit
        
        ts1  = np.asarray([s.ts1().value for s in self.slices])
        ts2  = np.asarray([s.ts2().value for s in self.slices])
        
        # Find in which slices are the start/stop of the visibility
        (id1, id2) = (self.find(tstart), self.find(tstop))
        
        if self.dbg : print("Visibility start / stop in slices ",id1,"/",id2)
        
        # Convert visibility values to slice units - keep original unit
        tvis_unit = tstart.unit
        tstart    = tstart.to(t_unit).value
        tstop     = tstop.to(t_unit).value
        
        # Select core, complete slices strictly within visibility window
        # Aplly name
        newslot = self.slices[(ts1>tstart) * (ts2<tstop)]

        if self.dbg :
            if (len(newslot)):
                print(" Slices retrieved : ")
                for s in newslot : print("*    core     *", s)
            else: print(" Core is empty")
        
        if (id1 == id2): # Visibility window is within a single slice
            newslice = [Slice(id1, 
                              tstart*t_unit, 
                              tstop*t_unit, 
                              site=site)
                        ]
            newslot  = np.hstack( (newslice, newslot) )

        else: # Visibility window spans over more than one slice
            
            if (id1 >= 0): # Add partial slice piece BEFORE the core
                tmin = max(self.slices[id1].ts1(),tstart*t_unit) # if tstart < start
                newslice = [Slice(id1,
                                  tmin,  
                                  self.slices[id1].ts2(),
                                  site=site)
                            ]
                newslot  = np.hstack( (newslice, newslot) )
                if self.dbg : print("* first added *", newslot[0])
            
            if (id2 >=0): # Add partial slice piece AFTER the core
                newslice = [Slice(id2,
                                  self.slices[id2].ts1(), 
                                  tstop*t_unit,
                                  site=site)
                            ]
                newslot = np.hstack( (newslot, newslice) )
                if self.dbg : print("* last added  *", newslot[-1])
                  
        self.site     = site
        self.slices   = newslot
        self.tstart    = tstart*tvis_unit
        self.tstop     = tstop*tvis_unit
        
        # Renumber the slices - redefine the site
        for i, s in enumerate(self.slices): 
            s.set_id(i)
            s.set_site(site)

        return
    #------------------------------------------------------------
    def merge(self,slot,name="Merged"):
        """
        Merge the current slot (self) with "slot"
        Slots are ordered from their visibility start.
        
        If the GRB duration is shorter than both visbility wndows, then they 
        will have exactly the same time slices (Slot) in both site. 
        The Slot is unchanged and the slices get two IRF attributed.
        
        In the other cases, i.e. the GRB duration goes accross the visibility 
        windows in some way, the two slots are not identical (at least it is 
        hard to imagine and no protection against this is implemented).
        
        The time refernce is the GRB trigger time.
        
        Several cases can occur :
        1. The visibility of slot2 starts after the end of slot1
        slots do not overlap they are simply concatenated.
        
        1.  : +------------------+ S1
            :                        +-------------------+ S2
            : +------------------+   +-------------------+ merged
            :        S1                        S2   
        
        2a  : +------------------+ S1
            :       +---------------------+ S2
            : +-----+------------+--------+ merged
            :   S1       S1,S2       S2
        
        2b  : +------------------+ S1
            :        +-----+ S2
            : +------+-----+-----+
            :    S1    S1,S2  S1
        
        2. The visibility of slot2 starts within the one of slot1
            2a. slot2 ends after slot1 ends
            There are 2 perios to merge :
                - a slot1 sloct before the overlap
                - a common slot1, slot2 slot that overlaps.
            2b. slot2 ends before slot1 (cannot happen for N,S merging)
                - compared to previous case, there an extra slot from 
                slot2 to be added
        
        Parameters
        ----------
        slot : a Sliceslot class instantiation
            A slot to be merged with the current simulation

        Returns
        -------
        A modified slot instantiation containing the merged simulations

        """
        
        # Check that the slots are attributed to an authorised site keyword
        if not (self.site in self.grb.site_keys and 
                slot.site in self.grb.site_keys) :
            print("merge : site required -> ",self.grb.site_keys)
            print("Slot unchanged")
            return  
        if (self.site == slot.site): 
            newsite = slot.site
        else:
            newsite = "Both"
            
        # Check whether the two slots have all slices identical
        # This happens in particular for prompt simulations that fall 
        # within a single slice
        if self.tstart == slot.tstart and len(self.slices)==len(slot.slices):
                if self.dbg:
                    warning("The two slots start at same moment and "
                        + "have same number of slices : considered IDENTICAL")         
                self.name   = "Doubled"
                self.site   = newsite
                self.delay  = -1*u.s # Loose the delay information after merging 
                for i, s in enumerate(self.slices): s.set_site(self.site)
                return      
        
        # Order slots with time, slot1 first, selot2 second
        if (self.tstart <= slot.tstart): 
            slot1 = self
            slot2 = slot
            if self.dbg: 
                print(slot1.name," start is before/coincident with",slot2.name)

        else:
            slot2 = self
            slot1 = slot
            if self.dbg: 
                print(slot1.name," starts after ",slot2.name)

        
        if (slot2.tstart > slot1.tstop):
        # 1. The visibility of slot2 starts after the end of slot1
        # Slots do not overlap they are simply concatenated. 
        # Slice sites are unchanged
            self.slices = np.hstack( (slot1.slices,slot2.slices) )      
        # End of case 1
        else:
        # 2. The visibility of slot2 starts within the one of slot1
        #First create the slotbefore overlap if any
            if self.dbg: print( "2nd slot start within first -> Overlapping")       
            
            if (slot1.tstart != slot2.tstart): 
                # There is a slot before the overlap
                slotbef = slot1.copy(name="Beg S1")
                slotbef.__mask(slot1.tstart, 
                               slot2.tstart, 
                               site = slot1.site)
                if self.dbg: print(" * 1st part slot1 -->",slotbef)
                before = True
            else:
                # The slots start at the same moment, nothing before overlap
                before = False
                if self.dbg: print(" The 2 slots start at same time")
                
            if (slot2.tstop > slot1.tstop):
            # 2a. slot2 ends after slot1 ends
                # A common slot1, slot2 slot that overlaps.
                slotoverlap = slot1.copy(name="Overlap")
                slotoverlap.__mask(slot2.tstart,
                                   slot1.tstop,
                                   site=newsite)
                
                if self.dbg: print(" * Overlap -->",slotoverlap)            
                
                # A remaining piece from slot2 slot.
                slotafter = slot2.copy(name="End S2")
                slotafter.__mask(slot1.tstop,
                                 slot2.tstop,
                                 site=slot2.site)
                after = True
                if self.dbg: print(" * Last part slot2 -->",slotafter)            

            elif (slot1.tstop == slot2.tstop):
                if self.dbg: print(" The 2 slots stop at same time")
                after = False
                
            else: # (slot2.tstop < slot1.tstart)
            # 2b. slot2 ends before slot1 (cannot happen for N,S merging)
                # A common slot1, slot2 slot that overlaps.
                slotoverlap = slot1.copy(name="overlap")
                slotoverlap.__mask(slot2.tstart,
                                 slot2.tstop,
                                 site=newsite)
                if self.dbg: print(" * First part overlap -->",slotoverlap)            
                
                # A remaing piece from slot1 slot.
                slotafter = slot1.copy(name="end S1")
                slotafter.__mask(slot2.tstop,
                               slot1.tstop,
                               site=slot1.site)
                after = True
                if self.dbg: print(" * Last part S1 -->",slotafter)                        
                
            # Now stacking - overlap always exist
            self.slices = np.hstack((slotbef.slices,slotoverlap.slices)) \
                          if before else slotoverlap.slices
            
            self.slices =  np.hstack((self.slices,slotafter.slices)) \
                          if after else self.slices
        
        self.name   = "Merged"
        self.tstart = min(slot1.tstart,slot2.tstart)
        self.tstop  = max(slot1.tstop,slot2.tstop) 
        self.site   = newsite
        self.delay  = -1*u.s # Loose the delay information after merging 

        # Renumber the slices
        for i, s in enumerate(self.slices): s.set_id(i)
            
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
        time = time.to(self.slices[0].ts1().unit).value
        
        tstart = np.asarray([s.ts1().value for s in self.slices])
        tstop  = np.asarray([s.ts2().value for s in self.slices])        
        
        if   (time < tstart[0]) : idx = -1
        elif (time > tstop[-1]): idx = -2
        else : idx = np.max(np.where(time>=tstart))
            
        return idx

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
        
        # GRB lightcurve near Eref
        iref = np.argmin(np.abs(Eref-self.grb.Eval))
        ax.plot(self.grb.tval,self.grb.fluxval[:,iref],
                marker="+",
                ls=":",
                label="grb",**kwargs)

        
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
        
        tref = self.grb.t_trig
        if (site == "Both"): sitelist = ["North","South"]
        else: sitelist = [site]
        
        for loc in sitelist:
            for elt in self.grb.t_true[loc]:
                ax.axvline( (elt[0]- tref).sec,
                           ls=":",color="tab:green",label="Start " + loc)
                ax.axvline( (elt[1]- tref).sec,
                           ls=":",color="tab:red",label="Stop " + loc)  
            
            if(self.delay.value >= 0): 
                ax.axvline( (self.grb.t_true[loc][0][0] - tref + self.delay).sec,
                           ls=":",color="tab:blue",
                           label="Delay " + loc)       
        ax.grid("both")
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
            if (show_header): 
                s.print_header()
                show_header = False
            print(s)
            
        return " "
    #------------------------------------------------------------
    def single_site(self, site="None", delay=0*u.s, debug=False):
        """
        Compute the time slices from the GRB originals, taking into account
        the visbility windows and the detection delays.

        Parameters
        ----------
        site : TYPE, optional
            DESCRIPTION. The default is "None".
        day_after : TYPE, optional
            DESCRIPTION. The default is 0.
        delay : TYPE, optional
            DESCRIPTION. The default is 0*u.s.
        debug : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
                                       
        # Should check that the visibility window list is in the right order
        # athough it had been build in that order
        first = True
        i = 0
        for elt in self.grb.t_true[site]:
            
            tvis1 = (elt[0]-self.grb.t_trig).sec*u.s 
            tvis2 = (elt[1]-self.grb.t_trig).sec*u.s   
            
            if (first):
                slot = self.copy(name=self.grb.name+"-"+site+"-naked") 
                tvis1 = tvis1 + delay
                slot.apply_visibility(tvis1,tvis2,site=site)  
                first = False
            else:
                i+=1
                slot1 = self.copy(name=self.grb.name+"-"+site+str(i)) 
                slot1.apply_visibility(tvis1,tvis2,site=site)  
                slot.merge(slot1,name=site)
        #print(slot)
            
        slot.delay = delay
        slot.dress(name=self.grb.name+"-"+site) # Dress the slot with physics
                                
        return slot
    
    #------------------------------------------------------------    
    def both_sites(self, delay = 0*u.s, debug = False):
        """
        

        Parameters
        ----------
        delay : TYPE, optional
            DESCRIPTION. The default is 0*u.s.
        debug : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """            
        
        # Compute the delays to be applied to the sites
        delta = (self.grb.t_true["North"][0][0] 
                 -  self.grb.t_true["South"][0][0] ).sec*u.s
            
        # Add delay to the first one plus the difference to the second 
        if (delta.value < 0):
            if debug: print("North is before South by ",-delta)
            delay_n = delay
            delay_s = max(0*u.s,delay-abs(delta))
        elif (delta.value > 0):
            if debug: print("South is before North by ",delta)
            delay_n = max(0*u.s,delay-delta)
            delay_s = delay
        elif (delta.value == 0):
            if debug: print("North and South are simultaneaous ",delta)
            delay_n = delay
            delay_s = delay
       
        # Get merged slots for both sites
        slot_n = self.single_site(site="North",delay = delay_n)
        slot_s = self.single_site(site="South",delay = delay_s)

        # Merge the two site slots
        slot = slot_n.copy()
        slot.merge(slot_s)
        # print(slot_b)
            
        # Dress with physics
        slot.dress(name="Both")
        if (debug): print(slot)   
                    
        return slot
        
 
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
    
    import warnings
    import gammapy
    from SoHAPPy import get_grb_fromfile, init, get_delay
    import grb_plot as gplt
    
    if (gammapy.__version__ == "0.12"):
        from   gammapy.spectrum.models import Absorption
    if (gammapy.__version__ == "0.16"):
        from gammapy.modeling.models import Absorption
    
    import ana_config as cf # Steering parameters
    from utilities import Log
    
    
    logfilename = "timeslot.log"  # Log file
    log = Log(name  = logfilename,talk=True)
    
    
    cf.dbg_level  = 0
    
    ifirst= [64]
    ngrb = 1
    # GRB list to be analysed
    if type(ifirst)!=list:
        grblist = list(range(ifirst,ifirst+ngrb))
    else:
        grblist = ifirst    
    
    # Loop over GRB list accordng to the config. file
    import visibility as v
    for i in grblist:

        # Get GRBs
        grb = get_grb_fromfile(i)
        print(grb)
        # Create a slot 
        #print(slot_original)
        slot_original = Slot(grb,opt="end")

        # Compute visibility - update GRB information
        for loc in ["North","South"]:
            vis = v.Visibility(grb,loc=loc) # Constructor    
            vis.compute(altmin=cf.altmin,depth=10)
            grb.update_visibility(vis)         
            # grb.show_visibility(loc=loc,log=log)
            # gplt.visibility_plot(grb,loc=loc)
            
            # if grb.vis_tonight[loc]:
            #     slot = slot_original.single_site(site=loc)
            #     slot.plot()
            #     print(slot)
        
        if grb.vis_tonight["North"] and grb.vis_tonight["South"]:
            slot_b = slot_original.both_sites(delay = get_delay(), debug =False)
            slot_b.plot()
            print(slot_b)
                
        ### TO BE REWRITTEN
    
        # slot_original = Slot(grb,opt="end")
        # # print(slot_original)
        # # slot_original.plot(ax=plt.subplots()[1])   
        
        # if (grb.vis_tonight["North"]):
            
        #     slot_n = slot_original.copy(name="North naked") 
        #     #print(slot_n)
        #     with warnings.catch_warnings():
        #         warnings.filterwarnings("ignore")
        #         tn_1 = (grb.t_true["North"][0] - grb.t_trig).sec * u.s
        #         tn_2 = (grb.t_true["North"][1] - grb.t_trig).sec * u.s
        #     slot_n.apply_visibility(tn_1,tn_2,site="North",debug=False)
        #     #print(slot_n)
        #     slot_n.dress(debug=False)
        #     slot_n.plot(ax=plt.subplots()[1])
        #     print(slot_n)
            
        # if (grb.vis_tonight["South"]):
        #     slot_s = slot_original.copy(name="South naked")
        #     with warnings.catch_warnings():
        #         warnings.filterwarnings("ignore")
        #         ts_1 = (grb.t_true["South"][0] - grb.t_trig).sec * u.s
        #         ts_2 = (grb.t_true["South"][1] - grb.t_trig).sec * u.s
        #     slot_s.apply_visibility(ts_1,ts_2,site="South",debug=False)
        #     # print(slot_s)
        #     slot_s.dress(debug=False)
        #     print(slot_s)
        #     slot_s.plot(ax=plt.subplots()[1])    
            
        # if (grb.vis_tonight["South"]) and (grb.vis_tonight["North"]):
        #     slot_b = slot_n.copy()
        #     # print(slot_b)
        #     slot_b.merge(slot_s,debug=False)
        #     # print(slot_b)
        #     slot_b.dress(name="Both")
        #     print(slot_b)
        #     slot_b.plot(ax=plt.subplots()[1])
                
    print("c'est fini")
    # tn_1 = 9*u.h
    # tn_2 = 32*u.h
    # slot_n.mask(tn_1,tn_2)
    # slot_n.dress()
    #slot_original.dress()
    #print(slot_original)
    
    # test find
    # slot_original.find(21*u.h)
    
    #test_flux(slot_original)
    
    #test_visibility(slot_original,debug=False)
    # myslot = slot_original.copy()
    # myslot.mask(12*u.h,42*u.h)
    # print(myslot)
    # myslot.plot(ax=plt.subplots()[1])
    # myslot.dress()
    # print(myobs)
    # myslot.plot(ax=plt.subplots()[1])    
     
    # Test merging
    
    # # North
    # slot_n = slot_original.copy(name="North",site="North")
    # # print(slot_n)
    # tn_1 = 9*u.h
    # tn_2 = 32*u.h
    # slot_n.mask(tn_1,tn_2)
    # # slot_n.dress()
    # print(slot_n)
    
    # # South
    # slot_s = slot_original.copy(name="South",site="South")
    # # print(slot_s)
    # ts_1 = 15*u.h
    # ts_2 = 42*u.h
    # slot_s.mask(ts_1,ts_2,)
    # # slot_s.dress()
    # print(slot_s)
    
    # # Plot South and NOrth
    # fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(10,5))
    # slot_n.plot(ax=ax[0])
    # slot_s.plot(ax=ax[1])
    # plt.tight_layout()
    # plt.show()
    
    # Merged
    # slot_b = slot_n.copy(name="North")
    # slot_b.merge(slot_s,debug=True)
    # print(slot_b)
    # slot_b.dress()
    # print(slot_b)
    # slot_b.plot(ax=plt.subplots()[1])
    
    
      
    
    
    