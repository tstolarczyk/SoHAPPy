# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:04:35 2020

@author: Stolar
"""
import numpy as np

import astropy.units as u
from astropy.utils import iers

from irf import IRF
from irf_onaxis import CTAPerf_onaxis

###----------------------------------------------------------------------------
class Slice():
    """
    A (time) slice has an observation point fixed in time, associated to a 
    physical flux.
    The observation can be simply set at the beginning or at the end of the 
    slice or at a date defined in a more elabrated way.
    Variations due to these various hypothesis highly depends on
    the slice length and how the flux variates in the slice.
    Note that a slice can have several performance information (IRF) since the
    same time slice can belong to 2 different sites (overlapping 
    observation windows).
    """
    #--------------------------------------------------------------------------
    def __init__(self,idt, t1, t2,site="?"):
        """
        
        This defines a so called naked slice. It has no physical information.
        
        Parameters
        ----------
        idt : TYPE
            DESCRIPTION.
        t1 : TYPE
            DESCRIPTION.
        t2 : TYPE
            DESCRIPTION.
        tobs : TYPE, optional
            DESCRIPTION. The default is 0.
        flux : TYPE, optional
            DESCRIPTION. The default is [0].
        site : TYPE, optional
            DESCRIPTION. The default is "?".
        perf : TYPE, optional
            DESCRIPTION. The default is [].

        Returns
        -------
        None.
        
        """
        
        self.__idt    = idt    # Slice identifier
        self.__ts1    = t1     # Slice start
        self.__ts2    = t2     # Slice stop
        self.__tobs   = 0*u.s  # Observation point
        self.__f_id   = -1     # Spectrum slice id in the grb object 
        self.__perf   = []     # Detector performance(s) irf at obs. point
        self.__site   = site   # Sites North, South, Both
            
        return
        
    def idt(self): return self.__idt
    def ts1(self): return self.__ts1
    def ts2(self): return self.__ts2
    def tobs(self): return self.__tobs
    def fid(self): return self.__f_id
    def site(self): return self.__site
    def perf(self): return self.__perf
    
    def set_id(self,idt): self.__idt = idt 
    def set_site(self,site="?"):
        """
        Possible site are the following : North, South, Both 
        More could be added.
        Todo; if more possibilities are adeed might be worth to create a list

        Parameters
        ----------
        site : TYPE, optional
            DESCRIPTION. The default is "?".

        Returns
        -------
        None.

        """
        self.__site = site
        
        return
    #--------------------------------------------------------------------------
    def dress(self,grb,opt="end",debug=False):
        """
        Add physical information to the slice

        Parameters
        ----------
        grb : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        self.obs_point(opt)
        self.get_flux(grb,self.__tobs)    
        self.get_perf(grb,debug)
        return

    #------------------------------------------------------------
    def obs_point(self,opt):
        """
        

        Parameters
        ----------
        Get the observation point time depending on where the flux will be
        taken

        Parameters
        ----------
        tmin : Quantity (date)
            Slice start date
        tmax : Quantity (date)
            Slice stop date
        opt : TYPE
            DESCRIPTION.

        Returns
        -------
        tobs : TYPE
            DESCRIPTION.

        """    
        if   (opt =="end")    : self.__tobs = self.__ts2
        elif (opt =="start")  : self.__tobs = self.__ts1
        elif (opt =="mean")   : self.__tobs = 0.5*(self.__ts1+self.__ts2)
        elif (opt =="interp") : self.__tobs = self.__ts2
        else: self.__tobs=0*(self.__ts1.unit)
        
        return
    
    #------------------------------------------------------------
    def get_flux(self,grb,t):
        """
        Compute the flux to be associated to a given  time
        So far get the flux at the closest measurement point.
        """
        t = t.to(grb.tval[0].unit).value
        self.__f_id = np.argmin(np.abs(t - grb.tval.value))
        
        return

    #--------------------------------------------------------------------------
    def get_perf(self,grb,debug=False):
        """
        
        Obtained the best performance for a given slice, indpendently of
        the observation point that was chosen.
        The altitue and azimuth are computed for at the slice edges.
        The performance is obatined at the beginning of the slice where
        the flux is higher.
        The observation time is the length of the slice.
        
        """

        # Observation window 
        obstime = self.__ts2 - self.__ts1 # Duration (not cumulated !)
       

        (self.__altaz , self.__perf, self.__Emin, self.__Emax) = ([],[],[],[])
        

        site_list = []
        if (self.__site == "North") or (self.__site == "South"): 
            site_list=[self.__site]
        elif (self.__site == "Both"):
            site_list = grb.site_keys
            
        for site in site_list:
            
            altaz =  grb.altaz(dt=self.__ts1,loc=site)
                
            irf_file = IRF.get_irf_file(theta   = 90*u.degree-altaz.alt,
                                        phi     = altaz.az,
                                        obstime = obstime,
                                        khem    = site)
            if (debug>1): print(self.__idt,"/",site,"-->",irf_file)
            
            perfi        = CTAPerf_onaxis.read(irf_file)

            self.__perf.append(perfi)

      
        if (debug):
            print("                       ===> {} IRF found"
              .format(len(self.__perf)),end="\n")
    
        return
    
    #--------------------------------------------------------------------------
    def print_header(self):
        print("  {:>2s} {:>10s} {:>10s} {:>10s} {:>3s} {:5s} {:4s}"
              .format("#",
                      "Start ("+str(self.__ts1.unit)+")",
                      "Stop ("+str(self.__ts2.unit)+")",
                      "Obs.("+str(self.__tobs.unit)+")",
                      "Fid",
                      "Site",
                      "Perf"))
        print("  {:>2s} {:>10s} {:>10s} {:>10s} {:>3s} {:>5s} {:>4s}"
              .format("--",
                      "----------",
                      "----------",
                      "----------",
                      "---",
                      "-----",
                      "----"))    
        return            
        
    #--------------------------------------------------------------------------
    def __str__ (self):
        txt = "  {:2d} {:10.2f} {:10.2f} {:10.2f} {:3d} {:>5s} {:4d}" \
        .format(self.__idt,
                self.__ts1.value,
                self.__ts2.value,
                self.__tobs.value,
                self.__f_id,
                self.__site,
                len(self.__perf))
        return txt