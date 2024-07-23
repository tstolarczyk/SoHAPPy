# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:04:35 2020

@author: Stolar
"""
from pathlib import Path

import numpy as np
import astropy.units as u

from irf import IRF

__all__ = ["Slice"]

###----------------------------------------------------------------------------
class Slice():
    """
    A (time) slice is a time window associated to a physical flux, an
    observation time within the window, a list of sites and the corresponding
    IRF as for combined analysis a slice can be associated to more than one
    site.
    The observation time can be simply set at the beginning or at the end of
    the slice or at a date defined in a more elaborated way. Variations due to
    these various hypothesis highly depends on the slice length and how the
    flux varies in the slice.
    """
    #--------------------------------------------------------------------------
    def __init__(self,idt, t1, t2, site="?"):
        """
        This defines a so called naked slice. It has no physical information
        defined.

        Parameters
        ----------
        idt : integer
            Slice identtifier
        t1 : Quantity Time
            Slice start time
        t2 : Quantity Time
            Slice stop time.
        tobs : Quantitiy Time, optional
            Observation time for that slice. The default is 0.
        flux : Float, optional
            Flux value at observation point. The default is [0].
        site : String, optional
            The site(s) associated to that slice. The default is "?".
        perf : CTAIrf, optional
            A list of CTA Irf instances. The default is [].

        Returns
        -------
        None.

        """

        self.__idt    = idt       # Slice identifier
        self.__ts1    = t1        # Slice start
        self.__ts2    = t2        # Slice stop
        self.__tobs   = 0*u.s     # Observation point
        self.__f_id   = -1        # Spectrum slice id in the grb object
        self.__irf    = []        # Detector IRFs at obs. point
        self.__site   = site      # Sites North, South, Both

    #--------------------------------------------------------------------------
    def idt(self) :
        """
        Get the Slice identifier.

        Returns
        -------
        integer
            Slice identifier.

        """
        return self.__idt
    #--------------------------------------------------------------------------
    def ts1(self) :
        """
        Get the slice start time.

        Returns
        -------
        Astropy Quantity
            Start time.

        """
        return self.__ts1
    #--------------------------------------------------------------------------
    def ts2(self) :
        """
        Get the slice stop time.

        Returns
        -------
        Astropy Quantity
            Stop time.

        """

        return self.__ts2
    #--------------------------------------------------------------------------
    def tobs(self):
        """
        Get the slice observation point.

        Returns
        -------
        Astropy Quantity
            Observation time.

        """
        return self.__tobs

    #--------------------------------------------------------------------------
    def fid(self) :
        """
        Get the flux identifier.

        Returns
        -------
        integer
            Flux identifuier in the GRB object.

        """
        return self.__f_id

    #--------------------------------------------------------------------------
    def site(self):
        """
        Get Site location.

        Returns
        -------
        String
            A String describing the site(s) attached to this Slice instance.
        """

        return self.__site
    #--------------------------------------------------------------------------
    def irf(self):
        """
        get detector IRF(s) at obs. point

        Returns
        -------
        IRF
            List of IRF attached to this Slice instance.

        """
        return self.__irf

    #--------------------------------------------------------------------------
    def set_id(self,idt):
        """
        Set the slice identifier

        Parameters
        ----------
        idt : integer
            Slice identifier.

        Returns
        -------
        None.

        """
        self.__idt = idt

    #--------------------------------------------------------------------------
    def set_site(self,site="?"):
        """
        Possible sites are the following : North, South, Both
        More could be added.
        Todo. if more possibilities are added might be worth to create a list

        Parameters
        ----------
        site : String, optional
            Site(s) associated to the site. The default is "?".

        Returns
        -------
        None.

        """
        self.__site = site

    #--------------------------------------------------------------------------
    def dress(self,grb,
                   irf_dir = Path("./"),
                   arrays  = None,
                   opt     = "end",
                   zenith  = None):
        """
        Add physical information to the slice.

        Parameters
        ----------
        grb : GammaRayBurst
            A GRB instance
        irf_dir : Path
            IRF top directory. The default is the local directory.
        arrays : string, optional
            String defining the telescope (sub)array. The defualt is None.
        opt : string, optional
            Define the position at which the flux is estimated. The default is "end".
        zenith : astropy Quantity, optional
            Value of the forced zenith angle. The default is None.

        Returns
        -------
        None.

        """
        self.obs_point(opt)
        self.get_flux(grb,self.__tobs)
        self.get_perf(grb, irf_dir= irf_dir, arrays=arrays, zenith=zenith)

    #------------------------------------------------------------
    def merge(self, next_slice):
        # Set end time of current slice to end time of another slce
        # Leaves the rest unchanged
        self.__ts2 = next_slice.ts2()
        self.__tobs = next_slice.tobs()

    #------------------------------------------------------------
    def obs_point(self, opt):
        """
        Get the observation point time depending on where the flux will be
        considered.

        Parameters
        ----------
        opt : string
            String defining where the flux is estimated in the slice.
            "end", "start", "mean", etc.

        Returns
        -------
        None.

        """

        if   opt == "end":
            self.__tobs = self.__ts2
        elif opt == "start":
            self.__tobs = self.__ts1
        elif opt == "mean":
            self.__tobs = 0.5*(self.__ts1+self.__ts2)
        elif opt == "interp":
            self.__tobs = self.__ts2
        else:
            self.__tobs=0*(self.__ts1.unit)

    #------------------------------------------------------------
    def get_flux(self,grb,t):
        """
        Compute the flux to be associated to a given  time.
        So far get the flux at the closest measurement point.

        Parameters
        ----------
        grb : GammaRayBurst
            Current GRB instance
        t : Quantity Time
            The time at which the flux should be estimated

        Returns
        -------
        None.

        """

        t = t.to(grb.tval[0].unit).value
        self.__f_id = np.argmin(np.abs(t - grb.tval.value))

    #--------------------------------------------------------------------------
    def get_perf(self, grb,
                       irf_dir = Path("./"),
                       arrays  = None,
                       zenith  = None,
                       debug   = False):
        """
        Obtain the best performance for a given slice, indpendently of
        the observation point that was chosen.
        The altitude and azimuth are computed for the slice start since
        the performance is obatined at the beginning of the slice where
        the flux is higher. The observation time is the length of the slice.

        Parameters
        ----------
        grb : GammaRayBurst
            The current GRB instance.
        debug : boolean, optional
            If true, display chosen IRF. The default is False.

        Returns
        -------
        None.

        """
        # debug=2

        # Observation window duration
        obstime = self.__ts2 - self.__ts1 # Duration (not cumulated !)

        site_list = []
        if self.__site in ("North", "South"):
            site_list=[self.__site]
        elif self.__site == "Both":
            site_list = ["North","South"]

        for site in site_list:
            # Altitude at slice start where the flux is expected larger
            altaz =  grb.altaz(dt=self.__ts1,loc=site)

            # Zenith is obtained from altaz except if fixed beforehand
            if zenith is None or zenith is False:
                zenith   = 90*u.degree-altaz.alt
            else:
                zenith = u.Quantity(zenith) # Was string so far

            irf = IRF.from_observation(zenith   = zenith,
                                       azimuth  = None,
                                       obstime  = obstime,
                                       loc      = site,
                                       irf_dir  = irf_dir,
                                       subarray = arrays[site])
            if debug>1:
                print(self.__idt,"/",site,"-->",irf)

            self.__irf.append(irf)

        if debug:
            print(f"                       ===> {len(self.__irf)} IRF found")

    #--------------------------------------------------------------------------
    def print(self,header=False):
        """
        Allow displaying

        Returns
        -------
        txt : String
            Slice information

        """
        if header:
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
        print("  {:2d} {:10.2f} {:10.2f} {:10.2f} {:3d} {:>5s} {:4d}"
              .format(self.__idt,
                self.__ts1.value,
                self.__ts2.value,
                self.__tobs.value,
                self.__f_id,
                self.__site,
                len(self.__irf)),end="")
        for perf in self.__irf:
            print(" ",perf.filename.parts[-2],end="")
        print()
