# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:17:27 2020

This modules compute is organise around the :class:`visibility`class.
It combines the rise and set, and night (twilight) windows to the GRB data window (defined by the two extreme GRB time of the data points ) to produce the visibility windows for a given GRB on a given site.
It can also be used to check that the ones given by default are in agreement for the default horizon value (10 degrees).
The windows found are later passed to the GammaRayBurst class for further processing.

@author: Stolar
"""
import math
import astropy.units as u
from astropy.time import Time
from astroplan import Observer, FixedTarget

# Refresh...
# from astroplan import download_IERS_A
# download_IERS_A

__all__ = ["Visibility"]
###------------------------------------------------------------------------
def Df(x):
    """
    Returns a time in Julian day format if the argument is finite
    """
    if math.isfinite(x):
        return Time(x,format="jd")
    else:
        return x
###------------------------------------------------------------------------
def Dp(x):
    """
    Returns a time in ISO format if the argument is an astropy Time, and if
    not returns the argument (unchanged).
    """
    if isinstance(x,Time):
        return x.iso
    else:
        return x

###############################################################################
class Visibility():
    """
    This class defines the visibility periods for a given GRB and a certain
    site (North or South).
    It can be called to update the default visibilities given with the GRB data
    file, in particular to modify the minimum altitude to declare a source
    being visible.
    One of the method, :method:`check` is used to compare the visibility windows with the one
    given by default in the GRB file.
    """

    ###------------------------------------------------------------------------
    def __init__(self,grb,loc=None, npt=150):
        """
        Visibility constructor
        I did not find how to have this displayed with automodapi
        Parameters
        ----------
        grb : GammaRayBurst
            A GammarayBurt instance
        loc : string, optional
            Site poisiton, either North or South. The default is None.

        Returns
        -------
        ok : TYPE
            DESCRIPTION.

        """
        self.grb     = grb
        self.loc     = loc
        self.depth   = 0
        self.vis     = False
        self.prompt  = False
        self.tonight = False
        self.t_grb   = []
        self.t_night = []
        self.t_above = []
        self.t_vis   = []

        return

    ###------------------------------------------------------------------------
    def compute(self,altmin = 10*u.degree,
                     depth  = 3*u.day,
                     end    = 0,
                     skip   = 0,
                     npt    = 150, debug=False):

        """
        Compute the visibility periods for a given GRB and site.
        The constructor takes as arguments a grb (GammaRayBurst) and a
        location (string).
        The algorithm is the following:

            1.	Build periods for:
                - Nights (night), at least two consecutive nights (more is possible).
                - Above horizon (above), at least two consecutive rise/set.
                - Bright (bright) from trigger time and to last GRB data point in time.

            2.	Put all start and stop of all these periods in a list, sort the
            list, resulting in a series of ordered “ticks”.

            3.	For each tick pairs [t1,t2]:
                - Compute the mean time $$0.5*(t1+t2)$$
                - Check if that time belongs to one of the (night), (above) or (bright)
                intervals.
                - If so get True, otherwise False.

            4.	For each of the tick pairs compute the Boolean:
                (visible) = (bright)*(above)*(night),
                and get the visibility windows when True.

            5.	Restrict the visible=True window to :
            stop time < depth*u.dayafter trigger.

        Parameters
        ----------
        altmin : float, optional
            Minimum altitude. The default is 10*u.deg.
        depth : Quantity Time, optional
            Depth up to which time windows are search for (at least one day)
            The default is 3 days, and the reference time is the start of the
            visibility window (i.e. the window end can be beyond the depth)
        end : integer, optional
            Defines the crietria to accept a visibility window within the
            depth. If 0, the window has to start before depth, if 1 it has to
            stop before depth.
        npt : integer, optional
            Number of grid points for horizon crossing. The default is 150.
        debug : bool, optional
            Print additional comments at excecution time if True .
            The default is False.

        """

        ###---------------------------------------------------
        def valid(t0,tslices):
            ok = False
            for slices in tslices:
                if (t0 >= slices[0] and t0 <= slices[1]):
                    if (ok == False): ok = True
            return ok
        ###---------------------------------------------------

        self.depth  = depth
        self.altmin = altmin

        obs  = Observer(location  = self.grb.pos_site[self.loc],
                        name = self.loc, timezone ="utc")
        target = FixedTarget(coord=self.grb.radec, name=self.grb.name)

        # Temporrary local variables in Julian days (to use sort)
        # Will be copied as Time objects later into the class variable
        t_grb   = []
        t_night = []
        t_above = []
        t_vis   = []

        ### GRB ---
        t_grb.append([self.grb.t_trig.jd, (self.grb.t_trig
                                           + self.grb.tval[-1]).jd])
        ### NIGHT ---
        night  = obs.is_night(self.grb.t_trig, horizon = -18*u.deg)

        # Get the first night : can be the current night
        # Omitit if requested
        if (night): search="previous"
        else: search = "next"

        inight = 0
        t_dusk = obs.twilight_evening_astronomical(self.grb.t_trig,
                                                   which = search,
                                                   n_grid_points = npt)
        t_dawn = obs.twilight_morning_astronomical(t_dusk,
                                                   which="next",
                                                   n_grid_points = npt)
        if (skip ==0): t_night.append([t_dusk.jd, t_dawn.jd])

        # Add subsequent nights until reaching the end of GRB data
        while (t_dusk.jd < t_grb[0][1]):
            inight +=1
            t_dusk = obs.twilight_evening_astronomical(t_dawn,
                                                        which = "next",
                                                        n_grid_points = npt)
            t_dawn = obs.twilight_morning_astronomical(t_dusk,
                                                        which = "next",
                                                        n_grid_points = npt)
            if (inight >= skip): t_night.append([t_dusk.jd, t_dawn.jd])

        # For test, add the previous night
        # t_dawn0 = obs.twilight_morning_astronomical(self.grb.t_trig,
        #                                            which="previous",
        #                                            n_grid_points = npt)
        # t_dusk0 = obs.twilight_evening_astronomical(t_dawn0,
        #                                             which="previous",
        #                                             n_grid_points = npt)
        # t_night.append([t_dusk0.jd, t_dawn0.jd])


        ### HORIZON ---
        # Get first period above horizon : can be the present period...
        high = obs.target_is_up(self.grb.t_trig, target,
                                horizon = self.altmin)
        if (high): search="previous"
        else: search = "next"
        t_rise = obs.target_rise_time(self.grb.t_trig,
                                      target,
                                      which   = search,
                                      horizon = self.altmin,
                                      n_grid_points = npt)

        # If rise time is undefined, this means that the GRB is always above
        # or below the horizon - Otherwise the set time can be found.
        if (math.isnan(t_rise.jd)):
            if (high):
                #t_above.append([-math.inf, math.inf])
                t_above = t_grb
                self.vis = True
            else:
                self.vis = False
                return
        else:
            self.vis = True
            t_set = obs.target_set_time(t_rise,
                                        target,
                                        which="next",
                                        horizon = self.altmin,
                                        n_grid_points = npt)
            t_above.append([t_rise.jd,t_set.jd])


            # Add a subsequent above-horizon periods until GRB end of data
            while (t_rise.jd < t_grb[0][1]):
                t_rise = obs.target_rise_time(t_set,
                                          target,
                                          which="next",
                                          horizon = self.altmin,
                                          n_grid_points = npt)
                t_set = obs.target_set_time(t_rise,
                                            target,
                                            which="next",
                                            horizon = self.altmin,
                                            n_grid_points = npt)
                t_above.append([t_rise.jd,t_set.jd])

        # Prompt appears above horiozn during night
        if (high and night): self.prompt = True

        # Now prepare the ticks from all the intervals
        ticks = []
        for elt in t_grb:   ticks.extend(elt)
        for elt in t_night: ticks.extend(elt)
        for elt in t_above: ticks.extend(elt)

        ticks.sort()
        if (debug): print("Number of ticks = ",len(ticks),ticks)

        # Loop over slices and check visibility
        if (debug):
            print(" {:<23s}   {:<23s} {:>10s} {:>6s} {:>6s} {:>6s}"
                  .format("T1", "T2", "bright", "dark", "above", "vis."))

        for i in range(len(ticks)-1):
            t1 = ticks[i]
            t2 = ticks[i+1]
            tmid = 0.5*(t1+t2)
            bright  = valid(tmid,t_grb)
            dark    = valid(tmid,t_night)
            above   = valid(tmid,t_above)
            visible = bright and dark and above
            if (visible):
                t_vis.append([t1, t2])
                self.tonight = True
                self.vis = True

            if (debug):
                if math.isinf(t1): t1 = "--"
                else : t1 = Df(t1).iso

                if math.isinf(t2): t2 = "--"
                else: t2 = Df(t2).iso
                print(" {:>23}   {:>23} {:>10} {:>6} {:>6} {:>6}"
                      .format(t1, t2,
                              bright, dark, above, visible),end="")
                if (visible): print(" *")
                else: print()

        # Write back all intervals into Time
        for elt in t_grb:
            self.t_grb.append( [Df(elt[0]), Df(elt[1])] )
        for elt in t_night:
            self.t_night.append( [Df(elt[0]), Df(elt[1])] )
        for elt in t_above:
            self.t_above.append( [Df(elt[0]), Df(elt[1])] )

        # Finalise visibility wondows, taking into account the depth
        # If no visibility window is left, re-assign vis_tonight
        for elt in t_vis:
            if (elt[end] - self.grb.t_trig.jd < self.depth.to(u.d).value):
                self.t_vis.append( [Df(elt[0]), Df(elt[1])] )

        if len(self.t_vis)== 0: self.tonight=False

        return

    ###------------------------------------------------------------------------
    def check(self, delta_max = 5., log=None):
        """
        Checks whether the visibility found in this code is
        identical to the one given by default for altmin=10°.
        The computation was done by maria Grazia Bernardini (MGB).
        See README file for more information.

        Parameters
        ----------
        delta_max : float, optional
            Max tolerance in seconds on window time limits.
            The default is 5..
        log : file pointer, optional
            Allows having the output both on screen and in a file. The default is None.

        """

        #----------------------------------------------------------------
        def status(tnew,torg,case=None):
            """
            This function is part of the check method only, and factorize a
            debugging printout

            Parameters
            ----------
            tnew : a time object pair in a list
                The visibility windows in this code.
            torg : a time object pair in a list
                The visibility windows in the original data.
            case : string, optional
                Some words describingthe data sets. The default is None.

            Returns
            -------
            bool
                Return True if both time set are identical within the tolerance
                delta_max.

            """
            delta1 = abs(tnew[0] -torg[0])
            delta2 = abs(tnew[1] -torg[1])
            if (delta1 > delta_max or delta2 > delta_max):

                if (log != None):
                    log.prt("   ==> {:>10s}: {} * {}"
                            .format(case,tnew[0].datetime,tnew[1].datetime))
                    log.prt("   ==> {:>10s}: {} * {}"
                          .format("ORG",torg[0].datetime,torg[1].datetime))
                    log.prt("   ==> {:>10s}: {:26.2f} * {:26.2f}"
                          .format("Delta",delta1,delta2))
                    log.prt("   ==> Not compatible with original values")
                    log.prt(" {:10s} {:15s} {:8.2f} {:8.2f} "
                          .format(self.name,case,delta1,delta2))
                return False
            else:
                return True
        #----------------------------------------------------------------

        matching = True
        log.prt(" *** Check {} {}".format(self.grb.name,self.loc),end=" ")

        # Check that general visibilities agree - if not return False
        if (self.vis != self.grb.vis[loc]):
            log.prt(" ==> Vis wrong ", end="")
            matching = False
        if (self.tonight != self.grb.vis_tonight[loc]):
            log.prt(" ==> Vis_tonight wrong ",end="")
            matching = False
        if (self.prompt != self.grb.vis_prompt[loc]):
            log.prt(" ==> Vis_prompt wrong ",end="")
            matching = False
        if len(self.t_vis) > 1:
            log.prt(" ==> >1 window",end="")
            matching = False

        if (not matching):
            print()
            return matching

        # If visibilities agree and there is nothing return True
        if (self.vis == False):
            log.prt ("--- not visible")
            return True
        if (self.tonight == False):
            log.prt ("--- not tonight")
            return True

        # both "tonight" are OK, compare delta-time

        if (self.tonight and self.grb.vis_tonight[loc]):
            # status(self.t_above,self.grb.t_event[loc],case = "Above")
            # status(self.t_night,self.grb.t_twilight[loc],case = "Night")
            matching = status(self.t_vis[0],self.grb.t_true[loc][0],case = "Vis")
            if matching: log.prt("--- ok")
            else: log.prt("--- DOES NOT MATCH")
        return matching

###---------------------------------------------------------------------
if __name__ == "__main__":
    """
    Compute the visibilities for some GRB and compare to the default stored in
    the GRB class. Update the GRB visibility windows.
    """

    from SoHAPPy import get_grb_fromfile
    from utilities import Log
    #from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib.pyplot as plt
    import grb_plot as gplt

    altmin     = 10*u.deg
    npt        = 150 # Change the default number of grid points in astroplan

    debug      = False
    showplot   = True # Plot visibility for all events

    check      = False # perform visibility window check within tolerance
    delta_max  = 5. # Tolerance in seconds
    depth      = 10*u.day # Keep windows up to that duration since the trigger (in days)

    logfilename = "visibility.log"  # Log file
    log = Log(name  = logfilename,talk=True)

    ### GRB list to be analysed ---
    # Events with Astroplan problem (risetime too close from tref)
    astroplan_north = [522, 644]
    astroplan_south = [233, 769]
    # Events for which MG misses a short window
    short_missed_North = [129, 174, 197, 214, 309, 380, 381, 479, 488, 502, 609,
                         673, 933, 965]
    short_missed_South = [24, 30, 182, 186, 319, 338, 444, 475, 734, 772, 826,
                         879, 917, 986]
    # Events with > 1 visibility window
    vis2_North = [75, 137, 139, 151, 153, 239, 270, 353, 507, 727, 859, 906,
                  923, 948]
    vis2_South = [85, 115, 414, 754, 755]

    # Consecutive list
    ifirst = 85
    ngrb = 1
    sites = ["South"]
    #sites   = ["North","South"]

    out = None
    if (check == True) and (ngrb > 10): ### If check done on many  GRB
        showplot = False
        # Enable checks if the minimal altiutde is the default
        out = log #open("deltas.txt", 'w')

    if type(ifirst)!=list: grblist = list(range(ifirst,ifirst+ngrb))
    else: grblist = ifirst

    # Let's go !
    for igrb in grblist:

        grb = get_grb_fromfile(igrb,log=None)

        # Check the update visibility process
        for loc in sites:

            # Compute new visibilities
            vis = Visibility(grb,loc=loc) # Constructor
            vis.compute(altmin=altmin,
                        end=1,
                        depth = depth,
                        debug=debug,npt=npt)

            # If requested check compatibility
            if (check):
                matching = vis.check(log=log)
            else:
                matching = True

            if (not matching) or (showplot):

                fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1,figsize=(12,12))

                print(" ORIGINAL")
                grb.print_visibility(loc=loc,log=log) # original visibility
                if (grb.vis[loc]):
                    gplt.visibility_plot(grb,ax=ax1,loc=loc)

                print(" NEW !!!")
                grb.print_visibility(loc=loc,log=log, alt=vis) # alt. visibility

                print(" UPDATED ")
                grb.update_visibility(vis)
                grb.print_visibility(loc=loc,log=log) # New visibility
                gplt.visibility_plot(grb,ax=ax2, loc=loc)

                # To be corrected
                # if (vis.t_vis[0][1] - vis.t_vis[0][0]).sec < 1000:
                #     t0 = 100*u.s
                #     t1 = +900*u.s
                #     axx =  inset_axes(ax, width="25%", height=1.2,
                #                       loc="upper right")
                #     vis.plot(axx,dt_before=t0, dt_after = t1, inset=True)
                plt.tight_layout()
                plt.show()
                if (check): fig.savefig("vis" + "_" + str(igrb) +"_"+loc+".pdf")

    if out: out.close()
