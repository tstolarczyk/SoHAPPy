# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:17:27 2020

This module is organised around the :class:`visibility`class.
It combines the rise, set and night (twilight) windows, and apply the moon veto
window to the GRB data window (defined by the two extreme GRB time of the
data points) to produce the visibility windows for a given GRB on a given site.
It can also be used to check that the ones given by default are in agreement
for the default horizon value (10 degrees) and for the nights found within
24 hours.

@author: Stolar
"""
import math
import astropy.units as u
from   astropy.time import Time
from   astropy.coordinates import AltAz, get_moon
from   astropy.table         import Table

from   astroplan import Observer, FixedTarget, moon_illumination

import warnings

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
    file, in particular to modify the minimum altitude for a source to be
    decently detectable.
    The method, :method:`check` is used to compare the visibility
    windows with the one given by default in the GRB file or any other
    pre-loaded visibility.
    """

    ###------------------------------------------------------------------------
    def __init__(self,target=None, site=None, tstart=0, tstop=0, name="Unknown"):
        """
        Visibility constructor
        I did not find how to have this displayed with automodapi

        Excerpt from the astroplan documentation:
        * Moonrise is defined as the instant when, in the eastern sky, under
        ideal meteorological conditions, with standard refraction of the
        Moon's rays, the upper edge of the Moon's disk is coincident with an
        ideal horizon.

        * Moonset is defined as the instant when, in the western sky, under
        ideal meteorological conditions, with standard refraction of the
        Moon's rays, the upper edge of the Moon's disk is coincident with
        an ideal horizon.

        The Moon angular diameter varies between 29.3 to 34.1 arcminutes
        (Wikipedia), namely 0.488 and 0.568 degree. The rise time is when
        the moon is below the horizon at a negative angle half of
        these values, respectively -0.244 and 0.284 degree

        Parameters
        ----------
        target : TYPE, optional
            DESCRIPTION. The default is None.
        site : TYPE, optional
            DESCRIPTION. The default is None.
        tstart : TYPE, optional
            DESCRIPTION. The default is 0.
        tstop : TYPE, optional
            DESCRIPTION. The default is 0.
        name : TYPE, optional
            DESCRIPTION. The default is "Unknown".

        Returns
        -------
        None.

        """


        self.status  = "init"

        self.site    = site
        self.target  = target
        self.tstart  = tstart
        self.tstop   = tstop
        self.name    = name
        self.depth   = 0

        self.altmin = 0*u.deg # GRB Minimal allowed altitude

        # Visible any moment of the year from the site
        self.vis         = True
        # Visible any moment within the 24h after the trigger from the site
        self.vis_tonight = True
        # Visible at the moment of the GRB trigger
        self.vis_prompt  = True

        # Visible above min. alt. - list of pair
        self.t_true     = [[Time('2000-01-01 00:00:00', scale='utc'),
                           Time('2000-01-01 08:00:00', scale='utc')]]

        # Astronomical twilight - list of pair
        self.t_twilight = [[Time('2000-01-01 00:00:00', scale='utc'),
                           Time('2000-01-01 8:00:00', scale='utc')]]

        # Rise and set of the target - list of pair
        self.t_event    = [[Time('2000-01-01 00:00:00', scale='utc'),
                           Time('2000-01-01 08:00:00', scale='utc')]]

        # Nights
        self.t_night  = [[]]

        # Moon period veto
        self.moon_maxalt     = -0.25*u.deg # Maximal allowed altitude
        self.moon_mindist    =  0*u.deg # Minimum distance to source
        self.moon_maxlight   =  1       # Maximum allowed brightness

        self.t_moon_up       = [[]] # Altitude exceeds moon_maxalt
        self.moon_too_bright = [] # Moon brigthness above threshold
        self.moon_too_close  = [] # Moon distance too small

        return

    ###-----------------------------------------------------------------------
    def from_fits(self, grb, hdr, hdul, hdu=1, loc="None"):
        """
        Visibility from input file
        The start and stop dates are searched during 24h after the trigger
        and correspond to the first visibility interval.
        Does not report a second visibility interval during 24h,
        that should be possible only if the target is promptly visible
        (about 15% of the targets)
        By default seen in the North.
        Default times are defined for the first GRB sample with no
        visibikity given. The North visibility is chosen to have the
        trigger therein. The South visibility starts after the trigger.
        This allows having the two conditions explored.

        Parameters
        ----------
        grb : TYPE
            DESCRIPTION.
        hdr : TYPE
            DESCRIPTION.
        hdul : TYPE
            DESCRIPTION.
        hdu : TYPE, optional
            DESCRIPTION. The default is 1.
        loc : TYPE, optional
            DESCRIPTION. The default is "None".

        Returns
        -------
        None.

        """

        self.status  = "Default"
        self.site    = grb.pos_site[loc]
        self.target  = FixedTarget(coord=grb.radec, name=grb.name)
        self.tstart  = grb.t_trig
        self.tstop   = grb.t_trig + grb.tval[-1]
        self.name    = grb.name+"_"+loc

        vis = Table.read(hdul,hdu=hdu)

        def f(t,loc):
            t = Time(t.data,format="jd",scale="utc")
            if (loc == "North"): return [ [ t[0:2][0], t[0:2][1] ] ]
            if (loc == "South"): return [ [ t[2:4][0], t[2:4][1] ] ]

         # Visibility has been computed with this minimum altitude
        if (loc == "North"):
            self.vis          = hdr['V_N_ANYT']
            self.vis_tonight  = hdr['V_N_TNG']
            self.vis_prompt   = hdr['V_N_PR']

        if (loc == "South"):
            self.vis          = hdr['V_S_ANYT']
            self.vis_tonight  = hdr['V_S_TNG']
            self.vis_prompt   = hdr['V_S_PR']

        self.altmin     = vis.meta["MIN_ALT"]*u.deg # Minimum altitude
        self.t_true     = f(vis["True"],loc)
        self.t_twilight = f(vis["Twilight"],loc)
        self.t_event    = f(vis["Event"],loc)

        #print(self.grb.name," : visibility read for ",loc)

        return

    ###-----------------------------------------------------------------------
    def write(self, folder=".",debug=False):
        """
        Write current visibility instance to disk

        Parameters
        ----------
        folder : TYPE, optional
            DESCRIPTION. The default is ".".
        debug : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """

        import pickle
        from pathlib import Path

        filename = Path(folder,self.name+"_vis.bin")
        outfile  = open(filename,"wb")
        pickle.dump(self,outfile)
        outfile.close()
        if (debug): print(" >>> Visibility written to : ",filename)
        return

    ###-----------------------------------------------------------------------
    @classmethod
    def read(cls,name,folder=".",debug=False):
        """
        Read a visibility instance from disk


        Parameters
        ----------
        cls : TYPE
            DESCRIPTION.
        name : TYPE
            DESCRIPTION.
        folder : TYPE, optional
            DESCRIPTION. The default is ".".
        debug : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        vis : TYPE
            DESCRIPTION.

        """

        import pickle
        from pathlib import Path

        filename = Path(folder,name)
        infile  = open(filename,"rb")
        vis =  pickle.load(infile)
        infile.close()
        if (debug): print(" <<<< Visibility read from : ",filename)

        return vis
    ###-----------------------------------------------------------------------
    def moonlight_veto(self,dt,debug=False):
        """
        Check if the Moon period defined by the rise and set time correspond
        to a situation where the moon is too bright or too close from the
        source.
        If this is the case (too bright or too close), returns True
        (the veto is confirmed)
        Criteria used by the GW paper group
        (F. Schussler, April 14th 2021, private communication) :
            max moon phase to 60%
            min moon separation 30deg
            max moon altitude 50 deg


        """
        too_bright = False
        too_close  = False

        # Check moon illumination
        for t in dt:
            moonlight = moon_illumination(Df(t))
            if (moonlight >= self.moon_maxlight):
                too_bright = True
                if (debug):
                    print("Moonlight :",moonlight," too bright ! -> confirmed")
                break

        # Check distance to source at rise and set
        for t in dt:
            moon_radec = get_moon(Df(t), self.site)
            dist = moon_radec.separation(self.target.coord)
            if dist <= self.moon_mindist:
                too_close = True
                if (debug): print(" Moon Distance : ",dist,"too close !")
                break

        return (too_bright, too_close)

    ###----------------------------------------------------------------------------
    def moon_halo(x, r0 = 0.5, rc=30, epsilon=0.1, q0 = 1):
        """
        Returns the moon light intensity (hand-made model) as a function of distance r
        It is assumed that the intensity is constant (q0) along the moon radius r0, and reach
        the fraction epsilon at the critical distance rc, following a 1/r2 function:
        1/(a*(x-r0)**2 +c)
        """
        f = 1/ ( (1 -epsilon)/(q0*epsilon)/(rc-r0)**2 * (x-r0)**2 + 1/q0)
        #f = 0.5*(1 - abs(r0-x)/(r0-x))* f + 0.5*(1 - abs(x-r0)/(x-r0))*q0
        return f

    ###----------------------------------------------------------------------------
    def moon_halo_veto(q0, threshold=0.1, r0 = 0.5, rc=30, epsilon=0.1):
        """
        Returns the distance at which the moon intensity reaches the threshold
        This is the analytical inverse of the moon_halo function
        """
        a = (1 -epsilon)/epsilon/q0/(rc-r0)**2
        return math.sqrt((1/threshold - 1/q0)/a)

    ###-----------------------------------------------------------------------
    def moon_alt_veto(self, obs, npt=150):
        """
        The first veto is the moon alitude.
        First find windows where the Moon is too high

        Parameters
        ----------
        obs : TYPE
            DESCRIPTION.
        twindow : TYPE
            DESCRIPTION.
        altmin : TYPE, optional
            DESCRIPTION. The default is 10*u.deg.
        npt : TYPE, optional
            DESCRIPTION. The default is 150.

        Returns
        -------
        tmoons : TYPE
            DESCRIPTION.

        """

        tmoons = []

        # Is the Moon there at trigger tigger time ?
        radec = get_moon(self.tstart,obs.location)
        altaz = radec.transform_to(AltAz(location=obs.location,
                                         obstime=self.tstart))

        # Search next rise except if Moon is already here
        search="next"
        if  (altaz.alt > self.moon_maxalt): search="previous"

        t_rise = obs.moon_rise_time(self.tstart,
                                    which = search, horizon=self.moon_maxalt,
                                    n_grid_points = npt)
        if (math.isnan(t_rise.jd)):
            # Moon will never rise
            return [[]] # No veto period

        t_set  = obs.moon_set_time(t_rise,
                                   which="next", horizon=self.moon_maxalt,
                                   n_grid_points = npt)

        tmoons.append([t_rise.jd, t_set.jd])

        # Add subsequent nights until reaching the end of GRB data
        while (t_set.jd < self.tstop.jd):
            t_rise = obs.moon_rise_time(t_set,
                                        which = "next", horizon=self.moon_maxalt,
                                        n_grid_points = npt)
            t_set  = obs.moon_set_time(t_rise,
                                       which = "next",  horizon=self.moon_maxalt,
                                       n_grid_points = npt)

            tmoons.append([t_rise.jd, t_set.jd])

        if len(tmoons): return tmoons
        else: return [[]]

        return tmoons
    ###-----------------------------------------------------------------------
    def nights(self,obs, skip=0, npt=150):
        """


        Parameters
        ----------
        obs : TYPE
            DESCRIPTION.
        is_night : TYPE
            DESCRIPTION.
        skip : TYPE, optional
            DESCRIPTION. The default is 0.
        tstop : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        tnights : TYPE
            DESCRIPTION.

        """

        tnights = []

        # Get the first night : can be the current night
        is_night = obs.is_night(self.tstart, horizon = -18*u.deg)

        if (is_night):
            search="previous"
        else:
            search = "next"

        inight = 0
        t_dusk = obs.twilight_evening_astronomical(self.tstart,
                                                    which = search,
                                                    n_grid_points = npt)
        t_dawn = obs.twilight_morning_astronomical(t_dusk,
                                                    which="next",
                                                    n_grid_points = npt)
        # Omit first night if requested
        if (skip ==0): tnights.append([t_dusk.jd, t_dawn.jd])

        # Add subsequent nights until reaching the end of GRB data
        while (t_dusk < self.tstop):
            inight +=1
            t_dusk = obs.twilight_evening_astronomical(t_dawn,
                                                        which = "next",
                                                        n_grid_points = npt)
            t_dawn = obs.twilight_morning_astronomical(t_dusk,
                                                        which = "next",
                                                        n_grid_points = npt)
            if (inight >= skip): tnights.append([t_dusk.jd, t_dawn.jd])

        # For test, add the previous night
        # t_dawn0 = obs.twilight_morning_astronomical(self.grb.t_trig,
        #                                            which="previous",
        #                                            n_grid_points = npt)
        # t_dusk0 = obs.twilight_evening_astronomical(t_dawn0,
        #                                             which="previous",
        #                                             n_grid_points = npt)
        # tnights.append([t_dusk0.jd, t_dawn0.jd])

        return (is_night, tnights)


    ###-----------------------------------------------------------------------
    def horizon(self, obs, npt=150):

        t_above = []

        # Get first period above horizon : can be the present period...
        high = obs.target_is_up(self.tstart, self.target,
                                horizon = self.altmin)

        if (high): search="previous"
        else: search = "next"

        t_rise = obs.target_rise_time(self.tstart,
                                      self.target,
                                      which   = search,
                                      horizon = self.altmin,
                                      n_grid_points = npt)

        # If rise time is undefined, this means that the GRB is always above
        # or below the horizon - Otherwise the set time can be found.
        if (math.isnan(t_rise.jd)):
            if (high):
                self.vis = True
                return high, [[self.tstart.jd,self.tstop.jd]]
            else:
                self.vis = False
                return high, [[]]
        else:
            self.vis = True
            t_set = obs.target_set_time(t_rise,
                                        self.target,
                                        which="next",
                                        horizon = self.altmin,
                                        n_grid_points = npt)
            t_above.append([t_rise.jd,t_set.jd])


            # Add a subsequent above-horizon periods until GRB end of data
            while (t_rise < self.tstop):
                t_rise = obs.target_rise_time(t_set,
                                          self.target,
                                          which="next",
                                          horizon = self.altmin,
                                          n_grid_points = npt)
                t_set = obs.target_set_time(t_rise,
                                            self.target,
                                            which="next",
                                            horizon = self.altmin,
                                            n_grid_points = npt)
                t_above.append([t_rise.jd,t_set.jd])

        return (high, t_above)

    ###-----------------------------------------------------------------------
    def compute(self,altmin     = 10*u.degree,
                     altmoon    = 0*u.degree,
                     moondist   = 0*u.degree,
                     moonlight  = 1,
                     depth      = 3*u.day,
                     end_of_day = 0,
                     skip       = 0,
                     npt        = 150,
                     debug=False):

        """
        Compute the visibility periods for a given GRB and site.
        The constructor takes as arguments a grb (GammaRayBurst) and a
        location (string).
        # Temporary local variables in Julian days (to use sort)
        # Will be copied as Time objects later into the class variable
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
            if len(tslices[0]) == 0 : return False # No slice !
            ok = False
            for slices in tslices:
                if (t0 >= slices[0] and t0 <= slices[1]):
                    if (ok == False): ok = True
            return ok
        ###---------------------------------------------------

        self.depth     = depth
        self.altmin    = altmin

        self.moon_maxalt   = altmoon
        self.moon_mindist  = moondist
        self.moon_maxlight = moonlight

        self.status    = "Updated"

        obs  = Observer(location  = self.site, name = self.name, timezone ="utc")

        ### NIGHT ---
        is_night, t_night  = self.nights(obs, skip=skip, npt=npt)

        ### MOON VETOES (high enough, close enough, bright enough) ---
        t_moon_alt_veto    = self.moon_alt_veto(obs, npt=npt)

        t_moon_veto = []
        for dt in t_moon_alt_veto:
            (too_bright, too_close) = self.moonlight_veto(dt)
            self.moon_too_bright.append(too_bright)
            self.moon_too_close.append(too_close)
            if too_bright or too_close: t_moon_veto.append(dt)
        if len(t_moon_veto) == 0: t_moon_veto = [[]]

        ### HORIZON ---
        (high, t_above) = self.horizon(obs)

        # Prompt appears above horizon during night
        if (high and is_night): self.vis_prompt = True

        # Now prepare the ticks from all the intervals
        ticks = [self.tstart.jd,self.tstop.jd]
        for elt in t_night       : ticks.extend(elt)
        for elt in t_above       : ticks.extend(elt)
        for elt in t_moon_veto   : ticks.extend(elt)

        ticks.sort() # Requires numerical values -> all times are in jd
        if (debug): print("Number of ticks = ",len(ticks),ticks)

        # Loop over slices and check visibility
        if (debug):
            print(" {:<23s}   {:<23s} {:>10s} {:>6s} {:>6s} {:>6s} {:>6s}"
                  .format("T1", "T2", "bright", "dark", "above", "moon.","vis."))

        t_vis   = []
        for i in range(len(ticks)-1):
            t1 = ticks[i]
            t2 = ticks[i+1]
            tmid = 0.5*(t1+t2)
            bright  = valid(tmid,[[self.tstart.jd,self.tstop.jd]])
            dark    = valid(tmid,t_night)
            above   = valid(tmid,t_above)
            moon    = not valid(tmid,t_moon_veto) # Moon acts as a veto
            visible = bright and dark and above and moon
            if (visible):
                t_vis.append([t1, t2])
                self.vis_tonight = True
                self.vis = True

            if (debug):
                if math.isinf(t1): t1 = "--"
                else : t1 = Df(t1).iso

                if math.isinf(t2): t2 = "--"
                else: t2 = Df(t2).iso
                print(" {:>23}   {:>23} {:>10} {:>6} {:>6} {:>6} {:>6}"
                      .format(t1, t2,
                              bright, dark, above, not moon, visible),end="")
                if (visible): print(" *")
                else: print()

        # Write back all intervals into Time and into Class
        if len(t_night[0])==0 :
            self.t_twilight  = [[]]
        else:
            self.t_twilight  = []
            for elt in t_night:
                self.t_twilight.append( [Df(elt[0]), Df(elt[1])] )

        if len(t_above[0])==0 :
            self.t_event  = [[]]
        else:
            self.t_event  = []
            for elt in t_above:
                self.t_event.append( [Df(elt[0]), Df(elt[1])] )

        if len(t_moon_alt_veto[0])==0 :
            self.t_moon_up  = [[]]
        else:
            self.t_moon_up  = []
            for elt in t_moon_alt_veto:
                self.t_moon_up.append( [Df(elt[0]), Df(elt[1])] )

        # Finalise visibility wondows, taking into account the depth
        # and additionnal moon vetoes.
        # If no visibility window is left, re-assign vis_tonight
        if len(t_vis) == 0 :
            self.t_true  = [[]]
            self.vis_tonight = False
            self.vis         = False
        else:
            ### At least one visibility period is found
            # Will it survive moon distance and brigthness cut ?
            self.vis_tonight = True
            self.t_true  = []
            for elt in t_vis:
#                if not moonlight_veto(t_vis):
                #self.t_moon_alt.append( [Df(elt[0]), Df(elt[1])] )
                if (elt[end_of_day] - self.tstart.jd <= self.depth.to(u.d).value):
                    self.t_true.append( [Df(elt[0]), Df(elt[1])] )

        return

    ###------------------------------------------------------------------------
    def print(self,log=None, alt=None):
        """
        Print out the GRB night, above-the-horizon periods, and default
        visibility window.

        Parameters
        ----------
        loc : String, optional
            Either "North" or "South". The default is None.
        log : TextIO, optional
            Log file. The default is None.
        alt : Visibility instance
            Get the values from the Visibility object (used for comparisons)

        Returns
        -------
        None.

        """

        log.prt("=================   {:10s}   {:10s}   ================"
          .format(self.name,self.status))
        log.prt(' Visible : {} - tonight, prompt : {}, {})'
              .format(self.vis, self.vis_tonight,self.vis_prompt))
        log.prt(' Altitude : Horizon > {:3.1f} - Moon > {:4.2f}'
              .format(self.altmin, self.moon_maxalt))
        #log.prt(" Trigger: {}".format(self.t_trig.datetime))

        if self.vis_tonight or self.status =="Updated": # Seen within 24hr after the trigger
            log.prt("+----------------------------------------------------------------+")

            ###------------------
            def show(t,case="Dummy"):
                if len(t[0]) == 0:
                    log.prt(" {:6s} : {:26s} * {:26s}".format(case,"--","--"))
                    return
                for i, elt in enumerate(t):
                    t1 = elt[0]
                    t2 = elt[1]
                    log.prt(" {:6s} : {} * {}".format(case,
                                                      t1.datetime,
                                                      t2.datetime),end="")
                    # t1  = (t1-self.grb.t_trig).sec*u.s
                    # t2  = (t2-self.grb.t_trig).sec*u.s
                    # log.prt("        : {:7.2f} {:6.2f} * {:7.2f} {:6.2f}"
                    #       .format(t1,self.grb.altaz(loc=loc,dt=t1).alt,
                    #               t2,self.grb.altaz(loc=loc,dt=t2).alt))
                    if case=="Moon":
                        log.prt(" B:{} D:{}"
                                .format(str(self.moon_too_bright[i])[0],
                                        str(self.moon_too_close[i])[0]))
                    else:
                        log.prt("")
                return
            #-------------------

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                show(self.t_event,case="Event") # Event - above altmin
                show(self.t_twilight,case="Twil.")  # Twilight - Dark time
                show(self.t_moon_up,case="Moon")  # Moon altitude veto
                show(self.t_true,case="True")  # True : dark time + altmin + triggered
            log.prt("+----------------------------------------------------------------+")

        return
    ###------------------------------------------------------------------------
    def check(self, view, loc=None, delta_max = 5*u.s, log=None):
        """
        Checks whether the visibility found in this code is
        identical to the one in argument.
        In particular, allow checking the default visbility for altmin=10°
        and one night (The computation was done by maria Grazia Bernardini (MGB).
        See README file for more information.

        Parameters
        ----------
        vis: Visibility
            The visibility to be compared wiht
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
                          .format("was",torg[0].datetime,torg[1].datetime))
                    log.prt("   ==> {:>10s}: {:26.2f} * {:26.2f}"
                          .format("Delta",delta1.sec,delta2.sec))
                    log.prt("   ==> Not compatible with original values")

                return False
            else:
                return True
        #----------------------------------------------------------------

        matching = True
        log.prt(" *** Check {} {}, '{}' and '{}' with tolerance : {}"
                .format(self.grb.name,self.loc,self.name,view.name,delta_max),end=" ")

        # Check that general visibilities agree - if not return False
        if (self.vis != view.vis):
            log.prt(" ==> Vis wrong ", end="")
            matching = False
        if (self.vis_tonight != view.vis_tonight):
            log.prt(" ==> Vis_tonight wrong ",end="")
            matching = False
        if (self.vis_prompt != view.vis_prompt):
            log.prt(" ==> Vis_prompt wrong ",end="")
            matching = False
        if len(self.t_true) > 1:
            log.prt(" ==> >1 window",end="")
            matching = False

        if (not matching):
            print()
            return matching

        # If visibilities agree and there is nothing to compare with, return True
        if (self.vis == False):
             log.prt ("--- not visible")
             return True
        if (self.tonight == False):
             log.prt ("--- not tonight")
             return True

        # both "tonight" are OK, compare delta-time

        if (self.tonight and view.vis_tonight):
            print()
            # status(self.t_above,self.grb.t_event,case = "Above")
            # status(self.t_night,self.grb.t_twilight,case = "Night")
            matching = status(self.t_true[0],view.t_true[0],case = "Vis")
            if matching: log.prt("--- ok")
            else: log.prt("--- DOES NOT MATCH")
        return matching

###---------------------------------------------------------------------
if __name__ == "__main__":
    """
    Compute the visibilities for some GRB and compare to the default stored in
    the GRB class. Update the GRB visibility windows.

    * Events with Astroplan problem (risetime too close from tref)
    astroplan_north = [522, 644]
    astroplan_south = [233, 769]

    * Events for which MG misses a short window
    short_missed_North = [129, 174, 197, 214, 309, 380, 381, 479, 488, 502, 609,
                         673, 933, 965]
    short_missed_South = [24, 30, 182, 186, 319, 338, 444, 475, 734, 772, 826,
                         879, 917, 986]

    * Events with > 1 visibility window
    vis2_North = [75, 137, 139, 151, 153, 239, 270, 353, 507, 727, 859, 906,
                  923, 948]
    vis2_South = [85, 115, 414, 754, 755]

    """

    from   utilities import Log
    from   SoHAPPy import get_grb_fromfile
    import grb_plot as gplt

    newvis    = True
    altmin    = 24*u.deg
    altmoon   = -0.25*u.deg
    moondist  = 30*u.deg
    moonlight = 0.6
    depth     = 3*u.day
    skip      = 0

    # GRB list to be analysed
    # 751 and 815 are killed by the Moon brigthness
    ifirst = [54, 85, 815, 751]
    ifirst = 51
    ngrb   = 50

    dbg   = 1
    show  = 1

    res_dir = './'

    log_filename    = res_dir  + "visibility.log"
    log = Log(name  = log_filename, talk=True)

    check = False # If True,run comparison between old and new visibilities
    if (check):
        altmin    = 10*u.deg
        moondist  = 0*u.deg
        moonlight = 1.
        depth     = 1*u.day
        skip      = 0

    if type(ifirst)!=list:
        grblist = list(range(ifirst,ifirst+ngrb))
    else: grblist = ifirst

    # Loop over GRB list accordng to the config. file
    for i in grblist:

        grb = get_grb_fromfile(i,log=log)
        print(grb)

        for loc in ["North","South"]:

            # Current visibility - not exciting
            log.prt("\n--------------- D E F A U L T  -----------------\n")
            grb.vis[loc].print(log=log)
            if (check and ngrb<30 and show>0):
                gplt.visibility_plot(grb, loc=loc)
            # Save old visibilities
            import copy
            vis_def = copy.deepcopy(grb.vis[loc])

            if (newvis):
            # Compute visibility with new altitude
                log.prt("\n--------------- New visibility  -----------------\n")
                grb.vis[loc].compute(altmin    = altmin,
                                     altmoon   = altmoon,
                                     moondist  = moondist,
                                     moonlight = moonlight,
                                     depth     = depth,
                                     skip      = skip,
                                     debug     = False)

                grb.vis[loc].print(log=log)

                if ngrb<30 and show>0:
                    gplt.visibility_plot(grb, loc=loc)

            if (check): # If requested, perform check
                log.prt()("\n--------------- C H E C K  ----------------\n")
                grb.vis[loc].check(vis_def,log=log,delta_max=1*u.s)

    log.close()

    print(" All done !")