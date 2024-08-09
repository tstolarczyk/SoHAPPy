# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:17:27 2020

This module is organised around the :class:`visibility` class.
It combines the rise, set and night (twilight) windows, and apply the Moon veto
window to the GRB data window (defined by the two extreme GRB time of the
data points) to produce the visibility windows for a given GRB on a given site.
It can also be used to check that the ones given by default are in agreement
for the default horizon value (10 degrees) and for the nights found within
24 hours.

@author: Stolar

"""

import warnings

import sys
import json
import yaml
from yaml import SafeLoader

import numpy as np
from pathlib import Path

import astropy.units as u
from   astropy.time import Time
from   astropy.coordinates import AltAz, SkyCoord, get_moon, EarthLocation
from   astropy.table import Table
from   astropy.io import fits

from observatory import xyz as obs_loc
from niceprint import Log
from utilities import get_filename, Df, Dp

from astroplan import Observer, FixedTarget, moon_illumination
from moon import moon_alt_plot, moonlight_plot, moon_dist_plot  #, moonphase_plot

__all__ = ["Visibility"]

###############################################################################
class Visibility():

    """
    This class defines the visibility periods for a given GRB and a certain
    site (North or South).
    The method, :class:`check` is used to compare the visibility
    windows with the one given by default in the GRB file or any other
    pre-loaded visibility.

    Excerpt from the astroplan documentation:

    * Moonrise is defined as the instant when, in the eastern sky, under ideal
      meteorological conditions, with standard refraction of the Moon's rays,
      the upper edge of the Moon's disk is coincident with an ideal horizon.
    * Moonset is defined as the instant when, in the western sky, under ideal
      meteorological conditions, with standard refraction of the Moon's rays,
      the upper edge of the Moon's disk is coincident with an ideal horizon.

    The Moon angular diameter varies between 29.3 to 34.1 arcminutes
    (Wikipedia), namely 0.488 and 0.568 degree. The rise time is when the moon
    is below the horizon at a negative angle half of these values, respectively
    -0.244 and -0.284 degree

    """
    ignore = []
    """ Members to be ignored when exported to Json. """

    def_vis_dicts  = "visibility.yaml"
    """ default visibility parameter file """

    ###------------------------------------------------------------------------
    def __init__(self, pos    = SkyCoord(0*u.deg,0*u.deg, frame='icrs'),
                       site   = None,
                       window = [0,0],
                       status = "Empty", name="Dummy"):
        """
        Create a Visibility instance with default values.

        Parameters
        ----------
        pos : Astropy SkyCoord, optional
            Source postion in the sky.
            The default is SkyCoord(0*u.deg,0*u.deg, frame='icrs').
        site : string, optional
            A keyword describing the site (e.g. `North`). The default is None.
        window : python list of float, optional
            A list with the start and stop time of the source.
            The default is [0,0].
        status : string, optional
            A keyword cahracterising this visibility. The default is "Empty".
        name : string, optional
            A name for this visibility. The default is "Dummy".

        Returns
        -------
        None.

        """

        self.status  = status # Status of the instance, e.g. `recomputed`
        self.site    = site
        self.target  = FixedTarget(coord=pos, name="Source")
        self.tstart  = window[0]
        self.tstop   = window[1]
        self.name    = name

        # GRB Minimal allowed altitude
        self.altmin  = 10*u.deg

        # Moon period veto - Default is no veto
        self.moon_maxalt     = 90*u.deg # Maximal allowed altitude
        self.moon_mindist    =  0*u.deg # Minimum distance to source
        self.moon_maxlight   =  1       # Maximum allowed brightness

        self.vis         = True
        # These three variables were defined by M.G.Bernardini in the first
        # version of the visibilty, without moon veto, and with a visibility
        # computed for the first night only. They are kept for consistency but
        # have a somewhat smaller importance with the presence of moon vetoes
        # and a visibility that can be computed on more than one night.
        # However the vis_prompt variable still indicates the chance to have
        # the prompt emission detected if the alert and slewing times are not
        # too long.
        # Here below, better definitions of these variables are given.

        # Originally: Visible any moment of the year from the site
        # This means that the GRB appears above the defined horizon whathever
        # the other conditions are. In particular, mitigation against the
        # Moonlight would make it visible from the first or next nights.

        # Originally: Visible any moment within the 24h after the trigger from
        # the site, names vis_tonight in the GRB file.
        # Changed to vis_night and True if the GRB is visible during any night
        # This the variable on which the decision of simulation is given.
        self.vis_night = True

        # Originally: Visible at the moment of the GRB trigger
        # Keeps the original meaning but add the moon veto, which can differ
        # with the assumprions made with the affordable Moonlight.
        self.vis_prompt  = True

        # Visible above min. alt. - list of pair
        self.t_true     = [[self.tstart, self.tstop]]

        # Astronomical twilight - list of pair
        self.t_twilight = [[self.tstart, self.tstop]]

        # Rise and set of the target - list of pair
        self.t_event    = [[self.tstart,self.tstop]]

        self.t_moon_up       = [[]] # Altitude exceeds moon_maxalt
        self.moon_too_bright = [] # Moon brigthness above threshold
        self.moon_too_close  = [] # Moon distance too small

        self.depth = 3 # Number of nights to be considered
        self.skip  = 0 # Number of first nights to be skipped

        return

    ###-----------------------------------------------------------------------
    def compute(self, param = None,
                      npt   = 150,
                      debug = False):

        """
        Compute the visibility periods until the end of the last night within
        the data window. This function takes as arguments the parameters
        required for the computation obtained from the `visibility.yaml` local
        file through the :func:`Configuration.decode_visibility_keyword`
        function.

        The algorithm is the following:

        1. Build periods for:
            * Nights (*t_night*), and check if trigger occurs during night
              (*is_night*);
            * Potential moon vetoes (*t_moon_alt_veto*) due to the Moon above
              the defined horizon (*altmoon*);
            * Check (*True*, *False*) whether these Moon veto periods are kept
              by the moon brightness and distance to the source
              (*moonlight_veto*). If this is the case, keep the corresponding
              moon veto period in *t_moon_veto* otherwise discard the period.
              In order to have all Moon periods kept (i.e. veto is valid as soon
              as Moon is above *altmoon*), *moondist* should be maximised
              (veto at any distance) and *moonlight* should be minimised (even
              if no moonligt, veto is valid);
            * Above horizon (*t_above*) and whether the object is above the
              horizon (*altmin*) at trigger time.

        2. Put all *start* and *stop* of all these periods in a list, sort the
           list, resulting in a series of ordered “ticks”.

        3. For each consecutive tick pairs in the list:
            * Compute the mean time (middle of the a tick pair)
            * Check if that time belongs to one of the *night*, *above*, or \
              *bright* (i.e. in the available data interval) intervals and not
              to any  *moon* interval. If so get True, otherwise False.

        4. For each of the tick pairs compute the Boolean `visible = bright
           and dark and above and moon` and get the visibility windows when True.

        Parameters
        ----------
        self : Visibility Object
            The present instance.
        param : Dictionnary, optional
            A dictionnary of parameters to compute the visibility.
            The default is None.
        npt : integer, optional
            Number of grid points for horizon crossing. The default is 150.
        debug : bool, optional
            Print additional comments at excecution time if True .
            The default is False.

        Returns
        -------
        Visibility Object
            The updated instance.

        """

        # Decode the parameter dictionnay - keep default otherwise
        if param is not None:
            # observatory = param["where"]
            self.altmin        = u.Quantity(param["altmin"])
            self.moon_maxalt   = u.Quantity(param["altmoon"])
            self.moon_mindist  = u.Quantity(param["moondist"])
            self.moon_maxlight = param["moonlight"]
            self.depth         = param["depth"]
            self.skip          = param["skip"]

        self.status = "Computed"

        obs  = Observer(location= self.site, name= self.name, timezone="utc")

        self.vis = False
        self.vis_night = False
        self.vis_prompt  = False

        ###---------------------
        ### Find the nights  ---
        ###---------------------
        is_night, self.t_twilight  = self.nights(obs, npt=npt)
        if len(self.t_twilight) ==0:
            self.t_true  = [[]]
            return self

        ###---------------------
        ### MOON VETOES (high enough, close enough, bright enough) ---
        ###---------------------
        # These are the periods when the Moon is above horizon
        self.t_moon_up    = self.moon_alt_veto(obs, npt=npt)

        # When Moon is up, check if moonlight is affordable
        t_moon_veto = [] # This variable is internal !!!
        for dt in self.t_moon_up:
            (too_bright, too_close) = self.moonlight_veto(dt)
            self.moon_too_bright.append(too_bright)
            self.moon_too_close.append(too_close)
            # If the Moon being above the horizon it gives too much
            # light, add the corresponding period to the Moon veto
            if too_bright or too_close: t_moon_veto.append(dt)
        # In case there is no Moon veto period, then conclude with
        # an empty array of [t1,t2] arrays
        if len(t_moon_veto) == 0: t_moon_veto = [[]]

        ###---------------------
        ### HORIZON ---
        ###---------------------
        (high, self.t_event) = self.horizon(obs)

        ###---------------------
        ### Collect the ticks from all periods (authorised and forbidden)
        ###---------------------
        # Note : ticks are in mjd (float) to use the sort function
        # First tick is the start of data
        ticks = [self.tstart.mjd]

        # Restrict the windows to the last night end or the GRB data length
        if self.tstop < self.t_twilight[-1][1]:
            # The GRB data stop before the end of last night
            # print(" >>>> Analysis shortened by lack of data")
            ticks.extend([self.tstop.mjd])
        else:
            # Set the end at last night for convenience (Time)
            self.tstop = self.t_twilight[-1][1]

        # Add the ticks of all previously computed windows
        for elt in self.t_twilight :
            ticks.extend([t.mjd for t in elt])
        for elt in self.t_event :
            ticks.extend([t.mjd for t in elt])
        for elt in t_moon_veto   : # And not only t_moon_alt_veto!
            ticks.extend([t.mjd for t in elt])

        # Sort in time
        ticks.sort() # Requires numerical values -> all times are in mjd

        if debug:
            print("Ticks : ",len(ticks))
            for t in ticks:
                print("{:10s} {:<23s} ".format(self.name, Df(t).iso))

        ###---------------------
        ### Check the visibility within all tick intervals
        ###---------------------
        # Loop over slices and check visibility
        if debug: # Header
            print(" {:<23s}   {:<23s} {:>10s} {:>6s} {:>6s} {:>6s} {:>6s}"
                  .format("T1", "T2", "bright", "dark", "above", "moon.","vis."))

        # t_vis   = []
        self.t_true   = []
        for i in range(len(ticks)-1):
            t1 = ticks[i]
            t2 = ticks[i+1]
            tmid = 0.5*(t1+t2)
            # Check if the current tick interval corresponds to an
            # authorised or forbidden window

            # The GRB shines (in the limit of the available data
            bright = self.valid(tmid,[[self.tstart,self.tstop]])

            # It is night
            dark  = self.valid(tmid,self.t_twilight)

            # It is above the horizon
            above = self.valid(tmid,self.t_event)

            # The moon authorises the observation
            not_moon = not self.valid(tmid,t_moon_veto) # Moon vetoes

            # In the end the GRB is visible if all conditions are fulfilled
            visible = bright and dark and above and not_moon
            # if visible: t_vis.append([t1, t2])
            if visible: self.t_true.append([ Time(t1,format="mjd"),
                                             Time(t2,format="mjd") ] )
            if debug:
                if np.isinf(t1):
                    t1 = "--"
                else :
                    t1 = Df(t1).iso

                if np.isinf(t2):
                    t2 = "--"
                else:
                    t2 = Df(t2).iso
                print(" {:>23}   {:>23} {:>10} {:>6} {:>6} {:>6} {:>6}"
                      .format(t1, t2,
                              bright, dark, above, not not_moon, visible),
                      end="")
                if visible:
                    print(" *")
                else:
                    print()

        ###---------------------------------------
        ### Finalise the visibility windows and flags
        ###---------------------------------------

        # At least a period above the horizon
        # Note that this function could stop as soon as no window above
        # the horizon is found, which is not the case in this implementation
        # if len(t_above) > 0: self.vis = True
        if len(self.t_event) > 0:
            self.vis = True

        # At least a visibility period for observation
        # if len(t_vis) > 0:
        if len(self.t_true) > 0:
            self.vis_night = True

            # In this first window the prompt is visible
            # Note that tstart is considered to be grb.t_trig
            if  self.valid(self.tstart.mjd,[self.t_true[0]]):
                self.vis_prompt=True

        # There is no visibility window at all - sad World
        else:
            self.t_true  = [[]]

        return self

    ###----------------------------------------------------------------------------
    @staticmethod
    def valid(t0,tslices):

        """
        Check it t0 in MJD is within the boundaries of tslices given as an
        array of two Time objects.

        Parameters
        ----------
        t0 : float
            A time in days.
        tslices : list
            A list of Astropy Time object pair.

        Returns
        -------
        ok : boolean
            return True if t0 within a slice.

        """
        if len(tslices[0]) == 0 : return False # No slice !
        ok = False
        for slices in tslices:
            if t0 >= slices[0].mjd and t0 <= slices[1].mjd:
                if ok == False: ok = True
        return ok

    ###------------------------------------------------------------------------
    def force_night(self, altmin = 24*u.deg, depth = 3, npt = 150):

        """
        From an existing instance, compute the visibility assuming the night \
        is of infinite length (i.e. from the trigger time to the end of the \
        data window with a safety margin of 1/10 of a day). Keep the vetoes \
        from the Moon and the altitude (above horizon).

        Parameters
        ----------
        altmin : Qauntity degree, optional
            Minimal altitude. The default is 24*u.deg.
        depth : integer, optional
            The number of nights considered. The default is 3.
        npt : integer, optional
            `astroplan` number of search points. The default is 150.
        debug : boolean, optional
            If True, display information. The default is False.

        Returns
        -------
        Visibility instance
            Updated visibility instance.

        """

        self.status    = "Forced"
        self.depth     = depth
        self.altmin    = altmin

        # Horizon (defines the visibiltiy)
        obs  = Observer(location  = self.site,
                        name      = self.name,
                        timezone  ="utc")

        (high, self.t_event) = self.horizon(obs)
        if len(self.t_event[0])==0 : # Above horizon period
            self.t_true = [[]]
            self.vis_night = False # Not visibile even during night
        else:
            self.vis = True
            self.vis_prompt = True
            self.t_event = [ [t[0], t[1]] for t in self.t_event]
            self.t_true = self.t_event
            self.vis_night = True

        # Infinite nights - starts at trigger (with margins)
        self.t_twilight  = [[self.tstart - 0.1*u.d,
                             self.tstop  + 0.1*u.d]] # Infini

        return self

    ###-----------------------------------------------------------------------
    def nights(self, obs, npt=150):

        """
        Return night periods withinthe observation window.

        Parameters
        ----------
        obs : Astroplan Observer instance
            Current observation window.
        npt: integer
            Number of points for window searches (Astroplan).
            The default is 150.

        Returns
        -------
        tnights : List of Time intervals
            List of nights.
        is_night: Boolean
            True if the observation starts at night.

        """

        tnights = []
        inight  = 0 # night (after trigger) counter

        # Get the first night : can be the current night
        is_night = obs.is_night(self.tstart, horizon = -18*u.deg)

        if is_night:
            search="previous"
        else:
            search = "next"

        t_dusk = obs.twilight_evening_astronomical(self.tstart,
                                                    which = search,
                                                    n_grid_points = npt)
        t_dawn = obs.twilight_morning_astronomical(t_dusk,
                                                    which="next",
                                                    n_grid_points = npt)
        # Omit first night if requested
        inight = 1 # night counter
        if self.skip == 0: tnights.append([t_dusk, t_dawn])

        # Add subsequent nights until reaching the end of GRB data
        while (t_dusk < self.tstop) and (inight < self.depth):
            t_dusk = obs.twilight_evening_astronomical(t_dawn,
                                                        which = "next",
                                                        n_grid_points = npt)
            t_dawn = obs.twilight_morning_astronomical(t_dusk,
                                                        which = "next",
                                                        n_grid_points = npt)
            if inight >= self.skip:
                tnights.append([t_dusk, t_dawn])

            inight +=1

        # For test, add the previous night
        # t_dawn0 = obs.twilight_morning_astronomical(self.grb.t_trig,
        #                                            which="previous",
        #                                            n_grid_points = npt)
        # t_dusk0 = obs.twilight_evening_astronomical(t_dawn0,
        #                                             which="previous",
        #                                             n_grid_points = npt)
        # tnights.append([t_dusk0.mjd, t_dawn0.mjd])

        if len(tnights) == 0:
            print(">>> No night found, please check your input parameters")

        return (is_night, tnights)

    ###-----------------------------------------------------------------------
    def horizon(self, obs, npt=150):

        """
        Compute periods above horizon within the observation.

        Parameters
        ----------
        obs : Astroplan Observer instance
            Current observation window.
        npt: integer
            Number of points for window searches (Astroplan).
            The default is 150.

        Returns
        -------
        t_above : List of Time intervals
            List of periods above horizon.
        high: Boolean
            True if the observation starts above horizon.

        """

        t_above = []

        # Get first period above horizon : can be the present period...
        high = obs.target_is_up(self.tstart, self.target,
                                horizon = self.altmin)

        if high:
            search="previous"
        else:
            search = "next"

        t_rise = obs.target_rise_time(self.tstart,
                                      self.target,
                                      which   = search,
                                      horizon = self.altmin,
                                      n_grid_points = npt)

        # If rise time is undefined, this means that the GRB is always above
        # or below the horizon - Otherwise the set time can be found.
        # Note that : if np.isnan(t_rise.mjd):
        # through a Warning. After investigation the test implemented below
        # seems correct (In numpy, Masked arrays, numpy.na, are arrays that may
        # have missing or invalid entries.)
        # if isinstance(t_rise.mjd ,np.ma.core.MaskedArray):
        # With gammapy1.2 installation, slight change in the test:
        if np.ma.is_masked(t_rise.mjd):
            if high:
                self.vis = True
                return high, [[self.tstart,self.tstop]]
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

            t_above.append([t_rise,t_set])


            # Add a subsequent above-horizon periods until GRB end of data
            while t_set < self.tstop:
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
                t_above.append([t_rise,t_set])

        return (high, t_above)

    ###-----------------------------------------------------------------------
    def moonlight_veto(self, dt, debug=False):

        """
        Check if the Moon period defined by the rise and set time correspond
        to a situation where the moon is too bright or too close from the
        source.
        If this is the case (too bright or too close), returns True
        (the veto is confirmed).

        Parameters
        ----------
        dt : Time interval
            Current time interval.
        debug : Boolean, optional
            If True, displays information. The default is False.

        Returns
        -------
        too_bright : Boolean
            True if Moon too bright.
        too_close : Boolean
            True if Moon too close.

        """

        too_bright = False
        too_close  = False

        # Check moon illumination
        for t in dt:
            moonlight = moon_illumination(t)
            if moonlight >= self.moon_maxlight:
                too_bright = True
                if debug:
                    print("Moonlight :",moonlight," too bright ! -> confirmed")
                break

        # Check distance to source at rise and set
        for t in dt:
            moon_radec = get_moon(t, self.site)
            dist = moon_radec.separation(self.target.coord)
            if dist <= self.moon_mindist:
                too_close = True
                if debug:
                    print(" Moon Distance : ",dist,"too close !")
                break

        return (too_bright, too_close)

    ###----------------------------------------------------------------------------
    def moon_halo(x, r0 = 0.5, rc=30, epsilon=0.1, q0 = 1):

        """
        Returns the moon light intensity (hand-made model) as a function of
        distance r. It is assumed that the intensity is constant (q0) along the
        moon radius r0, and reach the fraction epsilon at the critical
        distance rc, following a 1/r2 function 1/(a*(x-r0)**2 +c).

        Parameters
        ----------
        x : float
            Distance from the Moon.
        r0 : float, optional
            Moon radius. The default is 0.5.
        rc : float, optional
            Moon halo critical radius. The default is 30.
        epsilon : float, optional
            Affordable luminosity fraction. The default is 0.1.
        q0 : float, optional
            Moon light intensity (arbitrary unit). The default is 1.

        Returns
        -------
        f : float
            Moon halo intensity.

        """

        f = 1/ ( (1 -epsilon)/(q0*epsilon)/(rc-r0)**2 * (x-r0)**2 + 1/q0)
        #f = 0.5*(1 - abs(r0-x)/(r0-x))* f + 0.5*(1 - abs(x-r0)/(x-r0))*q0
        return f

    ###----------------------------------------------------------------------------
    def moon_halo_veto(q0, threshold=0.1, r0 = 0.5, rc=30, epsilon=0.1):

        """
        Returns the distance at which the moon intensity reaches the threshold
        This is the analytical inverse of the moon_halo function

        Parameters
        ----------
        q0 : float, optional
            Moon light intensity (arbitrary unit). The default is 1.
        threshold: float
            Affordable Moon light liminosity. The default is O.1.
        r0 : float, optional
            Moon radius. The default is 0.5.
        rc : float, optional
            Moon halo critical radius. The default is 30.
        epsilon : float, optional
            Affordable luminosity fraction. The default is 0.1.


        Returns
        -------
        float
            Radius with affordbale Moon light intensity.

        """

        a = (1 -epsilon)/epsilon/q0/(rc-r0)**2
        return np.sqrt((1/threshold - 1/q0)/a)

    ###-----------------------------------------------------------------------
    def moon_alt_veto(self, obs, npt=150):

        """
        The first veto is the Moon alitude.
        First find windows where the Moon is too high

        Parameters
        ----------
        obs : Astroplan Observer instance
            Current observation window.
        npt: integer
            Number of points for window searches (Astroplan).
            The default is 150.

        Returns
        -------
        tmoons : List of Time intervals
            Periods where the Moon light is not affordable for the observation.

        """

        tmoons = []

        # Is the Moon there at trigger tigger time ?
        radec = get_moon(self.tstart, obs.location)
        altaz = radec.transform_to(AltAz(location = obs.location,
                                         obstime  = self.tstart))

        # Search next rise except if Moon is already here
        search="next"
        if altaz.alt > self.moon_maxalt:
            search="previous"

        t_rise = obs.moon_rise_time(self.tstart,
                                    which = search, horizon=self.moon_maxalt,
                                    n_grid_points = npt)
        # if np.isnan(t_rise.mjd):
        if isinstance(t_rise.mjd ,np.ma.core.MaskedArray):
            # Moon will never rise
            print(" >>>>> Moon will never rise above ",self.moon_maxalt)
            return [[]] # No veto period

        t_set  = obs.moon_set_time(t_rise,
                                   which="next", horizon=self.moon_maxalt,
                                   n_grid_points = npt)

        tmoons.append([t_rise, t_set])

        # Add subsequent nights until reaching the end of GRB data
        while t_set < self.tstop:
            t_rise = obs.moon_rise_time(t_set,
                                        which = "next", horizon=self.moon_maxalt,
                                        n_grid_points = npt)
            t_set  = obs.moon_set_time(t_rise,
                                       which = "next",  horizon=self.moon_maxalt,
                                       n_grid_points = npt)

            tmoons.append([t_rise, t_set])

        if len(tmoons):
            return tmoons
        else:
            return [[]]

    ###------------------------------------------------------------------------
    @classmethod
    def from_fits(cls, grb, loc="None"):

        """
        Default visibility from input file.

        * The start and stop dates are searched during 24h after the trigger \
        and correspond to the first visibility interval.
        * Does not report a second visibility interval during 24h, that \
        should be possible only if the target is promptly visible (about 15% of the targets)

        Parameters
        ----------
        grb : GammaRayBurst instance
            The current source.
        loc : string, optional
            A string defining the site position. The default is "None".

        Returns
        -------
        None.

        """
        hdul   = fits.open(get_filename(grb.filename))
        hdr    = hdul[0].header
        hdu = hdul["VISIBILITY"]

        inst = cls(pos    = grb.radec,
                   site   = obs_loc["CTA"][loc],
                   window = [grb.tstart, grb.tstop],
                   name   = grb.id+"_"+loc)

        vis = Table.read(hdul,hdu=hdu)
        inst.status = "built-in"

        # ----------------------------------------------------
        def f(t,loc):
            t = Time(t.data,format="jd",scale="utc")
            if loc == "North":
                dt= [ [ t[0:2][0], t[0:2][1] ] ]
            if loc == "South":
                dt= [ [ t[2:4][0], t[2:4][1] ] ]

            # Check undefined intervals
            if dt[0][0].value == -9 or dt[0][1].value == -9:
                return [[]]
            else:
                return dt
        # ----------------------------------------------------

         # Visibility has been computed with this minimum altitude
        if loc == "North":
            inst.vis          = hdr['V_N_ANYT']
            inst.vis_night    = hdr['V_N_TNG'] # formerly vis_tonight
            inst.vis_prompt   = hdr['V_N_PR']

        if loc == "South":
            inst.vis          = hdr['V_S_ANYT']
            inst.vis_night    = hdr['V_S_TNG'] # formerly vis_tonight
            inst.vis_prompt   = hdr['V_S_PR']

        inst.altmin     = vis.meta["MIN_ALT"]*u.deg # Minimum altitude
        inst.t_true     = f(vis["True"],loc)
        inst.t_twilight = f(vis["Twilight"],loc)
        inst.t_event    = f(vis["Event"],loc)

        return inst

    ###-----------------------------------------------------------------------
    @classmethod
    def from_dict(cls, d):

        """
        Reads a Visibility instance from a dictionnary.

        Parameters
        ----------
        d : Dictionnary
            Dictionnary with the Visibility parameters.

        Returns
        -------
        inst : Visibility instance
            Instance.
       """

        inst = cls()
        inst.status = "From_dict"
        inst.site    = EarthLocation.from_geocentric(x=u.Quantity(d["site"][0]),
                                                    y=u.Quantity(d["site"][1]),
                                                    z=u.Quantity(d["site"][2]))
        inst.target = FixedTarget(SkyCoord(ra=u.Quantity(d["target"][0]),
                                           dec=u.Quantity(d["target"][1])),
                                           name="Source")

        inst.name            = d["name"]
        inst.moon_maxlight   = d["moon_maxlight"]
        inst.vis             = d["vis"]
        inst.vis_night       = d["vis_night"]
        inst.vis_prompt      = d["vis_prompt"]
        inst.moon_too_bright = d["moon_too_bright"]
        inst.moon_too_close  = d["moon_too_close"]

        for key in ["altmin","moon_maxalt", "moon_mindist"]:
            inst.__dict__[key] = u.Quantity(d[key])

        for key in ["tstart","tstop","t_true","t_twilight","t_event","t_moon_up"]:
            inst.__dict__[key] = Time(d[key],format="mjd")

        return inst

    ###-----------------------------------------------------------------------
    def to_json(self,file=None, debug=False):

        """
        Dump individual visibility instance to json.
        This is here essentially for didactical purpose as it is preferred to
        dump a dictionnary of visibilities.

        Parameters
        ----------
        file : pointer to file, optional
            Output file. The default is None.
        debug : boolean, optional
            If True, say something. The default is False.

        Returns
        -------
        None.

        """

        if debug:
            print("Write {} to Json output".format(self.name))
            record = json.dumps( self,
                                default= Visibility.object_to_serializable,
                                indent=None)
            print(record)

        json.dump(self, file,
                  default = Visibility.object_to_serializable, indent=None)

    ###------------------------------------------------------------------------
    @staticmethod
    def object_to_serializable(obj):

        """
        Turn any non-standard object into a serializable type so that
        JSON can write it.
        Note that this function explicitely returns a dictionnary.
        Written by K. Kosack, September 2022
        """
        from visibility import Visibility
        from astropy.time import Time
        import astropy.units.quantity
        from   astroplan import FixedTarget

        if isinstance(obj, Visibility):
            return {k:v for k,v in obj.__dict__.items()
                    if k not in Visibility.ignore}

        if isinstance(obj, Time): return obj.mjd
        if isinstance(obj, FixedTarget):
            return (str(u.Quantity(obj.ra)),
                    str(u.Quantity(obj.dec)) )
        if isinstance(obj, astropy.coordinates.earth.EarthLocation):
            return [str(obj.x), str(obj.y), str(obj.z)]
        if isinstance(obj, astropy.units.quantity.Quantity): return str(obj)

        return
    ###------------------------------------------------------------------------
    @staticmethod
    def params_from_key(keyword, parfile = None, debug=False):

        """
        Get the parameters to compute the visibility from the existing
        dictionnaries. Could be moved as a utility in the Configuration module.

        Parameters
        ----------
        keyword : string
            One of the key in the dictionnary file.
        parfile : string, optional
            The file name holding the dictionnaries. The default is
            `None`.
        debug : boolean, optional
            If True display information. The default is False.

        Returns
        -------
        Dictionnary
            Parameters used to compute the visibility.

        """
        if parfile is None:
            parfile = Visibility.def_vis_dicts

        try:
            with open(Path(Path(__file__).parent, parfile)) as file:
                visdict  = yaml.load(file, Loader=SafeLoader)
                if keyword in visdict.keys():

                    if debug:
                        print("   Vis. computed up to   : {} night(s)"
                                .format(keyword["depth"]))
                        print("   Skip up to            : {} night(s)"
                                .format(keyword["skip"]))
                        print("   Minimum altitude      : {}"
                                .format(keyword["altmin"]))
                        print("   Moon max. altitude    : {}"
                                .format(keyword["altmoon"]))
                        print("   Moon min. distance    : {}"
                                .format(keyword["moondist"]))
                        print("   Moon max. brightness  : {}"
                                .format(keyword["moonlight"]))

                    return visdict[keyword]
                else:
                    if debug:
                        print("{}.py: visibility keyword not referenced"
                             .format(__name__))
                    return None

        except IOError:
            sys.exit("{}.py: {} not found"
                     .format(__name__, parfile))

    ###------------------------------------------------------------------------
    def print(self, log=None):

        """
        Print out the GRB night, above-the-horizon periods, and default
        visibility window.

        Parameters
        ----------
        log : TextIO, optional
            Log file. The default is None.

        Returns
        -------
        None.

        """
        if log is None:
            log=Log()

        log.prt("=================   {:10s}   {:10s}   ================"
          .format(self.name,self.status))
        log.prt(' Visible  : {} - tonight, prompt : {}, {})'
              .format(self.vis, self.vis_night,self.vis_prompt))
        log.prt(' Altitude : Horizon > {:3.1f} - Moon veto > {:4.2f}'
              .format(self.altmin, self.moon_maxalt))
        #log.prt(" Trigger: {}".format(self.t_trig.datetime))

        #   WARNING : might be wrong since vis_tonight is not vis_night
        if self.vis_night or self.status =="Updated":
            log.prt("+----------------------------------------------------------------+")

            ###------------------
            def show(t,case="Dummy"):
                if len(t[0]) == 0:
                    log.prt(" {:6s} : {:26s} * {:26s}".format(case,"--","--"))
                    return
                for i, elt in enumerate(t):
                    t1 = elt[0]
                    t2 = elt[1]
                    log.prt(" {:6s} : {} * {}"
                            .format(case, Dp(t1), Dp(t2)),end="")
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
        In particular, allow checking the default visibility for altmin=10°
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
            if delta1 > delta_max or delta2 > delta_max:

                if log is not None:
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
        if self.vis != view.vis:
            log.prt(" ==> Vis wrong ", end="")
            matching = False
        if self.vis_tonight != view.vis_tonight:
            log.prt(" ==> Vis_tonight wrong ",end="")
            matching = False
        if self.vis_prompt != view.vis_prompt:
            log.prt(" ==> Vis_prompt wrong ",end="")
            matching = False
        if len(self.t_true) > 1:
            log.prt(" ==> >1 window",end="")
            matching = False

        if not matching:
            print()
            return matching

        # If visibilities agree and there is nothing to compare with, return True
        if self.vis == False:
            log.prt ("--- not visible")
            return True
        if self.tonight == False:
            log.prt ("--- not tonight")
            return True

        # both "tonight" are OK, compare delta-time

        if self.tonight and view.vis_tonight:
            print()
            # status(self.t_above,self.grb.t_event,case = "Above")
            # status(self.t_night,self.grb.t_twilight,case = "Night")
            matching = status(self.t_true[0],view.t_true[0],case = "Vis")
            if matching: log.prt("--- ok")
            else: log.prt("--- DOES NOT MATCH")
        return matching

    ###------------------------------------------------------------------------
    def plot(self,  grb,
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

        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from astropy.visualization import quantity_support

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

            duration = (self.tstop - self.tstart + dt_after + dt_before).to(u.d)
            dt   = np.linspace(0,duration.value, nalt)
            tobs = self.tstart  - dt_before + dt*duration.unit

            # Minimal altitude
            ax[0].axhline(y =self.altmin.value,
                       ls="--",color="tab:blue",label="Min. Alt.")

            ### GRB (main plot)
            grb.plot_altitude_and_flux(tobs, site=self.site, ax=ax[0])
            ax[0]. set_title(self.name)
            ax[0].grid("both",ls="--",alpha=0.5)

            ## GRB above horizon
            minalt, maxalt = ax[0].get_ylim()
            ymax = self.altmin.value/(maxalt - minalt)
            self.period_plot(self.t_event,
                        ax = ax[0],color="tab:blue",alpha=0.2,
                        tag="Above horizon",
                        ymin=0.,  ymax= ymax)

            ### Moon if requested
            if moon_alt or moon_dist:
                radec = get_moon(tobs, self.site) # Moon position along time

                if moon_alt:
                    ### Moon veto periods
                    self.period_plot(self.t_moon_up,
                                ax =ax[1],color="tab:orange",alpha=0.2,
                                tag="Moon veto")
                    ax[1].grid("both",ls="--",alpha=0.5)

                    ### Moon altitude
                    alt = radec.transform_to(AltAz(location=self.site,
                                                   obstime=tobs)).alt
                    moon_alt_plot(tobs, alt, ax = ax[1], alpha=0.5)
                    ax[1].axhline(y=self.moon_maxalt ,
                               color= "darkblue", alpha=1, ls=":")

                    if moon_dist:
                        ### Moon distance and brigthness
                        moon_dist_plot(grb.radec,tobs, radec,site=self.site,
                                        ax = ax[2],alpha=0.5) # Distance
                        self.period_plot(self.t_moon_up,
                        ax =ax[2],color="tab:orange",alpha=0.2, tag="Moon veto")
                        ax[2].grid("both",ls="--",alpha=0.5)
                        ax[2].axhline(y=self.moon_mindist ,
                                      color= "purple", ls=":",
                                      alpha=1, label="Min. distance")

                        axx = ax[2].twinx()
                        moonlight_plot(tobs,ax=axx) # Brightness
                        #moonphase_plot(tobs,ax=axx) # Phase (correlated to Brightness)

                        axx.axhline(y=self.moon_maxlight ,
                                    color= "tab:orange", ls=":",
                                    label="Max. illumination")

            ### Nights on all plots
            for axis in ax:
                if len(self.t_twilight) !=0:
                    self.period_plot(self.t_twilight,
                                ax =axis,color="black",alpha=0.1, tag="Night")

            ### Visibility windows - if the GRB is visible
            if self.vis:
                for axis in ax:
                    if len(self.t_true) !=0:
                        self.period_plot(self.t_true,
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
            axis.set_xlim([min(tobs).datetime,max(tobs).datetime])

        fig.tight_layout(h_pad=0, w_pad=0)

        return fig
    ###------------------------------------------------------------------------
    @staticmethod
    def period_plot(twindow,
                    ax=None,
                    alpha = 0.2, color="black", color2="black",tag="?",
                    tshift=0*u.s,
                    **kwargs):

        import matplotlib.pyplot as plt
        ax = plt.gca() if ax is None else ax

        if len(twindow[0]) == 0: return ax
        for elt in twindow:
            if elt[0] != -9 and elt[1] != -9:  # Check if -9, i.e. undefined
                t1  = (elt[0] - tshift).datetime
                t2  = (elt[1] - tshift).datetime
                if isinstance(tag, list):
                    ax.axvline(t1,color=color,label=tag[0])
                    ax.axvline(t2,color=color2,label=tag[1])
                else:
                    ax.axvspan(t1,t2,
                               alpha = alpha, color=color,label=tag,**kwargs)

        return ax

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

    # Complies with gamapy 1.2 installation
    import os
    os.environ["HAPPY_IN"] = "D:\\CTAO\SoHAPPy\input"
    os.environ["HAPPY_OUT"] = "D:\\CTAO\SoHAPPy\output"

    from configuration import Configuration
    from grb import GammaRayBurst

    # Read default configuration
    cf = Configuration()
    cf.find_and_read("config.yaml")

    # Supersede some parameters
    cf.prompt_dir = None
    cf.nsrc     = 1 # 250
    cf.ifirst   = [2] # ["190829A"]

    # Get visibility information
    visinfo = cf.decode_visibility_keyword()

    # Log fle
    log_filename    = Path(cf.out_dir,"/visibility.log")
    log = Log(log_name  = log_filename,talk=True)

    # Print configuration with possible superseded values
    cf.print(log)

    # Loop over GRB list
    data_path   = Path( Path(os.environ["HAPPY_IN"]),cf.data_dir) # Input data
    grblist     = cf.source_ids() # GRB list to be analysed

    for item in grblist:
        fname = cf.prefix+str(item)+cf.suffix

        grb = GammaRayBurst.from_fits(Path(data_path,fname),
                                      prompt  = cf.prompt_dir,
                                      ebl     = "dominguez")

        print(grb)
        grb.plot_energy_spectra()
        grb.plot_time_spectra()

        for loc in ["North","South"]:
            grb.set_visibility(item,loc,info=visinfo)
            grb.vis[loc].print(log)
            grb.vis[loc].plot(grb)

    log.close()

    print(" All done !")
