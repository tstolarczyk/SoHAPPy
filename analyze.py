# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 12:50:28 2022.

@author: Stolar
"""
import sys
import pickle
import numpy as np
import astropy
import astropy.units as u
from astropy.coordinates import AltAz
from astropy.visualization import quantity_support

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from gammapy.stats import WStatCountsStatistic

from niceprint import t_fmt

# Defaults colors for the plots
col3 = "tab:orange"  # 3 sigma
col5 = "tab:green"   # 5 sigma
colmx = "black"      # Sigma max

__all__ = ["Analysis"]


###############################################################################
class Analysis():
    """Handle Monte Carlo, detection results and analysis outputs."""

    ignore = ["slot",
              "det_level",
              "alpha",
              "nstat",
              "id_smax_list",
              "smax_list",
              "nex_mx_list",
              "nb_mx_list",
              "id_3s_list",
              "nex_3s_list",
              "nb_3s_list",
              "id_5s_list",
              "nex_5s_list",
              "nb_5s_list",
              "sigma_mean",
              "sigma_std",
              "abort"]
    """
    Members of the class not to be dumped into the population file
    """

    # ##-----------------------------------------------------------------------
    def __init__(self, slot, nstat=None, alpha=0., cl=0., loc=None):
        """
        Define the default instance.

        Note that the site can be obtained from the slot instance. But
        it is useful to attach a site to an analysis even if the slot
        (and therefore the slot site) is not defined.
        A subset of all member values are dumped to the output file as a table
        with columns having the member name (better to keep them short for the
        analysis). Some members are not written out (internal usage) and are
        defined in the 'ignore' list. Some members from external classes have
        to be dumped out : as the analysis instance can be created from a
        dummy empty slot, some default values are defined here.

        Parameters
        ----------
        slot : A Slot instance
            The slot, a set of time slices, to be analyzed.
        nstat : integer, optional
            The number of trials in the simulation. The default is None.
        alpha : float, optional
            The on to off ration number. The default is 0.
        cl : float, optional
            The fraction of Monte Carlo trials that have lead to a detection
            at a certain (3 sigma) siginificance. The default is 0.
        loc : string, optional
            The site considered for the analysis. The default is None.

        Returns
        -------
        None.

        """
        self.slot = slot  # Will be updated after the slot is dressed

        # Allows superseding the name in particular for an empty slot
        if loc is None:
            loc = slot.loc

        self.name = slot.grb.id  # Source name

        # Done here to keep order in the ouput file
        self.loca = loc  # Site nickname - loc is a pandas function !

        # slot related information to be dumped to the output - dummy so far
        self.nt = len(slot.slices)
        self.t1 = -1.
        self.t2 = -1.
        self.alt1 = self.alt2 = self.az1 = self.az2 = -1*u.deg
        self.prpt = -1  # O or 1 if prompt if visible (not valid for "Both")

        self.det_level = cl  # Fraction declaring a detection above 3,5 sig.
        self.alpha = alpha   # 5 zones (1/0.2) for estimating B in off-regions

        if nstat is not None:
            self.nstat = nstat
        else:
            sys.exit("{}.py: nstat should be defined ", format(__name__))

        # ##-------------------------
        # Significances over simulations
        # ##-------------------------
        self.id_smax_list = []  # Slice indices to get back the time/ altaz
        self.smax_list = []  # List of max. significances along trials
        self.nex_mx_list = []  # List of on counts at max. signif.
        self.nb_mx_list = []  # List of off counts at max. signif.

        self.id_3s_list = []  # Slice indices ot get back the time/altaz
        self.nex_3s_list = []  # List of on counts
        self.nb_3s_list = []  # List of off counts
        self.d3s = 0   # Number of trials 3 sigma was reached

        self.id_5s_list = []  # Slice indices to get back the time/altaz
        self.nex_5s_list = []  # List of on counts
        self.nb_5s_list = []  # List of off counts
        self.d5s = 0   # Number of trials 5 sigma was reached

        # Mean sigma versus time - one value per time slice
        self.sigma_mean = []  # The length is unknow as slot.dress can modify
        self.sigma_std = []

        # ##-------------------------
        # Results
        # ##-------------------------

        # ## Maximal significance
        self.sigmx = -1.  # Mean maximum significance reached
        self.esigmx = -1.  # Error on mean maximum significance reached
        self.nexmx = -1.  # Mean number of excess counts at sigmx
        self.enexmx = -1.  # Error on mean number of excess counts at sigmx
        self.nbmx = -1.  # Mean number of background counts at sigmx
        self.enbmx = -1.  # Error on mean number of background counts at sigmx

        self.tmx = -1*u.s  # Mean time at sigmx (related to time binning)
        self.etmx = -1*u.s  # Error on mean time at sigmx
        self.altmx = -1*u.deg  # Mean altitude at sigmx
        self.ealtmx = -1*u.deg  # Error on mean altitude at sigmx
        self.azmx = -1*u.deg  # Mean azimuth at sigmx
        self.eazmx = -1*u.deg  # Error on mean azimuth at sigmx (deg)

        # ## 3 sigma
        self.nex3s = -1.  # Mean number of excess counts at 3 sigma
        self.enex3s = -1.  # Error on mean number of excess counts at 3 sigma
        self.nb3s = -1.  # Mean number of background counts at 3 sigma
        self.enb3s = -1.  # Error on mean number of bckg. counts at 3 sigma

        self.t3s = -1*u.s  # Mean time at 3 sigma (related tot time binning)
        self.et3s = -1*u.s  # Error on mean time at 3 sigma
        self.alt3s = -1*u.deg  # Mean altitude at 3 sigma
        self.ealt3s = -1*u.deg  # Error on mean altitude at 3 sigma
        self.az3s = -1*u.deg  # Mean azimuth at 3 sigma
        self.eaz3s = -1*u.deg  # Error on mean azimuth at 3 sigma

        self.d3s = 0  # Number of trials reaching 3 sigma

        # ## 5 sigma
        self.nex5s = -1.  # Mean number of excess counts at 5 sigma
        self.enex5s = -1.  # Error on mean number of excess counts at 5 sigma
        self.nb5s = -1.  # Mean number of background counts at 5 sigma
        self.enb5s = -1.  # Error on mean number of bckg. counts at 5 sigma

        self.t5s = -1*u.s  # Mean time at 5 sigma (related tot time binning)
        self.et5s = -1*u.s  # Error on mean time at 5 sigma
        self.alt5s = -1*u.deg  # Mean altitude at 5 sigma
        self.ealt5s = -1*u.deg  # Error on mean altitude at 5 sigma
        self.az5s = -1*u.deg  # Mean azimuth at 5 sigma
        self.eaz5s = -1*u.deg  # Error on mean azimuth at 5 sigma

        self.d5s = 0  # Number of trials reaching 5 sigma

        # ## Error variables
        # * If positive, give the time slice at which the simulation stopped.
        # A complete simulation has err = number of trials requested.
        # * If negative:
        # - n : simulations completed but errors reported in n trials
        # - 999 : simulation/analysis not done (GRB not visible)
        self.err = -999
        self.abort = False  # True if analysis was aborted

    # ##-----------------------------------------------------------------------
    def run(self):
        """
        Analyse data cumlulated by the :func:`analyze.Analysis.fill` function.

        Modify the current instance with new values.
        Compute and store :
            * Significance values for all slices in the trials;
            * Maximum significance reached for the GRB;
            * Time of maximum significance, number of excess and bckg. events;
            * Time to reach 3 or 5 sigma, number of excess and bckg. events.

        Returns
        -------
        None.

        """
        # Add information on the slot  - now that we know that it exists
        self.nt = len(self.slot.slices)
        self.t1 = self.slot.tstart
        self.t2 = self.slot.tstop

        if self.loca in ["North", "South"]:
            altaz = [self.slot.grb.altaz(dt=self.t1, loc=self.loca),
                     self.slot.grb.altaz(dt=self.t2, loc=self.loca)]
            self.alt1 = altaz[0].alt
            self.alt2 = altaz[1].alt
            self.az1 = altaz[0].az
            self.az2 = altaz[1].az

        if self.loca in ["North", "South"]:
            self.prpt = int(self.slot.grb.vis[self.loca].vis_prompt)

        self.err = self.nstat  # GRB visible, simul. ought to be complete

        # From the cumlated values compute the mean and standard values
        # for each time slices
        self.sigma_mean = self.sigma_mean/self.nstat
        sigma2_mean = self.sigma_std/self.nstat
        self.sigma_std = np.sqrt(sigma2_mean-self.sigma_mean**2)

        # Max siginificance - mean values and standard deviations
        self.sigmx = np.mean(self.smax_list)
        self.esigmx = np.std(self.smax_list)
        self.nexmx = np.mean(self.nex_mx_list)
        self.enexmx = np.std(self.nex_mx_list)
        self.nbmx = np.mean(self.nb_mx_list)
        self.enbmx = np.std(self.nb_mx_list)

        self.tmx, self.etmx, \
            self.altmx, self.ealtmx, \
            self.azmx, self.eazmx \
            = self.mean_values(self.id_smax_list)

        # ## 3 sigma
        if len(self.id_3s_list) != 0:
            nex = [n for n in self.nex_3s_list if n >= 0]
            nb = [n for n in self.nb_3s_list if n >= 0]
            self.nex3s = np.mean(nex)
            self.enex3s = np.std(nex)
            self.nb3s = np.mean(nb)
            self.enb3s = np.std(nb)

        self.t3s, self.et3s, \
            self.alt3s, self.ealt3s, \
            self.az3s, self.eaz3s = self.mean_values(self.id_3s_list)

        # ## 5 sigma
        if len(self.id_5s_list) != 0:
            nex = [n for n in self.nex_5s_list if n >= 0]
            nb = [n for n in self.nb_5s_list if n >= 0]
            self.nex5s = np.mean(nex)
            self.enex5s = np.std(nex)
            self.nb5s = np.mean(nb)
            self.enb5s = np.std(nb)

        self.t5s,  self.et5s, \
            self.alt5s, self.ealt5s, \
            self.az5s, self.eaz5s = self.mean_values(self.id_5s_list)

    # ##-----------------------------------------------------------------------
    def mean_values(self, idlist):
        """
        Extract mean values and RMS from the index list.

        Parameters
        ----------
        idlist : list of indices
            List of slice indices to be analysed.

        Returns
        -------
        t_det : Astropy Quantity
            Detection time.
        e_t_det : Astropy Quantity
            Error on detection time.
        alt_det : Astropy Angle
            Altitude of detection.
        e_alt_det : Astropy Angle
            Error on the altiude of detection.
        az_det : Astropy Angle
            Azimut angle of the detection.
        e_az_det : Astropy Angle
            Error on the azimut angle of detection.

        """
        t_det, e_t_det = 2*(-1*u.s,)
        alt_det, e_alt_det, az_det, e_az_det = 4*(-1*u.deg,)

        if len(idlist) != 0:
            t_unit = self.slot.slices[0].tobs().unit

            t = [self.slot.slices[i].tobs().value for i in idlist
                 if self.slot.slices[i].tobs().value >= 0]
            t_det = np.mean(t) * t_unit
            e_t_det = np.std(t) * t_unit

            if self.slot.loc != "Both":
                altaz = self.slot.grb.altaz(dt=t*t_unit,
                                            loc=self.slot.loc)
                alt = altaz.alt
                alt_det = np.mean(alt)
                e_alt_det = np.std(alt)
                az = altaz.az
                az_det = np.mean(az)
                e_az_det = np.std(az)

        return t_det, e_t_det, alt_det, e_alt_det, az_det, e_az_det

    # ##-----------------------------------------------------------------------
    def fill(self, itrial, non, noff, dbg=0):
        """
        Fill arrays with information from one Monte Carlo iteration.

        Keep track of slice indices where a 3 sigma or 5 sigma detection
        occured, and where the max siginificance was reached and at which
        value.

        Also cumulate the sigma and sigma**2 values in each time slices in
        order to later compute the mean and standard deviation of the
        significance for each time slice.

        The siginificance is computed under the assmption of a measured
        background, i.e. with possble fluctuation.
        Note that WStatCountsStatistic returns a significance with the
        sign of the excess (negative excess give negative significance).
        If the count number is not fluctuated, the excess can only be
        positive since it is obtained from a physical flux. Excess at zero
        gives a significance at zeo.

        More on the various statistics `onthis page
        <https://docs.gammapy.org/0.17/stats/fit_statistics.html>`_ :

        Parameters
        ----------
        itrial: integer
            The Mone Carlo simulation current number of trials. Used to abort
            the simulation if 3s detection is not reachable.
        non : list of float
            Cumulated count numbers in the signal plus background region along
            time slices.
        noff : list of float
            Cumulated count numbers in the background regions along time
            slices.
        dbg : integer, optional
            Debugging flag. The default is 0.

        Returns
        -------
        None.

        """
        self.err = itrial

        # Compute the significances of all slices for the current trail
        wstat = WStatCountsStatistic(n_on=non,
                                     n_off=noff,
                                     alpha=self.alpha)
        sigma = wstat.sqrt_ts  # ? check

        # Cumulate sigma and sigma squared values for each slice
        if len(self.sigma_mean) == 0:
            self.sigma_mean = np.zeros(len(self.slot.slices))
            self.sigma_std = np.zeros(len(self.slot.slices))
        self.sigma_mean += sigma    # To be averaged later
        self.sigma_std += sigma**2  # To be averaged later

        # Go back to excess and background counts
        nb = noff * self.alpha
        nex = non - nb

        # Find slice with maximum - cannot be a slice with non or noff below
        # limit
        sigmx = np.nanmax(sigma)  # Returns Nan only if all are Nan
        if np.isnan(sigmx):
            print(" All sigma values are Nan !!! ")
            maxidx = -1
            sigmx = -999
            nex_mx = -1
            nb_mx = -1
        else:
            maxidx = np.nanargmax(sigma)  # If nan, it would be the max !
            nex_mx = nex[maxidx]
            nb_mx = nb[maxidx]

        self.id_smax_list.append(maxidx)
        self.smax_list.append(sigmx)
        self.nex_mx_list.append(nex_mx)
        self.nb_mx_list.append(nb_mx)

        # Find where 3 sigma is reached
        nex_3s = -1
        nb_3s = -1
        id_3s = np.where(sigma >= 3)[0]  # This ignore Nan values

        if np.size(id_3s) != 0:
            id_3s = id_3s[0]  # First one
            nex_3s = nex[id_3s]
            nb_3s = nb[id_3s]
            self.id_3s_list.append(id_3s)
            self.d3s += 1

        self.nex_3s_list.append(nex_3s)
        self.nb_3s_list.append(nb_3s)

        # Find where 5 sigma is reached
        nex_5s = -1
        nb_5s = -1
        id_5s = np.where(sigma >= 5)[0]  # This ignore Nan values
        if np.size(id_5s) != 0:
            id_5s = id_5s[0]  # First one
            nex_5s = nex[id_5s]
            nb_5s = nb[id_5s]
            self.id_5s_list.append(id_5s)
            self.d5s += 1

        self.nex_5s_list.append(nex_5s)
        self.nb_5s_list.append(nb_5s)

        if dbg > 2:
            print(" >>> t_smax = {:8.2f}, ({:5.2f} sigma at slice {:2d})"
                  .format(self.slot.slices[maxidx].tobs(),
                          sigma[maxidx],
                          maxidx), end="")
            if np.size(id_3s) != 0:
                print(" - 3s at t3s = {:8.2f} at slice {:2d})"
                      .format(self.slot.slices[id_3s].tobs(), id_3s))
            else:
                print(" - 3s not reached")

        # If the fraction of trial exceeds 1-CL, check is 1 detection is found
        if itrial/self.nstat > 1 - self.det_level:
            if self.d3s == 0:
                self.abort = True

        return sigma

    # ##-----------------------------------------------------------------------
    def summary(self, log=None):
        """Print out a summary."""
        log.prt(f" GRB {self.slot.grb.name[5:]:<4} {self.slot.loc:<s}")
        log.prt(f"  z    = {self.slot.grb.z:6.2f}")
        log.prt(f"  Eiso = {self.slot.grb.Eiso:6.2e}")
        log.prt("  Detection")
        log.prt(f"   - Visible : {t_fmt(self.slot.tstart):6.2f} -"
                " {t_fmt(self.slot.tstop)):4.2f}")
        log.prt(f"   - Slices  : {len(self.slot.slices):d}")
        log.prt(f"   - Delay   : {t_fmt(self.slot.delay):6.2f} ")
        log.prt(f"   - sigmax  : {self.sigmx:<6.2f} @ {self.tmx:<6.2f}")
        log.prt(f"   - 5 sigma @ {self.t5s:<6.2f}")

    # ##-----------------------------------------------------------------------
    def print(self, log=None):
        """Dump results."""
        # ---------------------------------------
        def time_convert(time, etime):
            """
            Simplify the time printing.

            Parameters
            ----------
            time : TYPE
                DESCRIPTION.
            etime : TYPE
                DESCRIPTION.

            Returns
            -------
            time : TYPE
                DESCRIPTION.
            etime : TYPE
                DESCRIPTION.
            utime : TYPE
                DESCRIPTION.

            """
            time = t_fmt(time)
            utime = time.unit
            etime = etime.to(time.unit)

            return time, etime, utime
        # ---------------------------------------

        if self.d5s/self.nstat >= self.det_level:
            message = "5 sigma detected at " \
                      + str(round(100*self.d5s/self.nstat, 1))+"%"
        elif self.d3s/self.nstat >= self.det_level:
            message = "3 sigma detected at " \
                      + str(round(100*self.d3s/self.nstat, 1))+"%"
        else:
            message = "NOT detected    "

        log.prt("| Analysis     |  ===> ", end="")
        log.highlight(f"{message:42s}", end="")
        log.prt("|")
        # .....................................................................
        log.prt("+" + 14*"-" + "+" + 49*"-" + "+")
        log.prt(" Window : {:6.3f} - {:6.3f} * Delay: {:6.2f} * {:3d} slices"
                .format(t_fmt(self.slot.tstart), t_fmt(self.slot.tstop),
                        t_fmt(self.slot.delay), len(self.slot.slices)))
        log.prt("+" + 64*"-" + "+")
        # .....................................................................
        message = "  RESULTS      : {}{:<10s}"\
            .format(14*" ", self.name+"_"+self.loca)
        log.highlight(f"{message:66s}")
        log.prt("+" + 64*"-" + "+")

        # ### Max significance and related values
        log.prt(f"            Max. sigma = {self.sigmx:>8.1f}  "
                f"+/- {self.esigmx:<8.1f}")

        tmx, etmx, utmx = time_convert(self.tmx, self.etmx)
        log.prt(f"                  time = {tmx.value:>8.2f}  "
                f"+/- {etmx.value:<8.2f} {utmx:5s}")

        log.prt(f"                  alt. = {self.altmx.value:>8.2f}  "
                f"+/- {self.ealtmx.value:<8.2f} {self.altmx.unit:5s}")
        log.prt(f"                   az. = {self.azmx.value:>8.2f}  "
                f"+/- {self.eazmx.value:<8.2f} {self.azmx.unit:5s}")
        log.prt(f"           S, B counts = {self.nexmx:>8.1f} "
                f"{self.nbmx:0.1f}")
        log.prt(f"        on, off counts = {self.nexmx + self.nbmx:>8.1f} "
                f"{self.nbmx/self.alpha:0.1f}")
        log.prt("+" + 64*"-" + "+")
        log.prt("                           3 sig                   5 sig")
        log.prt("+" + 64*"-" + "+")

        # ## 3 and 5 sigma times and other values
        t3s, et3, ut3 = time_convert(self.t3s, self.et3s)
        t5s, et5, ut5 = time_convert(self.t5s, self.et5s)
        log.prt(f" Time ({ut3:3s})     = "
                f"{t3s.value:>8.2f} +/- {et3.value:<8.2f} "
                f" {t5s.to(ut3).value:>8.2f} "
                f" +/- {et5.to(ut3).value:<8.2f}"
                )

        log.prt(f" Alt. ({self.alt3s.unit:3s})     = "
                f"{self.alt3s.value:>8.2f}  +/- {self.ealt3s.value:<8.2f}  "
                f"{self.alt5s.value:>8.2f}  +/- {self.ealt5s.value:<8.2f}"
                )
        log.prt(f" Az.  ({self.az3s.unit:3s})     = "
                f"{self.az3s.value:>8.2f}  +/- {self.eaz3s.value:<8.2f}  "
                f"{self.az5s.value:>8.2f}  +/- {self.eaz5s.value:<8.2f}"
                )
        log.prt(f" Det. level (%) =  {100*self.d3s/self.nstat:>8.2f}"
                f"                {100*self.d5s/self.nstat:>8.2f}"
                )
        log.prt(f" S, B           = {self.nex3s:>8.1f}, {self.nb3s:<8.1f}"
                f"      {self.nex5s:>8.1f}, {self.nb5s:<8.1f}"
                )
        log.prt(f" On, Off        = {self.nex3s + self.nb3s:>8.1f}, "
                f"{self.nb3s/self.alpha:<8.1f}      "
                f"{self.nex5s + self.nb5s:>8.1f}, "
                f"{self.nb5s/self.alpha:<8.1f}"
                )
        log.prt(66*"+")
        log.prt(f"{'SoHAPPy!':>}")

    # ##------------------------------------------------------------------------
    def dump_to_file(self, grb, pop, header=False, debug=False):
        """
        Dump the class content into a file with an appropriate format.

        Parameters
        ----------
        pop : pointer to file
            The output file.
        grb : GammaRayBurst instance
            The current GRB.
        header : boolean, optional
            If True, write the header to file before the grb data.
            The default is False. It shall be True for the first GRB to have a
            header on the first line of the file.

        Returns
        -------
        Boolean
            Return False to ensure that the header is not printed out twice.

        """
        # ## ------------------------------------------------------------
        def head_fmt(kh, vh):
            if kh == "name":
                return ">25s"
            if kh == "t_trig":
                return ">23s"
            if kh in ("radec", '   ra   dec'):
                return ">14s"
            if isinstance(vh, str):
                return ">10s"
            if isinstance(vh, int):
                return ">4s"
            if isinstance(vh, (float, np.float32, np.float64)):
                return ">10s"
            if isinstance(vh, (astropy.time.core.Time)):
                return ">23s"
            if isinstance(vh, (astropy.units.quantity.Quantity)):
                return ">23s"

            return ""

        def val_fmt(k, v):
            if k == "name":
                return ">25s"
            if k == "t_trig":  # "str" returns to few spaces
                return ">23s"
            if k in ("radec", '   ra   dec'):
                return ">13s"
            if isinstance(v, str):
                return ">10s"
            if isinstance(v, int):
                return ">4d"
            if isinstance(v, astropy.time.core.Time):
                return ">23s"
                # return ">12.4f"
            if isinstance(v, (float, np.float32, np.float64)):
                if abs(v) < 3.6e2 and abs(v) > 1e-2:
                    return ">10.4f"
                else:
                    return ">10.4e"
            return ""

        def values(var):
            if isinstance(var, astropy.coordinates.sky_coordinate.SkyCoord):
                var = [y for y in [var.ra.value, var.dec.value]]
                return f"{var[0]:>6.2f} {var[1]:>6.2f}"
            if isinstance(var, astropy.units.quantity.Quantity):
                return var.value
            if isinstance(var, astropy.time.core.Time):
                return var.isot
            return var

        # ##------------------------------------------------------------

        # Table header
        if header:

            for data in (self.__dict__.items(), grb.__dict__.items()):
                for k, v in data:
                    # Pay attention to identical member names
                    if k in self.ignore or k in grb.ignore:
                        continue
                    if k == "radec":
                        k = f"{'ra':>6s} {'dec':>6s}"
                    if debug:
                        print(f"> {k}: {v} -> {values(v)} "
                              f"head_fmt(v) = {head_fmt(k, values(v))}")
                    print("{val:{fmt}}".format(val=k,
                                               fmt=head_fmt(k, values(v))),
                          end=" ", file=pop)
            print("", file=pop)

        # Table data
        for data in (self.__dict__.items(), grb.__dict__.items()):
            for key, val in data:
                # Pay attention to identical member names
                if key in self.ignore or key in grb.ignore:
                    continue
                if debug:
                    print(f"> {key}: {val}, "
                          f"values(val) = {val_fmt(key, values(val))}")
                print(f"{values(val):{val_fmt(key, values(val))}}",
                      end=" ", file=pop)
        print("", file=pop)

        return False

    # ##-----------------------------------------------------------------------
    def write(self, filename=None):
        """
        Save the present class to disk for further use.

        Parameters
        ----------
        filename : string, optional
            Output file name. The default is "None".

        Returns
        -------
        None.

        """
        if filename is None:
            sys.exit("Output file not defined)")

        outfile = open(filename, "wb")
        pickle.dump(self, outfile)
        outfile.close()

        print(f" Saving analysis to file : {filename}")

    # ##-----------------------------------------------------------------------
    def show(self, pdf=None):
        """Display analysis relevant plots."""
        fig_sig = self.plot_sigma_vs_time()

        if self.nstat > 1:
            fig_n = None  # self.plot_non_vs_noff()
        else:
            fig_n = None  # Otherwise plot one single point

        fig_vis = self.plot_story(ref="VIS")
        fig_grb = self.plot_story(ref="GRB")

        if pdf is not None:
            pdf.savefig(fig_sig)
            if fig_n is not None:
                pdf.savefig(fig_n)
            pdf.savefig(fig_vis)
            pdf.savefig(fig_grb)

        # used to pause plot display in interactive mode on a shell script.
        # In the abscence of a call to that function figures will stack on the
        # screen during the run and all disappear at the end of the run.
        # Using this, figures will be stacked on screen and displayed for each
        # event until they are closed by the user.
        plt.show(block=True)

        return

    # ##-----------------------------------------------------------------------
    def plot_sigma_vs_time(self):
        """
        Display mean siginificance versus time.

        Returns
        -------
        fig : matplolib figure
            Current figure.

        """
        # If only one MC trial, no max significance distribution
        if self.nstat > 1:
            fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 5),
                                           gridspec_kw={'width_ratios':
                                                        [2, 1]})
        else:
            fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))

        # Measurement points and units
        base_unit = self.slot.slices[0].tobs().unit
        t_s = np.asarray([s.tobs().value for s in self.slot.slices])*base_unit
        t_unit = t_fmt(self.tmx).unit

        t_s = t_s.to(t_unit)
        tmax = self.tmx.to(t_unit)
        t3s = self.t3s.to(t_unit)
        t5s = self.t5s.to(t_unit)

        with quantity_support():

            # ## Mean significance at each slice
            ax1.errorbar(t_s, self.sigma_mean, yerr=self.sigma_std, fmt='o')

            for t, sig, errs, c, tag in zip(
                                [tmax, t3s, t5s],
                                [self.sigmx, 3, 5],
                                [self.esigmx, 0, 0],
                                [colmx, col3, col5],
                                [r"$\sigma_{max}$", r"$3\sigma$", r"$5\sigma$"]
                                              ):
                import matplotlib
                if matplotlib.__version__ < "3.4.2":
                    # Does not support Quantity errors
                    ax1.errorbar(x=np.mean(t).value, y=sig,
                                 xerr=np.std(t).value, yerr=errs,
                                 fmt="o", color=c, label=tag)
                    ax1.set_xlabel('Obs. duration (' + str(t_unit) + ")")

                else:
                    ax1.errorbar(x=np.mean(t), y=sig,
                                 xerr=np.std(t), yerr=errs,
                                 fmt="o", color=c, label=tag)
                    ax1.set_xlabel('Obs. duration ('+ax1.get_xlabel()+")")

                xmin, xmax = ax1.get_xlim()
                ymin, ymax = ax1.get_ylim()
                ax1.vlines(np.mean(t), ymin=ymin, ymax=self.sigmx,
                           alpha=0.5, ls=":", color=c)
                ax1.hlines(sig, xmin=xmin, ls=":", xmax=np.mean(t),
                           alpha=0.5, color=c)

            ax1.set_ylabel(r"Significance $\sigma$")
            ax1.legend()

            import matplotlib
            if matplotlib.__version__ >= "3.5.0":
                ax1.set_xscale("log", nonpositive='clip')  # 3.5.0
            else:
                ax1.set_xscale("log", nonposx='clip')  # 3.1.1

            ax1.grid(which='both', alpha=0.2)
            ttl = self.name+"_" + self.loca
            ttl += ' (' + str(self.nstat) + ' iter.)'
            ax1.set_title(ttl, loc="right")

            # Sigma max distribution (if niter > 1)
            if self.nstat > 1:
                ax2.hist(self.smax_list,
                         bins=max(int(self.nstat/2), 1),  # Cannot be < 1
                         range=[self.sigmx - 3*self.esigmx,
                                self.sigmx + 3*self.esigmx],
                         alpha=0.5,
                         color="grey",
                         label=r" {:5.1} $\pm$ {:5.1}"
                         .format(self.sigmx, self.esigmx))
                ax2.set_xlabel(r"$\sigma_{max}$")
                ax2.set_ylabel('Trials')

        plt.tight_layout()

        return fig

    # ##-----------------------------------------------------------------------
    def plot_non_vs_noff(self, ax=None, logscale=True):
        """
        Plot Non versus Noff from the background and excess counts.

        Draw the error bars from the variances.

        Parameters
        ----------
        ax : Matplolib axes, optional
            Current axis. The default is None.
        logscale : String, optional
            If True, logarithmic scales. The default is True.

        Returns
        -------
        fig : Matplotlib figure
            Current figure.

        """
        if ax is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))
        else:
            fig = ax.get_figure()

        nexlist = [self.nex_3s_list, self.nex_5s_list, self.nex_mx_list]
        nblist = [self.nb_3s_list, self.nb_5s_list, self.nb_mx_list]
        collist = [col3, col5, colmx]
        taglist = [r"$3\sigma$", r"$5\sigma$", r"$\sigma_{max}$"]

        for nex, nb, col, tag in zip(nexlist, nblist, collist, taglist):

            non = np.add(nex, nb)
            noff = np.true_divide(nb, self.alpha)
            non_mean = np.mean(non)
            non_std = np.std(non)
            noff_mean = np.mean(noff)
            noff_std = np.std(noff)
            xlabel = "$N_{off}$"
            ylabel = "$N_{on}$"
            if logscale:
                non = np.log10(non)
                non_mean = np.log10(non_mean)
                non_std = np.log10(non_std)
                noff = np.log10(noff)
                noff_mean = np.log10(noff_mean)
                noff_std = np.log10(noff_std)
                xlabel = "$log_{10}$" + xlabel
                ylabel = "$log_{10}$" + ylabel

            ax.plot(noff, non,
                    alpha=0.3, marker=".", ls="", color=col, label=tag)
            ax.errorbar(x=noff_mean, y=non_mean,
                        xerr=np.std(noff), yerr=np.std(non), color=col)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            title = self.name+"_" + self.loca \
                + ' ('+str(self.nstat) + ' iter.)'
            ax.set_title(title, loc="right")
            ax.grid(which='both', alpha=0.2)

        plt.tight_layout()

        return fig

    # ##----------------------------------------------------------------------------
    def plot_story(self, ref="VIS"):
        """Show altitude versus time and points used for computation."""
        if (self.loca != "North" and self.loca != "South"):
            return None

        # Plot sigma versus time on the full
        # GRB lifetime together with visibility

        coord = self.slot.grb.vis[self.loca].site
        dtplot = 500*u.s  # Add some space for the plot

        # These are the visibilities in the sense of the Dark time
        # Create points from the start of the first night to the end of the
        # last night
        tvis1 = self.slot.grb.vis[self.loca].t_twilight[0][0]
        tvis2 = self.slot.grb.vis[self.loca].t_twilight[-1][1]
        dtvis = (tvis2-tvis1).to(u.s)
        tvis = np.linspace(0, dtvis.value, 100)*dtvis.unit

        # Set the reference time
        event_name = self.slot.grb.id+" - " + self.loca
        if ref == "VIS":
            tref = tvis1
            # text_ref =event_name + " (Ref: vis. start)"
            text_ref = event_name + " (Ref: Dark time start)"
        else:
            tref = self.slot.grb.t_trig
            text_ref = event_name + " (Ref: GRB trigger)"

        t_trig_rel = self.slot.grb.t_trig - tref

        # Recompute time sampling in the two reference systems
        tvis_abs = tvis1 + tvis
        tvis_rel = tvis_abs - tref

        # Compute alt-az points for visibility time sampling
        altazvis = self.slot.grb.radec.transform_to(AltAz(obstime=tvis_abs,
                                                          location=coord))

        # Compute GRB measurement points in the two refrence systems
        tgrb_abs = self.slot.grb.t_trig + self.slot.grb.tval  # Absolute time
        tgrb_rel = tgrb_abs - tref
        altazgrb = self.slot.grb.radec.transform_to(AltAz(obstime=tgrb_abs,
                                                          location=coord))

        # Compute GRB observation points
        t_s = [s.tobs() for s in self.slot.slices] \
            * self.slot.slices[0].tobs().unit
        tgrb_obs_abs = self.slot.grb.t_trig + t_s  # Absolute time
        tgrb_obs_rel = tgrb_obs_abs - tref
        altazobs = self.slot.grb.radec.transform_to(AltAz(obstime=tgrb_obs_abs,
                                                          location=coord))

        with quantity_support():
            fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

            # Plot limits
            xmin = -dtplot
            xmax = (tvis2 - tref + dtplot).sec
            ymin = 0*u.deg

            # Relative to reference
            ax1.plot(tvis_rel.sec, altazvis.alt,
                     linestyle="--", lw=1, color="tab:blue", label="Altitude")

            ax1.plot(tgrb_rel.sec, altazgrb.alt,
                     linestyle="", marker="*", markerfacecolor="b",
                     label="End of intervals")

            ax1.plot(tgrb_obs_rel.sec, altazobs.alt,
                     linestyle="",
                     marker="o",
                     markersize=10,
                     markerfacecolor="none",
                     markeredgewidth=1.5,
                     markeredgecolor="r",
                     label="Measurements")

            # Trigger
            ax1.axvline(x=t_trig_rel.sec,
                        ls=(0, (1, 1)),  # densely dotted
                        color="grey",
                        label="Trigger")

            # Dark time
            first = True

            for elt in self.slot.grb.vis[self.loca].t_twilight:
                if first:
                    label = "Night"
                    first = False
                else:
                    label = None
                ax1.axvspan((elt[0] - tref).sec, (elt[1] - tref).sec,
                            alpha=0.2, color="black", label=label)

            # minimal altitude in simulation
            ax1.axhline(y=self.slot.grb.vis[self.loca].altmin, ls="dashdot",
                        color="grey")

            # Sig max
            ax1.errorbar(self.tmx.to(u.s).value + t_trig_rel.sec,
                         self.altmx.value,
                         xerr=self.etmx.to(u.s).value,
                         yerr=self.ealtmx.value,
                         fmt='o', color=colmx, label=r"$\sigma_{max}$")
            ax1.vlines(self.tmx.to(u.s).value + t_trig_rel.sec,
                       ymin=0,
                       ymax=self.altmx.value,
                       alpha=0.5, lw=1, color=colmx)

            ax1.hlines(self.altmx.value,
                       xmin=xmin,
                       xmax=self.tmx.to(u.s).value + t_trig_rel.sec,
                       alpha=0.5, lw=1, color=colmx)

            # 3 sigma - if reached
            if self.t3s.to(u.s).value > 0:
                ax1.errorbar(self.t3s.to(u.s).value + t_trig_rel.sec,
                             self.alt3s.value,
                             xerr=self.et3s.to(u.s).value,
                             yerr=self.ealt3s.value,
                             fmt='o', color=col3,
                             label=r"$3\sigma$")

                ax1.vlines(self.t3s.to(u.s).value + t_trig_rel.sec,
                           ymin=0,
                           ymax=self.alt3s.value, lw=1, color=col3)

                ax1.hlines(self.alt3s.value,
                           xmin=xmin,
                           xmax=self.t3s.to(u.s).value + t_trig_rel.sec,
                           alpha=0.5, lw=1, color=col3)

            # 5 sigma - if reached
            if self.t5s.to(u.s).value > 0:
                ax1.errorbar(self.t5s.to(u.s).value + t_trig_rel.sec,
                             self.alt5s.value,
                             xerr=self.et5s.to(u.s).value,
                             yerr=self.ealt5s.value,
                             fmt='o', color=col5,
                             label=r"$5\sigma$")

                ax1.vlines(self.t5s.to(u.s).value + t_trig_rel.sec,
                           ymin=0,
                           ymax=self.alt5s.value, lw=1, color=col5)

                ax1.hlines(self.alt5s.value,
                           xmin=xmin,
                           xmax=self.t5s.to(u.s).value + t_trig_rel.sec,
                           alpha=0.5, lw=1, color=col5)

            # Title and limits
            ax1.set_title(text_ref, fontsize=12, loc="right")
            ax1.set_xlim(xmin=xmin, xmax=xmax)
            ax1.set_ylim(ymin=ymin)

            if ref == "VIS":
                ax1.set_xlabel("Elapsed time since visible (s)")
            else:
                ax1.set_xlabel("Elapsed time since Trigger (s)")
            ax1.set_ylabel("Altitude")

            # ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax1.legend(fontsize=10)

            # Display second axis
            ax = ax1.twiny()
            ax.plot(tvis_abs.datetime, altazvis.alt,
                    linestyle="",
                    color="b")
            ax.tick_params(axis='x', rotation=70)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

            if ref == "VIS":
                ax.set_xlabel("Date since twilight (hh:mm)")
            else:
                ax.set_xlabel("Date since Trigger (hh:mm)")

        plt.tight_layout()

        return fig
