# -*- coding: utf-8 -*-
"""
This code handles the observation windows in CTA North and South sites and
combine the observations common to both sites.

Created on Mon Mar 30 10:36:09 2020

@author: Stolar
"""
import warnings
from copy import deepcopy
from pathlib import Path

import numpy as np

import matplotlib.pyplot as plt
import astropy.units as u

from niceprint import failure, heading, Log
from niceplot import show_noticeable_times

from obs_slice import Slice

__all__ = ["Slot"]


# ##---------------------------------------------------------------------------
class Slot():
    """
    A set of time slices, called `slot` is related to a GRB
    lightcurve given as a set of measurement points, and a visibility,
    namely an observation site.
    Slots are essentially a list of Slices (see the :class:`obs_slice.Slice`
    class).

    The slice number `i` starts at zero, and a slice `i` covers the observation
    from the sampled flux at `i` to the flux at `i+1`.

    Time slots can result from the merging of initial slots, allowing for
    managing overlapping or consecutive observations.

    Slots get a name for clarity, and can be associated to a site, or
    a combination of sites explicitely. They follow the site naming convention
    of slices, which is mandatory for slot merging.

    """

    # -------------------------------------------------------------------------
    def __init__(self, grb, name="naked", loc="?",
                 delay=-1*u.s, opt="end", debug=False):
        """

        Create a naked observation slot consisting of a GRB, naked slices,
        and a flag mentionning if the slot was dressed with physics
        information. The time reference is the GRB trigger time.

        Parameters
        ----------
        grb :  GammaRayBurst instance
            a GRB class instantiation
        name : String, optional
            The slot name. The default is "naked".
        loc : string, optional
            The site to which this slot is associated.Can be 'North', 'South'
            or 'Both'. The default is '?'.
        delay : astropy.time, optional
            The observation delay of the slot. The default is -1*u.s.
        opt : string, optional
            Defines how the observation points are chosen.
            The default is "end".
        debug : Boolean, optional
            If True, will be verbosy. The default is False.

        Returns
        -------
        None.
        """

        self.grb = grb
        self.name = name
        self.loc = loc
        self.delay = delay  # Has so far no sense if merged sites
        self.opt = opt
        self.tstart = -np.inf*u.h  # Negative infinite
        self.tstop = np.inf*u.h    # Positive infinite
        self.phys = False          # Initially not associated to flux values
        self.dbg = debug

        slices = []
        for i in range(1, len(grb.tval)):  # Starts at 1, but uses i-1
            time1 = grb.tval[i-1]
            time2 = grb.tval[i]
            slc = Slice(i-1, time1, time2)
            slices.append(slc)

        self.slices = np.asarray(slices)  # Works also it already an array

    # -------------------------------------------------------------------------
    def dress(self, irf_dir=Path("./"),
              arrays=None,
              zenith=None):
        """
        Dress the slot slices with physics information.

        Parameters
        ----------
        irf_dir : Path, optional
            The IRF folder path name. The default is Path("./").
        arrays : String, optional
            The subarray to which belong the slot. The default is None.
        zenith : astropy.coordinates.Angle, optional
            Observation zenith angle. The default is None.

        Returns
        -------
        None.

        """
        self.phys = True
        for sls in self.slices:
            sls.dress(self.grb,
                      irf_dir=irf_dir, arrays=arrays,
                      opt=self.opt, zenith=zenith)

        # After dressing, two consecutives slices can be associated to the
        # same spectrum id and should be merged ONLY if the IRFs have not
        # been associated yet.
        to_be_compacted = True

        while to_be_compacted:

            fidlist = [slc.fid() for slc in self.slices]
            if len(fidlist) != len(set(fidlist)):
                to_be_compacted = self.compact()  # Returns True or False
            else:
                to_be_compacted = False

    # -------------------------------------------------------------------------
    def compact(self):
        """
        Scan the slot slices, merge those with the same physical flux.

        Returns
        -------
        None.

        """
        prev_fid = -999
        prev_irf_files = [Path("")]
        newslices = []
        idt = 0

        for i, slc in enumerate(self.slices):

            fid = slc.fid()
            irf_files = [resp.filename for resp in slc.irf()]

            # Check if prevois IRFs are identical
            if prev_irf_files == irf_files:  # Same IRFs -> can be merge
                tobemerged = True
            else:  # Not same IRF -> CANNOT be merged
                tobemerged = False

            if fid == prev_fid and tobemerged:
                # Change end-time of prev-spec to this spec
                # print(" Slice ",i," merged with slice ",i-1)
                self.slices[i-1].merge(slc)
            else:
                slc.set_id(idt)
                newslices.append(slc)
                idt += 1
            prev_fid = fid
            prev_irf_files = irf_files

        if len(self.slices) == len(newslices):
            # Nothing was compacted
            return False
        else:
            self.slices = newslices
            return True

    # -------------------------------------------------------------------------
    def copy(self, name=None):
        """
        Copy a slot, change the initial name.

        Parameters
        ----------
        name : Sring, optional
            The name of the slot copy. The default is None.

        Returns
        -------
        slot_copy : Slot
            The slot copy.

        """
        if name is not None:
            self.name = name
        slot_copy = deepcopy(self)

        return slot_copy

    # -------------------------------------------------------------------------
    def apply_visibility(self, site="?", delay=0*u.s):
        """

        Apply visibility window and delay to the slot.

        This changes the time boundaries.

        Parameters
        ----------
        site : String, optional
            The site identifier. The default is "?".
        delay : Quantity, optional
            The delay before detection can start. The default is 0*u.s.

        Returns
        -------
        False if no slice survive to the reorganisation, True otherwise.

        """
        self.loc = site
        unit = self.grb.tval[0].unit

        # ------------------------------------------
        # ## Check if t is within a visibility window
        def visible(time):
            is_ok = False
            for time1, time2 in zip(tstart, tstop):
                if time1 <= time <= time2:
                    is_ok = True
                    continue
            return is_ok
        # ------------------------------------------

        # Get the list of all slices edges
        ticks = []

        for slc in self.slices:  # Original slices
            ticks.append(slc.ts1().value)
            ticks.append(slc.ts2().value)

        # Create visibility start and stop time and store
        shift = True
        tstart = []
        tstop = []
        for tvis in self.grb.vis[self.loc].t_true:
            if len(tvis) == 0:
                continue
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                time0 = (tvis[0] - self.grb.t_trig).sec\
                    + delay.to(unit).value*shift
                time1 = (tvis[1] - self.grb.t_trig).sec
            if time1 < time0:
                print(f" Delay made this slice vanished, "
                      f"duration ={time1-time0}")
            else:
                tstart.append(time0)
                tstop.append(time1)
                ticks.append(time0)
                ticks.append(time1)
                if shift:
                    shift = False

        if len(tstart) == 0:
            # If delays are such that all slices disappear, then
            # the GRB become not visible
            return False

        ticks = list(set(ticks))  # FIRST remove duplicated entries
        ticks.sort()  # THEN Sort (because 'set' change the order)

        # Now loop over all intervals and add to slot
        newslices = []
        idx = 0
        for inb in range(len(ticks)-1):
            time1 = ticks[inb]
            time2 = ticks[inb+1]
            if visible(0.5*(time1 + time2)):
                slc = Slice(idx, time1*unit, time2*unit, site=self.loc)
                idx += 1
                newslices.append(slc)

        # Create a slot from existing and replace the contents
        self.slices = newslices
        self.name = "Visible"
        self.tstart = self.slices[0].ts1()
        self.tstop = self.slices[-1].ts2()
        self.delay = delay  # Loose the delay information after merging

        return True

    # -------------------------------------------------------------------------
    def merge(self, slot):
        """
        Merge current slot (`self`) with another slot.

        Parameters
        ----------
        slot : Slot instance
            The slot to be merged with the current slot

        Returns
        -------
        None.

        """
        # Get the list of all slices edges
        ticks = []
        unit = self.slices[0].ts1().unit

        for obj in [self, slot]:
            for slc in obj.slices:
                ticks.append(slc.ts1().value)
                ticks.append(slc.ts2().value)
        ticks = list(set(ticks))  # FIRST remove duplicated entries
        ticks.sort()  # THEN Sort (because 'set' change the order)

        newslices = []
        icount = 0
        for i in range(len(ticks)-1):
            time1 = ticks[i]
            time2 = ticks[i+1]
            tmid = 0.5*(time1 + time2)
            idxa = self.find(tmid*unit)
            idxb = slot.find(tmid*unit)
            if idxa >= 0:
                if idxb >= 0:
                    loc = "Both"
                else:
                    loc = self.loc
            else:
                if idxb >= 0:
                    loc = slot.loc
                else:
                    loc = None

            # Append only if associated to one or two sites
            if loc is not None:
                slc = Slice(icount, time1*unit, time2*unit, site=loc)
                slc.set_id(icount)
                icount += 1
                newslices.append(slc)

        # Create a slot from existing and replace content
        self.slices = newslices
        self.tstart = min(ticks)*unit
        self.tstop = max(ticks)*unit
        self.loc = "Both"
        self.delay = -1*u.s  # Loose the delay information after merging

    # -------------------------------------------------------------------------
    def find(self, time):
        """
        Find the slice containing the given time.

        Parameters
        ----------
        time : astropy.time (time)
            An observation time

        Returns
        -------
        Index or error code:
            Index of the slice or -1, -2 if not within the slices

        """
        time = time.to(self.slices[0].ts1().unit)
        for slc in self.slices:
            if slc.ts1() <= time <= slc.ts2():
                return slc.idt()

        return -1

    # -------------------------------------------------------------------------
    def plot(self, axa=None, eref=100*u.GeV, **kwargs):
        """
        Display the slices of the time slot.

        Including the measurement points if the slot was dressed.

        Parameters
        ----------
        axa : matplotlib.axes
            Axes on which the plots are shown. The default is None.
        eref : Quantity (energy), optional
            The reference energy used to display the lightcurve.
        **kwargs : keyword arguments
            Extra arguents from matplotlib.

        Returns
        -------
        matplotlib figure
            Current Matplotlib figure
        """
        if axa is None:
            fig, axa = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
        else:
            fig = axa.get_figure()

        # ## GRB lightcurve neareEref ---
        iref = np.argmin(np.abs(eref - self.grb.Eval))
        axa.plot(self.grb.tval, self.grb.fluxval[:, iref],
                 marker="o", ls=":", lw=1, color="grey",
                 markeredgecolor="black",
                 markeredgewidth=1.,
                 markerfacecolor="white",
                 alpha=0.5,
                 label="grb", **kwargs)

        # Compute display interval - takes delay into accont
        # Because a log scale does not accept zero
        tmin = max(1.0, self.slices[0].ts1().value - self.delay.to(u.s).value)
        tmax = self.slices[-1].ts2().value
        axa.set_xlim(0.9*tmin, 1.10*tmax)

        # ## Show delay window
        if self.delay.value >= 0:
            label = f"Delay {self.delay}"
            axa.axvspan(tmin, tmin + self.delay.to(u.s).value,
                        alpha=0.5, color="grey", label=label)

        # ## Slices and flux points
        if self.phys:  # If dressed
            color = {"North": "blue", "South": "red", "Both": "purple"}
            for slc in self.slices:
                tobs = slc.tobs().value
                flux = self.grb.fluxval[slc.fid(), iref].value  # Not absorbed
                ts1 = slc.ts1().value
                ts2 = slc.ts2().value
                site = slc.site()
                axa.plot(tobs, flux,
                         marker="o",
                         markerfacecolor=color[site],
                         markeredgecolor=color[site],
                         markersize=8,
                         ls="",
                         label=site, **kwargs,)

                axa.hlines(xmin=ts1, xmax=ts2,
                           y=flux, ls="-", color="green", **kwargs)

        show_noticeable_times(axa, vpos=0.1*axa.get_ylim()[1],
                              tmax=max(self.grb.tval), size=12)
        axa.set_xlabel("Time since trigger (" + str(self.grb.tval[0].unit)+")")
        axa.set_ylabel("Flux at " + str(eref) + " - "
                       + str(self.grb.fluxval[0][0].unit))
        axa.set_xscale("log")
        axa.set_yscale("log")

        # ## Visibility windows
        tref = self.grb.t_trig

        if self.loc in ("North", "Both"):
            for elt in self.grb.vis["North"].t_true:
                axa.axvspan((elt[0] - tref).sec,
                            (elt[1] - tref).sec,
                            alpha=0.2, color="tab:blue", label="vis. North")

        if self.loc in ("South", "Both"):
            for elt in self.grb.vis["South"].t_true:
                axa.axvspan((elt[0] - tref).sec,
                            (elt[1] - tref).sec,
                            alpha=0.2, color="tab:red", label="vis. South")

        axa.grid("both", linestyle=':')
        handles, labels = axa.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        axa.legend(by_label.values(), by_label.keys())
        axa.set_title(self.grb.id + " - " + self.loc + " - " + self.name)

        return fig

    # -------------------------------------------------------------------------
    def __str__(self):
        """Print out the content of a slice set."""
        title = "Slot "+self.name
        if self.loc not in ("North", "South", "Both"):
            heading(title+" Original")
        else:
            heading(title+" "+self.loc)

        print(f" Observation set {self.name:15s} - site={self.loc:5s}")
        print(f"   Visibility : {self.tstart:8.2f} to {self.tstop:8.2f}")
        print(f"    including : {self.delay:8.2f} of delay")

        show_header = True

        for slc in self.slices:
            slc.print(header=show_header)
            if show_header:
                show_header = False

        return " "

    # -------------------------------------------------------------------------
    def both_sites(self, delay=None, debug=False):
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
        A slot instance or None.
            Resulting Slot or `None`  if no time slices left

        """
        if delay is None:
            delay = {"North": 0*u.s, "South": 0*u.s}

        # Compute the delays to be applied to the sites
        delta = (self.grb.vis["North"].t_true[0][0]
                 - self.grb.vis["South"].t_true[0][0]).sec*u.s

        # Add delay to the first one plus the difference to the second
        if delta.value < 0:
            if debug:
                print("North is before South by ", -delta)
            delay_n = delay["North"]
            delay_s = max(0*u.s, delay["South"] - abs(delta))

        elif delta.value > 0:
            if debug:
                print("South is before North by ", delta)
            delay_n = max(0*u.s, delay["North"] - delta)
            delay_s = delay["South"]

        elif delta.value == 0:
            if debug:
                print("North and South are simultaneaous ", delta)
            delay_n = delay["North"]
            delay_s = delay["South"]

        # Get slots for both sites and merge
        slot_n = self.copy()
        still_vis_n = slot_n.apply_visibility(site="North", delay=delay_n)

        slot_s = self.copy()
        still_vis_s = slot_s.apply_visibility(site="South", delay=delay_s)

        if still_vis_s and still_vis_n:
            slot_n.merge(slot_s)
            return slot_n

        failure(" N and/or S vanished because of delays")
        failure(" Combined analysis not possible")

        return None


# ## --------------------------------------------------------------------------
# ## TESTS
# ##---------------------------------------------------------------------------
def test_visibility(slotet):
    """
    Test visibility.

    Parameters
    ----------
    slotet : Slot instance
        Current time Slot.

    Returns
    -------
    None.

    """
    tvis = [[18*u.h, 40*u.h],
            [3*u.h, 60*u.h],
            [3*u.h, 44*u.h],
            [3*u.h, 39*u.h],
            [16*u.h, 60*u.h],
            [7*u.h, 13*u.h]]
    title = ["Inside",
             "Outside",
             "Before1",
             "Before2",
             "After",
             "Within"]

    for time, label in zip(tvis, title):
        print(f" Apply visibility {time[0]} - {time[1]}")

        _, axa = plt.subplots(nrows=1, ncols=1)

        myset = slotet.copy(name=label)
        myset.mask(time[0], time[1])

        print(myset)

        myset.plot(axa, name=label)
        axa.axvline(time[0].value, ls=":", color="tab:green")
        axa.axvline(time[1].value, ls=":", color="tab:red")

        plt.show()


# ##---------------------------------------------------------------------------
def test_flux(slotet, time=23*u.h):
    """
    Test flux retrieval.

    Parameters
    ----------
    slotet : Slot instance
        Current time Slot.
    time: Astropy Quantity
        Time of evaluation

    Returns
    -------
    None.

    """
    flux = slotet.get_flux(time)
    print(f" flux at {time:8.2f} is {flux:8.3f}")


# ##---------------------------------------------------------------------------
# ## MAIN
# ##---------------------------------------------------------------------------
if __name__ == "__main__":

    # Create time slices collection for a series of GRB.
    # Print them, plot them, test visibility application, merging etc.

    import os
    from grb import GammaRayBurst
    from configuration import Configuration

    # This is required to have the EBL models read from gammapy
    os.environ['GAMMAPY_DATA'] = str(Path(Path(__file__).absolute().parent,
                                          "data"))

    log = Log(talk=True)

    # ## -------------------
    # ## Configuration
    # ## -------------------
    cf = Configuration.command_line()

    # Supersede some parameters
    # 85 : not visible at all
    cf.dbg_level = 2
    cf.ifirst = 343  # 1, 85, 204
    cf.nsrc = 1

    cf.print(log)

    # Input data
    grblist = cf.source_ids(cf.infolder,)
    visinfo = cf.decode_visibility_keyword(cf.create_output_folder())

    # ##---------------------------
    # ## Loop over sources
    # ##---------------------------
    for item, fname in enumerate(grblist):

        grb = GammaRayBurst.from_fits(Path(fname),
                                      prompt=cf.prompt_dir,
                                      ebl=cf.ebl_model)

        for loc in ["North", "South"]:
            grb.set_visibility(item, loc, info=visinfo)

        # Printout grb and visibility windows
        if cf.niter <= 1 or cf.dbg > 0 or cf.ngrb == 1:
            heading(grb.id)
            log.prt(grb)
            grb.vis["North"].print()
            grb.vis["South"].print()

        # Create the original GRB slot
        origin = Slot(grb,
                      opt=cf.obs_point,
                      name=grb.id,
                      debug=bool(cf.dbg > 1))

        print(origin)
        origin.plot(axa=plt.subplots()[1])  # Empty dots (not dressed)

        # ## Compute slots
        delay = cf.get_delay()

        for loc in ["North", "South", "Both"]:

            name = grb.id + "-" + loc
            still_vis = False  # Assumed not visible

            # ## Both sites - create a slot
            if loc == "Both":
                if grb.vis["North"].vis_night and grb.vis["South"].vis_night:
                    slot = origin.both_sites(delay=delay, debug=(cf.dbg > 1))
                    if slot is not None:
                        still_vis = True

            # ## North or South, create a slot
            else:
                if grb.vis[loc].vis_night:  # Apply delays to original slot
                    slot = origin.copy(name="loc")
                    still_vis = slot.apply_visibility(delay=delay[loc],
                                                      site=loc)

            if still_vis:
                # Add IRF feature and run - Note that this can
                # modify the number of slices (merging)
                slot.dress(irf_dir=Path(cf.infolder, cf.irf_dir),
                           arrays=cf.arrays,
                           zenith=cf.fixed_zenith)

                print(slot)
                slot.plot()

    log.close()
    print("c'est fini")
