# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

This module handles the Monte carlo simulation of the detector response and
prepare the data for analysis.

@author: Stolar
"""
import warnings

import sys
import os
import time
import pickle
from pathlib import Path
import itertools

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion
import matplotlib.pyplot as plt

import mcsim_config as mcf

from dataset_tools import check_dataset

from gammapy.utils.random import get_random_state
from gammapy.maps import MapAxis
from gammapy.maps import RegionGeom

# Avoid deprecation Astropy warnings in gammapy.maps
with warnings.catch_warnings():
    from gammapy.data import Observation

    from gammapy.datasets import SpectrumDataset
    from gammapy.makers import SpectrumDatasetMaker

__all__ = ['MonteCarlo']


###############################################################################
class MonteCarlo():
    """
    This class handles a GRB simulation either on one of the two
    sites, or simultaneously on the two sites.
    It makes an aperture photometry count simulation of a GRB (on-off method)
    along the GRB time slices.
    """

    # ##-----------------------------------------------------------------------
    def __init__(self,
                 niter=1,
                 fluctuate=True,
                 nosignal=False,
                 seed='random-seed',
                 name="Unknown",
                 emin=None,
                 emax=None,
                 edense=False,
                 tmin=None,
                 tmax=None,
                 dbg=0):
        """
        Initialize class members to default values

        Parameters
        ----------
        niter : Integer, optional
            Number of Monte carlo iterations. The default is 1.
        fluctuate : Boolean, optional
            If False, generate one simulation with no fluctuation.
            The default is True.
        nosignal: Boolean, optional
            If True force signal to stricly zero. Default is False.
        seed : String or integer, optional
            The value of the seed to obtain the random state. Using a fix
            number will have as consequence to have, for all GRB, the same
            fluctuations generated along the trials.
            This can systematically bias the fraction of iterations reaching 3
            sigma in the first trials, and make the acceleration option
            inoperant. It is also very dangerous if the number of iteration
            is low. The default is 'random-seed'.
        name : String, optional
            Name of the simulation (usually related to the GRB name and the
            site). The default is "Unknown".
        dbg: integer
            Debugging flag. The default is O.

        Returns
        -------
        None.

        """

        self.niter = niter     # Number of trials
        self.fluctuate = fluctuate  # Poisson fluctuate the count numbers
        self.nosignal = nosignal  # Force signal count to zero
        self.slot = None      # The time slot (slices) of this simulation
        self.name = name
        self.dbg = dbg

        # The random state is reinitialised here and would lead to the
        # same sequence for all GRB if the seed is a fixed number
        self.rnd_state = get_random_state(seed)

        # Energy boundaries for specific studies
        self.emin = emin if emin is not None else -np.inf
        self.emax = emax if emax is not None else np.inf
        self.edense = edense

        # Time boundaries for specific studies
        self.tmin = tmin if tmin is not None else 0*u.s
        self.tmax = tmax if tmax is not None else 30*u.d

        # Data set list
        # Not a `Datasets` object, just my own list so far
        # Contains nslice arrays of datasets for each site
        # (i.e. for N or S, 2 for Both)
        self.dset_list = [[]]

        self.mcerr = -999  # Only for the status messages !
        self.mctime = 0.00

    # ##-----------------------------------------------------------------------
    @staticmethod
    def welcome(subarray, log=None):
        """
        Welcome and informative messages.

        Parameters
        ----------
        subarray : dictionnary or string
            Sub-array names or name.
        log : file, optional
            pointer to log file. The default is None.

        Returns
        -------
        None.

        """

        log.banner(f"+{78*'=':78s}+")
        log.banner(f"+{'LAUNCHING SIMULATION':^78s}+")
        log.banner(f"+{78*'=':78s}+")

        log.prt(f"   On-region size             : "
                f"N: {mcf.on_size[subarray['North']]:5}"
                f" -  S: {mcf.on_size[subarray['South']]:5}")

        log.prt(f"   Offset from center         : "
                f"N: {mcf.offset[subarray['North']]:5}"
                f" -  S: {mcf.offset[subarray['South']]:5}")

        log.prt(f"   Eff. area cont.            : {mcf.containment}")
        log.prt(f"   Min 'on' and 'off' counts  : {mcf.nLiMamin}")

    # ##-----------------------------------------------------------------------
    def run(self, slot, ana,
            boost=True,
            dump_dir=None):
        """
        Run simulations of the current grb for a given site.
        A simulation is based on a series of observations related to the GRB
        time slices (slot from class :class:`Slot`).
        It computes the `on` and `off` counts and pass the information to the
        :class:`Analysis` for analysis.

        Parameters
        ----------
        slot : Slot object
            A colection of time slices to be analysed.
        ana : Analyze instance
            A :class:`Analysis` instance
        boost : Boolen, optional
            If True, skip simulations if the firsts are not detected.
            The default is True.
        dump_dir : String, optional
            If defined, will dump information on problematic slices to a text
            file. The default is None.

        Returns
        -------
        None.

        """

        self.slot = slot
        self.mctime = time.time()  # Starts chronometer

        # Prepare to dump the slices if requested and an anomaly was found
        if dump_dir is not None:
            (fslice, dump_name) = self.dump_slices(phase="open", dir=dump_dir)
            dump = False

        # ## Create list of Datasets,and get all counts, once for all
        self.dset_list = self.create_dataset_list()

        # ##############################################
        # ## Monte Carlo iterations
        # ##############################################

        iMC = 1
        if self.dbg > 0:
            print("\n", self.name, ": ", end="")

        while iMC <= self.niter:
            if (iMC <= 10 or np.mod(iMC, 10) == 0) and self.dbg > 0:
                print("#", iMC, " ", end="")

            # Get running cumulated counts
            (non_t, noff_t) = self.onoff_counts_along_time(ana.alpha)

            # Update analysis data and get running significance for dumping
            sigma = ana.fill(iMC, non_t, noff_t)

            if ana.abort and boost:
                print(" Aborted")
                break

            # Dump slice stat if requested
            if dump_dir is not None:  # dump slices to track problems
                status = self.dump_slices(iMC=iMC,
                                          data=[non_t, noff_t, sigma],
                                          file=fslice)
                if not dump:
                    dump = status  # If True, dump unchanged

            if self.dbg > 2 and self.niter < 10:
                self.plot_onetrial(iMC)

            iMC += 1  # End of MC loop

        # Close special file for slice dumping
        if dump_dir is not None:
            self.dump_slices(phase="close",
                             dump=dump, name=dump_name, file=fslice)

        if self.dbg > 0:
            print()  # terminate the iteration counting line

        # Timing
        self.mctime = (time.time() - self.mctime)/self.niter
        self.mcerr = ana.err  # Only for the status messages !

    # ##------------------------------------------------------------------------
    def onoff_counts_along_time(self, alpha):
        """
        Compute on and off counts in the field of view for all slices using
        Gammapy functions and return the cumulated counts in time.

        The counts are summed-up over enegy from the :obj:`SpectrumDataset`
        object list (The :obj:`SpectrumDatasetOnOff` objects are not created).

        Parameters
        ----------
        alpha : float
            1 over the total number of on-off regions.

        Returns
        -------
        non_vs_time : list of float
            Cumulated counts along time slices.
        noff_vs_time : list of float
            Cumulated off-counts along time slices.

        """

        non_vs_time = []
        noff_vs_time = []
        non = noff = ns = nb = 0

        header = True  # for debuging

        # cumulate on and off counts along slices
        for ds_site in self.dset_list:

            # Cumulate on and off counts from potentially several sites
            for ds in ds_site:
                ns = ds.npred_signal().data[ds.mask_safe].sum()
                nb = ds.npred_background().data[ds.mask_safe].sum()

                if self.nosignal:
                    ns = 0

                non += ns+nb
                noff += nb/alpha

                if self.dbg > 2:
                    if header:
                        print()
                    header = check_dataset(ds,
                                           deeper=(self.dbg > 3),
                                           masked=True,
                                           show_header=header)

            # Fluctuate (summed site if required)
            if self.fluctuate:
                non = self.rnd_state.poisson(non)
                noff = self.rnd_state.poisson(noff)

            # Much faster to append list element than array
            # elemenst(because it creates new objects)
            non_vs_time.append(non)
            noff_vs_time.append(noff)

            # End of loop over slices / datasets

        # Access to arrays is much faster than access to lists
        non_vs_time = np.array(non_vs_time)
        noff_vs_time = np.array(noff_vs_time)

        return (non_vs_time, noff_vs_time)

    # ##-----------------------------------------------------------------------
    def create_dataset_list(self):
        """
        Create the :obj:`dataset` list as a list of list to handle more than
        one irf per slice (GRB seen on both sites).

        The pointing is the center of the field of view, and it is the same for
        all time slices as it is considered that the GRB is constantly followed
        within the pre-computed visibililty window. The pointing is displaced
        compared to the GRB position by a distance :code:`mc.offset`, a value
        choosen to give enough space to have :code:`1/mcf.alpha` off-regions of
        radius :code:`mcf.on_size`, the detection area.

        The detection areas, on and off regions, have a size independent of
        the energy, but the effective area and background are mofified to take
        into account the PSF evolution with energy (see the :obj:`irf.py`
        module).

        The set of `dataset` are defined with the same true and reconstructed
        enegy axes so that counts can be added bin per bin. Differences in
        thresholds due to different IRF are taken into account by masking some
        bins (the mask at a certain energy value cancel the whole bin it
        belongs to, see the :obj:`irf.py` module for more explanations)

        * A word on masking
            Direct masking is not made available to users in `Gammapy 0.18.2`.
            It is possible to mask on effective area criteria (minimum value),
            or data counts (:obj:`counts.geom.energy_mask`).
            The direct masking on energy is made as shown below, as a result of
            discussions with Gammapy developpers on the Gammapy slack channel
            (Nov. 2:sup:nd, 2020).
            Note that the :obj:`mask_safe` is not considered in the printout or
            peek functions and has to be obtained by hand (Could be corrected
            in further versions): this is implemented in the
            :func:`dataset_tools.check_dataset`  function.

            *Discussion on masking binning with A. Donath, dec. 3rd, 2020*

            *Applying a mask (Emin, Emax) is simply cancelling the data of the
            energy axis bins outside the mask. Even partially covered bins are
            lost (e.g. bins (E1,E2) with E1<Emin<E2 ).*

            *There could have been a partial recovery of the data (i.e. keeping
            the counts or flux between Emin and E2, but this is not
            what has been chosen.*

            *The influence of this becomes small if Energy bin sizes are small
            enough. However small bins size can give low statititics*

            *A. Donath says it is not a problem as long as the statitical
            computation is done properly : this is in fact approaching the
            unbinned likelihood case. Only for flux points computation, e.g.
            to get a significant point at high energies, bin should have enough
            statistics. This can be done since v0.18.2 using
            :obj:`SpectrumDatasetOnOff.resample_energy_axis()` to re-group the
            data before model fitting*

        * A word on stacking
            In order to be stacked, datasets should be based on the same
            energy axis, which is indeed implemented in SoHAPPy.
            Two observation have usually different IRfs, and in particular
            energy thresholds. The thresholds are take into account by masking,
            as mentionned above.
            Stacking in gammapy is intended to merge observations that were to
            large in data volue to be handled together, and or to sum-up
            consecutive observations of the same object. Stacking 2
            observations results in a longer observation with the same model,
            in particular the same spectral model, with an observation livetime
            being the sum of the two initial observation livetimes.
            It is not intended to sum-up two observations done at the same
            time, with the same livetime.
            As a consequence the gammapy stacking cannot be used to "merge"
            the simultaneous observations of two sites. This would result in a
            livetime larger than the simultaneous observations and would
            therefore shift the observation times in the consecutive time
            slices. There was a tentative to modify this in reducing the
            livetime ds_stacked.livetime /= 2, and multiplying the attached
            spectrum by 2 to keep the predicted count number coherent. But this
            causes problems later in the analysis.
            As a consequence, it is chosen to keep track of each observation
            identity, having a dataset for each time slice on a given site, and
            two datasets for slice in common on two sites (or more).
            The aperture photometry nalysis then simply adds-up the cont
            numbers from the two sites, whereas the spectrum extraction handles
            the situation outside the simulation code, in a separate notebook.

        * Saving datasets
            The choice made here is to save the simulation to disk as a bin
            file which has the advantage to have all information saved.
            An alternative is to save the dataset list only.

            *Discussion on the  slack channel with A. Donath (Nov. 27th, 2020)*

            *For writing datasets in the present use-case there is a bit
            more bookkeeping to do, so that you can read back the model and
            datasets correctly. Here is a minimal example:*

        .. code-block:: python

            from gammapy.datasets import SpectrumDatasetOnOff, Datasets
            from gammapy.modeling.models import PowerLawSpectralModel, SkyModel

            path = "$GAMMAPY_DATA/joint-crab/spectra/hess/"

            obs_1 = SpectrumDatasetOnOff
                    .from_ogip_files(path + "pha_obs23523.fits")
            obs_2 = SpectrumDatasetOnOff
                    .from_ogip_files(path + "pha_obs23592.fits")

            model_1 = SkyModel(PowerLawSpectralModel(), name="model-1",
                               datasets_names=[obs_1.name])
            model_2 = SkyModel(PowerLawSpectralModel(), name="model-2",
                               datasets_names=[obs_2.name])
            obs_1.models = [model_1]
            obs_2.models = [model_2]
            datasets = Datasets([obs_1, obs_2])

            datasets.write("test", prefix="test", overwrite=True)

        *This was something that we started to refactor in v0.17 and are
        still improving. So you have the possibility to “assign” models to
        datasets in a declarative way using the dataset_names keyword.
        The resulting YAML file looks like:*

        .. code-block:: python

            -   name: model-1
                type: SkyModel
                spectral:
                    type: PowerLawSpectralModel
                    parameters:
                    - {name: index, value: 2.0, unit: '', min: .nan, max: .nan,
                       frozen: false, error: 0}
                    - {name: amplitude, value: 1.0e-12, unit: cm-2 s-1 TeV-1,
                       min: .nan, max: .nan, frozen: false, error: 0}
                    - {name: reference, value: 1.0, unit: TeV, min: .nan,
                       max: .nan, frozen: true, error: 0}
                datasets_names: ['23523']
            -   name: model-2
                type: SkyModel
                spectral:
                    type: PowerLawSpectralModel
                    parameters:
                    - {name: index, value: 2.0, unit: '',
                       min: .nan, max: .nan, frozen: false, error: 0}
                    - {name: amplitude, value: 1.0e-12, unit: cm-2 s-1 TeV-1,
                       min: .nan, max: .nan, frozen: false, error: 0}
                    - {name: reference, value: 1.0, unit: TeV, min: .nan,
                       max: .nan, frozen: true, error: 0}
                datasets_names: ['23592']
            covariance: test/test_models_covariance.dat

        *So explicitly mentions for each model component the dataset it is
        assigned to. I think without defining dataset_names the information
        will not be present in the output YAML file and when reading back the
        data and model, the assignment is not correct anymore.
        We will add more documentation on this “global” model handling soon.*

        Returns
        -------
        dset_list : Numpy array of Dataset objects
            A list of datasets (not a collection like a Datasets.

        """

        dset_list = []

        for aslice in self.slot.slices:

            # Note: the spectrum is related to the slice, not the site
            model = self.slot.grb.models[aslice.fid()]

            # The reference time and duration of the observation
            # Note that this is the start of the slice
            # Not necessarilty the point at which the flux and altitude
            # are evaluated
            tref = self.slot.grb.t_trig + aslice.ts1()  # start of slice
            dt = aslice.ts2() - aslice.ts1()

            # Each slice can have more than one IRF since it can be observed
            # from both sites in some cases.
            # Two independent datasets are created
            dset_site = []

            for ip, perf in enumerate(aslice.irf()):

                array = perf.subarray
                kzen = perf.kzen
                on_size = mcf.on_size[array]
                offset = mcf.offset[array]

                # The on-region is on the GRB
                on_region = CircleSkyRegion(center=self.slot.grb.radec,
                                            radius=on_size)

                # Create the observation - The pointing is not on the GRB
                on_ptg = SkyCoord(self.slot.grb.radec.ra + offset,
                                  self.slot.grb.radec.dec, frame="icrs")

                # with warnings.catch_warnings(): # because of t_trig
                #     warnings.filterwarnings("ignore")
                obs = Observation.create(obs_id=aslice.idt(),
                                         pointing=on_ptg,
                                         livetime=dt,
                                         irfs=perf.irf,
                                         deadtime_fraction=0,
                                         reference_time=tref)

                # Create dataset - correct for containment - add model
                ds_name = (aslice.site() + "-"
                           + str(aslice.idt()) + "-" + str(ip))

                # Reconstructed energy axis
                # Stacking requires same original binning -> uses largest
                # interval
                # Ensure that all possible edges are in, later apply masking
                # Use the Optimised binning
                # There is a bug in 0.17 (unit not taken into account
                # correctly) preventig from simply writing
                # erec_axis = MapAxis.from_edges(erec_edges, name="energy")
                if self.edense:
                    e_rec = mcf.erec_spectral.to("TeV").value
                else:
                    e_rec = mcf.erec_sparse.to("TeV").value

                e_reco = MapAxis.from_edges(e_rec,
                                            unit="TeV",
                                            name="energy",
                                            interp="log")
                e_true = perf.etrue

                # gammapy.__version__ < "1.2":
                # ds_empty = SpectrumDataset.create(e_reco=e_reco,
                #                                   e_true=e_true,
                #                                   region=on_region,
                #                                   name=ds_name)

                geom = RegionGeom.create(region=on_region, axes=[e_reco])
                ds_empty = SpectrumDataset.create(geom=geom,
                                                  energy_axis_true=e_true,
                                                  name=ds_name)

                maker = SpectrumDatasetMaker(
                        selection=["exposure", "background", "edisp"])

                ds = maker.run(ds_empty, obs)

                # Compute containment factor
                # Energy approximation:
                # The PSF being given versus true energy, using the
                # reconstructed energy axis to compute the factor assumes that
                # the reconstructed energy is strictly equals to the true
                # energy, which is certainly not the case at the lowest
                # energies. Maybe this could desserve a specific study.

                # gammapy.__version__ < "1.2":
                # This returns a list of infividual quantities,
                # [1°, 2.2°, 4°], initially within a list (i.e. [[a,b,c]],
                # thus justifying taking the first element)
                #     radii = perf.irf['psf']\
                #             .containment_radius(energy=e_reco.center,
                #                                 theta=mcf.offset[array],
                #                                 fraction=mcf.containment)[0]

                radii = perf.irf['psf']\
                            .containment_radius(energy_true=e_reco.center,
                                                offset=mcf.offset[array],
                                                fraction=mcf.containment)

                # Angle computation in numpy:
                # The containment radius return an astropy Quantity in degree.
                # It can be checked that numpy handles correctly the conversion
                # to radian (i.e. np.cos(180*u.deg ) returns the dimensionless
                # Quantity -1). If this would not have been the case, the
                # angles in that formula are small enough so that even not
                # converted, and still in degrees, the result would still be
                # approximately valid. If the angle theta is small, then
                # cos(theta) is close to 1 -theta²/2. The factor simplifies as
                # the ratios of the theta²/2, thus the conversion factor 180/pi
                # does not count.
                factor = (1-np.cos(radii))/(1 - np.cos(mcf.on_size[array]))

                # If factor is too large above threshold, error
                idx = np.where((e_reco.center)[np.where(factor > 1)]
                               >= mcf.erec_min[array][kzen])

                if np.size(idx) != 0:
                    # Get guilty energies
                    print(f" E = {e_reco.center[idx]:}")
                    print(f" R = {radii[idx].value:} - "
                          f"max = {mcf.on_size[array].value}")
                    print(f" F = {factor[idx]}")
                    sys.exit(f"{__name__}.py: IRF, Initial region too small")

                # In Gammapy 0.17, it is mandatory to change the effective
                # area before the model is set, since once the model is set it
                # cannot be changed anymore.
                # This feature (bug) was discussed in Gammapy issue #3016 on
                # Sept 2020.

                # ds.exposure.data   *= mcf.containment
                ds.exposure *= mcf.containment
                ds.background.data *= factor.value.reshape((-1, 1, 1))
                ds.models = model
                mask = ds.mask_safe.geom\
                    .energy_mask(energy_min=max(mcf.erec_min[array][kzen],
                                                min(e_reco.edges)),
                                 energy_max=min(mcf.erec_max[kzen],
                                                max(e_reco.edges))
                                 )

                mask = mask & ds.mask_safe.data
                # gammapy.__version__ <= "0.18.2":
                #    ds.mask_safe = RegionNDMap(ds.mask_safe.geom, data=mask)

                ds.mask_safe = mask
                dset_site.append(ds)

            dset_list.append(dset_site)

        # return np.asarray(dset_list) crashes as there are lists of lists with
        # different size. The error is (found when going to gammapy 1.2):
        # setting an array element with a sequence. The requested array has an
        # inhomogeneous shape after 1 dimensions. The detected shape
        # was (6,) + inhomogeneous part.

        return dset_list

    # ##------------------------------------------------------------------------
    def status(self, log=None):
        """
        Display the simulation status

        Parameters
        ----------
        log : Log instance, optional
            Pointer to the log file. The default is None.

        Returns
        -------
        None.

        """

        if self.mcerr < 0:
            message = "Not possible (GRB not visible)"
        elif self.mcerr > 0 and self.mcerr != self.niter:
            message = f"Aborted at trial {self.mcerr:>3d}"
        else:
            message = f"Successful ({self.mctime:5.3f} s)"

        log.prt("+" + 14*"-" + "+" + 49*"-" + "+")
        log.prt("| Simulation   |  ===> ", end="")
        log.failure(f"{message:42s}", end="")
        log.prt("|")
        log.prt("+" + 14*"-" + "+" + 49*"-" + "+")

    # ##------------------------------------------------------------------------
    def dump_slices(self, phase=None, **kwargs):
        """
        Print out the list of on and off counts and the corresponding
        Li & Ma value. If `non` and `off` are less than a minimum value
        return a True status flag for further actions.

        This function can be used to dump all slices (just changing the
        `nLiMamin` value to a very large value)

        Parameters
        ----------
        phase : String, optional
            Select the phase among "open", "close", or by default writing.
            The default is None.

        **kwargs : TYPE
            Extra arguments for each of the phases.

        Returns
        -------
        Boolean
            File status.

        """

        if phase == "open":
            # Open output file
            name = self.slot.grb.id + "-" + self.slot.site + "_slices.txt"
            name = Path(Path(kwargs["dir"]), name)
            with open(name, "w", encoding="utf-8") as fslice:

                print(f"{'iMC / dt':9s}", end="", file=fslice)

                dt_cumul = 0*u.s
                for s in self.slot.slices:
                    dt_cumul += s.ts2() - s.ts1()
                    print(f"{dt_cumul:12.1f}", end="", file=fslice)
                print(file=fslice)  # Line feed

            print("   ---", name, " opened - header written")
            return (fslice, name)

        elif phase == "close":
            # Close output file - delete file if no anomalues were found
            kwargs["file"].close()
            print("   --- Slice file closed")
            if kwargs["dump"] is False:
                # If no anaomaly detected ever, delete file
                os.remove(kwargs["name"])
                print("   --- No anomaly - file deleted")
            return None

        else:  # Print slice features
            status = False
            fslice = kwargs["file"]
            [non, noff, _] = kwargs["data"]

            # Check if non or noff below limit
            badon = np.where(non < mcf.nLiMamin)
            badoff = np.where(noff < mcf.nLiMamin)

            # If an anomaly is found, the status become true, the file will
            # not be deleted, and the anomalous event is written out
            if len(non[badon]) + len(noff[badoff]):
                status = True
            else:
                return False

            # Dump
            for name in ["non", "noff", "sigma"]:
                d = globals()[name]  # Gives value for the name
                print(f"{kwargs['iMC']:3d} {name:5s}", end="", file=fslice)
                for x in d:
                    print(f"{x.item():14.1f}", end="", file=fslice)
                print(file=fslice)

            return status

    # ##-----------------------------------------------------------------------
    def write(self, filename=None):
        """
        Save the present class to disk for further use in particular a
        spectral analysis

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

        with open(filename, "wb") as outfile:
            pickle.dump(self, outfile)

        print(f" Saving simulation to file : {filename}")

    # ##-----------------------------------------------------------------------
    def plot_onetrial(self, itrial):
        """
        Display information obtained after one Monte Carlo trial.
        To be updated

        Parameters
        ----------
        itrial : integer
            Monte Carlo trial counter

        Returns
        -------
        None.

        """

        # sys.exit("{}.py: plot_onetrial to be reimplemented",__name__)

        print(" One trial plots")

        nplots = len(self.dset_list)
        ncols = min(5, nplots)
        nrows = 1 if ncols <= nplots else int((nplots)/ncols) + 1

        # ## ----------------
        # ## Predicted counts
        # ## ----------------
        fig, ax = plt.subplots(ncols=ncols, nrows=nrows,
                               figsize=(4*ncols, 5*nrows),
                               sharey=True)
        iplot = 0

        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            if nplots != 1:
                ax0 = ax[jrow][icol] if (nrows > 1) else ax[icol]
            else:
                ax0 = ax

            if iplot < nplots:
                # self.dset_list[iplot].plot_counts(ax=ax0)
                for ds in self.dset_list[iplot]:
                    ds.npred().plot(ax=ax0, label=ds.name)

            ax0.grid(which="both", axis="both", alpha=0.5)
            ax0.legend()
            ax0.set_ylabel("Predicted counts")

            # Remove y-axis label except if first column
            if icol != 0:
                # ax0.axes.get_yaxis().set_visible(False)
                ax0.set_ylabel(None)

            # Remove x-axis, label except if last row
            if jrow < nrows-1:
                # ax0.axes.get_xaxis().set_visible(False)
                ax0.set_xlabel(None)
            iplot += 1

        fig.suptitle("Trial " + str(itrial))
        fig.tight_layout(w_pad=0)
