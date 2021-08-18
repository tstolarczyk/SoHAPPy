# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""

import sys
import numpy as np
import time
import pickle

from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion

import mcsim_config as mcf

from dataset_tools import check_dataset

# Avoid deprecation Astropy warnings in gammapy.maps
from gammapy.modeling.models import SkyModel
from gammapy.utils.random import get_random_state
from gammapy.stats import WStatCountsStatistic
from gammapy.maps import RegionNDMap

import warnings

import gammapy

with warnings.catch_warnings():
    from gammapy.data import Observation

    from gammapy.datasets import SpectrumDataset
    from gammapy.makers   import SpectrumDatasetMaker

__all__ = ['MonteCarlo']

###############################################################################
class MonteCarlo():
    """
    This class handles a GRB simulation and analysis, either on one of the two
    sites, or silmutaneously on the two sites.
    It requires the following additionnal modules :
    - mcsim_config, that contains a certain number of analysis parameters
    (e.g. detection level)
    - mcsim_res : dump the results on screen and in a csv file.
    - mcsim_plot : display plots related to the analysis for each event.
    It makes an aperture phtometry analysis of a GRB and is particulary suited
    to a population studies.
    The extraction of the energy spectrum is done through standalone notebooks
    from the save simualtion for individual GRBs.
    """
    ###------------------------------------------------------------------------
    def __init__(self,
                 niter  = 1,
                 method = 0,
                 debug  = 0,
                 fluctuate = True,
                 seed = 'random-seed',
                 name = "Unknown"):
        """
        Initialize class members to default values

        Parameters
        ----------
        niter : Integer, optional
            Number of Monte carlo iterations. The default is 1.
        method : integer, optional
            Aperture photometry if 0, energy on-off if 1. The default is 0.
        debug : Boolean, optional
            If True, verbosy mode. The default is 0.
        fluctuate : Boolean, optional
            If false, generate one simulation with no fluctuation.
            The default is True.
        seed : String or integer, optional
            The value of the seed to obtain the random state. Using a fix
            number will have as consequence to have, for all GRB, the same
            fluctuations generated along the trials.
            This can systematically bias the fractio of iteration reaching 3
            sigma in the first trials, and make the acceleration option
            inoperant. It is also very dangerous if he number of iteration
            is low.The default is 'random-seed'
        name : String, optional
            Name of the simulation (usually related to the GRB name and the
            site). The default is "Unknown".

        Returns
        -------
        None.

        """
        self.dbg       = debug

        # Input parameters and objects
        self.niter     = niter     # Number of trials
        self.method    = method    # Analysis method
        self.fluctuate = fluctuate # Poisson fluctuate the count numbers
        self.slot      = None      # The time slot (slices) of this simulation
        self.name      = name

        # The random state is reinitialised here and would lead to the
        # same sequence for all GRB if the seed is a fixed number
        self.rnd_state =  get_random_state(seed)

        # Data set list
        self.dset_list = [] # Not a Datasets object, just my own list so far

        # list of simulations (one per slice)
        self.simulations = None # For gammapy 0.12 compatibility

        # Significance over simulations
        self.id_smax_list  = [] # Slice indices to get back the time/ altaz
        self.smax_list     = [] # List of max. significances along trials
        self.nex_smax_list = [] # List of excess counts at max. signif.
        self.nb_smax_list  = [] # List of background counts at max. signif.

        self.id_3s_list    = [] # Slice indices ot get back the time/altaz
        self.nex_3s_list   = [] # List of excess
        self.nb_3s_list    = [] # List of background
        self.detect_3s     = 0  # Number of trials 3 sigma was reached

        self.id_5s_list    = [] # Slice indices to get back the time/altaz
        self.nex_5s_list   = [] # List of excess
        self.nb_5s_list    = [] # List of background
        self.detect_5s     = 0  # Number of trials 5 sigma was reached

        # Mean sigma versus time - one value per time slice
        self.sigma_mean = []
        self.sigma_std  = []

        # Slice number with error or warning
        self.err_slice  = []  # useful ?"

        self.mctime = 0.00
        self.err    = -999 # error code : default, simulation is not completed

        return

    ###------------------------------------------------------------------------
    def run(self, slot, boost    = True,
                        savedset = False,
                        dump_dir = None):
        """
        Run simulations of the current grb, for a given hemisphere.
        A simulation correspond to a series of observations corresponding
        to the GRB time slices.
        This performs an aperure photometry analysis, and is particularly
        suitabel for a popualtion analyis.
        This s not intene
        Compute and store :
            - Significance values for all slices in the trials;
            - Maximum significance reached for the GRB;
            - Time of maximum siginificance, number of excess and background
            events
            - Time to reach 3 or 5 sigma, number of excess and background
            events

        Parameters
        ----------
        slot : Slot object
            A colection of time slices to be analysed
        boost : Boolen, optional
            If True, skip simulations if the first are not detected.
            The default is True.
        savedset : Boolean, optional
            If True, save the dataset to disk. The default is False.
        dump_dir : String, optional
            If defined, will dump information on problematic slices to a text file.
            The default is None.
        debug : Boolean, optional
            If True, verbose mode. The default is False.

        Returns
        -------
        None.

        """

        self.slot   = slot
        self.mctime = time.time() # Starts chronometer
        self.err    = self.niter # GRB visible, simul. ought to be complete

        if (boost): abort_test = False # Abortion test not yet performed
        else: abort_test  = True   # Abortion test supposed already performed

        # Prepare to dump the slices if requested and an anomaly was found
        if dump_dir != None:
            (fslice, dump_name) = self.dump_slices(phase="open",dir=dump_dir)
            dump = False

        ### Create list of Datasets,and get counts,one for all
        self.dset_list = self.create_dataset_list()

        ###############################################
        ### Monte Carlo iterations
        ###############################################
        sigma_sum   = 0           # Compute mean and std of sig for slices
        sigma2_sum  = 0           #     "          "

        iMC=1
        while(iMC <= self.niter):
            if (iMC <= 10) or (np.mod(iMC,10) == 0):
                if (iMC ==1): print("\n",self.name,": ",end="")
                print("#",iMC," ",end="")

            if (self.method == 0): ### Aperture photometry
                (sigma, nxs, nbck) = self.aperture_photometry()
            else:
                sys.exit("Unrecognised analysis method")

            # Acummulate sum and sum**2 for mean / error in each slice
            sigma_sum  += sigma
            sigma2_sum += sigma**2

            # Dump slice stat if requested
            if (dump_dir != None): # dump slices to track problems
                status = self.dump_slices(iMC  = iMC,
                                          data = [nxs,nbck,sigma],
                                          file = fslice)
                if (dump == False): dump= status # If True, dump unchanged

            # Update statistics
            self.fill_stat(sigma, nxs, nbck) # Update stat list

            # In case 3 signma is not reached in the first 10% of the trials
            # then the 90% CL can not be reached.

            if (abort_test == False):
                 if (iMC/self.niter > 1 - mcf.det_level):
                     abort_test = True
                     # If 3 sigma is not reached in the first trials,
                     # the CL will not be reached - stop simulation
                     if (self.detect_3s == 0):
                         self.err = iMC
                         break

            iMC += 1 # End of MC loop

        # Close special file for slice dumping
        if (dump_dir !=None) :
            self.dump_slices(phase= "close",
                             dump = dump, name=dump_name, file=fslice)

        # Timing
        self.mctime = time.time() - self.mctime
        self.mctime /= self.niter

        ### Mean values
        self.sigma_mean = sigma_sum/self.niter
        sigma2_mean     = sigma2_sum/self.niter
        self.sigma_std  = np.sqrt(sigma2_mean-self.sigma_mean**2)

        return

    ###------------------------------------------------------------------------
    def aperture_photometry(self):
        """
        Perform an aperture photometry simulation and analysis.
        The counts are summed-up over enegy from the SpectrumDataset object
        list (The SpectrumDatasetOnOff objects are not created).
        The siginificance is computed under the assmption of a measured
        background, i.e. with possble fluctuation.
        Note that WStatCountsStatistic returns a significance with the
        sign of the excess (negative excess give negative significance).
        If the count number is not fluctuated, the excess can only be
        positive since it is obtained from a physical flux. Excess at zero
        gives a significance at zeo.

        More on the various statistics here :
            https://docs.gammapy.org/0.17/stats/fit_statistics.html

        Parameters
        ----------
        debug : Boolean, optional
            If True, verbose mode. The default is False.

        Returns
        -------
        sigma : numpy array, float
            One significance per time slice
        nxs : numpy array, float
            One excess count per time slice
        nbck : numy array, float
            One background count per time slice

        """

        sigma = []
        nxs   = []
        nbck  = []

        non  = 0
        noff = 0
        sig  = 0
        ns   = 0
        nb   = 0

        header = True # for debuging

        # cumulate on and off counts along slices
        for ds_site in self.dset_list:

            # Cumulate on and off counts from potentially several sites
            for ds in ds_site:
                if gammapy.__version__ == "0.17":
                    ns   = ds.npred_sig().data[ds.mask_safe].sum()
                    nb   = ds.background.data[ds.mask_safe].sum()
                else: # 0.18.2
                    ns   = ds.npred_signal().data[ds.mask_safe].sum()
                    nb   = ds.npred_background().data[ds.mask_safe].sum()
                non  += (ns+nb)
                noff += (nb/mcf.alpha)
                
                if (self.dbg>2):
                    if header: print()
                    header = check_dataset(ds,
                                           deeper      = True,
                                           masked      = True,
                                           show_header = header)

            # Fluctuate (summed site if required)
            if (self.fluctuate):
                non   = self.rnd_state.poisson(non)
                noff  = self.rnd_state.poisson(noff)

            # Compute significance and background/excess counts

            wstat = WStatCountsStatistic(n_on  = non,
                                         n_off = noff,
                                         alpha = mcf.alpha)
            if gammapy.__version__ == "0.17":
                sig =  wstat.significance[0] # This is sigma*sign(nxs)
            else: # 0.18.2
                sig =  wstat.sqrt_ts # ? check
            nb  = mcf.alpha*noff
            ns  = non - nb

            # Much faster to append list than arrays
            # Array appending create a new object
            sigma.append(sig)
            nxs.append(ns)
            nbck.append(nb)

            # End of loop over slices / datasets

        # Access to arrays is much faster than access to lists
        sigma = np.array(sigma)
        nxs   = np.array(nxs)
        nbck  = np.array(nbck)

        return (sigma, nxs, nbck)

    ###------------------------------------------------------------------------
    def create_dataset_list(self):
        """
        Create the dataset list as a list of list to handle more than one
        irf per slice (GRB seen on both sites).

        The pointing, center of the field of view, is the same for all time
        slices as it is considered that the GRB is constantly followed within
        the pre-computed visibililty window. The pointing is displaced compared
        to the GRB position by a distance mc.offset, a value
        choosen to give enough space to have 1/mcf.alpha off regions of radius
        mcf.on_size, the detection area.
        The detection areas, on and off regions, have a size independent of
        the energy, but the effective area and background are mofified to take
        into account the PSF evolution with energy (see the irf module).

        The datasets are defined with the same true and reconstructed enegy
        axes so that counts can be added bin per bin. Differences in threshold
        due to different IRF are taken into account by masking some bins
        (the mask at a certain energy value cancel the whole bin it belongs to,
        see the irf module for more explanations)

        A word on masking
        =================
        Direct masking is not made available to users in Gammapy 0.17:
        It is possible to mask on effective area criteria (minimum value),
        or data counts (counts.geom.energy_mask)
        The direct masking on energy is made as shown below, as a result of
        discussions with Gammapy developpers on the Gammapy slack channel
        (Nov. 2nd, 2020).
        Note that the mask_safe is not considered in the printout or peek
        functions and has to be obtained by hand (Could be corrected in
        further versions): this is implemented in the check_dataset SoHAPPy
        function.

        Discussion on masking binning with A. Donath, dec. 3rd, 2020
        ------------------------------------------------------------
        Applying a mask (Emin, Emax) is simply cancelling the data of the
        energy axis bins outside the mask. Even partially covered bins are
        lost (e.g. bins (E1,E2) with E1<Emin<E2 ).
        There could have been a partial recovery of the data (i.e. keeping
        the counts or flux between Emin and E2, but this is not
        what has been chosen.
        The influence of this becomes small if Energy bin sizes are small
        enough. However small bins size can give low statititics
        A. Donath says it is not a problem as long as the statitical
        computation is done properly : this is in fact approaching the
        unbinned likelihood case. Only for flux points computation, e.g.
        to get a significant point at high energies, bin should have enough
        statistics. This can be done since v0.18.2 using
        SpectrumDatasetOnOff.resample_energy_axis() to re-group the data
        before model fitting

        A word on stacking
        ==================
        In order to be stacked, datasets should be based on the same
        energy axis, which is indeed implemented in SoHAPPy.
        Two observation have usually different IRfs, and in particular energy
        thresholds. The thresholds are take into account by masking, as
        mentionned aboce.
        Stacking in gammapy is intended to merge observations that were to
        large in data volue to be handled together, and or to sum-up
        consecutive observations of the same object. Stacking 2 observations
        results in a longer observation with the same model, in particular the
        same zspectral model, with an observation livetime being the sum of the
        two initial observation livetime. It is not intended to sum-ip two
        observations done at the same time, with the same livetime.
        As a consequence the gammapy stacking cannot be used to "merge"
        the simultaneous observations of two sites. This would result in a
        livetime larger than the simultaneous observations and would therefore
        shift the observation times in the consecutive time slices.
        There was a tentative to modify this in reducing the livetime
        ds_stacked.livetime /= 2, and multiplying the attached spectrum by 2
        to keep the predicted count number coherent. But this causes problems
        later in the analysis.
        As a consequence, it is chosen to keep track of each observation
        identity, having a dataset for each time slice on a given site, and
        two datasets for slice in common on two sites (or more).
        The aperture photometry nalysis then simply adds-up the cont numbers
        from the two sites, whereas the spectrum extraction handles the
        situation outside the simulation code, in a separate notebook.

        Saving datasets
        ===============
        The choice made here is to save the simulation to disk as a bin file
        which has the advantage to have all information saved.
        An alternative is to save the dataset list only.

        Discussion on the  slack channel with A. Donath (Nov. 27th, 2020)
        -----------------------------------------------------------------
        For writing datasets in the present use-case there is a bit
        more bookkeeping to do, so that you can read back the model and
        datasets correctly. Here is a minimal example:
        <code>
        from gammapy.datasets import SpectrumDatasetOnOff, Datasets
        from gammapy.modeling.models import PowerLawSpectralModel, SkyModel

        path = "$GAMMAPY_DATA/joint-crab/spectra/hess/"

        obs_1 = SpectrumDatasetOnOff.from_ogip_files(path + "pha_obs23523.fits")
        obs_2 = SpectrumDatasetOnOff.from_ogip_files(path + "pha_obs23592.fits")

        model_1 = SkyModel(PowerLawSpectralModel(), name="model-1", datasets_names=[obs_1.name])
        model_2 = SkyModel(PowerLawSpectralModel(), name="model-2", datasets_names=[obs_2.name])
        obs_1.models = [model_1]
        obs_2.models = [model_2]
        datasets = Datasets([obs_1, obs_2])

        datasets.write("test", prefix="test", overwrite=True)
        </code>

        This was something that we started to refactor in v0.17 and are
        still improving. So you have the possibility to “assign” models to
        datasets in a declarative way using the dataset_names keyword.
        The resulting YAML file looks like:
        components:

        <code>
        -   name: model-1
            type: SkyModel
            spectral:
                type: PowerLawSpectralModel
                parameters:
                - {name: index, value: 2.0, unit: '', min: .nan, max: .nan, frozen: false,
                    error: 0}
                - {name: amplitude, value: 1.0e-12, unit: cm-2 s-1 TeV-1, min: .nan, max: .nan,
                    frozen: false, error: 0}
                - {name: reference, value: 1.0, unit: TeV, min: .nan, max: .nan, frozen: true,
                    error: 0}
            datasets_names: ['23523']
        -   name: model-2
            type: SkyModel
            spectral:
                type: PowerLawSpectralModel
                parameters:
                - {name: index, value: 2.0, unit: '', min: .nan, max: .nan, frozen: false,
                    error: 0}
                - {name: amplitude, value: 1.0e-12, unit: cm-2 s-1 TeV-1, min: .nan, max: .nan,
                    frozen: false, error: 0}
                - {name: reference, value: 1.0, unit: TeV, min: .nan, max: .nan, frozen: true,
                    error: 0}
            datasets_names: ['23592']
        covariance: test/test_models_covariance.dat
        </code>

        So explicitly mentions for each model component the dataset it is
        assigned to. I think without defining dataset_names the information
        will not be present in the output YAML file and when reading back the
        data and model, the assignment is not correct anymore.
        We will add more documentation on this “global” model handling soon.

        Parameters
        ----------
        debug : Boolean, optional
            If True, verbose mode. The default is False.

        Returns
        -------
        dset_list : Numpy array of Dataset objects
            A list of datasets (not a collection like a Datasets.

        """

        dset_list = []
        
        for aslice in self.slot.slices:

            # Note: the spectrum is related to the slice, not the site
            spec  = self.slot.grb.spectra[aslice.fid()]
            name  = "Spec"+"-"+str(aslice.fid())
            model = SkyModel(spectral_model = spec, name = name)

            # The reference time and duration of the observation
            # Note that this is the start of the slice
            # Not necessarilty the point at which the flux and altitude
            # are evaluated
            tref = self.slot.grb.t_trig + aslice.ts1() # start of slice
            dt   = aslice.ts2() - aslice.ts1()

            # Each slice can have more than one IRF since it can be observed
            # from both sites in some cases.
            # Two independent datasets are created
            dset_site = []
            for ip, perf in enumerate(aslice.irf()):

                on_size = mcf.on_size[perf.subarray]
                offset  = mcf.offset[perf.subarray]
                # The on-region is on the GRB
                on_region = CircleSkyRegion(center = self.slot.grb.radec,
                                            radius = on_size)
                # Create the observation - The pointing is not on the GRB
                on_ptg = SkyCoord(self.slot.grb.radec.ra + offset,
                                  self.slot.grb.radec.dec, frame="icrs")

                with warnings.catch_warnings(): # because of t_trig
                    warnings.filterwarnings("ignore")
                    obs = Observation.create(obs_id   = aslice.idt(),
                                             pointing = on_ptg,
                                             livetime = dt,
                                             irfs     = perf.irf,
                                             deadtime_fraction = 0,
                                             reference_time = tref)

                # Create dataset - correct for containment - add model
                ds_name  = aslice.site()+"-"+str(aslice.idt())+"-"+str(ip)
                
                if gammapy.__version__== "0.17":
                    e_reco = perf.ereco.edges
                    e_true = perf.etrue.edges
                else: # 0.18.2
                    e_reco = perf.ereco
                    e_true = perf.etrue
                    
                ds_empty = SpectrumDataset.create(e_reco = e_reco,
                                                  e_true = e_true,
                                                  region = on_region,
                                                  name   = ds_name)

                if gammapy.__version__ == "0.17":
                    maker = SpectrumDatasetMaker(
                            selection=["aeff", "edisp", "background"])
                else: #0.18.2
                    maker = SpectrumDatasetMaker(
                            selection=["exposure", "background","edisp"])  
                    
                ds = maker.run(ds_empty, obs)

                # In Gammapy 0.17, it is mandatory to change the effective
                # area before the mpdel is set, since once the model is set it
                # cannot be changed anymore.
                # This feature (bug) was discussed in Gammapy issue #3016 on
                # Sept 2020.
                if gammapy.__version__== "0.17":
                    ds.aeff.data.data  *= mcf.containment
                    ds.background.data *= perf.factor
                    ds.models           = model
                    mask = ds.mask_safe.geom.energy_mask(emin = perf.ereco_min,
                                                         emax = perf.ereco_max)
                else: #0.18.2
                    # ds.exposure.data   *= mcf.containment
                    ds.exposure   *= mcf.containment
                    ds.background.data *= perf.factor.reshape((-1, 1, 1))
                    ds.models           = model
                    mask = ds.mask_safe.geom.energy_mask(energy_min = perf.ereco_min,
                                                         energy_max = perf.ereco_max)
                    
                mask = mask & ds.mask_safe.data
                ds.mask_safe = RegionNDMap(ds.mask_safe.geom,data=mask)

                dset_site.append(ds)

            dset_list.append(dset_site)

        return np.asarray(dset_list)

    ###------------------------------------------------------------------------
    def fill_stat(self,sigma, nex, nb):
        """
        Get the statistics and handle the exceptions of the current MC
        simulation

        Parameters
        ----------
        sigma : numpy array, float
            One siginificance per time slice
        nxs : numpy array, float
            One excess count per time slice
        nbck : numy array, float
            One background counter per time slice

        Returns
        -------
        None.

        """

        ### Find maximum - cannot be a slice with non or noff below limit
        sigmax = np.nanmax(sigma) # Returns Nan only if all are Nan
        if np.isnan(sigmax):
            print(" All sigma values are Nan !!! ")
            maxidx = -1
            sigmax = -999
            nexmax = -1
            nbmax  = -1
        else:
            maxidx = np.nanargmax(sigma) # If sigma is nan, it would be the max !
            nexmax = nex[maxidx]
            nbmax  = nb[maxidx]

        self.id_smax_list.append(maxidx)
        self.smax_list.append(sigmax)
        self.nex_smax_list.append(nexmax)
        self.nb_smax_list.append(nbmax)

        # Find where 3 sigma is reached
        nex_3s = -1
        nb_3s  = -1
        id_3s = np.where(sigma>=3)[0] # This ignore Nan values
        if (np.size(id_3s) != 0):
            id_3s  = id_3s[0] # First one
            nex_3s = nex[id_3s]
            nb_3s  = nb[id_3s]
            self.id_3s_list.append(id_3s)
            self.detect_3s += 1

        self.nex_3s_list.append(nex_3s)
        self.nb_3s_list.append(nb_3s)

        # Find where 5 sigma is reached
        nex_5s = -1
        nb_5s  = -1
        id_5s = np.where(sigma>=5)[0]  # This ignore Nan values
        if (np.size(id_5s) != 0):
            id_5s  = id_5s[0] # First one
            nex_5s = nex[id_5s]
            nb_5s  = nb[id_5s]
            self.id_5s_list.append(id_5s)
            self.detect_5s += 1

        self.nex_5s_list.append(nex_5s)
        self.nb_5s_list.append(nb_5s)

        if (self.dbg>2):
            print(" >>> t_smax = {:8.2f}, ({:5.2f} sigma at slice {:2d})"
                  .format(self.slot.slices[maxidx].tobs(),
                          sigma[maxidx],
                          maxidx),end="")
            if (np.size(id_3s) != 0):
                print(" - 3s at t3s = {:8.2f} at slice {:2d})"
                      .format(self.slot.slices[id_3s].tobs(),id_3s))
            else:
                print(" - 3s not reached")
        return

    ###------------------------------------------------------------------------
    def dump_slices(self, phase=None,**kwargs):
        """
        Print out the list of on and off counts and the corresponding
        Li & Ma value. If non and off are less than a minimum value
        return a True status flag for further actions.
        This function ca nbe used to dump all slices (just changing the
        nLiMamin value to a very large value)

        Parameters
        ----------
        phase : String, optional
            Select the phase among "open", "close", or by default writing.
            The default is None.

        **kwargs : TYPE
            Extra arguments for each of the phases.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        if phase== "open":
            # Open output file
            name  = "/" + self.slot.grb.name+"-"+ self.slot.site+"_slices.txt"
            name  = kwargs["dir"]+name
            fslice  = open(name,"w")

            print("{:9s}".format("iMC / dt"),end="",file = fslice)
            import astropy. units as u
            dt_cumul = 0*u.s
            for s in self.slot.slices:
                dt_cumul += s.ts2()-s.ts1()
                print("{:12.1f}".format(dt_cumul), end="",file = fslice)
            print(file=fslice) # Line feed

            print("   ---",name," opened - header written")
            return (fslice, name)

        elif phase=="close":
            # Close output file - delete file if no anomalues were found
            kwargs["file"].close()
            print("   --- Slice file closed")
            if (kwargs["dump"] == False):
                # If no anaomaly detected ever, delete file
                import os
                os.remove(kwargs["name"])
                print("   --- No anomaly - file deleted")
            return

        else: # Print slice features
            status = False
            fslice = kwargs["file"]
            [nxs,nb,sigma] = kwargs["data"]

            non  = nxs + nb
            noff = nb/mcf.alpha

            # Check if non or noff below limit
            badon  = np.where( non < mcf.nLiMamin)
            badoff = np.where(noff < mcf.nLiMamin)

            # If an anomaly is found, the status become true, the file will
            # not be deleted, and the anomalous event is written out
            if len(non[badon]) + len(noff[badoff]) :
                status = True # will dump
                self.err_slice += list(badon[0])+list(badoff[0])
            else:
                return False

            for name in ["non", "noff","sigma"]:
                d = eval(name)
                print("{:3d} {:5s}".format(kwargs["iMC"],name),end="",
                      file = fslice)
                for x in d:
                    print("{:14.1f}".format(x.item()), end="",file = fslice)
                print(file=fslice)

            return status

    ###------------------------------------------------------------------------
    def write(self,filename="None"):
        """
        Save the present class to disk for further use in particular a
        E-dependt on-off analysis

        Parameters
        ----------
        filename : string, optional
            Output file name. The default is "None".

        Returns
        -------
        None.

        """

        if (filename == "None"):
            sys.exit("Output file not defined)")

        outfile  = open(filename,"wb")
        pickle.dump(self,outfile)
        outfile.close()

        from utilities import success
        success(" Saving simulation to file : {}".format(filename))
        return