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

from gammapy.utils.random import get_random_state
from gammapy.maps import RegionNDMap, MapAxis

# Avoid deprecation Astropy warnings in gammapy.maps
import warnings
with warnings.catch_warnings():
    from gammapy.data import Observation

    from gammapy.datasets import SpectrumDataset
    from gammapy.makers   import SpectrumDatasetMaker
    
__all__ = ['mc_welcome', 'MonteCarlo']

###########################################################################
def mc_welcome(subarray,log=None):
    """
    Welcome and informative message

    Parameters
    ----------
    subarray : string
        Sub-array name.
    log : file, optional
        pointer to log file. The default is None.

    Returns
    -------
    None.

    """
    log.banner("+================================================================+")
    log.banner("||                     LAUNCHING SIMULATION                     ||")
    log.banner("+================================================================+")
    log.prt("   On-region size     : N: {:5} -  S: {:5}"
            .format(mcf.on_size[subarray["North"]],
            mcf.on_size[subarray["South"]]))
    log.prt("   Offset from center : N: {:5} -  S: {:5}"
            .format(mcf.offset[subarray["North"]],
            mcf.offset[subarray["South"]]))
    log.prt("   Eff. area cont.    : {}".format(mcf.containment))
    log.prt("   Min on/off counts  : {}".format(mcf.nLiMamin))

    return

###############################################################################
class MonteCarlo():
    """
    This class handles a GRB simulation and analysis, either on one of the two
    sites, or simultaneously on the two sites.
    It requires the following additionnal modules :
    - mcsim_config, that contains a certain number of analysis parameters
    (e.g. detection level)
    - mcsim_res : dump the results on screen and in a csv file.
    - mcsim_plot : display plots related to the analysis for each event.
    It makes an aperture photometry analysis of a GRB and is particulary suited
    to a population studies.
    The extraction of the energy spectrum is done through standalone notebooks
    from the save simulation for individual GRBs.
    Note that this class also contains the analysis steps consisting in finding 
    where the 3 sigma, 5 sigma levels are found, and where the maximal 
    signiificance is reached. 
    """
    ###------------------------------------------------------------------------
    def __init__(self,
                 niter  = 1,
                 debug  = 0,
                 fluctuate = True,
                 nosignal  = False,
                 seed = 'random-seed',
                 name = "Unknown"):
        """
        Initialize class members to default values

        Parameters
        ----------
        niter : Integer, optional
            Number of Monte carlo iterations. The default is 1.
        debug : Boolean, optional
            If True, verbosy mode. The default is 0.
        fluctuate : Boolean, optional
            If false, generate one simulation with no fluctuation.
            The default is True.
        nosignal: Boolean, optional
            If True force signal to stricly zero. Default is False.
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
        self.niter     = niter     # Number of trials
        self.fluctuate = fluctuate # Poisson fluctuate the count numbers
        self.nosignal  = nosignal  # Force signal count to zero
        self.slot      = None      # The time slot (slices) of this simulation
        self.name      = name

        # The random state is reinitialised here and would lead to the
        # same sequence for all GRB if the seed is a fixed number
        self.rnd_state =  get_random_state(seed)

        # Data set list
        # Not a Datasets object, just my own list so far
        # Contains nslice arrays of datasets for each site 
        # (i.e. for N or S, 2 for Both)
        self.dset_list = [[]] 

        self.mcerr = -999 # Only for the status messages !
        self.mctime = 0.00

        return

    ###------------------------------------------------------------------------
    def run(self, slot, ana,
                        boost    = True,
                        savedset = False,
                        dump_dir = None):
        """
        Run simulations of the current grb for a given hemisphere.
        A simulation correspond to a series of observations corresponding
        to the GRB time slices.
        This performs an aperure photometry analysis, and is particularly
        suitable for a population analyis.

        Compute and store :
            * Significance values for all slices in the trials;
            * Maximum significance reached for the GRB;
            * Time of maximum siginificance, number of excess and background events;
            * Time to reach 3 or 5 sigma, number of excess and background events.

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

        Returns
        -------
        None.

        """
        
        self.slot   = slot
        self.mctime = time.time() # Starts chronometer
        
        # Prepare to dump the slices if requested and an anomaly was found
        if dump_dir != None:
            (fslice, dump_name) = self.dump_slices(phase="open",dir=dump_dir)
            dump = False

        ### Create list of Datasets,and get all counts, once for all
        self.dset_list = self.create_dataset_list()

        ###############################################
        ### Monte Carlo iterations
        ###############################################

        iMC=1
        print("\n",self.name,": ",end="")

        while(iMC <= self.niter):
            if (iMC <= 10) or (np.mod(iMC,10) == 0):
                print("#",iMC," ",end="")

            # Get running cumulated counts and signifcance for this iteration            
            (non_t, noff_t) = self.aperture_photometry(ana.alpha)
            sigma = ana.fill(iMC, non_t, noff_t) # Update analysis data
            
            if ana.abort and boost:
                print(" Aborted")
                break

            # Dump slice stat if requested
            if (dump_dir != None): # dump slices to track problems
                status = self.dump_slices(iMC  = iMC,
                                          data = [non_t,noff_t, sigma],
                                          file = fslice)
                if not dump: dump= status # If True, dump unchanged

            if self.dbg> 2: self.plot_onetrial(iMC+1)

            iMC += 1 # End of MC loop

        # Close special file for slice dumping
        if (dump_dir !=None) :
            self.dump_slices(phase= "close",
                             dump = dump, name=dump_name, file=fslice)
            
        print() # terminate the iteration counting line

        # Timing
        self.mctime = (time.time() - self.mctime)/self.niter
        self.mcerr  = ana.err # Only for the status messages !
        
        return 

    ###------------------------------------------------------------------------
    def status(self,log=None):

        if self.mcerr < 0: message = "Not possible (GRB not visible)"
        elif self.mcerr> 0 and self.mcerr != self.niter:
            message = "Aborted at trial {:>3d}".format(self.mcerr)
        else: message = "Successful ({:5.3f} s)".format(self.mctime)
         
        log.prt("+" + 14*"-" + "+" + 49*"-" + "+")
        log.prt("| Simulation   |  ===> ",end="")
        log.failure("{:42s}".format(message),end="")
        log.prt("|")
        log.prt("+" + 14*"-" + "+" + 49*"-" + "+")
          
        return        
    
    ###------------------------------------------------------------------------
    def aperture_photometry(self,alpha):
        
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

        Returns
        -------
        sigma : numpy array, float
            One significance per time slice
        nxs : numpy array, float
            One excess count per time slice
        nbck : numy array, float
            One background count per time slice

        """

        non_vs_time   = []
        noff_vs_time  = []
        non = noff = ns = nb = 0

        header = True # for debuging

        # cumulate on and off counts along slices
        for ds_site in self.dset_list:

            # Cumulate on and off counts from potentially several sites
            for ds in ds_site:
                ns   = ds.npred_signal().data[ds.mask_safe].sum()
                nb   = ds.npred_background().data[ds.mask_safe].sum()

                if self.nosignal: ns = 0
                                    
                non  += (ns+nb)
                noff += (nb/alpha)
                
                if (self.dbg>2):
                    if header: print()
                    header = check_dataset(ds,
                                           deeper      = (True if self.dbg>3 else False),
                                           masked      = True,
                                           show_header = header)

            # Fluctuate (summed site if required)
            if (self.fluctuate):
                non   = self.rnd_state.poisson(non)
                noff  = self.rnd_state.poisson(noff)
                
            # Much faster to append list element than array
            # elemenst(because it creates new objects)
            non_vs_time.append(non)
            noff_vs_time.append(noff)

            # End of loop over slices / datasets

        # Access to arrays is much faster than access to lists
        non_vs_time     = np.array(non_vs_time)
        noff_vs_time    = np.array(noff_vs_time)

        return (non_vs_time, noff_vs_time)

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
            model  = self.slot.grb.models[aslice.fid()]
            
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

                array   = perf.subarray
                kzen    = perf.kzen
                on_size = mcf.on_size[array]
                offset  = mcf.offset[array]
                # The on-region is on the GRB
                on_region = CircleSkyRegion(center = self.slot.grb.radec,
                                            radius = on_size)
                # Create the observation - The pointing is not on the GRB
                on_ptg = SkyCoord(self.slot.grb.radec.ra + offset,
                                  self.slot.grb.radec.dec, frame="icrs")

                # with warnings.catch_warnings(): # because of t_trig
                #     warnings.filterwarnings("ignore")
                obs = Observation.create(obs_id   = aslice.idt(),
                                         pointing = on_ptg,
                                         livetime = dt,
                                         irfs     = perf.irf,
                                         deadtime_fraction = 0,
                                         reference_time = tref)

                # Create dataset - correct for containment - add model
                ds_name  = aslice.site()+"-"+str(aslice.idt())+"-"+str(ip)
                
                # Reconstructed energy axis
                # Stacking requires same original binning -> uses largest interval
                # Ensure that all possible edges are in, later apply masking
                # Use the Optimised binning
                # There is a bug in 0.17 (unit not taken into account
                # correctly)preventig from simply writing
                # erec_axis = MapAxis.from_edges(erec_edges,name="energy")
                e_reco = MapAxis.from_edges(mcf.erec_edges[array].to("TeV").value,
                                            unit="TeV",
                                            name="energy",
                                            interp="log")
                e_true = perf.etrue
                
                ds_empty = SpectrumDataset.create(e_reco = e_reco,
                                                  e_true = e_true,
                                                  region = on_region,
                                                  name   = ds_name)

                maker = SpectrumDatasetMaker(
                        selection=["exposure", "background","edisp"])  
                    
                ds = maker.run(ds_empty, obs)

                # Compute containment factor  
                # The PSF being given versus true energy, using the reconstructed energy
                # axis to compute the factor assumes that the reconstructed energy is
                # strictly equals to the true energy, which is certainly not the case at
                # the lowest energies.
                # Maybe this could desserve a specific study.                       
                radii = perf.irf['psf'].containment_radius(energy   = e_reco.center,
                                                           theta    = mcf.offset[array],
                                                           fraction = mcf.containment)[0]
                factor  = (1-np.cos(radii))/(1 - np.cos(mcf.on_size[array]))
                
                # If factor is too large above threshold, error
                idx = np.where((e_reco.center)[np.where(factor>1 )] >= mcf.erec_min[array][kzen])
                if np.size(idx):
                    # Get guilty energies
                    print(" E = ",e_reco.center[idx])
                    print(" R = ",radii[idx].value," max = ",
                                  mcf.on_size[array].value)
                    print(" F = ",factor[idx])
                    sys.exit("{}.py: IRF, Initial region too small"
                             .format(__name__))

                # In Gammapy 0.17, it is mandatory to change the effective
                # area before the model is set, since once the model is set it
                # cannot be changed anymore.
                # This feature (bug) was discussed in Gammapy issue #3016 on
                # Sept 2020.

                # ds.exposure.data   *= mcf.containment
                ds.exposure        *= mcf.containment
                ds.background.data *= factor.value.reshape((-1, 1, 1))
                ds.models           = model
                mask = ds.mask_safe.geom.energy_mask(energy_min = mcf.erec_min[array][kzen],
                                                     energy_max = mcf.erec_max[kzen])
                    
                mask = mask & ds.mask_safe.data
                ds.mask_safe = RegionNDMap(ds.mask_safe.geom,data=mask)

                dset_site.append(ds)

            dset_list.append(dset_site)

        return np.asarray(dset_list)

    
    ###------------------------------------------------------------------------
    def dump_slices(self, phase=None,**kwargs):
        """
        Print out the list of on and off counts and the corresponding
        Li & Ma value. If non and off are less than a minimum value
        return a True status flag for further actions.
        This function can be used to dump all slices (just changing the
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
            from pathlib import Path
            name  = self.slot.grb.id+"-"+ self.slot.site+"_slices.txt"
            name  = Path(Path(kwargs["dir"]),name)
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
            [non,noff,sigma] = kwargs["data"]

            # Check if non or noff below limit
            badon  = np.where( non < mcf.nLiMamin)
            badoff = np.where(noff < mcf.nLiMamin)

            # If an anomaly is found, the status become true, the file will
            # not be deleted, and the anomalous event is written out
            if len(non[badon]) + len(noff[badoff]) : status = True
            else: return False

            # Dump
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

        print(" Saving simulation to file : {}".format(filename))
        return


    ###------------------------------------------------------------------------
    def plot_onetrial(self, itrial):
        """
        
    
        Parameters
        ----------
        self.dset_list : TYPE
            DESCRIPTION.
        slot : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        """
        import matplotlib.pyplot as plt
        from gammapy.estimators import FluxPointsEstimator
        import astropy.units as u
    
        # import sys
        # sys.exit("{}.py: plot_onetrial to be reimplemented",__name__)
        
        print(" One trial plots")

        nplots = len(self.dset_list)
        ncols  = 5
        nrows  = int((nplots)/ncols)+1
        
    
        # Predicted counts
        fig, ax = plt.subplots(ncols=ncols,nrows=nrows,
                               figsize=(20,5*nrows), sharey=True)
        iplot = 0
        import itertools
        
        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            if nplots !=1: ax0 = ax[jrow][icol] if (nrows>1) else ax[icol]
            else: ax0=ax

            if iplot<nplots:    
                # self.dset_list[iplot].plot_counts(ax=ax0)
                for ds in self.dset_list[iplot]:
                    ds.npred().plot(ax=ax0, label=ds.name)
                ax0.legend()
            
            if (jrow<nrows-1):
                ax0.axes.get_xaxis().set_visible(False)
            iplot+=1
            
        fig.suptitle("Trial "+str(itrial))
        fig.tight_layout(h_pad=0)
        plt.show()
    
        fig, ax = plt.subplots(ncols=ncols,nrows=nrows,figsize=(15,5*nrows))
        iplot = 0
        import itertools
        
        for jrow, icol in itertools.product(range(nrows), range(ncols)):

            if nplots !=1: ax0 = ax[jrow][icol] if (nrows>1) else ax[icol]
            else: ax0=ax
                
            if iplot<nplots: 
                for ds in self.dset_list[iplot]:               
                    slc = self.slot.slices[iplot]    
                    e_true = slc.irf()[0].irf["aeff"].data.axes["energy_true"].edges
                    fpe = FluxPointsEstimator(energy_edges=e_true)
                    flux_points = fpe.run(datasets=ds)
        #             flux_points.table["is_ul"] = flux_points.table["ts"] < 4
        #             flux_points.plot(energy_power=2,
        #                       flux_unit="erg-1 cm-2 s-1",
        #                       color="tab:blue",ax=ax0)
        #             spectrum = self.slot.grb.spectra[slc.fid()]
        #             t = slc.tobs().to(u.s)
        #             if (t.value > 3600):
        #                 t = t.to(u.h)
        #             if (t.value > 3600*24):
        #                 t = t.to(u.day)
        #             tobs = str( round(t.value,2)) + str(t.unit)
        #             spectrum.plot(energy_range=dset.energy_range,
        #                               flux_unit='cm-2 s-1 erg-1',
        #                               energy_power=2,
        #                               energy_unit='TeV',
        #                               n_points=10,
        #                               ax=ax0,ls=":",color="red",marker="",
        #                               label=tobs)
        #             ax0.legend()
        #             if (jrow<nrows-1):
        #                 ax0.axes.get_xaxis().set_visible(False)
        #             iplot+=1
        #             if (icol !=0): ax0.set_ylabel(None)
        # fig.tight_layout(h_pad=0)
        # plt.show()
        
        # for dset, slice in zip(self.dset_list, slot.slices):
        #     fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols= 3, figsize=(16,5))
    
        #     # Extract spectrum
        #     # flux_points.to_sed_type("e2dnde").plot_ts_profiles(ax=ax3,
        #     #                                                    cmap="Blues",
        #     #                                                    alpha=0.5)
        #     spectrum = slot.grb.spectra[slice.fid()]
        #     spectrum.plot(energy_range=dset.energy_range,
        #                               flux_unit='cm-2 s-1 erg-1',
        #                               energy_power=2,
        #                               energy_unit='TeV',
        #                               n_points=10,
        #                               ax=ax3,ls=":",color="red",marker="",
        #                               label="Theory")
        #     ax3.legend()
        
 ### This seems very old       
 #    @staticmethod
 #    def plot(simu, target):
 #        """Plot some simulation results"""
 #
 #        import matplotlib.pyplot as plt
 #        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
 #                                       figsize=(10, 5))
 #
 #        # Spectrum plot
 #        energy_range = [0.01 * u.TeV, 100 * u.TeV]
 #        target.model.plot(ax=ax1, energy_range=energy_range,
 #                          label='Model')
 #        plt.text(0.55, 0.65, target.__str__(),
 #                 style='italic', transform=ax1.transAxes, fontsize=7,
 #                 bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
 #        ax1.set_xlim([energy_range[0].value, energy_range[1].value])
 #        ax1.set_ylim(1.e-17, 1.e-5)
 #        ax1.grid(which='both')
 #        ax1.legend(loc=0)
 #
 #        # Counts plot
 #        on_off = simu.on_vector.data.data.value
 #        off = 1. / simu.off_vector.backscal * simu.off_vector.data.data.value
 #        excess = on_off - off
 #        bins = simu.on_vector.energy.lo.value
 #        x = simu.on_vector.energy.nodes.value
 #        ax2.hist(x, bins=bins, weights=on_off,
 #                 facecolor='blue', alpha=1, label='ON')
 #        ax2.hist(x, bins=bins, weights=off,
 #                 facecolor='green', alpha=1, label='OFF')
 #        ax2.hist(x, bins=bins, weights=excess,
 #                 facecolor='red', alpha=1, label='EXCESS')
 #        ax2.legend(loc='best')
 #        ax2.set_xscale('log')
 #        ax2.set_xlabel('Energy [TeV]')
 #        ax2.set_ylabel('Expected counts')
 #        ax2.set_xlim([energy_range[0].value, energy_range[1].value])
 #        ax2.set_ylim([0.0001, on_off.max() * (1 + 0.05)])
 #        ax2.vlines(simu.lo_threshold.value, 0, 1.1 * on_off.max(),
 #                   linestyles='dashed')
 #        ax2.grid(which='both')
 #        ax2.text(0.55, 0.05, simu.__str__(),
 #                 style='italic', transform=ax2.transAxes, fontsize=7,
 #                 bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
 #        plt.tight_layout()   
        return