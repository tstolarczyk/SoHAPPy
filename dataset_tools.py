# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:51:18 2020

This modules contains the functions to manipulate the Gammapy `dataset` and
`datasets` objects as used in `SoHAPPy`, both for the population analysis
(see `<population_analysis.rst>`_) or the individual source spectra analysis
(see `<spectral_analysis.rst>`_).

@author: Stolar
"""
import sys

import pickle
from pathlib import Path
import warnings

import numpy as np

import astropy.units as u
from astropy.time          import Time
# This causes a pylint crash. Issue opened.
from   astropy.coordinates import SkyCoord

from regions import CircleSkyRegion

from gammapy.utils.random import get_random_state
from gammapy.stats import WStatCountsStatistic

from gammapy.datasets import SpectrumDataset, Datasets, SpectrumDatasetOnOff
from gammapy.makers   import SpectrumDatasetMaker

import gammapy
if gammapy.__version__ < "1.2":
    from gammapy.irf import load_cta_irfs
else:
    from gammapy.irf import load_irf_dict_from_file

from gammapy.data import Observation
from gammapy.maps import RegionNDMap, MapAxis

from gammapy.modeling.models import TemplateSpectralModel
from gammapy.modeling.models import SkyModel

from niceprint import t_fmt
from configuration import Configuration

#warnings.filterwarnings('error')
warnings.filterwarnings('ignore')

sys.path.append("../SoHAPPy")

__all__ = ["generate_dataset","stacked_model","compare_stacked_dataset","stacking",
           "get_masked_dataset","compactify","createonoff_from_simulation",
           "read_from_ogip","get_axis","check_dataset","check_datasets","sigmax"]

#------------------------------------------------------------------------------
def generate_dataset(Eflux, flux, Erange = None,
                     tstart   = Time('2000-01-01 02:00:00',scale='utc'),
                     tobs     = 100*u.s,
                     irf_file = None,
                     alpha    = 1/5,
                     name     = None,
                     fake     = True,
                     onoff    = True,
                     seed     = 'random-seed',
                     debug    = False):
    """
    Generate a dataset from a list of energies and flux points either as
    a SpectrumDataset or a SpectrumDatasetOnOff

    Note :
    - in SpectrumDataset, the backgound counts are assumed precisely known and
    are not fluctuated.
    - in SpectrumDatasetOnOff, the background counts (off counts) are
    fluctuated from the IRF known values.

    Parameters
    ----------
    Eflux : Quantity
        Energies at which the flux is given.
    flux : Quantity
        Flux corresponding to the given energies.
    Erange : List, optional
        The energy boundaries within which the flux is defined, if not over all
        energies. The default is None.
    tstart : Time object, optional
        Start date of the dataset.
        The default is Time('2000-01-01 02:00:00',scale='utc').
    tobs : Quantity, optional
        Duration of the observation. The default is 100*u.s.
    irf_file : String, optional
        The IRf file name. The default is None.
    alpha : Float, optional
        The on over off surface ratio for the On-Off analysis.
        The default is 1/5.
    name : String, optional
        The dataset name, also used to name the spectrum. The default is None.
    fake : Boolean, optional
        If True, the dataset counts are fluctuated. The default is True.
    onoff : Boolean, optional
        If True, use SpectrumDatasetOnOff, otherwise SpectrumDataSet.
        The default is True.
    seed : String, optional
        The seed for the random generator; An integer will generate the
        same random series at each run. The default is 'random-seed'.
    debug: Boolean
        If True, let's talk a bit. The default is False.

    Returns
    -------
    ds : Dataset object
        The dataset.

    """
    random_state = get_random_state(seed)

    ### Define on observation region
    on_pointing = SkyCoord(ra=0*u.deg, dec=0*u.deg, frame="icrs")
    on_region   = CircleSkyRegion(center = on_pointing, radius = 0.5*u.deg)

    # Define energy axis (see spectrum analysis notebook)
    # edges for SpectrumDataset - all dataset should have the same axes
    # Note that linear spacing is clearly problematic for powerlaw fluxes
    # Axes can also be defined using MapAxis
    unit = u.GeV
    E1v  = min(Eflux).to(unit).value
    E2v  = max(Eflux).to(unit).value

    ereco_axis  = MapAxis.from_energy_bounds(1.1*E1v*unit,
                                             0.9*E2v*unit,
                                             nbin = 4,
                                             per_decade=True,
                                             name="energy")

    etrue_axis  = MapAxis.from_energy_bounds(E1v*unit,
                                             E2v*unit,
                                             nbin = 4,
                                             per_decade=True,
                                             name="energy_true")
    if debug:
        print("Dataset ",name)
        print("  Etrue : ", etrue_axis.edges)
        print("  Ereco : ", ereco_axis.edges)

    # Load IRF
    if gammapy.__version__ < "1.2":
        irf   = load_cta_irfs(irf_file)
    else:
        irf   = load_irf_dict_from_file(irf_file)

    spec = TemplateSpectralModel(energy = Eflux,
                                 values = flux,
                                 interp_kwargs={"values_scale": "log"})

    model = SkyModel(spectral_model = spec, name  = "Spec" + str(name))
    obs   = Observation.create(obs_id   = 1,
                               pointing = on_pointing,
                               livetime = tobs,
                               irfs     = irf,
                               deadtime_fraction = 0,
                               reference_time = tstart)

    ds_empty  = SpectrumDataset.create(e_reco = ereco_axis, # Ereco.edges,
                                       e_true = etrue_axis, #Etrue.edges,
                                       region = on_region,
                                       name   = name)
    maker = SpectrumDatasetMaker(containment_correction = False,
                                  selection=["exposure", "background","edisp"])
    ds = maker.run(ds_empty, obs)
    ds.models = model
    mask = ds.mask_safe.geom.energy_mask(energy_min = Erange[0],
                                         energy_max = Erange[1])

    mask = mask & ds.mask_safe.data
    ds.mask_safe = RegionNDMap(ds.mask_safe.geom,data=mask)

    ds.fake(random_state = random_state)   # Fake is mandatory ?

    # Transform SpectrumDataset into SpectrumDatasetOnOff if needed
    if onoff:
        ds = SpectrumDatasetOnOff.from_spectrum_dataset(
                                                        dataset=ds,
                                                        acceptance=1,
                                                        acceptance_off=1/alpha)
        print("Dataset was transformed in ONOFF")

    if fake:
        print(" Fluctuations : seed = ",seed)
        if onoff:
            ds.fake(ds.npred_background())
        else:
            ds.fake(random_state = random_state)

    print("ds.energy_range = ",ds.energy_range)

    return ds

#------------------------------------------------------------------------------
def stacked_model(ds, ds_stack = None, first = False, debug = False):

    """
    This function is not finalised and should be checked. It is intended to
    recompute an effective spectral model from the exiisting already stacked
    model -from the previous call- and the current model in a dataset list. It
    also extract the masked dataset of the first dataset in the list
    (when first == True).

    The reconstructed binning is required to apply the masking
    The theoretical spectrum is therefore evaluated on reconstrcuted
    energies which assumes that the reconstructed energy is not too
    different from the true one.


    Parameters
    ----------
    ds : Dataset
        The current dataset.
    ds_stack : Dataset, optional
        The current stacked dataset if first is False. The default is None.
    first : Boolean, optional
        If True, extract the masked dataset - Valid for the first dataset of
        the list. The default is False.
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    model : Skymodel
        The current stcake model.

    """

    # Get unmasked Reconstructed E bining from current dataset
    # (but they are all identical)
    e_axis = ds.background.geom.axes[0] # E reco sampling from the IRF

    if first:
        # The reconstructed binning is required to apply the masking
        # The theoretical spectrum is therefore evaluated on reconstrcuted
        # energies which assumes that the reconstructed energy is not too
        # different from the true one.

        # e_axis = dsets[0].aeff.data.axes[0] # E true sampling from the IRF
        flux_org = ds.models[0].spectral_model(e_axis.center)
        mask_org = ds.mask_safe.data.flatten()

        spec = TemplateSpectralModel(energy = e_axis.center,
                                     values = flux_org*mask_org,
                                     interp_kwargs ={"values_scale": "log"})
        model = SkyModel(spectral_model = spec, name  = "Stack"+"-"+ds.name)

    else:
        flux_stacked = ds_stack.models[0].spectral_model(e_axis.center)
        dt_stack     = ds_stack.gti.time_sum # Duration on present stack

        flux_org = ds.models[0].spectral_model(e_axis.center)
        mask_org = ds.mask_safe.data.flatten()
        dt       = ds.gti.time_sum # Duration on present stack

        # Create ad-hoc new flux and model
        dt_new    = dt+ dt_stack
        flux_new  = (dt_stack.value*flux_stacked + dt.value*flux_org*mask_org)
        flux_new  = flux_new/(dt_stack.value+dt.value)

        # Create a new SkyModel from the flux template model
        spec = TemplateSpectralModel(energy = e_axis.center,
                                     values = flux_new,
                                     interp_kwargs ={"values_scale": "log"})
        model = SkyModel(spectral_model = spec, name  = "Stack"+"-"+ds.name)

        if debug:
            livetime     = ds.gti.time_sum # Duration on present stack
            print(72*"-")
            print(f" Current dataset dt={livetime:10.2f} -" \
                  f" with model {ds.name,ds.models[0].name:}")
            print("    On stack  : dt=",dt_stack," F0 = ",flux_stacked[0])
            print("    To add    : dt=",dt," F0 = ",flux_org[0])
            print("    To stack  : dt=",dt_new," F0 = ",flux_new[0])
            print("")

    return model
#------------------------------------------------------------------------------
def compare_stacked_dataset(dsets,dsets_stacked,count="pred"):

    """
    Compare the stacked counts of a Dataset list with the counts of a stacked
    Dataset list (they should match).

    Parameters
    ----------
    dsets : Dataset List
        The original Dataset list.
    dsets_stacked : Dataset
        The final stacked dataset.
    count : String, optional
        Defines which cou_nt should be checked among "pred" (from the model),
        "excess" and "background". The default is "pred".

    """

    print(" Check cumulated "+count)
    tot_counts  = 0
    duration    = 0*u.s

    i=0
    print(f"{'Id':>3} {'dt':>8s} {count:>10s} "\
          f"{'dtsum':>8s} {count+'sum':>10s} --> "\
          f"{'dtnew':>8s} {count+'new':>10s}")

    for ds,dss in zip(dsets,dsets_stacked):
        if count=="pred":
            masked_counts  = ds.npred_signal().data[ds.mask_safe].sum()
            stacked_counts = dss.npred_signal().data[dss.mask_safe].sum()
        elif count=="excess":
            masked_counts  = ds.excess.data[ds.mask_safe].sum()
            stacked_counts = dss.excess.data[dss.mask_safe].sum()
        elif count == "background":
            masked_counts  = ds.background.data[ds.mask_safe].sum()
            stacked_counts = dss.background.data[dss.mask_safe].sum()
        tot_counts   += masked_counts

        duration += ds.gti.time_sum
        print(f"{i:3} {ds.gti.time_sum.value:8.2f} "\
              f"{masked_counts:10.2f} {duration.value:8.2f} "\
              f"{tot_counts:10.2f} --> {dss.gti.time_sum.value:8.2f} "\
              f"{stacked_counts:10.2f}")
        i+=1

#------------------------------------------------------------------------------
def stacking(dsets, tryflux=False, debug=False):

    """
    Create a new dataset collection (Datasets) with (consecutively) stacked
    datasets.
    When datasets are stacked they are automatically masked, whihc means that
    the energy bins outside the valid energy range are hidden and not
    considered in the following use of this dataset.

    Note that the first dataset has to be masked explicitely as the stacked
    ones are. Note that by default the model carried by the stacked dataset
    in a stacked list is the model of the first dataset.
    An option exists in gammapy to supersede the models in each consecutively
    stacked dataset but it essentially works only for dataset with the same
    `mask_safe` (See documentation for more explanations).

    Parameters
    ----------
    dsets : List of Dataset objects
        A collection of original datasets.
    tryflux : boolean
        If true a mean spectrum is recomputed an assigned to the stacked models.
        See documentation for caveat.
    debug : boolean
        If True, let's talk a bit. Default is False.

    Returns
    -------
    dsets_stacked : Datasets object
        A collection of stacked datasets.

    """

    #------------------
    def info(ds0,dsstk0):
        print(f" * Current dataset : "\
              f"{ds0.name:20} dt= {t_fmt(ds0.gti.time_sum):10} "\
              f"- Model : {ds0.models[0].name:}")
        print(f" *         Stacked : {dsstk0.name:20} "\
              f"dt= {t_fmt(dsstk0.gti.time_sum):10} "\
              f"- Model : {dsstk0.models[0].name:}")
        print(f"                  -> t1={dsstk0.gti.time_start:} "\
                                   f"t2={dsstk0.gti.time_stop} "\
                                   f"dt={dsstk0.gti.time_delta} ")
    #------------------

    # Put first MASKED dataset on stack - masking is mandatory
    ds = dsets[0]
    stacked = get_masked_dataset(ds.copy(name="1st unmasked"))

    if tryflux: # Change model
        stacked.models = stacked_model(ds, stacked, first=True, debug=debug)
    dsets_stacked = Datasets(ds.copy(name="1st unmasked"))

    if debug:
        info(ds,ds)

    # Stack following ones
    for ds in dsets[1:]:

        stacked.stack(ds) # Stack applies the mask !
        dss = stacked.copy(name="Stacked_"+ ds.name)
        if tryflux:
            dss.models = stacked_model(ds, dsets_stacked[-1], debug=debug)

        # Add the dataset to the stack
        dsets_stacked.append(dss)
        if debug:
            info(ds,dss)

    print(" Initial dataset collection stacked")

    return dsets_stacked
#------------------------------------------------------------------------------
def get_masked_dataset(ds0):

    """
    This is a trick to get a masked dataset from a dataset.

    Parameters
    ----------
    ds0 : Dataset object
        The original complete unmasked dataset.

    Returns
    -------
    masked_dataset : Dataset object
        The original dataset with masking applied.

    """
    # print("dataset_tools/get_maked_dataset(ds) : masking disabled")
    # return ds0
    e_true = ds0.exposure.geom.axes["energy_true"]
    e_reco = ds0.counts.geom.axes[0]
    region = ds0.counts.geom.region
    masked_dataset = SpectrumDataset.create(
        e_true=e_true, e_reco=e_reco, region=region, name=ds0.name
    )
    masked_dataset.models = ds0.models
    masked_dataset.stack(ds0)

    return masked_dataset
#------------------------------------------------------------------------------
def compactify(dsets,dtmin=1*u.h,debug=False):

    """
    Returns a list of stacked Dataset having a minimal total duration from an
    original unstacked Dataset list.
    Note that the model stacking is not applied.

    Parameters
    ----------
    dsets : Dataset List
        The initial list of Datasets.
    dtmin : astropy.time, optional
        The stacked Dataset minimal duration. The default is 1*u.h.
    debug : Boolean, optional
        If True, let's talk a bit. The default is False.

    Returns
    -------
    ds_compacted : Dataset List
        The compacted Dataset list.

    """

    duration     = 0*u.s
    tmp_stack    = Datasets()
    ds_compacted = Datasets()

    iset = 0
    # Loop over the `dataset` objetcs
    for ds in dsets:

        tmp_stack.append(ds)
        duration += ds.gti.time_delta[0]

        if debug:
            print("  ",ds.name," : ",ds.livetime," appended")

        # If max duration reached, stack
        if duration >= dtmin :

            dset_stacked = stacking(tmp_stack, tryflux=False, debug=False)
            name = "Compacted-"+str(iset)
            ds_compacted.append(dset_stacked[-1].copy(name=name))

            if debug:
                print("  Dt exceeded - stack",len(tmp_stack)," datasets")
                print(tmp_stack)
                print("   Duration and stack reset")
                print(dset_stacked)
                print(dset_stacked[-1].name," should be kept as ",name)

            # Reset stack and duration
            duration = 0*u.s
            tmp_stack = Datasets()
            iset+=1

    return ds_compacted
#------------------------------------------------------------------------------
def createonoff_from_simulation(mc, random_state='random-seed', alpha=0.2,
                                debug=False):

    """
    Transform the stored SpectrumDataset list in the MonteCarlo class of
    SoHAPPy into a list of Datasets of type SpectrumDatasetOnOff.
    Note that the On-Off data set is randomised from the initial contents
    ('fake' function).

    Parameters
    ----------
    mc : MonteCarlo instantiation
        A SoHAPPy MonteCarlo class instantiation
    random_state : randome seed
        The default is 'randome-seed'
    debug : boolean, optional
        If True, verbosy. The default is False.

    Returns
    -------
    dlist_onoff : Datasets
        A list of on-off datasets from the originals.

    """

    if mc.slot is None:
        print(" Simulation was not runable - no dataset available")
        return None

    if debug:
        mc.slot.plot()
        print(mc.slot)

    dset_list = mc.dset_list

    # Create on-off datasets fom the original list
    # which is a list of list to account for multiple sites
    dlist_onoff = Datasets()

    for ds_site in dset_list:

        for ds in ds_site:
            # Remove NAN background (corrected in 0.19)
            ds.background.data = np.nan_to_num(ds.background.data )

            # It seems that above 0.17, dataset need to be simulated in order
            # to go to on-off
            ds.fake(random_state = random_state)

            ds_onoff = SpectrumDatasetOnOff.from_spectrum_dataset(
                                            dataset        = ds,
                                            acceptance     = 1,
                                            acceptance_off = 1/alpha)

            ds_onoff.fake(npred_background=ds_onoff.npred_background())

            # if debug:
            #     print(ds_onoff)

            dlist_onoff.append(ds_onoff)

    return dlist_onoff

#------------------------------------------------------------------------------
def read_from_ogip(file=None):

    """
    https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007.pdf
    """
    print(" Not implemented yet, bugs to be solved")

        ### Work on this :
        # # Read the data set - this does not work
        # import os
        # path = "../SoHAPPy/Event85data/grb"

        # filedata = Path(path + "_datasets.yaml")
        # print(" Datastes yaml = ",filedata,end="")
        # if (os.path.isfile(filedata)): print(" exists")
        # else: print(" Not available")

        # filemodel = Path(path + "_models.yaml")
        # print(" Model yaml = ",filemodel,end="")
        # if (os.path.isfile(filemodel)): print(" exists")
        # else: print(" Not available")
        # datasets = Datasets.read(filedata=filedata, filemodel=filemodel)

#------------------------------------------------------------------------------
### UTILITIES
#------------------------------------------------------------------------------
def get_axis(dsets, tunit = u.s, Eunit=u.GeV):

    """
    Returns the list of observing time and the reconstructed energy axis from
    a list od Dataset

    Parameters
    ----------
    dsets : Dataset List
        List of Dataset.
    tunit : astropy.unit, optional
        The unit to express time. The default is u.s.
    Eunit : astropy.unit, optional
        The unit ti express energies. The default is u.GeV.

    Returns
    -------
    Quantity List
        list of oberving times.
    Quantity numpy array
        Numpy array ot the reconstructed energy axis.

    """

    #--- Build t axis
    dt = []
    t = 0
    for ds in dsets:
        t += ds.gti.time_delta[0].to(tunit).value
        dt.append(t)
    dt = np.asarray(dt)

    #--- Get E axis
    Erec = dsets[0].excess.geom.axes[0].center
    Erec = np.asarray(np.log10(Erec.to(Eunit).value))

    return dt*tunit, Erec*Eunit


# ##--------------------------------------------------------------------------------------
def check_dataset(ds, tag="?", e_unit="GeV",
                  masked=False, show_header=False, deeper=False):
    """
    Printout some content of a dataset for verification purposes

    Parameters
    ----------
    ds : Dataset object
        The dataset to be scrutinized.
    e_unit: Quantity Unit
        Energy unit
    tag : String, optional
        A descritpion tag. The default is "?".
    masked : Boolean, optional
        If True, display the contents of the masked energy bins.
        The default is False.
    header : Boolean, optional
        If False, print the table header that was not printed yet.
        The default is True.
    deeper : Boolean, optional
        If True, ,describe the dataset content bin per bin.
        The default is False.

    Returns
    -------
    header : Boolean
        Retruns True if the header was displayed during the call.

    """

    if masked is True:
        mask = ds.mask_safe.data
    else:
        mask = np.asarray(ds.mask_safe.data.size*[[[True]]])

    if show_header is True:
        print(90*"=")
        print(" id         name  livetime           on        bck",end="")
        print("       excs       pred    Model     Flux")
        if deeper is True:
            print("Ecenter (",e_unit,")         Ebin         ",end="")
            print("                                   Mask ",end="")
            print(ds.models[0].spectral_model(1*u.TeV).unit)
        show_header=False # Since it was shown already

#   if ds.counts undefined, set ncounts to np.nan etc.
    if ds.counts is not None: # ds.counts != None crahses in Jupyter
        ncounts     = ds.counts.data[mask]
        ncounts_sum = ncounts.sum()
        nxs         = ds.excess.data[mask]
        nxs_sum     = nxs.sum()
    else:
        ncounts     = np.array(ds.background.data[mask].size*[[[np.nan]]])
        ncounts_sum = np.nan
        nxs         = np.array(ds.background.data[mask].size*[[[np.nan]]])
        nxs_sum     = np.nan

    nbck        = ds.background.data[mask]

    if ds.models is None:
        sys.exit("Please, define a Dataset model")

    npred   = ds.npred_signal().data[mask]
    print(90*"=")
    print(f" {tag:3s}  {ds.name:>8s}    {ds.gti.time_sum:6.2f} "
          f"{ncounts_sum:10.2f} {nbck.sum():10.2f} {nxs_sum:10.2f} "
          f"{npred.sum():10.2f}   {ds.models[0].name:6s} ")

    # If requested go into the energy bins
    if deeper:
        e_center = ds.background.geom.axes[0].center
        e_edges = ds.background.geom.axes[0].edges
        print(90*"=")

        for i,_ in enumerate(ds.background.data[mask]):

            energy = e_center[mask.flatten()][i].to(e_unit) # Current ecenter
            emax   = e_edges[e_edges>energy][0].to(e_unit)
            emin   = e_edges[e_edges<energy][-1].to(e_unit)
            print(f" {energy.value:7.1f}  "\
                  f"[{emin.value:7.1f}, {emax.value:7.1f}]",end="")
#            if ds.counts != None:
            print(f" {ncounts.flatten()[i]:10.2f} {nbck.flatten()[i]:10.2f} "\
                  f"{nxs.flatten()[i]:10.2f} {npred.flatten()[i]:10.2f}" \
                  ,end="")
#            else:
#
#                print(" {:>10s} {:10.2f} {:>10s} {:10.2f}"
#                      .format("-",
#                              nbck.flatten()[i],
#                              "-",
#                              npred.flatten()[i]),end="")
            print(f"  {ds.mask_safe.data[mask].flatten()[i]:^6}  "\
                  f"{ds.models[0].spectral_model(energy).value:8.2e} ")

    return show_header

###---------------------------------------------------------------------------------------
def check_datasets(dsets,masked=False,deeper=False,header=True):

    """
    Printout the content of a dataset collection (Datasets), based on the
    check_dataset function.

    Parameters
    ----------
    dsets : Datasets object
        A collection of datasets to be scrutinized.
    masked : Boolean, optional
        If True, apply mask to the bins. The default is False.
    deeper : Boolean, optional
        If True display the content bin per bin. The default is False.

    Returns
    -------
    None.

    """
    for i,ds in enumerate(dsets):
        header = check_dataset(ds,tag=str(i),
                               show_header = header,
                               deeper=deeper,
                               masked=masked)

#------------------------------------------------------------------------------
def sigmax(dsets, alpha = 0.2, stacked = False):

    """
    Compute the max significance from the cumulated counts in a fluctuated
    dataset. Note: the value can change depending on the random
    generation (this is not comparable to the mean maximal siginificance from
    a full simulation).

    Parameters
    ----------
    dsets : list of Dataset objects
        Input datasets

    Returns
    -------
    sigmax : float
        Significance

    """

    non    = 0
    noff   = 0
    sigmax = -999

    for _ ,ds in enumerate(dsets):
        if stacked : # Data have already been masked
            ns = ds.excess.data.flatten().sum()
            nb = ds.background.data.flatten().sum()
        else: # Data should be masked
            ns = ds.excess.data[ds.mask_safe].flatten().sum()
            nb = ds.background.data[ds.mask_safe].flatten().sum()

        # ns   = ds.npred_signal().data.sum() # don't use this for merged ds !!!!
        # if stacked: # mask applied, and existing mask if for 1st element
        #     nb   = ds.npred_background().data.sum()
        # else:
        #     nb   = ds.npred_background().data[ds.mask_safe].sum()
        on  = ns+nb
        off = nb/alpha
        non  += on
        noff += off
        wstat = WStatCountsStatistic(n_on= non, n_off= noff, alpha= alpha)
        if wstat.sqrt_ts > sigmax:
            sigmax = wstat.sqrt_ts

    return sigmax

#------------------------------------------------------------------------------
if __name__ == "__main__":

    # Test the main function

    debug = True

    # Data folder - get configuration
    folder = Path(r"D:/CTA/SoHAPPy/output/long_1_1000/test_omega/strictmoonveto\strictmoonveto_343")

    # Get configuration file
    cf = Configuration()
    cf.read_from_yaml(Path(folder,"config_ref.yaml"))
    if debug:
        cf.print()
    file = "Event343-South_sim.bin" # 4 slices

    # Open archived data
    infile  = open(Path(folder,file),"rb")
    mc      = pickle.load(infile)
    infile.close()

    # Create a ON-Off list of datasets from the stored simulation
    dlist_onoff = createonoff_from_simulation(mc, alpha=cf.alpha, debug=True)

    print("Done")
