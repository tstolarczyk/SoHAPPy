# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:51:18 2020

@author: Stolar
"""
import numpy
import numpy as np

import astropy.units as u
from   astropy.time          import Time
from   astropy.coordinates import SkyCoord

from regions import CircleSkyRegion

from utilities import t_fmt

import gammapy

from gammapy.utils.random import get_random_state

from gammapy.datasets import SpectrumDataset, Datasets, SpectrumDatasetOnOff
from gammapy.makers   import SpectrumDatasetMaker

from gammapy.irf import load_cta_irfs
from gammapy.data import Observation
from gammapy.maps import RegionNDMap, MapAxis

from gammapy.modeling.models import TemplateSpectralModel
from gammapy.modeling.models import SkyModel

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
                     debug=False):
    """
    Generate a dataset from a list of energies and flux points either as
    a SpectrumDataset or a SpectrumDatasetOnOff

    Note :
    - in SpectrumDataset, the backgound counts are assumed precisely know and
    are not fluctuated
    - in SpectrumDatasetOnOff, the background counts (off counts) are
    fluctuated from the IRF known value

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
        The dataset name, also used to name tthe spectrum. The default is None.
    fake : Boolean, optional
        If True, the dataset counts are fluctuated. The default is True.
    onoff : Boolean, optional
        If True, use SpectrumDatasetOnOff, otherwise SpectrumDataSet.
        The default is True.
    seed : String, optional
        The sedd for the randome generator; If an integer will generate the
        same random series at each run. The default is 'random-seed'.

    Returns
    -------
    ds : Dataset object
        The dataset.

    """
    random_state = get_random_state(seed)

    ### Define on region
    on_pointing = SkyCoord(ra=0*u.deg,dec=0*u.deg,frame="icrs") # Observing region
    on_region   = CircleSkyRegion(center = on_pointing, radius = 0.5*u.deg)

    # Define energy axis (see spectrum analysis notebook)
    # edges for SpectrumDataset - all dataset should have the same axes
    # Note that linear spacing is clearly problematic for powerlaw fluxes
    # Axes can also be defined using MapAxis
    unit = u.GeV
    E1v  = min(Eflux).to(unit).value
    E2v  = max(Eflux).to(unit).value
#     ereco = np.logspace(np.log10(1.1*E1v), np.log10(0.9*E2v), 20) * unit
#     ereco_axis = MapAxis.from_edges(ereco.to("TeV").value,
#                                    unit="TeV",
#                                    name="energy",
#                                    interp="log")

    ereco_axis  = MapAxis.from_energy_bounds(1.1*E1v*unit,
                                             0.9*E2v*unit,
                                             nbin = 4,
                                             per_decade=True,
                                             name="energy")


#     etrue = np.logspace(np.log10(    E1v), np.log10(    E2v), 50) * unit
#     etrue_axis = MapAxis.from_edges(etrue.to("TeV").value,
#                                    unit="TeV",
#                                    name="energy_true",
#                                    interp="log")
    etrue_axis  = MapAxis.from_energy_bounds(E1v*unit,
                                             E2v*unit,
                                             nbin = 4,
                                             per_decade=True,
                                             name="energy_true")
    if (debug):
        print("Dataset ",name)
        print("Etrue : ", etrue_axis.edges)
        print("Ereco : ", ereco_axis.edges)

    # Load IRF
    irf  = load_cta_irfs(irf_file)

    if gammapy.__version__ == "0.17":
        spec = TemplateSpectralModel(energy = Eflux, 
                                     values= flux, 
                                     norm = 1.,
                                     interp_kwargs={"values_scale": "log"})
    elif gammapy.__version__ == "0.18.2":
         spec = TemplateSpectralModel(energy = Eflux, 
                                      values= flux, 
                                      interp_kwargs={"values_scale": "log"})        

    model = SkyModel(spectral_model = spec, name  = "Spec"+str(name))
    obs   = Observation.create(obs_id   = 1, 
                               pointing = on_pointing,
                               livetime = tobs,  
                               irfs     = irf,
                               deadtime_fraction = 0,
                               reference_time = tstart)
    
    if gammapy.__version__ == "0.17":
        ds_empty = SpectrumDataset.create(e_reco = ereco_axis.edges,
                                          e_true = etrue_axis.edges,
                                          region = on_region,
                                          name   = name)
        maker = SpectrumDatasetMaker(selection=["aeff", "edisp", "background"])
        ds = maker.run(ds_empty, obs)
        ds.models = model
        mask = ds.mask_safe.geom.energy_mask(emin = Erange[0],
                                            emax = Erange[1])

    elif gammapy.__version__ == "0.18.2":
        ds_empty  = SpectrumDataset.create(e_reco = ereco_axis, # Ereco.edges,
                                       e_true = etrue_axis, #Etrue.edges,
                                       region = on_region,
                                       name   = name)
        maker = SpectrumDatasetMaker(containment_correction=False,
                                      selection=["exposure", "background","edisp"]) 
        ds = maker.run(ds_empty, obs)
        ds.models = model
        mask = ds.mask_safe.geom.energy_mask(energy_min = Erange[0],
                                            energy_max = Erange[1])

    mask = mask & ds.mask_safe.data
    ds.mask_safe = RegionNDMap(ds.mask_safe.geom,data=mask)

    if gammapy.__version__ == "0.18.2": # Fake is mandatory ?
        ds.fake(random_state = random_state)    
    
    # Transform SpectrumDataset into SpectrumDatasetOnOff if needed
    if (onoff):


        ds = SpectrumDatasetOnOff.from_spectrum_dataset(
                                    dataset=ds,
                                    acceptance=1,
                                    acceptance_off=1/alpha)
        print("Transformed in ONOFF")

    if fake:
        print(" Fluctuations : seed = ",seed)
        if (onoff):
            if gammapy.__version__ == "0.17":
                ds.fake(random_state=random_state,
                        background_model=ds.background)
            elif gammapy.__version__ == "0.18.2":
                ds.fake(npred_background=ds.npred_background())
                
        else:
            ds.fake(random_state = random_state)

    print("ds.energy_range = ",ds.energy_range)

    return ds

#------------------------------------------------------------------------------
def stacked_model(ds, ds_stack=None, first=False, debug=False):

    # Get unmasked Reconstructed E bining from current dataset (but they are all identical)
    e_axis = ds.background.geom.axes[0] # E reco sampling from the IRF

    if first:
        # The reconstructed binning is required to apply the masking
        # The theoretical spectrum is therefore evalauted on reconstrcuted energies which
        # assumes that the reconstructed energy is not too different from the true one.
        # e_axis = dsets[0].aeff.data.axes[0] # E true sampling from the IRF
        flux_org = ds.models[0].spectral_model(e_axis.center)
        mask_org = ds.mask_safe.data.flatten()
        if gammapy.__version__ == "0.17":
            spec = TemplateSpectralModel(energy = e_axis.center,
                                         values = flux_org*mask_org,
                                         norm   = 1.,
                                         interp_kwargs ={"values_scale": "log"})
        else:#0.18.2
            spec = TemplateSpectralModel(energy = e_axis.center,
                                         values = flux_org*mask_org,
                                         interp_kwargs ={"values_scale": "log"})
        model = SkyModel(spectral_model = spec, name  = "Stack"+"-"+ds.name)

    else:
        flux_stacked = ds_stack.models[0].spectral_model(e_axis.center)
        if gammapy.__version__ == "0.17":
            dt_stack     = ds_stack.livetime # Duration on present stack
        else: # 0.18.2
            dt_stack     = ds_stack.gti.time_sum # Duration on present stack

        flux_org = ds.models[0].spectral_model(e_axis.center)
        mask_org = ds.mask_safe.data.flatten()
        dt       = ds.gti.time_sum # Duration on present stack

        # Create ad-hoc new flux and model
        dt_new    = dt+ dt_stack
        flux_new  = (dt_stack.value*flux_stacked + dt.value*flux_org*mask_org)/(dt_stack.value+dt.value)

        # Create a new SkyModel from the flux template model
        if gammapy.__version__ == "0.17":
            spec = TemplateSpectralModel(energy = e_axis.center,
                                     values = flux_new,
                                     norm   = 1.,
                                     interp_kwargs ={"values_scale": "log"})
        else: #0.18.2
            spec = TemplateSpectralModel(energy = e_axis.center,
                                     values = flux_new,
                                     interp_kwargs ={"values_scale": "log"})
        model = SkyModel(spectral_model = spec, name  = "Stack"+"-"+ds.name)

        if debug:

            livetime     = ds.gti.time_sum # Duration on present stack                  print(72*"-")
            print(" Current dataset dt={:10.2f} - {} with model {}"
                  .format(livetime, ds.name,ds.models[0].name))
            print("    On stack  : dt=",dt_stack," F0 = ",flux_stacked[0])
            print("    To add    : dt=",dt," F0 = ",flux_org[0])
            print("    To stack  : dt=",dt_new," F0 = ",flux_new[0])
            print("")

    return model
#------------------------------------------------------------------------------
def compare_stacked_dataset(dsets,dsets_stacked,count="pred"):

    print(" Check cumulated "+count)
    tot_counts  = 0
    duration = 0*u.s

    i=0
    print("{:>3} {:>8s} {:>10s} {:>8s} {:>10s} --> {:>8s} {:>10s}"
          .format("Id","dt",count,"dtsum",count+"sum","dtnew",count+"new"))
    for ds,dss in zip(dsets,dsets_stacked):
        if count=="pred":
            if gammapy.__version__ == "0.17":
                masked_counts  = ds.npred_sig().data[ds.mask_safe].sum()
                stacked_counts = dss.npred_sig().data[dss.mask_safe].sum()
            else: #0.18.2
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
        print("{:3} {:8.2f} {:10.2f} {:8.2f} {:10.2f} --> {:8.2f} {:10.2f}"
              .format(i,
                      ds.gti.time_sum.value,
                      masked_counts,
                      duration.value,
                      tot_counts,
                      dss.gti.time_sum.value,
                      stacked_counts,
                          ))
        i+=1
    return
#------------------------------------------------------------------------------
def stacking(dsets, tryflux=False, debug=False):
    """
    Create a new dataset collection (Datasets) with (consecutively) stacked datasets.

    Note that the first dataset has to be masked explicitely as the stacked
    ones are. Note that by default the model carried by the stacked dataset
    is the model of the first dataset.
    An option exists to supersede the models in each consecutively stacked
    dataset but it essentially work only for dataset with the same mask_safe
    (See documentation for more explanations)

    Parameters
    ----------
    dsets : Datasets object
        A collection of original datasets.
    tryflux : boolean
        If true a mean spectrum is recomputed an assigned to the stacked models. See documentation for caveat.
    debug : boolean
        If True, verbosy

    Returns
    -------
    dsets_stacked : Datasets object
        A collection of stacked datasets.

    """

    #------------------
    def info(ds0,dss0):
        print(" * Current dataset : {:20} dt= {:10} - Model : {}"
              .format(ds0.name, t_fmt(ds0.gti.time_sum),ds0.models[0].name))
        print(" *         Stacked : {:20} dt= {:10} - Model : {}"
              .format(dss0.name, t_fmt(dss0.gti.time_sum),dss0.models[0].name))
        print("                  -> t1={} t2={} dt={} "
              .format(dss0.gti.time_start,dss0.gti.time_stop,dss0.gti.time_delta))
#         print("        stacked    ",dplt.t_fmt(dss0.gti.time_sum)," - ",dss0.name,
#               " Model",dss0.models[0].name)
    #------------------

    # Put first MASKED dataset on stack
    ds = dsets[0]
    stacked = get_masked_dataset(ds.copy(name="1st unmasked"))
    if tryflux: # Change model
        stacked.models = stacked_model(ds, stacked, first=True, debug=debug)
#     dsets_stacked = Datasets(stacked.copy(name="1st unmasked"))
    dsets_stacked = Datasets(ds.copy(name="1st unmasked"))

    if (debug): info(ds,ds)

    # Stack following ones
    for ds in dsets[1:]:

        stacked.stack(ds) # Stack applies the mask !
        dss = stacked.copy(name="Stacked_"+ ds.name)
        if tryflux:
            dss.models = stacked_model(ds, dsets_stacked[-1], debug=debug)

        # Add the dataset to the stack
        dsets_stacked.append(dss)
        if debug: info(ds,dss)

#     print(" Initial dataset collection stacked")

    return dsets_stacked
#------------------------------------------------------------------------------
def get_masked_dataset(ds0):
    """
    This is a trick to get a masked datset from a dataset.

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

    if gammapy.__version__ == "0.17":
        e_true = ds0.aeff.energy.edges
        e_reco = ds0.counts.geom.axes[0].edges
    else: # 0.18.2
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

    duration = 0*u.s
    tmp_stack    = Datasets()
    ds_compacted = Datasets()
    iset = 0

    for ds in dsets:
        tmp_stack.append(ds)
        duration += ds.gti.time_delta[0]

        if debug: print("  ",ds.name," : ",ds.livetime," appended")

        " If max duration reached, stack"
        if (duration >= dtmin):
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
def createonoff_from_simulation(mc, random_state='rendom-seed', debug=False):
    """


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

    import sys
    sys.path.append("../SoHAPPy")
    import mcsim_res as mcres
    import mcsim_config as mcf
    from utilities    import Log

    if (mc.slot == None):
        print(" Simulation was not runable - no dataset available")
        return None

    if debug:
        # story(mc, ref="VIS",saveplots="False",outfile="")
        # story(mc, ref="notVIS",saveplots="False",outfile="")
        # stat(mc,saveplots="False",outfile="")
        log = Log(name  = "out.log", talk=True)
        mcres.result(mc,mc.slot.grb,log=log)
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
            if gammapy.__version__ != "0.17": # 0.18.2
                ds.fake(random_state = random_state)

            ds_onoff = SpectrumDatasetOnOff.from_spectrum_dataset(
                                            dataset        = ds,
                                            acceptance     = 1,
                                            acceptance_off = 1/mcf.alpha)

            if gammapy.__version__ == "0.17":
                ds_onoff.fake(background_model=ds_onoff.background,
                              random_state = random_state)
            else: # 0.18.2
                ds_onoff.fake(npred_background=ds_onoff.npred_background())

            if (debug): print(ds_onoff)

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

    return
#------------------------------------------------------------------------------
### UTILITIES
#------------------------------------------------------------------------------
def get_axis(dsets, tunit = u.s, Eunit=u.GeV, debug=False):

    #--- Build t axis
    dt = []
    t = 0
    for ds in dsets:
        if gammapy.__version__ == "0.17":
            t += ds.livetime.to(tunit).value
        else: #0.18.2
            t += ds.gti.time_delta[0].to(tunit).value

        dt.append(t)
    dt = np.asarray(dt)

    #--- Get E axis
    Erec = dsets[0].excess.geom.axes[0].center
    Erec = np.asarray(np.log10(Erec.to(Eunit).value))

    return dt*tunit, Erec*Eunit

###---------------------------------------------------------------------------------------
def check_dataset(ds, tag="?", e_unit="GeV",
                  masked = False, show_header=False, deeper=False):
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

    if (masked==True):
        mask = ds.mask_safe.data
    else:
        mask = np.asarray(ds.mask_safe.data.size*[[[True]]])

    if (show_header==True):
        print(90*"=")
        print(" id         name  livetime           on        bck",end="")
        print("       excs       pred    Model     Flux")
        if (deeper==True):
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

    if (ds.models == None):
        import sys
        sys.exit("Please, define a Dataset model")
    if gammapy.__version__ == "0.17":

        npred   = ds.npred_sig().data[mask]
        print(90*"=")
        print(" {:3s}  {:>8s}    {:6.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}   {:6s} "
          .format(tag, ds.name, ds.livetime,
                  ncounts_sum, nbck.sum(), nxs_sum,
                  npred.sum(), ds.models[0].name
                  ))

    else: #0.18.2
        npred   = ds.npred_signal().data[mask]
        print(90*"=")
        print(" {:3s}  {:>8s}    {:6.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}   {:6s} "
              .format(tag, ds.name, ds.gti.time_sum,
                      ncounts_sum, nbck.sum(), nxs_sum,
                      npred.sum(), ds.models[0].name
                      ))

    # If requested go into the energy bins
    if (deeper):
        e_center = ds.background.geom.axes[0].center
        e_edges  = ds.background.geom.axes[0].edges
        print(90*"=")

        for i,_ in enumerate(ds.background.data[mask]):

            energy = e_center[mask.flatten()][i].to(e_unit) # Current ecenter
            emax   = e_edges[e_edges>energy][0].to(e_unit)
            emin   = e_edges[e_edges<energy][-1].to(e_unit)
            print(" {:7.1f}  [{:7.1f}, {:7.1f}]"
                  .format(energy.value, emin.value, emax.value),end="")
#            if ds.counts != None:
            print(" {:10.2f} {:10.2f} {:10.2f} {:10.2f}"
                  .format(ncounts.flatten()[i],
                          nbck.flatten()[i],
                          nxs.flatten()[i],
                          npred.flatten()[i]),end="")
#            else:
#
#                print(" {:>10s} {:10.2f} {:>10s} {:10.2f}"
#                      .format("-",
#                              nbck.flatten()[i],
#                              "-",
#                              npred.flatten()[i]),end="")
            print("  {:^6}  {:8.2e} "
                  .format(ds.mask_safe.data[mask].flatten()[i],
                          ds.models[0].spectral_model(energy).value))

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
    return



#------------------------------------------------------------------------------
if __name__ == "__main__":

    import pickle
    from pathlib import Path

    import warnings
    #warnings.filterwarnings('error')
    warnings.filterwarnings('ignore')

    # Read data from a simulation on disk
    folder = r"D:\000_Today\testnew"
    #file = "Event3-North.bin" # 4 slices
    #file = "Event3-South.bin" # 4 slices
    file = "Event85-South.bin" # 4 slices
    #file = "Event815-North.bin" # 24 slices

    infile  = open(Path(folder,file),"rb")
    mc      = pickle.load(infile)
    infile.close()

    dlist_onoff = createonoff_from_simulation(mc,debug=False)

    # Actions tests

    print("Done")