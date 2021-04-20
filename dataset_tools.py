# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:51:18 2020

@author: Stolar
"""


import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from   astropy.time          import Time
from   astropy.coordinates import SkyCoord
from   astropy.visualization import quantity_support

from regions import CircleSkyRegion

from gammapy.utils.random import get_random_state

from gammapy.datasets import SpectrumDataset, Datasets, SpectrumDatasetOnOff
from gammapy.makers   import SpectrumDatasetMaker

from gammapy.irf import load_cta_irfs
from gammapy.data import Observation
from gammapy.maps import RegionNDMap, MapAxis

from gammapy.modeling.models import TemplateSpectralModel
from gammapy.modeling.models import SkyModel

#------------------------------------------------------------------------------
def generate_dataset(Eflux,flux,Erange=None,
                     tstart   = Time('2000-01-01 02:00:00',scale='utc'),
                     tobs     = 100*u.s,
                     irf_file =None,
                     alpha    = 1/5,
                     name     = None,
                     fake     = True,
                     onoff    = True,
                     seed     = 'random-seed'):
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

    ### Define on region
    on_pointing = SkyCoord(ra=0*u.deg,dec=0*u.deg,frame="icrs") # Observing region
    on_region   = CircleSkyRegion(center = on_pointing, radius = 0.5*u.deg)

    # Define energy axis
    # edges for SpectrumDataset - all dataset should have the same axes
    Etrue    = MapAxis.from_edges(np.linspace(20,260,13),unit=u.GeV,name="Etrue") # Eff. area table
    Ereco    = MapAxis.from_edges(np.linspace(30,300,10),unit=u.GeV,name="energy") # Count vector
    # Etrue    = MapAxis.from_edges(np.linspace(20,260,39)*u.GeV)
    # Ereco    = MapAxis.from_edges(np.linspace(30,300,100)*u.GeV)

    print("Dataset ",name)
    print(" E reco edges : ",Ereco.edges)
    print(" E true edges : ",Etrue.edges)

    # Load IRF
    irf = load_cta_irfs(irf_file)

    spec  = TemplateSpectralModel(energy = Eflux, values= flux, norm = 1.,
                                    interp_kwargs={"values_scale": "linear"})
    model = SkyModel(spectral_model = spec, name  = "Spec"+str(name))
    obs   = Observation.create(obs_id = 1, pointing = on_pointing,
                               livetime = tobs,  irfs = irf,
                               deadtime_fraction = 0,
                               reference_time = tstart)

    ds    = SpectrumDataset.create(e_reco = Ereco.edges,
                                   e_true = Etrue.edges,
                                   region = on_region,
                                   name   = name)
    maker = SpectrumDatasetMaker(selection=["aeff", "edisp", "background"])
    ds = maker.run(ds, obs)
    ds.models = model

    if (Erange != None): # Create an energy mask
        mask = ds.mask_safe.geom.energy_mask(emin = Erange[0],
                                             emax = Erange[1])
        mask = mask & ds.mask_safe.data
        ds.mask_safe = RegionNDMap(ds.mask_safe.geom,data=mask)

    if (onoff):
        ds = SpectrumDatasetOnOff.from_spectrum_dataset(dataset=ds,
                                                        acceptance=1,
                                                        acceptance_off=1/alpha)

    if fake:
        print(" Fluctuations : seed = ",seed)
        if (onoff):
            ds.fake(random_state=get_random_state(seed),
                    background_model=ds.background)
        else:
            ds.fake(random_state = get_random_state(seed))

    return ds



#------------------------------------------------------------------------------
def stacking(dsets):
    """
    Create a new dataset collection (Datasets) with consecutively stacked datasets
    Note that the first dataset has to be masked explicitely as the stacked
    ones are. Note that by default the model carried by the stacked dataset
    is the model of the first dataset.

    Parameters
    ----------
    dsets : Datasets object
        A collection of original datasets.

    Returns
    -------
    dsets_stacked : Datasets object
        A collection of stacked datasets.

    """
    # Initialise with first one
    # Warning : if the trick in not applied, the first dataset is not masked
    stacked = get_masked_dataset(dsets[0].copy(name="1st unmasked"))
    dsets_stacked=Datasets(stacked.copy(name="1st unmasked"))

    # Stack following ones
    for ds in dsets[1:]:
        stacked.stack(ds) # Stack applies the mask !
        dsets_stacked.append(stacked.copy(name="Stacked_"+ds.name))

    return dsets_stacked

# def cumulate(dsets, debug=False):
#     """
#     Returns a collection of stacked Gammapy datasets (Datasets) from an
#     initial colection.
#     In the original Gammapy 0.17 code, the models are not stacked.
#     Done by hand here.

#     Parameters
#     ----------
#     dsets : Datasets
#         A collection of dataset objects

#     Returns
#     -------
#     ds_stacked : Datasets
#         A collection of stacked dataset objects

#     """

#     # Copy first unchanged dataset
#     stacked    = dsets[0].copy()

#     # Get E bining from first dataset (but they are all identical)
#     ecenter    = dsets[0].aeff.data.axes[0].center # Original E sampling

#     # Init. new Dataset collection from the first dataset
#     ds_stacked = Datasets([dsets[0].copy(name=dsets[0].name)])


#     if debug:
#         check_dataset(ds_stacked[0],"First")

#     # Loop over subsequent datasets and stack
#     flux = np.zeros(len(ecenter))

#     for i,ds in enumerate(dsets[1:]):

#         check_dataset(ds,str(i+1))

#         #---------------------
#         # Create stacked model
#         #---------------------
#         # Original flux points to be modified for stacking
#         flux_org = ds.models[0].spectral_model(ecenter)
#         dt       = ds.livetime      # Duration of the element tobe stacked

#         # Flux of current stacked datasets
#         flux_stacked = stacked.models[0].spectral_model(ecenter)
#         dt_stack = stacked.livetime # Duration on present stack

#         # Build a flux weighted by livetime
#         flux_new  = (dt_stack.value*flux_stacked
#                      + dt.value*flux_org)/(dt_stack.value+dt.value)

#         print(" Stacked so far",flux_stacked[0])
#         print(" Current       ",flux_org[0])
#         print(" Weigted sum   ",flux_new[0])
#         # Create a new SkyModel from the flux templat emodel
#         spec = TemplateSpectralModel(energy = ecenter,
#                                      values = flux_new,
#                                      norm= 1.,
#                                      interp_kwargs={"values_scale": "linear"})
#         stacked_model = SkyModel(spectral_model = spec,
#                                   name  = "Stacked"+"-"+ds.name)
#         #---------------------
#         # Add to current stack
#         #---------------------
#         stacked.stack(ds) # Stack current ds to stacked list-model unchanged
#         stacked.models = stacked_model # Change model

#         # Check dataset to be added to the Datasets
#         check_dataset(stacked,"New")


#         if debug:
#             print("Dt= {:8.2f}    Model: {:6s} --- Stacked --> Dt= {:8.2f}    Model: {:6s}"
#               .format(ds.livetime,ds.models[0].name,stacked.livetime,stacked.models[0].name))

#         ds_stacked.append(stacked.copy(name="Stacked-"+ds.name))


#     if (debug):
#         print(" Initial dataset ")
#         tot_npred = 0
#         duration = 0*u.s
#         i=0
#         for ds,dss in zip(dsets,ds_stacked):
#             npred_masked = ds.npred_sig().data[ds.mask_safe].sum()
#             tot_npred   += npred_masked
#             duration += ds.livetime
#             print("{:3} dt={:8.2f} DT={:8.2f} n={:10.2f} N={:10.2f} -- dt={:8.2f} N={:10.2f}"
#                   .format(i,
#                           ds.livetime,
#                           duration,
#                           npred_masked,
#                           tot_npred,
#                           dss.livetime,
#                           dss.npred_sig().data[dss.mask_safe].sum()
#                           ))
#             i+=1
#         print("  ")

#     return  ds_stacked

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

    e_true = ds0.aeff.energy.edges
    e_reco = ds0.counts.geom.axes[0].edges
    region = ds0.counts.geom.region
    masked_dataset = SpectrumDataset.create(
        e_true=e_true, e_reco=e_reco, region=region, name=ds0.name
    )
    masked_dataset.models = ds0.models
    masked_dataset.stack(ds0)


    return masked_dataset

#------------------------------------------------------------------------------
def createonoff_from_simulation(mc, fake=True, debug=False):
    """


    Parameters
    ----------
    mc : TYPE
        DESCRIPTION.
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    dlist_onoff : TYPE
        DESCRIPTION.

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

    dlist_onoff = Datasets()

    for ds_site in dset_list:

        for ds in ds_site:

            ds_onoff = SpectrumDatasetOnOff.from_spectrum_dataset(
                                            dataset=ds,
                                            acceptance=1,
                                            acceptance_off=1/mcf.alpha)
            if (fake): ds_onoff.fake(background_model=ds_onoff.background)
            if (debug): print(ds_onoff)

            dlist_onoff.append(ds_onoff)

    return dlist_onoff

#------------------------------------------------------------------------------
def read_from_ogip(file=None):
    """
    https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007.pdf
    """
    print(" Not implemented yet, bugs to be solved")

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
        t += ds.livetime.to(tunit).value
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
        print("MASKED")
        mask = ds.mask_safe
    else:
        mask = np.ones(ds.mask_safe.data.size,dtype=bool)

    if (show_header==True):
        print(90*"=")
        print(" id         name  livetime           on        bck",end="")
        print("       excs       pred   Mask       Flux")
        if (deeper==True):
                print("Ecenter (",e_unit,")         Ebin                                              ",
                      ds.models[0].spectral_model(1*u.TeV).unit)
        show_header=False # Since it was shown already

    if (ds.counts!= None):
        ncounts     = ds.counts.data[mask]
        ncounts_sum = ncounts.sum()
        nxs         = ds.excess.data[mask]
        nxs_sum     = nxs.sum()
    else:
        ncounts     = np.nan
        ncounts_sum = np.nan
        nxs         = np.nan
        nxs_sum     = np.nan

    nbck    = ds.background.data[mask]

    if (ds.models == None):
        import sys
        sys.exit("Please, define a Dataset model")
    npred   = ds.npred_sig().data[mask]

    print(90*"=")
    print(" {:3s}  {:>8s}    {:6.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}   {:6s} "
          .format(tag, ds.name, ds.livetime,
                  ncounts_sum, nbck.sum(), nxs_sum,
                  npred.sum(), ds.models[0].name
                  ))

    if (deeper):
        axes    = ds.background.geom.axes[0]
        print(90*"=")
        for i,_ in enumerate(ds.background.data[mask]):
            ecenter = axes.center.flatten()[i].to(e_unit)
            emin    = axes.edges.flatten()[i].to(e_unit)
            emax    = axes.edges.flatten()[i+1].to(e_unit)
            print(" {:7.1f}  [{:7.1f}, {:7.1f}]"
                  .format(ecenter.value, emin.value, emax.value),end="")
            if (ds.counts != None):
                print(" {:10.2f} {:10.2f} {:10.2f} {:10.2f}"
                      .format(ncounts.flatten()[i],
                              nbck.flatten()[i],
                              nxs.flatten()[i],
                              npred.flatten()[i]),end="")
            else:

                print(" {:>10s} {:10.2f} {:>10s} {:10.2f}"
                      .format("-",
                              nbck.flatten()[i],
                              "-",
                              npred.flatten()[i]),end="")
            print("  {:^6}  {:8.2e} "
                  .format(ds.mask_safe.data[mask].flatten()[i],
                          ds.models[0].spectral_model(ecenter).value))

    return show_header
###---------------------------------------------------------------------------------------
def check_datasets(dsets,masked=False,deeper=False):
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
    header=False
    for i,ds in enumerate(dsets):
        header=check_dataset(ds,tag=str(i),
                             show_header=header,
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




