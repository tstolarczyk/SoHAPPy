
Analyzing a single objects
==========================

Introduction
------------

Introduce the concept

stacking, masking and associated models
---------------------------------------


 A dataset ds can carry models, in particular a spectral model and a mask_safe.
    
    The mask_safe define the energy region on which the dataset data are considered valid.
    It is in particular used during a fit and when several datasets are stacked together.
    
    The model attached to a dataset can be used to simulate the count content. It can also be used for fitting. A model can be an analytical expression liek a powerlaw or a TemplateSpectraModel obtained from point measurement. In that case the model is an interpolation through the points.
    
    Not that fitting the data with a model different from the one used for simulation is possible.
    However the initial model should be saved before hand if it has to be reused:
    
    .. code-block:: python
    
        model_init = ds.models
        model_fit = SkyModel(spectral_model=PowerLawSpectralModel(...))
        ds.models = model_fit
        Fit(ds)
        ds.models = model_init # Go back to initial model
        
    The predcited counts in a dataset are computed using the model (and the IRF), which is a continous 
    function, at the center of the bins of the reconstrcuted axis (although the model is defined upon true energies).
    When the mask is applied, the predicted counts in the masked bins are set to zero. But the model is 
    unchanged.
    
    Stacking is usually intended to group several observations of a same phenomemon.
    Each observation can have a different masking (e.g. varying energy thresholds), but all observations are intended
    to be interpreted from the same fit model.
    
    When simulatiing a source with a varying energy spectrum in time, this is the case with GRB simulations, each time
    slice carries its own spectral model (in the most general case a template spectral model).
    
    In Gammapy the model in a stacked model dsstack obtained from two datasets ds1 and ds2 is the ds1 model 
    (as it is supposed to be a fit model for a set of observations). 
    
    However it is possible to compute an equivalent mean model as the sum of the original models weighted by the 
    dataset livetime using a template spectral model.
    
    .. code-block:: python
    
        flux1 = ds1.models[0].spectral_model(e_axis.center)
        dt1   = ds1.livetime 
        flux2 = ds2.models[0].spectral_model(e_axis.center)
        dt2   = ds2.livetime 
        
        flux_new  = (dt1*flux1 + dt2*flux2)/(dt1+dt2)
        spec      = TemplateSpectralModel(energy = e_axis.center,
                                 values = flux_org*mask_org,
                                 norm   = 1.,
                                 interp_kwargs ={"values_scale": "log"})
        dsstack.models = SkyModel(spectral_model = spec, name  = "Stack"+"-"+ds1.name)         
    
    This works perfectly as long as the two datasets are defined on the same energy range, 
    i.e. have the same mask_safe attribute.
    With this, the predicted counts, that computed from the model, can still be compared to the generated counts 
    in a simulation
    
    If the two datasets have NOT the same mask_safe attribute, one could set to zero the flux point in the definition
    of the template spectral model as follows. Note that this requires to define the model on the reconstructed 
    energy axis (axis on whih the mask _safe is defined):
    
    .. code-block:: python
    
        e_axis = ds1.background.geom.axes[0] # E reco sampling from the IRF    
        mask1 = ds1.mask_safe.data.flatten()
        flux1 = ds1.models[0].spectral_model(e_axis.center)*mask1 # masked bins at zero 
        dt1   = ds1.livetime     
    
    However the resulting interpolated spectrum can/will have edge effects that will not guarantee that the flux is 
    stricly at zero.
    
    As a consequence it is not strictly possible to propagate a mean effective model when stackig two datasets.
    
    The solution to check that the simulated count numbers follow the original modelling is to cumulate by hand the 
    prdicted counts from the original individual datasets.

#.. automodapi:: singlesource_analysis
#   :include-all-objects:
#   :no-inheritance-diagram: