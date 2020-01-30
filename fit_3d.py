# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
#import numpy as np
import matplotlib.pyplot as plt
from gammapy.cube import MapDataset, MapEvaluator
from gammapy.utils.fitting import Fit
from gammapy.cube.models import BackgroundModel

# Suppress RunTimeWarning in mc_3d background fitting
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
###############################################################################
def fit_3d(mc, counts_map, irf, debug=0):

    """
    Simulate 3D likelihood observation

    the :meth:`Fit.run` method minimizes the likelihood() :math:`L` function of the MapDataset.
    By default the Cash statistic is used, and the returned value is :math:`FCN = -2*log(L)`.
    The Cash statistics simply computes the Poisson probability that the data map ressemble the model map at each pixel, and compute the likelihood as the product of all the pixel probabilities.

    The likelihood is computed two times, first with a source and background model, :math:`L_{S+B}`, and then with a background only, :math:`L_B` (see Mattox et al. for details).
    The test statistics

    :math:`TS=FCN_{B}-FCN_{S+B}`

    behaves as :math:`chi^2` with a number of degree of freedom equal to the differences between the two hypothesis  (See Wilks).
    In the case where only the source amplitude is the free parameter, the number of degree of freedom is 1. If a spectral index is fitted, the number of degree of freedom is 2. Adding the fit of the source position would increase this number to 4.
    Because the amplitude (or any other parameter) is restricted to positive (physical) values, the fit process explores only half of the phase space. As a consequence the TS correspond numerically to :math:`chi^2/2`.
    The :math:`p_{value}` can be obtained from the :math:`chi^2/2` probability distribution.

    References :
        - Mattox et al., Astrophys. J. 461 (1996) 396-407, for the TS defintion
        - S.S.Wilks (1938), 'The large sample distribution of the likelihood ratio for testing composite hypotheses', in the Annals of Mathematical Statistics, http://www.jstor.org

    """
    expo, bckg, psf, edisp = irf


    ### Background to be fitted
    bckg_model = BackgroundModel(bckg)
    bckg_model.parameters["norm"].value  = 1.0
    bckg_model.parameters["norm"].frozen = False
    bckg_model.parameters["tilt"].frozen = True

    ### Fit signal + background
    fitdataset = MapDataset(model     = mc.source_model,
                            exposure  = expo,
                            counts    = counts_map,
                            background_model = bckg_model,
                            psf       = psf,
                            edisp     = edisp,)
    fit    = Fit(fitdataset)
    result = fit.run(optimize_opts={"print_level": 0})
    # print(result)

    if (debug >1): DataSetShow(fitdataset,label = "S+B fit ", plot=False)

    # Results
    fcn_fit = result.total_stat
    ns_fit  = MapEvaluator(fitdataset.model,expo,psf,edisp)
    ns_fit  = ns_fit.compute_npred().data.sum()
    nb_fit  = fitdataset.background_model.evaluate()
    nb_fit  = nb_fit.data.sum()

    ### Fit signal + background
    bckset = MapDataset(model    = mc.null_model,
                        exposure = expo,
                        counts   = counts_map,
                        background_model = bckg_model,
                        psf      = psf,
                        edisp    = edisp,)
    fit_bckg = Fit(bckset) # Putting dataset here lead to an error
    result_bckg = fit_bckg.run(optimize_opts={"print_level": 0})

    if (debug>1): DataSetShow(bckset,label = "B fit   ", plot=False)

    # results
    fcn_bonly = result_bckg.total_stat
    ns_bonly  = MapEvaluator(bckset.model,expo,psf,edisp)
    ns_bonly  = ns_bonly.compute_npred().data.sum()
    nb_bonly  = bckset.background_model.evaluate()
    nb_bonly  = nb_bonly.data.sum()

    if (debug>1): print_stat(fcn_fit,ns_fit,nb_fit,fcn_bonly,ns_bonly,nb_bonly)

    # Returns the simulation results
    return fcn_fit, fcn_bonly, ns_fit, nb_fit, ns_bonly, nb_bonly

###############################################################################
def print_stat(fcn_fit,ns_fit,nb_fit,fcn_bonly,ns_bonly,nb_bonly):

    print(" FIT RESULTS (simulate_3d)")
    print(" ----------------------------")
    print("* S+B:  FCN = {:5.2f} (S={:5.2f} B={:5.2f})"
          .format(fcn_fit,ns_fit,nb_fit))
    #print("       amplitude -> {:5.2e}"
    #      .format(result.parameters["amplitude"].value))
    print("* B:    FCN = {:5.2f} (S={:5.2f} B={:5.2f})"
          .format(fcn_bonly,ns_bonly,nb_bonly))
    print(" ----------------------------\n")

    return

###############################################################################
def DataSetShow(dataset, label="", debug=0, plot=True):
    """Plot data set maps and print statistics"""
    # plot = True
    map_model = dataset.npred()
    if (map_model != None):
        # print(" dbg : dataset contains a model")
        map_sig     = MapEvaluator(dataset.model,
                                   dataset.exposure,dataset.psf,
                                   dataset.edisp).compute_npred()
        map_bkg     = dataset.background_model.evaluate()
        mtot_counts = map_model.data.sum()
        msig_counts = map_sig.data.sum() # To be checked
        mbkg_counts = map_bkg.data.sum()

        map_counts = dataset.counts
        if (map_counts != None):
            ngen = map_counts.data.sum()
        else:
            ngen = 0

        if (debug):
            print("* ",label, "Tot={:10.2f} (gen={:10.2f}), B={:5.2f} S={:5.2f}"
              .format(mtot_counts,ngen, mbkg_counts,msig_counts ))

#    print(dataset.model)
#    print(dataset.background_model)

#    map_initial = dataset.counts  # given to fit
    map_final   = dataset.npred() # after fit from the model evaluation
#    fit = (map_initial != None)
#
#    # Count map
    count_bckg    = dataset.background_model.evaluate()
    count_src     = MapEvaluator(dataset.model,
                                 dataset.exposure,dataset.psf,
                                 dataset.edisp).compute_npred()
    count_final   = map_final.sum_over_axes()
#    if (fit):
#        count_initial = map_initial.sum_over_axes()
#        count_res     = count_initial - count_final
#
#    print("\n DATASET statistics : ",label)
#    print("   Initial counts :")
#    print("         - Total         = {:5.2f}".format(count_final.data.sum()))
#    print("                   rate -> {:5.2f}".format(count_final.data.sum()/duration))
#    print("         - Background    = {:5.2f}".format(count_bckg.data.sum()) )
#    print("         - Signal        = {:5.2f}".format(count_src.data.sum()) )
#    if (fit):
#        print("   Generated counts :")
#        print("     - Generated     = {:5.2f}".format(count_initial.data.sum()))
#        print("     - Difference    = {:5.2f}".format(count_initial.data.sum()
#                                                     -count_src.data.sum()
#                                                     -count_bckg.data.sum()))
#        print("     - Residual      = {:5.2f}",count_res.data.sum())
#
#    if (plot and fit):
#        fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(10,4))
#
#        count_initial.plot(ax= ax[0],cmap='viridis',add_cbar=True)
#        ax[0].set_title("Initial count map")
#
#        count_res.plot(ax=ax[1],cmap='viridis',add_cbar=True)
#        ax[1].set_title("Residuals")
#        plt.show()
#
    if (plot):
        fig, ax = plt.subplots(ncols=3,nrows=1,figsize=(14,4))

        c0 = ax[0].imshow(count_src.data.sum(axis=0),cmap='viridis')
        ax[0].set_title("Signal model (sum over E)")
        plt.colorbar(c0 , ax = ax[0])

        c1 = ax[1].imshow(count_bckg.data.sum(axis=0),cmap='viridis')
        ax[1].set_title("Bckg model (sum over E)")
        plt.colorbar(c1,ax=ax[1])

        count_final.plot(ax=ax[2],cmap='viridis',add_cbar=True)
        ax[2].set_title("Resulting count map")

        plt.show()
        print(" ------------------------------------------------------------------")

    return