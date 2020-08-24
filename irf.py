# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 09:01:41 2019

@author: Stolar
"""
import os
import sys
import numpy as np
from pathlib import Path

import astropy.units as u
from gammapy.cube import make_map_exposure_true_energy, make_map_background_irf
from gammapy.cube import PSFKernel

import ana_config as cf
from utilities import warning, failure, success, highlight,banner

__all__=['IRF']


###############################################################################
class IRF(object):
    """
    Class to handle time variable IRF
    
    Information from G. Maier obtained by I. Sadeh
    Preamble : the latest prod versions are prod3v2

    Former IRf can be used they essentially differ in the angular resolution at 
    high energies. Include both full array, threshold array and LST-only IRFs. 
    See here :
    https://forge.in2p3.fr/projects/cta_analysis-and-simulations/wiki/Prod3b_based_instrument_response_functions

    
    High NSB : 
        The earlier prod3b productions can be used.
        See here : https://forge.in2p3.fr/attachments/download/55554
        They cover 20 and 40 degrees and not 60 degrees.
        They can be used as an indicator for the systematics, 
        not as a precise prediction of performance.
        They have 5x or 30x the nominal NSB. 
        High NSB IRFs are only available for 20 deg zenith angle. 

    Interpolation : 
        
        ZENITH 
        Gernot recommends to define the zenith  windows by 1/cos(theta). 
        This means (Iftach says):
            zen < 33        —> use IRF 20deg
            33 <= zen < 54  —> use IRF 40deg
            zen > 54        —> use IRF 60deg
        
        AZIMUTH
        We should use the average azimuth (not the step function):
        
        ENERGY
        We should generate events from the minimum allowed value for 
        a given IRF, but then we should only use for the analysis events above 
        a slightly higher threshold, in order to avoid cutoff effects. 
        
        The IRFs are defined with these minimum energy thresholds:
            20deg —> E_gen >= 12.6 GeV
            40deg —> E_gen >= 26.4 GeV
            60deg —> E_gen >= 105.3 GeV
        Iftach suggests (does not says if generated ot reconstrcuted):
            20deg —> E_grb >= 30 GeV
            40deg —> E_grb >= 40 GeV
            60deg —> E_grb >= 110 GeV
        That corresponds roughly to [E_gen + deltaE], which is the 
        expected energy resolution at low energies (deltaE/E (68%) = 0.25)

    Note also that the CTA performance plots have a cutoff of E>20 GeV 
    since we can’t trust the sims below this range. 
    In terms of requirements, we actually only have a guarantee above 30 GeV, 
    again, due to the large uncertainties at these energies. 
    -> We can therefore not use <30GeV energies for spectra/detection in any 
    case for a consortium publication
    
    """

    ###########################################################################
    def __init__(self,
                 folder = ".",
                 kchain = "Reco1",
                 khem   = "North",
                 array  = "FullArray",
                 kzen   = "20deg",
                 kaz    = "average",
                 kdt    = "100s",
                 nsb=False):
        """
        tbw

        """

        self.theta_list =  {"20deg":20*u.deg,"40deg":40*u.deg,"60deg":60*u.deg}
        self.dt_list    = {"100s": (100*u.s),
                           "30m" : (30*u.min).to(u.s),
                           "05h" : (5*u.h).to(u.s),
                           "50h" : (50*u.h).to(u.s)}

        self.folder = folder # path
        self.kchain = kchain # key
        self.array  = array
        self.khem   = khem   # key
        self.kzen   = kzen   # key
        self.kaz    = kaz    # key
        self.kdt    = kdt    # key
        self.nsb    = nsb    # bool

        return

    ###########################################################################
    def print(self):

        print("========================================================")
        if (os.path.exists(self.folder)):
            print("==== IRF folder  : ",self.folder,end="")
        else:
            print(" NOT EXISTING !!!")
        print()
        print("========================================================")
        print("  Analysis chain : ",self.kchain)
        print("  Hemisphere     : ",self.khem)
        print("  Array          : ",self.array)
        print("  Zenith         : ",self.kzen)
        print("  Azimuth        : ",self.kaz)
        print("  Duration       : ",self.kdt)
        print("  High NSB       : ",self.nsb)

        return

    ###########################################################################
    @classmethod
    def get_irf_file(cls, theta, phi, obstime, khem="North",**kwargs):

        # This calls the default constructor
        # Don't forget to pass the location key otheriwise, North!
        irf = cls(khem=khem,**kwargs)

        # First find the best keys
        def find_best_key(mydict,x):
            dict_values    = np.asarray([x.value for k, x in mydict.items() ])
            closest_val = dict_values [  (np.abs(dict_values - x.value)).argmin() ]
            closest_key = [k for k, x in mydict.items() if x.value == closest_val]
            return closest_key[0]


        irf.kzen   = find_best_key(irf.theta_list,theta)
        irf.kaz    = "average"
        irf.kdt    = find_best_key(irf.dt_list,obstime.to(u.s))

        # print(" Best theta : ",irf.kzen," -> ",irf.theta_list[irf.kzen])
        # print(" Best phi = ",irf.kaz)
        # print(" Best obs : ",irf.kdt," -> ",irf.dt_list[irf.kdt])

        return irf.get_file_fromkeys()

###########################################################################
    def get_file_fromkeys(self,missing=0):
        """
        Find folder name and filename of an IRf file based on a series of
        keywords. In case no corresponding file is found, go back to
        a default set of keywords.
        This has to be tuned with the existing list of available IRF
        """

        if (cf.fixed_zenith != False):
            self.kzen = cf.fixed_zenith
            
        folder = Path(self.folder,self.kchain,self.khem,self.array,self.kzen)
        if folder.exists() != True:
            self.khem = "North"
            warning(" Unable to find folder {} -> Trying North"
                    .format(folder))
            folder = Path(self.folder,self.kchain,self.khem,self.array,self.kzen)
            if folder.exists() != True:
                failure(" Impossible to find a suitable IRF folder")
                sys.exit("Stopped")
            
        irf_file = []
        
        for file in os.listdir(folder):
            if (file.find(self.khem) != -1) and \
               (file.find(self.kzen) != -1) and \
               (file.find(self.kaz)  != -1) and \
               (file.find(self.kdt)  != -1):
                   irf_file.append(file)

        if (len(irf_file)>1):  sys.exit(" More than one IRF file found - ERROR")
        if (len(irf_file)==0):
            print(" Missing: ",self.khem,self.kzen,self.kaz,self.kdt,end=" ")
            if (missing != 0):
                failure(" Unable to find a suitable file")
                sys.exit("IRF file search aborted")
            self.khem = "North"
            print(" ==> ",self.khem,self.kzen,self.kaz,self.kdt)

            def_file  = self.get_file_fromkeys(missing=1)
            irf_file  = [def_file.name]
            folder    = def_file.parents[0]

        return Path(folder, irf_file[0])

    ###############################################################################
    def get_irf_components(self,obs_param, pointing, irfs, debug=False):

        """Get the IRf components"""

        axis = obs_param.axis
        duration = obs_param.livetime
        geom = obs_param.geom
        fov = obs_param.fov

        offset = 0*u.deg

        ### Compute the exposure map
        exposure = make_map_exposure_true_energy(
                    pointing = pointing.galactic, # does not accept SkyCoord altaz
                    livetime = duration,
                    aeff     = irfs["aeff"],
                    geom     = geom)

        ### Compute the background map
        background = make_map_background_irf(
                        pointing = pointing.galactic, # does not accept SkyCoord altaz
                        ontime   = duration,
                        bkg=irfs["bkg"],
                        geom=geom
                        )

        ### Point spread function
        psf = irfs["psf"].to_energy_dependent_table_psf(theta=offset)
        psf_kernel = PSFKernel.from_table_psf(psf, geom, max_radius=fov/2) # Why this max_radius ?

        ### Energy dispersion
        edisp = irfs["edisp"].to_energy_dispersion(offset,
                    e_reco=axis.edges,
                    e_true=axis.edges)

        if (debug>2):
            # Maybe this could go into a irf_plot module
            import matplotlib.pyplot as plt
            islice = 3

            exposure.slice_by_idx({"energy": islice-1}).plot(add_cbar=True)
            plt.title("Exposure map")

            background.slice_by_idx({"energy": islice}).plot(add_cbar=True)
            plt.title("Background map")

            psf_kernel.psf_kernel_map.sum_over_axes().plot(stretch="log")
            plt.title("point spread function")

            edisp.plot_matrix()
            plt.title("Energy dispersion matrix")
            plt.show()

        return exposure, background, psf_kernel, edisp


###############################################################################
if __name__ == "__main__":

    # Create an IRF object - get the default values
#    myirf = IRF(folder=cf.irf_dir,array="LST")
#    myirf.print()
#    irf_file = myirf.get_file_fromkeys()
#    print(irf_file)
#
#    # Create an IRF object - get a  values
#    myirf = IRF(kzen="40deg")
#    myirf.print()
#    irf_file = myirf.get_file_fromkeys()
#    print(irf_file)


 
    # Get file form specific values
    print("=====================================================")
    irffile = IRF.get_irf_file(folder =cf.irf_dir,
                               khem="South",
                               theta=55*u.deg, 
                               phi=123*u.deg,
                               obstime=23*u.min)
    print(irffile)

#    myirf = IRF(khem="North", kaz="average",kzen="20deg",kdt="100s")
#    myirf.print()
#    irf_file = myirf.get_file_fromkeys()
#    print(irf_file)
