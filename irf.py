# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 09:01:41 2019

@author: Stolar
"""
import os
import sys
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt

import astropy.units as u
from gammapy.cube import make_map_exposure_true_energy, make_map_background_irf
from gammapy.cube import PSFKernel

__all__=['IRF']


###############################################################################
class IRF(object):
    """
    Class to handle time variable IRF

    """

    ###########################################################################
    def __init__(self,
                 kchain="Reco1",folder="../input/irf",
                 khem="North",kzen="20deg",kaz="average",kdt="100s",
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
        self.khem   = khem   # key
        self.kzen   = kzen   # key
        self.kaz    = kaz    # key
        self.kdt    = kdt    # key
        self.nsb    = nsb    # bool

        return

    ###########################################################################
    def print(self):

        print("========================================================")
        print("==== IRF folder  : ",self.folder)
        print("========================================================")
        print("  Analysis chain : ",self.kchain)
        print("  Hemisphere     : ",self.khem)
        print("  Zenith         : ",self.kzen)
        print("  Azimuth        : ",self.kaz)
        print("  Duration       : ",self.kdt)
        print("  High NSB       : ",self.nsb)

        return

    ###########################################################################
    @classmethod
    def get_irf_file(cls, theta, phi, obstime, khem="North"):

        # This calls the default constructor
        # Don't firget to pass the location key otheriwise, North!
        irf = cls(khem=khem)

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
    def get_file_fromkeys(self):
        """
        Find folder name and filename of an IRf file based on a seried of
        keywords. In case no corresponding file is found, go back to
        a default set of keywords.
        This has to be tuned with the existing list of available IRF
        """

        folder = Path(self.folder,self.kchain,self.khem)
        irf_file = []
        for file in os.listdir(folder):
            if (file.find(self.khem) != -1) and \
               (file.find(self.kzen) != -1) and \
               (file.find(self.kaz)   != -1) and \
               (file.find(self.kdt)  != -1):
                   irf_file.append(file)

        if (len(irf_file)>1):  sys.exit(" More than one IRF file found - ERROR")
        if (len(irf_file)==0):
            print(" No IRF found for keys : ",
                  self.khem,self.kzen,self.kaz,self.kdt,end=" ")
            self.khem = "North"
            #self.kzen = "20deg"
            #self.kaz  = "average"
            #self.kdt  = "100s"
            print("   ==> :",
                  self.khem,self.kzen,self.kaz,self.kdt)

            def_file  = self.get_file_fromkeys()
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
            islice = 3

            exposure.slice_by_idx({"energy": islice-1}).plot(add_cbar=True);
            plt.title("Exposure map")

            background.slice_by_idx({"energy": islice}).plot(add_cbar=True);
            plt.title("Background map")

            psf_kernel.psf_kernel_map.sum_over_axes().plot(stretch="log");
            plt.title("point spread function")

            edisp.plot_matrix();
            plt.title("Energy dispersion matrix")
            plt.show()

        return exposure, background, psf_kernel, edisp


###############################################################################
if __name__ == "__main__":

    # Create an IRF object - get the default values
#    myirf = IRF()
#    myirf.print()
#    irf_file = myirf.get_file_fromkeys()
#    print(irf_file)
#
#    # Create an IRF object - get a  values
#    myirf = IRF(kzen="40deg")
#    myirf.print()
#    irf_file = myirf.get_file_fromkeys()
#    print(irf_file)


    print("............................................\n")

    # Get file form specific values
    print("=====================================================")
    irf_file = IRF.get_irf_file(theta=45*u.deg, phi=123*u.deg,obstime=23*u.min)
    print(irf_file)

#    myirf = IRF(khem="North", kaz="average",kzen="20deg",kdt="100s")
#    myirf.print()
#    irf_file = myirf.get_file_fromkeys()
#    print(irf_file)
