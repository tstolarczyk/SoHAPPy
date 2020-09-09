# -*- coding: utf-8 -*-

"""
Created on Thu Dec 12 09:01:41 2019

@author: Stolar

This module hosts the functionnalities to choose the best IRF from the available
files and the recipes on their validity domains in observation duration and
angles.

The recommandations for prod3v2 IRF are the following (Information from
G. Maier obtained by I. Sadeh on January 30th, 2020).

Angular "interpolation"
----------------------

* Zenith
The zenith windows should be defined following a 1/cos(theta) law, which means
(Iftach says) that depending on the zenith angle the following IRF should be
used :
    - zen < 33        : "20deg"
    - 33 <= zen < 54  : "40deg"
    - zen > 54        : "60deg"

On top of this, it is foreseen that the instrument performance will not be
reliable above 66 degrees as implicitly mentionned in A-PERF-0720
(*Observable Sky: The system as a whole must be able to target any astrophysical
object in the sky which has an elevation >24 degrees*).

* Azimuth:
It is recommended to use the average azimuth (not a step function going when
using the North and South IRF).

Energy thresholds
-----------------
Although events should be generated from the minimum allowed value for a given
IRF, the analysis should restrict to a higher threshold, in order to avoid
cutoff effects.

Minimum *generated* energies are the following:
    - 20deg —> E_gen >= 12.6 GeV
    - 40deg —> E_gen >= 26.4 GeV
    - 60deg —> E_gen >= 105.3 GeV

Minimum *reconstructed* energies are the following (Iftach suggests these
values but does not say it this is his conclusion or more generally accepted
values):
    - 20deg —> E_grb >= 30 GeV
    - 40deg —> E_grb >= 40 GeV
    - 60deg —> E_grb >= 110 GeV
This corresponds approximately to the generated energies to which one standard
deviation has been added (the expected energy resolution at low energies is
deltaE/E (68%) = 0.25.

Note that the CTA performance plots have a cutoff of E>20 GeV.
Requirements expect the IRF to be reliable above 30 GeV.


older IRF files
---------------
IRF older than prodv3b2 can be used. They essentially differ in the angular
resolution at high energies. They include both full arrays, threshold arrays
and LST-only IRFs. They cover 20 and 40 degrees but not 60 degrees.
Some 20 degree IRF have 5x 5North and South) or 30x the nominal NSB.

See here :
https://forge.in2p3.fr/projects/cta_analysis-and-simulations/wiki/Prod3b_based_instrument_response_functions
And download there : https://forge.in2p3.fr/attachments/download/55554

"""
import os
import sys
import numpy as np
from pathlib import Path

import astropy.units as u

# # Avoid deprecation Astropy warnings in gammapy.maps
import warnings
with warnings.catch_warnings():
    from gammapy.cube import make_map_exposure_true_energy, make_map_background_irf
    from gammapy.cube import PSFKernel

from utilities import warning, failure

# Keyword lists and true values, also in log for time
zenith_list =  {"20deg": 20*u.deg,
                "40deg": 40*u.deg,
                "60deg": 60*u.deg}

dt_list    =  {"100s" : (100*u.s),
               "30m"  : (30*u.min).to(u.s),
               "05h"  : (5*u.h).to(u.s),
               "50h"  : (50*u.h).to(u.s)}

dtl    =  {"100s" : np.log10((100*u.s).value),
           "30m"  : np.log10((30*u.min).to(u.s).value),
           "05h"  : np.log10((5*u.h).to(u.s).value),
           "50h"  : np.log10((50*u.h).to(u.s).value)
                   }

# Validity range of IRF in zenith (see documentation of this module).
# The 60deg IRF is allowed to be used down to alitude zero for tests
# Its use if ofreseen tobe limitee bythe altmin variable
theta_valid = {"20deg": [0*u.deg, 33*u.deg],
               "40deg": [33*u.deg, 54*u.deg],
               "60deg": [54*u.deg, 90*u.deg] # Should be limited by altmin
               }

# Validity range of IRF in time, taking into account that the validity
# intervals are somehow in log. scale
# The edge values are the following "
# 0, 424.26 s, 5692.01 s (94.9'), 56921.0 s (15.8h)

dt_log_valid = {"100s": [0,
                         10**(0.5*( dtl["100s"] + dtl["30m"] )) ],
                "30m" : [10**(0.5*( dtl["100s"] + dtl["30m"] )),
                         10**(0.5*( dtl["30m"]  + dtl["05h"] )) ],
                "5h"  : [10**(0.5*( dtl["30m"]  + dtl["05h"] )),
                         10**(0.5*( dtl["05h"]  + dtl["50h"] )) ],
                "50h" : [10**(0.5*( dtl["05h"]  + dtl["50h"] )),
                         np.Inf]
                }

# Minimal acceptable energies depending on the IRF
# (resp. generated and reconstructed energies)
egen_min = {"20deg": 12.6*u.GeV,
            "40deg": 26.4*u.GeV,
            "60deg": 105.3*u.GeV
            }
erec_min = {"20deg": 40*u.GeV,
            "40deg": 26.4*u.GeV,
            "60deg": 110*u.GeV
            }


__all__=['IRF']

###############################################################################
class IRF(object):
    """
        Class to handle time variable IRF.

    """

    ###########################################################################
    def __init__(self,
                 folder = "../input/irf/OnAxis/",
                 kchain = "Reco1",
                 khem   = "North",
                 array  = "FullArray",
                 kzen   = "20deg",
                 kaz    = "average",
                 kdt    = "100s",
                 nsb    = False):
        """
        Define default values for getting IRF data from an IRF file folder
        and keywords.

        """
        self.folder = folder # path
        self.name   = "undef"# irf filename
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
        """
        Print out the IRF class content

        """

        print("========================================================")
        if (os.path.exists(self.folder)):
            print("==== IRF folder  : ",self.folder,end="")
        else:
            print(" NOT EXISTING !!!")
        print()
        print("========================================================")
        print("  File name      : ",self.name)
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
    def get_irf_file(cls, zenith, azimuth, obstime, loc = None,
                     closest = False):
        """
        From the current observation features, zenith and azimuth angles
        and observation time, find the best IRF available among the available
        list of samples angles and duration. Define the keys allowing to
        retrieve the corresponding IRF file on disk.

        Parameters
        ----------
        cls : IRF class
            Current IRF instance
        theta :
            DESCRIPTION.
        phi : TYPE
            DESCRIPTION.
        obstime : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        self=cls()

        if (loc == None):
            sys.exit("Get_Irf_file : location should be defined")

        ###--------------------------
        def find_best_key(mydict,x):
            dict_values    = np.asarray([x.value for k, x in mydict.items() ])
            closest_val = dict_values [  (np.abs(dict_values - x.value)).argmin() ]
            closest_key = [k for k, x in mydict.items() if x.value == closest_val]
            return closest_key[0]
        ###--------------------------

        # Zenith
        if (closest):
            self.kzen = find_best_key(zenith_list, obstime)
        else:
            found = False
            for k,v in theta_valid.items():
                if (zenith >=v[0] and zenith < v[1] ):
                    self.kzen = k
                    found = True
                    continue
            if (not found):
                sys.exit("get_irf_file: zenith= {} => range not found"
                         .format(zenith))

        # Azimuth - implement here N, S choice
        self.kaz = "average"

        # Observation time
        if (closest):
            self.kdt =find_best_key(dt_list, obstime.to(u.s))
        else:
            found = False
            for k,v in dt_log_valid.items():
                if obstime.to(u.s).value >= v[0] \
                   and obstime.to(u.s).value < v[1] :
                    self.kdt = k
                    found = True
                    continue
            if (not found):
                sys.exit("get_irf_file: obstime range not found")

        return self.get_file_fromkeys()

    ###########################################################################
    def get_file_fromkeys(self,fixed_zenith = None, missing=0):
        """
        Find folder name and filename of an IRf file based on a series of
        keywords. In case no corresponding file is found, go back to
        a default set of keywords.
        This has to be tuned with the existing list of available IRF.
        """

        if (fixed_zenith != None):
            self.kzen = fixed_zenith

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
            self.irf_file  = [def_file.name]
            self.folder    = def_file.parents[0]

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
    """
    Code example to use the IRF class
    """

    import ana_config as cf

    # Get file form specific values
    print("=====================================================")
    irf = IRF()
    irffile = irf.get_irf_file(loc="South",
                               zenith  =55*u.deg,
                               azimuth =123*u.deg,
                               obstime =95*u.h)
    irf.print()
    print(irffile)
    print(dt_log_valid)
    # Create an IRF object - get the default values and the corresponding file
    case1 = False
    if (case1):
        myirf = IRF(folder=cf.irf_dir,array="LST")
        myirf.print()
        irf_file = myirf.get_file_fromkeys()
        print(irf_file)

    case2 = False
    if (case2):
        # Create an IRF object - get a  values
        myirf = IRF(kzen="40deg")
        myirf.print()
        irf_file = myirf.get_file_fromkeys()
        print(irf_file)