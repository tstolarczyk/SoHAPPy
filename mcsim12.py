# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:22:11 2019

@author: Stolar
"""

import sys
import numpy as np
import os
import time
import pickle

from   astropy.coordinates import SkyCoord

from regions import CircleSkyRegion
from gammapy.utils.random import get_random_state

import gammapy
if (gammapy.__version__=="0.12"):
    from gammapy.spectrum import SpectrumDatasetOnOffStacker
    from gammapy.data import Observations

from fit_onoff12 import mc_onoff, cumulative_OnOffstats

import mcsim_config as mcf

from utilities import success

__all__ = ['MonteCarlo']
###############################################################################
class MonteCarlo():
    """
    This class handles a GRB simulation, either on one of the two sites, and
    in a next version combining both sites.
    It uses several additionnal modules :
        - mcsim_config : contains variables and flags with some useful options
        that are not intended to be changed by the user (suggested values are
        mentionned between brackets):
            - det_level (0.9):
                Fraction declaring a detection above the 3,5 sigma threshold
            - alpha (0.2) : related to the number of off N zones in the
            on-off analysis, alpha = 1/N
            - fov (5.*u.deg) : Full Field-of-view for the sky map
            - pxsize (0.125*u.deg) : square pixel length for the sky map
            - dtslew (30*u.s) : slewing time (supposed identival for all
            telescope types)
            - fixslew (True) : if True teh slewing time if fixed to dtslew,
            or randomly generated between 0 and dtslew if False.
            - n_lima_min = 10 # Count number below which Li&Ma cannot be trusted anymore
write_slices = True # Store detailed information on slices if True


from   astropy.coordinates   import EarthLocation
# The following information can be obtained from EarthLocation.of_site
# but this requires to have an internet connection
# Get all available site names :
# sitelist = EarthLocation.get_site_names()

# Get coordinattes of Paranal and LaPalma

# xyz_north = EarthLocation.of_site('Roque de los Muchachos')
xyz_north = EarthLocation.from_geocentric( 5327448.9957829,
                                          -1718665.73869569,
                                           3051566.90295403,
                                           unit="m")

# xyz_south = EarthLocation.of_site('Paranal Observatory')
xyz_south = EarthLocation.from_geocentric( 1946404.34103884,
                                          -5467644.29079852,
                                          -2642728.20144425,
                                          unit="m")

pos_site ={"North":xyz_north, "South":xyz_south}
        - mcsim_plot : generate plot from the simulation results.
    """
    ###------------------------------------------------------------------------
    def __init__(self,
                 niter  =1,
                 method   = 0,
                 fluctuate = True,
                 seed= 'random-seed',
                 debug  = 0,
                 name="Unknown"):
        """
        Initialize class members to default values

        """
        self.dbg          = debug

        # Input parameters and objects
        self.niter        = niter     # Number of trials
        self.method       = method   # Deg. of freedom in likelikood, 0 if on-off
        self.fluctuate    = fluctuate
        self.on_region    = None
        self.slot         = None   # This is initialized at run time
        self.name         = name
        self.rnd_state    = get_random_state(seed)
        # Initialise model associated to the MC
        #self.source_model = self.source_spectrum()
        #self.null_model   = self.source_spectrum(norm=0/(u.cm)**2/u.s/u.TeV)

        # list of simulations (one per slice)
        self.simulations = None

        # Significance over simulations
        self.id_smax_list  = [] # Slice number to get back the time/ altaz
        self.smax_list     = [] # List of max. significances along trials
        self.nex_smax_list = [] # List of excess counts at max. signif.
        self.nb_smax_list  = [] # List of background counts at max. signif.

        self.id_3s_list    = [] # Slice number ot get back the time/altaz
        self.nex_3s_list   = [] # List of excess
        self.nb_3s_list    = [] # List of background
        self.detect_3s     = 0  # Number of trials 3 sigma was reached

        self.id_5s_list    = [] # Slice number to get back the time/altaz
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
    def run(self,slot,boost=True,savedset="notused",dump_dir=None):
        """
        - Run simulations of the current grb, for a given hemisphere.
          A simulation correspond to a series of observations corresponding
          to the GRB time slices.

        - Compute and store :
            - Significance values for all slices in the trials;
            - Maximum significance reached for the GRB;
            - Time of maximum siginificance, number of excess and background
            events
            - Time to reach 3 or 5 sigma, number of excess and background
            events

        """

        self.slot = slot

        on_pointing = SkyCoord(self.slot.grb.radec.ra + mcf.offset,
                                    self.slot.grb.radec.dec,
                                    frame="icrs")

        self.mc_on_region = CircleSkyRegion(center = on_pointing,
                                            radius = mcf.on_size)
        self.mctime = time.time() # Starts chronometer
        self.err    = self.niter # GRB is visible, simul. ought to be complete

        if (boost):
            abort_test  = False # Abortion test not yet performed
        else:
            abort_test  = True   # Abortion test supposed already performed

        dump = False
        if (dump_dir != None): # File for problem tracking - deleted if none
            sdump = dump_dir +"/"+ self.slot.grb.name+"-"+self.slot.site+"_slices.txt"
            fslice  = open(sdump,"w")
            print("   ---",sdump," opened")
            dump = True
            to_be_dumped = False
            self.dump_slices(header=True, file=fslice)

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

            # Loop on GRB time slices - get the signal
            ###
            # Gammpay 0.12 - JLF implementation
            ###
            obs_stat = self.loop_slices_onoff(slot)
            sigma    = obs_stat['sigma']
            non      = obs_stat['n_on']
            noff     = obs_stat['n_off']
            nex      = non - mcf.alpha*noff
            nb       = mcf.alpha*noff

            if (dump): # dump slices to track problems
                status = self.dump_slices(iMC=iMC,
                                          non=non,noff=noff,sigma=sigma,
                                          file=fslice)
                if (to_be_dumped == False): to_be_dumped = status

            if (self.dbg > 2):
                from mcsim_plot import onetrial12
                from fit_onoff12 import cumulative_OnOffstats

                onetrial12(cumulative_OnOffstats(self.simulations))

            # Acummulate sum and sum**2 for mean / error in each slice
            sigma_sum  += sigma
            sigma2_sum += sigma**2

            self.fill_stat(sigma, nex, nb) # Update stat list

            # In case 3 signma is not reached in the first 10% of the trials
            # then the 90% CL can not be reached.
            if (abort_test==False):
                if (iMC/self.niter > 1 - mcf.det_level):
                    abort_test = True
                    if (self.detect_3s ==0):
                        # print("\n *** Will never reach ",
                        #       100*mcf.det_level,"% of CL",end="")
                        # print(" ===> Simulation stopped")
                        self.err = iMC
                        break

            iMC+=1 # End of MC loop

        if (dump) :
            print()
            fslice.close()
            print("   ---",sdump," closed")
            if (not to_be_dumped):
                print("   ---",sdump," deleted")
                os.remove(sdump)

        self.mctime = time.time() - self.mctime

        ### Mean values
        self.mctime /= self.niter
        self.sigma_mean = sigma_sum/self.niter
        sigma2_mean     = sigma2_sum/self.niter
        self.sigma_std  = np.sqrt(sigma2_mean-self.sigma_mean**2)

        return

    ###------------------------------------------------------------------------
    def fill_stat(self,sigma, nex, nb):
        """
        Get the statistics and handle the exceptions of the current MC
        simulation
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
            print(" >>> t_smax = {:5.2f}, ({:5.2f} sigma at slice {:2d})"
                  .format(self.slot.slices[maxidx].tobs(),sigma[maxidx],maxidx+1))

        return

    ###------------------------------------------------------------------------
    def loop_slices_onoff(self,slot):
        """
        Get counts for each time slice.
        Stack observations.

        Loop over time interval and create a target, a set of observation parameters.
        Make the Monte Carlo simulations

        """


        self.simulations = []
        for s in slot.slices:
            #print("*===> Simulating interval ",s.idt())
            simu = mc_onoff(s, mc=self, alpha = mcf.alpha)
            self.simulations.append(simu)
            self.stack_obs = SpectrumDatasetOnOffStacker(Observations(simu))
            #stack.run()
        obs_stat = cumulative_OnOffstats(self.simulations,
                                      n_min = mcf.nLiMamin,
                                      alpha = mcf.alpha)
        return obs_stat

    ###------------------------------------------------------------------------
    def dump_slices(self,iMC=0,non=0,noff=0,sigma=0,file=None,header=False):
        """
        Print ou the list of on and off counts and the corresponding
        Li & Ma value. If non and off are less than a minimul value
        return a True status flag for further actions
        """
        status = False

        if (header):
            print("     ",file=file)
            for s in self.slot.slices:
                print(s)
                print("{:8.1f}".format(s.ts2()-s.ts1()), end="    ",file = file)
            print(file=file)

        else:
            print("{:4d}".format(iMC),end="",file = file)
            for i in range(len(non)):
                if(non[i]<=mcf.n_lima_min):
                    sep="*   "
                    status = True
                    self.err_slice += [i]
                else:
                    sep="    "
                print("{:8.1f}".format(non[i]), end=sep,file = file)
            print(file = file)

            print("{:4d}".format(iMC),end="",file = file)
            for i in range(len(noff)):
                if(noff[i]<=mcf.n_lima_min):
                    sep="*   "
                    status = True
                    self.err_slice += [i]
                else:
                    sep="    "
                print("{:8.1f}".format(noff[i]), end=sep,file = file)
            print(file = file)

            print("{:4d}".format(iMC),end="",file = file)
            for i in range(len(sigma)):
                if (sigma[i]>=3):
                    sep = "*** "
                else:
                    sep= "    "
                print("{:8.1f}".format(sigma[i]), end=sep,file = file)
            print(file = file)

        return status

    ###------------------------------------------------------------------------
    def save_to_disk(self,filename="None"):

        if (filename == "None"):
            sys.exit("Output file not defined)")

        outfile  = open(filename,"wb")
        pickle.dump(self,outfile)
        outfile.close()
        success(" Saving simulation to file : {}".format(filename))
        return