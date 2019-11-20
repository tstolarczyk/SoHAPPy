# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import yaml
from yaml import CLoader as Loader
#from glob import glob

import numpy as np
import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt

from gammapy.maps import WcsGeom, MapAxis
from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel
from gammapy.image.models import SkyPointSource
#from cta_grb_observation import (GRBObservationSimulation, GRBTarget,
#                                       GRBObservationParameters)

from gammapy.spectrum import SpectrumDatasetOnOffStacker
from gammapy.data import Observations
from gammapy.stats.poisson import significance_on_off


__all__=['GammaRayBurst','GRBTarget']

###############################################################################
class GammaRayBurst(object):
    """
    Class to store GRB properties

    This class includes the calls for simulated observations.
    The GRB information is read from file.
    A GRB is composed of :
        - a name
        - a redshift
        - a position in the sky (temporarily given outside the file)
        - a list of time-intervals with correspoding :
            - energy spectrum.
            - observations
            - simulations

    Simulations are done with appropriate functions
    ~run_OnOffsimulation
    ~run_3Dsimulation
    """

    ###########################################################################
    def __init__(self, filepath      = None,
                       name          = None,
                       z             = None,
                       time_interval = None,
                       spectral_model= None,
                       pointing      = None):

        # Path of the GRB properties/model
        self.filepath = filepath

        # GRB properties
        self.name     = name
        self.z        = z
        self.fluence  = 0.
        self.pointing = pointing

        self.time_interval  = time_interval # Time intervals
        #self.t_s            = 0.5*(time_interval[:,1]-time_interval[:,0]).to(u.s).value
        self.t_s            = time_interval[:,1].to(u.s).value # End of det. interval
        self.spectral_model = spectral_model # Gammapy models (one per t slice)
        self.spatial_model  = SkyPointSource(lon_0=pointing.galactic.l,
                                             lat_0=pointing.galactic.b)

        # Stacked obs
        self.stack_obs = None

    ###########################################################################
    @classmethod
    def read(cls, filepath, absorption,reduction_factor=1, pointing = None, dbg=False):
        """Read GRB charcateristics from an ascii file.

        Read the spectrum characterisitics and create the spectrum model list.
        Initialise the GRB object content.
        Temporary :
            - The GRB position is not in the file
            - The time intervals are set by hand because of a file format that do not permit to retrieve them.
        """
        data = cls.get_grb_properties(filepath + '/grb_properties.txt')

        name          = data['name']
        z             = data['redshift']
        time_interval = data['time intervals']

        # HACK, waiting for Lara
        time_interval = [[30., 50.], [50., 80.], [80., 125.], [125., 200.],
                         [200., 315.], [315., 500.], [500., 800.], [800., 1250.],
                         [1250., 2000.], [2000., 3150.], [3150., 5000.],
                         [5000., 10000.]]
        # Attempt to associate the files with the time interval fails because the time are rounded
        #folder_file = os.listdir(filepath)[0:] # First file is grb_properties.txt
        #print(" ThS - ",folder_file)
        #interval name = []
        spectral_model = []

        #----------------------------------------------------------------------
        if (dbg > 1):
            print("======= DEBUG - Flux on Earth")
            print("        Time intervals from 'properties' file : ", len(time_interval))

            fig = plt.figure(figsize=(18.,8.))
            flux_min = +float('inf')
            flux_max = -float('inf')
            a1 = fig.add_subplot(121)
            a1.set_xscale("log")
            a1.set_yscale("log")
            a1.set_xlabel("Energy (TeV)")
            a1.set_title(filepath)

            done = False # Will get flux unit once

        #----------------------------------------------------------------------
        # Reads spectral shape from each time interval
        for t in time_interval:
            filename = '{}_{:.0f}-{:.0f}.txt'.format(name,t[0],t[1])
            table = Table.read(filepath + '/' + filename, format='ascii')

            energy = table['col1'] * u.TeV
            flux = (table['col2'] / reduction_factor) * u.Unit('1 / (cm2 s TeV)')

            table_model = TableModel(energy=energy,
                                     values=flux,
                                     norm=1.,
                                     values_scale='lin')
            spectral_model.append(
                AbsorbedSpectralModel(spectral_model=table_model,
                                      absorption=absorption,
                                      parameter=z))

            #------------------------------------------------------------------
            if (dbg > 1):
                print(" {0:6.1f} - {1:7.1f} ({2:6.1f}) - (nE = {3:3d})"
                      .format(t[0],t[1],t[1]-t[0],len(energy)))
                # Set flux unit once
                if (not done):
                   a1.set_ylabel(str(flux[0].unit))
                   done = True
                # Get min and max flux for further plots
                if flux[0].value > flux_max :  flux_max = flux[0].value
                if flux[len(flux)-1].value < flux_min: flux_min = flux[len(flux)-1].value
                # Plot flux and absorbed flux

                a1.plot(energy.value,flux.value,label=str(t[0])+"-"+str(t[1]))
                # another way to get it
                # table_model.plot([min(energy),max(energy)],a1)

        #------------------------------------------------------------------------
        if (dbg > 1):
            print("        Number of absorbed spectral model generated =",len(spectral_model))
            a1.legend()
            a2 = fig.add_subplot(122)
            a2.set_xscale("log")
            a2.set_yscale("log")
            a2.set_ylim([flux_min,flux_max])
            a2.set_title(filepath)
            for absorbed_model in spectral_model:
                absorbed_model.plot([min(energy),max(energy)],a2)
            plt.show()
            print("======= END DEBUG - Flux on Earth\n")
        #------------------------------------------------------------------------

        return cls(
            filepath=filepath,
            name=name,
            z=z,
            time_interval=time_interval * u.s,
            spectral_model=spectral_model,
            pointing=pointing)

    ###########################################################################

    def __str__(self):
        """ Printout the GRB properties """
        txt = ''
        txt += '========================== GRB summary\n'.format()
        txt += '  Name: {}\n'.format(self.name)
        txt += '  Redshift: {}\n'.format(self.z)
        txt += '  Pointing  {} :\n'.format(self.pointing)
        txt += '  Time intervals ({}):\n'.format(self.time_interval[0].unit)
        for t in self.time_interval:
            txt += '[{},{}], '.format(t[0].value, t[1].value)
        txt+= '\n'

        if self.stack_obs is not None:
            txt += str(self.stack_obs.total_stats_safe_range)

        return txt


    ###########################################################################
    def add_stack_OnOffobs(self,simulations):
        """Stack observations"""
        stack = SpectrumDatasetOnOffStacker(Observations(simulations))
        #stack.run()
        return stack.stacked_obs

    ###########################################################################
    def get_cumulative_stats(self,simulations):
        """Get cumulative statistics"""
        # Init vectors
        cumulative_time = np.zeros(len(simulations))
        cumulative_n_on = np.zeros(len(simulations))
        cumulative_n_off = np.zeros(len(simulations))
        cumulative_alpha = np.zeros(len(simulations))
        delta_t = np.zeros(len(simulations))

        # Loop on observations
        for idx, obs in enumerate(simulations):

            alpha = obs.alpha

            if idx == 0:
                cumulative_time[idx] = obs.total_stats_safe_range.livetime.value
                cumulative_n_on[idx] = obs.total_stats_safe_range.n_on
                cumulative_n_off[idx] = obs.total_stats_safe_range.n_off
                cumulative_alpha[idx] = obs.total_stats_safe_range.alpha
            else:
                cumulative_time[idx] += cumulative_time[idx-1] + obs.total_stats_safe_range.livetime.value
                cumulative_n_on[idx] += cumulative_n_on[idx-1] + obs.total_stats_safe_range.n_on
                cumulative_n_off[idx] += cumulative_n_off[idx-1] + obs.total_stats_safe_range.n_off
                cumulative_alpha[idx] = obs.total_stats_safe_range.alpha

            delta_t[idx] =  obs.total_stats_safe_range.livetime.value
            cumulative_excess = cumulative_n_on - alpha * cumulative_n_off
            cumulative_bkg = alpha * cumulative_n_off
            cumulative_significance = significance_on_off(cumulative_n_on,
                                                          cumulative_n_off,
                                                          cumulative_alpha)

        return dict(livetime=cumulative_time,
                    excess=cumulative_excess,
                    bkg=cumulative_bkg,
                    sigma=cumulative_significance,
                    n_on=cumulative_n_on,
                    n_off=cumulative_n_off,
                    alpha=cumulative_alpha,
                    delta_t=delta_t)

    ###########################################################################
    def quicklook(self,plot=False):

        for idx, obs in enumerate(self.simulations):
            print("GRB ",self.name," - Observation: ",idx)
            #obs.peek() # n_on, alpha*n_off Espec. eff.area,, Ematrix and stat
            print("    - lifetime:",obs.total_stats_safe_range.livetime.value)
            print("    - excess vector:",obs.excess_vector.data)
            if (plot):
                obs.excess_vector.peek()
            plt.show()

        return


    @staticmethod
    ###########################################################################
    def get_grb_properties(filename):
        """Get GRB properties"""
        with open(filename,'r') as stream:
            data = yaml.load(stream,Loader=Loader)
        delta_t = data['time intervals'].split()
        intervals = []
        for idx in range(0,24,2):
            # in seconds
            intervals.append([float(delta_t[idx]), float(delta_t[idx + 1])] * u.s)
        data['time intervals'] = intervals
        return data

###############################################################################
#
###############################################################################
class GRBTarget(object):
    """
    Observation target information.

    Parameters
    ----------
    name : `str`
        Name of the source
    model : `~gammapy.spectrum.models.SpectralModel`
        Model of the source
    """

    ###########################################################################
    def __init__(self, name=None, model=None,pointing = None):
        self.name = name
        self.model = model
        self.pointing = pointing
    ###########################################################################
    def __str__(self):
        """Target report (`str`)."""
        ss = '*** Target parameters ***\n'
        ss += 'Name={}\n'.format(self.name)
        for par in self.model.parameters.parameters:
            ss += '{}={} {}\n'.format(par.name, str(par.value), par.unit)
        return ss
    ###########################################################################
    def from_fermi_lat_catalogue(name):
        raise NotImplementedError

###############################################################################
#
###############################################################################
class GRBObservationParameters(object):
    """
    Container for observation parameters.

    Parameters:

    alpha : `~astropy.units.Quantity`
        Normalisation between ON and OFF regions

    livetime :  `~astropy.units.Quantity`
        Observation time

    emin :  `~astropy.units.Quantity`
        Minimal energy for simulation

    emax : `~astropy.units.Quantity`
        Maximal energy for simulation
    """

    ###########################################################################
    def __init__(self, alpha=None, livetime=None,
                 emin=None, emax=None, pointing=None, fov=None, binsize=None):
        self.alpha = alpha
        self.livetime = livetime
        self.emin = emin
        self.emax = emax
        self.fov = fov
        self.binsize = binsize

        # Create sky map frame - Emin, Emax, and the bin size could be extracted from the IRF
        nedges = 10
        logEedges = np.logspace(np.log10(self.emin.value), np.log10(self.emax.value), nedges)
        if (pointing):
            self.axis = MapAxis.from_edges(logEedges, unit="TeV", name="energy", interp="log")
            self.geom = WcsGeom.create(skydir   = pointing,
                                       binsz    = binsize,
                                       width    = (fov, fov),
                                       coordsys ="GAL",
                                       axes=[self.axis]) # Width in degrees

#        print(" - Number of Energy log. bins : ",nedges-1," from ",self.emin," to ",self.emax)
#        print(" - Number of spatial bins : ",self.fov/self.binsize," x ",self.fov/self.binsize)
#        print(" *** Sky map")

    ###########################################################################
    def __str__(self):
        """Observation summary report (`str`)."""
        ss = '*** Observation parameters summary ***\n'
        ss += 'alpha={} [{}]\n'.format(self.alpha.value, self.alpha.unit)
        ss += 'livetime={} [{}]\n'.format(self.livetime.value,
                                          self.livetime.unit)
        ss += 'emin={} [{}]\n'.format(self.emin.value, self.emin.unit)
        ss += 'emax={} [{}]\n'.format(self.emax.value, self.emax.unit)
        return ss






