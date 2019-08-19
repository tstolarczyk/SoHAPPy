import yaml
from glob import glob
import numpy as np
import astropy.units as u
from astropy.table import Table
from gammapy.spectrum.models import TableModel, AbsorbedSpectralModel
#from gammapy.scripts.cta_utils import (CTAObservationSimulation, Target,
#                                       ObservationParameters)
from cta_utils_imported import (CTAObservationSimulation, Target,
                                       ObservationParameters)
from gammapy.spectrum import (SpectrumObservationList,
                              SpectrumObservationStacker)
from gammapy.stats.poisson import (significance_on_off,
                                   excess_error,
                                   background_error)

#### ThS
import matplotlib.pyplot as plt
import os

#------------------------------------------------------------------------------
class GammaRayBurst(object):
    """
    Class to store GRB properties.

    Simulations are done with appropriate functions
    """

    #--------------------------------------------------------------------------
    def __init__(self, filepath=None,
                 name=None,
                 z=None,
                 time_interval=None,
                 spectral_model=None):

        # Path of the GRB properties/model
        self.filepath=filepath

        # GRB properties
        self.name=name
        self.z=z

        # Time intervals
        self.time_interval=time_interval
        # Gammapy models
        self.spectral_model=spectral_model

        # list of simulations
        self.simulations = None
        # Stacked obs
        self.stack_obs = None
        
   #--------------------------------------------------------------------------
    @classmethod
    def from_file(cls, filepath, absorption,dbg=False):
        """Read from an ascii file."""
        data = cls.get_grb_properties(filepath + '/grb_properties.txt')
        
        name = data['name']
        z = data['redshift']
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
            #for ti in time_interval:
            #print(ti.value)
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
        for interval in time_interval:
            filename = '{}_{:.0f}-{:.0f}.txt'.format(name,
                                                     interval[0],
                                                     interval[1])
            table = Table.read(filepath + '/' + filename, format='ascii')

            energy = table['col1'] * u.TeV
            flux = table['col2'] * u.Unit('1 / (cm2 s TeV)')
            
            
                                
            # Deprecated 
            # table_model = TableModel(energy=energy,
            #                                     values=flux,
            #                                     scale=1.,
            #                                     scale_logy=False)
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
                print("       Interval :",interval[0],interval[1]," - E points = ",len(energy))
                # Set flux unit once
                if (not done):
                   a1.set_ylabel(str(flux[0].unit))
                   done = True
                # Get min and max flux for further plots
                if flux[0].value > flux_max :  flux_max = flux[0].value
                if flux[len(flux)-1].value < flux_min: flux_min = flux[len(flux)-1].value
                # Plot flux and absorbed flux
                
                a1.plot(energy.value,flux.value,label=str(interval[0])+"-"+str(interval[1]))
                # another way to get it 
                # table_model.plot([min(energy),max(energy)],a1)
            #------------------------------------------------------------------
            
            
                  
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
            spectral_model=spectral_model
        )

    #--------------------------------------------------------------------------
    def __str__(self):
        txt = ''
        txt += '========================== GRB summary\n'.format()
        txt += '  Name: {}\n'.format(self.name)
        txt += '  Redshift: {}\n'.format(self.z)
        txt += '  Time intervals:\n'
        for t in self.time_interval:
            txt += '     {} -- {}\n'.format(t[0], t[1])
        txt+= '\n'

        if self.stack_obs is not None:
            txt += str(self.stack_obs.total_stats_safe_range)
            
        return txt
    

#### ---------------------------------------------------------------------
# Note from gammapy.scripts.cta_utils import (CTAObservationSimulation, Target, ObservationParameters)
    def run_simulation(self, perf,
                       emin=0.03 * u.TeV,
                       alpha=0.2,
                       random_state='random-seed'):
        
        self.simulations = []
        for idx, interval in enumerate(self.time_interval):
            # print(" == Simulating interval ",idx)
            target = Target(name=self.name, model=self.spectral_model[idx])
            
            livetime = interval[1] - interval[0]
            obs_param = ObservationParameters(
                alpha=alpha * u.Unit(''), livetime=livetime,
                emin=emin, emax=1000 * u.TeV
            )

            simu =  CTAObservationSimulation.simulate_obs(
                perf=perf,
                target=target,
                obs_param=obs_param,
                random_state=random_state
            )

            self.simulations.append(simu)
            self.stack_obs = self.add_stack_obs()
            
    def add_stack_obs(self):
        """
        Stack observations
        """
        stack = SpectrumObservationStacker(
            SpectrumObservationList(self.simulations)
        ) 
        stack.run()
        return stack.stacked_obs

    def get_cumulative_stats(self):
        """
        Get cumulative statistics
        """
        # Init vectors
        cumulative_time = np.zeros(len(self.simulations))
        cumulative_n_on = np.zeros(len(self.simulations))
        cumulative_n_off = np.zeros(len(self.simulations))
        cumulative_alpha = np.zeros(len(self.simulations))
        delta_t = np.zeros(len(self.simulations))
        
        # Loop on observations
        for idx, obs in enumerate(self.simulations):

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
        
    def quicklook(self,plot=False):
        
        for idx, obs in enumerate(self.simulations):
            print("GRB ",self.name," - Observation: ",idx)
            #obs.peek() # Plots n_on, alpha*n_off energy spectrum, effective area,, dispersion matrix and stat
            print("    - lifetime:",obs.total_stats_safe_range.livetime.value)
            print("    - excess vector:",obs.excess_vector.data)
            if (plot): 
                obs.excess_vector.peek()
            plt.show()
        
        return


    @staticmethod
    def get_grb_properties(filename):
        """
        Get GRB properties
        """
        with open(filename,'r') as stream:
            data = yaml.load(stream)
        delta_t = data['time intervals'].split()
        intervals = []
        for idx in range(0,24,2):
            # in seconds
            intervals.append([float(delta_t[idx]), float(delta_t[idx + 1])] * u.s)
        data['time intervals'] = intervals
        return data

#------------------------------------------------------------------------------
class GammaRayBurstPop(object):
    """Class to store a GRB sample

    Run simulations and store results for a GRB population

    Parameters
    ----------
    filepath: `str`
        Path where GRB folders are stored
    absorption: `~gammapy.spectrum.models.Absorption`
        Utility handling EBL models
    """

    def __init__(self, filepath, absorption):
        # Load GRBs
        self.grb_list = []
        file_list = glob(filepath + '/*')
        for ifile in file_list:
            print(ifile)
            self.grb_list.append(GammaRayBurst.from_file(ifile, absorption))

    def run_simulations(self, perf,
                        emin=0.03 * u.TeV,
                        alpha=0.2,
                        random_state = 1): # Fixed random state to get same results
#                        random_state='random-seed'):
        """Run simulations"""
        for grb in self.grb_list:
            grb.run_simulation(perf, emin, alpha, random_state)




