# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
import matplotlib.pyplot as plt
import sys

from grb_utils import GammaRayBurst
import numpy as np

class GRBMonteCarlo():
    
    def __init__(self,
                 niter,
                 alpha,
                 cta_perf,
                 grb
                 ):
        
        # Input parameters and object
        self.niter       = niter
        self.alpha       = alpha
        self.cta_perf    = cta_perf
        self.grb         = grb
        
        # List, mean and rms related to trials
        # NOTE : list could be dropped and become local variables
        self.sigmax_list    = []
        self.sigmax = 0
        self.err_sigmax = 0
        self.nex_sigmax_list  = []
        self.nb_sigmax_list  = []
        self.t_sigmax_list  = []
        self.nex_sigmax = -1
        self.nb_sigmax = -1
        self.t_sigmax = 0
        self.err_t_sigmax = 0
        
        self.nex_3sigma_list = []
        self.nb_3sigma_list = []
        self.t_3sigma_list  = []
        self.nex_3sigma = -1
        self.nb_3sigma = -1
        self.t_3sigma = -1e9
        self.err_t_3sigma = -1e9
        
        self.nex_5sigma_list = []
        self.nb_5sigma_list = []
        self.t_5sigma_list  = []
        self.t_5sigma = -1
        self.nex_5sigma = -1
        self.nb_5sigma = -1
        self.err_t_5sigma = -1e9
        self.sigfull_list   = []
        

        # Related to time and time bins
        self.detect_3sigma = 0  # Fraction of time the 3 sigma level was reached
        self.detect_5sigma = 0  # Fraction of time the 3 sigma level was reached
        
        # Mean sigma versus time plot
        self.sigma_mean = []
        self.sigma_rms  = []
        
        return
        

    
    def result(self,out=sys.stdout):

        print("\n========================== RESULTS",file=out)
        print("  GRB name               : ",self.grb.name,file=out)
        print("  Redshift             z : ",self.grb.z,file=out)       
        print("==========================  ",file=out)
        print("  Max. significance       : ",round(self.sigmax,1),
              "+/-",round(self.err_sigmax,1),file=out)
        print("                  at time : ",round(self.t_sigmax,1),
              "+/-",round(self.err_t_sigmax,1),file=out)
        print("  3 sigma reached at time : ",round(self.t_3sigma,1),
              "+/-",round(self.err_t_3sigma,1),file=out)              
        print("                  Fraction:",100*self.detect_3sigma,"%",file=out)
        print("  5 sigma reached at time : ",round(self.t_5sigma,1),
              "+/-",round(self.err_t_5sigma,1),file=out)              
        print("                  Fraction:",100*self.detect_5sigma,"%",file=out)
        print("==========================\n",file=out)

        return


#------------------------------------------------------------------------------
# Run simulations of the current grb
# Compute significance on the full observation, 
# significance max on intervals,
# Time of siginificance max, time of 3 sigma reached
    def run(self, dbg_level):
        
        iMC=1
        prt_frequency = 10
        
        n_3sigma = 0
        n_5sigma = 0
        
        sigma_sum = 0
        sigma2_sum = 0 
        
        while(iMC <= self.niter):
            if (self.niter <= 10) or (np.mod(iMC,prt_frequency) == 0): 
                print(iMC," ",end="")
                
            self.grb.run_simulation(self.cta_perf,alpha=self.alpha)
            obs_stat = self.grb.get_cumulative_stats()
        
            # Best sensitivity for the trial
            sigmax     = max(obs_stat['sigma'])
            t_sigmax   = obs_stat['livetime'][ (np.where( obs_stat['sigma']==sigmax)[0][0]) ]  
            non        = obs_stat['n_on'][ (np.where( obs_stat['sigma']==sigmax)[0][0]) ]
            noff       = obs_stat['n_off'][ (np.where( obs_stat['sigma']==sigmax)[0][0]) ]
            self.nex_sigmax_list.append(non -self.alpha*noff)
            self.nb_sigmax_list.append(self.alpha*noff)
            self.sigmax_list.append(sigmax)
            self.t_sigmax_list.append(t_sigmax)
            if (dbg_level) : print("\n --- Max significance = ",sigmax," at tmax=",t_sigmax)
            

            # 3 sigma reached for the trial
            mask3s = obs_stat['sigma']>=3
            if (mask3s.any()):
                t_3sigma   = obs_stat['livetime'][mask3s][0]
                non = obs_stat['n_on'][mask3s][0]
                noff = obs_stat['n_off'][mask3s][0]
                nex_3sigma = non - self.alpha*noff 
                nb_3sigma  = self.alpha*noff
                n_3sigma +=1
            else:
                t_3sigma = -1e9
                nex_3sigma = -1
                nb_3sigma = -1
            self.nex_3sigma_list.append(nex_3sigma)
            self.nb_3sigma_list.append(nb_3sigma)
            self.t_3sigma_list.append(t_3sigma)
            if (dbg_level) : print(" --- Alert time (3 sigma) =",t_3sigma)
            

            # 5 sigma reached for the trial
            mask5s = obs_stat['sigma']>=5
            if (mask5s.any()):
                t_5sigma = obs_stat['livetime'][mask5s][0]
                non = obs_stat['n_on'][mask5s][0]
                noff = obs_stat['n_off'][mask5s][0]
                nex_5sigma = non - self.alpha*noff 
                nb_5sigma  = self.alpha*noff
                n_5sigma +=1
            else:
                t_5sigma = -1e9
                nex_5sigma = -1
                nb_5sigma = -1

            self.nex_5sigma_list.append(nex_5sigma)
            self.nb_5sigma_list.append(nb_5sigma)
            self.t_5sigma_list.append(t_5sigma)
            if (dbg_level) : print(" --- Alert time (5 sigma) =",t_5sigma)

            # Overall sensitivity for the trial - not exciting
            sig_full = obs_stat['sigma'][len(obs_stat['sigma'])-1]
            self.sigfull_list.append(sig_full)
            if (dbg_level) : print(" --- Overal signif. = ",sig_full)
            
            
            # Acummulate sum and sum**2 for each observation (time bins)
            sigma_sum += obs_stat['sigma']
            sigma2_sum += obs_stat['sigma']**2
            
        
            iMC+=1
            
        ### Compute mean sigma and rms for each observation
        self.sigma_mean = sigma_sum/self.niter
        sigma2_mean = sigma2_sum/self.niter
        self.sigma_rms = np.sqrt(sigma2_mean-self.sigma_mean**2)
        # print("Mean sigma =",sigma_mean,"rms=",sigma_rms)
        
        ### Compute mean values and erros
        self.sigmax         = round(np.mean(self.sigmax_list),3)
        self.err_sigmax     = round(np.std( self.sigmax_list),3)
        self.nex_sigmax     = round(np.std(self.nex_sigmax_list),3)
        self.nb_sigmax      = round(np.std(self.nb_sigmax_list),3)
        self.t_sigmax       = round(np.mean(self.t_sigmax_list),3)
        self.err_t_sigmax   = round(np.std( self.t_sigmax_list),3)
        
        t_detected   = [time for time in self.t_3sigma_list   if time >= 0]
        nex_detected = [nex  for nex  in self.nex_3sigma_list if nex>=0]
        nb_detected  = [nb   for nb   in self.nb_3sigma_list  if nb>=0]
        if (len(t_detected)):
            self.nex_3sigma     = round(np.std(nex_detected),3)
            self.nb_3sigma      = round(np.std(nb_detected),3)
            self.t_3sigma       = np.mean(t_detected)
            self.err_t_3sigma   = np.std( t_detected)
        
        t_detected =   [time for time in self.t_5sigma_list if time >= 0]
        nex_detected = [nex  for nex  in self.nex_5sigma_list if nex>=0]
        nb_detected  = [nb   for nb   in self.nb_5sigma_list  if nb>=0]

        if (len(t_detected)):
            self.nex_5sigma     = round(np.std(self.nex_5sigma_list),3)
            self.nb_5sigma      = round(np.std(self.nb_5sigma_list),3)
            self.t_5sigma       = round(np.mean(t_detected),3)
            self.err_t_5sigma   = round(np.std( t_detected),3)
        
        ### Compte detection level
        self.detect_3sigma = n_3sigma/self.niter
        self.detect_5sigma = n_5sigma/self.niter
          
        return
                
                
#------------------------------------------------------------------------------
# plot simulations of the current grb    
    def plot(self, saveplots, out_filename):

        fig = plt.figure(figsize=(12,18))
        #fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)
    
        
        a1 = plt.subplot(321)
        a1.set_xlabel('Observation duration (s)')
        a1.set_ylabel('Mean significance')
        a1.set_title('Significance for '+str(self.niter)+' realisations')
        a1.set_xscale("log", nonposx='clip')
        a1.grid(which='both')
        plt.errorbar(self.grb.get_cumulative_stats()['livetime'],self.sigma_mean,yerr=self.sigma_rms,fmt='o') # Use last simulation time interval - they are all the same
        
        # Plot max significance
        a2 = plt.subplot(322)
        a2.set_xlabel('Standard deviation (Li&Ma)')
        a2.set_ylabel('#Realisation')
        a2.set_title('Max. significance- Mean = '
                     + str( round(self.sigmax,1))
                     + ' +/- '
                     + str( round(self.err_sigmax,1)))
        a2.grid(which='both')
        plt.hist(self.sigmax_list,color="grey")
        #y,binEdges = np.histogram(sigmax_list)
        #bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        #plt.bar(bincenters, y, yerr=np.sqrt(y))
        
        # Plot excess counts at  max significance
        a3 = plt.subplot(323)
        a3.set_xlabel('Counts')
        a3.set_ylabel('#Realisation')
        a3.set_title('S & B at Max. significance')
        a3.grid(which='both')
        plt.hist(np.add(self.nex_sigmax_list,self.nb_sigmax_list),bins=25,color="grey",alpha=0.5,label="Excess+Bckgd")
        plt.hist(self.nb_sigmax_list, bins=25,color="pink",alpha=0.5,label="Bckgd")
        plt.hist(self.nex_sigmax_list,bins=25,color="green",alpha=0.5,label="Excess")
        plt.legend()
#        print(np.add(self.nex_sigmax_list,self.nb_sigmax_list))
#        print(self.nb_sigmax_list)
#        print(self.nex_sigmax_list)

        # Plot time to get 3 sigma for detected cases
        a5 = plt.subplot(325)
        a5.set_xlabel('Observation duration (s)')
        a5.set_ylabel('#Realisation')
        a5.set_title('Time 3s - Mean = '
                     + str( round(self.t_3sigma,1))
                     + ' +/- '
                     + str( round(self.err_t_3sigma,1))
                     + ' ('
                     + str(round(100*self.detect_3sigma,1))
                     +'%)')
        a5.grid(which='both')
        t_detected = [time for time in self.t_3sigma_list if time >= 0]
        plt.hist(t_detected, color="grey")
 
        # Plot time to get 5 sigma for detected cases
        a6 = plt.subplot(326)
        a6.set_xlabel('Observation duration (s)')
        a6.set_ylabel('#Realisation')
        a6.set_title('Time 5s - Mean = '
                     + str( round(self.t_5sigma,1))
                     + ' +/- '
                     + str( round(self.err_t_5sigma,1))
                     + ' ('
                     + str(round(100*self.detect_5sigma,1))
                     +'%)')
        a6.grid(which='both')
        t_detected = [time for time in self.t_5sigma_list if time >= 0]
        plt.hist(t_detected, color="grey")       
        
#        # Plot overall significance - not exciting
#        a4 = plt.subplot(224)
#        a4.set_xlabel('Standard deviation (Li&Ma)')
#        a4.set_ylabel('#Realisation')
#        a4.set_title('Overall significance - mean = '
#                     + str( round(np.mean(self.sigfull_list),1))
#                     + ' +/- '
#                     + str( round(np.std(self.sigfull_list),1)))
#        a4.grid(which='both')
#        plt.hist(self.sigfull_list,color="grey")
        
        
        
        
        
        if (saveplots>0): plt.show()
        if (saveplots>1):
            #fig.savefig(out_filename+'.eps',transparent=True)
            fig.savefig(out_filename+'.jpg')
            
        return