# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:45:41 2022

@author: Stolar
"""

# List of source trigger dates between two dates

import numpy as np
from astropy.time import Time
import sys
from pathlib import Path

###############################################################################
class trigger_dates():
    
    #--------------------------------------------------------------------------
    def __init__(self, tstart, tstop, folder, seed=2022, fext=[".gz",".fits"]):
        
        if not Path(folder).is_dir():
            sys.exit("{} is not a valid folder".format(folder) )
                
        self.tstart   = tstart
        self.tstop    = tstop
        self.folder   = folder
        self.seed     = seed
        self.filelist = [f.name for f in Path(folder).iterdir() \
                         if f.is_file() and (f.suffix in fext)]
        
        print(" Found {:d} files in {}"
              .format(len(self.filelist),self.folder))
       
        return
    
    #--------------------------------------------------------------------------
    def generate(self):
        
        delta  = (self.tstop-self.tstart).jd

        np.random.seed(self.seed)
        days = np.random.random(len(self.filelist))*delta
        
        self.dates = self.tstart+days   
        
        print(40*"-")
        print(" Source positions read from files in : ",self.folder)
        print(" Generated dates ")
        print("  - from     : ", self.tstart)
        print("  - to       : ", self.tstop)
        print("  - Duration : {:5.2f} days ({:3.2f} yrs)".format(delta,delta/365.25))
        print(40*"-")
        
        return
        
    #--------------------------------------------------------------------------
    def dump2yaml(self):
        
        import datetime

        filename = "Trigger_{:d}-{:%Y_%m_%d_%H%M%S}-{:%Y_%m_%d_%H%M%S}.yaml".format(len(self.filelist),self.tstart.datetime,self.tstop.datetime)
        out = open(filename,"w")
            
        print(" Now dumping generated dates into :",filename)

        print("sources: {}".format(self.folder),file=out)
        print("created: {}".format(datetime.datetime.now()),file=out)
        print("seed: {:d}".format(self.seed),file=out)
        print("start: {}".format(self.tstart),file=out)
        print("stop: {}".format(self.tstop),file=out)
        for fname,d in zip(self.filelist, self.dates):
            print("{}: {:20.10f} # {}"
                  .format(fname,d.jd,Time(d.jd,format="jd",scale="utc").isot),file=out)
        out.close()
        print("Done!")    
            
        return filename
    
    #--------------------------------------------------------------------------
    def plot(self, nbin=100):
            
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,5))
        ax.hist(Time(self.dates,format="jd",scale="utc").datetime,bins=nbin,alpha=0.8)
        ax.set_xlabel("Date")
        ax.grid(which="both")
        plt.show()    
        
        return
        
###############################################################################
def read_from_yaml(folder, filename, nmax=10,  fext=[".gz",".fits"]):
    
    import yaml
    from yaml.loader import SafeLoader
    
    filelist = [f.name for f in Path(folder).iterdir() \
                     if f.is_file() and (f.suffix in fext)]
    infile = open(filename,"r")
    data =  yaml.load(infile, Loader=SafeLoader)
    
    print(40*"-")
    print(" Read dates from ",filename)    
    for f in filelist[: min(nmax,len(filelist))]:
        print(f," found: ",data[f])
    
    infile.close()

    return    

#------------------------------------------------------------------------------
def get_from_yaml(filename):
    
    import yaml
    from yaml.loader import SafeLoader
    infile = open(filename,"r")
    data =  yaml.load(infile, Loader=SafeLoader)
    
    return data # That's a dictionnary

###############################################################################
def get_trigger_dates(trigger):
    """
    

    Parameters
    ----------
    trigger : Path or float
        Either a time shift in days for all initial trigger dates or a yaml 
        file with a number of Julian days (mjd) for each file names .

    Returns
    -------
    trig_data : float or list
        DESCRIPTION.
    trig_abs : TYPE
        DESCRIPTION.

    """
    if type(trigger) == float: # The triger shift is a fixed number in days
        trig_data =  trigger #np.array(len(grblist)*[cf.trigger])
        trig_abs   = False    
    elif trigger.is_file():
            
        from trigger_dates import get_from_yaml
        trig_data  = get_from_yaml(trigger)
        trig_abs = True
    else:    
        sys.exit("Trigger keyword is not a float nor a file name")

    return trig_data, trig_abs

###############################################################################
if __name__ == "__main__":

    """
    A standalone function to generate trigger times
    """
        
    # Get this from the command line    
    tstart    = Time('2028-01-01T00:00:00', format='isot', scale='utc')
    tstop     = Time('2034-12-31T23:59:59', format='isot', scale='utc')
    seed      = 2022 #'random-seed'
    nmax      = 10 
    datafiles = "D:\CTAA\SoHAPPy\input\lightcurves\LONG_FITS"
    
    dates = trigger_dates(tstart, tstop, datafiles)
    
    dates.generate()
    
    filename = dates.dump2yaml()
    
    dates.plot()    
    
    # read_from_yaml(dates.folder, filename, nmax=1000000)    # For test 
    
    get_from_yaml(filename)
    
    print(" *** CONSIDER RECOMPUTING THE VISIBILITIES ***")