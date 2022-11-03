# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:24:09 2021

@author: Stolar

A standalone function to compute and dump visibilities to disk.
See visibility.yaml for possible parameters.
"""

import warnings
#warnings.filterwarnings('error') # Transform warnings into errors 
warnings.filterwarnings('ignore')

class Config():
    def __init__(self,filename):
        with open(filename,"r") as f:
            mydict = yaml.load(f)
            
        for k,v in mydict.items():
           # print(k,v)
           setattr(self, k, v)
        return
###############################################################################
if __name__ == "__main__":
    """
    
    This script can be used to save/read viibility classes in binary format
    using pickle. just use, for each individual instance, where appropriate :
        method = "binary"
        
        For saving the visibility `vis`
        import pickle
        pickle.dump(vis,outfile)
        
        For reading data into the visibility `vis`:
         vis =  pickle.load(infile)
    """

    import time
    import json, yaml
    from pathlib import Path
    from utilities import  source_ids
    from trigger_dates import get_trigger_dates
    from grb import GammaRayBurst
    from visibility import Visibility, params_from_key, object_to_serializable
    
    ###---------------------------
    ### Parameters
    ###---------------------------
    cf = Config("config_visibility.yaml")

    # GRB data
    grblist    = source_ids(cf.ifirst, cf.ngrb)
    grb_folder = Path(cf.infolder,cf.data_dir)
    trigger = Path(cf.infolder,cf.trigger) if isinstance(cf.trigger,str) \
                                           else float(cf.trigger)
    dt, dt_abs = get_trigger_dates(trigger)  
    
    # Visibility conditions
    vispar =  (True, params_from_key(cf.visibility))
    
    # Output file name
    Path(cf.resfolder).mkdir(parents=True, exist_ok=True)    
    name = "vis"+"_"+str(grblist[0])+"_"+str(grblist[-1])+"_"+cf.visibility+".json"
   
    print("Visibilities to/from '{}/{}' ".format(cf.resfolder, name)) 
    f_all = open(Path(cf.resfolder,name),"w" if cf.save else "r")
    
    ###---------------------------
    ### Let's go
    ###---------------------------
    start = time.time() # Starts chronometer
    
    if cf.save: vis_list = []
    if cf.read: data_vis = json.load(f_all)

    # Loop over sources
    for i, item in enumerate(grblist):
        
        print(item)
        name = "Event"+str(item)   

        # Get GRB
        grb = GammaRayBurst.from_fits(Path(grb_folder,name+".fits.gz"),
                                  ebl     = None,
                                  prompt  = None, 
                                  dt      = dt[name+".fits.gz"] if dt_abs else dt,
                                  dt_abs  = dt_abs)    
        
        ### Save visibilities into a list
        if cf.save: 
            for loc in ["North", "South"]:
                grb.set_visibility(loc, info=vispar, name="Json") # computed or from a binary
                vis_list.append(grb.vis[loc])

         ### Read visibiities
        if cf.read: 
            for loc in ["North", "South"]:
                grb.vis[loc] = Visibility.from_dict(data_vis[name+"_"+loc])
                
        ### Show visibiities
        if cf.dbg:
            for loc in ["North", "South"]:
                grb.vis[loc].print()
                import grb_plot as gplt
                gplt.visibility_plot(grb, loc=loc)    
                    
    ### Save all visibilities in a single dictionnary
    if cf.save:
        json.dump({v.name:v for v in vis_list}, f_all, 
                    default=object_to_serializable, indent=None)
        f_all.close()
    
    ### End
    stop = time.time() # Starts chronometer
    print("Completed in {:8.2f} s ({:4.2f} s per source)"
            .format(stop-start,(stop-start)/cf.ngrb))