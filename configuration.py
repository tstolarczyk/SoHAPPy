import os
import sys, getopt

import yaml
from yaml.loader import SafeLoader

import astropy.units as u
from pathlib import Path

from utilities    import Log

def_conf = "config.yaml"
def_vis  = "visibility.yaml"

###############################################################################
class Configuration(object):
    """
    This class handles all the external parameters required to run a simulation
    and the corresponding analysis.
    """
    
    ###------------------------------------------------------------------------    
    def __init__(self, argv, conf_file = def_conf, 
                             vis_file  = def_vis, copy=True, debug=False):
        """
        The defaut configuration file is found in the code repository, but it 
        can be changed for tests by changing the def_conf variable in the 
        configuration module. Same for the visibility parameter file.        

        Parameters
        ----------
        argv :  list
            DESCRIPTION.
        conf_file : String, optional
            Default configuration file name. The default is def_conf.
        vis_file : String, optional
            Default visibility parameter file name. The default is def_vis.
        copy : Boolean, optional
            If True, copy the configuration file to the output folder. 
            The default is True.
        debug : Boolean, optional
            If True, verbosy. The default is False.

        Returns
        -------
        None.

        """
        
        ### --------------------------
        ### Get command line arguments
        ### --------------------------
        try:
            opts, args = getopt.getopt(argv,
                                       "hc:n:f:N:d:o:i:D",
                                       ["config=",
                                        "ngrb=",
                                        "first=",
                                        "niter=",
                                        "dbg=",
                                        "output=",
                                        "input=",
                                        "days="])
        except getopt.GetoptError:
            print("No line arguments passed: Using default values")              

        ### --------------------------
        ### Check if help is needed
        ### --------------------------
        for opt, arg in opts:
            if opt == '-h':
                print(" SoHAPPy.py "
                      + "-c <config.file> "
                      + "-N <ngrb> "
                      + "-f <1st grb or list> "
                      + "-n <MC iterations> "
                      + "-o <Output folder> "
                      + "-i <Input folder> "
                      + "-d <debug> "
                      + "-D <date shift or file>")
                sys.exit()
                
        ### --------------------------
        ### Check if configuration file is given
        ### --------------------------            
        self.filename = None # Assume not given
        
        for opt, arg in opts:
            if opt in ("-c","--config"):
                if os.path.isfile(arg):
                    self.filename =  Path(arg)
                else:
                    sys.exit("Error: configuration file does not exist")

        if self.filename == None:
            # Default parameter read from an existing default file
            # Where the code is
            basefolder = Path(__file__).parent
            self.filename = Path(basefolder,conf_file)
            
        # Change filename to absolute Path
        self.filename = self.filename.resolve()
        
        ### --------------------------
        ### Read configuration file
        ### --------------------------            
        if (debug): print(" Now read config file ", self.filename)
        self.read()
        
        ### --------------------------
        ### Supersed config. parameters by command line
        ### --------------------------         
        for opt, arg in opts:
            if opt in ("-N", "--ngrb"):
                self.ngrb =  int(arg)
            elif opt in ("-f", "--first"):
                self.ifirst = int(arg)
            elif opt in ("-n", "--niter"):
                self.niter = int(arg)
            elif opt in ("-o", "--output"):
                self.res_dir = Path(arg)
            elif opt in ("-d", "--debg"):
                dbg = int(arg)
                self.dbg = dbg
            elif opt in ("-D", "--days"):
                self.trigger = arg 
                
        ### --------------------------
        ### deduce additionnal parameters
        ### --------------------------          
        
        # Define the input folder and files from the base and the subfolders
        self.data_dir   = Path(self.infolder,self.data_dir)
        self.swiftfile  = Path(self.infolder,self.swiftfile)
        self.irf_dir    = Path(self.infolder,self.irf_dir)
        if type(self.trigger) is str:
            self.trigger = Path(self.infolder,self.trigger)
            if not self.trigger.is_file():
                sys.exit("{} does not point to a valid file "
                         .format(self.trigger))
        else:
            self.trigger = float(self.trigger)
        
        # The visibility variable can be either:
        # - A subfolder where pre-computed visibility files can be found
        # - A keyword for computing the visibility from predefined parameters
        # - "None" in order to use the visibility in the data files
        if self.visibility != None:
                      
            if Path(self.infolder,self.visibility).is_dir():
                self.visibility = Path(self.infolder,self.visibility)
                print(" Visibilities read from disk, folder ",self.visibility)            
            elif self.visibility != "built-in":
                self.read_vis_param(visfilename = vis_file)
        else:
            sys.exit(" Visibility information not defined")
           
        # Create the show debugging flag from the general debug flag
        if (self.dbg < 0): 
            self.show = 0
            self.dbg = -self.dbg
        else: self.show = abs(self.dbg)

        # If debugging is requested, cannot be silent        
        if (self.dbg>0): self.silent = False

        # If the simulation is saved, it is not fluctuated
        if self.save_dataset or self.save_simu:
            self.do_fluctuate = False 
       
        # If no fluctuation, one simulation is enough !
        if self.do_fluctuate == False: 
            self.niter = 1
        
        return
        
    ###------------------------------------------------------------------------    
    def create_output_folder(self, log):
    
        # Get output folder absolute repository 
        # (from where the code is launched, not the code repository)
        self.res_dir = Path(self.res_dir).resolve()
        
        # Check that the output folder exists, otherwise try to create it
        if not self.res_dir.exists():
            log.warning("Creating {}".format(self.res_dir))
            try:
                Path.mkdir(self.res_dir)
            except:
                sys.exit(" Could not create {}".format(self.res_dir))
        else:
            log.warning(" Already exists :",self.res_dir)
        return
    
    ###------------------------------------------------------------------------    
    def read_vis_param(self, visfilename=None):
                
        """
        Read the visibility parameters from a yaml file or try to 
        use the default visibility in the input data files.
        
        Parameters
        ----------
        visfilename : String, optional
            The visibility yaml file where the visibility parameters are 
            stored along keywords. The default is None.

        Returns
        -------
        None.

        """
        
        # Try to get the parameters from a visibility yaml file
        print(" Visibilities defined by the keyword ",self.visibility)
        print(" Get visibility dictionnaries from ",visfilename)

        try:
            with open(visfilename) as f:
                visdict  = yaml.load(f, Loader=SafeLoader)
                # Check if the visibility variable correspond to an entry
                # Get the corresponding dictionnary values
                if self.visibility in visdict.keys():
                    self.visibility = visdict[self.visibility]           
                else:
                    sys.exit("{} found but keyword {} not referencec"
                             .format(visfilename,self.visibility))
        except IOError:
            sys.exit("{} not found".format(visfilename))
                
        return
    
    ###------------------------------------------------------------------------    
    def read(self):
        """
        Read configuration file

        Returns
        -------
        None.

        """

        #---------------------------------------------------
        def obj_dic(d):
            # top = type('new', (object,), d)
            seqs = tuple, list, set, frozenset
            for i, j in d.items():
                if isinstance(j, dict):
                    setattr(self, i, obj_dic(j))
                elif isinstance(j, seqs):
                    setattr(self, i, 
                        type(j)(obj_dic(sj) if isinstance(sj, dict) else sj for sj in j))
                else:
                    # if (j== "None"):
                    #     setattr(top,i,None)
                    # Not necessary, parsed correclty    
                    # elif (j== "False"):
                    #     setattr(top,i,False)
                    # elif (j== "True"):
                    #     setattr(top,i,True)
                    # else:
                    setattr(self, i, j)
                    
            return
        #---------------------------------------------------
        
        print(">>> Read configuration from ",self.filename)
        with open(self.filename) as f:
            data = yaml.load(f, Loader=SafeLoader)
            obj_dic(data)
            self.dtslew_North = u.Quantity(self.dtslew_North)
            self.dtslew_South = u.Quantity(self.dtslew_South )
            self.dtswift = u.Quantity(self.dtswift)
            self.arrays = {"North": self.array_North, "South":self.array_South}
            self.dtslew = {"North": self.dtslew_North, "South":self.dtslew_South}
            if self.Emax != None: self.Emax = u.Quantity(self.Emax)
 
        return
    ###------------------------------------------------------------------------ 
    def print(self,log):
        """
        

        Parameters
        ----------
        log : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        log.prt(" Configuration file      : {} ".format(self.filename))
        log.prt(" Simulation:")
        if self.save_simu:
             log.prt("     Save simulation     : {}".format(self.save_simu))
           
        if type(self.ifirst)!=list and not isinstance(self.ifirst, str):
            log.prt("     Number of sources*  : {:>5d}".format(self.ngrb))
            log.prt("     First source*       : {:>5d}".format(self.ifirst))
        else:
            log.prt("     Source list         : {}".format(self.ifirst))
        if type(self.trigger) == float:
            log.prt("     Date shift (days)*  : {:<10.2f}".format(self.trigger))
        elif self.trigger.is_file():
            log.prt("     Date shift file  *  : {}".format(self.trigger))
        else:
            sys.exit("Trigger keyword is not file nor float")
        log.prt("     Number of trials*   : {:>5d}".format(self.niter))
        log.prt(" EBL model               : {}".format(self.EBLmodel))
        log.prt(" Input/output :")
        log.prt("      Input folder*      : {}".format(self.infolder))
        # log.prt("      Output folder*     : {}".format(self.outfolder))
        log.prt("      Debug mode*        : {:>5d}".format(self.dbg))
        log.prt("      Show plots         : {:>5d}".format(self.show))
        log.prt("      Analysing files in : {}".format(self.data_dir))
        log.prt("      IRF files in       : {}".format(self.irf_dir))
        log.prt("      Results in*        : {}".format(self.res_dir))
        log.prt(" Site sub-arrays         : N:{} S:{}"
                .format(self.arrays["North"],self.arrays["South"]))
        log.prt(" Slewing time            : N:{} S:{}"
                .format(self.dtslew["North"],self.dtslew["South"]))
        log.prt("      Fixed              : {}".format(self.fixslew))
        log.prt("   SWIFT delay fixed     : {}".format(self.fixswift))
        if (self.fixswift):
            log.prt("                   value : {}".format(self.dtswift))
        else:
            log.prt("            Read from : {}".format(self.swiftfile))
            
        # if (self.method ==0): method = "On-Off Aperture photometry"
        # elif (self.method == 1): method = "On-off Energy dependent"
        log.prt(" Analysis (ndof)         : {}".format(self.method))
        if not self.forced_visible:
            if self.visibility == "built-in":
                print(" Default visibilities read from the GRB input files")
            elif isinstance(self.visibility, dict):
                print(" Visibilities recomputed on the fly :")
                # print(self.visibility)
                log.prt("   Vis. computed up to   : {} night(s)"
                        .format(self.visibility["depth"]))
                log.prt("   Skip up to            : {} night(s)"
                        .format(self.visibility["skip"]))
                log.prt("   Minimum altitude      : {}"
                        .format(self.visibility["altmin"]))
                log.prt("   Moon max. altitude    : {}"
                        .format(self.visibility["altmoon"]))
                log.prt("   Moon min. distance    : {}"
                        .format(self.visibility["moondist"]))
                log.prt("   Moon max. brightness  : {}"
                        .format(self.visibility["moonlight"]))               
        
            else:
                if not Path(self.visibility).is_dir():
                    print(self.visibility," folder does not exist")
                    sys.exit("Check the visibility attribute")
                else:    
                    print(" Visibilities read from disk, folder ",self.visibility)           
                
    
        log.prt("+----------------------------------------------------------------+")
        log.prt("|                 *: can be changed with command line (use -h)   |")
        log.prt("+----------------------------------------------------------------+")
        log.prt(" Developments:")
        if (self.save_grb == True):
            log.highlight("{:60s}"
            .format("Simulation saved to disk save_grb (save_grb = True)"))

        if (self.save_fig == True):
            log.highlight("{:60s}"
            .format("Plots saved as pdf booklet (fig_save = True)"))
    
        if (self.write_slices == True):
            log.highlight("{:60s}"
            .format("Slice information saved to disk (write_slices=True)"))
    
        if (self.signal_to_zero == True):
            log.warning("{:60s}"
            .format("Signal set to zero (signal_to_zero==True)"))
    
        if (self.do_fluctuate == False):
            log.warning("{:60s}"
            .format("No fluctuation in simulation (do_fluctuate==False)"))
    
        if (self.do_accelerate  == False):
            log.warning("{:60s}"
            .format("No abortion if first 10% undetected (do_accelarate==False)"))
        
        if (self.n_night != None):
             log.warning("GRB data limited to the first {} nights"
                    .format(self.n_night))
        if (self.Emax != None):
             log.warning("GRB energy bins limited to {}"
                    .format(self.Emax))
        if (self.fixed_zenith != None):
             log.warning("Zenith angle requested to be fixed at value '{:5s}'     "
                    .format(self.fixed_zenith))
        if (self.magnify !=1):
            log.warning("GRB flux values are multiplied by {}"
                    .format(self.magnify))        
        if (self.forced_visible):
            log.warning("{:60s]}".format("GRB always visible (infinite nights, always above horizon, no Moon) "))
    
        if (self.niter == 0):
            log.failure(" Cannot run simulation with ZERO trials")
            log.warning("{:60s]}".format(" Use other main specific scripts for tests"))
            sys.exit( " At least one trial is requested")
        if (self.test_prompt):
            log.warning("{:60s]}".format("Test prompt simulation"))
    
        log.prt("")
        
        return
    ###------------------------------------------------------------------------
    def write(self):
        name = Path(self.res_dir, os.path.basename(self.filename))
        print("<<< Writing configuration in",name)
        file = open(name,"w")
        
        strlist = [ "datafile",
                   "logfile",
                   "EBLmodel",
                   "swiftfile",
                   "altmin",
                  "altmoon",
                   "moondist",
                   "dtswift",      
                   "array_North",
                   "array_South",
                   "dtslew_North",
                   "dtslew_South"]
        dirlist = ["filename", "res_dir",  "grb_dir",   "irf_dir"]
        
        ## Special actions to be taken for vis_dir that can be None
        for k in vars(self):
            val = vars(self)[k]
            if k in (strlist): val = '"'+str(val)+'"'
            if k in (dirlist): val = '"'+Path(val).as_posix()+'"'
            if val == None: 
                val = "Null"
            if k != "arrays" and k!= "dtslew":
                print("{:16s}: {}".format(k,str(val)),file=file)
                  
        file.close()
        return

###############################################################################
if __name__ == "__main__":

    """
    A standalone function to read a configuration
    """
    
    cf = Configuration(sys.argv[1:])
    log_filename    = Path(cf.res_dir,cf.logfile)
    log = Log(name  = log_filename, talk=not cf.silent)    
    cf.print(log)
    cf.write()