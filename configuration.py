

import os
import sys, getopt

import yaml
from yaml.loader import SafeLoader

import astropy.units as u
from pathlib import Path

from utilities    import Log, warning
###############################################################################
class Configuration(object):
    """
    """
    ###------------------------------------------------------------------------    
    def __init__(self, argv, def_file="config.yaml", debug=False):
        """
        The defalut configuartion file is found in the code repository

        Parameters
        ----------
        argv : list
            Command line arguments
        def_file: string
            Default configuration file name without path 
            (Assumed to be sored where SoHAPPy scripts are)
        Returns
        -------
        None.

        """
        
        ### --------------------------
        ### Get command line arguments
        ### --------------------------
        try:
            opts, args = getopt.getopt(argv,
                                       "hc:n:f:N:d:o:",
                                       ["config=",
                                        "ngrb=",
                                        "first=",
                                        "niter=",
                                        "dbg=",
                                        "output=",])
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
                      + "-d <debug> ")
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
            self.filename = Path(basefolder,def_file)
            
        # Change filename to absolute Path
        self.filename = self.filename.resolve()
        self.read()
        if (debug): print(" Now read config file ", self.filename)
        
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
                self.dbg = abs(dbg)
    
        ### --------------------------
        ### deduce additionnal parameters
        ### --------------------------          
        
        # Create the show debugging flag from the general debug flag
        if (self.dbg < 0): self.show = 0
        else: self.show = abs(self.dbg)
    
        if (self.do_fluctuate == False): self.niter = 1
        #if (self.niter == 1): self.do_fluctuate=False
        if (self.dbg>0): self.silent = False
    
        # Avoid writing mutliple datasets when iteration number > 1
        if (self.niter > 1):
            self.save_dataset = False

        return
        
    ###------------------------------------------------------------------------    
    def create_output(self):
    
        # Get output folder absolute repository 
        # (from where the code is launched,not the code repository)
        self.res_dir = Path(self.res_dir).resolve()
        
        # Check that the output folder exists, otherwise try to create it
        if not self.res_dir.exists():
            warning("Creating {}".format(self.res_dir))
            try:
                Path.mkdir(self.res_dir)
            except:
                sys.exit(" Could not create {}".format(self.res_dir))
        else:
            warning(" Already exists :",self.res_dir)
        return
    ###------------------------------------------------------------------------    
    def read(self):

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
            self.altmin = u.Quantity(self.altmin)
            self.altmoon = u.Quantity(self.altmoon)
            self.moondist = u.Quantity(self.moondist)
            self.dtslew_North = u.Quantity(self.dtslew_North)
            self.dtslew_South = u.Quantity(self.dtslew_South )
            self.dtswift = u.Quantity(self.dtswift)
            self.arrays = {"North": self.array_North, "South":self.array_South}
            self.dtslew = {"North": self.dtslew_North, "South":self.dtslew_South}

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
        if type(self.ifirst)!=list and not isinstance(self.ifirst, str):
            log.prt("     Number of sources*  : {:>5d}".format(self.ngrb))
            log.prt("     First source*       : {:>5d}".format(self.ifirst))
        else:
            log.prt("     Source list         : {}".format(self.ifirst))
        log.prt("     Number of trials*   : {:>5d}".format(self.niter))
        log.prt(" EBL model               : {}".format(self.EBLmodel))
        log.prt(" Input/output :")
        log.prt("      Debug mode*        : {:>5d}".format(self.dbg))
        log.prt("      Show plots*        : {:>5d}".format(self.show))
        log.prt("      Analysing files in : {}".format(self.grb_dir))
        log.prt("      IRF files in       : {}".format(self.irf_dir))
        log.prt("      Result folder*     : {}".format(self.res_dir))
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
        # if (self.method ==0):
        #     method = "On-Off Aperture photometry"
        # elif (self.method == 1):
        #     method = "On-off Energy dependent"
        log.prt(" Analysis (ndof)         : {}".format(self.method))
        if (self.vis_dir != None):
            log.prt(" Vis. read from          : {}".format(self.vis_dir))
        elif (self.vis_cmp):
            # vis_self.print(log)
            log.prt(" Vis. computed up to     : {} night(s)".format(self.depth))
            log.prt(" Skip up to              : {} night(s)".format(self.skip))
            log.prt(" Minimum altitude        : {}".format(self.altmin))
            log.prt(" Moon max. altitude      : {}".format(self.altmoon))
            log.prt(" Moon min. distance      : {}".format(self.moondist))
            log.prt(" Moon max. brightness    : {}".format(self.moonlight))
    
        else:
            log.prt(" Visibility              : default")
    
        log.prt("+----------------------------------------------------------------+")
        log.prt("|                 *: can be changed with command line (use -h)   |")
        log.prt("+----------------------------------------------------------------+")
        log.prt(" Developments:")
        if (self.save_grb == True):
            log.highlight("Simulation saved to disk save_grb (save_grb = True)       ")
    
        if (self.write_slices == True):
            log.highlight("Slice information saved to disk (write_slices=True)       ")
    
        if (self.signal_to_zero == True):
            log.warning(  "Signal set to zero (signal_to_zero==True)                 ")
    
        if (self.do_fluctuate == False):
            log.warning(  "No fluctuation in simulation (do_fluctuate==False)        ")
    
        if (self.do_accelerate  == False):
            log.warning(  "No abortion if first 10% undetected (do_accelarate==False)")
    
        if (self.fixed_zenith != None):
             log.warning(  "Zenith angle requested to be fixed at value '{:5s}'     "
                    .format(self.fixed_zenith))
        if (self.magnify !=1):
            log.warning(  "GRB flux values are multiplied by {}"
                    .format(self.magnify))        
        if (self.forced_visible):
            log.warning(  "GRB always visible (infinite nights, always above horizon, no Moon) ")
    
        if (self.niter == 0):
            log.failure(  " Cannot run simulation with ZERO trials")
            log.warning(  " Use other main specific scripts for tests")
            sys.exit( " At least one trial is requested")
        if (self.test_prompt):
            log.warning(  "Test prompt simulation")
    
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
# arrays          : {'North': 'FullArray', 'South': 'FullArray'}
# dtslew          : {'North': <Quantity 30. s>, 'South': <Quantity 30. s>}

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

if __name__ == "__main__":


    """
    A standalone function to read a configuration
    """
    
    cf = Configuration(sys.argv[1:])
    log_filename    = Path(cf.res_dir,cf.logfile)
    log = Log(name  = log_filename, talk=not cf.silent)    
    cf.print(log)
    cf.write()