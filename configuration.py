import sys, getopt
import numpy as np

import yaml
from yaml.loader import SafeLoader

import astropy.units as u
from pathlib import Path

from niceprint    import Log

def_conf = "config.yaml"

__all__=['Configuration']

###############################################################################
class Configuration(object):
    """
    This class handles all the external parameters required to run a simulation
    and the corresponding analysis.
    """
    
    # When dumpingthe Configuration class instance to a fileon disk
    # Keyword to be ignored
    ignore = ["filename","dtslew", "arrays", "dbg"]
    
    # Member of the class that should be handled as quantities and not strings
    quantities = ["dtslew_South","dtslew_North","dtswift"]

    ###------------------------------------------------------------------------    
    def __init__(self):
        """
        Parameter defintion and default values

        Returns
        -------
        None.

        """
        
        self.filename = "config.yaml"
        self.dbg = 0
        
        ### -----------------
        ### PHYSICS PARAMETERS
        ### -----------------
        self.ifirst  = 1  # id to start with, or a list (then nsrc is not used)
        self.nsrc    = 1  # integer, number of sources to be read
        self.trigger = 0. #real (days) or a file with a list of days 
      
        # visibility can be :
        # - "built-in" (read from the data file if it exists);
        # - a subfolder where to find the pre-computed visibility files;
        # - a json file with a colection of visibilities; 
        # - a keyword corresponding to a dictionnay in visibility.yaml
        #   to compute the visivility on the fly.
        self.visibility = "stricmoonveto"

        # Possible EBL models are from gammapy or gilmore data on disk. 
        self.EBLmodel   = "dominguez" 
                
        ### -----------------
        ### INPUT/OUPUT
        ### -----------------
        self.infolder   = Path("../input")  # Base inoput folder 
        self.resfolder  = Path("../output") # Base output main folder
        self.data_dir   = "lightcurves/LONG_FITS/" #Afterglow subfolder
        self.out_dir    = "test" # Result subfolder
        self.prompt_dir = None # Prompt subfolder (if None, not considered)        
        self.irf_dir    = "irf/Full/prod3-v2" # IRF subfolder
        self.datafile   = "data.txt"   # Population main output file name
         
        ### -----------------
        ### SIMULATION PARAMETERS
        ### -----------------
        self.niter  = 1 # Number of Monte Carlo trials 
        self.seed   = 2021  # Choose ‘random-seed’ to randomize   
       
        # If True Statistical fluctuations are enabled, 
        # If False niter forced to one
        self.do_fluctuate = False 
        
        # When True, the simulation is stopped if none of the first 
        # trials in the limit of 1 - det_level have reached the minimal
        # significance (3 sigma).        
        self.do_accelerate = True   

        # Observation position in the time slice
        self.obs_point       = "end" 

        ### -----------------
        ### DETECTION PARAMETERS
        ### -----------------
        
        # Fraction of trials above significance threshold to declare a 
        # detection at that threshold
        self.det_level    = 0.9
        
        # Number of on over off regions for detection
        self.alpha        = 0.2
        
        self.array_North  = "FullArray"  # IRF subarrays in North
        self.array_South  = "FullArray"  #IRF subarrays in South 
        self.dtslew_North = 30*u.s # Maximum slewing time delay in North
        self.dtslew_South = 30*u.s # Maximum slewing time delay in South
    
        self.dtslew = {"North": self.dtslew_North,
                       "South": self.dtslew_South}
        self.arrays = {"North": self.array_North,
                       "South": self.array_South}

    
        # If True use max. slewing value, otherwise a random delay < dtslew
        self.fixslew  = True   
        
        # Alert latency (e.g. SWIFT latency, with amean value of 77 s)
        self.dtswift  = 77*u.s     
        
        # If True, use above value.If False, latency generated from Swift data 
        self.fixswift = True 
        
        # SWIFT Latency distribution data file 
        self.swiftfile  = "swift/Swift_delay_times.txt"  
        
        # Gammapy datasets folder
        self.extra_dir = "gammapy-extra-master/datasets"

        ### -----------------
        ### DEBUGGING / BOOKKEEPING 
        ### -----------------

        # From 0 to 3, increasingly talkative
        # if negative or zero : no plot
        # 0: evt counting, 1: + some results, 2: + event details
        # 3: details for each trials
        self.dbg        = 1
        self.show       = 1
        self.logfile    = "analysis.log" # log file (results, status, warnings)

        self.save_simu  = False  # If True, Simulation saved to file 
        self.save_grb   = False  # If True, GRB class saved to disk 
        self.save_fig   = False  # If True, plots saved to pdf file 
        self.remove_tar = False  # If True, remove tarred files

        ### -----------------
        ### EXPERTS/DEVELOPPERS ONLY
        ### -----------------
        self.test_prompt     = False # If True test prompt alone (experimental)
        
        # Prompt characteristics from the afterglow with same id.
        self.use_afterglow   = False 
        
        self.fixed_zenith    = None  # If a value ("20 deg") freezes zenith
        self.magnify         = 1 # Multiplicative factor of the input flux
        
        # If True, nothing on screen (output to log (if dbg =0))
        self.silent          = True 
        
        # Store detailed information on slices if True        
        self.write_slices    = False 
        
        # self.forced_visible  = False # If True, the GRB is always visible (infinite nights)
        self.n_night         = None # Limit data to a maximal number of nights
        self.Emax            = None # Limit data energy bins to Emax                
            
        return   
    
    ###------------------------------------------------------------------------    
    @classmethod
    def build(cls, argv, conf_file = def_conf, copy=True, debug=False):
        """
        
        This acts as a construtor.
        The defaut configuration file is found in the code repository, but it 
        can be changed.        

        Parameters
        ----------
        argv :  list
            DESCRIPTION.
        conf_file : String, optional
            Configuration file name. The default is def_conf.
        copy : Boolean, optional
            If True, copy the configuration file to the output folder. 
            The default is True.
        debug : Boolean, optional
            If True, verbosy. The default is False.

        Returns
        -------
        None.

        """
        cls = Configuration()
        
        ### --------------------------
        ### Get command line arguments
        ### --------------------------
        try:
            opts, args = getopt.getopt(argv,
                                       "hc:n:f:N:d:o:i:D:V:",
                                       ["config=",
                                        "nsrc=",
                                        "first=",
                                        "niter=",
                                        "dbg=",
                                        "output=",
                                        "input=",
                                        "days=",
                                        "visibility="])
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
                      + "-D <date shift or file>"
                      + "-V <visibility keyword>")
                sys.exit()
                
        ### --------------------------
        ### Check if configuration file is given on command line
        ### --------------------------            
        cls.filename = None # Assume not given
        
        for opt, arg in opts:
            if opt in ("-c","--config"):
                if Path(arg).is_file():
                    cls.filename =  Path(arg)
                else:
                    sys.exit("{}.py: '{} configuration file does not exist"
                             .format(arg, __name__))
                    
        # Not on command line - get default config file 
        # Can be obtained from a file name or from an archive
        import tarfile
        if not type(conf_file) == tarfile.ExFileObject:
            if cls.filename == None:
                # Default parameter read from an existing default file
                basefolder = Path(__file__).parent  # Where the code is
                cls.filename = Path(basefolder,conf_file)
                
            cls.filename = cls.filename.resolve() #  absolute Path
        else:
            cls.filename = conf_file # In memory
        
        ### --------------------------
        ### Read configuration file from file name
        ### --------------------------            
        if debug: print(" Now read config file ", cls.filename)
        cls.read_from_yaml(filename = cls.filename)
        
        ### --------------------------
        ### Supersed config. parameters by command line
        ### --------------------------   
        for opt, arg in opts:
            if opt in ("-N", "--ngrb"):
                cls.nsrc =  int(arg)
            elif opt in ("-f", "--first"):
                cls.ifirst = int(arg)
            elif opt in ("-n", "--niter"):
                cls.niter = int(arg)
            elif opt in ("-o", "--output"):
                cls.out_dir = Path(arg)
            elif opt in ("-d", "--debg"):
                dbg = int(arg)
                cls.dbg = dbg
            elif opt in ("-D", "--days"):
                try:
                    cls.trigger = float(arg)
                except ValueError:
                    cls.trigger = arg
            elif opt in ("-V", "--visibility"):
                cls.visibility = arg 
                
        # Show plots or not        
        cls.show = 0 if cls.dbg  < 0 else abs(cls.dbg)
          
        # If debugging is requested, cannot be silent        
        if cls.dbg>0: cls.silent = False

        # If the simulation is saved, it is not fluctuated
        if cls.save_simu:  cls.do_fluctuate = False 
       
        # If no fluctuation, one simulation is enough !
        if cls.do_fluctuate == False:  cls.niter = 1
        
        return cls   
    
    ###------------------------------------------------------------------------    
    def read_from_yaml(self, filename=None):
        """
        Read configuration file from disk.
        
        Returns
        -------
        None.
        
        """
                     
        #---------------------------------------------------
        def obj_dic(d):
            
        
            # top = type('new', (object,), d)
            seqs = tuple, list, set, frozenset
            for k, v in d.items():
                
                # Do not read wrong/deprecated keywords
                if k not in self.__dict__.keys():
                    print(" warning: ",k," ignored")
                    continue
                    
                if isinstance(v, dict):
                    setattr(self, k, obj_dic(v))
                elif isinstance(v, seqs):
                    setattr(self, k, 
                        type(v)(obj_dic(sv) if isinstance(sv, dict) else sv for sv in v))
                else:
                    if k in self.quantities:    
                        setattr(self, k, u.Quantity(v))
                    else:
                        setattr(self, k, v)
                    
            return self
    #---------------------------------------------------
        if filename != None: self.filename = filename 
        print(">>> Read configuration from ",self.filename)
        
        try:
            f = open(filename, "r")
        except IOError:
            sys.exit("{}.py : {} does not exist."
                     .format(__name__,self.filename))
        
        data = yaml.load(f, Loader=SafeLoader)
        obj_dic(data)
        
        return
    
    ###------------------------------------------------------------------------    
    def write(self, filename = None):
        """
        Write current configuration to a yaml file

        Parameters
        ----------
        filename : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        if filename == None: filename = "config_backup.yaml"
        
        # Create a dict with unnecessary keywords dropped and no quantities
        newdict = {}
        for k,v in self.__dict__.items():
        
            if k in self.ignore: continue
        
            if isinstance(v,u.Quantity) or isinstance(v, Path):
                newdict[k] = str(v)
                
        f = open(filename, 'w')
        import yaml
        yaml.dump(newdict, f, sort_keys=False)   
        
        return
    

        
    ###------------------------------------------------------------------------ 
    def print(self,log):
        """
        Print confuguartion class instance to screen and log file

        Parameters
        ----------
        log : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        log.prt(" Configuration file*     : {} ".format(self.filename))
            
        #----------------------------------------------------
        def title(s):
            log.prt("")
            log.prt("{:<27s} {}".format("=== "+s,40*"="))
            return
        #----------------------------------------------------
            
        ### -----------------
        ### PHYSICS PARAMETERS
        ### -----------------
        title("Physics parameters")
        if type(self.ifirst)!=list and not isinstance(self.ifirst, str):
            log.prt(" Number of sources*         : {:>5d}".format(self.nsrc))
            log.prt(" First source*              : {:>5d}".format(self.ifirst))
        else:
            log.prt(" Source list*               : {}".format(self.ifirst))        

        if isinstance(self.trigger, float) or isinstance(self.trigger, int) :
            log.prt(" Date shift (days)*         : {:<10.2f}"
                    .format(self.trigger))
        elif Path(self.infolder,self.trigger).is_file():
            log.prt(" Date shift file*           : {}"
                    .format(self.trigger))
        else:
            sys.exit("{}.py: Trigger keyword is not file nor float"
                     .format(__name__))        
            
        log.prt(" Visibility*                : {}".format(self.visibility))
        log.prt(" EBL model                  : {}".format(self.EBLmodel))
        
        ### -----------------
        ### INPUT/OUPUT
        ### -----------------  
        title("Input/output")
        log.prt(" Input folder*              : {}".format(self.infolder))
        log.prt(" Output folder*             : {}".format(self.resfolder))
        log.prt(" Data subfolder             : {}".format(self.data_dir))
        log.prt(" IRF subfolder              : {}".format(self.irf_dir))
        
        if self.prompt_dir != None:
            if Path(self.infolder,self.prompt_dir).is_dir():
                log.prt(" Prompt data subfolder      : {}"
                        .format(self.prompt_dir))
            else:
                sys.exit("Prompt data : {} not a valid subfolder"
                         .format(self.prompt_dir))        
        else:
            log.prt(" Prompt data                : not considered")
        log.prt(" Output population          : {}".format(self.datafile))
            
        ### -----------------
        ### SIMULATION PARAMETERS
        ### -----------------        
        title("Simulation")
        if not self.niter:
            sys.exit( "{}.py: At least 1 trial is requested".format(__name__))        
        else:
            log.prt(" Number of trials*          : {:>5d}".format(self.niter))
        log.prt(" Seed                       : {}".format(self.seed))

        log.prt(" Counts randomised          :{} ".format(self.do_fluctuate),end="")        
        if not self.do_fluctuate: log.warning("-> Only 1 trial")
        else: log.prt("")

        log.prt(" Stop if cannot be detected : {}".format(self.do_accelerate))
    
        ### -----------------
        ### DETECTION PARAMETERS
        ### -----------------
        title("Detection")
        log.prt(" Detection level            : {}".format(self.det_level))
        log.prt(" Alpha (1/n)                : {}".format(self.alpha))
        log.prt(" Site sub-arrays            : N:{} S:{}"
                .format(self.array_North,self.array_South))
        log.prt(" Slewing time               : N:{} S:{}"
                .format(self.dtslew_North,self.dtslew_South))
        log.prt("    Fixed                   : {}".format(self.fixslew))
        log.prt(" SWIFT delay fixed          : {}".format(self.fixswift))
        if (self.fixswift):
            log.prt("    Value                   : {}".format(self.dtswift))
        else:
            log.prt("            Read from : {}".format(self.swiftfile))

        ### -----------------
        ### DEBUGGING / BOOKKEEPING 
        ### -----------------
        title("Debugging / bookkeeping")
        log.prt(" Debug mode*                : {:>5d}".format(self.dbg))
        log.prt(" Log file                   : {}".format(self.logfile))
        log.prt(" Show plots                 : {:>5d}".format(self.show))
        log.prt(" Save simulation            : {} ".format(self.save_simu),end="")
        if self.save_simu: log.warning("-> Only 1 trial")
        else: log.prt("")
        log.prt(" Save source class          : {}".format(self.save_grb))
        log.prt(" Save figures to pdf        : {}".format(self.save_fig))
        log.prt(" Remove tarred file         : {}".format(self.remove_tar))

        log.prt("+"+66*"-"+"+")
        log.prt(" *: can be changed with command line (use -h)")
        log.prt("+"+66*"-"+"+")
        log.prt(" Developments:")

        ### -----------------
        ### EXPERTS/DEVELOPPERS ONLY
        ### -----------------    
        if self.test_prompt:
            log.warning("{:60s]}".format("Test prompt simulation"))
            if self.use_afterglow:
                log.warning("{:60s]}".format("-> use afterglow information"))               
        if self.fixed_zenith != None:
             log.warning("Zenith angle requested to be fixed at value '{}'"
                    .format(self.fixed_zenith))
        if self.magnify !=1:
            log.warning("GRB flux values are multiplied by {}"
                    .format(self.magnify))                
        if not self.silent:
            log.warning("Maximise screen output")
            
        if self.write_slices == True:
            log.highlight("{:60s}"
            .format("Slice information saved to disk (write_slices=True)"))

        if self.n_night != None:
             log.warning("GRB data limited to the first {} nights"
                    .format(self.n_night))
        if self.Emax != None:
             log.warning("GRB energy bins limited to {}"
                    .format(self.Emax))
    
        log.prt("    --end\n")
        
        return
    
    ###------------------------------------------------------------------------    
    def create_output_folder(self, log):
        """
        Create the code outptut folder

        Parameters
        ----------
        log : TYPE
            DESCRIPTION.

        Returns
        -------
        res_dir : TYPE
            DESCRIPTION.

        """
    
        res_dir = Path(self.resfolder, self.out_dir).resolve()
        
        # Check that the output folder exists, otherwise try to create it
        if not res_dir.exists():
            log.warning("Creating {}".format(res_dir))
            try:
                Path.mkdir(res_dir)
            except:
                sys.exit("{}.py: Could not create {}".format(__name__,res_dir))
        else:
            log.warning(" Already exists :",res_dir)
            
        return res_dir   

    ###------------------------------------------------------------------------
    def source_ids(self):
        """
        Obtain the source list from the input parameters.
        `first`is either an integer (identifier of the first source to be analysed)
        or list with integer (identifiers) or string (source name).
    
        Returns
        -------
        A list of source identifiers to be analysed.
    
        """
        if type(self.ifirst)!=list:
            if isinstance(self.ifirst, str):
                srclist = [self.ifirst]
            elif isinstance(self.ifirst, int):
                srclist = list(range(self.ifirst,self.ifirst+self.nsrc))
        else:
            srclist = self.ifirst   
            
        return srclist    

    ###------------------------------------------------------------------------
    def decode_keyword(self, debug=False):
        
        """
        Decode the `visibility`keyword value and return an information
        to be handled at running time.
        
        The keyword parameter can be:
        * a json file name, where all visibility members of a data subset is 
        stored.
        * 'force' : visibility forced to be maximal (infinite nights)
        And for backward compatibility:
        * 'built-in' : will use the default keyword stored in the source
        files, if any.
        * a folder name, where class instances produced beforehand
        can be reloaded as binary files from the folder on disk.
        * a keyword giving access to information in a yaml file, as a python
        dictionnary. 

        Parameters
        ----------
        debug : boolean, optional
            If True, talk a bit. The default is False.

        Returns
        -------
        Various type of information in a single variable
            The information necessary to get the visibility at running time.
            Can be a string or a dictionnary.

        """
    
        # Test whether the keyword corresponds to a file or a folder
        test = Path(self.infolder,self.visibility)
        
        # print(" Visibilities read from json file : {}".format(test))
        if test.is_file() and test.suffix == ".json":

            if isinstance(self.ifirst, list): 
                range = [min(self.ifirst), max(self.ifirst)]
            elif isinstance(self.ifirst,int):
                range= [self.ifirst, self.ifirst+self.nsrc-1 ]
            else:
                sys.exit("{}.py : check visibility keyword. Recomputation seems required")
                
            # Crude check to see if the visibility file matches the grb list
            # as in the configuration file
            bnds = [int(s) for s in test.name.split("_") if s.isdigit()]
            if len(bnds)!=2:
                sys.exit(" {}.py : {} does not follow naming scheme"
                .format(__name__, test.name))
            elif range[0] < bnds[0] or  range[1] > bnds[1]:
                    sys.exit("{}.py : {} does not match the event list ({}-{})"
                             .format(__name__, test.name, range[0], range[1]))

            # Get visibility from file given in the configuration file
            with open(test,"r") as f:
                import json
                return json.load(f)
            
        # print(" Each keyword is read from disk as a binary, folder : {}".format(test))
        elif test.is_dir(): return test     
         
        else: # Not a folder nor a file
            # print(" Visibilities from input data file")
            # print(" Visibilities forced to be maximal")
            if self.visibility in ["built-in","forced"]: 
                return self.visibility           
                      
            # print(" Visibilities computed from keyword : '{}' "
            #   .format(self.visibility))
            else:   
                from visibility import params_from_key
                vispar =  params_from_key(self.visibility)
                return vispar

        return            

    ###------------------------------------------------------------------------    
    def get_delay(self):

        """
        Compute the overall delay to be applied to the start of detection
        (satellite and telescope slewing), according to the user parameters.

        Returns
        -------
        delay : Quantity (time)
                Delay before the detection can start.

        """
        delay = {"North":0*u.s, "South":0*u.s}

        for loc in ["North", "South"]:

            dt = 0*u.s
            if (self.fixslew): dt = self.dtslew[loc]
            else:              dt = self.dtslew[loc]*np.random.random()

            if (self.fixswift): dt = dt + self.dtswift # don't do += !!!
            else: sys.exit("{}.py: Variable SWIFT delay not implemented)"
                           .format(__name__))
            delay[loc] = dt.to(u.s)

        return delay
    
###############################################################################
if __name__ == "__main__":

    """
    A standalone function to read a configuration
    """
    log=Log("dummy.log")
    cf1 = Configuration() # Default
    # cf1.print(log)
    
    cf1.read_from_yaml(filename="config.yaml")
    # cf1.print(log)
    cf1.write(filename = "config_bck.yaml")
    
    cf2 = Configuration()
    cf2.read_from_yaml(filename="config_bck.yaml")
    cf2.print(log)

    cf4 = Configuration.build(sys.argv[1:])
    cf4.print(log)
