import sys
import argparse        
import numpy as np

import yaml
from yaml.loader import SafeLoader

import astropy.units as u
from pathlib import Path

from niceprint import warning

__all__=['Configuration']

###############################################################################
class Configuration(object):
    """
    This class handles all the external parameters required to run a simulation
    and the corresponding analysis. It also handles the tools to build the 
    input and output folder names, and create the necessary output folders 
    if they do not exist. 
    """

    # Keyword to be ignored when dumping the Configuration class instance to a 
    # file on disk
    ignore    = ["filename","dtslew", "arrays", "dbg", "srclist"]

    # Member of the class that should be handled as quantities and not strings
    quantities = ["dtslew_South","dtslew_North","dtswift"]
    
    # Default configuration file name
    def_conf   = "config.yaml" 

    ###-----------------------------------------------------------------------
    def __init__(self):
        """
        Parameter definition and default values

        Returns
        -------
        None.

        """

        self.filename = None
        self.dbg = 0

        ### -----------------
        ### PHYSICS PARAMETERS
        ### -----------------
        self.ifirst  = 1  # id to start with, or a list (then nsrc is not used)
        self.nsrc    = 1  # integer, number of sources to be read
        
        # visibility can be :
        # - "built-in" (read from the data file if it exists);
        # - a subfolder where to find the pre-computed visibility files;
        # - a json file with a colection of visibilities;
        # - a keyword corresponding to a dictionnay in visibility.yaml
        #   to compute the visivility on the fly.
        self.visibility = "strictmoonveto"

        # Possible EBL models are from gammapy or Gilmore data on disk.
        self.EBLmodel   = "dominguez"

        ### -----------------
        ### INPUT/OUPUT
        ### -----------------
        self.infolder   = Path("../input")  # Base input folder
        self.resfolder  = Path("../output") # Base output main folder
        self.data_dir   = "lightcurves/LONG_FITS/" #Afterglow subfolder
        self.out_dir    = "test" # Result subfolder
        self.prompt_dir = None # Prompt subfolder (if None, not considered)
        self.irf_dir    = "irf/Full/prod3-v2" # IRF subfolder
        self.datafile   = "data.txt"   # Population main output file name
        self.prefix     = "Event" # Prefix file name, followed by a number
        self.suffix     = "" # Suffix file name, following a number

        ### -----------------
        ### SIMULATION PARAMETERS
        ### -----------------
        self.niter = 1 # Number of Monte Carlo trials
        self.seed  = 2021  # Choose ‘random-seed’ to randomize

        # If True Statistical fluctuations are enabled,
        # If False niter forced to one
        self.do_fluctuate = False

        # When True, the simulation is stopped if none of the first
        # trials in the limit of 1 - det_level have reached the minimal
        # significance (3 sigma).
        self.do_accelerate = True

        # Observation position in the time slice
        self.obs_point = "end"

        ### -----------------
        ### DETECTION PARAMETERS
        ### -----------------

        # Fraction of trials above significance threshold to declare a
        # detection at that threshold
        self.det_level = 0.9

        # Number of on over off regions for detection
        self.alpha = 0.2

        self.array_North  = "FullArray"  # IRF subarrays in North
        self.array_South  = "FullArray"  # IRF subarrays in South
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

        # If True, use above value. If False, latency generated from Swift data
        self.fixswift = True

        # SWIFT Latency distribution data file in input base folder
        self.swiftfile = "swift/Swift_delay_times.txt"

        # Gammapy datasets folder
        self.extra_dir = "gammapy-extra-master/datasets"

        ### -----------------
        ### DEBUGGING / BOOKKEEPING
        ### -----------------

        # From 0 to 3, increasingly talkative
        # if negative or zero : no plot (show=0)
        # 0: evt counting, 1: + some results, 2: + event details
        # 3: details for each trials
        self.dbg        = 1
        self.show       = 1
        self.logfile    = "analysis.log" # log file

        self.save_simu  = False  # If True, Simulation saved to file
        self.save_grb   = False  # If True, GRB class saved to disk
        self.save_vis   = False  # If True, svae computed visibility
        self.save_fig   = False  # If True, plots saved to pdf file
        self.remove_tar = False  # If True, remove tarred output files

        ### -----------------
        ### EXPERTS/DEVELOPPERS ONLY
        ### -----------------
        self.test_prompt     = False # If True test prompt alone (experimental)

        # Prompt characteristics from the afterglow with same id.
        self.use_afterglow   = False

        # Shitf in time - Useful for tests
        self.tshift          = 0.    # float (days) 
        
        # Force fixed zenith
        self.fixed_zenith    = None  # If a value ("20 deg"), freezes zenith
        
        # Multiplicative factor of the input flux
        self.magnify         = 1 

        # If True, nothing on screen (output to log (if dbg =0))
        self.silent          = True

        # Store detailed information on slices if True
        self.write_slices    = False

        self.n_night         = None # Limit data to a maximal number of nights
        self.Emax            = None # Limit data energy bins to Emax
        
        # Bulld the source list to be processed
        self.srclist = self.source_ids()

    ###------------------------------------------------------------------------
    @classmethod
    def command_line(cls, debug=False):
        """
        Update an instance from the command line parameters        

        Parameters
        ----------

        debug : Boolean, optional
            If True, talk a bit. The default is False.

        Returns
        -------
        Configuration object
            Current instance.

        """
    
        cls = Configuration() # Initialize default
    
        # Define command line arguments - default values will come from
        # the configuraion file and are not known at this stage
        # When using getopt, it was possible to treat the reading
        # the configuraion file separately.
        # As a consequence, teh default values from the __init__ constructor
        # are not used here.
        parser = argparse.ArgumentParser(description="SoHAPPy", epilog="---")
        
        parser.add_argument('-f', '--first',  
                            help ="First source id",
                            type = int,
                            default = None)  
                
        parser.add_argument('-N', '--nsrc',   
                            help ="Number of source files",
                            type = int)  

        parser.add_argument('-n', '--niter',   
                            help ="Number of Monte Carlo iteration",
                            type = int)  
        
        parser.add_argument('-o', '--output',
                            help ="Output base folder (path)")
        
        parser.add_argument('-i', '--input',
                            help ="Input base folder (path)")  

        parser.add_argument('-c', '--config',
                            help ="Configuration file name")  
        
        parser.add_argument('-V', '--visibility',
                            help ="Visibility keyword")  
        
        parser.add_argument('-d', '--debug',
                            help ="Debugging flag", type = bool)  
        
        args, extra_args = parser.parse_known_args()
        
        # Find a configuration file, load the data    
        # This replaces the constructor defaults
        cls.find_and_read(args.config)
        
        # Supersede parameters if given
        if args.first is not None:
            cls.ifirst     = args.first
        if args.nsrc is not None:
            cls.nsrc       = args.nsrc
        if args.niter is not None:
            cls.niter      = args.niter
        if args.output is not None:
            cls.resfolder  = args.output
        if args.input is not None:
            cls.infolder   = args.input
        if args.visibility is not None:
            cls.visibility = args.visibility
        if args.debug is not None:
            cls.dbg        = args.debug
        
        # Show plots or not
        cls.show = 0 if cls.dbg  < 0 else abs(cls.dbg)

        # If debugging is requested, cannot be silent
        if cls.dbg>0: cls.silent = False

        # If the simulation is saved, it is not fluctuated
        if cls.save_simu:  cls.do_fluctuate = False

        # If no fluctuation, one simulation is enough !
        if cls.do_fluctuate == False:  cls.niter = 1
        
        # Bulld the source list to be processed
        cls.srclist = cls.source_ids()

        return cls
    
    ###------------------------------------------------------------------------
    def find_and_read(self, testname, debug=True):
        """
        Check existence of given file name, and read data if it exists.

        Parameters
        ----------
        testname : None, String or Path
            A tentative configuration file name.
        debug : boolean, optional
            If True, talk a bit. The default is True.

        Returns
        -------
        None.

        """
        
        self.filename = None # Assume not given
                
        if testname is None: # No argument was given on the command line
            basefolder   = Path(__file__).parent  # Where the code is
            self.filename = Path(basefolder,Configuration.def_conf)
            if debug: print(" Using default config file")
            
        else: # A name was given
            if  Path(testname).is_file(): # The argument is an existing file
                self.filename =  Path(testname)
            else:
                sys.exit(f"{__name__:}.py/build: "\
                         f"{testname:} configuration file does not exist")
                    
        ### Read configuration file from file name
        self.read_from_yaml()     
          
    ###------------------------------------------------------------------------
    def read_from_yaml(self, filename=None):
        """
        Read configuration file from disk.

        Parameters
        ----------
        filename : string, optional
            Configuration file name. The default is None.

        Returns
        -------
        None.

        """
        if filename is not None: 
            self.filename = Path(filename)
            
        #---------------------------------------------------
        def obj_dic(d):

            # top = type('new', (object,), d)
            seqs = tuple, list, set, frozenset
            for k, v in d.items():

                # Do not read wrong/deprecated keywords
                if k not in self.__dict__.keys():
                    print(f"{__name__:}.py: warning: {k:} ignored")
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

    #---------------------------------------------------
        print(">>> Read configuration from ",self.filename)

        try:
            f = open(self.filename, "r")
        except IOError:
            sys.exit("{}.py : {} does not exist."
                     .format(__name__,self.filename))

        data = yaml.load(f, Loader=SafeLoader)
        obj_dic(data)
        
        # Bulld the source list to be processed
        self.srclist = self.source_ids()

    ###------------------------------------------------------------------------
    def write(self, out_name = None):
        """
        Write current configuration to a yaml file for further use. 
        Note that this does not copy the possible comments of the original 
        file.

        Parameters
        ----------
        filename : String, optional
            Configuration file name. The default is None.

        Returns
        -------
        None.

        """
        
        if out_name is None: out_name = Path("config_backup.yaml")

        # Create a dict with unnecessary keywords dropped and no quantities
        newdict = {}
        for k,v in self.__dict__.items():

            if k in self.ignore:
                continue
            if isinstance(v,u.Quantity) or isinstance(v, Path):
                newdict[k] = str(v)
            else:
                newdict[k] = v

        f = open(out_name, 'w')
        yaml.dump(newdict, f, sort_keys=False)

    ###------------------------------------------------------------------------
    def print(self,log):
        
        """
        Print confuguration class instance to screen and log file

        Parameters
        ----------
        log : Log class instance
            Logbook output.

        Returns
        -------
        None.

        """

        log.prt(" Configuration file*     : {} ".format(self.filename))

        #----------------------------------------------------
        def title(s):
            log.prt("")
            log.prt("{:<27s} {}".format("=== "+s,40*"="))
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

        if self.prompt_dir is not None:
            if Path(self.infolder,self.prompt_dir).is_dir():
                log.prt(" Prompt data subfolder      : {}"
                        .format(self.prompt_dir))
            else:
                sys.exit("Prompt data : {} not a valid subfolder"
                         .format(self.prompt_dir))
        else:
            log.prt(" Prompt data                : not considered")
        log.prt(" Output population          : {}".format(self.datafile))
        log.prt(" Source file prefix         : {}".format(self.prefix))
        log.prt("             suffix         : {}".format(self.suffix))

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
        log.prt(" Save computed visibility   : {}".format(self.save_vis))
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
        if self.fixed_zenith is not None:
            log.warning("Zenith angle requested to be fixed at value '{}'"
                    .format(self.fixed_zenith))
        if self.magnify !=1:
            log.warning("GRB flux values are multiplied by {}"
                    .format(self.magnify))        
        if self.tshift !=0:
            log.warning("Times shifted by {} days"
                        .format(self.tshift))
        if not self.silent:
            log.warning("Maximise screen output")

        if self.write_slices == True:
            log.highlight("{:60s}"
            .format("Slice information saved to disk (write_slices=True)"))

        if self.n_night is not None:
             log.warning("GRB data limited to the first {} nights"
                    .format(self.n_night))
        if self.Emax is not None:
             log.warning("GRB energy bins limited to {}"
                    .format(self.Emax))

        log.prt("    --end\n")

    ###------------------------------------------------------------------------
    def create_output_folder(self):
        """
        Create the code outptut folder.
        The convention is that the output folder refers to the folder name
        containing the lightcurves, a subfolder refering to the visibility 
        which will be used.
        
        Parameters
        ----------
        log : Log instance
            pointer to the logbook file.

        Returns
        -------
        res_dir : pathlib Path
            output folder.

        """

        if len(self.srclist) > 1:
            ext = "_" + str(min(self.srclist)) + "_"+str(max(self.srclist))
        else:
            ext = "_" + str(self.srclist[0])
        
        res_dir = Path(self.resfolder, 
                       Path(self.data_dir).name, 
                       self.visibility,
                       self.out_dir,
                       self.visibility+ext).resolve()

        # Check that the output folder exists, otherwise try to create it
        if not res_dir.exists():
            warning("Creating {}".format(res_dir))
            try:
                Path.mkdir(res_dir,parents=True)
            except:
                sys.exit(f"{__name__:}.py: Could not create {res_dir:}")
        else:
            warning(f" Already exists :{res_dir:}")

        return res_dir

    ###------------------------------------------------------------------------
    def source_ids(self):
        """
        Obtain the source list from the input parameters.
        `first`is either an integer (identifier of the first source to be 
        analysed) or a list with integers (identifiers) or string (source 
        names).

        Returns
        -------
        A list of source identifiers to be analysed.

        """
        if type(self.ifirst) != list:
            if isinstance(self.ifirst, str):
                srclist = [self.ifirst]
            elif isinstance(self.ifirst, int):
                srclist = list(range(self.ifirst,self.ifirst+self.nsrc))
        else:
            srclist = self.ifirst

        return srclist

    ###------------------------------------------------------------------------
    def decode_visibility_keyword(self, folder=None, debug=False):

        """
        Decode the `visibility` keyword value and return an information
        to be handled at running time.

        The visibility keyword can be:
            
        * a visibility folder that should contain a suitable visibility json 
        file where all visibility instances of a data subset is stored.
        By definition the folder structure is the following:
        - output
        - population name
        - a subfolder name as described in skygen, i.e
            keyword_year1_Nyeasr_version

        * 'force' : visibility forced to be maximal (infinite nights)
        
        * 'built-in' : will use the default keyword stored in the source
        files, if any (for backward compatibility).
        
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
        
        if self.visibility in ["built-in","forced"]:
            print(" Visibilities from keyword: '{}' ".format(self.visibility))
            return self.visibility
                    
        else: # Seems to be a keyword to be found in visibility.yaml    
            from visibility import params_from_key
            vispar =  params_from_key(self.visibility)
            return vispar   

        # else: # A folder is given
        
        #     test = Path(folder.parent.parent, "visibility")
        
        #     if not test.is_dir():
        #         sys.exit(f"{__name__}.py : Keyword supposed to be a folder and is not")               
               
        #     # Find file name pattern from source identifiers
        #     idmin = min(self.srclist)
        #     idmax = max(self.srclist)
           
        #     # Complete version that find at least a valid file for the source 
        #     # id range
            
        #     candidates = list(test.glob("*.json"))
        #     visfile = None
        #     for fname in candidates:
        #         # Get ids in the filename
        #         bnds = [int(s) for s in fname.stem.split("_") if s.isdigit()]
                
        #         if len(bnds) == 4: # Range
        #             bndmin = bnds[-2]
        #             bndmax = bnds[-1]
        #         elif len(bnds) == 3: # One source only
        #             bndmin = bndmax = bnds[-1]
        #         # print(fname.name,": ",bnds,bndmin, bndmax )
                
        #         if bndmin<= idmin and idmax<= bndmax:
        #             visfile = fname
        #             break
            
        #     # Load data from the visibility file
        #     if visfile is not None:
        #         print(" Visibility data :",visfile)                     
        #         with open(visfile,"r") as f:
        #             import json
        #             return json.load(f)       
        #     else:
        #         sys.exit(f"{__name__}.py : Visibility: Recomputation seems required")               

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
    A standalone function to read a configuration froma configuration file
    with wome parameters possibly ovsersed by the command line.
    * Directly from the fileon disk, e.g?
    
    """
    
    from niceprint import Log
    log = Log()
    
    testfile = Path("myConfigs/config-LongFinalTest-omega.yaml")
    
    # Read the confi
    
    # testfile = Path("toto.yaml")
    # testfile = Path("D:\CTA\SoHAPPy\output\long\omega_prod3\strictmoonveto\pop_vis24_strictmoonveto-100iter-noacc/pop_vis24_strictmoonveto-100iter-noacc_20220316_183213.tar.gz")
    # sys.argv = ["", "-f", "343", "-N", "10","-n", "100", 
    #             "-V", "nomoonveto"]
    cf1 = Configuration() # Default
    cf1.read_from_yaml(filename=testfile)
    res_dir  = cf1.create_output_folder()  # Create output folder
    print(" >>> ",res_dir)
    cf1.decode_visibility_keyword(res_dir)
    cf1.print(log)



    # cf1.read_from_yaml(filename="config.yaml")
    # # cf1.print(log)
    # cf1.write(out_name = "config_bck.yaml")

    # cf2 = Configuration()
    # cf2.read_from_yaml(filename="config_bck.yaml")
    # cf2.print(log)

    # cf4 = Configuration.build(sys.argv[1:])
    # cf4.print(log)

    # test automatic output file naming
    
    # data_sub = Path(cf1.data_dir).name
    # if cf1.trigger ==0:
    #     trigpos_sub = "TP_default"
    # else:
    #     trigpos_sub = "TP_tobedefined"
        
    # print(" output subfolder for this run :")
    # out_sub = Path(data_sub,trigpos_sub)
    # print(Path(out_sub))
    # print(Path(cf1.resfolder,out_sub))
    