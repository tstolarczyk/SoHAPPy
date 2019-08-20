res_folder = "Result" # Folder for results
res_folder = "Today" # Folder for results

#------------------------------------------------------------------------------
### Code Initialisations
#------------------------------------------------------------------------------
dbg_level = 0 # 0 : event counting, 1: event counting + some results, 2: + details for each event, 3: + plots 
showplots = 2 # 0 : not shown, not written, 1: shown on screen, 2 : shown on screen and written out

#------------------------------------------------------------------------------
### Simulation Initialisations
#------------------------------------------------------------------------------
niter = 10 # Number of iterations for the Monte Carlo
alpha = 0.2 # 5 zones (1/0.2) to define the off region for background estimation
#------------------------------------------------------------------------------
### Physics Initialisations
#------------------------------------------------------------------------------
det_level       = 0.9 # Will be detected if n=3,5 sigma reached in det_level of the population
#redfactor       = [1, 2, 3, 5, 10, 20, 50, 100] # Initial flux is divided by this value, for test
#redfactor       = [0.5,1, 2] # Initial flux is divided by this value, for test
redfactor = [1]

grb_repository  = 'data/long_grb_test/' # to be changed into a repository to be scanned
#grblist         = ["LGRB1","LGRB2","LGRB3","LGRB4","LGRB5","LGRB6","LGRB7","LGRB8","LGRB9","LGRB10"]
grblist         = ["LGRB1"]

EBLmodel        = "dominguez"

#------------------------------------------------------------------------------
### Instrument response functions
#------------------------------------------------------------------------------
irf_repository = "irf"
chain       = ["Reco1"] # or Reco2
hemisphere  = ["North"] # or South
#hemisphere  = ["South"] # or South
#zenith      = ["20","40","60"] # or 40 deg...
zenith      = ["20"] # or 40 deg...
azimuth     = ["average"]  # or N, or S
#duration    = ["100s","30m"] # or 30m, 05h 50h
duration    = ["100s"] # or 30m, 05h 50h
NSB         = False # if True try to find a high NSB file
