"""
This module contains the classes and tools to handle a GRB object, i.e.
a collection of lightcurve bins for each of which is given an energy spectrum,
and information on site location.

Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""

import sys, os
import numpy as np
from pathlib import Path

import astropy
import astropy.units         as u
from   astropy.table         import Table, QTable
from   astropy.io            import fits
from   astropy.time          import Time
from   astropy.coordinates   import SkyCoord
from   astropy.coordinates   import AltAz

from visibility import Visibility

from utilities import warning, failure

from gammapy.modeling.models import EBLAbsorptionNormSpectralModel     
from gammapy.modeling.models import TemplateSpectralModel
from gammapy.modeling.models import PointSpatialModel

# Transform warnings into errors - useful to find who is guilty !
import warnings
#warnings.filterwarnings('error')
warnings.filterwarnings('ignore')

__all__ = ["GammaRayBurst","GRBDummy"]

###############################################################################
class GRBDummy():
    """
    A list of obaservation times and an integrated flux at each point
    Useful for slice merging tests or prototype for e.g. Fermi extrapolated
    spectra.

    """

    def __init__(self,tmin = 6*u.h, tmax = 54*u.h, nbin = 7, beta= 1):
        self.tmeas = tmeas = np.linspace(tmin,tmax,7)
        self.flux  = u.Unit("1/(GeV cm2 s)")*(tmeas/(1*u.h))**(-beta)

    def __str__(self):
        with np.printoptions(precision=3, suppress=False):
            print("Points measured at  : ",self.tmeas)
            print("Flux                : ",self.flux)
        return ""

    def plot(self,ax):
        ax.plot(self.tmeas,self.flux,marker="*")

###############################################################################
class GammaRayBurst(object):
    """
    Class to store GRB properties

    The GRB information is read from files.
    A GRB is composed of several parameters among which, a name, a redshift,
    a position in the sky, measurement points at which
    an energy spectrum is given and in some cases a default visibility window.
    Each energy spectrum is given as a series of flux measurements but is then
    stored as an interpolated spectrum (i.e. the energy bin width has a modest
    importance).
    The GRB properties like z, Eiso are generated in the population synthesis
    (G. Ghirlanda) and they connect the afterglow and prompt emission together.

    """

    ###------------------------------------------------------------------------
    def __init__(self):
        """
        This initializes a default GRB.

        Returns
        -------
        None.

        """
        self.name = 'dummy'
        self.z    = 0

        # GRB properties - Dummy default values
        self.radec    = SkyCoord(100*u.degree,-15*u.degree, frame='icrs')
        self.Eiso     = 0.*u.erg
        self.Epeak    = 0.*u.keV
        self.t90      = 0.*u.s
        self.G0H      = 0.*u.dimensionless_unscaled
        self.G0W      = 0.*u.dimensionless_unscaled
        self.FluxPeak = 0.*u.erg
        self.gamma_le = 0.*u.dimensionless_unscaled
        self.gamma_he = 0.*u.dimensionless_unscaled

        # GRB detection positions
        self.site_keys = ["North","South"] # Put it somewhere else !
        self.site      = {"North": 'Roque de los Muchachos',
                          "South": 'Paranal'}
  
        # GRB alert received
        self.t_trig   = Time('2000-01-01 02:00:00', scale='utc')                

        # Afterglow Flux table - Flux at a series of points
        self.Eval           = [0]*u.GeV
        self.tval           = [0]*u.s
        self.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")
        self.spec_afterglow = [] # Interpolated, non attenuated

        # Prompt energy spectrum
        self.prompt      = False # Default : no prompt
        self.id90        = -1 
        self.E_prompt    = [0]*u.GeV
        self.flux_prompt = [0]*u.Unit("1 / (cm2 GeV s)")
        self.spec_prompt = None # One Interpolated, non attenuated E-spectrum

        # Visibility (requires GRB points interval)
        self.vis  = { "North": Visibility(self,"North"),
                      "South": Visibility(self,"South")
                                          }
        # GRB cumulated spectra and spatial model
        self.spectra = [] # Gammapy models (one per t slice)
        self.spatial = PointSpatialModel(lon_0=0*u.deg,lat_0=0*u.deg)
        
        return
   
    ###------------------------------------------------------------------------   
    def EBLabsorbed(self, tab, model, debug=False):
        """
        Return the EBL-abosbrbed model.

        Parameters
        ----------
        tab :  TemplateSpectralModel
            A flux versus energy as an interpolated table.
        model : String
            An EBL model name among those available.
        debug : Boolean, optional
            If True let's talk a bit. The default is False.

        Returns
        -------
        attflux : TemplateSpectralModel
            An attenuated flux versus energy as an interpolated table.

        """
        
        attflux = tab
        if (model == "gilmore"):
            from ebl import EBL_from_file
            eblabs = EBL_from_file("EBL/data/ebl_gilmore12-10GeV.dat")
            attflux.values = tab.values*eblabs(self.Eval,self.z)
            if debug: print(" EBL : Gilmore")
        elif model != None and model != 'in-file':
            eblabs = EBLAbsorptionNormSpectralModel.read_builtin(model, 
                                                              redshift=self.z)
            attflux = tab*eblabs
            if debug: print(" EBL : ",model)
           
        return attflux
        
    ###------------------------------------------------------------------------
    @classmethod
    def from_fits(cls, filename, 
                  vis     = None, 
                  prompt  = False, 
                  ebl = None, 
                  t0 = Time('1800-01-01T00:00:00', format='isot', scale='utc'),
                  magnify = 1, 
                  forced_visible = False, 
                  dbg = 0):

        """
        Fluxes are given for a series of (t,E) values
        So far no flux exists beyond the last point.

        The spectrum is stored as a table, TemplateSpectralModel, that takes as 
        an input a series of flux values as a function of the energy.
        The model will return values interpolated in log-space with
        math::`scipy.interpolate.interp1d`, returning zero for energies
        outside of the limits of the provided energy array.
        
        The trigger time, time of the GRB explosion, is given either as an 
        absolute date or is arbitrary. In that case, an arbitrary start date is
        added to the elapsed times.
        
        The spectra is modifief for an EBL absorption.
        
        The afterglow flux can be adjusted by a multiplicatice facot for 
        various tests.
        
        If requested, the associated prompt spectrum is read. 

        Note that TemplateSpectralModel makes an interpolation.
        Following a question on the Slack gammapy channel on
        November 27th, and the answer by Axel Donath:
        The following statement later in the code gave an error
        (dlist_onoff is a collection of Dataset)
        dlist_onoff.write(datapath,prefix="cls",overwrite=True)
        gives:
        ...\gammapy\modeling\models\spectral.py", line 989, in to_dict
        "data": self.energy.data.tolist(),
        NotImplementedError: memoryview: unsupported format >f
        This error comes from the fact that the energy list as to be
        explicitely passed as a float as done below:
            cls.Eval.astype(float)
        (A Quantity is passed as requested but the underlying numpy
        dtype is not supported by energy.data.tolist())


        Parameters
        ----------
        cls : GammaRayBurst class
            GammaRayBurst class instance
        filename : Path
            GRB input file name
        prompt: Boolean
            If True, read information from the associated prompt component
        ebl : STRING, optional
            The EBL absorption model considered. If ebl is "built-in", uses
            the absorbed specrum from the data.
        t0: Astropy Time
            A start date to be added to the slice elapsed time when no GRB 
            explosion date is given.
        magnify : float, optional
            Flux multiplicative factor for tests to the afterglow model. 
            Default is 1
        forced_visible: Boolean
            If True, the GRB is always visible, in particular the obervation
            is always during night. Default is False.
        dbg: Integer
            A debugging level. Default is zero, no debugging messages.
            
        Returns
        -------
        A GammaRayBurst instance.

        """
        
        # Check if file was compressed and/or if it exists
        # Strangely path.is_file() does not give False but an error
        if not os.path.exists(filename): # Current file does not exist
            if filename.suffix == ".gz": # If this is a gz file, try the non gz file
                filename = filename.with_suffix("")
            else: # If this not a gz file, try the gz file
                filename = filename.with_suffix(filename.suffix+".gz")
        if not os.path.exists(filename): # Check that the new file exist   
            sys.exit("grp.py: {} not found".format(filename))    

        # Open files and get header and data
        hdul = fits.open(filename)
        hdr  = hdul[0].header
        
        keys_0 = list(hdul[0].header.keys()) # Keys in header
       
        cls = GammaRayBurst() # constructor
        
        cls.z    = hdr['Z']
        # Remove all extensions
        cls.name = str(filename.name).rstrip(''.join(filename.suffixes)) 
        
        cls.radec    = SkyCoord(hdr['RA']*u.degree, hdr['DEC']*u.degree,
                                frame='icrs')
        cls.Eiso     = hdr['EISO']*u.erg
        cls.Epeak    = hdr['EPEAK']*u.keV
        cls.t90      = hdr['Duration']*u.s
        cls.G0H      = hdr['G0H']
        if "G0W" in keys_0: # Not in SHORTFITS
            cls.G0W      = hdr['G0W']
        cls.Fluxpeak = hdr['PHFLUX']*u.Unit("ph.cm-2.s-1")
        cls.gamma_le = hdr['LOWSP']
        cls.gamma_he = hdr['HIGHSP']

        # GRB trigger time
        if "GRBJD" in keys_0: # Not in SHORTFITS
            cls.t_trig = Time(hdr['GRBJD']*u.day,format="jd",scale="utc")
        elif "GRBTIME" in keys_0: # Add start date
            cls.t_trig = Time(hdr['GRBTIME'] + t0.jd,format="jd",scale="utc")           

        ###--------------------------
        ### Time slices - so far common to afterglow and prompt
        ###--------------------------
        tval  = QTable.read(hdul["TIMES (AFTERGLOW)"])
        
        # Temporary : in shot GRB fitrs file, unit is omitted
        if isinstance(tval[0][0], astropy.units.quantity.Quantity):
            cls.tval = tval[tval.colnames[0]]
        else: # In SHORT GRB, the time bin is given, [t1, t2].
            cls.tval = tval["col1"]*u.s
            
        ###--------------------------
        ### Visibilities --- includes tval span
        ###--------------------------
        # Either read the default visibility from the GRB fits file (None)
        # or recompute it(a keyword has been given and a dictionnary retrieved)
        # or read it from the specified folder
        for loc in ["North","South"]:
            if vis == None:
                cls.vis[loc] = cls.vis[loc].from_fits(cls, hdr, hdul,
                                                           hdu=1,loc=loc)
            elif isinstance(vis,dict):
                cls.vis[loc] = Visibility.compute(cls,
                                                  loc,
                                                  param     = vis,
                                                  force_vis = forced_visible,
                                                  debug     = bool(dbg>2))                
            else:
                name = Path(vis,cls.name+"_"+loc+"_vis.bin")
                cls.vis[loc] = Visibility.read(name)       
 
        ###--------------------------
        ### Afterglow
        ###--------------------------
        
        ### Read and store Energies, times and fluxes ---
        k = "Energies (afterglow)"
        cls.Eval  = Table.read(hdul[k])[k].quantity  
            
        ### Get flux,possibly already absorbed
        if (ebl == "in-file"): # Read default absorbed model
            flux = QTable.read(hdul["EBL-ABS. SPECTRA (AFTERGLOW)"])
        else:
            flux = QTable.read(hdul["SPECTRA (AFTERGLOW)"])
                
        ### Store the flux
        flux_unit    = u.Unit(flux.meta["UNITS"])/u.Unit("ph") # Removes ph
        # flux_unit = u.Unit("1/(cm2 GeV s)")
        icol_t    = len(flux.colnames)          # column number - time
        jrow_E     = len(flux[flux.colnames[0]]) # row number

        # Note the transposition from flux to fluxval
        cls.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                cls.fluxval[i][j] = magnify* flux[j][i]*flux_unit # transp!
                
        # Build time series of interpolated spectra        
        for i,t in enumerate(cls.tval):
            glow = TemplateSpectralModel(energy = cls.Eval.astype(float),
                                         values = cls.fluxval[i],
                                         interp_kwargs={"values_scale": "log"})                 
            cls.spec_afterglow.append(glow)
            
        ###--------------------------
        ### Prompt - a unique energy spectrum 
        ###--------------------------
        ### Get the prompt if potentially visible and if requested

        cls.prompt = False # No prompt component was found
        if prompt:
            if cls.vis["North"].vis_prompt or cls.vis["South"].vis_prompt:
                # Deduce prompt folder from GRB path name
                folder = Path(filename.absolute().parents[0], "../prompt/")
                cls.spec_prompt = cls.get_prompt(grb_id = cls.name[5:],
                                                 folder = folder,
                                                 debug  = bool(dbg))
                if cls.spec_prompt != None: 
                    cls.prompt = True
                
        ###--------------------------
        ### Total attenuated spectra
        ###--------------------------
        ### Build the total attenuated model
        for i,t in enumerate(cls.tval):
            spec_tot = cls.spec_afterglow[i] 
            if cls.prompt:
                if i<=cls.id90: spec_tot += cls.spec_prompt[i]
            model = cls.EBLabsorbed(spec_tot,ebl)
            cls.spectra.append(model)
 
        hdul.close()

        return cls
    
    ###------------------------------------------------------------------------
    @classmethod
    def from_yaml(cls, data, ebl= None,magnify=1):

        """

        Returns
        -------
        A GammaRayBurst instance.

        """

        cls = GammaRayBurst() # This calls the constructor
        
        cls.z        = data["z"]
        cls.name     = data["name"]
        cls.radec    = SkyCoord(data["ra"], data["dec"], frame='icrs')
        cls.Eiso     = u.Quantity(data["Eiso"])
        cls.Epeak    = u.Quantity(data['Epeak'])
        cls.t90      = u.Quantity(data['t90'])
        cls.G0H      = data['G0H']
        cls.G0W      = data['G0W']
        cls.Fluxpeak = u.Quantity(data['Fluxpeak'])
        cls.gamma_le = data['gamma_le']
        cls.gamma_he = data['gamma_he']
        cls.t_trig = Time(data["t_trig"], format="datetime",scale="utc")

        ### Energies, times and fluxes ---
        Emin     = u.Quantity(data["Emin"])
        Emax     = u.Quantity(data["Emax"])
        tmin     = u.Quantity(data["tmin"])
        tmax     = u.Quantity(data["tmax"])
        ntbin    =  data["ntbin"]
        cls.Eval = np.asarray([Emin.value, Emax.to(Emin.unit).value])*Emin.unit

        if ntbin != 1:
            cls.tval = np.logspace(np.log10(tmin.to(u.s).value),
                                   np.log10(tmax.to(u.s).value),
                                   ntbin)*u.s
        else: # A single wtime window
            cls.tval = np.array([tmin.value,tmax.to(tmin.unit).value])*tmin.unit

        flux_unit = u.Unit("1/(cm2 GeV s)")
        cls.fluxval = np.zeros( (len(cls.tval),len(cls.Eval)) )*flux_unit

        for i,t in enumerate(cls.tval):
            for j,E in enumerate(cls.Eval):
                dnde = (u.Quantity(data["K"])*(E/data["E0"])**-data["gamma"]
                                 *(t/data["t0"])**-data["beta"])
                #print(i,j,dnde)
                cls.fluxval[i][j] = magnify* dnde.to(flux_unit)

        ### Visibilities --- includes tval span
        for loc in ["North","South"]:
            cls.vis[loc]  = Visibility(cls,loc)
            # Recomputing done in the main
            # cls.vis[loc].compute(debug=False)
            
        ### No prompt component foreseen in this case
        cls.prompt = False

        for i,t in enumerate(cls.tval):
            # Note that TableModel makes an interpolation
            # Following a question on the Slack gammapy channel on
            # November 27th, and the answer by Axel Donath:
            # The foloowing statement later in the code gave an error
            # (dlist_onoff is a collection of Dataset)
            # dlist_onoff.write(datapath,prefix="cls",overwrite=True)
            # gives:
            # ...\gammapy\modeling\models\spectral.py", line 989, in to_dict
            # "data": self.energy.data.tolist(),
            # NotImplementedError: memoryview: unsupported format >f
            # This error comes from the fact that the energy list as to be
            # explicitely passed as a float as done below:
            #    cls.Eval.astype(float)
            # (A Quantity is passed as requested but the underlying numpy
            # dtype is not supported by energy.data.tolist()
            tab = TemplateSpectralModel(energy = cls.Eval.astype(float),
                                          values = cls.fluxval[i],
                                          interp_kwargs={"values_scale": "log"})  
                
            model = cls.EBLabsorbed(tab, ebl)
            cls.spectra.append(model)

        return cls

    ###------------------------------------------------------------------------
    def get_prompt(self, grb_id = None, folder = None,
                   subfolder= "ctagrbs_spop", debug = False):
        """
        Get the prompt component associated to the afterglow.
        The prompt spectra are produced so that they correspond to the values 
        provided by Giancarlo (Lorentz Factor, peak energy, photon flux, and 
        redshift). 
        A single set of spectra is given spanning over the T90 GRB duration. 
        The flux is an "average" prompt spectrum over T90. 

        Parameters
        ----------
        grb_id : integer, optional
            The GRB identifier. The default is None.
        folder : String, optional
            Where to find the prompt data, below the 'lightcurves/prompt' 
            folder. The default is None.
        debug : Boolean, optional
            If True, let's talk a bit. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        
        ###----------------------------
        ### Check everything goes well        
        ###----------------------------
        
       # Define the prompt data folder from the afterglow parent folder
        folder = Path(folder,subfolder)
        if not folder.exists(): sys.exit(" Wrong prompt data folder")

        if grb_id == None: sys.exit(" Provide a GRB identifier")
        
        # if grb_id in [6, 30 ,191]:
        #     failure(" GRB prompt",grb_id," simulated with a Wind profile")
        #     failure(" It cannot be associated to the current afterglow !")
        #     return -1, None
        
        ###----------------------------
        ### Mask times beyong t90         
        ###----------------------------
        
        # Find the closest time bin at t90 in the Afterglow
        # id90 is the index of the time following the t90 value
        # t90 is therefore between the index id90-1 and id90
        self.id90 = np.where(self.tval >= self.t90)[0][0]
        if self.id90 == 0: 
            warning("{} : t90 = {} is before the first afterglow time bin"
                    .format(self.name, self.t90))
            return None
        
        # Define a reduction for each time slice, either 0 or 1 except at t90.
        fraction = (self.t90 - self.tval[self.id90-1])
        fraction = fraction/(self.tval[self.id90]-self.tval[self.id90-1])
        self.weight = np.concatenate(  (np.ones(self.id90),  
                                        [fraction], 
                                        np.zeros(len(self.tval)-self.id90-1)))

        ###----------------------------
        ### Read prompt data    
        ###----------------------------       
        # Open file - read the lines
        filename = Path(folder,grb_id+"-spec.dat")
        if not filename.exists(): 
            failure("{} does no exist".format(filename))
            return None
        
        file  = open(filename,"r")
        lines = file.readlines()
        
        # get Gamma and z
        data         = lines[0].split()
        gamma_prompt = float(data[2])
        redshift     = float(data[5])
        
        # Consistency check - Zeljka data are rounded at 2 digits
        if np.abs(self.z - redshift)>0.01 or \
           np.abs(self.G0H - gamma_prompt)>0.01:
              failure(" {:10s}: Afterglow / prompt : z= {:4.2f} / {:4.2f}  G0H/gamma= {:5.2f} / {:5.2f} (G0W= {:5.2f})"
                    .format(self.name, self.z,redshift, self.G0H, gamma_prompt, self.G0W ))
        
        # Get units - omit "ph" in flux unit
        data     = lines[1].split()
        unit_e   = data[0][1:-1]
        unit_flx = ''.join([ x + " " for x in data[2:]])[1:-2]
        
        # Get flux points - assign units to data, use afterglow units
        data = lines
        Eval = []
        fluxval = []
        for line in lines[4:]:
            data = line.split()
            Eval.append(float(data[0]))
            fluxval.append(float(data[1]))
        
        self.E_prompt    = Eval*u.Unit(unit_e)
        self.flux_prompt = fluxval*u.Unit(unit_flx)
        
        if (unit_e != self.Eval.unit):
            if debug:
                warning("Converting energy units from {} to {}"
                        .format(unit_e,self.Eval[0].unit))
            self.E_prompt = self.E_prompt.to(self.Eval[0].unit)
            
        if (unit_flx != self.fluxval[0].unit):
            if debug:
                warning("Converting flux units from {} to {}"
                        .format(unit_flx, self.fluxval[0][0].unit))
            self.flux_prompt = self.flux_prompt.to(self.fluxval[0][0].unit) 
            
        ###----------------------------
        ### Create a list of weighted models for non-zero weights    
        ###----------------------------
        models = []
        for i in range(self.id90+1):
            flux_w = self.flux_prompt*self.weight[i]
            models.append(TemplateSpectralModel(energy = self.E_prompt,
                                        values = flux_w,
                                        interp_kwargs={"values_scale": "log"}))
            
        if debug:
            print("Prompt associated to ",self.name)
            print("Gamma = ",gamma_prompt," redshift =",redshift)
            print("  {:8} : {:10}".format(unit_e, unit_flx))
            for E, flux in zip(self.E_prompt,self.flux_prompt):
                print("  {:8.2f} : {:10.2e} "
                      .format(E.value,flux.value))  
            islice=0
            print(" {:>8} - {:>5} ".format("Time","Weight"),end="")
            for t, w in zip(self.tval, self.weight):
                print(" {:8.2f} - {:5.2f} ".format(t,w),end="")
                print("*") if islice== self.id90 else print("")
                islice+=1
                
            print("Prompt is between: ",self.tval[self.id90-1],self.tval[self.id90])
            
        return models
    
    ###------------------------------------------------------------------------
    @classmethod
    def read_prompt(cls, filename, glow= None, ebl= None,
                    z=0*u.dimensionless_unscaled, magnify=1):
        """
        This function is for tests using ttime-resolved spectra.
        Read prompt data from a file and associate it to the afterglow
        information if requested (or keep the default from the constructor
        otherwise)

        Parameters
        ----------
        cls : GammaRayBurst
            Class intance
        filename : String
            DESCRIPTION.
        glow : boolean, optional
            If defined, the Prompt characteritics are obatined from the
            afterglow with same id. The default is None.
        ebl : string, optional
            The name of an EBL absoprtion model in the Gammapy data.
            The default is None.
        z : float, optional
            GRB redshift. The default is 0*u.dimensionless_unscaled.
        magnify : float, optional
            Flux multiplicative factor for tests. Default is 1

        Returns
        -------
        grb : GammaRayBurst
            A GammaRayBurst instance

        """

        if (glow != None):
            grb = glow # Copy afterglow parameters
            #grb.z = 0
            # Flux table - Flux at a series of points
            grb.Eval           = [0]*u.GeV
            grb.tval           = [0]*u.s
            grb.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")
            grb.spectra = [] # Gammapy models (one per t slice)
        else:
            grb = cls() # This calls the constructor

        hdul = fits.open(filename)

        grb.name = Path(filename.name).stem

        # Energies, times and fluxes
        grb.Eval     = Table.read(hdul,hdu=1)["energy"].quantity*u.TeV
        grb.tval     = Table.read(hdul,hdu=2)["time"].quantity*u.s
        flux         = Table.read(hdul,hdu=3)
        flux_unit = u.Unit("1/(cm2 TeV s)")

        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Note the transposition from flux to fluxval
        grb.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                f = flux[j][i]
                if (f>1):
                    # print(i,j,f)
                    f =0 # Correct a bug in event #172 - to be removed
                grb.fluxval[i][j] = magnify*f*flux_unit # transp!

        for i,t in enumerate(grb.tval):
            # Note that TableModel makes an interpolation
            tab = TemplateSpectralModel(energy  = grb.Eval.astype(float),
                                         values = grb.fluxval[i],
                                         norm   = 1.,
                                         interp_kwargs={"values_scale": "log"})
            model = cls.EBLabosrbed(tab, ebl)
            grb.spectra.append(model)

        hdul.close()
        return grb
    ###------------------------------------------------------------------------
    def write(self, folder, debug=True):

        import pickle

        filename = Path(folder,self.name+".bin")
        outfile  = open(filename,"wb")
        pickle.dump(self,outfile)
        outfile.close()
        if (debug): print(" GRB saved to : {}".format(filename))

        return

    ###------------------------------------------------------------------------
    def altaz(self,loc="",dt=0*u.s):
        """
        Get altitude azimuth for the GRB at a given site at GRB time t (s)

        Parameters
        ----------
        location : string
            Either North or South
        dt : Quantity (time)
            The time elapsed since the trigger time

        Returns
        -------
        altaz : astropy.coordinates.AltAz
            The GRB poistion in the sky at the given time and site

        """
        if (type(dt) is not astropy.units.quantity.Quantity):
            sys.exit("Time has to be a quantity (weird behaviour otherwise")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            altaz = self.radec.transform_to(AltAz(obstime  = dt + self.t_trig,
                                                  location = self.vis[loc].site))

        return altaz

    ###------------------------------------------------------------------------
    def __str__(self):
        """ Printout the GRB properties """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ttrig = self.t_trig.datetime

        # After glow duration
        txt = '\n'
        txt += "==================================================================\n"
        txt += "===                             {:10s}                     ===\n"\
        .format(self.name)
        txt += "==================================================================\n"
        txt += '  RA, DEC       : {:6.2f} {:6.2f}\n'\
        .format(self.radec.ra, self.radec.dec)
        txt += '  Redshift      : {:6.2f}\n'.format(self.z)
        txt += '  Eiso          : {:6.2e}\n'.format(self.Eiso)
        txt += '  Epeak         : {:6.2f}\n'.format(self.Epeak)
        txt += '  t90           : {:6.2f}\n'.format(self.t90)
        txt += '  G0H / G0W     : {:6.2f} / {:6.2f}\n' \
        .format(self.G0H,self.G0W)
        txt += '  Flux peak     : {:6.2f}\n'.format(self.Fluxpeak)
        txt += '  gamma LE / HE : {:6.2f} / {:6.2f}\n' \
        .format(self.gamma_le,self.gamma_he)

        txt += '  t_trig        : {}\n'.format(ttrig)
        txt += '  Duration      : {:6.2f} {:10.2f}\n' \
            .format( (self.tval[-1]-self.tval[0]).to(u.d),
                     (self.tval[-1]-self.tval[0]).to(u.s))
        # txt += '  Bins : E, t, dt    : {:4d} {:4d} {:4d} \n' \
        # .format(len(self.Eval), len(self.tval), len(self.time_interval))
        txt += '  Bins : E, t   : {:4d} {:4d} \n' \
        .format(len(self.Eval), len(self.tval))
        if self.prompt:
            txt += ' Prompt component : \n'
            txt += '  Up to slice   : {:3d}\n'.format(self.id90)
            txt += "  Bins : E      : {:3d}\n".format(len(self.E_prompt))
        else: 
            txt += ' Prompt component not considered'
        return txt
    
###############################################################################
if __name__ == "__main__":

    """
    A standalone function to read a GRB and make various tests
    """
    os.environ['GAMMAPY_DATA'] =r'../input/gammapy-extra-master/datasets'
    # Transform warnings into errors - useful to find who is guilty !
    
    #warnings.filterwarnings('error')
    warnings.filterwarnings('ignore')
    
    from   utilities import Log
    import yaml
    from yaml.loader import SafeLoader
    import grb_plot as gplt


    ###------------------
    ### GRB to read
    ###------------------
    ifirst     = 1 #  ["190829A"]
    ngrb       = 1 # 250
    ebl        = "dominguez"
    if type(ifirst)!=list: grblist = list(range(ifirst, ifirst + ngrb))
    else: grblist = ifirst
    
    ###------------------
    ### Visibility
    ###------------------
    # visibility      : "D:/CTAA/SoHAPPy/input/visibility/short/vis_24_strictmoonveto"      
    visibility = "D:/CTAA/SoHAPPy/input/visibility/long/vis_24_strictmoonveto"      
    # visibility = "../input/visibility/vis_24_strictmoonveto"
    # visibility = None # Default in file if it exists
    visibility = "strictmoonveto"
    with open("visibility.yaml") as f:
        visibility  = yaml.load(f, Loader=SafeLoader)[visibility]
        
    
    prompt     = False
    save_grb   = False # (False) GRB saved to disk -> use grb.py main
    dbg        = 1
    # grb_folder = "../input/lightcurves/LONG_FITS"
    # grb_folder = "../input/lightcurves/SHORT_FITS"
    # grb_folder = "D:/CTAA/SoHAPPy/input/lightcurves/SHORT_FITS/"    
    grb_folder = "D:/CTAA/SoHAPPy/input/lightcurves/LONG_FITS/"
    res_dir    = "."    
    log_filename    = res_dir  + "/grb.log"
    log = Log(name = log_filename,talk=True)
    
    ## Test dummy GRB
    # grb = GammaRayBurst()
    # print(grb)
    # for loc in ["North", "South"]:
    #     print(grb.vis[loc].print(log=log))
    
    # Loop over GRB list
    for item in grblist:
        filename = Path(grb_folder,"Event"+str(item)+".fits")
        grb = GammaRayBurst.from_fits(filename,
                                  ebl     = ebl,
                                  prompt  = prompt, 
                                  vis     = visibility,
                                  dbg     = dbg)

        print(grb)
        gplt.time_spectra(grb)                
        gplt.energy_spectra(grb)                
        gplt.visibility_plot(grb, loc ="North")
        gplt.visibility_plot(grb, loc ="South")

        if (save_grb): grb.write(res_dir) # Save GRB if requested
