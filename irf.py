# -*- coding: utf-8 -*-

"""
Created on Thu Dec 12 09:01:41 2019

@author: Stolar

"""
import sys
import numpy as np
from pathlib import Path
import astropy.units as u
from astropy.visualization import quantity_support

from gammapy.maps import MapAxis
from gammapy.irf import load_cta_irfs

# Keyword lists and true values, also in log for time
zenith_list =  {"20deg": 20*u.deg,
                "40deg": 40*u.deg,
                "60deg": 60*u.deg }

dt_list     =  {"prod3":
                    {"100s" : (100*u.s),
                     "30m"  : (30*u.min).to(u.s),
                     "5h"   : (5*u.h).to(u.s),
                     "50h"  : (50*u.h).to(u.s) },
                 "prod5":
                     {"1800s"  : (1800*u.s),   # 30 min
                      "18000s" : (18000*u.s),  # 5h
                      "180000s": (180000*u.s)} # 50h
                }

dtl         =  {"prod3":
                    {"100s" : np.log10((100*u.s).value),
                     "30m"  : np.log10((30*u.min).to(u.s).value),
                     "5h"  : np.log10((5*u.h).to(u.s).value),
                     "50h"  : np.log10((50*u.h).to(u.s).value) },
                "prod5":
                    {"1800s"   : np.log10((1800*u.s).value),    # 30 min
                     "18000s"  : np.log10((18000*u.s).value),   # 5h
                     "180000s" : np.log10((180000*u.s).value) } # 50h
                }

# Validity range of IRF in zenith (see documentation of this module).
# The 60deg IRF is allowed to be used down to alitude zero for tests
# Its use is foreseen to be limited by the altmin variable
zenith_valid = {"20deg": [0*u.deg, 33*u.deg],
                "40deg": [33*u.deg, 54*u.deg],
                "60deg": [54*u.deg, 90*u.deg] # Should be limited by altmin
                }

# Validity range of IRF in time, taking into account that the validity
# intervals are somehow in logarithmic scale.
# The edge values are the following :
# 0, 424.26 s, 5692.01 s (94.9'), 56921.0 s (15.8h)

dt_log_valid = \
{"prod3":
    {"100s": [0,
             10**(0.5*( dtl["prod3"]["100s"] + dtl["prod3"]["30m"] )) ],
    "30m" : [10**(0.5*( dtl["prod3"]["100s"] + dtl["prod3"]["30m"] )),
             10**(0.5*( dtl["prod3"]["30m"]  + dtl["prod3"]["5h"] )) ],
    "5h"  : [10**(0.5*( dtl["prod3"]["30m"]  + dtl["prod3"]["5h"] )),
             10**(0.5*( dtl["prod3"]["5h"]  + dtl["prod3"]["50h"] )) ],
    "50h" : [10**(0.5*( dtl["prod3"]["5h"]  + dtl["prod3"]["50h"] )),
             np.Inf]
    },
 "prod5":
    {"1800s" : [0,
                10**(0.5*( dtl["prod5"]["1800s"]  + dtl["prod5"]["18000s"] ))],
    "18000s" : [10**(0.5*( dtl["prod5"]["1800s"]  + dtl["prod5"]["18000s"] )),
                10**(0.5*( dtl["prod5"]["18000s"] + dtl["prod5"]["180000s"] ))],
    "180000s": [10**(0.5*( dtl["prod5"]["18000s"] + dtl["prod5"]["18000s"])),
                 np.inf]
    }
}

# Minimal acceptable energies depending on the IRF
# (resp. generated and reconstructed energies)
#
# The masking later will remove the bin containing the E value.
# If the E value is an edge, the subsequent bin is lost.
# The minimal and maximal energies need therefore to be slighlty
# before, resp. after the first, resp last edge of interest.

egen_min = {"20deg": 12.6*u.GeV,
            "40deg": 26.4*u.GeV,
            "60deg": 105.3*u.GeV
            }

etrue_max = {"20deg": 17*u.TeV,
             "40deg": 17*u.TeV,
             "60deg": 17*u.TeV}

nbin_per_decade = 4


__all__=['IRF']

###############################################################################
class IRF():
    """
    This class handles the Instrument response function information and
    utilities
    """

    ###------------------------------------------------------------------------
    def __init__(self,filename  = None,
                      irf       = None,
                      subarray  = None,
                      etrue     = None,
                      kzen      = None,
                      kaz       = None,
                      kdt       = None):
        """
        Initialise the object.

        Parameters
        ----------
        name : String, optional
            The IRF file name on disk. The default is "none".
        irf : TYPE, optional
            DESCRIPTION. The default is None.
        array : String
            The array configuration. default is None
        etrue : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.


        """
        self.filename  = filename
        self.irf       = irf
        self.subarray  = subarray
        self.etrue     = etrue # idem
        self.kzen      = kzen
        self.kaz       = kaz
        self.kdt       = kdt
  
        return

    ###------------------------------------------------------------------------
    def find_best_keys(zenith, azimuth, obstime,
                       closest=False, fixed_zenith=None, prod=None):
        """
        Find the best keys to later identify the IRF data file.

        Parameters
        ----------
        zenith : Quantity, angle
            zenith angle.
        azimuth : Quantity, angle
            Azimuth angle.
        obstime : Quantity, time
            Observation duration.
        closest : Boolean, optional
            If True, obtain zenith and observation time from the closest
            available sampling, instead of using the predefined validity value.
            The default is False.
        fixed_zenith : Quantity angle, optional
            If defined, the zenith is fixed at the given value, wathever the
            observation (for tests). The default is None.
        prod: S

        Returns
        -------
        String, String, String
            The three strings defining the IRF file and/or folder.

        """
        ###--------------------------
        def find_closest_key(mydict,x):
            dict_values = np.asarray([x.value for k, x in mydict.items() ])
            closest_val = dict_values [(np.abs(dict_values - x.value)).argmin() ]
            closest_key = [k for k, x in mydict.items() if x.value == closest_val]
            return closest_key[0]
        ###--------------------------

        # Zenith

        if (closest):
            kzen = find_closest_key(zenith_list, zenith)
        else:
            found = False
            for k,v in zenith_valid.items():
                if (zenith >=v[0] and zenith < v[1] ):
                    kzen = k
                    found = True
                    continue
            if (not found):
                sys.exit("get_irf_file: zenith= {} => range not found"
                         .format(zenith))

        # Azimuth - implement here N, S choice
        kaz = "average"

        # Observation time
        if (closest):
            kdt =find_closest_key(dt_list[prod], obstime.to(u.s))
        else:
            found = False
            for k,v in dt_log_valid[prod].items():
                if obstime.to(u.s).value >= v[0] \
                   and obstime.to(u.s).value < v[1] :
                    kdt = k
                    found = True
                    continue
            if (not found):
                sys.exit("get_irf_file: obstime range not found")

        return kzen, kaz, kdt

    ###------------------------------------------------------------------------
    @classmethod
    def from_observation(cls,
                         zenith   = 0*u.deg,
                         azimuth  = 0*u.deg,
                         obstime  = 0*u.h,
                         subarray = None,
                         loc      = None,
                         nsb      = None,
                         irf_dir  = None,
                         closest  = False,
                         debug    = False):
        """
        Get the IRF data from the chaarcteristics of an observation.
        Note that MapAxis accepts ony a normalised list of axis type as
        described `here <https://docs.gammapy.org/dev/irf/index.html#irf-axis-naming>`_.

        Parameters
        ----------
        cls : IRF class instance
            Current object.
        zenith : Quantity angle, optional
            Observation zenith angle. The default is 0*u.deg.
        azimuth : Quantity angle, optional
            Observation azimuth angle. The default is 0*u.deg.
        obstime : Quantity time, optional
            Observation time. The default is 0*u.h.
        kchain : String, optional
            IRF version folder. The default is "prod3-v2".
        array : String, optional
            Sub-array folder name. The default is "FullArray".
        loc : String, optional
            Site location. The default is None.
        nsb : TYPE, optional
            Handle special NSB cases. The default is None.
        irf_dir : string, optional
            IRF folder. The default is None.
        closest : Boolean, optional
            Choose IRF closes to IRF range bounds. The default is False.
        debug : Boolean, optional
            If True, verbose mode. The default is False.

        Returns
        -------
        IRF ojbect
            Initialise an IRF object with the obtained values.

        """

        if (loc == None):
            sys.exit("from_observation : location should be defined")

        # Find best keys (IRF interpolation)
        idx = irf_dir.find("prod")
        prod = irf_dir[idx:idx+5]
        kzen, kaz, kdt = cls.find_best_keys(zenith, azimuth, obstime,
                                                closest = closest,
                                                prod=prod)
            
        # Find base folder from the keys
        if prod=="prod3":
            folder = Path(irf_dir,subarray,loc,kzen)
                    # Build the filename, get the IRF
            if (subarray != "FullArray"):
                kaz = kaz+"_"+subarray
            subfolder = loc+"_z"+kzen[:2]+"_"+kaz+"_"+kdt
            irf_file = Path(folder,subfolder,"irf_file.fits.gz")

        elif prod=="prod5":
            folder = Path(irf_dir,subarray)
            if kaz=="average":
                filename = "Prod5-"+loc+"-"+kzen+"-AverageAz-"\
                    +subarray+"."+kdt+"-v0.1.fits.gz"
            else:
                sys.exit("Azimuth, only average is implemented")
                
            irf_file = Path(folder,filename)
            
        else:
            sys;exit("Not implemented")

        irf = load_cta_irfs(irf_file)

        eirf_min = min(irf["aeff"].data.axes["energy_true"].edges)
                
        etrue_axis  = MapAxis.from_energy_bounds(eirf_min,
                                                 etrue_max[kzen],
                                                 nbin = nbin_per_decade,
                                                 per_decade=True,
                                                 name="energy_true")


        # Alternatively, this is not optimal for masking except if the
        # number of bins is very large:
        # erec_axis  = MapAxis.from_energy_bounds(min(erec_min.values()),
        #                                         max(erec_max.values()),
        #                                         nbin = nbin_per_decade,
        #                                         per_decade=True,
        #                                         name="Rec. energy")

        if (irf_file.exists() != True):
            sys.exit(" This file does not exist :",irf_file)

        return cls(filename  = irf_file,
                   irf       = irf,
                   subarray  = subarray,
                   etrue     = etrue_axis,
                   kzen      = kzen,
                   kaz       = kaz,
                   kdt       = kdt
                   )
    
    ###------------------------------------------------------------------------
    def print(self):
        """
        Print out some IRF class contents

        Returns
        -------
        None.

        """
        print("IRF             : ",self.filename)
        print(" Etrue          : ",self.etrue) 
        print("          edges : ",self.etrue.edges)
        print(" Sub-array      : ",self.subarray) 
        print(" Zenith         : ",self.kzen) 
        print(" Azimuth        : ",self.kaz) 
        print(" Duration       : ",self.fr) 
        return

###############################################################################
### Utilities and check plots
###############################################################################
###------------------------------------------------------------------------
import mcsim_config as mcf
def containment_plot(irf,
                     eunit="GeV",erec_min =10*u.GeV, erec_max = 100*u.TeV,
                     subarray=None,tag=None, ax=None):

    import matplotlib.pyplot as plt
    plt.style.use('seaborn-poster') # Bug with normal x marker !!!

    if (ax == None):
        ax = plt.subplots()[1]
    # irfname =  self.filename.parts[-2]

    # Add lower bin to exisiting array
    unit = "GeV"
    e_edges = np.append(np.array([10.,20.]), 
                        mcf.erec_edges[subarray].to(unit).value)
    e2plot = MapAxis.from_edges(e_edges,unit=unit,name="energy",interp="log")
    radii = irf.irf['psf'].containment_radius(energy = e2plot.edges,
                                             theta    = mcf.offset[subarray],
                                             fraction = mcf.containment)[0]

    ax.plot(e2plot.edges.to(eunit).value,radii.value,
            marker="o",alpha=0.5, label= tag)
    ax.axvline(erec_min.to(eunit).value,ls=":")
    ax.axvline(erec_max.to(eunit).value,ls=":")
    ax.set_xlabel("Energy ("+eunit+")")
    ax.set_ylabel("Containment radius (°)")
    ax.set_xscale("log")
    ax.legend()

    return
def onoff_sketch_plot(irf,Emin=30*u.GeV,subarray=None, tag=None,
                      nmaxcol=4, debug=False):
    """
    On-off geometry skeches as a function of reconstrcuted energy edges
    for the current IRF

    Parameters
    ----------
    eunit : String, optional
        Energy unit used. The default is "GeV".
    debug : Boolean, optional
        Let's' say a word on what is done. The default is False.

    Returns
    -------
    None.

    """

    # Retrieve radii at EDGES
    # Assumes that Erec = Etrue as the PSF is given versus Etrue
    radii = irf.irf['psf'].containment_radius(energy    = mcf.erec_edges[subarray],
                                               theta    = mcf.offset[subarray],
                                               fraction = mcf.containment)[0]
    if (debug):
        print(72*"=")
        print(radii.value)

    ###-----------------------------
    def onoffsketch(radius,energy,ax=None):
        """
        Plot the on-off sketch n the field of view from a containment
        radius at a given energy (used for labelling only)
        """

        if (ax == None): ax = plt.subplots()

        # On region
        x_on = mcf.offset[subarray].value
        y_on = 0

        # Default on region (useless)
        on_reg = plt.Circle( (x_on, y_on), mcf.on_size[subarray].value ,
                             fill = True, color="tab:blue", alpha =0.1)
        # 68% Aeff contained region
        on68 = plt.Circle( (x_on,y_on), radius.value,
                            fill = True, color="green", alpha =0.5)
        # ax.text(x_on,y_on,s="on")

        # Build label
        txt = str(round( energy.value,1) ) + " " + str(Emin.unit)
        # txt+= " - "+ str(round(100*mcf.containment,0)) + "%"
        # ax.set_title('Field of view -'+txt )
        ax.add_artist( on_reg )
        ax.add_artist( on68)
        ax.legend([on68], [txt] )

        # Equal acceptance circle
        accept =  plt.Circle( (0, 0), mcf.offset[subarray].value ,
                            fill = False, color="black",ls="--", alpha =1)
        ax.add_artist( accept )

        # Off regions
        for i in range(1, noff_reg):
            theta = i*dtheta
            x = x_on*np.cos(theta)
            y = x_on*np.sin(theta)
            off_reg = plt.Circle( (x,y ), radius.value ,
                                 fill = True, color="red",alpha = 0.5)
            ax.add_artist( off_reg )

        # Mark center of FOV
        ax.axvline(x=0,ls=":")
        ax.axhline(y=0,ls=":")

        # Set field of view, and aspect ratio
        view = mcf.fov.value/2
        ax.set_xlim([-view,view])
        ax.set_ylim([-view,view])

        ax.set_aspect( 1 ) # Nice but conflicts with grid spacing

        return
    ###-----------------------------

    noff_reg = int(1/mcf.alpha)+1
    dtheta   = 360*u.degree/noff_reg
    #print(noff_reg," regions separated by ", dtheta)

    xsize    = 3
    ysize    = 3
    erec = mcf.erec_edges[subarray]
    mask = mcf.erec_edges[subarray]>=Emin
    nplots   = len(erec[mask]) -1 #
    ifirst   = np.where(mask)[0][0]
    ncols    = min(nmaxcol,nplots) # If nplots < nmaxcol, take nplots
    nrows    = int(nplots/ncols)+ 1*(nplots%ncols != 0)

    f_adj = 0.92 # adjusts the size when plot scale is forced to squre
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows,
                           figsize=(f_adj*xsize*ncols,ysize*nrows),
                           sharex=True, sharey=True)

    iplot = 0
    import itertools

    for jrow, icol in itertools.product(range(nrows), range(ncols)):

        ax0 = ax[jrow][icol] if (nrows>1) else ax[icol]

        if iplot < nplots:
            radius = radii[iplot+ifirst]
            energy = mcf.erec_edges[subarray][iplot+ifirst].to(Emin.unit)
            onoffsketch(radius,energy,ax=ax0)
        else:
            ax0.axis('off')
            continue # Next plot

        # Compactify
        if (jrow+1 != nrows): ax0.set_xlabel(None)
        if (icol !=0): ax0.set_ylabel(None)
        ax0.tick_params(which='major', length=10, width=2, direction='in')
        ax0.tick_params(which='minor', length=5, width=2, direction='in')

        iplot+=1

    # Figure title
    fig.suptitle(irf.filename.parts[-2]+tag,fontsize=18, y = 1.02)
    
    # know feature/bug : does not take supttile into account !
    # fig.suptitle("Title centered above all subplots", fontsize=14)
    fig.tight_layout(h_pad=0,w_pad=0)

    # Adjust AFTER tight_layout
    #plt.subplots_adjust(top=0.4,hspace=None,wspace=None)
    plt.subplots_adjust(hspace=0,wspace=0,top=0.95)

    return

###------------------------------------------------------------------------
def  aeff_plot(irf,min_fraction = 0.05, unit="GeV"):

    # Effective area
    effarea = irf.irf["aeff"].data
    e_edges = effarea.axes[0].center
    for j, off in enumerate([0*u.deg, 0.5*u.deg, 1*u.deg]):
        axij = axi[j]
        effoff = irf.irf["aeff"].to_effective_area_table(off)
        effmax = effoff.max_area # or max( effarea.evaluate(energy_true=e_edges,offset=off) )
        emin    = effoff.find_energy(effmax*min_fraction)
        print(" {:5.1f} ".format(emin[0].to(unit).value),end="")
        with quantity_support():
            
            
            label = str(off.value)+"° " + "-" + (irf.filename.parts[-2])[6:]
            
            # effarea.evaluate(energy_true = e_edges,offset =off))
            p = axij.plot(e_edges,
                          effoff.data.evaluate(energy_true = e_edges)/1e6,
                          label=label)
            axij.axhline(y=min_fraction*effmax/1e6,ls=":",color=p[0].get_color())
            # axij.axvline(x=irf.ereco_min,
            #              ls=":",color=p[0].get_color(),
            #              label="Emin")
            # axij.axvline(x=irf.ereco_max,
            #              ls=":",color=p[0].get_color(),
            #              label="Emax")
            axij.set_xscale("log")
            axij.set_yscale("log")
            axij.legend()
            if j>0 : axij.set_ylabel(None)
            if i>1 : axij.set_xlabel(None)
            axij.grid("both",which="major",alpha=0.8)
            axij.grid("both",which="minor",alpha=0.5)
    print()
            #axis.legend()
    plt.tight_layout(h_pad=0, w_pad=0)
    
    return

###############################################################################
if __name__ == "__main__":
    """
    Code example to use the IRF class
    """
    
    # Suppress invalid unit warning when reading IRF files
    import logging
    logging.basicConfig()
    log = logging.getLogger("gammapy.irf")
    log.setLevel(logging.ERROR)

    from utilities import t_str
    import matplotlib.pyplot as plt
    # plt.style.use('seaborn-talk') # Make the labels readable
    plt.style.use('seaborn-poster') # Bigger - bug with normal x marker !!!
    
    import gammapy
    print(" Running with gammapy ",gammapy.__version__)

    # irf_dir = r"D:\CTA\00-Data\IRF-SoHAPPy\prod3-v2"
    # array   = {"North":"FullArray", "South":"FullArray"} # "FullArray", "LST",...
    # array   = {"North":"LST", "South":"LST"} # "FullArray", "LST",...
    # array   = {"North":"MST", "South":"MST"} # "FullArray", "LST",...

    irf_dir = r"D:\CTA\00-Data\IRF-SoHAPPy\prod5-v0.1"
    array   = {"North":"4LSTs09MSTs", "South":"14MSTs37SSTs"}

    idx = irf_dir.find("prod")
    prod = irf_dir[idx:idx+5]
    print(" Processing ",prod," files")
    show = []
    # show = ["containment", "onoff", "effearea","generic"]
    # show = ["effarea", "containment"]
    show = ["containment","onoff"]
        
    ###-------------------------
    ### Generic    
    ###-------------------------
    if "generic" in show:
        for loc in ["North","South"]:
            irf = IRF.from_observation(loc=loc,
                                       subarray   = array[loc],
                                       zenith  = 20*u.deg,
                                       azimuth = 123*u.deg,
                                       obstime = 100*u.s,
                                       irf_dir = irf_dir )
            irf.print()
            irf.irf["psf"].peek()
            irf.irf["edisp"].peek()
    
        print(dt_log_valid)
    
    ###-------------------------
    ### Show containment radii    
    ###-------------------------
    if "containment" in show:
        for loc in ["North","South"]:
            fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(15,5), 
                                   sharex=True, sharey=True)
            iplot=0
            for ztag, axi in zip(zenith_list,ax):
                for dt in dt_list[prod].values():           
                    irf = IRF.from_observation(loc       = loc,
                                               subarray = array[loc],
                                               zenith   = zenith_list[ztag],
                                               azimuth  = 123*u.deg,
                                               obstime  = dt,
                                               irf_dir  = irf_dir )    
                    print("---->",irf.filename,dt)
                    containment_plot(irf, subarray=array[loc],
                                     erec_min = mcf.erec_min[array[loc]][ztag],
                                     erec_max = mcf.erec_max[ztag],
                                     ax=axi,
                                     tag=ztag[:2]+"° - "+t_str(dt))
                    
                    axi.grid("both",which="major",alpha=0.8)
                    axi.grid("both",which="minor",alpha=0.5)
                    axi.set_ylim(ymax=0.5)
                    
                    if (iplot>0): axi.set_ylabel(None)
                    iplot +=1
            
            fig.suptitle(array[loc]+" "+loc,fontsize=18,y=0.95)
            plt.tight_layout(w_pad=0)
            # plt.subplots_adjust(top=0.95)
            
    ###-------------------------
    ### Show on-off sketch / containment radii
    ###-------------------------   
    nmaxcol = 6
    if "onoff" in show: 
        for loc in ["North","South"]:
            for ztag in zenith_list:
                for dt in dt_list[prod].values():
    
                    irf = IRF.from_observation(loc       = loc,
                                               subarray = array[loc],
                                               zenith   = zenith_list[ztag],
                                               azimuth  = 123*u.deg,
                                               obstime  = dt,
                                               irf_dir  = irf_dir )
                    print("---->",irf.filename)
                    Emin = mcf.erec_min[array[loc]][ztag]
                    onoff_sketch_plot(irf,Emin=Emin,nmaxcol=nmaxcol,
                                      subarray=array[loc],
                                      tag=" - "+ztag[:2]+"° - "+t_str(dt))
                                       
                    plt.show()

    ###-------------------------
    ### Effectiva area plots    
    ###-------------------------
    if "effarea" in show: 
        unit = "GeV"
        min_fraction = 0.05 # 0.05, 0.1
    
        print(" *** Threshold for {:2.0f}% eff.max ({:3s})  *** "
              .format(100*min_fraction,unit))
    
        for loc in ["North","South"]:
            fig, ax = plt.subplots(nrows=3,ncols=3, figsize=(15,15),
                                        sharex=True, sharey=True)
            fig.suptitle(array[loc] + "-" + loc,fontsize=30)
    
            print("{:5s} {:10s} {:7s} {:7s} {:7s}"
                  .format(loc,array[loc],"0°","0.5°","1.0°"))
    
            i=0
            for z, axi in zip([20, 40,57],ax):
        #       print(axi)
                for dt in [100*u.s,0.5*u.h, 5*u.h, 50*u.h]:
                    irf = IRF.from_observation(loc      = loc,
                                               subarray = array[loc],
                                               zenith   = z*u.deg,
                                               azimuth  = 123*u.deg,
                                               obstime  = dt,
                                               irf_dir  = irf_dir )
                    #print(" Found : ",irf.filename)
                    print("{:13s} :"
                          .format(str(z)+"° "+str(dt.value)+" "+str(dt.unit)),end="")
    
                    aeff_plot(irf,unit=unit,min_fraction=min_fraction)
                    i+=1
            plt.subplots_adjust(top=0.95)

