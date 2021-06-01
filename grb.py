# -*- coding: utf-8 -*-
"""
This module contains the classes and tools to handle a GRB object, i.e.
a collection of lightcurve bins for each of which is given an energy spectrum,
and information on site location.

Created on Tue Jan 22 11:41:34 2019

@author: Stolar
"""
# import yaml
# from yaml import CLoader as Loader
#from glob import glob

import sys
import warnings

from pathlib import Path

import numpy as np

import astropy
import astropy.units         as u
from   astropy.table         import Table
from   astropy.io            import fits
from   astropy.time          import Time
from   astropy.coordinates   import SkyCoord
from   astropy.coordinates   import EarthLocation
from   astropy.coordinates   import AltAz

from visibility import Visibility

from gammapy.modeling.models import Absorption, AbsorbedSpectralModel
from gammapy.modeling.models import TemplateSpectralModel
from gammapy.modeling.models import PointSpatialModel


# The following information can be obtained from EarthLocation.of_site
# but this requires to have an internet connection
# Get all available site names :
# sitelist = EarthLocation.get_site_names()

# Get coordinattes of Paranal and LaPalma
# Note : with the Gammapy-0.12 environment,this requires to
#        downgrade numpy to version 1.16.2 (instead on >1.18)
#        See discussion here :
#        https://github.com/astropy/astropy/issues/9306

# xyz_north = EarthLocation.of_site('Roque de los Muchachos')
# xyz_south = EarthLocation.of_site('Paranal Observatory')

xyz_north = EarthLocation.from_geocentric( float(5327448.9957829),
                                          float(-1718665.73869569),
                                            float(3051566.90295403),
                                            unit="m")
xyz_south = EarthLocation.from_geocentric( float(1946404.34103884),
                                          float(-5467644.29079852),
                                          float(-2642728.20144425),
                                            unit="m")

# Values taken from Maria Grazia Bernardini - See Readme, July 28, 2020
# xyz_north = EarthLocation.from_geocentric( 5327285.09211954,
#                                           -1718777.11250295,
#                                           3051786.7327476,
#                                           unit="m")
# xyz_south = EarthLocation.from_geocentric(1946635.7979987,
#                                           -5467633.94561753,
#                                           -2642498.5212285,
#                                           unit="m")

# xyz_north = EarthLocation.from_geodetic('342.1184','28.7606',2326.* u.meter)
# xyz_south = EarthLocation.from_geodetic('289.5972','-24.6253',2635.* u.meter)

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
    a position in the sky, a visibility window, measurement points at which
    an energy spectrum is given.
    Each energy spectrum is given as a series of flux measurements but is then
    stored as an interpolated spectrum (i.e. the energy bin width has a modest
    importance).
    The GRB properties like z, Eiso are generated in the population synthesis
    (G. Ghirlanda) and they connect the afterglow and prompt emission together
    (Other parameters like the magnetic field are independent in the two
    models).

    """

    ###########################################################################
    def __init__(self, ebl_model=None,z=0):
        """
        This initializes a default GRB.
        The default visibility has been put by hand for backward compatibility
        with the old file format and definitions where it was assumed that the
        GRB was seen without any delay. The dummy GRB is visible in North and
        South (inspired form GRB 289 in Lara's 1000 GRB files)

        Parameters
        ----------
        ebl_model : string, optional
            A normalized name of an EBL abosrtpion model (e.g. "dominguez")
            from the list available in Gammapy. The default is None.
        z : Dimensionnles quantity, optional
            The GRB default redshift. The default is 0*u.dimensionless_unscaled.

        Returns
        -------
        None.

        """
        self.name      = 'dummy'


        # Initialise absorption
        if (ebl_model != None) and (ebl_model != "built-in"):
            self.eblabs = Absorption.read_builtin(ebl_model)
        else:
            self.eblabs = None
            #print(" EBL absorption not defined or built-in")

        # GRB properties - Dummy default values
        self.z        = z
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
        self.pos_site  = {"North":xyz_north, "South":xyz_south}

        # GRB alert received
        self.t_trig   = Time('2000-01-01 02:00:00', scale='utc')

        # Flux table - Flux at a series of points
        self.Eval           = [0]*u.GeV
        self.tval           = [0]*u.s
        self.fluxval        = [0]*u.Unit("1 / (cm2 GeV s)")

        # Visibility (requires GRB points interavl)
        self.vis  = { "North": Visibility(self.radec,
                                          self.pos_site["North"]),
                      "South": Visibility(self.radec,
                                          self.pos_site["South"])
                                          }

        # GRB spectral and spatial model
        self.spectra = [] # Gammapy models (one per t slice)

        self.spatial = PointSpatialModel(lon_0=0*u.deg,lat_0=0*u.deg)

        return

    ###########################################################################
    @classmethod
    def read(cls, filename, ebl= None, newis=False, magnify=1):

        """
        Fluxes are given for a series of (t,E) values
        So far no flux is given beyond the last point.

        The spectrum is stored as a table, TableModel, that takes as an input
        a series of flux values as a function of the energy.
        The model will return values interpolated in log-space with
        math::`scipy.interpolate.interp1d`, returning zero for energies
        outside of the limits of the provided energy array.

        Parameters
        ----------
        cls : GammaRayBurst class
            GammaRayBurst class instance
        filename : STRING
            GRB input file name
        ebl : STRING, optional
            The EBL absorption model considered. If ebl is "built-in", uses
            the absorbed specrum from the data.
        magnify : float, optional
            Flux multiplicative factor for tests. Default is 1

        Returns
        -------
        A GammaRayBurst instance.

        """

        grb = cls(ebl_model=ebl) # This calls the constructor

        hdul = fits.open(filename)

        grb.name = Path(filename.name).stem

        ### GRB properties ---
        hdr = hdul[0].header
        grb.radec    = SkyCoord(hdr['RA']*u.degree,
                                hdr['DEC']*u.degree,
                                frame='icrs')
        grb.z        = hdr['Z']
        grb.Eiso     = hdr['EISO']*u.erg
        grb.Epeak    = hdr['EPEAK']*u.keV
        grb.t90      = hdr['Duration']*u.s
        grb.G0H      = hdr['G0H']
        grb.G0W      = hdr['G0W']
        grb.Fluxpeak = hdr['EISO']*u.erg
        grb.gamma_le = hdr['LOWSP']
        grb.gamma_he = hdr['HIGHSP']

        # GRB trigger time
        grb.t_trig   = Time(hdr['GRBJD']*u.day,format="jd",scale="utc")

        ### Energies, times and fluxes ---
        grb.Eval     = Table.read(hdul,hdu=2)["Energies (afterglow)"].quantity
        grb.tval     = Table.read(hdul,hdu=3)["Times (afterglow)"].quantity
        if (ebl == "built-in"):
            flux = Table.read(hdul,hdu=5)
        else:
            flux = Table.read(hdul,hdu=4)

        ### Visibilities --- includes tval span
        # One should consider not reading the default and go directly to new
        # visibilities
        grb.vis["North"].from_fits(grb, hdr, hdul,hdu=1,loc="North")
        grb.vis["South"].from_fits(grb, hdr, hdul,hdu=1,loc="South")

        #flux_unit    = u.Unit(flux.meta["UNITS"])/u.Unit("ph") # Removes ph
        # flux_unit    = u.Unit(flux.meta["UNITS"]) # Removes ph
        flux_unit = u.Unit("1/(cm2 GeV s)")

        icol_t       = len(flux.colnames)          # column number - time
        jrow_E       = len(flux[flux.colnames[0]]) # row number

        # Note the transposition from flux to fluxval
        grb.fluxval = np.zeros( (icol_t,jrow_E) )*flux_unit
        for i in range(0,icol_t):
            for j in range(0,jrow_E):
                grb.fluxval[i][j] = magnify* flux[j][i]*flux_unit # transp!

        for i,t in enumerate(grb.tval):
            # Note that TableModel makes an interpolation
            # Following a question on the Slack gammapy channel on
            # November 27th, and the answer by Axel Donath:
            # The foloowing statemetn later in the code gave an error
            # (dlist_onoff is a collection of Dataset)
            # dlist_onoff.write(datapath,prefix="grb",overwrite=True)
            # gives:
            # ...\gammapy\modeling\models\spectral.py", line 989, in to_dict
            # "data": self.energy.data.tolist(),
            # NotImplementedError: memoryview: unsupported format >f
            # This error comes from the fact that the energy list as to be
            # explicitely passed as a float as done below:
            #    grb.Eval.astype(float)
            # (A Quantity is passed as requested but the underlying numpy
            # dtype is not supported by energy.data.tolist()
            tab = TemplateSpectralModel(energy = grb.Eval.astype(float),
                                        values = grb.fluxval[i],
                                        norm   = 1.,
                                        interp_kwargs={"values_scale": "log"})

            if (grb.eblabs != None) and (ebl != "built-in"):
                model = AbsorbedSpectralModel(tab,grb.eblabs,grb.z)
            else:
                model = tab

            grb.spectra.append(model)

        hdul.close()

        return grb

    ###########################################################################
    @classmethod
    def read_prompt(cls, filename, glow= None, ebl= None,
                    z=0*u.dimensionless_unscaled, magnify=1):
        """
        Read prompt data froma file and associate it to the afterglow
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
            grb = cls(ebl_model = ebl,z=z) # This calls the constructor

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

            if (grb.eblabs != None) and (ebl != "built-in"):
                absmodel = AbsorbedSpectralModel(spectral_model = tab,
                                                 absorption    = grb.eblabs,
                                                 parameter     = grb.z)
                grb.spectra.append(absmodel)

        hdul.close()
        return grb
    ###------------------------------------------------------------------------
    def write(self, folder, debug=True):

        import pickle

        filename = Path(folder,self.name+".bin")
        outfile  = open(filename,"wb")
        pickle.dump(grb,outfile)
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
        time : Quantity (time)
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
                                                  location = self.pos_site[loc]))

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
        txt += '  Flux peak     : {:6.2f}\n'.format(self.FluxPeak)
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

        return txt

#------------------------------------------------------------------------------
if __name__ == "__main__":
    """
    A standalone function to read a GRB and make various tests
    """
    from   utilities import Log

    from SoHAPPy import get_grb_fromfile
    import grb_plot as gplt

    dbg      = 0

    ngrb     = 1 # 250
    ifirst   = 1
    save_grb = False # (False) GRB saved to disk -> use grb.py main
    res_dir  = "."

    log_filename    = res_dir  + "/grb.log"
    log = Log(name  = log_filename)

    # GRB list to be analysed
    if type(ifirst)!=list:
        grblist = list(range(ifirst, ifirst + ngrb))
    else:
        grblist = ifirst

    ## Test dummy GRB
    # grb = GammaRayBurst()
    # print(grb)
    # for loc in ["North", "South"]:
    #     print(grb.vis[loc].print(log=log))

    # Loop over GRB list
    for i in grblist:

        grb = get_grb_fromfile(i,log=log)
        print(grb)
        gplt.spectra(grb,opt="Packed")
        gplt.visibility_plot(grb, loc ="North")
        gplt.visibility_plot(grb, loc ="South")

        if (save_grb): grb.write(res_dir) # Save GRB if requested
