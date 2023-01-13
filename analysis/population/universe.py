# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 17:22:01 2023

@author: Stolar
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource

from mpl_toolkits import mplot3d


import astropy.coordinates as coord
import astropy.units as u
from astropy.cosmology import Planck13 as cosmo # astropy v5
from astropy.coordinates.distances import Distance
from   astropy.visualization import quantity_support

from niceplot import col_size, vals_legend

###----------------------------------------------------------------------------
def draw_sphere(radius=1, colormap=plt.cm.viridis,ax=None, **kwargs):

    if ax==None:
        fig = plt.figure(figsize=(8,8), dpi=300)
        ax = fig.add_subplot(111, projection='3d')
    
    ax.set_box_aspect(aspect = (1,1,1))

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = radius* np.outer(np.cos(u), np.sin(v))
    y = radius* np.outer(np.sin(u), np.sin(v))
    z = radius* np.outer(np.ones(np.size(u)), np.cos(v))
    
    
    ls = LightSource(azdeg=0, altdeg=65)
    rgb = ls.shade(z, colormap)

    ax.plot_surface(x, y, z,  rstride=1, cstride=1, 
#                     color=color,   
                    facecolors=rgb, linewidth=0, **kwargs)

    return
  
###----------------------------------------------------------------------------
def universe_coverage(sub_pop):
    
    import matplotlib
    if matplotlib.__version__ < "3.5":
        print(" Does not work with matplotlib ",matplotlib.__version__)
        return
    
    fig = plt.figure(figsize=(10,10)) 
    ax=fig.add_subplot(projection='3d') 
    ax.set_box_aspect(aspect = (1,1,1))

    colors, sizes= col_size(sub_pop.sigmx)

    distance = cosmo.luminosity_distance(sub_pop.z).to('Gpc')
    c1 = coord.SkyCoord(ra=sub_pop.ra*u.deg, dec=sub_pop.dec*u.deg,
                        distance=distance)

    gc1 = c1.transform_to(coord.Galactocentric)
    x = gc1.x.to('Gpc')
    y = gc1.y.to('Gpc')
    z = gc1.z.to('Gpc')
    dmax = max(max(x),max(y),max(z))
    print(dmax)

    with quantity_support():

        ax.scatter(gc1.x, gc1.y, gc1.z, zdir='z',
                   marker="o",s=sizes,color=colors)
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_zlabel(None)
        ax.set_xlim(xmin=-dmax, xmax=dmax)
        ax.set_ylim(ymin=-dmax, ymax=dmax)
        ax.set_zlim(zmin=-dmax, zmax=dmax)

    draw_sphere(radius=dmax.value,ax=ax,alpha=0.1)
    draw_sphere(radius=cosmo.luminosity_distance(2.5).to('Gpc').value,ax=ax,alpha=0.1)


    patches=vals_legend(ax)
    fig.legend(title="$\sigma_{max}$",handles=patches,bbox_to_anchor=[1.15, 0.800],ncol=1)

    plt.tight_layout()
################################################################################################
if __name__ == "__main__":

    # A standalone function to read a GRB and make various tests

    from population import Pop
    from pop_io import create_csv

    import sys
    codefolder = "../../"
    sys.path.append(codefolder)

    nyears, file, _ = create_csv(file="parameter.yaml",debug=True)
    pop = Pop(filename=file, nyrs= nyears)
    
    # draw_sphere(radius=2.3)  
    universe_coverage(pop.g_tot[pop.g_tot.d5s>pop.eff_lvl])
