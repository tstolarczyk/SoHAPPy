# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 17:22:01 2023

@author: Stolar
"""

import matplotlib
import matplotlib.pyplot as plt

import astropy.coordinates as coord
import astropy.units as u
from astropy.cosmology import Planck13 as cosmo # astropy v5
# from astropy.coordinates.distances import Distance
from   astropy.visualization import quantity_support

from niceplot import col_size, vals_legend, draw_sphere

__all__ = ["universe_coverage"]
###----------------------------------------------------------------------------
def universe_coverage(sub_pop):
    """


    Parameters
    ----------
    sub_pop : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """


    # if matplotlib.__version__ < "3.5":
    #     print(" Does not work with matplotlib ",matplotlib.__version__)
    #     return

    fig = plt.figure(figsize=(10,10))
    ax=fig.add_subplot(projection='3d')
    if matplotlib.__version__ >= "3.5":
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
    fig.legend(title=r"$\sigma_{max}$",handles=patches,bbox_to_anchor=[1.15, 0.800],ncol=1)

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
