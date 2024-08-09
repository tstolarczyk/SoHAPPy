# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:45:41 2024

@author: Stolar
"""
import sys

import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

from niceplot import single_legend
from population import Pop
from pop_io import get_data

# Bigger texts and labels
sns.set_context("notebook") # poster, talk, notebook, paper

codefolder = "."
sys.path.append(codefolder)

#------------------------------------------------------------------------------
def contour(ax, pop):
    """
    Under test

    Parameters
    ----------
    ax : TYPE
        DESCRIPTION.
    pop : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    H, xe, ye = np.histogram2d(pop.g_tot.z,np.log10(pop.g_tot.Eiso), 10)
    # extent = np.array([xe.min(), xe.max(), ye.min(), ye.max()])

    plt.pcolormesh(xe, ye, H, cmap="rainbow")

    midpoints = (xe[1:] + xe[:-1])/2, (ye[1:] + ye[:-1])/2
    # cs = plt.contourf(*midpoints, H, levels=2, extend='max', cmap='plasma', alpha=0.5, zorder=100)
    ax.contour(*midpoints, H,levels=2,colors='k',zorder=5)



#------------------------------------------------------------------------------
def lower_limit(pop, varx, vary, mask,
                ax=None, logx=False, logy=False,
                prt      = False,
                plot_lim = False,
                color    = "red",
                fitdeg   = 3,
                bins     = 25,
                cmap     = "magma_r"):


        plot = False if ax is None else True

        xdet = pop[varx][mask]
        if logx:
            xdet = np.log10(xdet)

        ydet = pop[vary][mask]
        if logy:
            ydet = np.log10(ydet)

        if plot:
            H, xe, ye, img = ax.hist2d(xdet, ydet, bins=bins, cmap=cmap)
        else:
            H, xe, ye = np.histogram2d(xdet, ydet, bins=bins)

        # Bin centre and heights
        xctr = (xe[1:] + xe[:-1])/2
        yctr = (ye[1:] + ye[:-1])/2
        dy   = ye[1:] - ye[:-1]

        ylow = []
        xlow = []

        for ibin in range(len(xctr)):

            # Display counts on each cell
            if prt:
                for jbin in range(len(yctr)):
                    ax.text(xctr[ibin],yctr[jbin],str(int(H[ibin, jbin])) )

            # Find position of minimal value at bottom of each column
            ycol = H[ibin,:]
            ipos = np.where(ycol>0)[0]

            # Print found positions on screen and draw a circle on the plot
            if len(ipos) != 0:
                xlow.append(xctr[ibin])
                ylow.append(yctr[ipos[0]] - 1.5*dy[ipos[0]])
                if prt:
                    print("bin  #",ibin," : ",ycol, ipos, yctr[ipos])
                    ax.text(xctr[ibin],yctr[ipos[0]],"O", size=20)
            else:
                if prt:
                    print("bin  # No data")

        # Fit low value points with a polynomial
        pfit = np.polyfit(xlow,ylow,fitdeg)
        # from scipy.interpolate import UnivariateSpline
        # spl = UnivariateSpline(xlow, ylow, k=3)

        ### Plot results if required
        if plot:
            if plot_lim:
                ax.plot(xlow,ylow,color="black",lw=1, ls=":",label="lower limits")

            ax.plot(xlow,np.polyval(pfit,xlow),color=color,lw=2,
                    label=" Fit order="+str(fitdeg))

            # ax.plot(xlow,spl(xlow),color="tab:green",lw=2,
            #         label=" Spline"+str(fitdeg))

            if logx:
                ax.set_xlabel("$Log_{10}( "+varx+")$")
            else:
                ax.set_xlabel("$"+varx+"$")

            if logy:
                ax.set_ylabel("$Log_{10}("+vary+")$")
            else:
                ax.set_ylabel("$"+vary+"$")

            ax.legend()

            fig = plt.gcf()
            cbar = fig.colorbar(img, ax=ax)
            cbar.set_label('Counts')

        return pfit

###############################################################################
if __name__ == "__main__":

    import os
    os.environ["HAPPY_IN"] = "D:\\CTAO\SoHAPPy\input"
    os.environ["HAPPY_OUT"] = "D:\\CTAO\SoHAPPy\output"

    # Use default as a demo
    # nyears, files, tag = get_data(parpath=None,debug=True)

    # Select a particular data set
    nyears, files, tag = get_data(parpath="../../data/samples/max_detection_parameter.yaml",debug=False)

    # Read population
    pop   = Pop(files, tag=tag, nyrs= nyears)


    poplist  = [pop.g_n,pop.g_s,pop.g_tot]
    taglist  = ["North", "South", "Combined"]
    masklist = [(pop.g_n.d5s   >= pop.eff_lvl) & (pop.g_n.z <5),
                (pop.g_s.d5s   >= pop.eff_lvl) & (pop.g_s.z <5),
                (pop.g_tot.d5s >= pop.eff_lvl) & (pop.g_tot.z <5)]

    # Eiso versus z
    varx    = "z"
    vary    = "Eiso"

    first = True
    for grbs, mask, tag in zip(poplist, masklist, taglist):

        fig, ax = plt.subplots(figsize=(10,7))

        if  first:
            pfit0 = lower_limit(grbs, varx, vary, mask,
                                ax = ax, logx=False, logy=True,
                                bins=25, fitdeg=2, color="red")
            first = False
        else:
            pfit = lower_limit(grbs, varx, vary, mask,
                               ax =ax, logx=False, logy=True,
                                bins=25, fitdeg=2, color="black")
            xlim = ax.get_xlim()
            xsample = np.linspace(xlim[0],xlim[1],20)
            ax.plot(xsample, np.polyval(pfit0, xsample), ls="--", alpha=0.5,
                    color= "red", label="Reference")

        ax.set_xlabel(varx)


        ax.set_title(tag)
        single_legend(ax)






































