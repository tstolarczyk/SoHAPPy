# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 11:58:24 2021

@author: Stolar
"""
import numpy as np
import matplotlib.pyplot as plt

from niceplot import MyLabel, single_legend, stamp, projected_scatter, col_size, vals_legend
from historical import plot_historical, historical

###-----------------------------------------------------------------------------
def distri(pop, grbs,
            var="z", varmin=0, varmax=4, varlog = False,
            varname=None, tag="dummy",
            ax=None,
            reference=True, ratio=True,
            nbin=20, color="tab:blue", color2="red",
            **kwargs):
    """
    

    Parameters
    ----------
    pop : TYPE
        DESCRIPTION.
    grbs : TYPE
        DESCRIPTION.
    var : TYPE, optional
        DESCRIPTION. The default is "z".
    varmin : TYPE, optional
        DESCRIPTION. The default is 0.
    varmax : TYPE, optional
        DESCRIPTION. The default is 4.
    varlog : TYPE, optional
        DESCRIPTION. The default is False.
    varname : TYPE, optional
        DESCRIPTION. The default is None.
    tag : TYPE, optional
        DESCRIPTION. The default is "dummy".
    ax : TYPE, optional
        DESCRIPTION. The default is None.
    reference : TYPE, optional
        DESCRIPTION. The default is True.
    ratio : TYPE, optional
        DESCRIPTION. The default is True.
    nbin : TYPE, optional
        DESCRIPTION. The default is 20.
    color : TYPE, optional
        DESCRIPTION. The default is "tab:blue".
    color2 : TYPE, optional
        DESCRIPTION. The default is "red".
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    ax : TYPE
        DESCRIPTION.
    axx : TYPE
        DESCRIPTION.

    """
    

    ax = plt.gca() if ax is None else ax

    # If varaibale namenot explcitely given, use the column name
    if varname is None:
        varname = var

    # Plot the reference "1000" GRB population
    if reference == True:
        mask = (pop.ref[var] <= varmax) & (pop.ref[var] >= varmin)

        if varlog:
            x = [np.log10(v) if v>0 else 0 for v in pop.ref[mask][var]]
        else:
            x = pop.ref[mask][var]

        nref, bins, _  = ax.hist(x, bins=nbin,
                                 facecolor="none",edgecolor="black")
    else:
        bins = nbin
        nref = 0

    # Plot the requested population
    mask = (grbs[var] <= varmax) & (grbs[var] >= varmin)
    if varlog:
        x  =[ np.log10(v) if v>0 else 0 for v in grbs[mask][var]]
    else:
        x = grbs[mask][var]
    n, bins, _ = ax.hist(x,
                         bins=bins, color=color,
                         label=MyLabel(x,label=tag),
                         **kwargs)
    ax.legend()
    ax.set_xlabel(varname)
    ax.set_yscale("log")

    # Plot the ratio
    ratio = [ (x/xtot if xtot!=0 else 0) for x, xtot in zip(n, nref) ]
    axx = ax.twinx()
    axx.plot(bins[:-1] + 0.5*(bins[1:]-bins[:-1]),ratio,
             alpha=1,color=color2,label="Ratio")

    axx.grid(ls=":")
    axx.legend(loc="lower right")

    return ax, axx

###----------------------------------------------------------------------------
def coverage(varx, vary,
             ref = None, pop = None, mask = None,
             xrange = [None,None], yrange = [None, None],
             lblx   = "",          lbly   ="",
             xscale = "log",       yscale ="log",
             title  = "dummy",
             nbin=25, alpha=0.5):
    """
    

    Parameters
    ----------
    varx : TYPE
        DESCRIPTION.
    vary : TYPE
        DESCRIPTION.
    ref : TYPE, optional
        DESCRIPTION. The default is None.
    pop : TYPE, optional
        DESCRIPTION. The default is None.
    mask : TYPE, optional
        DESCRIPTION. The default is None.
    xrange : TYPE, optional
        DESCRIPTION. The default is [None,None].
    yrange : TYPE, optional
        DESCRIPTION. The default is [None, None].
    lblx : TYPE, optional
        DESCRIPTION. The default is "".
    lbly : TYPE, optional
        DESCRIPTION. The default is "".
    xscale : TYPE, optional
        DESCRIPTION. The default is "log".
    yscale : TYPE, optional
        DESCRIPTION. The default is "log".
    title : TYPE, optional
        DESCRIPTION. The default is "dummy".
    nbin : TYPE, optional
        DESCRIPTION. The default is 25.
    alpha : TYPE, optional
        DESCRIPTION. The default is 0.5.

    Returns
    -------
    None.

    """

    fig, ax, axh, axv = projected_scatter()

    ### --------------------
    ### Central scatter plot
    ### --------------------
    # Generated population

    ax.scatter(np.log10(ref[varx]) if xscale=="log" else ref[varx],
               np.log10(ref[vary]) if yscale=="log" else ref[vary],
               marker=".", color="black",s=10, alpha=0.5, label="All")

    # Detected population
    colors, sizes = col_size(pop[mask].sigmx)
    ax.scatter(np.log10(pop[mask][varx]) if xscale=="log" else pop[mask][varx],
               np.log10(pop[mask][vary]) if yscale=="log" else pop[mask][vary],
               marker="o", s= sizes, c=colors, alpha = alpha)



    plot_historical(ax,historical(), obs=["H.E.S.S"])

    # Decoration
    if xscale=="log":
        ax.set_xlim(xmin=np.log10(xrange[0]), xmax=np.log10(xrange[1]))
        ax.set_xlabel("$log_{10} \ $"+lblx)
    else:
        ax.set_xlim(xmin=xrange[0],xmax=xrange[1])
        ax.set_xlabel(lblx)

    if yscale == "log":
        ax.set_ylim(ymin=np.log10(yrange[0]),ymax=np.log10(yrange[1]))
        ax.set_ylabel("$log_{10} / $" +lbly)
    else:
        ax.set_ylim(yrange)
        ax.set_ylabel(lbly)

    ax.grid("both",ls="--")

    patches = vals_legend(ax, alpha = alpha)
    fig.legend(title="$\sigma_{max}$",handles=patches, bbox_to_anchor=[1.04, 1.01],ncol=2)

    single_legend(ax)

    ### --------------------
    ### horizontal data
    ### --------------------
    hist_mask = (ref[varx] >= xrange[0]) & (ref[varx] <= xrange[1])

    xref = np.log10(ref[hist_mask][varx]) if xscale=="log" else ref[hist_mask][varx]
    n, bins, _ = axh.hist(xref, bins = nbin, facecolor="none",edgecolor="black")

    x = np.log10(pop[mask][varx]) if xscale=="log" else pop[mask][varx]
    axh.hist(x, bins = bins, color="purple",alpha=0.3)

    if xscale == "log":
        axh.set_xlim(xmin = np.log10(xrange[0]), xmax=np.log10(xrange[1]))
    else:
        axh.set_xlim(xmin = xrange[0], xmax=xrange[1])

    axh.set_title(tag+" - " + title)
    axh.set_yscale("log")

    ### --------------------
    ### Vertical data
    ### --------------------
    hist_mask = (ref[vary] >= yrange[0]) & (ref[vary] <= yrange[1])

    yref = np.log10(ref[hist_mask][vary]) if yscale=="log" else ref[hist_mask][vary]
    n, bins, _ = axv.hist(yref, bins = nbin,
                          facecolor="none",edgecolor="black",orientation="horizontal",
                          label=MyLabel(yref ,label="All"))

    y = np.log10(pop[mask][vary]) if yscale=="log" else pop[mask][vary]
    axv.hist(y,bins = bins,
             color="purple",alpha=0.3,orientation="horizontal",
             label = MyLabel(y ,label="Detected"))
    axv.set_xscale(yscale)

    return

###############################################################################
if __name__ == "__main__":

    # A standalone function to read a GRB and make various tests

    from population import Pop
    from pop_io import create_csv

    import sys
    codefolder = "../../"
    sys.path.append(codefolder)

    nyears, file, _ = create_csv(file="parameter.yaml",debug=True)

    pop = Pop(filename=file, nyrs= nyears)
    popNS = pop.grb[(pop.grb.loca=="North") | (pop.grb.loca=="South")]


    # Plot one-dim varaibale coverage
    fig, ax = plt.subplots(nrows=3,ncols=2,figsize=(15,10))

    for icol, grbs, tag in zip(range(2),
                               [popNS, pop.g_tot],
                               [r"North & South 5$\sigma$","Combined 5$\sigma$"]):

        print(icol)
        # Selection
        grbs = grbs[grbs.d5s >= pop.eff_lvl]

        distri(pop, grbs, var="z", tag=tag, varmax=5, ax=ax[0][icol],
               varname="Redshift (z)")

        distri(pop, grbs ,var="Eiso", tag=tag, varmin=5.5e50, varmax=5.5e55,
               varlog=True, ax=ax[1][icol], varname=r"$log_{}10 E_{iso}$")

        distri(pop, grbs ,var="Epeak", tag=tag, varmin=10, varmax=1e5,
               varlog=True, ax=ax[2][icol])

        plt.tight_layout()

    ### Plot coverage in 2-dim space
    poplist  = [pop.g_n,pop.g_s,pop.g_tot]
    taglist  = ["North", "South", "Combined"]
    masklist = [pop.g_n.d5s   >= pop.eff_lvl,
                pop.g_s.d5s   >= pop.eff_lvl,
                pop.g_tot.d5s >= pop.eff_lvl]
    title = r"Detected at $5\sigma \ (90\% \ C.L.)$"

    # Eiso versus z
    varx    = "z"
    vary    = "Eiso"
    lblx    = "$z$"
    lbly    = "$E_{iso} \ (erg)$"
    xrange  = [0,5]
    yrange  = [1e50, 5e55]

    alpha = 0.5

    for grbs, mask, tag in zip(poplist, masklist, taglist):

        coverage(varx,vary,ref=pop.ref, pop=grbs,
                     xrange=xrange,yrange=yrange,
                     lblx = lblx, lbly = lbly,
                     xscale="linear",
                     mask=mask,
                     title=title)
        stamp(pop.tag, axis=fig, x=1.04, y = 0.5, rotation=270)

    # Eiso versus Epeak
    varx    = "Epeak"
    vary    = "Eiso"
    lblx    = r"$E_{peak} \ (keV)$"
    lbly    = r"$E_{iso} \ (erg)$"
    xrange  = [10,1e5]
    yrange  = [1e50, 5e55]

    for grbs, mask, tag in zip(poplist, masklist, taglist):

        coverage(varx,vary,ref=pop.ref, pop=grbs,
                     xrange=xrange,yrange=yrange,
                     lblx = lblx, lbly = lbly,
                     mask=mask,
                     title=title)
        stamp(pop.tag, axis=fig, x=1.04, y = 0.5, rotation=270)