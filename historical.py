# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 10:38:19 2022

@author: Stolar
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

history = {"190114C":
           {"Observatory": "MAGIC",
            "z": 0.4245,
            "Eiso": 3.00E+53*u.erg,
            "t90": 25*u.s,
            "marker": "^",
            "col": "black"},

           "201216C":
           {"Observatory": "MAGIC",
            "z": 1.1,
            "Eiso": 5.00E+53*u.erg,
            "t90": 30*u.s,
            "marker": ">",
            "col": "black"},

           "180720B":
           {"Observatory": "H.E.S.S.",
            "z":  0.654,
            "Eiso": 6.00E+53*u.erg,
            "t90": 49*u.s,
            "marker": "v",
            "col": "red"},

           "190829A":
           {"Observatory": "H.E.S.S.",
            "z": 0.0785,
            "Eiso": 2.00E+50*u.erg,
            "t90": 63*u.s,
            "marker": "<",
            "col": "red"},

           "080916C":
           {"Observatory": "Fermi/LAT",
            "z": 4.3,
            "Eiso": 8.80E+54*u.erg,
            "t90": -1*u.s,
            "marker": "*",
            "col": "tab:green"},

           "090902B":
           {"Observatory": "Fermi/LAT",
            "z": 1.822,
            "Eiso": 2.20E+52*u.erg,
            "t90": -1*u.s,
            "marker": "*",
            "col": "tab:green"},

           "130427A":
           {"Observatory": "Fermi/LAT",
            "z": 0.34,
            "Eiso": 9.60E+53*u.erg,
            "t90": -1*u.s,
            "marker": "*",
            "col": "tab:green"}
           }

__all__ = []


# ##----------------------------------------------------------------------------
def historical():
    return history


# ##----------------------------------------------------------------------------
def plot_historical(ax, ddict, obs=None):
    """
    Add z, Eiso from already detected GRBs

    Parameters
    ----------
    ax : matplotlib axis
        Current axis
    ddict : Dictionnay
        A dictionnary with predefined values
    obs : String, optional
        Observatories to be shown from the dictionnary. The default is [].

    Returns
    -------
    ax : matpotlib axis
        Current axis

    """

    ax = plt.gca() if ax is None else ax

    for target, _ in ddict.items():
        data = ddict[target]
        if data["Observatory"] in obs or obs == "all":
            print(target)
            ax.scatter(data["z"],
                       np.log10(data["Eiso"].value),
                       marker=data["marker"],
                       color=data["col"],
                       label=target)
    return ax


###############################################################################
if __name__ == "__main__":

    fig, ax = plt.subplots()
    plot_historical(ax, history, obs="all")
    ax.legend()
