# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 17:40:49 2022

@author: Stolar
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm

import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from pop_io import get_config_data
from niceprint import t_fmt, heading, warning
from niceplot import MyLabel, stamp

plt.style.use('seaborn-talk') # Make the labels readable

__all__ = ["Pop"]
###############################################################################
class Pop():
    """
    This class handles poulation data from a csv file on disk.
    """

    ###------------------------------------------------------------------------
    def __init__(self, filename="data.csv", nyrs=1., nmaxsrc=100000,
                 debug=False):
        """
        Get data from a csv SoHAPPy output file
        Compute combinatory of site visibilities.

        Parameters
        ----------
        filename : string, optional
            Input file name. The default is "data.csv".
        nyrs : float, optional
            Number of years in the simulation or data sample. The default is 1.
        nmaxsrc : integer, optional
            Maximum number of sources to be read. The default is 100000.
        debug : boolean, optional
            A general flag for debugging. The default is False.

        Returns
        -------
        None.

        """

        self.dbg      = debug
        self.filename = filename
        self.tag      = filename.parent.name
        self.nyears   = nyrs # Rate normalisation

        # The full GRB population (Entries North, South and Both).
        self.grb = pd.read_csv(filename)

        # Extract "not visible" flag and iteration from the data
        self.unvis     = min(set(self.grb.err))
        self.niter     = max(set(self.grb.err))
        self.niter_3s  = max(self.grb.d3s)
        self.niter_5s  = max(self.grb.d5s)

        # If "loca" column does not exist, this is an old file using "site"
        if "loca" not in self.grb.columns:
            if "site" in self.grb.columns:
                self.grb = self.grb.rename(columns={'site': 'loca'})
                warning(" Deprecated column name 'site' changed to 'loca'")
            else:
                sys.exit(f"{__name__}.py: Missing column for location")

        # If "loca" column does not exist, this is an old file using "site"
        if "prpt" not in self.grb.columns:
            if "vis" in self.grb.columns:
                self.grb = self.grb.rename(columns={'vis': 'prpt'})
                warning(" Deprecated column name 'vis' changed to 'prpt'")
            else:
                sys.exit(f"{__name__}.py: Missing column for prompt")

        # If "loca" column does not exist, this is an old file using "site"
        if "t_trig" not in self.grb.columns:
            if "ttrig" in self.grb.columns:
                self.grb = self.grb.rename(columns={'ttrig': 't_trig'})
                warning(" Deprecated column name 'ttrig' changed to 't_trig'")
            else:
                sys.exit(f"{__name__}.py: Missing column for prompt")
        # Site populations - even if not analysed
        g_s = self.grb[self.grb.loca=="South"][:nmaxsrc]
        g_n = self.grb[self.grb.loca=="North"][:nmaxsrc]
        g_b = self.grb[self.grb.loca=="Both" ][:nmaxsrc]

        # Keep one of them for simulated GRB parameters (e.g. z)
        self.ref = g_n # Reference, no selection applied

        # Get GRB names in the data
        if  len(g_s.name) != len(g_n.name) or \
            len(g_s.name) != len(g_b.name) :
            sys.exit(f"{__name__:}.py: Inconstistency in GRB name list")

        self.names = g_s.name

        # Warn user in case the maximal GRB limit is below the data content
        if len(self.names) > nmaxsrc:
            print(" WARNING : limited statistics (",len(self.names),") >",nmaxsrc)
        self.ngrb     = min(nmaxsrc,len(self.names))

        # Add combinatory to data frame for the name list
        self.add_combinatory()

        # The iteration number of the simulation can be guessed
        # from the error code (except if all simulations failed!)
        # All iterations were processed when niter == err
        self.g_ana = self.grb[ self.grb.err == self.niter]

        # Population of all sites that could be fully analysed
        self.g_n0 = self.g_ana[ (self.g_ana.loca =="North") & (self.g_ana.N==1)] # North only
        self.g_s0 = self.g_ana[ (self.g_ana.loca =="South") & (self.g_ana.S==1)] # South only
        self.g_n = self.g_ana[  self.g_ana.loca =="North"] # North and maybe elsewhere
        self.g_s = self.g_ana[  self.g_ana.loca =="South"] # South and maybe elsewhere
        self.g_b = self.g_ana[ (self.g_ana.loca =="Both")  & (self.g_ana.B==1)] # Seen both
        # The total unique population
        self.g_tot = pd.concat([self.g_n0,self.g_s0,self.g_b],axis=0)

        # Get configuration file data and some parameters
        conf = get_config_data(self.filename, debug=False)

        if conf["niter"] != self.niter:
            sys.exit(f"{__name__:}.py: niter mismatch wrt configuration")

        if "det_level" not in conf.keys():
            warning(f"{__name__}.py: old config. file, 'det_level' set to 0.9")
            conf["det_level"] = 0.9

        self.eff_lvl = conf["det_level"]*conf["niter"]
        self.dtswift = u.Quantity(conf["dtswift"])
        self.dtslew  = {"North":u.Quantity(conf["dtslew_North"]),
                        "South":u.Quantity(conf["dtslew_South"])}

    ###-------------------------------------------------------------------
    def add_combinatory(self):
        """
        Add combinatory to data frame
        """

        # If columns do not exist, create them
        if "N" not in self.grb:
            self.grb.insert(1,"N",0) # North only
        if "S" not in self.grb:
            self.grb.insert(1,"S",0) # South only
        if "B" not in self.grb:
            self.grb.insert(1,"B",0) # North and South

        if self.dbg:
            print(f"{'name':>10s} {'N':>3s} {'S':>3s} {'B':>3s}" \
                  f" {'No':>3s} {'So':>3s}")

        for name in self.names:
            grb = self.grb[self.grb.name==name]
            seen_n = (grb[grb.loca=="North"].err!=self.unvis).bool()
            seen_s = (grb[grb.loca=="South"].err!=self.unvis).bool()
            seen_b = (grb[grb.loca=="Both" ].err!=self.unvis).bool()
            seen_sonly = seen_s & ~seen_n
            seen_nonly = seen_n & ~seen_s
            if self.dbg:
                print(f"{name:>10s} {seen_n:3d} {seen_s:3d} {seen_b:3d}"\
                      f"{seen_nonly:3d} {seen_sonly:3d}")
            #print(g.index)
            for idx in grb.index:
                # Not in pandas 1.0.3
                # self.grb.set_value(idx,"N",int(seen_nonly))
                # self.grb.set_value(idx,"S",int(seen_sonly))
                # self.grb.set_value(idx,"B",int(seen_b))
                self.grb.at[idx,"N"] = int(seen_nonly)
                self.grb.at[idx,"S"] = int(seen_sonly)
                self.grb.at[idx,"B"] = int(seen_b)

    ###------------------------------------------------------------------------
    def print(self):
        """
        Printout the class content.

        """
        heading(self.tag) #======================

        print(" DATA READING from ",self.filename)
        if ("N" in self.grb) and ("S" in self.grb) and ("B" in self.grb):
            print("Supplementary information is present")
            print(" grb.N==1 seen North only")
            print(" gbr.S==1 seen South only")
            print(" grb.B==1 seen on both")
            print("")
        print("+-------------------------- Flags ---------------------------+")
        print(" Flags:")
        print("   No visible flag, unvis           = ",self.unvis)
        print("   Iteration # from error code, 3s and 5s counts : ",
              self.niter, self.niter_3s, self.niter_5s)

        print("+----------------------- Statistics -------------------------+")
        print(f" {'Not visible':^15s} {'Fully analyzed':^15s} {'Aborted':^15s}")
        print(" {:^15d} {:^15d} {:^15d}"
          .format(len(self.grb[  self.grb.err == self.unvis]),
                  len(self.grb[  self.grb.err == self.niter]),
                  len(self.grb[ (self.grb.err != self.niter)
                              & (self.grb.err!=self.unvis) ])))
        print()
        print(" Raw statistics - max per site =",self.ngrb)
        print("  - total      : ",len(self.grb))
        print("  - analyzable : ",len(self.g_ana))
        print("  - North      : ",len(self.g_n))
        print("  - South      : ",len(self.g_s))
        print("  - Both sites : ",len(self.g_b),
              "-> total = ",len(self.g_n)+len(self.g_s)+len(self.g_b))

        print("  - North only : ",len(self.g_n0))
        print("  - South only : ",len(self.g_s0))

        print("+------------------------------------------------------------+")

    ###-------------------------------------------------------------------
    def slice_stat(self, gname,ax=None, **kwargs):
        """

        Slice number distribution

        Parameters
        ----------
        gname : string
            Subpopulation name.
        ax : Matplotlib axes, optional
            Current plots axis. The default is None.
        **kwargs : Dictionnary
            Extra arguments.

        Returns
        -------
        ax : Matplotlib axes
            Current matplotlib axis.

        """

        bins = np.linspace(0,30,31)
        pop = self.__dict__[gname]

        if ax==None:
            fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(5,5))

        n, bins, _ = ax.hist(pop.nt,
                             bins=bins,facecolor="none",edgecolor="black",
                             label=MyLabel(pop.nt))
        ax.legend()
        ax.set_xlabel("Number of slices in the GRB")
        ax.xaxis.set_major_locator(MultipleLocator(2.000))
        ax.set_title(gname)

        stamp(self.tag,axis=fig,where="bottom")

        return ax

    ###-------------------------------------------------------------------
    def trust_CL(self, gname, sigma=3):
        """
        Since the detection is at 3/5 sigma,
        It is expected that as a mean, the background fluctuate at 3/5 sigma
        in a fraction of the time. This function try to check it.

        Parameters
        ----------
        gname : string
            Subpopulation name.
        sigma : float, optional
            Siginificance thresold value. The default is 3.

        Returns
        -------
        None.

        """

        heading("Trust confidence level") #======================

        print(" Sub population ",gname)
        pop = self.__dict__[gname]

        print(" Err code < 0          = ",len(pop[pop.err<0]))
        nSL   = sum(pop.nt)

        if sigma == 3:
            nok = sum(pop.d3s)
            var = pop.d3s
            tag = r"$3\sigma$"
            print(" d3s code < 0          = ",len(pop[pop.d3s<0]))

        if sigma == 5:
            nok = sum(pop.d5s)
            var = pop.d5s
            tag = r"$5\sigma$"
            print(" d5s code < 0          = ",len(pop[pop.d5s<0]))

        print("Total numbers : ")
        print(f" nMC: trials  = {self.niter:10d}")
        print(f" nSL: slices  = {nSL:10d} (mean = {np.mean(pop.nt):3.2f})")
        print(f" -> MC trials (nMC x nSL) = {self.niter*nSL:10d}")

        p_value = 1-norm.cdf(sigma)
        nexpect = self.niter*nSL*p_value
        print()
        print(f" Probability to get a signal (p_value for {sigma:2d} sigma"\
              f"   = {100*p_value:7.6f}%)")
        print(f"    => {sigma:2d} sigma occurence in a no signal scenario"\
              f"    = {nexpect:10.2f}")
        print(f" Observed {sigma:2d} sigma trials in the data"\
              f"                = {nok:10.2f}")

    ###-------------------------------------------------------------------
    def sanity_check(self):
        """
        Perform some sanity checks on the data

        Returns
        -------
        None.

        """
        heading("Sanity checks") #======================

        # Min altitude
        _ , ax = plt.subplots(nrows=1, ncols=3, figsize=(10,2))
        for axi, gpop, tag in zip(ax,
                                  [self.g_n,self.g_s,self.g_b],
                                  ["North","South","Both"]):
            if tag != "Both":
                print(" Estimated min altitude in ",tag," :",min(gpop.altmx))
            axi.hist(gpop.altmx,bins=100,label=tag)
            axi.set_title(r"Altitude at max $\sigma$")
            axi.legend()
        plt.show()

        _ , ax0 = plt.subplots(nrows=1, ncols=2, figsize = (10,2))

        for axi, gpop, loc in zip(ax0,[self.g_n,self.g_s],["North","South"]):

            axi.hist(gpop[gpop.t3s>=0].t3s,bins=100,label=loc)
            delay = self.dtslew[loc]

            axi.axvline(x = delay.value,
                        color="red",ls=":",
                        label=str(delay.value)+" "+str(delay.unit))

            delay = delay+ self.dtswift
            axi.axvline(x= delay.value,
                       color="green",ls=":",
                       label=str(delay.value)+" "+str(delay.unit))

            axi.set_title(r"Time delay to $3 \sigma$  - "+loc)
            axi.set_xlim(xmin=0,xmax=10+1.3*delay.value)
            #axi.set_yscale("log")
            axi.legend()
            print(f" Estimated total delay in {loc:5s}")
            print(f" - From visibility start : {min(gpop[gpop.t1>=0].t1):5.1f}")
            print(f" - From 3s detection     : {min(gpop[gpop.t3s>=0].t3s):5.1f}")

    ###-------------------------------------------------------------------
    def negative_significance(self):
        """
        Negative significances
        The significance is returned by the WStatCountsStatistic function in
        gammapy, and has a sign which corresponds to the sign of the excess.
        As a consequence it can be negative. In case the simulation is done
        without fluctuation and one trial, the excess can not be negative,
        and all significances are zero (no excess count) or more.

        Returns
        -------
        None.

        """

        heading("Negative siginificances")

        for g,txt in zip([self.g_n, self.g_s, self.g_n0, self.g_s0, self.g_b],
                         ["North", "South", "North only","South only","Both"]):

            n_pos = len(g[g.sigmx>0])
            n0    = len(g[g.sigmx==0])
            n_neg = len(g[g.sigmx<0])

            if n_neg == 0:
                print(" This simualtion was probably without fluctuations")

            print(" {:10s} <0: {:<5d} ==0: {:<5d} >0: {:<5d}"
                  .format(txt,n_neg,n0,n_pos))

    ###-------------------------------------------------------------------
    def detectable_prompt(self):
        """
        Based on t90, prpt, avalaible in the most recent file.

        Note that counting for both sites has no sense since the delays are
        different on the two sites (the first has the intrinsic delay whereas
        the other is most likely having no delay if the time delay between the
        visibilities of the two sites is larger than the first site delay).

        Returns
        -------
        None.

        """
        heading(" GRB with prompt potentially visible")

        print(f" {'t90(s)>':8s} {'N':4s} {'S':4s}")
        for delay in [0,22,107]:
            prS = self.g_s[ (self.g_s.prpt==1) & (self.g_s.t90>delay) ]
            prN = self.g_n[ (self.g_n.prpt==1) & (self.g_n.t90>delay) ]
            n_S   = len(prS)
            n_N   = len(prN)
            print(f" {delay:<8.0f} {n_N:>4d} {n_S:>4d} ")

        print("In North : ",end="")
        print(" ".join([n[5:] for n  in prN.name]))
        print("In South : ",end="")
        print(" ".join([n[5:] for n  in prS.name]))

        fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(8,4),sharey=True)

        for ax, g, tag in zip(ax, [self.g_n,self.g_s],
                              ["North", "South"]):

            n, bins, _ = ax.hist(g.t90, bins=50,alpha=0.5)

            mask = (g.prpt==1)
            ax.hist(g[mask].t90,
                    bins=bins,alpha=0.5,
                    label=MyLabel(g[mask].t90,label="Prompt visible"))

            mask1 = (g.prpt==1) & (g[mask].t90>107)
            ax.hist(g[mask1].t90,
                    bins=bins,alpha=1,color="red",
                    label=MyLabel(g[mask1].t90,label="Prompt detectable >107s"))

            ax.axvline(x=107,color="green",ls="--")
            ax.set_title(tag)
            ax.set_yscale("log")
            ax.set_xlabel("$t_{90}$")
            ax.legend()

        plt.tight_layout()
        plt.show()

    ###-------------------------------------------------------------------
    def detection_above_sigma(self, sigma=5):
        """
        Display information on source detected above 5 sigma at 90%CL

        Parameters
        ----------
        sigma : float, optional
            Significance value. The default is 5.

        Returns
        -------
        None.

        """

        heading("Detection above a certain significance") #================

        for grb, tag in zip([self.g_n0, self.g_s0, self.g_b],
                            ["N only","S only", "Both"]):

            conf_level = 100*self.eff_lvl/self.niter

            print("###---------------------")
            print(f"###------- {tag:7s} -----")
            print("###---------------------")

            if sigma == 5:
                glist = grb[(grb.d5s >= self.eff_lvl) ]
                print(f" >= 5 sig. at {conf_level:2.0f}% CL :")

            elif sigma==3:
                glist = grb[(grb.d3s >= self.eff_lvl) ]
                print(f" >= 3 sig. at {conf_level:2.0f}% CL :")

            else:
                glist = grb[ (grb.d3s>=self.eff_lvl) & (grb.d5s<self.eff_lvl)]
                print(" >= 3 and <5 sigma at {conf_level:2.0f}% CL :")

            print(f"{' ':<5s} {'prt':>3s} {'z':>4s} {'t5s':>18s} "\
                  f"{'sigmax':>15s} {'tmax':>19s} {'date':>23s} "\
                  f"{'%3s':>3s} {'%5s':>3s}")

            for _, gg in glist.iterrows():
                t5s  = t_fmt(gg.t5s*u.s,digit=2)
                et5s = t_fmt((gg.et5s*u.s),digit=2).to(t5s.unit)

                tmx = t_fmt(gg.tmx*u.s,digit=2)
                etmx = t_fmt( (gg.etmx*u.s),digit=2).to(tmx.unit)

                t_trig = Time(gg.t_trig,format="mjd").isot
                if tag=="Both":
                    vis = str(int(self.g_n[self.g_n["name"]==gg["name"]].prpt))
                    vis+= "/"+str(int(self.g_s[self.g_s["name"]==gg["name"]].prpt))
                else:
                    vis = gg.prpt

                print(f"{gg['name'][5:]:<5s} {str(vis):3s} {gg.z:4.1f} " \
                      f"{t5s.value:5.1f} +/- {et5s.value:<5.1f} {t5s.unit:3s} "\
                      f"{gg.sigmx:>6.1f} +/- {gg.esigmx:>3.1f} "\
                      f"{tmx.value:5.1f} +/- {etmx.value:<5.1f} {tmx.unit:3s} "\
                      f"{t_trig:23s} {gg.d3s:3.0f} {gg.d5s:3.0f}")



###############################################################################
if __name__ == "__main__":

    # A standalone function to read a GRB and make various tests

    from pop_io import create_csv

    codefolder = "../../"
    sys.path.append(codefolder)

    nyears, file, _ = create_csv(file="parameter.yaml",debug=True)

    pop = Pop(filename=file, nyrs= nyears)
    pop.print()

    pop.sanity_check()
    pop.detection_above_sigma(sigma=5)
    pop.negative_significance()
    pop.detectable_prompt()
    pop.slice_stat("g_n")
    pop.trust_CL("g_tot")
