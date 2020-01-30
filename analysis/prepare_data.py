# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 18:17:09 2020

@author: Stolar
"""
import pandas as pd

# Productions where the altitude/azimuth what the one at the start of lasty observation slice cumulated
# Prod1 has only visible GRB, complicating the analysis
#file = "Prod1/"Pop_1000GRB_0dof_100iter.csv" # First production 
#file = "Prod2/Pop_1000GRB_0dof_100iter.csv"  # Second production with more columns : alt3s, alt5s, abort
#file = "Prod3_varslew/Pop_1000GRB_0dof_100iter.csv"  # Like the second one but with a varying slew time


# BUGGED Production where the alt-az is a tstart -fixed, the start of the observation
#file = "Prod-100sapcta/Pop_1-100GRB_0dof_100iter.csv"
#file = "Prod4_azcorr/Pop_1000GRB_0dof_100iter.csv"  # Like the second one but with a varying slew time

# First unbugged production
file = "Prod10_nobug/Pop_1-1000GRB_0dof_100iter.csv"
#file = "Prod10_nobug_noopt/Pop_1-1000GRB_0dof_100iter.csv"

filename = "../../output/"+file

grb = pd.read_csv(filename)
grb.columns

# Prepare subsamples
g_vis = grb[grb.abort>-1] # Visible
g_ana = grb[grb.abort==0] # Simulation complete (aborted value is at 10% of the iteration values)
g_abrt = grb[grb.abort>0]
#plt.hist(grb.abort,bins=100)

# initialize a few flags
i3s   = (grb.d3s >= 0.9) # 5s detected at 90% CL
i5s   = (grb.d5s >= 0.9) # 5s detected at 90% CL
south = grb.site=="South"
north = grb.site=="North"

g_3s  = g_ana[i3s]
g_5s  = g_ana[i5s]

print("Data where read from file : ",filename)
print("Data counts : ",len(grb))