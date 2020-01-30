# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019

@author: Stolar
"""
__all__ = [
    'GetNameNoExtension',
    'OutputPrefixName',
]

import os

###--------------------------------------------------------------------
def GetNameNoExtension(name):
    pos_dot = 1000000
    for pos,ch in enumerate(name):
        #print(ch,' at ',pos)
        if (ch == '/'): pos_slash = pos
        if (ch =='.'): 
            if (pos < pos_dot): pos_dot = pos
    #print(" last / =",pos_slash)
    #print(" last . =",pos_dot)
    return name[pos_slash+1:pos_dot]


###----------------------------------------------------------------------------------------------------

def OutputPrefixName(grbname, chain, hem, zenith, azimuth, duration,  niter, NSB=False):
    outfile = grbname +"-"+chain+"-"+hem+"-"+zenith+"-"+azimuth+"-"+duration
    if (NSB):
        outfile = outfile+"-highNSB"
        
    outfile = outfile+"-"+str(niter)+"iter-"
    return outfile