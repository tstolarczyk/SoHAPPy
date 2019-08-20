# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019

@author: Stolar
"""
__all__ = [
    'GetNameNoExtension',
    'GetIrfFile',
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

###--------------------------------------------------------------------

def GetIrfFile(irf_folder, chain, hemisphere, zenith, azimuth, duration, NSB=False):
    # This find the correct IRF from a series of keyword
    # Assumption
    # 1. The files atre stored in the irf_folder
    # 2. Subfolders contain the file for a certain analysis chain. The subfolder name contain any of authorised keywords
    # 3. North and South hemispheres are stored separately
    # high NSB file are stored in a sub folder
    
    # Remodify values totake into account current IRF faile format
    zenith   = zenith+"deg"

    irf_file = "Error"
    found = False
    
    # Check that IRF foler exists
    if (not os.path.exists(irf_folder)):
        print(" *** GetIrfFile ERROR : IRF folder ",irf_folder," not found ***")
        return found, irf_folder, irf_file
    
    # Get subfolder list, check presence of chain
    chainlist = os.listdir(irf_folder)
    # print(chainlist)
    
    for folder in chainlist:
        
        # print(" scan ",folder," res=",folder.find(chain))
        if (folder.find(chain)+1): # if not found return -1, i.e +1 = zero
            found = True
            irf_file = chain + " OK"
            # Check hemisphere
            irf_folder = irf_folder+'/'+folder+'/'+hemisphere
            if (not os.path.exists(irf_folder)):
                print(" *** GetIrfFile ERROR : ", hemisphere, " not found ***")
                return False, irf_folder, irf_file
            irf_file = irf_file + " " +hemisphere + " OK"
            # check NSB is required
            if (NSB): # If searching for high NSB file
                irf_folder = irf_folder + '/NSB'
                if (not os.path.exists(irf_folder)):
                    print(" *** GetIrfFile ERROR : no high NSB file found ***")
                    return False,irf_folder, irf_file

    # Redundant
    if (not found): 
        print(" *** GetIrfFile ERROR : ", chain, " not found ***")
        return found, irf_folder, irf_file


    # At this stage, the required folder exists - try to find a suite IRF file
    filelist = os.listdir(irf_folder)
    #print(filelist,"\n")
    nfound = 0
    for file in filelist:
              if (file.find(zenith)+1 and file.find(azimuth)+1 and file.find(duration)+1):
                nfound+=1
                irf_file = file

    if (nfound != 1):
        found = False
        if (nfound == 0): irf_file = " *** GetIrfFile ERROR : No file found for "
        if (nfound > 1): irf_file = " *** GetIrfFile ERROR : More than one file found for"   
        irf_file = irf_file \
        + chain +" " \
        + hemisphere + " " \
        + zenith + " " \
        + azimuth + " " \
        + duration + " ***"
    else:
        found = True
        
    return found, irf_folder, irf_file

###----------------------------------------------------------------------------------------------------

def OutputPrefixName(grbname, chain, hem, zenith, azimuth, duration,  niter, NSB=False):
    outfile = grbname +"-"+chain+"-"+hem+"-"+zenith+"-"+azimuth+"-"+duration
    if (NSB):
        outfile = outfile+"-highNSB"
        
    outfile = outfile+"-"+str(niter)+"iter-"
    return outfile