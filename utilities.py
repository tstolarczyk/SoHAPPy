# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019

@author: Stolar
"""
import sys
import numpy as np

###----------------------------------------------------------------------------
# A label for plots and histograms with statistics
def MyLabel(var,label="",stat="std"):
    if (len(label)!=0): label = label+"\n"
    legend="bad option"
    if (stat=="std"):
        legend = label  \
                + "$n$ : {:d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$\sigma$ : "+" {:5.3f}".format(np.std(var))
    elif (stat=="med"):
        legend = label  \
                + "$n$ : {:d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$Med.$ : "+" {:5.3f}".format(np.median(var))
    return legend

###----------------------------------------------------------------------------
# Annotation on the side of any plot referred from the axis
def stamp(ax,text):
    ax.text(x=1,
            y=0.,
            fontsize=8,
            horizontalalignment='left',
            verticalalignment='bottom',
            transform=ax.transAxes,
            rotation=270,
            s=text)
    return
###----------------------------------------------------------------------------
# Colored text
# see https://www.instructables.com/id/Printing-Colored-Text-in-Python-Without-Any-Module/
def textcol(text,t="black",b="white",s=None):
    color = {"black":"0",
             "red":"1",
             "green":"2",
             "yellow":"3",
             "blue":"4",
             "purple":"5",
             "cyan":"6",
             "white":"7"}
    style = {"normal":"0",
             "bold":"1",
             "underline":"2",
             "negative1":"3",
             "negative2":"5"}
    
    code = "\033["
    if (s != None): code = code + style[s] + ";"
    if  (t != None): code = code + "3"+color[t] + ";"
    if (b != None): code = code + "4"+color[b] + ";"
    code = code[:-1] + "m"
    endcode = "\033[m"
    return code+text+endcode
###
#def failure(text): print(textcol(text,t="white",b="red",s="bold"))
#def warning(text): print(textcol(text,t="white",b="purple",s="bold"))
#def success(text): print(textcol(text,t="black",b="green",s="bold"))
def failure(text,out=sys.stdout,**kwarg): 
    print(textcol(text,t="red",s="bold"),**kwarg)
def warning(text,out=sys.stdout,**kwarg): 
    print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
def success(text,out=sys.stdout,**kwarg): 
    print(textcol(text,t="green", b="black",s="bold"),**kwarg)
def highlight(text,out=sys.stdout,**kwarg):
    print(textcol(text,s="bold"),**kwarg)
def banner(text,out=sys.stdout,**kwarg):
    print(textcol(text,t="black",b="yellow",s="bold"),**kwarg)
    
###----------------------------------------------------------------------------
def backup_file(folder=None, dest=None, dbg=False):
    """
    Copy a file to a result folder
    If it already exists, make a backup with the date
    """
    import os
    import shutil
    import datetime

    # Create result folder if not exisitng
    if (not os.path.exists(folder)):
            if (dbg): print(" *** Creating output folder ",folder)
            os.makedirs(folder)
            
    conf_filename = folder+"/"+dest
    
    if (os.path.exists(conf_filename)):
        nw = datetime.datetime.now()
        newname = conf_filename + "_" \
                                + str(nw.hour) \
                                + str(nw.minute) \
                                + str(nw.second)
        
        os.rename(conf_filename, newname)
        if (dbg): print("     --- config.txt exists, renamed to ",newname)
        
    shutil.copy('ana_config.py',conf_filename)

    return conf_filename

###----------------------------------------------------------------------------
class Log():
    def __init__(self, name="default.log",talk=False):
        
        self.log_file = open(name,'w')  
        self.talk     = talk
        return

    def prt(self, text,**kwarg):
        if (self.talk): 
            print(text,**kwarg)
        if (self.log_file != None): print(text,**kwarg,file=self.log_file)
        
    def close(self):
        self.log_file.close()
        return
    
    def warning(self,text,**kwarg):
        if (self.talk): 
            print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** WARNING *** "+text,**kwarg,file=self.log_file)
    
    def failure(self,text,out=sys.stdout,**kwarg): 
        if (self.talk): 
            print(textcol(text,t="red",s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** FAILURE *** "+text,**kwarg,file=self.log_file)
        
    def success(self,text,out=sys.stdout,**kwarg): 
        if (self.talk): 
            print(textcol(text,t="green", b="black",s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** SUCCESS *** "+text,**kwarg,file=self.log_file)
    
    def highlight(self,text,out=sys.stdout,**kwarg):
        if (self.talk): 
            print(textcol(text,s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** LOOK ! *** "+text,**kwarg,file=self.log_file)
    
    def banner(self,text,out=sys.stdout,**kwarg):
        if (self.talk): 
            print(textcol(text,t="black",b="yellow",s="bold"),**kwarg)       
        if (self.log_file != None):
            print(" **** " + text + " *** ",**kwarg,file=self.log_file)  
###----------------------------------------------------------------------------
if __name__ == "__main__":
    failure("test failure",end="\n")
    failure("test 2 failure")