# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019

@author: Stolar
"""

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
def textcol(text,t="black",b="white",s="None"):
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
    if (s != "None"): code = code + style[s] + ";"
    if  (t != "None"): code = code + "3"+color[t] + ";"
    if (b != "None"): code = code + "4"+color[b] + ";"
    code = code[:-1] + "m"
    endcode = "\033[m"
    return code+text+endcode
###
#def failure(text): print(textcol(text,t="white",b="red",s="bold"))
#def warning(text): print(textcol(text,t="white",b="purple",s="bold"))
#def success(text): print(textcol(text,t="black",b="green",s="bold"))
def failure(text,**kwarg): 
    print(textcol(text,t="red",s="bold"),**kwarg)
def warning(text,**kwarg): 
    print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
def success(text,**kwarg): 
    print(textcol(text,t="green", b="black",s="bold"),**kwarg)
def highlight(text,**kwarg):
    print(textcol(text,s="bold"),**kwarg)
def banner(text,**kwarg):
    print(textcol(text,t="black",b="yellow",s="bold"),**kwarg)
    
###----------------------------------------------------------------------------
if __name__ == "__main__":
    failure("test failure",end="\n")
    failure("test 2 failure")