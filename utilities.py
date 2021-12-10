# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:41:11 2019

@author: Stolar
"""
import sys
import numpy as np
import astropy.units as u

###----------------------------------------------------------------------------
# A label for plots and histograms with statistics
def MyLabel(var,label="",stat="std"):
    """
    Add extra statistical information (counts, dispersions) to a classical 
    matplotlib label. 

    Parameters
    ----------
    var : A list or numpy array
        The plotted variable.
    label : String, optional
        The classical label text. The default is "".
    stat : String, optional
        A keyword defining the extra information to be displayed (on top 
        of counts) : "std" for standard deviation, "med" for median, None to 
        have counts only. The default is "std".

    Returns
    -------
    legend : String
        The modfied label text.
    
    """
    
    if (len(label)!=0): label = label+"\n"
    
    legend="bad option"
    
    if stat.findstr("std") !=-1:
        legend = label  \
                + "$n$ : {:d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$\sigma$ : "+" {:5.3f}".format(np.std(var))
    elif stat.findstr("med") !=-1 :
        legend = label  \
                + "$n$ : {:d} \n".format(len(var)) \
                + r"$\bar{n}$ : "+"{:5.3f}\n".format(np.mean(var)) \
                + r"$Med.$ : "+" {:5.3f}".format(np.median(var))
    elif stat == None:
        legend =  label  \
                + "$n$ : {:d}".format(len(var))
                
    return legend

 ###----------------------------------------------------------------------------
def single_legend(ax,**kwargs):
    """
    Remove duplicated labels in legend (e.g. a vertical and an horizontal 
    lines defining and intersection havingthe same labels).
    
    Parameters
    ----------
    ax : matplotlib.axes
        Current axis.

    Returns
    -------
    None.

    """
    from collections import OrderedDict
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),**kwargs)
    
    return

###----------------------------------------------------------------------------
def stamp(text, axis=None,
          where="right",x=None, y=None, rotation=0, **kwargs):
    """
    Annotation on the side of any plot referred from the axis, including
    the gammapy version.
    
    Parameters
    ----------
    text : String
        Text to be displayed
    where : String, optional
        Position of the text with respect to the axis. The default is "right".
    x : float, optional
        Text x position, supersedes the where variable. The default is None.
    y : float, optional
        Text x position, supersedes the where variable. The default is None.
    rotation : float, optional
        Text rotation - If not given use default. The default is None.
     axis : matplotlib.axes, optional
         Current plot axis. The default is None.
    **kwargs : 
        Any additionnal arguments for matplotlib.axes.text

    Returns
    -------
    None.
    """

    import gammapy
    text = text + " - " + gammapy.__version__
    
    if x==None or y== None:
        if   where =="right":  
            (x,y) = (  1, 0.5)
            rotation = 270
        elif where =="left":   
            (x,y) = (  0, 0.5)
            rotation = 90
        elif where =="top":    
            (x,y) = (0.5,   1)
            rotation = 0
        elif where =="bottom": 
            (x,y) = (0.5,   0)
            rotation = 0
    
    axis.text(x=x,y=y,s=text,
              horizontalalignment='center',
              verticalalignment="center",
              rotation=rotation)
    return
###----------------------------------------------------------------------------
def projected_scatter(xsize=12, ysize=8, 
                      left=0.1, width=0.7, bottom=0.1, height=0.7, 
                      spacing=0.02):
    """
    Matplotlib template to display a scatter plot and the horizontal and 
    vertical projections of the data. Adapted from :
    https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html#sphx-glr-gallery-lines-bars-and-markers-scatter-hist-py

    Parameters
    ----------
    xsize : float, optional
        Figure horizontal size. The default is 12.
    ysize : float, optional
        Figure vertical size. The default is 8.
    left : float, optional
        Scatter plot left position. The default is 0.1.
    width : float, optional
        Height of the histogram handling the y-projection. The default is 0.7.
    bottom : float, optional
        Scatter plot bottom position. The default is 0.1.
    height : float, optional
        Height of the histogram handling the x-projection. The default is 0.7.
    spacing : TYPE, optional
        Space between the subplots. The default is 0.02.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    ax : matplotlib.axes
        Scatter plot axis.
    axh : matplotlib.axes
        Horizontal projection axis.
    axv : matplotlib.axes
        Vertical projection axis.

    """
    import matplotlib
    
    rect_scatter = [left, bottom, width, height]
    rect_histx   = [left, bottom + height + spacing, width, 0.2]
    rect_histy   = [left + width + spacing, bottom, 0.2, height]

    fig = matplotlib.pyplot.figure(figsize=(xsize, ysize))
    ax  = fig.add_axes(rect_scatter)
    axh = fig.add_axes(rect_histx, sharex=ax)
    axv = fig.add_axes(rect_histy, sharey=ax)

    axh.tick_params(axis="x", labelbottom=False)
    axv.tick_params(axis="y", labelleft=False)    
    
    return fig, ax, axh, axv
#------------------------------------------------------------------------------
def t_str(t, digit=2):
    """
    Transform a time quantity into a string (used for plt labels)

    Parameters
    ----------
    t : Tilme Quantity
        Input time.
    digit : integer
        Siginificant digit to be printed.
    Returns
    -------
    String
        The formatted and rounded time as a string.

    """
    
    t = t_fmt(t)
    return str( round(t.value,digit)) +" "+ str(t.unit)
#------------------------------------------------------------------------------
def t_fmt(t):
    """
    A utility to have reasonable duration format displayed.

    Parameters
    ----------
    t : Quantity
        A time or duration

    Returns
    -------
    tobs : Quantity
        A time with an adapted unit and rounding for plotting

    """
    # Get and format livetimes
    t = t.to(u.s)
    if t.value > 1.5*3600*24:
        t = t.to(u.d)
    elif t.value > 1.5*3600:
        t = t.to(u.h)
    elif t.value > 2*60:
        t = t.to(u.min)

    return t

###----------------------------------------------------------------------------
def ColorMap(threshold,maxval):
    """
    Create a colormap based on a threshold and a maximal value 

    Parameters
    ----------
    threshold : float
        Starting values.
    maxval : float
        End value.

    Returns
    -------
    mymap : matplotlib.colors.LinearSegmentedColormap
        A matplotlib color map.

    """
    
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    color_cut = int(256*threshold/maxval)
    print("color_cut=",color_cut)
    colors1  = plt.cm.Greys(np.linspace(0, 1, color_cut))
    colors2  = plt.cm.plasma(np.linspace(0, 1, 256-color_cut))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    return mymap
###----------------------------------------------------------------------------
def backup_file(filename,folder=None, dbg=False):
    """
    Copy a file to a result folder.
    If it already exists, make a backup with the date.

    Parameters
    ----------
    file : TYPE
        DESCRIPTION.
    folder : String, optional
        Output folder. The default is None.
    dbg : Boolean, optional
        If True, display messages. The default is False.

    Returns
    -------
    None.

    """
    """

    """
    import os
    import shutil
    import datetime

    # Create result folder if not exisitng
    if (not os.path.exists(folder)):
            if (dbg): print(" *** Creating output folder ",folder)
            os.makedirs(folder)

    output_file = folder+"/"+filename

    if os.path.exists(output_file):
        nw = datetime.datetime.now()
        output_file = output_file + "_" \
                                  + str(nw.hour) \
                                  + str(nw.minute) \
                                  + str(nw.second)

    shutil.copy(filename, output_file)
    if (dbg): print("   ----",filename," copied to ",output_file())

    return 

###----------------------------------------------------------------------------
def textcol(text,t="black",b="white",s=None):
    """
    Change text color.
    See:
    https://www.instructables.com/id/Printing-Colored-Text-in-Python-Without-Any-Module/
    Add color to text in python : https://ozzmaker.com/add-colour-to-text-in-python/
    
    To make some of your text more readable, you can use ANSI escape codes to 
    change the colour of the text output in your python program. A good use 
    case for this is to highlight errors.
    The escape codes are entered right into the print statement.
    
    print("\033[1;32;40m Bright Green \n")
    
    The above ANSI escape code will set the text colour to bright green. 
    The format is; \033[ Escape code, this is always the same 1 = Style, 
    1 for normal. 32 = Text colour, 32 for bright green. 
    40m = Background colour, 40 is for black.
    
    This table shows some of the available formats; 
    TEXT COLOR CODE TEXT STYLE CODE BACKGROUND COLOR CODE 
    Black 30 No effect 0 Black 40 Red 31 Bold 1 Red 41 Green 32 Underline 
    2 Green 42 Yellow 33 Negative1 3 Yellow 43 Blue 34 Negative2 5 Blue 
    44 Purple 35 Purple 45 Cyan 36 Cyan 46 White 37 White 47

    Parameters
    ----------
    text : String
        Text to be displayed.
    t : String, optional
        Text color. The default is "black".
    b : String, optional
        Text background color. The default is "white".
    s : String, optional
        Text style keyword. The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

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
    if (t != None): code = code + "3"+color[t] + ";"
    if (b != None): code = code + "4"+color[b] + ";"
    code = code[:-1] + "m"
    endcode = "\033[m"
    return code+text+endcode

###----------------------------------------------------------------------------
def warning(text, **kwarg):
    print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
    return
def failure(text, **kwarg):
    print(textcol(text,t="red",s="bold"),**kwarg)
    return
def success(text,**kwarg):
    print(textcol(text,t="green", b="black",s="bold"),**kwarg)
    return
def highlight(text,**kwarg):
    print(textcol(text,s="bold"),**kwarg)
def banner(text,**kwarg):
    print(textcol(text,t="black",b="yellow",s="bold"),**kwarg)    
    return
    
###----------------------------------------------------------------------------

class Log():
    """
    A class to manage a logbook information
    
    """
    def __init__(self, name="default.log",talk=True):

        self.log_file = open(name,'w')
        self.name     = name
        self.talk     = talk
        return

    def prt(self, text,**kwarg):
        if (self.talk):
            print(text,**kwarg)
        if (self.log_file != None): print(text,**kwarg,file=self.log_file)

    def close(self, delete=False):
        self.log_file.close()
        if delete:
            import os
            os.remove(self.name)
        return

    def warning(self,text,**kwarg):
        if (self.talk):
            warning(text)
            # print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** WARNING *** "+text,**kwarg,file=self.log_file)

    def failure(self,text,out=sys.stdout,**kwarg):
        if (self.talk):
            failure(text, **kwarg)
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
    log = Log()
    log.failure("test failure",end="\n")
    log.failure("test 2 failure")
    log.close(delete = True)
