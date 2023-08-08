# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 09:55:35 2022

A bunch of functions to manipulate text information.

@author: Stolar
"""
import os
from pathlib import Path
import astropy.units as u

__all__ = ["t_str", "t_fmt","heading","textcol",
           "warning","failure","success","highlight","banner", "Log"]

#------------------------------------------------------------------------------
def t_str(t, digit=2):

    """
    Transform a time quantity into a string (used for plot labels).

    Parameters
    ----------
    t : astropy Time Quantity
        Input time.
    digit : integer
        Significant digits to be printed.

    Returns
    -------
    String
        The formatted and rounded time as a string.

    """

    t = t_fmt(t)
    return str( round(t.value,digit)) +" "+ str(t.unit)

#------------------------------------------------------------------------------
def t_fmt(t, digit=None):
    """
    A utility to have reasonable duration format displayed.

    Parameters
    ----------
    t : astropy Time Quantity
        A time or duration

    Returns
    -------
    tobs : Quantity
        A time with an adapted unit and rounding for plotting

    """

    # Get and format livetimes
    t = t.to(u.s)
    if   t.value > 1.5*3600*24:
        t = t.to(u.d)
    elif t.value > 1.5*3600:
        t = t.to(u.h)
    elif t.value > 2*60:
        t = t.to(u.min)

    if digit is not None:
        return round(t.value,digit)*t.unit

    return t

###----------------------------------------------------------------------------
def heading(title, deco="-"):

    """
    Display a centered heading with 2 decorated lines above and below a title

    Parameters
    ----------
    title : String
        Text to be displayed.
    deco : String
        Character to be repeated for decoration.

    Returns
    -------
    None.

    """
    print(f"+{78*deco:78s}+")
    print(f"+{title:^78s}+")
    print(f"+{78*deco:78s}+\n")

###----------------------------------------------------------------------------
def textcol(text,t="black",b="white",s=None):

    """
    Change text color.
    See:

    * https://www.instructables.com/id/Printing-Colored-Text-in-Python-Without-Any-Module/
    * Add color to text in python :
      https://ozzmaker.com/add-colour-to-text-in-python/

    To make some of your text more readable, you can use ANSI escape codes to
    change the colour of the text output in your python program. A good use
    case for this is to highlight errors.
    The escape codes are entered right into the print statement.

    This table shows some of the available formats:

    .. code-block::

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
    String
        Colored text.

    """

    color = {"black" :"0",
             "red"   :"1",
             "green" :"2",
             "yellow":"3",
             "blue"  :"4",
             "purple":"5",
             "cyan"  :"6",
             "white" :"7"}
    style = {"normal"   :"0",
             "bold"     :"1",
             "underline":"2",
             "negative1":"3",
             "negative2":"5"}

    code = "\033["
    if s is not None:
        code = code + style[s] + ";"
    if t is not None:
        code = code + "3"+color[t] + ";"
    if b is not None:
        code = code + "4"+color[b] + ";"
    code = code[:-1] + "m"
    endcode = "\033[m"

    return code+text+endcode

###----------------------------------------------------------------------------
def warning(text, **kwarg):
    """
    Display a warning message

    Parameters
    ----------
    text : String
        A text to be displayed.
    **kwarg :
        Extra arguments.

    Returns
    -------
    None.

    """
    print(textcol(text,t="purple",b="white",s="bold"),**kwarg)

###----------------------------------------------------------------------------
def failure(text, **kwarg):
    """
    Display a failure message

    Parameters
    ----------
    text : String
        A text to be displayed.
    **kwarg :
        Extra arguments.

    Returns
    -------
    None.

    """
    print(textcol(text,t="red",s="bold"),**kwarg)

###----------------------------------------------------------------------------
def success(text,**kwarg):
    """
    Display a success message.

    Parameters
    ----------
    text : String
        A text to be displayed.
    **kwarg :
        Extra arguments.

    Returns
    -------
    None.

    """
    print(textcol(text,t="green", b="black",s="bold"),**kwarg)

###----------------------------------------------------------------------------
def highlight(text,**kwarg):
    """
    Display a failure message

    Parameters
    ----------
    text : String
        A text to be displayed.
    **kwarg :
        Extra arguments.

    Returns
    -------
    None.

    """
    print(textcol(text,s="bold"),**kwarg)

###----------------------------------------------------------------------------
def banner(text,**kwarg):
    """
    Display a banner message

    Parameters
    ----------
    text : String
        A text to be displayed.
    **kwarg :
        Extra arguments.

    Returns
    -------
    None.

    """
    print(textcol(text,t="black",b="yellow",s="bold"),**kwarg)

###----------------------------------------------------------------------------
class Log():

    """
    A class to manage a logbook information.

    """

    ###------------------------------------------------------------------------
    def __init__(self, log_name = None, talk = True):

        """
        Initialize the log book, with display either on Screen or in a log
        file, or both

        Parameters
        ----------
        log_name : Path, optional
            Log file path or name (String). The default is None.
        talk : boolean, optional
            If True, display on screen. The default is True.

        Returns
        -------
        None.

        """

        self.write = False # Disable print out to file
        self.talk  = talk  # Enbale/disable print out on screen
        self.log_file = None


        if log_name is not None:

            log_name = Path(log_name) # In case this would not be a Path

            try:
                self.log_file = open(log_name,'w')
            except IOError:
                print(f"Failure opening {log_name}: locked?")

            self.filename = log_name
            self.write = True
            print(f"log information to file {self.filename}")
        else:
            if talk is None:
                print(" Minimal information displayed")
            else:
                print(" Information displayed on Screen")

        # if not name.absolute().parent.exists(): # Check folder exists
        #     name.absolute().parent.mkdir(parents=True, exist_ok=True)

    ###------------------------------------------------------------------------
    def close(self, delete=False):

        if self.write:
            self.log_file.close()
            if delete:
                os.remove(self.filename)
        else:
            warning("No log file opened - cannot close")

    ###------------------------------------------------------------------------
    def prt(self, text, func=print, **kwarg):

        if self.talk:
            func(text,**kwarg)
        if self.write is not None:
            func(text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def warning(self,text,**kwarg):
        if self.talk:
            warning(text)
            # print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
        if self.write is not  None:
            print("*** WARNING *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def failure(self,text,**kwarg):
        if self.talk:
            failure(text, **kwarg)
        if self.write is not  None:
            print("*** FAILURE *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def success(self,text,**kwarg):
        if self.talk:
            print(textcol(text,t="green", b="black",s="bold"),**kwarg)
        if self.write is not  None:
            print("*** SUCCESS *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def highlight(self,text,**kwarg):
        if self.talk:
            print(textcol(text,s="bold"),**kwarg)
        if self.write is not  None:
            print("*** LOOK ! *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def banner(self,text,**kwarg):
        if self.talk:
            print(textcol(text,t="black",b="yellow",s="bold"),**kwarg)
        if self.write is not  None:
            print(" **** " + text + " *** ",**kwarg,file=self.log_file)

###----------------------------------------------------------------------------
if __name__ == "__main__":
    log = Log()
    log.failure("test failure",end="\n")
    log.failure("test 2 failure")
    log.close(delete = True)
