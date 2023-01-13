# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 09:55:35 2022

@author: Stolar
"""

#------------------------------------------------------------------------------
def t_str(t, digit=2):
    """
    Transform a time quantity into a string (used for plot labels)

    Parameters
    ----------
    t : astropy Time Quantity
        Input time.
    digit : integer
        Siginificant digits to be printed.
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
    import astropy.units as u

    # Get and format livetimes
    t = t.to(u.s)
    if   t.value > 1.5*3600*24: t = t.to(u.d)
    elif t.value > 1.5*3600:    t = t.to(u.h)
    elif t.value > 2*60:        t = t.to(u.min)

    if digit != None: return round(t.value,digit)*t.unit
    else:             return t

###----------------------------------------------------------------------------
def heading(title):
    print(f"+{78*'-':78s}+")
    print(f"+{title:^78s}+")
    print(f"+{78*'-':78s}+\n")

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
    import sys
    ###------------------------------------------------------------------------
    def __init__(self, name="default.log",talk=True):

        from pathlib import Path
        name = Path(name) # In case this would not be a path

        if not name.absolute().parent.exists(): # Check folder exists
            name.absolute().parent.mkdir(parents=True, exist_ok=True)

        try:
            self.log_file = open(name,'w')
        except IOError:
            print("Failure opening {}: locked".format(name))

        self.name     = name
        self.talk     = talk
        return

    ###------------------------------------------------------------------------
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

    ###------------------------------------------------------------------------
    def warning(self,text,**kwarg):
        if (self.talk):
            warning(text)
            # print(textcol(text,t="purple",b="white",s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** WARNING *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def failure(self,text,out=sys.stdout,**kwarg):
        if (self.talk):
            failure(text, **kwarg)
        if (self.log_file != None):
            print("*** FAILURE *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def success(self,text,out=sys.stdout,**kwarg):
        if (self.talk):
            print(textcol(text,t="green", b="black",s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** SUCCESS *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
    def highlight(self,text,out=sys.stdout,**kwarg):
        if (self.talk):
            print(textcol(text,s="bold"),**kwarg)
        if (self.log_file != None):
            print("*** LOOK ! *** "+text,**kwarg,file=self.log_file)

    ###------------------------------------------------------------------------
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