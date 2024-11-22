# -*- coding: utf-8 -*-
"""
Running this script sets or unsets the default environment variables. To be
adaptaed to each installation.

Created on Wed Nov 13 14:37:26 2024

@author: Stolar

To be run to setup the developement environment
"""

import os


if "HAPPY_IN" in os.environ:
    del os.environ["HAPPY_IN"]
    print("HAPPY_IN deleted")
else:
    os.environ["HAPPY_IN"] = r"D:\\CTAO\SoHAPPy\input"
    print("Input folder  = ", os.environ["HAPPY_IN"])


if "HAPPY_OUT" in os.environ:
    del os.environ["HAPPY_OUT"]
    print("HAPPY_OUT deleted")
else:
    os.environ["HAPPY_OUT"] = r"D:\\CTAO\SoHAPPy\output"
    print("Output folder = ", os.environ["HAPPY_OUT"])
