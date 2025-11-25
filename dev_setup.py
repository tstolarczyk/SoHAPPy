# -*- coding: utf-8 -*-
"""
Running this script sets or unsets the default environment variables.

To be adaptaed to each installation.

Created on Wed Nov 13 14:37:26 2024

@author: Stolar

To be run to setup the developement environment.
"""

import os
import sys
from pathlib import Path
import seaborn as sns

style_list = ["poster", "talk", "notebook", "paper"]


# -----------------------------------------------------------------------------
def set_env(style="notebook"):
    """Set input/output folders from the environement varaiables and style."""
    # Bigger texts and labels
    if style in style_list:
        sns.set_context(style)
    else:
        sys.exit(f"dev_steup: Cannot set style. Use {style_list}")

    if "HAPPY_IN" in os.environ.keys():
        infolder = Path(os.environ["HAPPY_IN"])
    else:
        sys.exit("The HAPPY_IN environment variable should be defined")

    if "HAPPY_OUT" in os.environ.keys():
        resfolder = Path(os.environ["HAPPY_OUT"])
    else:
        sys.exit("The HAPPY_OUT environment variable should be defined")

    return infolder, resfolder


# #############################################################################
if __name__ == "__main__":

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
