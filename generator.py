# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 09:29:47 2023

@author: Stolar
"""
import sys
import argparse
from skygen import Skies
from configuration import Configuration
from utilities import subset_ids
from niceprint import heading

__all__ = []
###----------------------------------------------------------------------------
def slurm_cmd(job_cmd,
                nproc=1,
                mem_per_cpu="2G",
                duration="00:00:30",
                name="SoHAPPy"):
    """
    When run on Spyder, a 1000 source prodcution occupy 1.6 Gb of memory on
    a classical laptop.

    Parameters
    ----------
    job_cmd : TYPE
        DESCRIPTION.
    nproc : TYPE, optional
        DESCRIPTION. The default is 1.
    mem_per_cpu : TYPE, optional
        DESCRIPTION. The default is "1G".
    duration : TYPE, optional
        DESCRIPTION. The default is "00:00:30".
    name : TYPE, optional
        DESCRIPTION. The default is "SoHAPPy".

    Returns
    -------
    cmd : TYPE
        DESCRIPTION.

    """

    cmd  = "sbatch"
    cmd += " -c "  + str(nproc)
    cmd += " --mem-per-cpu "+mem_per_cpu
    cmd += " -t "+ duration
    cmd += " -J "+ name
    cmd += " --wrap '"+ job_cmd+"'"

    return cmd

###----------------------------------------------------------------------------
def decode_command_line():
    """


    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    parser = argparse.ArgumentParser(description="Generate batch scripts for SoHAPPy",
                                     epilog="---")
    parser.add_argument('-C', '--Code',
                        help ="Code to be run",
                        required=True)
    parser.add_argument('-P', '--Pop',
                        help ="Population statistics",
                        default=10,
                        type=int)
    parser.add_argument('-S', '--Sets',
                        help ="Number of sets",
                        default=2,
                        type=int)
    parser.add_argument('-b', '--batch',
                        help ="Batch flag (True, generates batch commands)",
                        default=True,
                        type=bool)

    return parser.parse_known_args()  # Separate known arguments from others

###############################################################################
if __name__ == '__main__':

    print(sys.argv)

    if len(sys.argv[1:]) <= 0:
        heading("Execute examples")
        sys.argv=["", "-C", "skygen.py","-V","strictmoonveto","-n","100", "-b", "True"]


    # Get command line arguments
    args, extra_args = decode_command_line()

    # Debrief
    cmd = ' '.join(f'--{k} {v}' for k, v in vars(args).items())
    print(" Equivalent command line:",cmd)
    print(" Extra arguments        :",extra_args)
    print(" To be passed to        :",args.Code)

    # Prepare job slicing
    dsets = subset_ids(args.Pop, args.Sets)
    print(" Batch commands to be defined for these sets:")
    print(dsets)

    print()
    print(" Generated commands: ")

    # Pass extra arguments to external code
    if args.Code == "skygen.py":
        for [id1, id2] in dsets:

            # Output folder name
            # resdir = "batch_"+str(dset[0])+"_"+str(dset[1])

            sys.argv = [args.Code]   + extra_args + \
                     [ "-f", str(id1), "-N", str(id2-id1+1) ]
            sky = Skies.from_command_line()
            job_cmd = sky.cmd_line
            print(" >",job_cmd)

            if args.batch:
                batch_cmd = slurm_cmd(job_cmd, name="tbd")
                print(" >> ",batch_cmd)

    elif args.Code == "SoHAPPy.py":
         cf = Configuration.build(extra_args)
    else:
        print(args.Code, "Not handled")

