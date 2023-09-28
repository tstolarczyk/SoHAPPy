# -*- coding: utf-8 -*-
"""
This script generates command lines for interactive running and batch
submissions of a series of simulation and analysis based on a large population.
It is written to use the
`slurm <https://slurm.schedmd.com/documentation.html>`_ batch system.

Created on Mon Mar  6 09:29:47 2023

@author: Stolar
"""
import sys
import argparse
import platform
from pathlib import Path

from skygen import Skies
from configuration import Configuration
from utilities import subset_ids
from niceprint import heading

__all__ = ["slurm_cmd", "decode_command_line"]

###----------------------------------------------------------------------------
def slurm_cmd(job_cmd,
              nproc       = 1,
              mem_per_cpu = "2G",
              duration    = "00:30:00",
              name        = "SoHAPPy"):
    """
    Note: When run on Spyder, a 1000 source production occupy 1.6 Gb of memory on
    a classical laptop, thus justifying the `mem_per_cpu` default value.

    Parameters
    ----------
    job_cmd : string
        The command to be submitted.
    nproc : integer, optional
        Number of processor per task (`srun` option). The default is 1.
    mem_per_cpu : string, optional
        Memory per cpu core. The default is "1G".
    duration : string, optional
        Estimated duration (optimise the choice of the bach queue).
        The default is "00:30:00" (30 minutes).
    name : string, optional
        Job name, will appear using the `squeue` command (list of running jobs).
        The default is "SoHAPPy".

    Returns
    -------
    cmd : string
        The `slurm` batch command.

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
    Decode the current command line.

    Returns
    -------
    parszer object
        List of arguments in the command line.

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

    parser.set_defaults(batch=True)

    parser.add_argument('--batch',
                        dest='batch',
                        action='store_true',
                        help = "Create batch commands")

    parser.add_argument('--nobatch',
                        dest='batch',
                        action='store_false',
                        help = "Does not create batch commands")

    return parser.parse_known_args()  # Separate known arguments from others

###############################################################################
if __name__ == '__main__':

    # print(sys.argv)

    # If no argument given, run an example - Used for tests
    if len(sys.argv[1:]) <= 0:
        heading("Execute examples")
        # sys.argv=["", "-C", r"./skygen.py",
        #               "-V", "strictmoonveto",
        #               "-P", "10",
        #               "-S", "2",
        #               "--nobatch",
        #               "-c", r"data/config_ref.yaml"]
        sys.argv=["", "-C", "SoHAPPy.py","-V","strictmoonveto","-P","10",
                      "-S", "3", "--nobatch", "-d", 0]

    # Get command line arguments and debrief
    args, extra_args = decode_command_line()
    cmd = ' '.join(f'--{k} {v}' for k, v in vars(args).items())

    # Debrief
    print(" Equivalent command line:",cmd)
    print(f" Arguments to be passed to {args.Code:} : {extra_args:}")

    # Prepare job slicing
    dsets = subset_ids(args.Pop, args.Sets)
    print(" Slicing:",dsets)

    print()
    print(" Generated commands: ")

    # Open output script file
    outname = Path(args.Code).stem+"_"+str(args.Pop)+"_"+str(args.Sets)
    extname = ".ps1" if platform.system() == "Windows" else ".sh"
    if args.batch:
        fbatch  = open("batch_"+outname+extname,"w")
    else:
        finter  = open("interactive_"+outname+extname,"w")

    # Special action to transform the config filename into a resolved path name
    if "-c" in extra_args:
        idx = extra_args.index("-c")
        extra_args[idx+1] = str(Path(extra_args[idx+1]).resolve())
    elif "--config" in extra_args:
        idx = extra_args.index("--config")
        extra_args[idx+1] = str(Path(extra_args[idx+1]).resolve())

    # Loop over sets
    for [id1, id2] in dsets:

        # Pass extra arguments to external code
        sys.argv = [args.Code]   + extra_args + \
                 [ "-f", str(id1), "-N", str(id2-id1+1) ]

        if args.Code.find("skygen") != -1:
            sky = Skies.command_line()
            job_cmd = sky.cmd_line

        elif args.Code.find("SoHAPPy") != -1:
            cf = Configuration.command_line()
            job_cmd = cf.cmd_line

        else:
            print("generator.py : ",args.Code, "Not handled")

        # Create batch command if requested - display and write commands
        if args.batch:
            batch_cmd = slurm_cmd(job_cmd, name="SoHAPPy")
            print(" > ",batch_cmd)
            print(batch_cmd, file = fbatch)
        else:
            print(" >",job_cmd)
            print("python "+job_cmd, file = finter)

    # Close the script file
    if args.batch:
        fbatch.close()
    else:
        finter.close()
