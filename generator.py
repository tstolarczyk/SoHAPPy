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


# ##---------------------------------------------------------------------------
def slurm_cmd(job_cmd,
              nproc=1,
              mem_per_cpu="2G",
              duration="30:00:00",
              name="SoHAPPy"):
    """
    Generate batch command.

    Note: When run on Spyder, a 1000 source production occupy 1.6 Gb of memory
    on a classical laptop, thus justifying the `mem_per_cpu` default value.

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
        The default is "10:00:00" (10 hours).
        Processing 300 GRBs in mode "infinite nights" requires 7h on the IRFU
        Feynman cluster if the ligh curves are limited to 3 days.
    name : string, optional
        Job name, will appear using the `squeue` command (list of running
        jobs). The default is "SoHAPPy".

    Returns
    -------
    cmd : string
        The `slurm` batch command.

    """
    cmd = "sbatch"
    cmd += " -c " + str(nproc)
    cmd += " --mem-per-cpu " + mem_per_cpu
    cmd += " -t " + duration
    cmd += " -J " + name
    cmd += " --wrap '" + job_cmd+"'"

    return cmd


# ##---------------------------------------------------------------------------
def decode_command_line():
    """
    Decode the current command line.

    Returns
    -------
    parser object
        List of arguments in the command line.

    """
    parser = argparse.ArgumentParser(description="Generate scripts",
                                     epilog="---")
    parser.add_argument('-C', '--Code',
                        help="Code to be run",
                        required=True)

    parser.add_argument('-P', '--Pop',
                        help="Population statistics",
                        default=10,
                        type=int)
    parser.add_argument('-S', '--Sets',
                        help="Number of sets",
                        default=2,
                        type=int)

    parser.add_argument('-B', '--Begin',
                        help="Starting identifier",
                        default=1,
                        type=int)

    parser.add_argument('-j', '--json',
                        help="Json file folder",
                        default=None)

    parser.set_defaults(batch=True)

    parser.add_argument('--batch',
                        dest='batch',
                        action='store_true',
                        help="Create batch commands")

    parser.add_argument('--nobatch',
                        dest='batch',
                        action='store_false',
                        help="Does not create batch commands")

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
        # sys.argv = ["", "-C", "SoHAPPy.py",
        #             "--config", "MyConfigs/prod5_std/config_prod5_std.yaml",
        #             "-P", "10000",
        #             "-S", "5", "--nobatch", "-d", "0"]
        sys.argv = ["", "-C", "SoHAPPy.py",
                   #"--config", r"MyConfigs/prod5_std/config_prod5_std.yaml",
                    "-P", "2000",
                    "-S", "1", "--nobatch", "-d", "0",
                    "--json", "data/det3s/test_detected"]

    # Get command line arguments and debrief
    args, extra_args = decode_command_line()
    cmd = ' '.join(f'--{k} {v}' for k, v in vars(args).items())

    # Debrief
    print(" Equivalent command line:", cmd)
    print(f" Arguments to be passed to {args.Code:} : {extra_args:}")

    # Prepare job slicing
    dsets = subset_ids(args.Begin, args.Pop, args.Sets)
    print(" Slicing:", dsets)

    print()
    print(" Generated commands: ")

    outname = Path(args.Code).stem + "_"\
        + str(args.Begin).zfill(5) + "_"\
        + str(args.Pop).zfill(5) + "_"\
        + str(args.Sets)

    extname = ".ps1" if platform.system() == "Windows" else ".sh"
    if args.batch:
        fbatch = open("batch_"+outname+extname, "w")
    else:
        finter = open("interactive_"+outname+extname, "w")

    # Special action to transform the config filename into a resolved path name
    cfg = Configuration()
    conf_file = cfg.def_conf

    if "-c" in extra_args:
        idx = extra_args.index("-c")
        conf_file = str(Path(extra_args[idx+1]).resolve())
        extra_args[idx+1] = conf_file
    elif "--config" in extra_args:
        idx = extra_args.index("--config")
        conf_file = str(Path(extra_args[idx+1]).resolve())
        extra_args[idx+1] = conf_file

    if args.json is not None:  # get some information from the config file
        cfg.read_from_yaml(filename=Path(conf_file))

    # Loop over sets
    for [id1, id2] in dsets:

        # Pass extra arguments to external code
        if args.json is None:
            sys.argv = [args.Code] + extra_args + \
                     ["-f", str(id1), "-N", str(id2 - id1 + 1)]
        else:
            # jname = str(Path(args.json+"_"+str(id1)+"_"+str(id2)).resolve())
            jname = (args.json + "_"
                     + str(id1).zfill(cfg.dgt) + "_"
                     + str(id2).zfill(cfg.dgt) + ".json")
            jname = "'"+str(Path(jname).resolve())+"'"
            sys.argv = [args.Code] + extra_args + ["-f", jname]

        if args.Code.find("skygen") != -1:
            sky = Skies.command_line()
            job_cmd = sky.cmd_line

        elif args.Code.find("SoHAPPy") != -1:
            # If a json file prefix and folder are given, replace 'id1'
            cf = Configuration.command_line()
            job_cmd = cf.cmd_line

        else:
            print("generator.py : ", args.Code, "Not handled")

        # Create batch command if requested - display and write commands
        if args.batch:
            batch_cmd = slurm_cmd(job_cmd, name="SoHAPPy")
            print(" > ", batch_cmd)
            print(batch_cmd, file=fbatch)
        else:
            print(" >", job_cmd)
            print(job_cmd, file=finter)

    # Close the script file
    if args.batch:
        fbatch.close()
    else:
        finter.close()
