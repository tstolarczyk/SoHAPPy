Introduction
############

This paragraph gives the first indications to install and run `SoHAPPy`. See
eselwhere how to `run large batch productions <production.html>`_ and perform
a `spectral analysis <spectral_analysis.html>`_.


Installation, dependencies
==========================

It is suggested to install the code inside a dedicated environment based on
`gammapy`. This version is compatible with the following dependencies:

    * `gammapy <https://gammapy.org/>`_: 0.18.2
	   Use :code:`conda env create -f '..\gammapy-0.18.2-environment.yml'` to
	   create a dedicated environement (the `yaml` file is found from the
	   `gammapy` site).
	   Activate the environment, :code:`conda activate gpy0.18.2`, pursue with
	   further installations.
    * `astroplan <https://pypi.org/project/astroplan/>`_ : version 0.8
      (use :code:`pip install astroplan`).
    * `seaborn <https://seaborn.pydata.org/>`_ : version 0.12.2

At this stage, check that the line command :code:`python SoHAPPy.py -h` gives
the list of possible options. Note that this version is compatible with
`spyder 5.2.2` with `spyder-kernels 2.2.1`. During the code analysis,
`pylint` is used, and a `.pylintrc` file is provided.

In order to get results two inputs are required:
    * **Instrument response functions** (IRFs) for CTA (or alternative
      instruments if they are implemented). For CTA, the official response
      functions for the alpha configuration ar stored on
      `Zenodo for prod5 <https://zenodo.org/record/5499840#.YUya5WYzbUI>`_.
      In `SoHAPPy` the IRF files are in folders with a fixed specific
      organisation (see `IRF <irf.rst>`_)
    * **Input object files** with a `dedicated format <file_format.rst>`_.

The code has been checked to work with `matplotlib 3.4.1`
(`3.2.1` requested by `gammapy 0.18.2`), a version that helps solving inconsistencies
in the installation.
Note that the `analysis/population/universe.py` module requires `matplotlib 3.5`
or higher to be used.

Repository structure
====================
The `SoHAPPy` repository has 3 subfolders:
    * ``analysis`` : contains the notebooks and scripts to analyse the
      `SoHAPPy` outputs, with the following subfolders:

      * ``population`` to analyze the population files and get the detection
        rates;
      * ``spectra`` to perform a spctral analysis of one Monte Carlo
        realization on one of the source in the population;
      * ``prompt`` specific preliminary tools to analyze the prompt component.

    * ``docs``: contains the `Sphinx` information to generate the present
      documentation (Version 4.5.0 was used for this document).
    * ``data``, contains the following subfolders:

       * ``historical``: a list of `yaml` files containing parameters to
         generate GRB light curves from an analytical function;
       * ``swift``: contains information on the Swift satellite latency time
         and in particular a distribution to generate realistic delays;
       * ``ebl``: this is a copy of the `ebl` subfolder of the
         `gammapy extras <https://github.com/gammapy/gammapy-extra>`_ and is
         necessary to have the `gammapy` models used in `SoHAPPy`. It also
         contains some additionnal EBL models for tests.

Input and output
================
`SoHAPPy` requires to have an input base folder and an output base folder
defined through environment variables:

* :code:`HAPPY_IN` defines the input base folder
* :code:`HAPPY_OUT` defines the output base folder

The variables are set as follows :

* `Windows powershell` : :code:`$Env:HAPPY_IN = "D:\\\\CTA\\SoHAPPy\\input"`
* `Interactive python` : :code:`os.environ["HAPPY_IN"] = "D:\\\\CTA\\SoHAPPy\\input"`
* `bash shell` : :code:`$ export HAPPY_IN="/home/user/SoHAPPY/input"`

Input folder
------------

The input folder has to contain the following
subfolders:

    * a folder handling the light curves and some subfolders with various
      subpopulation, e.g. `lightcurves/pop_name`. The name of the population
      folder `pop_name` is used in the output folder hierarchy (the last name
      in the path is used, see below).
    * a folder with the IRFs. The IRF file repository is expected to be
      organised in a strict way as the folder structure is used internally
      by `SoHAPPy` to get the suitable IRFs file (see
      `Instrument response functions <irf.html>`_ for details)

The graph below illustrates this organisation. For the light curves, none of
the names nor the organisation are mandatory (They will be passed by the user
in the configuration file or on the command line).

.. code-block:: python

 |  +input
 |    +-- lightcurves
 |           +-- long_grb
 |           +-- special_grbs
 |           +-- other_grbs
 |    +-- irf

Output folder
-------------

A `SoHAPPy` run output consists in 3 files:

* a log file, **analysis.log**, that keeps track of the processing;
* a result text file, **data.txt**, which contains the results of the
  simulation and analysis.
* a copy of the **configuration file** initially given (or the default) with
  some parameters changed through the command line (so that reusing this
  configuration file allows reproducing exactly the simulation/analysis).

These three files are stored into a `tar.gz` archive. If the code is run more
than once, then  the archive is copied and its name contains an information
on the date at which it was duplicated (results are never automatically
overwritten). The user can decide if both the files and the archive are kept,
or only the archive. The output files(s) are stored in a subfolder created
on-the-fly by `SoHAPPy`. The output folder position is free but its internal
structure is generated by SoHAPPy.


The output folder file structure follows this organisation:

.. code-block:: python

 | + output
 |         +-- pop_name_1
 |         |            +--  out_dir
 |         |                       +-- vis_name
 |         |                                  +--- vis_name_id1_id2
 |         |                                  +--- vis_name_id2_id3
 |         |                                  +--- vis_name_id3_id4
 |         |                                  +...
 |         +-- pop_name_2


`output` is the output base folder, `pop_name_1` is the population input
folder stored in the input base folder. The folder name,`out_dir`, is chosen by
the user and refer to his analysis (e.g. can be `test_omega` for results
testing the `omega` configuration). `vis_name` refers to the assumption
on the observation, including the minimal altitude for observation, the Moon
light veto etc. A collection of possible names is found in the `SoHAPPy`
`visibility.yaml` file (see `Visibility <visibility.html>`_) where more
visibility configurations can be added.

The last folder name use again the visibility keyword `vis_name` and add the
first and last source identifiers of the run. In case only one source is
analysed the names has only the first identifier (`vis_name_id0`).

Note that it is possible to run `SoHAPPy` on a limited number of sources and
obtain for each of them additionnal output files to be used for a spectral
analysis (save the simulation for this using the :code:`save_simu` keyword in
the `configuration file <configuration.html>`_ ).

Required data
-------------
The path to get access to these data are given in the configuration file or
on the command line. The necessary data files are the following:

    * **astrophysical object data files**, one per source, containing the
      energy spectra along time slices.
    * for each of these files, the **position in ra-dec** and the
      **explosion time** (referred often as the trigger time), generated
      independently from the :obj:`skygen <../../skygen.py>` application and
      stored in a `yaml` file (or a collection of `yaml` files). See the
      chapter on `productions <production.html>`_ for details). In some
      cases, this information can be inside the input astrophysical source
      data files.
    * a **visibility** file giving the rise and set time for the Sun, the Moon
      and the source itself. This information can be generated  from the
      :obj:`skygen` application. See the chapter on
      `productions <production.html>`_ for details). In some cases this
      visibility is encoded in the astrophysical source files or is computed
      on-the-fly from a given keyword referenced in the `visibility.yaml`
      file (see `visibilities <visibility.html>`_ for details).
    * The **instrument response function** set used (e.g. `prod3`, `omega`)
      and extra information on the array or subarry used, or specific flags
      used during the simulation and analysis (e.g. the slewing time)

Launching the code
==================

The steering parametres are obtained from the `yaml` configuration file.
A file is provided with the releas, `data/config_ref.yaml` and is used by
default.

    * :code:`python SoHAPPy.py` would simply run the code from the code folder
      with the `data/config_ref.yaml` parameters.
    * Some of the paramters in the configuration file can be superseded on
      the command line. :code:`python SoHAPPy.py -h` gives the list of
      accessible parameters. Here is the output:

.. code-block:: python

 | usage: SoHAPPy.py [-h] [-f FIRST] [-N NSRC] [-n NITER] [-o OUTPUT]
 |                   [-i INPUT] [-c CONFIG] [-V VISIBILITY] [-d DEBUG]
 |
 | SoHAPPy optional arguments:
 |  -h, --help            show this help message and exit
 |  -f FIRST, --first FIRST
 |                       First source id
 |  -N NSRC, --nsrc NSRC  Number of source files
 |  -n NITER, --niter NITER
 |                       Number of Monte Carlo iteration
 |  -o OUTPUT, --output OUTPUT
 |                       Output base folder (path)
 |  -i INPUT, --input INPUT
 |                       Input base folder (path)
 |  -c CONFIG, --config CONFIG
 |                        Configuration file name
 |  -m MAXNIGHT, --maxnight MAXNIGHT
 |                       Maximal number of nights
 |  -s SKIP, --skip SKIP
 |                       Number of nights to skip
 |  -V VISIBILITY, --visibility VISIBILITY
 |                       Visibility keyword
 |  -d DEBUG, --debug DEBUG
 |                       Debugging flag



In particular the default configuration file name can be superseded :
:code:`python SoHAPPy.py -c myconfig.yaml` and parameters in `myconfig.yaml`
can at their turn be superseded :code:`python SoHAPPy.py -c myconfig.yaml -N 1`

For large productions, it is useful to run `SoHAPPy` on subsets.
The parametesr accessible on the command line are intended for this purpose.
The base input and output base folders are considered installation dependent
whereas the subfolders can be changed to differentiate the runs.
See more on `productions and batch submissions <production.html>`_

