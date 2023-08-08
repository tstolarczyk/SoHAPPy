===========
Productions
===========

A production starts from a set of files with a collection of energy spectra
for a series of time slices, plus some physical information from the sources
(redshift, trigger time, sky position). The data are usually stored as `fits`
file with a specific format but other formats are possible (like `yaml` files
to reproduce simplistic lighcurve with two time and energy indices).

This section describes how to generate sky positions and/or trigger times
in files independent of the source file information in a way which is
manageable by `SoHAPPy`.

Generating visibilities
=======================

Introduction
------------

The module :obj:`skygen.py` generates trigger times, sky positions and the
related visibility from command line arguments. It is able to generate several
output files for a sequence of subsets of a given population. This subsets are
useful to launch `SoHAPPy` on a batch system
(see `batch submission`_).
The command line arguments are obtained from:


:code:`python skygen.py -h`

and the output is:

.. code-block::

    usage: skygen.py [-h] [-y YEAR1] [-n NYEARS] [-f FIRST] [-N NSRC] [-v VERSION]
                     [-D DAYS] [-V VISIBILITY] [-c CONFIG] [-o OUTPUT] [-s SEED]
                     [--debug] [--nodebug] [--trigger] [--notrigger] [--position]
                     [--noposition]

    Generate visibilities for SoHAPPy

    optional arguments:
      -h, --help            show this help message and exit
      -y YEAR1, --year1 YEAR1
                            First year
      -n NYEARS, --nyears NYEARS
                            Number of years
      -f FIRST, --first FIRST
                            First source identifier
      -N NSRC, --Nsrc NSRC  Number of sources
      -v VERSION, --version VERSION
                            version number
      -D DAYS, --days DAYS  Visibility range in days
      -V VISIBILITY, --visibility VISIBILITY
                            Visibility keywords
      -c CONFIG, --config CONFIG
                            Configuration file name
      -o OUTPUT, --output OUTPUT
                            Output base folder (path)
      -s SEED, --seed SEED  Seed from random generator
      --debug               Display processing details
      --nodebug             Does not display processing details
      --trigger             (re)generate dates
      --notrigger           Do not generate dates if already existing
      --position            Generate new positions in the sky
      --noposition          Do not generate new sky positions if already existing

    ---



Generating dates, positions and visibilities ex-nihilio
-------------------------------------------------------
Let's create dates, positions and visibilities for ten sources with date
randomly chosen between 2004 (January 1 :sup:`st`, 0h00'00) and end of 2013
(December 31 :sup:`st`, just before midnight,i.e. `23h59'59`):

``python skygen.py  -y 2004 -n 10 -f 1 -N 7``

The resulting file is in the folder:
``skygen_vis/`` **strictmoonveto_2004_10_1** ``/visibility``

It is intended to be used with a population having trigger dates along ten
years, starting in 2004. The last number indicates a version number (in case
other dates and positions would be randomly generated with the same
constrainsts). It can be changed with the `-v` option.

The file in the folder have the source indices in its name, from 1 to 7 :

 * ``DP_strictmoonveto_2004_10_1_1_7.yaml`` for the dates and positions.
 * ``strictmoonveto_2004_10_1_1_7.json`` for the visibilities.

Naturally, the next files in the series can be produced, e.g.
``python skygen.py  -y 2004 -N 10 -f 8 -N 5``

The code output is the following:

    .. code-block::

        --------------------------------------------------
         Generation number (version):  1
         Source identifiers from  8  to  12  ( 5 )
         Generated dates
          - from               :  2004
          - to                 :  2013
          - Duration (yr)      :  10
         Visibility:
          - Visibility keyword :  strictmoonveto
          -            range   :  3.0 d
          - Output folder      :  skygen_vis
         Debugging :  False

        Command line:
        skygen.py --year1 2004 --nyears 10 --first 8 --Nsrc 5 --version 1 --days 3.0 d --visibility strictmoonveto --output skygen_vis --seed 2022 --nodebug --notrigger --noposition
        --------------------------------------------------
        log information to file skygen.log
        +------------------------------------------------------------------------------+
        +                          Dates and positon - random                          +
        +------------------------------------------------------------------------------+

        +------------------------------------------------------------------------------+
        +                             Create output folder                             +
        +------------------------------------------------------------------------------+

        skygen_vis\strictmoonveto_2004_10_1\visibility Already exists
        +------------------------------------------------------------------------------+
        +                     Dumping generated dates and posiiton                     +
        +------------------------------------------------------------------------------+

         Output: skygen_vis\strictmoonveto_2004_10_1\visibility\DP_strictmoonveto_2004_10_1_8_12.yaml
        Done!
        +------------------------------------------------------------------------------+
        +                            Creating visibilities                             +
        +------------------------------------------------------------------------------+

        # 8  # 9  # 10  # 11  # 12   - Done
        +------------------------------------------------------------------------------+
        +                        Dumping generated visibilities                        +
        +------------------------------------------------------------------------------+

        Output: skygen_vis\strictmoonveto_2004_10_1\visibility\strictmoonveto_2004_10_1_8_12.json

        -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
         Duration      =    13.90 (s)
          per source   =     2.78 (s)
         ******* End of job - Total time =     0.23 min *****

        2023-08-08 16:27:04.568052

And, as an example, the date-position file,
`DP_strictmoonveto_2004_10_1_8_12.yaml`, in the same folder as the previous one,
has the following content:

.. code-block::

    created: 2023-08-08 16:26:50.677299
    id1: 8
    nsrc: 5
    seed: 2022
    start: 2004
    stop: 2013
    config: None
    basedir: skygen_vis
    key: strictmoonveto
    duration: 3.0 d
    version: 1
    ev8:     56041.9340364733             3.369101            -1.491226 # 2012-04-24T22:25:00.751
    ev9:     56028.1061458277           179.660812            52.684967 # 2012-04-11T02:32:51.000
    ev10:     56050.0662109165            40.818128            17.151789 # 2012-05-03T01:35:20.623
    ev11:     56501.0829470993            17.990647            52.553954 # 2013-07-28T01:59:26.629
    ev12:     54349.4663241195           246.746734            26.248797 # 2007-09-06T11:11:30.404

To get the files with other visibility conditions, pass the correct keyword to the command line :

  ``python skygen.py  -y 2004 -n 10 -f 1 -N 7 -V nomoonveto``

This will create the two files :
``skygen_vis\nomoonveto_2004_10_1\visibility\DP_nomoonveto_2004_10_1_8_12.yaml``

``skygen_vis\nomoonveto_2004_10_1\visibility\nomoonveto_2004_10_1_8_12.json``


Files with dates and positions already known
--------------------------------------------

In some cases, the source files have already trigger times and postions given
(and even sometimes a default visibility encoded). This is the case of the
first 1000 long afterglows primarliy studied with `SoHAPPy`.

The access to this information is done through the `config.yaml file` where
the file position will be read. Here is an example:

``python skygen.py  -y 2000 -n 5  -c config.yaml``

In this strict case, the dates and positions will not be recomputed since they
alreayd exist, despite a start year and a number of years are given.
`skygen` will generate the visibility for the default `strictmoonveto`
conditions.  Moreover, since no source range was  given only the first source is processed.
It correspond to the explicit command line:

``skygen.py --year1 2000 --nyears 5 --first 1 --Nsrc 1 --version 1 --days 3.0 d --visibility strictmoonveto --config config.yaml --output skygen_vis --seed 2022 --nodebug --notrigger --noposition``

The DP file contain the original information contained in the source file.
The visibility file contains the computed visibility.

More interesting is the case where it is requested to generate the trigger dates
(and not the source posiition) for a series of sources:

 ``python skygen.py  -y 2000 -n 5  -c config.yaml --trigger``

For illsutration, the start of the returned message is the following:

.. code-block::

    --------------------------------------------------
     Generation number (version):  1
     Source identifiers from  1  to  7  ( 7 )
     Generated dates
      - from               :  2000
      - to                 :  2004
      - Duration (yr)      :  5
     Reading dates and positions from surce files
      - Configuration file  :  config.yaml
      - New dates           :  True
      - New positions       :  False
     Visibility:
      - Visibility keyword :  strictmoonveto
      -            range   :  3.0 d
      - Output folder      :  skygen_vis
     Debugging :  False

    Command line:
    skygen.py --year1 2000 --nyears 5 --first 1 --Nsrc 7 --version 1 --days 3.0 d --visibility strictmoonveto --config config.yaml --output skygen_vis --seed 2022 --nodebug --trigger --noposition
    --------------------------------------------------

In final, the following files are
produced:
``strictmoonveto_2000_5_1/visibility/DP_strictmoonveto_2000_5_1_1_7.yaml``
``strictmoonveto_2000_5_1/visibility/strictmoonveto_2000_5_1_1_7.json``

Note, the following command where both the dates and positions are recomputed:

``python skygen.py  -y 2000 -n 5  -c config.yaml --trigger --position``

has the same effect as omiiting the configuration file with the drawback that
the source file are opened!

The `skygen` module
-------------------

.. autoclass:: skygen.Skies
    :special-members: __init__
    :members:


.. _batch submission:

Preparing batch submissions
===========================

The processing of 1000 GRB files with one iteration per trial and a classical
logarithmic spacing of time slices (ca. 40 slices in total) takes roughly 2
hours on a I5 -16GB laptop, including the visibility computation. With 100
iterations the computing time is only slightly increased.


The `generator`  module
-----------------------

.. automodapi:: generator
   :no-inheritance-diagram:
   :include-all-objects:
