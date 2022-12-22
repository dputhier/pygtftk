===================        =================
Pip package                |Pippackage|_
Bioconda package           |bioconda|_
License                    |license|_
Platforms                  |platform|_
Languages                  |lang|_
Build status               |build|_
Repository size            |size|_
Latest conda               |latestconda|
Downloads                  |downloads|_
Codacy                     |codacy|_
Contribution               |contrib|_
Issues                     |issues|_
Citing pygtftk             |citing_pygtftk|_
Citing OLOGRAM             |citing_ologram|_
Citing OLOGRAM-MODL        |citing_ologram_modl|_
Documentation              |documentation|_
===================        =================

.. |codacy| image:: https://app.codacy.com/project/badge/Grade/0e5b6b1af3a84af7a2d57e0465204e7f
.. _codacy: https://www.codacy.com/gh/dputhier/pygtftk/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=dputhier/pygtftk&amp;utm_campaign=Badge_Grade

.. |bioconda| image:: https://anaconda.org/bioconda/pygtftk/badges/version.svg
.. _bioconda: https://anaconda.org/bioconda/pygtftk

.. |license| image:: https://img.shields.io/github/license/dputhier/pygtftk.svg
.. _license: https://github.com/dputhier/pygtftk

.. |pippackage| image:: https://badge.fury.io/py/pygtftk.svg
.. _pippackage: https://badge.fury.io/py/pygtftk

.. |contrib| image::  https://img.shields.io/badge/contributions-welcome-brightgreen.svg
.. _contrib: https://github.com/dputhier/pygtftk/blob/master/CONTRIBUTING.rst

.. |build| image:: https://travis-ci.org/dputhier/pygtftk.svg?branch=master
.. _build: https://travis-ci.org/dputhier/pygtftk

.. |size| image:: https://img.shields.io/github/repo-size/badges/shields.svg
.. _size: https://travis-ci.org/dputhier/pygtftk

.. |platform| image:: https://anaconda.org/bioconda/pygtftk/badges/platforms.svg
.. _platform: https://anaconda.org/bioconda/pygtftk

.. |latestconda| image:: https://anaconda.org/bioconda/pygtftk/badges/latest_release_date.svg
.. _latestconda: https://anaconda.org/bioconda/pygtftk

.. |downloads| image:: https://anaconda.org/bioconda/pygtftk/badges/downloads.svg
.. _downloads: https://anaconda.org/bioconda/pygtftk

.. |reference| image:: https://img.shields.io/reference-yes-green.svg
.. _reference: http://hits.dwyl.io/dputhier/pygtftk

.. |issues| image:: https://img.shields.io/github/issues-raw/dputhier/pygtftk.svg
.. _issues: https://github.com/dputhier/pygtftk/issues

.. |citing_pygtftk| image:: https://img.shields.io/badge/citing__pygtftk-https%3A%2F%2Fdoi.org%2F10.1093%2Fbioinformatics%2Fbtz116-blue
.. _citing_pygtftk: https://doi.org/10.1093/bioinformatics/btz116

.. |documentation| image:: https://img.shields.io/badge/Documentation-https%3A%2F%2Fdputhier.github.io%2Fpygtftk%2F-blue.svg
.. _documentation: https://dputhier.github.io/pygtftk/

.. |lang| image:: https://img.shields.io/badge/Languages-Python%2C%20C%2C%20Cython%2C%20C++-blue.svg
.. _lang: https://github.com/dputhier/pygtftk

.. |citing_ologram| image:: https://img.shields.io/badge/citing__ologram-https%3A%2F%2Fdoi.org%2F10.1093%2Fbioinformatics%2Fbtz810-blue
.. _citing_ologram: https://doi.org/10.1093/bioinformatics/btz810

.. |citing_ologram_modl| image:: https://img.shields.io/badge/citing__ologram--modl-https%3A%2F%2Fdoi.org%2F10.1093%2Fnargab%2Flqab114--modl-blue
.. _citing_ologram_modl: https://doi.org/10.1093/nargab/lqab114


Python GTF toolkit (pygtftk)
=============================

The **Python GTF toolkit (pygtftk) package** is intended to ease handling of GTF/GFF2.0 files (Gene Transfer Format). It currently does not support GFF3 file format. The pygtftk package is compatible with Python 3.9 and relies on **libgtftk**, a library of functions **written in C**.

The package comes with a set of **UNIX commands** that can be accessed through the **gtftk  program**. The gtftk program proposes several atomic tools to **filter, convert, or extract data from GTF files**.

The newly released command, **OLOGRAM (OverLap Of Genomic Regions Analysis using Monte Carlo)** may be used to compute overlap statistics between user supplied regions (BED format) and annotation derived from :

- Gene centric features enclosed in a GTF (e.g. exons, promoters, terminators...).
- Regions in a GTFs flagged with built-in keys/values (e.g. check the 'gene_biotype' as provided by ensembl GTFs of the regions in which peaks fall).
- Same with custom keys/values through the gtftk CLI (e.g. adding a numeric value to a gene and discretizing this value to create gene classes).
- User supplied BED files.

With the most recent update, OLOGRAM is now also capable of computing the enrichment of **n-wise combinations** (ie. A+B, A+B+C, etc.) to find correlation groups of regions. Please see the documentation page of *OLOGRAM* for more details.

The gtftk set of Unix commands can be easily extended using a basic plugin architecture.

All these aspects are covered in the help sections ; please see the `documentation <https://dputhier.github.io/pygtftk/>`_.

While the gtftk Unix program comes with hundreds of unitary and functional tests, it is still in active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the GitHub repository.


Documentation
--------------

Documentation about the latest release is available as a `github page <https://dputhier.github.io/pygtftk/>`_.

Documentation about OLOGRAM (OverLap Of Genomic Regions Analysis using Monte Carlo) can be found in `the 'ologram' section of the documentation <https://dputhier.github.io/pygtftk/ologram.html>`_.

**NB:** The readthedoc version won't be maintained and will be closed in the near future. This choice was motivated by the impossibility to maintain a dynamic documentation (using sphinx/sphinxcontrib-programoutput) given the computing time provided by readthedoc server.

Note that example dataset are available to test the various subcommands (see `documentation page <https://dputhier.github.io/pygtftk/information.html#get-example>`_).

    gtftk get_example  -h
    # E.g. to get all file from the 'simple' dataset
    gtftk get_example -d simple -f "*"


System requirements
--------------------

Depending on the **size of the GTF file**, pygtftk and gtftk may require lot of memory to perform selected tasks. A computer with 16Go is recommended in order to be able to pipe several commands when working with human annotations from ensembl release (e.g. 91). When working with a cluster think about reserving sufficient memory.

At the moment, the gtftk program has been tested on:

  - Linux (Ubuntu 12.04 and 18.04)
  - OSX (Yosemite, El Capitan, Mojave).

Installation
-------------

Installation through conda package building
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installation through **conda** should be the **preferred install solution**. The pygtftk package and gtftk command line tool require external dependencies (bedtools ">v2.23.1", graphviz, unzip) with some version constrains.

If conda is not available on your system, first install miniconda from the official `web site <http://conda.pydata.org/miniconda.html>`_ and make sure you have bioconda and conda-forge channels set up in the order below. ::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Then you can simply install pygtftk in its own isolated environment and activate it. ::

    conda create -n pygtftk pygtftk
    conda activate pygtftk


Installation through setup.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is not the preferred way for installation. Choose conda whenever possible. We have observed several issues with dependencies that still need to be fixed. ::

    git clone http://git@github.com:dputhier/pygtftk.git pygtftk
    cd pygtftk
    # Check your Python version (>=3.8,<3.9)
    pip install -r requirements.txt
    python setup.py install


Installation through pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Prerequisites**


Again, this is not the preferred way for installation. Please choose conda whenever possible. We have observed several issues with dependencies that still need to be fixed.

**Running pip**


Installation through pip can be done as follow. ::

    pip install -r requirements.txt
    pip install pygtftk
    # It is important to call gtftk -h
    # to look for plugins and their
    # CLI in ~/.gtftk
    # before going further
    gtftk -h



Testing
--------

Running functional tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A lot of functional tests have been developed to ensure consistency with expected results. This does not rule out that bugs may hide throughout the code... In order to check that installation is functional you may be interested in running functional tests. The definition of all functional tests declared in  gtftk commands is accessible using the -p/--plugin-tests argument: ::

    gtftk -p


To run the tests, you will need to install `bats (Bash Automated Testing System) <https://github.com/sstephenson/bats>`_. Once bats is installed run the following commands: ::

    # The tests should be run in the pygtftk git
    # directory because several tests contains references (relative path)
    # to file enclosed in pygtftk/data directory.
    gtftk -p > gtftk_test.bats
    bats gtftk_test.bats


Note, alternatively you may directly call the tests using the Makefile. ::

    make clean
    make test


Or run tests in parallel using: ::

    make clean
    make test_para -j 10 # Using 10 cores



Running unitary tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several unitary tests have been implemented using doctests. You can run them using nose through the following command line: ::

    make nose
