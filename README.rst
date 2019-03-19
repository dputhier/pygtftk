.. image:: https://img.shields.io/github/license/mashape/apistatus.svg
    :alt: Licence
    :target: https://github.com/dputhier/pygtftk

.. image:: https://badge.fury.io/py/pygtftk.svg
    :alt: PyPI
    :target: https://badge.fury.io/py/pygtftk

.. image::  https://img.shields.io/badge/contributions-welcome-brightgreen.svg
    :alt: GitHub
    :target: https://github.com/dputhier/pygtftk/blob/master/CONTRIBUTING.rst

.. image:: https://readthedocs.org/projects/pygtftk/badge/?version=master
    :alt: Documentation Status
    :target: https://pygtftk.readthedocs.io/en/latest/

.. image:: https://travis-ci.org/dputhier/pygtftk.svg?branch=master
    :target: https://travis-ci.org/dputhier/pygtftk

.. image:: https://img.shields.io/github/repo-size/badges/shields.svg
    :target: https://travis-ci.org/dputhier/pygtftk

.. image:: https://anaconda.org/guillaumecharbonnier/pygtftk/badges/installer/conda.svg
    :target: https://anaconda.org/guillaumecharbonnier/pygtftk

.. image:: https://anaconda.org/guillaumecharbonnier/pygtftk/badges/platforms.svg
    :target: https://anaconda.org/guillaumecharbonnier/pygtftk

.. image:: https://anaconda.org/guillaumecharbonnier/pygtftk/badges/latest_release_date.svg
    :target: https://anaconda.org/guillaumecharbonnier/pygtftk

.. image:: https://anaconda.org/guillaumecharbonnier/pygtftk/badges/downloads.svg
    :target: https://anaconda.org/guillaumecharbonnier/pygtftk

.. highlight-language: shell



Python GTF toolkit (pygtftk)
=============================

The **Python GTF toolkit (pygtftk) package** is intented to ease handling of GTF/GFF2.0 files (Gene Transfer Format). It currently does not support GFF3 file format. The pygtftk package is compatible with Python  >=3.5,<3.7 and relies on **libgtftk**, a library of functions **written in C**.

The package comes with a set of **UNIX commands** that can be accessed through the **gtftk  program**. The gtftk program proposes several atomic tools to filter, convert, or extract data from GTF files. The gtftk set of Unix commands can be easily extended using a basic plugin architecture. All these aspects are covered in the help sections.

While the gtftk Unix program comes with hundreds of unitary and functional tests, it is still upon  active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the github repository. 

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

Installation through **conda** should be the **preferred install solution**. The pygtftk package and gtftk command line tool require external dependencies with some version constrains.

If conda is not available on your system, first install miniconda from the official `web site <http://conda.pydata.org/miniconda.html>`_ and make sure you have bioconda and conda-forge channels set up in the order below. ::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Then you can simply install pygtftk in its own isolated environment and activate it. ::

    conda create -n pygtftk -c guillaumecharbonnier pygtftk
    conda activate pygtftk


Installation through setup.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is not the preferred way for installation. Choose conda whenever possible. We have observed several issues with dependencies that still need to be fixed. ::

    git clone git@github.com:dputhier/pygtftk.git pygtftk
    cd pygtftk
    # Check your Python version (>=3.5,<3.7)
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



Documentation
--------------

Documentation about the latest release is dynamically produced and available at `readthedoc server <https://pygtftk.readthedocs.io/en/latest/>`_.

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


