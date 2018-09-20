[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/dputhier/pygtftk) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dputhier/pygtftk/issues) [![HitCount](http://hits.dwyl.com/dputhier/pygtftk.svg)](http://hits.dwyl.com/dputhier/pygtftk)


# Python GTF toolkit (pygtftk)


The **Python GTF toolkit (pygtftk) package** is intented to ease handling of GTF (Gene Transfer Format) files. The pygtftk package is compatible with Python 2.7 and Python >=3.5 and relies on **libgtftk**, a library of functions **written in C**. 

The package comes with a set of **UNIX commands** that can be accessed through the **gtftk  program**. The gtftk program proposes several atomic tools to filter, convert, or extract data from GTF files. The gtftk set of Unix commands can be easily extended using a basic plugin architecture. All these aspects are covered in the help sections.

While the gtftk Unix program comes with hundreds of unitary and functional tests, it is still upon  active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the github repository. 

## System requirements

Depending on the **size of the GTF file**, pygtftk and gtftk may require lot of memory to perform selected tasks. A computer with 16Go is recommended in order to be able to pipe several commands when working with human annotations from ensembl release (e.g. 91).

At the moment, the gtftk program has been tested on:

- Linux (Ubuntu 12.04 and 18.04)
- OSX (Yosemite, El Capitan).


## Installation through conda package building

Installation through **conda** should be the **prefered install solution**. The pygtftk package and gtftk command line tool require external dependencies with some version constrains (e.g. bedtools that, we observed, displays some back compatibility issues).

A conda package will be available in the near future. In the meantime, you can however create an environment with all prerequisites using the commands below.
If conda is not available on your system, first install miniconda from the official [web site](http://conda.pydata.org/miniconda.html).

    git clone git@github.com:dputhier/pygtftk.git pygtftk
    cd pygtftk
    conda env create -n pygtftk_py3k -f conda/env.yaml python=3.6
    source activate pygtftk_py3k
    make install
    # It is important to call gtftk -h
    # to find and dump plugin parsers
    # before going further
    gtftk -h 

## Installation through setup.py

This is not the prefered way for installation. Choose conda whenever possible. The gtftk Unix command line program has been tested with bedtools 2.27.1 (be aware that we have encountered some back compatibility issues with bedtools).

    git clone git@github.com:dputhier/pygtftk.git pygtftk
    cd pygtftk
    # Check your Python version before (2.7 or >=3.5)
    pip install -r requirements.txt
    python setup.py install

## Installation through pip (not functional at the moment)

### Prerequesites
 
Again, this is not the prefered way for installation. Please choose conda whenever possible. The gtftk Unix command line program has been tested with bedtools 2.27.1 (be aware that we have encountered some back compatibility issues with bedtools).

### Running pip 

Installation through pip can be done as follow.

    pip install -r requirements.txt
    pip install pygtftk
    # It is important to call gtftk -h
    # to find and dump plugin parsers
    # before going further
    gtftk -h     

## Running functional tests

A lot of functional tests have been developed to ensure consistency with expected results. This does not rule out that bugs may hide throughout the code... In order to check that installation is functional you may be interested in running functional tests. The definition of all functional tests declared in  gtftk commands is accessible using the -p/--plugin-tests argument:

    gtftk -p

To run the tests, you will need to install [bats (Bash Automated Testing System)](https://github.com/sstephenson/bats). Once bats is installed run the following commands:

    # The tests should be run in the pygtftk git
    # directory because several tests contains references (relative path)
    # to file enclosed in pygtftk/data directory.
    gtftk -p > gtftk_test.bats
    bats gtftk_test.bats

Note, alternatively you may directly call the tests using the Makefile
    
    make clean
    make test

Or run tests in parallel using:

    make clean
    make test_para -j 10 # Using 10 cores
        
## Running unitary tests

Several unitary tests have been implemented using doctests. You can run them using nose through the following command line:

    make nose
