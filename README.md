[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/dputhier/gtftk) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dputhier/gtftk/issues) [![HitCount](http://hits.dwyl.io/puthier/gtftk.svg)](http://hits.dwyl.io/puthier/gtftk)

# Python GTF toolkit (pygtftk)


The **Python GTF toolkit (pygtftk) package** is intented to ease handling of GTF (Gene Transfer Format) files. The pygtftk package is compatible with Python 2.7 and relies on **libgtftk**, a library of functions **written in C**. 

The package comes with a set of **UNIX commands** that can be accessed through the **gtftk  program**. The gtftk program proposes several atomic tools to filter, convert, or extract data from GTF files. The gtftk set of Unix command can be easily extended using a basic plugin architecture. All these aspects are covered in the help section.

While the gtftk Unix program comes with hundreds of unitary and functionnal tests, it is still upon  active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the github repository. 

## System requirements

Depending on the size of the GTF file, pygtftk may require lot of memory to perform selected tasks. A computer with 16Go is recommended in order to be able to pipe several commands when working with human annotations from ensembl release (e.g. 91).

At the moment, the gtftk program has been tested on:

- Linux (Ubuntu 12.04)
- OSX (Yosemite, El Capitan).

Windows users should try *Cygwin* or *Bash on Ubuntu on Windows* although it has not been tested in our hands.


## Installation through conda package building

Installation through **conda** should be the prefered install solution. Although the GTF interface of pygtftk should work properly after a pip install (see next section), the UNIX commands (gtftk program) require several external dependencies with some version constrains.

At the moment, there is no built conda package available. You can however create an environment with all prerequisites using the commands below.
If conda is not available on your system, first install miniconda from the official [web site](http://conda.pydata.org/miniconda.html).

    git clone git@github.com:dputhier/pygtftk.git pygtftk
    cd pygtftk
    conda env create -n pygtftk -f conda/env.yaml
    source activate pygtftk
    make install

Note for developpers: You can install the develop branch using the same approach.

    git clone git@github.com:dputhier/pygtftk.git pygtftk
    cd pygtftk
    git checkout develop
    conda env create -n pygtftk_dev -f conda/env.yaml
    source activate pygtftk_dev
    make install
    
## Installation through pip 

### Prerequesites
 
Again, this is not the prefered way for installation. Please choose conda whenever possible. The gtftk program requires at least bedtools (command line), R and some R packages. To get the list of required R packages, use:

    gtftk -r 

To get the list of all dependencies and associated versions, please have a look at conda/env.yaml file.

### Running pip 

Installation through pip can be done as follow.

    pip install -r requirements.txt
    pip install pygtftk
    

## Building doc files

The following commands can be used to build the help files:

    cd docs/manual/
    make html
    ls build/html/index.html
    
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

    make test

Or run tests in parallel using:

    make test_para -j 10 # Using 10 cores
        
## Running unitary tests

Several unitary tests have been implemented using doctests. You can run them using nose through the following command line:

    make nose
