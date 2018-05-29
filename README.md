[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://github.com/dputhier/gtftk) [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dputhier/gtftk/issues) [![HitCount](http://hits.dwyl.io/puthier/gtftk.svg)](http://hits.dwyl.io/puthier/gtftk)

# GTFtoolkit


The GTFtoolkit (gtftk) program is intented to ease handling of GTF (Gene Transfer Format) files. One of our objective is to propose a set of tools that perform common and simple tasks to limit the use of perl/awk onliners that tend to limit readability of bioinformatics workflows. The gtftk program implements several atomic tools than can be used and piped to filter, convert, or extract data from GTF files. It has more recently evolved to a more full-feature program that can be used to accomplish more sophisticated tasks ranging from annotation, functional enrichment or coverage analysis (*e.g* metaprofile around genomic features). The gtftk framework can be used directly through dedicated Unix commands or through a Python API that makes call to a dynamic library written in C. The gtftk set of Unix command, can be easily extended using a simple plugin architecture. All these aspects are covered in the help section.
While the gtftk program comes with hundreds of unitary and functionnal tests, it is still upon  active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the github repository. 

## System requirements

Depending on the GTF file, gtftk may require lot of memory to perform selected task. A computer with 16Go is recommended in order to be able to pipe several commands when working with human annotations from ensembl release (e.g. 91).

The gtftk program has been tested on:

- Linux (Ubuntu 12.04)
- OSX (Yosemite, El Capitan).

Windows user should try *Cygwin* or *Bash on Ubuntu on Windows* although it has not been tested in our hand.


## Installation through conda package building

As gtftk requires several external programs and libraries, installation through **conda should really be the preferred solution** in order to limit dependency issues.

At the moment, there is no built conda package available. You can however create an environment with all prerequisites using the commands below.
If conda is not available on your system, first install miniconda from the official [web site](http://conda.pydata.org/miniconda.html).

    git clone git@github.com:dputhier/gtftk.git gtftk
    cd gtftk
    conda env create -n gtftk -f conda/env.yaml
    source activate gtftk
    make install

Note for developpers: You can install the develop branch using the same approach.

    git clone git@github.com:dputhier/gtftk.git gtftk
    cd gtftk
    git checkout develop
    conda env create -n gtftk_dev -f conda/env.yaml
    source activate gtftk_dev
    make install
    
## Installation through setup.py

### Prerequesites
 
The gtftk program requires at least bedtools (command line), R and some R packages. To get the list of required R packages, use:


    gtftk -r 

To get the list of all dependencies and associated versions, have a look at conda/env.yaml file.

### Installation with setup.py

As any Python package, installation is normaly performed using:

    git clone git@github.com:dputhier/gtftk.git gtftk
    cd gtftk
    pip install -r requirements.txt
    python setup.py install
    # or for those without admin rights
    # python setup.py install --user

## Building doc files

The following commands can be used to build the help files:

    cd docs/manual/
    make html
    ls build/html/index.html
    
## Running functional tests

There a lot of functional tests that have been developped through gtftk project to ensure consistancy with expected results. This does not rule out that bugs may hide throughout the code... In order to check that installation is functional you may be interested in running functional tests. The definition of all functional tests declared in  gtftk commands/plugins is accessible using the -p/--plugin-tests argument:

    gtftk -p

To run these tests You will need to install [bats (Bash Automated Testing System)](https://github.com/sstephenson/bats). Once bats is installed you may run the tests using the following commands:

    # The tests should be run in the gtftk git
    # directory because several tests contains references (relative path)
    # to file enclosed in gtftk/data directory.
    gtftk -p > gtftk_test.bats
    bats gtftk_test.bats

Note, alternatively you may directly call the tests through the Makefile

    # make test
        
## Running unitary tests

Several unitary tests have been implemented using doctests. You can run them using nose through the following command line:

    make nose
