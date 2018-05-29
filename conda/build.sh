#!/bin/bash

git clone  git@github.com:dputhier/gtftk.git
cd gtftk
git checkout develop
make install

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
