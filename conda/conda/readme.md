# Installation of the development environment
With a working conda installation (gcc required), just run:
```
conda env create -f env.yaml
source activate gtftk_dev
cd ..
make install
```

# Building gtftk package
The current build for gtftk should work with the following channels:
```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels prometeia
```
Then use conda-build in the conda directory containing 'build.sh' and 'meta.yaml'.

