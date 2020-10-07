"""
The pygtfk package.

The Python GTF toolkit (pygtftk) package is intented to ease handling of GTF (Gene Transfer Format) files. The pygtftk package relies on libgtftk, a library of functions written in C.
The package comes with a set of UNIX commands that can be accessed through the gtftk program. The gtftk program proposes several atomic tools to filter, convert, or extract data from GTF files. The gtftk set of Unix commands can be easily extended using a basic plugin architecture. All these aspects are covered in the help section.

Authors: D. Puthier and F. Lopez
"""

# -------------------------------------------------------------------------
# A set of builtin packages
# -------------------------------------------------------------------------

import re
import shutil
import subprocess
import sys
from subprocess import DEVNULL

import glob
import hashlib
import numpy as np
import os
import platform
from Cython.Distutils import build_ext
from tempfile import NamedTemporaryFile

# -------------------------------------------------------------------------
# Python compiler version
# -------------------------------------------------------------------------
print("Python version and gcc version used for compilation : ")
print(sys.version)

# -------------------------------------------------------------------------
# Python Version
# -------------------------------------------------------------------------

PY2 = sys.version_info[0] == 2

if PY2:
    sys.stderr.write("\nStarting from version 0.9.8, gtftk does not support python 2 anymore.\n")
    sys.exit(1)

# -------------------------------------------------------------------------
# Informations about the project
# -------------------------------------------------------------------------


__author__ = 'Fabrice Lopez and Denis Puthier'
__email__ = 'denis.puthier@univ-amu.fr'
__description__ = 'The Python GTF toolkit (pygtftk) package: easy handling of GTF files'
__license__ = 'GPL-2'
__url__ = 'https://github.com/dputhier/pygtftk'
__url_source__ = 'https://github.com/dputhier/pygtftk'
__url_tracker__ = 'https://github.com/dputhier/pygtftk'
__keywords__ = 'genomics bioinformatics GTF BED'
__python_requires__ = '>=3.6,<=3.8'
__classifiers__ = ("License :: OSI Approved :: MIT License",
                   "Operating System :: MacOS",
                   "Operating System :: POSIX :: Linux",
                   "Development Status :: 4 - Beta",
                   "Environment :: Console",
                   "Programming Language :: Python :: 3.6",
                   "Programming Language :: Python :: 3.7",
                   "Programming Language :: Python :: 3.8",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Topic :: Scientific/Engineering :: Bio-Informatics",
                   "Operating System :: POSIX :: Linux",
                   "Operating System :: MacOS",
                   "Topic :: Documentation :: Sphinx"
                   )

# -------------------------------------------------------------------------
# Printing Python version
# -------------------------------------------------------------------------

sys.stderr.write('Python version : ' + str(sys.version_info) + '\n')
sys.stderr.write('Python path : ' + str(sys.prefix) + '\n')

# -------------------------------------------------------------------------
# Check setuptools is installed
# -------------------------------------------------------------------------

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    sys.stderr.write("Please install setuptools before installing pygtftk.\n")
    exit(1)

# -------------------------------------------------------------------------
# Check gtftk version
# -------------------------------------------------------------------------

version_fh = open("pygtftk/version.py")

for i in version_fh:
    if "__base_version__" in i:
        base_version = i.split("=")[1]
        base_version = re.sub("['\" \n\r]", "", base_version)

sha = ""

try:
    import git
    from git import InvalidGitRepositoryError

    try:
        repo = git.Repo(search_parent_directories=True)
        sha = repo.head.object.hexsha
        sha = repo.git.rev_parse(sha, short=4)

        if sha != "" and not os.path.exists("release_in_progress"):
            __version__ = base_version + ".dev0+" + sha
        else:
            __version__ = base_version
    except InvalidGitRepositoryError:
        __version__ = base_version

except ImportError:
    __version__ = base_version

version_file = open('pygtftk/version.py', "w")
version_file.write("__base_version__='" + base_version + "'\n")
version_file.write("__version__='" + __version__ + "'\n")
version_file.close()

# -------------------------------------------------------------------------
# Building C library
# -------------------------------------------------------------------------

cmd_src_list = glob.glob("pygtftk/src/libgtftk/*.c")
cmd_src_list += glob.glob("pygtftk/src/libgtftk/command/*.c")
cmd_src_list = list(set(cmd_src_list))

'''
if platform == "darwin":
    vars = sysconfig.get_config_vars()
    # vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')
    # dyn_lib_compil = ['-dynamiclib', '-shared']
    dyn_lib_compil = []
else:
    dyn_lib_compil = ['-shared']
'''

extra_compile_args = ['-Ipygtftk/src/libgtftk',
                      '-O3',
                      '-Wall',
                      '-fPIC',
                      '-MMD',
                      '-MP',
                      '-fmessage-length=0']  # + dyn_lib_compil

lib_pygtftk = Extension(name='pygtftk/lib/libgtftk',
                        include_dirs=[
                            'pygtftk/src/libgtftk'],
                        library_dirs=['/usr/lib'],
                        libraries=['z'],
                        extra_compile_args=extra_compile_args,
                        sources=cmd_src_list)

# ----------------------------------------------------------------------
#  Building Cython modules - mostly for OLOGRAM
# ----------------------------------------------------------------------

platform.system()

extra_comp_cython = ['-W', '-O3']
extra_link_cython = []

# Use OpenMP only on Linux, as clang by default does not support it on OSX
# TODO Make it a parameter
if platform.system() == 'Darwin':
    print("No openMP for you !")
if platform.system() == 'Linux':
    extra_comp_cython += ['-fopenmp']
    extra_link_cython += ['-fopenmp']

# Avoid Cython warning about NumPy API deprecation upon installation
if platform.system() == 'Darwin':
    extra_comp_cython += ['-Wno-#warnings']
if platform.system() == 'Linux':
    extra_comp_cython += ['-Wno-cpp']

# Avoid error "fatal error: 'complex' file not found" under OSX (Python 3.6)

if platform.system() == 'Darwin':
    if platform.python_version_tuple()[0:2] == ('3', '6'):
        extra_comp_cython += ["-stdlib=libc++"]
        extra_link_cython += ["-stdlib=libc++"]

# NOTE : the separation in several different modules was needed to make it
# work on MacOSX for some unfathomable reason.

cython_ologram_1 = Extension(name='pygtftk.stats.intersect.create_shuffles',
                             sources=["pygtftk/stats/intersect/create_shuffles.pyx"],
                             extra_compile_args=extra_comp_cython, extra_link_args=extra_link_cython,
                             language='c')

cython_ologram_2 = Extension(name='pygtftk.stats.intersect.overlap.overlap_regions',
                             sources=["pygtftk/stats/intersect/overlap/overlap_regions.pyx"],
                             extra_compile_args=extra_comp_cython, extra_link_args=extra_link_cython,
                             language='c')

cython_ologram_3 = Extension(name='pygtftk.stats.intersect.read_bed.read_bed_as_list',
                             sources=["pygtftk/stats/intersect/read_bed/read_bed_as_list.pyx",
                                      "pygtftk/stats/intersect/read_bed/exclude.cpp"],  # Include custom Cpp code
                             extra_compile_args=extra_comp_cython, extra_link_args=extra_link_cython,
                             include_dirs=[np.get_include()],
                             language='c++')

cython_ologram_4 = Extension(name='pygtftk.stats.multiprocessing.multiproc',
                             sources=["pygtftk/stats/multiprocessing/multiproc.pyx",
                                      "pygtftk/stats/multiprocessing/multiproc_structs.pxd",
                                      "pygtftk/stats/multiprocessing/multiproc.pxd"],
                             extra_compile_args=extra_comp_cython, extra_link_args=extra_link_cython,
                             language='c')

# ----------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------

long_description = open('README.rst', mode="r").read()

# ----------------------------------------------------------------------
# Update the docs/source/conf.py
# This will allow the doc to display the right version
# of the pygtftk library
# ----------------------------------------------------------------------

tmp_file_conf = NamedTemporaryFile(delete=False,
                                   mode="w",
                                   prefix="pygtftk_",
                                   suffix="_conf.py")

conf_file = open("docs/source/conf.py", "r")
for line in conf_file:
    if line.startswith("version = "):
        line = "version = " + "u'" + __version__ + "'\n"

    if line.startswith("release = "):
        line = "release = " + "u'" + __version__ + "'\n"

    tmp_file_conf.write(line)

tmp_file_conf.close()
conf_file.close()

# Using copy instead of move due to a bug both in
# os.rename and shutils.move...
# https://tinyurl.com/y8ghvgyq

os.remove(conf_file.name)
shutil.copy(tmp_file_conf.name, conf_file.name)
os.remove(tmp_file_conf.name)

# ----------------------------------------------------------------------
# Declare the setup function
# ----------------------------------------------------------------------

with open('requirements.txt') as f:
    pack_required = f.read().splitlines()

setup(name="pygtftk",
      include_dirs=[np.get_include()],
      version=__version__,
      author_email=__email__,
      author=__author__,
      description=__description__,
      url=__url__,
      zip_safe=False,
      project_urls={
          'Source': __url_source__,
          'Tracker': __url_tracker__
      },
      python_requires=__python_requires__,
      cmdclass={'build_ext': build_ext},
      keywords=__keywords__,
      packages=['pygtftk',
                'pygtftk/plugins',
                'docs',
                'docs/source',
                'pygtftk/bwig',
                'pygtftk/rtools',
                'pygtftk/stats',
                'pygtftk/stats/multiprocessing',
                'pygtftk/stats/intersect',
                'pygtftk/stats/intersect/overlap',
                'pygtftk/stats/intersect/read_bed',
                'pygtftk/stats/intersect/modl',
                'pygtftk/data',
                'pygtftk/data/simple',
                'pygtftk/data/simple_02',
                'pygtftk/data/simple_03',
                'pygtftk/data/simple_04',
                'pygtftk/data/simple_05',
                'pygtftk/data/simple_06',
                'pygtftk/data/simple_07',
                'pygtftk/data/mini_real',
                'pygtftk/data/mini_real_10M',
                'pygtftk/data/mini_real_noov_rnd_tx',
                'pygtftk/data/tiny_real',
                'pygtftk/data/hg38_chr1',
                'pygtftk/data/mini_real_ens',
                'pygtftk/data/control_list',
                'pygtftk/data/ologram_1',
                'pygtftk/data/ologram_2',
                'pygtftk/src/',
                'pygtftk/src/libgtftk',
                'pygtftk/src/libgtftk/command'],
      package_data={'pygtftk/data': ['*.*'],
                    'pygtftk/data/simple': ['*.*'],
                    'pygtftk/data/simple_02': ['*.*'],
                    'pygtftk/data/simple_03': ['*.*'],
                    'pygtftk/data/simple_04': ['*.*'],
                    'pygtftk/data/simple_05': ['*.*'],
                    'pygtftk/data/simple_06': ['*.*'],
                    'pygtftk/data/simple_07': ['*.*'],
                    'pygtftk/data/mini_real': ['*.*'],
                    'pygtftk/data/mini_real_10M': ['*.*'],
                    'pygtftk/data/mini_real_noov_rnd_tx': ['*.*'],
                    'pygtftk/data/tiny_real': ['*.*'],
                    'pygtftk/data/hg38_chr1': ['*.*'],
                    'pygtftk/data/control_list': ['*.*'],
                    'pygtftk/data/ologram_1': ['*.*'],
                    'pygtftk/data/ologram_2': ['*.*'],
                    'pygtftk/data/mini_real_ens': ['*.*'],
                    'pygtftk/plugins': ['*.*'],
                    'docs': ['Makefile'],
                    'docs/source': ['*.*'],
                    'pygtftk/stats': ['*.*'],
                    'pygtftk/stats/intersect': ['*.*'],
                    'pygtftk/stats/intersect/overlap': ['*.*'],
                    'pygtftk/stats/intersect/modl': ['*.*'],
                    'pygtftk/stats/multiprocessing': ['*.*'],
                    'pygtftk/stats/intersect/read_bed': ['*.*'],
                    'pygtftk/src': ['*.*'],
                    'pygtftk/src/libgtftk': ['*.*'],
                    'pygtftk/src/libgtftk/command': ['*.*']},
      scripts=['bin/gtftk'],
      license='LICENSE.txt',
      classifiers=__classifiers__,
      long_description=long_description,
      extras_require={
          'dev': ['pycodestyle >= 2.1.0',
                  'sphinx >=1.5.2',
                  'sphinxcontrib-programoutput >=0.8',
                  'sphinx_bootstrap_theme >=0.4.9',
                  'sphinxcontrib-googleanalytics'],
          'gffutils': ['gffutils']},
      install_requires=pack_required,
      ext_modules=[lib_pygtftk] + [cython_ologram_1, cython_ologram_2, cython_ologram_3, cython_ologram_4])

# ----------------------------------------------------------------------
# Update gtftk config directory
# ----------------------------------------------------------------------

current_install_path = subprocess.Popen(['which', 'gtftk'],
                                        stdout=subprocess.PIPE,
                                        stderr=DEVNULL
                                        ).stdout.read().rstrip()
if current_install_path == '':
    sys.stderr.write("Unable to find gtftk in the path...\n")

str_to_hash = current_install_path + __version__.encode()
config_dir_hash = hashlib.md5(str_to_hash).hexdigest()
config_dir = os.path.join(os.path.expanduser("~"),
                          ".gtftk",
                          config_dir_hash)

if os.path.exists(config_dir):
    shutil.rmtree(config_dir)
    sys.stderr.write("Deleting old configuration directory at: " + config_dir + "\n")

# ----------------------------------------------------------------------
# Print gtftk info (and load the plugins...)
# ----------------------------------------------------------------------
# Put this in a `try` as it could raise an exception if the program is not
# in the PATH.
try:
    gtftk_sys_config = subprocess.Popen(['gtftk', '-s'], stdout=subprocess.PIPE).stdout.read().rstrip()
    sys.stderr.write(gtftk_sys_config.decode())
except FileNotFoundError:
    pass

sys.stderr.write("\n\nInstallation complete.\n\n")
