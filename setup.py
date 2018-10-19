"""
The pygtfk package.

The Python GTF toolkit (pygtftk) package is intented to ease handling of GTF (Gene Transfer Format) files. The pygtftk package relies on libgtftk, a library of functions written in C.
The package comes with a set of UNIX commands that can be accessed through the gtftk program. The gtftk program proposes several atomic tools to filter, convert, or extract data from GTF files. The gtftk set of Unix commands can be easily extended using a basic plugin architecture. All these aspects are covered in the help section.

Authors: D. Puthier and F. Lopez
"""

# -------------------------------------------------------------------------
# A set of builtin packages
# -------------------------------------------------------------------------

import glob
import hashlib
import os
import re
import shutil
import sys
from distutils import sysconfig
from sys import platform
from tempfile import NamedTemporaryFile
import subprocess

#try:
#    out=os.popen("ls -R /home/docs/checkouts/readthedocs.org/user_builds/pygtftk/")
#    print(out.read())
#except:
#    print("Unable to print the test")

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
# Python Version
# -------------------------------------------------------------------------

PY3 = sys.version_info[0] == 3
PY2 = sys.version_info[0] == 2


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

    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    sha = repo.git.rev_parse(sha, short=4)

    if sha != "" and not os.path.exists("release_in_progress"):
        __version__ = base_version + ".dev0+" + sha
    else:
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

if platform == "darwin":
    vars = sysconfig.get_config_vars()
    # vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')
    # dyn_lib_compil = ['-dynamiclib', '-shared']
    dyn_lib_compil = []
else:
    dyn_lib_compil = ['-shared']

extra_compile_args = ['-Ipygtftk/src/libgtftk',
                      '-O3',
                      '-Wall',
                      '-fPIC',
                      '-MMD',
                      '-MP',
                      '-fmessage-length=0'] + dyn_lib_compil

lib_pygtftk = Extension(name='pygtftk/lib/libgtftk',
                        include_dirs=[
                            'pygtftk/src/libgtftk'],
                        library_dirs=['/usr/lib'],
                        libraries=['z'],
                        extra_compile_args=extra_compile_args,
                        sources=cmd_src_list)

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


setup(name="pygtftk",
      version=__version__,
      author_email='denis.puthier@univ-amu.fr',
      author="fabrice Lopez,Denis Puthier",
      description="The Python GTF toolkit (pygtftk) package: easy handling of GTF files",
      url="https://github.com/dputhier/pygtftk",
      zip_safe=False,
      project_urls={
          'Source': 'https://github.com/dputhier/pygtftk',
          'Tracker': 'https://github.com/dputhier/pygtftk/issues'
      },
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <3.7',
      keywords="genomics bioinformatics GTF BED",
      packages=['pygtftk',
                'pygtftk/plugins',
                'docs',
                'docs/source',
                'pygtftk/bwig',
                'pygtftk/rtools',
                'pygtftk/data',
                'pygtftk/data/simple',
                'pygtftk/data/simple_02',
                'pygtftk/data/simple_03',
                'pygtftk/data/simple_04',
                'pygtftk/data/simple_05',
                'pygtftk/data/mini_real',
                'pygtftk/data/mini_real_noov_rnd_tx',
                'pygtftk/data/control_list',
                'pygtftk/src/',
                'pygtftk/src/libgtftk',
                'pygtftk/src/libgtftk/command'],
      package_data={'pygtftk/data': ['*.*'],
                    'pygtftk/data/simple': ['*.*'],
                    'pygtftk/data/simple_02': ['*.*'],
                    'pygtftk/data/simple_03': ['*.*'],
                    'pygtftk/data/simple_04': ['*.*'],
                    'pygtftk/data/simple_05': ['*.*'],
                    'pygtftk/data/mini_real': ['*.*'],
                    'pygtftk/data/mini_real_noov_rnd_tx': ['*.*'],
                    'pygtftk/data/control_list': ['*.*'],
                    'pygtftk/plugins': ['*.*'],
                    'docs': ['Makefile'],
                    'docs/source': ['*.*'],
                    'pygtftk/src': ['*.*'],
                    'pygtftk/src/libgtftk': ['*.*'],
                    'pygtftk/src/libgtftk/command': ['*.*']},
      scripts=['bin/gtftk'],
      license='LICENSE.txt',

      classifiers=("License :: OSI Approved :: MIT License",
                   "Operating System :: MacOS",
                   "Operating System :: POSIX :: Linux",
                   "Development Status :: 4 - Beta",
                   "Environment :: Console",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3.5",
                   "Programming Language :: Python :: 3.6",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Topic :: Scientific/Engineering :: Bio-Informatics",
                   "Operating System :: POSIX :: Linux",
                   "Operating System :: MacOS",
                   "Topic :: Documentation :: Sphinx"
                   ),
      long_description=long_description,
      extras_require={
          'tests': [
              'nose',
              'pycodestyle >= 2.1.0'],
          'docs': [
              'sphinx >=1.5.2',
              'sphinxcontrib-programoutput >=0.8',
              'sphinx_bootstrap_theme >=0.4.9']},
      install_requires=['pyyaml >=3.12',
                        'argparse',
                        'cloudpickle >=0.5.6',
                        'ftputil >=3.3.1',
                        'pybedtools >=0.7.8',
                        'pandas >=0.23.3',
                        'requests >=2.13.0',
                        'pyBigWig >=0.3.12',
                        'cffi >=1.10.0',
                        'biopython >=1.69',
                        'pyparsing >=2.2.0',
                        'GitPython >=2.1.8',
                        'pyparsing',
                        'matplotlib >=2.0.2',
                        'plotnine >=0.4.0',
                        'future',
                        'setuptools'],
      ext_modules=[lib_pygtftk])


# ----------------------------------------------------------------------
# Update gtftk config directory
# ----------------------------------------------------------------------

current_install_path = subprocess.Popen(['which', 'gtftk'], stdout=subprocess.PIPE).stdout.read().rstrip()
if current_install_path == '':
    sys.stderr.write("Unable to find gtftk in the path...\n")
    sys.exit(1)

config_dir_hash = hashlib.md5(os.path.abspath(current_install_path)).hexdigest()
config_dir = os.path.join(os.path.expanduser("~"),
                                     ".gtftk",
                                     config_dir_hash)


if os.path.exists(config_dir):
    shutil.rmtree(config_dir)
    sys.stderr.write("Deleting old configuration directory at: " + config_dir + "\n")

# ----------------------------------------------------------------------
# Print gtftk info (and load the plugins...)
# ----------------------------------------------------------------------

gtftk_sys_config = subprocess.Popen(['gtftk', '-s'], stdout=subprocess.PIPE).stdout.read().rstrip()
sys.stderr.write(gtftk_sys_config.decode())
sys.stderr.write("\n\nInstallation complete.\n")