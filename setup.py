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
import os
import re
import shutil
import sys
from distutils import sysconfig
from sys import platform
from tempfile import NamedTemporaryFile
import subprocess




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
# Delete any existing .gtftk in home folder
# -------------------------------------------------------------------------

gtftk_cnf_dir = os.path.join(os.environ['HOME'], '.gtftk')
sys.stderr.write("Trying to find an existing .gtftk directory\n")
if os.path.exists(gtftk_cnf_dir):
    sys.stderr.write("Found an existing .gtftk directory\n")
    shutil.rmtree(gtftk_cnf_dir, ignore_errors=True)
    if os.path.exists(gtftk_cnf_dir):
        sys.stderr.write("Failed to delete .gtftk directory.\n")
else:
    sys.stderr.write("No .gtftk directory found.\n")
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
# Delete the first line from REAME.md
# and convert .md to .rst...
# This is required for Pypi that use rst format
# and will read the content of the long_description key
# ----------------------------------------------------------------------
# I just realized that github also support rst...
# This part should be deleted in the future.

long_description_file = open('README.md')
long_description = []
markup_char = {1: "=", 2: "-", 3: "~"}
markup_level = 0
past_line_len = 0

start_of_text = False

for pos, line in enumerate(long_description_file):
    if not line.startswith("    "):
        line = line.rstrip('\r\n')
    if line.startswith("# ") and start_of_text is False:
        start_of_text = True

    if start_of_text:
        # Replace title
        if markup_level > 0:
            title = markup_char[markup_level] * past_line_len
            long_description += [title]
        len_line = len(line)
        line = re.sub("^#+", "", line)
        markup_level = len_line - len(line)
        past_line_len = len(line)

        # replace URL
        for hit in re.finditer("\[(.*?)\]\((.*?)\)", line):
            line = line.replace("[" + hit.group(1) + "]", "`" + hit.group(1) + " <")
            line = line.replace("(" + hit.group(2) + ")", hit.group(2) + ">`_")

        if markup_level:
            long_description += [line.lstrip(" ")]
        else:
            long_description += [line]

long_description = "\n".join(long_description)
long_description_rst = open('README.rst', mode="w")
long_description_rst.write(long_description)

# ----------------------------------------------------------------------
# Update the docs/manual/source/conf.py
# This will allow the doc to display the right version
# of the pygtftk library
# ----------------------------------------------------------------------

tmp_file_conf = NamedTemporaryFile(delete=False,
                                   mode="w",
                                   prefix="pygtftk_",
                                   suffix="_conf.py")

conf_file = open("docs/manual/source/conf.py", "r")
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
                'docs/manual',
                'docs/manual/build',
                'docs/manual/build/html',
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
                    'pygtftk/src': ['*.*'],
                    'pygtftk/src/libgtftk': ['*.*'],
                    'pygtftk/src/libgtftk/command': ['*.*'],
                    'docs/manual/build/html': ['*.*']},
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

config_dir = os.path.join(os.path.expanduser("~"), ".gtftk")
dumped_plugin_path = os.path.join(config_dir, "plugin.pick")

sys.stderr.write("Removing old plugin dump.\n")

try:
    os.remove(dumped_plugin_path)
except OSError as e:
    pass

sys.stderr.write("Installation complete.\n")
result = subprocess.run(['gtftk', '-h'], stdout=subprocess.PIPE)
sys.stderr.write(result.stdout.decode() + "\n")
