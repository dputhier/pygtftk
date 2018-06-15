"""
The setup.py file of the gtfk package.
"""

import glob
import os
import sys
import re
from distutils import sysconfig
from pygtftk.utils import chomp
from sys import platform
from pygtftk.utils import  make_tmp_file
import git
from setuptools import setup, Extension

from pygtftk.version import __base_version__

try:
    from setuptools import setup
except:
    pass

# -------------------------------------------------------------------------
# Check gtftk version
# -------------------------------------------------------------------------

try:
    repo = git.Repo(search_parent_directories=True)
    branch = repo.active_branch
    sha = repo.head.object.hexsha
    sha = repo.git.rev_parse(sha, short=4)

except:
    sha = ""

if sha != "" and branch != "master" and not os.path.exists("pypi_release_in_progress"):
    __version__ = __base_version__ + ".dev0+" + sha
else:
    __version__ = __base_version__

version_file = open('pygtftk/version.py', "w")
version_file.write("__base_version__='" + __base_version__ + "'\n")
version_file.write("__version__='" + __version__ + "'\n")
version_file.close()

# -------------------------------------------------------------------------
# Building C library
# -------------------------------------------------------------------------

cmd_src_list = glob.glob("pygtftk/src/libgtftk/command/*.c")

if platform == "darwin":
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')
    dyn_lib_compil = ['-shared']
else:
    dyn_lib_compil = ['-shared']

extra_compile_args = ['-Ipygtftk/src/libgtftk',
                      '-O3',
                      '-Wall',
                      '-fPIC',
                      '-shared',
                      '-MMD',
                      '-MP',
                      '-fmessage-length=0'] + dyn_lib_compil

lib_pygtftk = Extension(name='pygtftk/lib/libgtftk',
                      include_dirs=[
                          'pygtftk/src/libgtftk'],
                      library_dirs=['/usr/lib'],
                      libraries=['z'],
                      extra_compile_args=extra_compile_args,
                      sources=cmd_src_list + ['pygtftk/src/libgtftk/column.c',
                                              'pygtftk/src/libgtftk/gtf_reader.c',
                                              'pygtftk/src/libgtftk/libgtftk.c'])


# ----------------------------------------------------------------------
# Delete the first line from REAME.md
# and convert .md to .rst...
# This is required for Pypi that use rst format
# and will read the content of the long_description key
# ----------------------------------------------------------------------



long_description_file = open('README.md')
long_description = []
markup_char = {1:"=", 2:"-", 3:"~"}
markup_level = 0
past_line_len = 0

for pos,line in enumerate(long_description_file):
    if not line.startswith("    "):
        line = chomp(line)

    if pos > 0:
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
            line= line.replace("[" + hit.group(1) + "]", "`" + hit.group(1) + " <")
            line = line.replace("(" + hit.group(2) + ")", hit.group(2) + ">`_")

        if markup_level:
            long_description += [line.lstrip(" ")]
        else:
            long_description += [line]

long_description = "\n".join(long_description)

# ----------------------------------------------------------------------
# Update the docs/manual/source/conf.py
# This will allow the doc to display the right version
# of the pygtftk library
# ----------------------------------------------------------------------





tmp_file_conf = make_tmp_file(prefix="pygtftk_", suffix="_conf.py")

conf_file = open("docs/manual/source/conf.py", "r")
for line in conf_file:
    if line.startswith("version = "):
        line = "version = " + "u'" + __version__ + "'\n"

    if line.startswith("release = "):
        line = "release = " + "u'" + __version__ + "'\n"


    tmp_file_conf.write(line)

tmp_file_conf.close()
conf_file.close()

os.remove(conf_file.name)
os.rename(tmp_file_conf.name, conf_file.name)

# ----------------------------------------------------------------------
# Declare the setup function
# ----------------------------------------------------------------------

setup(name="pygtftk",
      version=__version__,
      author_email='fabrice.lopez@inserm.fr,denis.puthier@univ-amu.fr',
      author="fabrice Lopez,Denis Puthier",
      description="The Python GTF toolkit (pygtftk) package: easy handling of GTF files",
      url="https://github.com/dputhier/pygtftk",
      zip_safe=False,
      project_urls={
        'Source': 'https://github.com/dputhier/pygtftk',
        'Tracker': 'https://github.com/dputhier/pygtftk/issues'
      },
      python_requires='~=2.7',
      keywords="genomics bioinformatics GTF BED",
      packages=['pygtftk',
                'pygtftk/plugins',
                'pygtftk/bwig',
                'pygtftk/rtools',
                'pygtftk/data',
                'pygtftk/data/simple',
                'pygtftk/data/simple_02',
                'pygtftk/data/simple_03',
                'pygtftk/data/simple_04',
                'pygtftk/data/simple_05',
                'pygtftk/data/mini_real',
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
                    'pygtftk/plugins': ['*.*'],
                    'pygtftk/src': ['*.*'],
                    'pygtftk/src/libgtftk': ['*.*'],
                    'pygtftk/src/libgtftk/command': ['*.*']},
      scripts=['bin/gtftk'],
      license='LICENSE.txt',

      classifiers=(
          "Programming Language :: Python :: 2.7",
          "License :: OSI Approved :: MIT License",
          "Operating System :: MacOS",
          "Operating System :: POSIX :: Linux",
          "Development Status :: 4 - Beta",
          "Environment :: Console"
      ),
      long_description=long_description,
      install_requires=['pyyaml', 'argparse', 'cloudpickle',
                        'ftputil', 'pybedtools', 'pandas', 'pyBigWig',
                        'requests', 'cffi', 'biopython', 'pyparsing'],
      ext_modules=[lib_pygtftk])

config_dir = os.path.join(os.path.expanduser("~"), ".gtftk")
dumped_plugin_path = os.path.join(config_dir, "plugin.pick")

sys.stderr.write("Removing old plugin dump.\n")

try:
    os.remove(dumped_plugin_path)
except OSError as e:
    pass

sys.stderr.write("Installation complete.\n")

