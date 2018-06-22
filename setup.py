"""
The setup.py file of the gtfk package.
#TODO
"""

import glob
import os
import sys

from sys import platform

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
    sha = repo.head.object.hexsha
    sha = repo.git.rev_parse(sha, short=4)

except:
    sha = ""

if sha != "":
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
    dyn_lib_compil = ['-dynamiclib', '-shared']
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

setup(name="pygtftk",
      version=__version__,
      author_email='puthier@gmail.com',
      author="Denis Puthier",
      description="Genomic tool suite",
      url="",
      zip_safe=False,
      keywords="genomics",
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
      long_description=open('README.md').read(),
      install_requires=['pyyaml', 'argparse', 'cloudpickle',
                        'ftputil', 'pybedtools', 'pandas', 'pyBigWig',
                        'requests', 'cffi', 'biopython', 'pyparsing'],
      classifiers=[
          ""
      ],
      ext_modules=[lib_pygtftk])

config_dir = os.path.join(os.path.expanduser("~"), ".gtftk")
dumped_plugin_path = os.path.join(config_dir, "plugin.pick")

sys.stderr.write("Removing old plugin dump.\n")

try:
    os.remove(dumped_plugin_path)
except OSError as e:
    pass

sys.stderr.write("Installation complete.\n")
