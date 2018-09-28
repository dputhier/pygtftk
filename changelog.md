# Changelog

## v0.9.7 (not yet released)

### Bug Fixes

- Fixed a bug in add_attr_from_file. The arg has_header led to empty result.
- Fixed a bug in select_by_reg_exp. The match method was used instead of search method...   	
	
### API Changes

### New Features

- Several changes in doc to make it compliant with readthedoc.


# Changelog

## v0.9.6

### Bug Fixes

### API Changes

### New Features

- Several changes in setup.py and requirements_develops.txt to makes it compliant with Pypi and buildable with manylinux.


## v0.9.5

### Bug Fixes

### API Changes

### New Features

- Several changes in setup.py and requirements_develops.txt to makes it compliant with Pypi and buildable with manylinux.

## v0.9.3

### Bug Fixes

- User-defined colors are now applied when calling profile command (bug in 0.9.2).
- Fixed test issues when using md5sum under Linux.
- Lots of typo fixed (their must be lot remaining unfortunatly...).

### API Changes


### New Features

- pygtftk is now compatible with both python 2.7 and >=3.5.
- libstdcxx is no more required in env.yaml.
- 'Makefile clean' also reset version.py and conf.py to repository version.
- Improved parser loading.
- 5p_3p_coord has been renamed get_5p_3p_coords.
- Added a pip requirement.txt file for developers.
- Added a specific conda environment file for py2k.
- The default conda environment targets Python 3.6.
- Improve garbage collector behavior upon exit.
- Added --system-info argument to gtftk.
