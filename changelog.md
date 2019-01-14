# Changelog

## v0.9.10

### Bug Fixes

- Argparse was part of the dependencies. However, argparse is part of Python 3. Thus, this caused pygtftk to come with an older version of argparse...
- Fixed gene sorting in tss_dict to ensure reproducible result.
- Fixed a problem with retrieve() when used from interpreter (#45).
- The load_gtf() function (C API) no raise an error if a GFF3 is passed with .gtf extension.


### API Changes

- Input BED file in bed3 format are now converted to bed6 automatically.
- The select_by_numeric() function has been renamed eval_numeric()
- It is now possible to use numpy array of boolean to index the GTF (i.e. using the indexing function).	
- the prepare_gffutils_db() function allows one to create a db for gffutils while selecting features and attributes.

### Code changes

- The argformatter module was refactored. Development of FormattedFile(argparse.FileType) that test for file extension and content (at least for bed).
- The BED convertion is now performed using the C function print_bed().

### New Features



## v0.9.9

### Bug Fixes

- Fix a critical bug in get_sequence that affected get_feat_seq and get_tx_seq.
- Select_by_key now throw an error when no key/val are available.
- No more function with mutable objects as default arguments.
- Fix temporary file deletion. 

### API Changes

- Refactored arg_formatter by  creating a single type (ranged_num) to test for numeric inputs.
- Refactored all plugins so that there is no more reference to unused arguments (tmp_dir, verbosity...).


### Code changes

- No more reference to PY2.
- Added several test to get_tx_seq and get_feat_seq.
- Added several script to manipulate fasta files (see 'tools' folder). For pygtftk dev.
- Added 'extra_require' slot in setup().
- The get-feature-seq program now relies on bedtools (not on internal C code). This may change in the future asa a more flexible C interface is available.

### New Features

- Added --list-bigwigs to profile (to display the content of a coverage file).
- Added a novel dataset	(mini_real_10M) derived from mini_real and containing 10 Mb of chr1.
- The configuration directory now supports several subdirectories named based on a hash string computed from path to the gtftk program.


## v0.9.8

### Bug Fixes

- Convergent and divergent return coordinate (e.g. dist_to_divergent) in integer not float...

### API Changes

- pygtftk is no more compatible with python 2. This decision aims at integrating last plotnine versions (starting from 0.5.1) that depends on matplotlib 3.0.0 which strictly depends on py3k.

### New Features

- Added an example for col_from_tab in presentation.rst
- Added a dataset mini_real_coding_pot.tab.
- Working on travis now.
- gtftk configuration directory now contains several subdirectories whose names are computed based on gtftk program location.
- Added a -d argument to gtftk program. This argument returns gtftk configuration directory. 
- All tests should be independent of the directory.
- Added test to count_key_values.

## v0.9.7

### Bug Fixes

- Fixed a bug in add_attr_from_file. The arg has_header led to empty result.
- Fixed a bug in select_by_reg_exp. The match method was used instead of search method...
- Fixed an error in count_key_values. 

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
