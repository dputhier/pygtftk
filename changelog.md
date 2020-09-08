# Changelog


## v1.1.5

This version introduces OLOGRAM-MODL, a new paradigm for OLOGRAM to find intersections between multiple sets of genomic regions at once and then compute their enrichment with OLOGRAM. An optional algorithm (MODL) to find interesting combinations with sparse dictionary learning and greedy submodular optimisation has also been added. Furthermore, it also contains major speedups to OLOGRAM itself.

### Bug Fixes

*   OLOGRAM - Fixed bug in multiple overlaps when using sets with no peak on a given chromosome.
*   Updated package requirements
*   OLOGRAM - The display graph has now the X labels at 90 degrees.
*   OLOGRAM - more-bed-labels should take and clean BED file names as default
*   OLOGRAM - Various graphical fixes
*   fix #124
*   fix BED to BED convertion in arg_formatted.FormattedFile(). BED6+ files were considered as BED6- files.
*   fix #136 although --show-group-number is no more supported with gtftk profile when plotnine > 0.6.0 is used.

### API Changes

*   Moved OLOGRAM-related commands to their own section in the documentation.
*   The MODL algorithm for combination mining can be accessed independantly.

### Code changes

*   Major speedups achieved in OLOGRAM by better typecasting in the Cython code.
*   Major speedup in OLOGRAM due to rewriting the pandas melt() function in C/Cython.
*   Added multithreading batch-by-batch for OLOGRAM
*   Renamed *merge_ologram_stats* to *ologram_merge_stats*.
*   Improved *ologram_merge_stats* visuals.
*   Added new *simple_07* and *ologram_2* example datasets to study multiple overlaps.
*   Added scikit-learn as a dependency.
*   Moved OLOGRAm functions to calculate enrichment to their own module.

### New Features

*   This version implements OLOGRAM-MODL to study the enrichment of intersections between multiple sets of genomic regions. Please see the documentation and code comments for more details.
*   The API contains a Modl class which is a dictionary-learning based itemset mining algorithm, used in OLOGRAM-MODL
*   Introduced a *treeify_ologram_modl* plugin to visualize n-wise enrichment results as a treee
*   Introduced a *ologram_merge_runs* command to merge several runs to save RAM, treating each as a superbatch.



## v1.1.4

### Bug Fixes

*   None

### API/CLI Changes

*   No more compatible with Python 3.5 (as BioPython).

### Code changes

*   None.

### New Features

*   None.

## v1.1.3


### Bug Fixes

*   fix #123, #122, #121, #120

### API/CLI Changes

*   None.

### Code changes

*   None.

### New Features

*   None.

## v1.1.2


### Bug Fixes

*   None

### API/CLI Changes

*   None.

### Code changes

*   None.

### New Features

*   The --more-bed-labels is now facultative in OLOGRAM.


## v1.1.1

### Bug Fixes

*   Fix #116 (pandas version issue)
*   Fix an issue related to pybedtool/bedtool version (naming of sequences that differs due to name/name+/nameOnly arguments).
*   Fix -n with integer values in get_5p_3p_coords.

### API/CLI Changes

*   None.

### Code changes

*   None.

### New Features

* md5sum-lite call have been replaced by "md5 -r" under darwin platforms.
* The tss_numbering command now allows to add the number of different TSSs to the gene feature.

## v1.1.0


### Bug Fixes

*   None.

### API/CLI Changes

*   None.

### Code changes

*   None.

### New Features

* Support for Python 3.7.
* The tss_numbering command now allows to add the number of different TSSs to the gene feature.


## v1.0.9


### Bug Fixes

*   None.

### API/CLI Changes

*   None.
*   bigwig_to_bed is know part of miscellaneous commands.

### Code changes

*   The select_by_key command now accept seq_name and seqname in addition to seqid and chrom as a key.

### New Features

*  This version contains the tss_numbering plugin. Annotate transcripts by computing their TSS position relative to the most five prime TSS of the corresponding gene.
*  The get_attr_value_list command now accepts a list of keys as input.
*  The get_attr_value_list command has an additional argument (--print-key-name).
*  Added a new miscellaneous command, get_ceas_record that extract records from CEAS (Cis-regulatory Element Annotation System).
*  Added a new miscellaneous command, great_reg_domains. This tool represents an attempt to process genomic annotations in GTF
format in order to produced a set of 'labeled' regions with the same rules as those described in GREAT (Genomic Regions
Enrichment of Annotations Tool) documentation. We can not warrant that the procedure is exactly the same. See the CLI
for more details.
* Added -y/--display-fit-quality to ologram


## v1.0.8

This version introduces *ologram_merge_stats* command that can be used to produce a heatmap from multiple OLOGRAM results.

### Bug Fixes

*   None.

### API Changes

*   None.

### Code changes

*   The select_by_key command now accept seq_name and seqname in addition to seqid and chrom as a key.

### New Features

*   This version implements ologram_merge_stats command that can be used to produce a heatmap from multiple OLOGRAM results.

## v1.0.7

This version contains some minor code refactoring. See 1.0.6 for recent major changes.

### Bug Fixes

*   None.

### API Changes

*   None.

### Code changes

*   Essentially some refactoring.

### New Features

*   None.


## v1.0.6

### Bug Fixes

*   *ologram* warnings are now clearer.
*   Fixed reproducibility issue when using --more_keys in *ologram*.
*   Fix #78 (Relaxing constraints on GTF format)

### API Changes

### Code changes

  * Added new tests to *overlapping*.
  * Excluding regions in *ologram* is now done in C++, bringing major speedups (400x).

### New Features

*   Added -b/--bool argument to *overlapping*.
*   Added -@/--annotate-all to *overlapping*.
*   Added bigwig_to_bed.
*   Added --bed-incl argument to *ologram*. (fix #73)
*   The --dpi has been removed from *ologram*.
*   The --pval-precision has been removed from *ologram*.
*   The --user-img-file argument in *ologram* has become pdf_file.
*   Added -W as command-wise argument. It force gtftk to print its
    message to a file. This may be handy when no access to stderr
    is available (e.g. through a scheduler).

## v1.0.5

This version provides several improvements and bug fix to ologram. The CLI of *mk_matrix*
and *join_multi_file* have slightly changed. An int ([0-4]) is now mandatory to control
verbosity level.

### Bug Fixes

*   Fix #76 (issue with chromosomes in ensembl format when using ologram).
*   Fix #72 (P-values equal to zero when nb_intersections_true and
    nb_intersections_expectation_shuffled are equal to zero).

### API Changes

### Code changes

*   Improved p-value precision in ologram.

### New Features

*   Added -j/--sort-features to ologram. Controls the feature sorting in
    barplot diagram. Changed doc accordingly.
*   The positional argument 'bigwiglist' in mk_matrix has been replaced by -y.
*   The verbosity argument now must take a value (0-4).
*   The positional argument in join_multi_file has been replaced by -m.

## v1.0.4

*   Updated changelog

### Bug Fixes

### API Changes

### Code changes

### New Features

## v1.0.3

*   Improved p-value precision in ologram.

### Bug Fixes

### API Changes

### Code changes

*   Improved p-value precision in ologram.
*   Added mpmath as dependency.

### New Features

*   Added tiny_real dataset.
*   Added -q argument (--quiet) to get_example.

## v1.0.2

*   No more using readthedoc as the doc is becoming too long to process.

### Bug Fixes

### API Changes

### Code changes

*   Many typo detected and fixed.

### New Features

*   Added tiny_real dataset.
*   Added -q argument (--quiet) to get_example.

## v1.0.0

### Bug Fixes

*   Fix --no-strandness in divergent.

### API Changes

### Code changes

*   Many typo detected and fixed.

### New Features

*   This version now integrates ologram (OverLap Of Genomic Regions Analysis using
    Monte Carlo). Ologram annotates peaks(in BED format) with region sets/features
    extracted from (i) GTF file features (e.g promoter, tts, gene body, UTR...) (ii)
    GTF file keys (e.g. gene_biotype, user defined keys...) (iii) or from a BED file.
*   The user can now use --chrom-info to provide the command with a file or a string.
    The string should be one of 'mm8', 'mm9', 'mm10', 'hg19', 'hg38', 'rn3' or 'rn4'.
    When a genome version is requested as a string, the conventional chromosomes
    are used (chrM is discarded together with alternative haplotypes,
    unlocalized regions...).

## v0.9.10

### Bug Fixes

*   Argparse was part of the dependencies. However, argparse is part of Python 3.
    Thus, this caused pygtftk to come with an older version of argparse...
*   Fixed gene sorting in tss_dict to ensure reproducible result.
*   Fixed a problem with retrieve() when used from interpreter (#45).
*   The load_gtf() function (C API) no raise an error if a GFF3 is passed
    with .gtf extension.

### API Changes

*   Input BED file in bed3 format are now converted to bed6 automatically.
*   The select_by_numeric() function has been renamed eval_numeric()
*   It is now possible to use numpy array of booleans to index the GTF (i.e.
    using the indexing function).
*   the prepare_gffutils_db() function allows one to create a db for
    gffutils while selecting features and attributes.

### Code changes

*   The argformatter module was refactored. Development of FormattedFile(argparse.FileType)
    that test for file extension and content (at least for bed).
*   The BED conversion is now performed using the print_bed() C function .

### New Features

## v0.9.9

### Bug Fixes

*   Fix a critical bug in get_sequence that affected get_feat_seq and get_tx_seq.
*   Select_by_key now throw an error when no key/val are available.
*   No more function with mutable objects as default arguments.
*   Fix temporary file deletion.

### API Changes

*   Refactored arg_formatter by  creating a single type (ranged_num)
to test for numeric inputs.
*   Refactored all plugins so that there is no more reference to
unused arguments (tmp_dir, verbosity...).

### Code changes

*   No more reference to PY2.
*   Added several test to get_tx_seq and get_feat_seq.
*   Added several script to manipulate fasta files (see 'tools' folder).
    For pygtftk dev.
*   Added 'extra_require' slot in setup().
*   The get-feature-seq program now relies on bedtools (not on internal C code).
    This may change in the future as a a more
    flexible C interface is available.

### New Features

*   Added --list-bigwigs to profile (to display the content of a coverage file).
*   Added a novel dataset	(mini_real_10M) derived from mini_real and containing
    10 Mb of chr1.
*   The configuration directory now supports several subdirectories named based
    on a hash string computed from path to the gtftk program.

## v0.9.8

### Bug Fixes

*   Convergent and divergent return coordinate (e.g. dist_to_divergent)
in integer not float...

### API Changes

*   pygtftk is no more compatible with python 2. This decision aims at
    integrating last plotnine versions (starting from 0.5.1) that
    depends on matplotlib 3.0.0 which strictly depends on py3k.

### New Features

*   Added an example for col_from_tab in presentation.rst
*   Added a dataset mini_real_coding_pot.tab.
*   Working on travis now.
*   gtftk configuration directory now contains several subdirectories whose
    names are computed based on gtftk program location.
*   Added a -d argument to gtftk program. This argument returns gtftk
    configuration directory.
*   All tests should be independent of the directory.
*   Added test to count_key_values.

## v0.9.7

### Bug Fixes

*   Fixed a bug in add_attr_from_file. The arg has_header led to empty result.
*   Fixed a bug in select_by_reg_exp. The match method was used instead of search method...
*   Fixed an error in count_key_values.

### API Changes

### New Features

*   Several changes in doc to make it compliant with readthedoc.

## v0.9.6

### Bug Fixes

### API Changes

### New Features

*   Several changes in setup.py and requirements_develops.txt to makes it
    compliant with Pypi and buildable with manylinux.

## v0.9.5

### Bug Fixes

### API Changes

### New Features

*   Several changes in setup.py and requirements_develops.txt to makes it compliant
    with Pypi and buildable with manylinux.

## v0.9.3

### Bug Fixes

*   User-defined colors are now applied when calling profile command (bug in 0.9.2).
*   Fixed test issues when using md5sum under Linux.
*   Lots of typo fixed (their must be lot remaining unfortunatly...).

### API Changes

### New Features

*   pygtftk is now compatible with both python 2.7 and >=3.5.
*   libstdcxx is no more required in env.yaml.
*   'Makefile clean' also reset version.py and conf.py to repository version.
*   Improved parser loading.
*   5p_3p_coord has been renamed get_5p_3p_coords.
*   Added a pip requirement.txt file for developers.
*   Added a specific conda environment file for py2k.
*   The default conda environment targets Python 3.6.
*   Improve garbage collector behavior upon exit.
*   Added --system-info argument to gtftk.
