#!/usr/bin/env python

import argparse
import os
import re
import sys
import time
import warnings
from functools import partial

import numpy as np
import pandas as pd
import pybedtools
from plotnine import (ggplot, aes, position_dodge,
                      geom_bar, ylab, theme, element_blank, element_text, geom_text, geom_errorbar)
from plotnine.ggplot import save_as_pdf_pages

from pygtftk import arg_formatter
from pygtftk.bedtool_extension import BedTool
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.stats.intersect import read_bed_as_list as read_bed  # Only used here for exclusions
from pygtftk.stats.intersect.overlap_stats_shuffling import \
    compute_overlap_stats  # Main function from the stats.intersect module
from pygtftk.utils import chrom_info_as_dict
from pygtftk.utils import close_properly
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import message

__updated__ = "2019-01-25"
__doc__ = """
 Annotate peaks (in bed format) with region sets/features computed on the
 fly from a GTF file  (e.g promoter, tts, gene body, UTR...). Custom features
 are supported.

 Each couple peak file/feature is randomly shuffled across the genome (inter-region
 lengths are considered). Then the probability of intersection under the null
 hypothesis (the peaks and this feature are independant) is deduced thanks to
 this Monte Carlo approach.

 Authors : Quentin Ferré <quentin.q.ferre@gmail.com> and Denis Puthier <denis.puthier@univ-amu.fr>
 """

__notes__ = """
 -- Genome size is computed from the provided chromInfo file (-c). It should thus only contain ordinary chromosomes.
 -- -\-chrom-info may also accept 'mm8', 'mm9', 'mm10', 'hg19', 'hg38', 'rn3' or 'rn4'. In this case the corresponding 
 size of conventional chromosomes are used. ChrM is not used.  
 -- The program produces a pdf file and a txt file ('_stats_') containing intersection statistics
 for the shuffled BEDs under H0 (peak_file and the considered genomic region are independant):
 number of intersections (= number of lines in the bed intersect) and total number of overlapping
 base pairs.
 The output figure gives, for both statistics, esperance and standard deviation (error bars)
 in the shuffles compared to the actual values.
    It also gives, under the 'fit' label for each statistic, the goodness of fit of the statistic under (H0)
 to a Negative Binomial assessed by a Cramer's V score (fit_quality gives 1-V ; as per Cramer (1948) a good fit should have a fit quality above 1 - 0.25 = 0.75).
    The p-value of the true intersection under the distribution characterized by the shuffles is also given, under 'p_value'.
    Finally, the log2 fold change between true and shuffles is also given.

 -- If -\-more-keys is used additional region sets will be tested based on the associated key value.
 As an example, if -\-more-keys is set to the 'gene_biotype' (a key generally found in ensembl GTF), the
 region related to 'protein_coding', 'lncRNA' or any other values for that key will be retrieved merged and tested
 for enrichment.

 -- Use -\-no-basic-feature if you want to perform enrichment analysis on focused annotations only (-\-more-bed or -\-more-key).


 -- BETA : The lists of region and inter-region lengths can be shuffled independantly, or by using two independant Markov models
 of order 2 respectively for each. This is not recommended in the general case and can *very* time-consuming (hours).

 -- The goal of the minibatch is to save RAM. Increase the number of minibatches, instead of their size.
 You may need to use very small minibatches if you have large sets of regions.

 -- You can exclude regions from the shuffling. This is done by shuffling across a concatenated "sub-genome" obtained by removing
 the excluded regions, but the same ones will be excluded from the peak_file and the GTF.
 This in Beta for now and will be very time-consuming (hours), especially if you have few CPU cores.
 Try using an exclusion file that is as small (around a thousand elements) as possible.

 -- Although peak_anno itself is not RAM-intensive, base pygtftk processing of a full human GTF can require upwards of 8Gb.
 It is recommended you do not run other programs in the meantime.

 -- If you are using the --no-basic-features argument *without* --more-keys, you can supply an empty file as the GTF, since it will be disregarded in the code.
 """


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputdir',
                            help='Output directory name.',
                            metavar="DIR",
                            default="peak_annotation",
                            type=str)

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Tabulated two-columns file. "
                                 "Chromosomes as column 1, sizes as column 2",
                            default=None,
                            metavar="TXT",
                            action=arg_formatter.CheckChromFile,
                            required=False)

    parser_grp.add_argument('-p', '--peak-file',
                            help='The file containing the peaks/regions to be annotated.'
                                 ' (bed format).',
                            default=None,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'),
                            required=True)

    parser_grp.add_argument('--more-bed',
                            help="A list of bed files to be considered as additional genomic annotations.",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'),
                            nargs='*',
                            required=False)

    parser_grp.add_argument('-l', '--more-bed-labels',
                            help="A comma separated list of labels (see --more-bed)",
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the TSS and TTS of in 5' by a given value.",
                            default=1000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-k', '--nb-threads',
                            help='Number of threads for multiprocessing.',
                            type=arg_formatter.ranged_num(0, None),
                            default=8,
                            required=False)

    parser_grp.add_argument('-s', '--seed',
                            help='Numpy random seed.',
                            type=arg_formatter.ranged_num(None, None),
                            default=42,
                            required=False)

    parser_grp.add_argument('-mn', '--minibatch-nb',
                            help='Number of minibatches of shuffles.',
                            type=arg_formatter.ranged_num(0, None),
                            default=10,
                            required=False)

    parser_grp.add_argument('-ms', '--minibatch-size',
                            help='Size of each minibatch, in number of shuffles.',
                            type=arg_formatter.ranged_num(0, None),
                            default=20,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the TSS and TTS of in  3' by a given value. ",
                            default=1000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-e', '--bed-excl',
                            help='Exclusion file. The chromosomes will be shortened by this much for the shuffles of peaks and features. Can take a long time.'
                                 ' (bed format).',
                            default=None,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'),
                            required=False)

    parser_grp.add_argument('-ma', '--use-markov',
                            help='Whether to use Markov shuffling instead of independant shuffles for respectively region lengths and inter-region lengths. Not recommended in the general case. Can take a *very* long time.',
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-pw', '--pdf-width',
                            help='Output pdf file width (inches).',
                            type=arg_formatter.ranged_num(0, None),
                            default=None,
                            required=False)

    parser_grp.add_argument('-ph', '--pdf-height',
                            help='Output pdf file height (inches).',
                            type=arg_formatter.ranged_num(0, None),
                            default=None,
                            required=False)

    parser_grp.add_argument('-m', '--more-keys',
                            help='A comma separated list of key used for labeling the genome. See Notes.',
                            type=str,
                            default=None,
                            required=False)

    parser_grp.add_argument('-n', '--no-basic-feature',
                            help="Don't compute statistics for genomic features but concentrates on --more-bed and --more-keys.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-if', '--user-img-file',
                            help="Provide an alternative path for the main image.",
                            default=None,
                            type=arg_formatter.FormattedFile(mode='w', file_ext='pdf'),
                            required=False)

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=arg_formatter.ranged_num(0, None),
                            default=300,
                            required=False)

    return parser


# -------------------------------------------------------------------------
# The command function
# -------------------------------------------------------------------------


def peak_anno(inputfile=None,
              outputdir=None,
              peak_file=None,
              chrom_info=None,

              more_bed=None,
              more_bed_labels=None,
              upstream=1000,
              more_keys=None,
              downstream=1000,
              no_basic_feature=False,
              bed_excl=None,
              use_markov=False,

              pdf_width=None,
              pdf_height=None,
              user_img_file=None,
              dpi=300,

              nb_threads=8,
              seed=42,
              minibatch_nb=8,
              minibatch_size=25,
              ):
    """
    This function is intended to perform statistics on peak intersection. It will compare your peaks to
    classical features (e.g promoter, tts, gene body, UTR,...) and to sets of user provided peaks.
    """

    # Set random seed
    np.random.seed(seed)

    # Load the peak file as pybedtools.BedTool object
    peak_file = pybedtools.BedTool(peak_file.name)

    # Just in case it was not, sort and merge the file.
    # In any case, it should be short compared to the expected total running time.
    peak_file = peak_file.sort().merge()

    # Are we using markov shuffling ?
    # If yes, send a warning to the user.
    if use_markov:
        message('Using Markov order 2 shuffling.', type='INFO')
        message(
            'Markov shuffling is still in beta at the moment and tends to biais the null hypothesis towards association.',
            type='WARNING')

    # -------------------------------------------------------------------------
    # If user wants no basic features (e.g prom, genes, exons) then he
    # needs to provide --more-keys or --more-bed
    # -------------------------------------------------------------------------

    if no_basic_feature:
        if more_keys is not None:
            if inputfile is None:
                message("If --more-keys is set you should provide a GTF",
                        type="ERROR")
        else:
            if more_bed is None:
                message("If --no-genomic-feature is set to True "
                        "provide --more-keys or --more-bed.",
                        type="ERROR")
    else:
        if inputfile is None:
            message("Please provide a GTF.",
                    type="ERROR")

        if chrom_info is None:
            message("Please provide a chromInfo file (--chrom-info)",
                    type="ERROR")

    # -------------------------------------------------------------------------
    # chrom_len will store the chromosome sizes.
    # -------------------------------------------------------------------------

    chrom_len = chrom_info_as_dict(chrom_info)

    # -------------------------------------------------------------------------
    # Region exclusion
    # -------------------------------------------------------------------------

    # If there is an exclusion of certain regions to be done, do it.
    # Here, we do exclusion on the peak file ('bedA') and the chrom sizes.
    # Exclusion on the other bed files or gtf extracted bed files ('bedB') is
    # done once we get to them.
    # overlap_stats_shuffling() will handle that, with the same condition : that bed_excl != None

    # WARNING : do not modify chrom_info or peak_file afterwards !
    # We can afford to modify chrom_len because we only modify its values, and the
    # rest of the peak_anno code relie otherwise only on its keys.

    if bed_excl is not None:
        # Treating bed_excl once and for all : turning it into a pybedtools file, merging it and sorting it.
        # NOTE This will prevent later conflicts, if two different pybedtools objects try to access it.
        bed_excl = pybedtools.BedTool(bed_excl)
        bed_excl = bed_excl.sort().merge()
        # Split in its constituent commands in case of very large files

        exclstart = time.time()
        message('Exclusion BED found, proceeding on the BED peaks file. This may take a few minutes.', type='INFO')

        chrom_len = read_bed.exclude_chromsizes(bed_excl, chrom_len)  # Shorten the chrom_len only once, and separately
        peak_file = read_bed.exclude_concatenate(pybedtools.BedTool(peak_file), bed_excl, nb_threads)

        exclstop = time.time()
        message('Exclusion completed for the BED PEAKS file in ' + str(exclstop - exclstart) + ' s', type='DEBUG')

    # -------------------------------------------------------------------------
    # Read the gtf file and discard any records corresponding to chr not declared
    # in ChromInfo file. This only needs to be done if one want basic feature
    # (default) or more-keys (e.g gene_biotype)
    # -------------------------------------------------------------------------

    if not no_basic_feature or more_keys:
        gtf = GTF(inputfile).select_by_key("seqid", ",".join(chrom_len.keys()))

        if len(gtf) == 0:
            message("The GTF file does not contain any genomic feature "
                    "falling in chromosomes declared in chromInfo file.",
                    type="ERROR")

        chrom_list = gtf.get_chroms(nr=True)

        # -------------------------------------------------------------------------
        # Check chromosomes are defined in the chrom-info file
        # -------------------------------------------------------------------------

        for i in chrom_list:
            if i not in chrom_len:
                message("Chromosome " + " i from GTF is undefined in --chrom-info file.",
                        type="ERROR")

    # -------------------------------------------------------------------------
    # Check user provided annotations
    # -------------------------------------------------------------------------

    if more_bed is not None:

        if more_bed_labels is not None:

            more_bed_labels = more_bed_labels.split(",")

            for elmt in more_bed_labels:
                if not re.search("^[A-Za-z0-9_]+$", elmt):
                    message(
                        "Only alphanumeric characters and '_' allowed for --more-bed-labels",
                        type="ERROR")
            if len(more_bed_labels) != len(more_bed):
                message("--more-bed-labels: the number of labels should be"
                        " the same as the number of bed files "
                        "("
                        "see --more-bed-labels).", type="ERROR")

            if len(more_bed_labels) != len(set(more_bed_labels)):
                message("Redundant labels not allowed.", type="ERROR")
        else:
            message(
                "--more-bed-labels should be set if --more-bed is used.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Preparing output files
    # -------------------------------------------------------------------------

    file_out_list = make_outdir_and_file(out_dir=outputdir,
                                         alist=["00_peak_anno_stats.tsv",
                                                "00_peak_anno_diagrams.pdf"
                                                ],
                                         force=True)

    data_file, pdf_file = file_out_list

    if user_img_file is not None:

        os.unlink(pdf_file.name)
        pdf_file = user_img_file

        test_path = os.path.abspath(pdf_file.name)
        test_path = os.path.dirname(test_path)

        if not os.path.exists(test_path):
            os.makedirs(test_path)

    # -------------------------------------------------------------------------
    # Check chromosomes for peaks are defined in the chrom-info file
    # -------------------------------------------------------------------------

    chrom_list = set()
    for i in pybedtools.BedTool(peak_file):
        chrom_list.add(i.chrom)

    for i in chrom_list:
        if i not in chrom_len:
            message("Chromosome " + " i from GTF is undefined in --chrom-info file.",
                    type="ERROR")

    # -------------------------------------------------------------------------
    # Fill the dict with info about basic features include in GTF
    # -------------------------------------------------------------------------

    # Prepare a partial call with all fixed parameters (ie. everything except
    # the two bed files) for code legibility.
    overlap_partial = partial(compute_overlap_stats, chrom_len=chrom_len,
                              minibatch_size=minibatch_size, minibatch_nb=minibatch_nb,
                              bed_excl=bed_excl, use_markov_shuffling=use_markov,
                              nb_threads=nb_threads)

    # Initialize result dict
    hits = dict()

    if not no_basic_feature:
        for feat_type in gtf.get_feature_list(nr=True):
            message("Processing " + str(feat_type), type="INFO")
            gtf_sub = gtf.select_by_key("feature", feat_type, 0)
            gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                               "gene_id",
                                               "exon_id"]).sort().merge()  # merging bed file !

            del gtf_sub

            hits[feat_type] = overlap_partial(bedA=peak_file, bedB=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the intergenic regions
        # -------------------------------------------------------------------------

        message("Processing intergenic regions", type="INFO")
        gtf_sub_bed = gtf.get_intergenic(chrom_info,
                                         0,
                                         0,
                                         chrom_len.keys()).merge()

        hits["Intergenic"] = overlap_partial(bedA=peak_file, bedB=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the intronic regions
        # -------------------------------------------------------------------------

        message("Processing on : Introns", type="INFO")
        gtf_sub_bed = gtf.get_introns()

        hits["Introns"] = overlap_partial(bedA=peak_file, bedB=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the promoter regions
        # -------------------------------------------------------------------------

        message("Processing promoters", type="INFO")
        gtf_sub_bed = gtf.get_tss().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits["Promoters"] = overlap_partial(bedA=peak_file, bedB=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the tts regions
        # -------------------------------------------------------------------------

        message("Processing terminator", type="INFO")
        gtf_sub_bed = gtf.get_tts().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits["Terminator"] = overlap_partial(bedA=peak_file, bedB=gtf_sub_bed)

    # -------------------------------------------------------------------------
    # if the user request --more-keys (e.g. gene_biotype)
    # -------------------------------------------------------------------------

    if more_keys is not None:

        more_keys_list = more_keys.split(",")

        if len(more_keys_list) > 50:
            message("The selected key in --more-keys should be "
                    "associated with less than 50 different values.",
                    type="ERROR")

        for user_key in more_keys_list:
            user_key_values = set(gtf.extract_data(user_key,
                                                   as_list=True,
                                                   hide_undef=True,
                                                   no_na=True,
                                                   nr=True))

            if len(user_key_values) > 50:
                message("The selected key in --more-keys "
                        "should be associated with less than 50 different values.",
                        type="ERROR")

            for val in user_key_values:

                gtf_sub = gtf.select_by_key(user_key, val, 0)

                if len(gtf_sub) > 0:
                    gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                                       "gene_id",
                                                       "exon_id"]).sort().merge()  # merging bed file !
                    del gtf_sub
                    ft_type = ":".join([user_key, val])  # Key for the dictionary
                    hits[ft_type] = overlap_partial(bedA=peak_file,
                                                    bedB=gtf_sub_bed)
                    message("Processing " + str(ft_type), type="INFO")

    # -------------------------------------------------------------------------
    # Process user defined annotations
    # -------------------------------------------------------------------------

    if more_bed is not None:
        message("Processing user-defined regions (bed format).")
        for bed_anno, bed_lab in zip(more_bed, more_bed_labels):
            message("Processing " + str(bed_lab), type="INFO")
            chrom_list = set()
            for i in BedTool(bed_anno.name):
                chrom_list.add(i.chrom)

            for i in chrom_list:
                if i not in chrom_len:
                    message("Chromosome " + " i from GTF is undefined in " + bed_anno.name + " file.",
                            type="ERROR")

            hits[bed_lab] = overlap_partial(bedA=peak_file,
                                            bedB=BedTool(bed_anno.name))

    # ------------------ Treating the 'hits' dictionary --------------------- #

    if len(hits) == 0:
        message("No feature found.", type="ERROR")

    ### Print the 'hits' dictionary into a tabulated file

    should_print_header = True

    for feature_type in hits.keys():

        current_dict = hits[feature_type]  # This is an ordered dict

        # First line should be a header
        if should_print_header:
            header = [str(s) for s in hits[feature_type].keys()]

            data_file.write("\t".join(['feature_type'] + header) + "\n")
            should_print_header = False

        values = []
        for k, v in current_dict.items():
            values = values + [str(v)]

        data_file.write("\t".join([feature_type] + values) + "\n")

    close_properly(data_file)

    # -------------------------------------------------------------------------
    # Read the data set and plot it
    # -------------------------------------------------------------------------

    d = pd.read_csv(data_file.name, sep="\t", header=0)

    plot_results(d, data_file, pdf_file, pdf_width, pdf_height, dpi)


def plot_results(d, data_file, pdf_file, pdf_width, pdf_height, dpi):
    """
    Main plotting function by Q. Ferré and D. Puthier
    """

    if d.shape[0] == 0:
        message("No lines found in input file.",
                type="ERROR")

    # Save the data file
    d.to_csv(open(data_file.name, 'w'), sep="\t", header=True, index=False)

    # -------------------------------------------------------------------------
    # Copy the data
    # -------------------------------------------------------------------------
    dm = d.copy()

    message('Adding bar plot.')

    # -------------------------------------------------------------------------
    # Create a new plot
    # -------------------------------------------------------------------------

    def plot_this(statname):

        # -------------- First plot : number of intersections ---------------- #

        # Collect true and shuffled number of intersections
        data_ni = dm[['feature_type', statname + '_esperance_shuffled', statname + '_true']]
        maximum = data_ni[[statname + '_esperance_shuffled', statname + '_true']].max(axis=1)

        data_ni.columns = ['Feature', 'Shuffled', 'True']  # Rename columns
        dmm = data_ni.melt(id_vars='Feature')
        dmm.columns = ['Feature', 'Type', statname]

        # Create plot
        p = ggplot(dmm)

        # Bar plot of shuffled vs true
        aes_plot = aes(x='Feature', y=statname, fill='Type')
        p += geom_bar(mapping=aes_plot, stat='identity', alpha=0.6, position='dodge', show_legend=True, width=.6)

        # Add error bars for the standard deviation of the shuffles
        errorbar_mins = dm[statname + '_esperance_shuffled'] - np.sqrt(dm[statname + '_variance_shuffled'])
        errorbar_maxs = dm[statname + '_esperance_shuffled'] + np.sqrt(dm[statname + '_variance_shuffled'])

        # True values have no error
        na_series = pd.Series([np.nan] * len(errorbar_mins))
        errorbar_mins = errorbar_mins.append(na_series)
        errorbar_mins.index = range(len(errorbar_mins))
        errorbar_maxs = errorbar_maxs.append(na_series)
        errorbar_maxs.index = range(len(errorbar_maxs))

        p += geom_errorbar(aes(x='Feature', ymin=errorbar_mins, ymax=errorbar_maxs, fill='Type'), width=.5,
                           position=position_dodge(.6))

        # Text for the p-value
        text = dm[statname + '_pvalue'].append(na_series)
        text.index = range(len(text))

        # Format the text
        def format_pvalue(x):
            if x < 1.12E-16:
                r = 'p<1.12E-16'  # If the p-value is under 1.12E-16 (log precision limit), say so
            else:
                r = 'p=' + '{0:.3g}'.format(x)  # Add 'p=' before and format the p value
            return r

        text = text.apply(format_pvalue)
        text_pos = (maximum + 0.05 * max(maximum)).append(na_series)
        text_pos.index = range(len(text_pos))
        aes_plot = aes(x='Feature', y=text_pos, label=text, fill='Type')
        p += geom_text(mapping=aes_plot, stat='identity', size=5)

        # Theme
        p += theme(legend_title=element_blank(),
                   legend_position="top",
                   legend_box_spacing=0.65,
                   legend_key_size=8,
                   legend_text=element_text(size=8),
                   legend_key=element_blank(),
                   axis_title_x=element_blank(),
                   axis_title_y=element_text(colour='#333333',
                                             size=8,
                                             hjust=4,
                                             angle=90,
                                             face="plain"),
                   axis_text_y=element_text(size=5,
                                            margin={'r': 0},
                                            angle=0),
                   axis_text_x=element_text(size=5,
                                            angle=45)
                   )

        return p

    # Compute the plots for both statistics
    p1 = plot_this('nb_intersections') + ylab("Number of intersections")
    p2 = plot_this('summed_bp_overlaps') + ylab("Nb. of overlapping base pairs")

    # -------------------------------------------------------------------------
    # Computing page size
    # -------------------------------------------------------------------------

    nb_ft = len(list(d['feature_type'].unique()))

    if pdf_width is None:
        panel_width = 0.5
        pdf_width = panel_width * nb_ft

        if panel_width > 25:
            panel_width = 25
            message("Setting --pdf-width to 25 (limit)")

    if pdf_height is None:
        pdf_height = 5

    message("Page width set to " + str(pdf_width))
    message("Page height set to " + str(pdf_height))

    # -------------------------------------------------------------------------
    # Turn warning off. Both pandas and plotnine use warnings for deprecated
    # functions. I need to turn they off although I'm not really satisfied with
    # this solution...
    # -------------------------------------------------------------------------

    def fxn():
        warnings.warn("deprecated", DeprecationWarning)

    # -------------------------------------------------------------------------
    # Saving
    # -------------------------------------------------------------------------

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()
        message("Saving diagram to file : " + pdf_file.name)
        message("Be patient. This may be long for large datasets.")

        save_as_pdf_pages(filename=pdf_file.name,
                          plots=[p1, p2],
                          width=pdf_width,
                          height=pdf_height,
                          dpi=dpi)

    close_properly(pdf_file, data_file)


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    peak_anno(**args)


if __name__ == '__main__':
    main()


else:

    # 'Bats' tests
    test = '''
        #peak_anno: get example files
        @test "peak_anno_0" {
             result=`gtftk get_example -d simple_02 -f '*'`
          [ "$result" = "" ]
        }

        #peak_anno: run on simple test file
        @test "peak_anno_1" {
             result=`rm -Rf peak_annotation; gtftk peak_anno -i simple_02.gtf -p simple_02_peaks.bed -c simple_02.chromInfo -u 2 -d 2 -K peak_annotation`
          [ "$result" = "" ]
        }

        #peak_anno: proper number of true intersections
        @test "peak_anno_2" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep gene | cut -f 5`
          [ "$result" = "16" ]
        }

        #peak_anno: proper number of shuffled intersections
        @test "peak_anno_3" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep gene | cut -f 2`
          [ "$result" = "14.97" ]
        }

        #peak_anno: overlapping bp
        @test "peak_anno_4" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep gene | cut -f 11`
          [ "$result" = "75" ]
        }

        #peak_anno: shuffled overlapping bp
        @test "peak_anno_5" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep gene | cut -f 8`
          [ "$result" = "61.35" ]
        }

        #peak_anno: shuffled overlapping bp variance
        @test "peak_anno_6" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep gene | cut -f 9`
          [ "$result" = "32.94" ]
        }

        #peak_anno: shuffled overlapping bp fitting
        @test "peak_anno_7" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep gene | cut -f 10`
          [ "$result" = "0.8227700000000001" ]
        }
        '''

    cmd = CmdObject(name="peak_anno",
                    message="Statistics on bed file intersections with genomic features.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="annotation",
                    notes=__notes__,
                    updated=__updated__,
                    test=test)
