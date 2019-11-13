#!/usr/bin/env python

import argparse
import math
import os
import re
import sys
import warnings
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd
from plotnine import (ggplot, aes, ggsave, ylab, theme, element_blank, element_text, facet_wrap, geom_tile,
                      xlab, scale_fill_gradient2)
from scipy.stats import binom_test

from pygtftk.arg_formatter import CheckChromFile
from pygtftk.arg_formatter import FormattedFile
from pygtftk.arg_formatter import ranged_num
from pygtftk.bedtool_extension import BedTool
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chrom_info_as_dict
from pygtftk.utils import close_properly
from pygtftk.utils import intervals
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

R_LIBS = "reshape2,ggplot2"

__updated__ = "2018-01-24"
__doc__ = """
 Annotate peaks (in bed format) with region sets computed on the
 fly from a GTF file  (e.g promoter, tts, gene body, UTR...). The midpoint of
 each peak is considered and intersected iteratively with region sets. A binomial
 p-value is computed based on hypothesized probability of success p (fraction of genome covered by the
 feature f), the number of trials (number of peaks) and the number of successes (number of intersections).
 """

__notes__ = """
 -- Genome size is computed from the provided chromInfo file (-c). It should thus only contain ordinary chromosomes.
 
 -- The program produces two pdf files and one txt file ('_stats_') containing intersection statistics.
 The two pdf files correspond to the statistics performed on the whole genome or at the level of the chromosomes.
 In the case of the chromosomes ('_by_chrom_' pdf file) the question is to test whether enrichments/depletions observed
 at a global level are also observed throughout chromosomes or whether some of them deviate from the general trend.
 
 -- If -\-more-keys is used additional region sets will be tested based on the associated key value.
 As an example, if -\-more-keys is set to the 'gene_biotype' (a key generally found in ensembl GTF), the
 region related to 'protein_coding', 'lncRNA' or any other values for that key will be retrieved merged and tested
 for enrichment.
 -- Use -\no-basic-feature if you want to perform enrichment analysis on focused annotations only (-\-more-bed or -\-more-key).
 -- TODO: This function does not support a mappability file at the moment...
 -- TODO: the png output by chromosomes is not functional at the moment. 
 """


def make_parser():
    """The main parser."""
    # parser = OptionParser()

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputdir',
                            help='Output directory name.',
                            metavar="DIR",
                            default="peak_annotation",
                            type=str)

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Tabulated two-columns file. Chromosomes as column 1, sizes as column 2",
                            default=None,
                            metavar="TXT",
                            action=CheckChromFile,
                            required=False)

    parser_grp.add_argument('-p', '--peak-file',
                            help='The file containing the peaks/regions to be annotated (bed format).',
                            default=None,
                            metavar="BED",
                            required=True,
                            type=FormattedFile(mode='r', file_ext=('bed')))

    parser_grp.add_argument('-b', '--more-bed',
                            help="A list of bed files to be considered as additional genomic annotations.",
                            type=FormattedFile(mode='r', file_ext='bed'),
                            nargs='*',
                            required=False)

    parser_grp.add_argument('-l', '--more-bed-labels',
                            help="A comma separated list of labels."
                                 "(see --more-bed)",
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the TSS and TTS of in 5' by a given value.",
                            default=1000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the TSS and TTS of in  3' by a given value. ",
                            default=1000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-q', '--nb-bin',
                            help="How many bins for the reference.",
                            default=10,
                            type=int,
                            required=False)

    parser_grp.add_argument('-e', '--nb-bin-peak',
                            help="How many bins for the peaks.",
                            default=10,
                            type=int,
                            required=False)

    parser_grp.add_argument('-pw', '--pdf-width',
                            help='Output pdf file width (inches).',
                            type=ranged_num(lowest=0, highest=None,
                                            val_type="float", linc=False),
                            default=7,
                            required=False)

    parser_grp.add_argument('-ph', '--pdf-height',
                            help='Output pdf file height (inches).',
                            type=ranged_num(lowest=0, highest=None,
                                            val_type="float", linc=False),
                            default=7,
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
                            type=argparse.FileType("w"),
                            required=False)

    parser_grp.add_argument('-pf', '--page-format',
                            help='Output file format.',
                            choices=['pdf', 'png'],
                            default='pdf',
                            required=False)

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=ranged_num(lowest=0, highest=None,
                                            val_type="float", linc=False),
                            default=300,
                            required=False)

    parser_grp.add_argument('-r', '--order-bar',
                            help='The way bar should be ordered.',
                            choices=['log2_ratio', 'pval_binom'],
                            default='log2_ratio',
                            required=False)

    parser_grp.add_argument('-z', '--no-gtf',
                            help="No GTF file is provided as input.",
                            action='store_true',
                            required=False)

    return parser


# -------------------------------------------------------------------------
# This function performs intersection between a bedTools object (peak_file/query)
# and another file (feature_file/reference) e.g promoter, intron, exon...
# Fill a dictionary of dictionary
# -------------------------------------------------------------------------


def _intersection_results(peak_file=None,
                          feature_bo=None,
                          my_dict=None,
                          ft_type=None,
                          qt=0,
                          qt_ref=1):
    """Get feature and peaks and compute intersections. Returns an updated
    dict.
    """

    feature_bo = BedTool(feature_bo)

    if len(feature_bo) == 0:
        return my_dict

    # -------------------------------------------------------------------------
    # Compute the number of regions of peaks
    # -------------------------------------------------------------------------

    peak_file_bo = BedTool(peak_file)
    nb_peaks = len(peak_file_bo)
    my_dict[ft_type][qt][qt_ref]["nb_peaks"] = nb_peaks

    # -------------------------------------------------------------------------
    # Compute the nucleotide size of all feature (promoter, exons, introns,...)
    # Save features (introns, intergenic,...) for tracability
    # -------------------------------------------------------------------------

    feature_bo = feature_bo.sort()
    feature_merge_bo = feature_bo.merge()

    peak_anno_tmp_file = make_tmp_file(prefix="peak_anno_" + ft_type + "_merged",
                                       suffix=".bed")

    feature_merge_bo.saveas(peak_anno_tmp_file.name)
    peak_anno_tmp_file.close()

    for i in feature_merge_bo:
        size = i.end - i.start
        my_dict[ft_type][qt][qt_ref]["coverage"] += size

    # -------------------------------------------------------------------------
    # Compute intersections
    # -------------------------------------------------------------------------

    intersections = peak_file_bo.intersect(feature_merge_bo, u=True)
    my_dict[ft_type][qt][qt_ref]["Observed"] = len(intersections)

    file_out_save = make_tmp_file(prefix="peak_anno_intersections_" +
                                         ft_type,
                                  suffix=".bed")

    intersections.saveas(file_out_save.name)

    return my_dict


# -------------------------------------------------------------------------
# The command function
# -------------------------------------------------------------------------


def pogos(inputfile=None,
          outputdir=None,
          peak_file=None,
          more_bed=None,
          more_bed_labels=None,
          upstream=1000,
          more_keys=None,
          downstream=1000,
          no_basic_feature=False,
          pdf_width=None,
          pdf_height=None,
          no_gtf=False,
          chrom_info=None,
          user_img_file=None,
          page_format=None,
          dpi=300,
          order_bar=None,
          nb_bin_peak=10,
          nb_bin=None):
    """
    This function is intended to perform statistics on peak intersection. It will compare your peaks to
    classical features (e.g promoter, tts, gene body, UTR,...) and to sets of user provided peaks.
    """

    # -------------------------------------------------------------------------
    # If the user wishes not to provide a GTF
    # -------------------------------------------------------------------------

    if no_gtf:
        if more_keys is not None:
            message("If --more-keys should be used with a GTF.",
                    type="ERROR")
        if no_basic_feature:
            message("If --no-basic-feature should be used with a GTF.",
                    type="ERROR")
        if more_bed is None:
            message("If --no-gtf is set to True provide --more-bed.",
                    type="ERROR")

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
    # The hits variable is a dict to store the results per ft_type
    #   - e.g my_dict[ft_type]["Observed"]
    #   - e.g my_dict[ft_type]["coverage"]
    #   - e.g my_dict[ft_type]["nb_peakss"]
    # -------------------------------------------------------------------------

    def nested_dict(n, type):
        """"http://stackoverflow.com/questions/29348345"""
        if n == 1:
            return defaultdict(type)
        else:
            return defaultdict(lambda: nested_dict(n - 1, type))

    hits = nested_dict(4, int)

    # -------------------------------------------------------------------------
    # Read the gtf file and discard any records corresponding to chr not declared
    # in ChromInfo file. This only needs to be done if one want basic feature
    # (default) or more-keys (e.g gene_biotype)
    # -------------------------------------------------------------------------

    if not no_gtf:
        if not no_basic_feature or more_keys:

            gtf = GTF(inputfile,
                      check_ensembl_format=False
                      ).select_by_key("seqid", ",".join(chrom_len.keys()))

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
                        "( see --bedAnnotationList).", type="ERROR")

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
                                         alist=["00_peak_anno_stats.txt",
                                                "00_peak_anno_diagrams." + page_format],
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

    for i in BedTool(peak_file):
        chrom_list.add(i.chrom)

    for i in chrom_list:
        if i not in chrom_len:
            message("Chromosome " + " i from BED is undefined in --chrom-info file.",
                    type="ERROR")

    # -------------------------------------------------------------------------
    # Get the quantiles of the peak region positions.
    # Store them in a list
    # -------------------------------------------------------------------------

    peak_bin_bo_list = []

    quantiles_list = np.linspace(0, 1, nb_bin_peak)

    for qt in quantiles_list:
        peak_bin_file = make_tmp_file("peaks_bin_" + str(round(qt, 10)), ".bed")

        # Loop through peaks

        peak_bin_bo = BedTool(peak_file.name).get_quantile_pos(quantile=qt)
        peak_bin_bo.saveas(peak_bin_file.name)
        peak_bin_bo_list += [peak_bin_bo]

    # -------------------------------------------------------------------------
    # Fill the dict with info about basic features include in GTF
    # -------------------------------------------------------------------------

    def binned_files_from_bed(input_bed, nb_bin, ft_type):

        input_bed = BedTool(input_bed)

        binned_files_bed_list = list()

        for bin in range(nb_bin):
            binned_files_bed_list += [make_tmp_file("reference_bin_" + ft_type + "_" + str(bin), ".bed")]

        for line in input_bed:

            my_intervals = intervals(range(line.start, line.end), nb_bin, silent=True)

            if my_intervals is None:
                continue

            if line.strand == "-":
                my_intervals = list(reversed(my_intervals))

            for pos, this_intervals in enumerate(my_intervals):
                binned_files_bed_list[pos].write("\t".join([str(line.chrom),
                                                            str(this_intervals[0]),
                                                            str(this_intervals[1])]) + "\n")
        for bed in binned_files_bed_list:
            bed.close()

        for pos, bed in enumerate(binned_files_bed_list):

            bed_bo = BedTool(bed.name)
            if len(bed_bo) > 0:
                binned_files_bed_list[pos] = bed_bo

        return binned_files_bed_list

    if not no_gtf:

        if not no_basic_feature:
            for feat_type in gtf.get_feature_list(nr=True):
                if feat_type not in ["start_codon", "stop_codon"]:

                    gtf_sub = gtf.select_by_key("feature", feat_type, 0)

                    gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                                       "gene_id",
                                                       "exon_id"]).sort().merge_by_strand()  # merging bed file !

                    a_bed_file_list = binned_files_from_bed(gtf_sub_bed, nb_bin, ft_type=feat_type)
                    for nb, a_bed_file in enumerate(a_bed_file_list):
                        for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                            hits = _intersection_results(peak_file=peak_bo.fn,
                                                         feature_bo=a_bed_file,
                                                         my_dict=hits,
                                                         ft_type=feat_type,
                                                         qt=round(qt, 10),
                                                         qt_ref=nb)

            # -------------------------------------------------------------------------
            # Fill the dict with info about  intergenic regions
            # -------------------------------------------------------------------------

            gtf_sub_bed = gtf.get_intergenic(chrom_info,
                                             0,
                                             0,
                                             chrom_len.keys()).merge()

            a_bed_file_list = binned_files_from_bed(gtf_sub_bed, nb_bin, ft_type="intergenic")

            for nb, a_bed_file in enumerate(a_bed_file_list):
                for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                    hits = _intersection_results(peak_file=peak_bo.fn,
                                                 feature_bo=a_bed_file,
                                                 my_dict=hits,
                                                 ft_type="Intergenic",
                                                 qt=round(qt, 10),
                                                 qt_ref=nb)

            # -------------------------------------------------------------------------
            # Fill the dict with info about  the intronic regions
            # -------------------------------------------------------------------------

            gtf_sub_bed = gtf.get_introns()

            a_bed_file_list = binned_files_from_bed(gtf_sub_bed, nb_bin, ft_type="intronic")

            for nb, a_bed_file in enumerate(a_bed_file_list):

                for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                    hits = _intersection_results(peak_file=peak_bo.fn,
                                                 feature_bo=a_bed_file,
                                                 my_dict=hits,
                                                 ft_type="Introns",
                                                 qt=round(qt, 10),
                                                 qt_ref=nb)

            # -------------------------------------------------------------------------
            # Fill the dict with info about promoter regions
            # -------------------------------------------------------------------------

            gtf_sub_bed = gtf.get_tss().slop(s=True,
                                             l=upstream,
                                             r=downstream,
                                             g=chrom_info.name).cut([0, 1, 2,
                                                                     3, 4, 5]).sort().merge_by_strand()

            a_bed_file_list = binned_files_from_bed(gtf_sub_bed, nb_bin, ft_type="promoter")

            for nb, a_bed_file in enumerate(a_bed_file_list):

                for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                    hits = _intersection_results(peak_file=peak_bo.fn,
                                                 feature_bo=a_bed_file,
                                                 my_dict=hits,
                                                 ft_type="Promoters",
                                                 qt=round(qt, 10),
                                                 qt_ref=nb)

            # -------------------------------------------------------------------------
            # Fill the dict with info about tts regions
            # -------------------------------------------------------------------------

            gtf_sub_bed = gtf.get_tts().slop(s=True,
                                             l=upstream,
                                             r=downstream,
                                             g=chrom_info.name).cut([0, 1, 2,
                                                                     3, 4, 5]).sort().merge_by_strand()

            a_bed_file_list = binned_files_from_bed(gtf_sub_bed, nb_bin, ft_type="tts")

            for nb, a_bed_file in enumerate(a_bed_file_list):

                for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                    hits = _intersection_results(peak_file=peak_bo.fn,
                                                 feature_bo=a_bed_file,
                                                 my_dict=hits,
                                                 ft_type="Terminator",
                                                 qt=round(qt, 10),
                                                 qt_ref=nb)

    # -------------------------------------------------------------------------
    # if the user request --more-keys (e.g. gene_biotype)
    # Fill the dict with corresponding info
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
                                                       "exon_id"]).sort().merge_by_strand()  # merging bed file !

                    a_bed_file_list = binned_files_from_bed(gtf_sub_bed, nb_bin, ft_type=val)

                    for nb, a_bed_file in enumerate(a_bed_file_list):

                        for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                            hits = _intersection_results(peak_file=peak_bo.fn,
                                                         feature_bo=a_bed_file,
                                                         my_dict=hits,
                                                         ft_type=":".join([user_key,
                                                                           val]),
                                                         qt=round(qt, 10),
                                                         qt_ref=nb)

    # -------------------------------------------------------------------------
    # Process user defined annotations
    # -------------------------------------------------------------------------

    if more_bed is not None:
        message("Processing user-defined regions (bed format).")
        for bed_anno, bed_lab in zip(more_bed, more_bed_labels):

            chrom_list = set()
            for i in BedTool(bed_anno.name):
                chrom_list.add(i.chrom)

            for i in chrom_list:
                if i not in chrom_len:
                    message("Chromosome " + " i from GTF is undefined in " + bed_anno.name + " file.",
                            type="ERROR")

            a_bed_file_list = binned_files_from_bed(bed_anno, nb_bin, ft_type=bed_lab)

            for nb, a_bed_file in enumerate(a_bed_file_list):

                for peak_bo, qt in zip(peak_bin_bo_list, quantiles_list):
                    hits = _intersection_results(peak_file=peak_bo.fn,
                                                 feature_bo=a_bed_file,
                                                 my_dict=hits,
                                                 ft_type=bed_lab,
                                                 qt=round(qt, 10),
                                                 qt_ref=nb)

    # Store the result into a file
    # before loading it into R
    message("Storing data in: " + data_file.name)

    # Write header
    data_file.write("\t".join(["ft_type",
                               "quantile",
                               "quantile_ref",
                               "reference_size",
                               "nb_peaks",
                               "coverage",
                               "Observed\n"]))

    if len(hits) == 0:
        message("No feature found.", type="ERROR")

    print(hits)
    # Write data
    for ft_type in hits:
        for bin in hits[ft_type]:
            for bin_ref in hits[ft_type][bin]:
                out_list = [ft_type,
                            str(bin),
                            str(bin_ref),
                            str(chrom_len['all_chrom']),  # genome size
                            str(hits[ft_type][bin][bin_ref]["nb_peaks"]),
                            str(hits[ft_type][bin][bin_ref]["coverage"]),
                            str(hits[ft_type][bin][bin_ref]["Observed"])]

                data_file.write("\t".join(out_list) + "\n")

    close_properly(data_file)

    # -------------------------------------------------------------------------
    # Read the data set
    # -------------------------------------------------------------------------
    message("Reading overlapping statistics")

    d = pd.read_csv(data_file.name, sep="\t", header=0)

    if d.shape[0] == 0:
        message("No lines found in input file.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Compute expected number of intersections
    # -------------------------------------------------------------------------
    message("Genome coverage")

    d['freq'] = d['coverage'] / d['reference_size']
    d['Expected'] = d['freq'] * d['nb_peaks']

    # -------------------------------------------------------------------------
    # Compute binomial p.val  (unilateral)
    # -------------------------------------------------------------------------

    message("Computing binomial test")

    for i, _ in d.iterrows():

        message("Computing log ratio")

        if d.loc[i, 'Expected'] > 0:
            if d.loc[i, 'Observed'] > 0:
                log2_ratio = np.log2(d.loc[i, 'Observed'] / d.loc[i, 'Expected'])
                d.loc[i, 'log2_ratio_str'] = "{0:0.2f}".format(log2_ratio)
                d.loc[i, 'log2_ratio'] = log2_ratio
            else:
                d.loc[i, 'log2_ratio_str'] = "NA"
                d.loc[i, 'log2_ratio'] = float('nan')
        else:
            d.loc[i, 'log2_ratio'] = float('nan')
            d.loc[i, 'log2_ratio_str'] = "NA"

        message("Computing binomial test.")

        if d.loc[i, 'nb_peaks'] > 0:

            pval = binom_test(n=int(d.loc[i, 'nb_peaks']),
                              x=int(d.loc[i, 'Observed']),
                              p=float(d.loc[i, 'freq']))

            d.loc[i, 'pval_binom'] = pval
            d.loc[i, 'pval_binom_str'] = "{0:0.3g}".format(pval)
            if pval == 0:
                pval = 1e-320
            d.loc[i, 'pval_binom_log10'] = - math.log10(pval)
        else:
            d.loc[i, 'pval_binom_str'] = "NA"
            d.loc[i, 'pval_binom_log10'] = np.nan

        message("Checking Enrichment.")

        d.loc[i, 'test_type'] = 'Unchanged'

        if d.loc[i, 'Observed'] > d.loc[i, 'Expected']:
            if d.loc[i, 'pval_binom'] < 0.05:
                d.loc[i, 'test_type'] = 'Enrichment'
        elif d.loc[i, 'Observed'] < d.loc[i, 'Expected']:
            if d.loc[i, 'pval_binom'] < 0.05:
                d.loc[i, 'test_type'] = 'Depletion'

    # -------------------------------------------------------------------------
    # ft_type for display
    # -------------------------------------------------------------------------

    for i, _ in d.iterrows():
        d.loc[i, 'ft_type'] = d.loc[i, 'ft_type'].replace(":", ":\n")

    # -------------------------------------------------------------------------
    # Save file
    # -------------------------------------------------------------------------

    d.to_csv(open(data_file.name, 'w'), sep="\t", header=True, index=False)

    # -------------------------------------------------------------------------
    # Prepare a function for drawing diagram
    # -------------------------------------------------------------------------

    def heatmap(d=None,
                x_ft_type=None):

        # -------------------------------------------------------------------------
        # Compute text position
        # -------------------------------------------------------------------------

        max_y = max(d['Expected'].tolist() + d['Observed'].tolist())

        offset = 10 / 100 * max_y

        for i, _ in d.iterrows():
            ## max_y to be use to plot pval
            max_y_col = max(d.loc[i, 'Observed'],
                            d.loc[i, 'Expected'])

            d.loc[i, 'y_log2_ratio'] = max_y_col + offset * 2.4
            d.loc[i, 'y_pval'] = max_y_col + offset * 1.8
            d.loc[i, 'y_text'] = max_y_col + offset / 10

            ## y_lim to set diagram limits in y coords.
            d.loc[i, 'y_lim'] = max_y + offset * 2.2

        # -------------------------------------------------------------------------
        # Melt the data frame
        # -------------------------------------------------------------------------

        message('Melting.')

        dm = d.melt(id_vars=[x for x in d.columns if x not in ['Expected',
                                                               'Observed']],
                    value_vars=['Observed', 'Expected'])

        dm['variable'] = pd.Categorical(dm['variable'])

        # -------------------------------------------------------------------------
        # Create a new plot
        # -------------------------------------------------------------------------

        p = ggplot()

        # -------------------------------------------------------------------------
        # Create facet
        # -------------------------------------------------------------------------

        if nb_bin_peak > 1:
            p += facet_wrap("~ ft_type")

        # -------------------------------------------------------------------------
        # Order features levels (Categories) based on binom test
        # -------------------------------------------------------------------------

        dm.loc[:, x_ft_type] = pd.Categorical(dm[x_ft_type].tolist())

        levels_ordered = [x for _, x in sorted(zip(dm[order_bar].tolist(),
                                                   dm[x_ft_type].tolist()),
                                               key=lambda x: x[0] if not math.isnan(x[0]) else 0,
                                               reverse=True)]
        unique = list(OrderedDict.fromkeys(levels_ordered))

        dm[x_ft_type].cat.reorder_categories(unique, inplace=True)

        # -------------------------------------------------------------------------
        # Display bars
        # -------------------------------------------------------------------------

        message('Adding bar plot.')

        if nb_bin_peak > 1:
            aes_plot = aes(x='quantile',
                           y='quantile_ref',
                           fill='log2_ratio')
        else:
            aes_plot = aes(x='ft_type',
                           y='quantile_ref',
                           fill='log2_ratio')

        p += geom_tile(data=dm,
                       mapping=aes_plot,
                       show_legend=True)

        p += ylab("5' <- Reference -> 3'")
        p += xlab("Peaks")

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

        p += scale_fill_gradient2(low="#005824",
                                  mid="#ffffff",
                                  high="#990000",
                                  midpoint=0)

        return (p)

    # -------------------------------------------------------------------------
    # Compute bar plot by feature
    # -------------------------------------------------------------------------

    message("Computing barplot.")

    p = heatmap(d=d,
                x_ft_type='ft_type')

    # -------------------------------------------------------------------------
    # Turn warning off. Both pandas and plotnine use warnings for deprecated
    # functions. I need to turn they off although I'm not really satisfied with
    # this solution...
    # -------------------------------------------------------------------------

    def fxn():
        warnings.warn("deprecated", DeprecationWarning)

    # -------------------------------------------------------------------------
    #
    # Saving
    #
    # -------------------------------------------------------------------------

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()
        message("Saving diagram to file : " + pdf_file.name)
        message("Be patient. This may be long for large datasets.")
        ggsave(filename=pdf_file.name,
               plot=p,
               width=pdf_width,
               height=pdf_height,
               dpi=dpi)
        # dm.to_csv(data_file, sep="\t", header=True, index=False)

    close_properly(pdf_file, data_file)


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    pogos(**args)


if __name__ == '__main__':

    main()


else:

    test = '''
            
        #pogos: chr2 len
        @test "pogos1" {
             result=`rm -Rf pogos; gtftk pogos -i pygtftk/data/simple_02/simple_02.gtf -p pygtftk/data/simple_02/simple_02_peaks.bed -c pygtftk/data/simple_02/simple_02.chromInfo -u 2 -d 2 -K peak_annotation`
          [ "$result" = "" ]
        }
        
        
        #pogos: all_chrom len
        @test "pogos2" {
         result=`cat pogos/00_pogosstats_* | grep chr2 | cut -f 3 | sort | uniq | perl -npe 's/\\n/,/'`
          [ "$result" = "1700,400," ]
        }
        
        
        #pogos: all_chrom len
        @test "pogos3" {
         result=`cat pogos/00_pogosstats_* | grep all_chrom | cut -f 3 | sort | uniq`
          [ "$result" = "1700" ]
        }
        
        #pogos: there are 11 guys (midpoints) intersecting CDS
        @test "pogos4" {
         result=`grep -w CDS pogos/00_pogosstats_* | cut -f 6 | perl -npe 's/\\n/,/'`
          [ "$result" = "0,0,0,11,11,0," ]
        }
        
        
        #pogos: there is 59 guys (midpoints) intersecting intergenic regions (and 19 on chr1).
        @test "pogos5" {
         result=`grep -i Interg pogos/00_pogosstats_* | cut -f 6 | perl -npe 's/\\n/,/'`
          [ "$result" = "0,0,40,19,59,0," ]
        }
        
        
        
        #pogos: size of intronic features is 21 bases
        @test "pogos7" {
         result=`grep -i Introns pogos/00_pogos_stats_* | cut -f 5 | sort | uniq | perl -npe 's/\\n/,/'`
          [ "$result" = "0,21," ]
        }
        
        
        #pogos: the size of all features (this has been checked)
        @test "pogos8" {
         result=`cat pogos/00_pogosstats_* | cut -f 5 | sort -n | uniq | grep -v cov | perl -npe 's/\\n/,/'`
          [ "$result" = "0,21,48,53,92,100,113,187,200,300,400,700,1587," ]
        }
        
        #pogos: 40 peaks tested for chromosome 2
        @test "pogos9" {
         result=`cut -f 1,2,4 pogos/00_pogos*txt | grep Promoter | grep chr2 | cut -f 3`
          [ "$result" -eq 40 ]
        }
        
        #pogos: 45 peaks on chromosome 1
        @test "pogos10" {
         result=`cut -f 1,2,4 pogos/00_pogos*txt | grep Promoter | grep chr1 | cut -f 3`
          [ "$result" -eq 45 ]
        }
        
        #pogos: Chromosomal coverage of promoter is 48 on chr1
        @test "pogos11" {
         result=`cut -f 1,2,5 pogos/00_pogos*txt | grep Promoter | grep chr1 | cut -f 3`
          [ "$result" -eq 48 ]
        }
        
        
        #pogos: Coverage of intergenic regions is 187 nuc
        @test "pogos12" {
         result=`cut -f 1,2,5 pogos/00_pogos*txt | grep Inter | grep chr1 | cut -f 3`
          [ "$result" -eq 187 ]
        }
        
        #pogos: Number of peaks midpoints intersecting intergenic is 19
        @test "pogos13" {
         result=`cut -f 1,2,6 pogos/00_pogos*txt | grep Inter | grep chr1 | cut -f 3`
          [ "$result" -eq 19 ]
        }
        
        
        #pogos: Coverage of intronic regions is 21 nuc
        @test "pogos14" {
         result=`cut -f 1,2,5 pogos/00_pogos*txt | grep Intro | grep chr1 | cut -f 3`
          [ "$result" -eq 21 ]
        }
        
        #pogos: The number of peaks midpoints falling in introns is 1
        @test "pogos15" {
         result=`cut -f 1,2,6 pogos/00_pogos*txt | grep Intro | grep chr1 | cut -f 3`
          [ "$result" -eq 1 ]
        }
        
        # This test does not testpogosand currently fail on some systems (issue 184). Its output is already git-controlled so skipping it allows to pursue on following tests.
        ##pogos: check more keys
        #@test "pogos15a" {
        # result=`gtftk nb_exons -i pygtftk/data/simple_02/simple_02.gtf -g > pygtftk/data/simple_02/simple_02_nbe.gtf`
        #  [ "$result" = "" ]
        #}
        
        #pogos: check more keys
        @test "pogos16" {
         result=`rm -Rf peak_annotation;  gtftkpogos -i pygtftk/data/simple_02/simple_02_nbe.gtf -p pygtftk/data/simple_02/simple_02_peaks.bed  -K peak_annotation --more-keys nb_exons -c pygtftk/data/simple_02/simple_02.chromInfo`
          [ "$result" = "" ]
        }
    
        #pogos: check more keys
        @test "pogos17" {
         result=`cat pogos/*txt| grep "nb_exons"| wc -l`
          [ "$result" -eq 18 ]
        }
        
        #pogos: check more keys
        @test "pogos18" {
         result=`cat pogos/*txt| grep "nb_exons"| wc -l`
          [ "$result" -eq 18 ]
        }
        
        #pogos: check no basic feature
        @test "pogos19" {
         result=`rm -Rf peak_annotation; gtftkpogos -i pygtftk/data/simple_02/simple_02_nbe.gtf -p pygtftk/data/simple_02/simple_02_peaks.bed  -K peak_annotation --more-keys nb_exons -c pygtftk/data/simple_02/simple_02.chromInfo -n`
          [ "$result" = "" ]
        }
        
        #pogos: check no basic feature
        @test "pogos20" {
         result=`cat pogos/*txt|  wc -l`
          [ "$result" -eq 24 ]
        }
    '''

    cmd = CmdObject(name="pogos",
                    message="Statistics on bed file intersections with genomic features.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="annotation",
                    notes=__notes__,
                    updated=__updated__,
                    test=test)
