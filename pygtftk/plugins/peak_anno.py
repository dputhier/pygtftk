#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import re
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from plotnine import (ggplot, aes, position_dodge,
                      geom_bar, ggsave, ylab, theme, element_blank, element_text, scale_fill_manual,
                      scale_color_manual, guides, guide_legend, geom_label, geom_text, scale_y_log10)
from scipy.stats import binom_test

from pygtftk.arg_formatter import FileWithExtension, bedFileList
from pygtftk.arg_formatter import bed6_or_bed3
from pygtftk.arg_formatter import checkChromFile
from pygtftk.arg_formatter import int_greater_than_null_or_None
from pygtftk.bedtool_extension import BedTool
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chrom_info_as_dict
from pygtftk.utils import chrom_info_to_bed_file
from pygtftk.utils import close_properly
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
                            type=FileWithExtension('r',
                                                   valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

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
                            action=checkChromFile,
                            required=False)

    parser_grp.add_argument('-p', '--peak-file',
                            help='The file containing the peaks/regions to be annotated.'
                                 ' (bed format).',
                            default=None,
                            metavar="BED",
                            action=bed6_or_bed3,
                            required=True)

    parser_grp.add_argument('-b', '--more-bed',
                            help="A comma separated list of bed files to be "
                                 "considered as additional genomic annotations.",
                            action=bedFileList,
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

    parser_grp.add_argument('-pw', '--pdf-width',
                            help='Output pdf file width (inches).',
                            type=int_greater_than_null_or_None,
                            default=None,
                            required=False)

    parser_grp.add_argument('-ph', '--pdf-height',
                            help='Output pdf file height (inches).',
                            type=int_greater_than_null_or_None,
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
                            type=argparse.FileType("w"),
                            required=False)

    parser_grp.add_argument('-pf', '--page-format',
                            help='Output file format.',
                            choices=['pdf', 'png'],
                            default='pdf',
                            required=False)

    return parser


# -------------------------------------------------------------------------
# This function performs intersection between a bedTools object
# an another file (feature_file) e.g promoter, intron, exon...
# Fill a dictionary of dictionary
# -------------------------------------------------------------------------


def _intersection_results(peak_file=None,
                          feature_bo=None,
                          my_dict=None,
                          chrom_len=None,
                          ft_type=None):
    """Get feature and peaks and compute intersections. Returns an updated
    dict.
    """

    # -------------------------------------------------------------------------
    # Compute the number of regions of interest falling in
    # each chromosome.
    # If ft_type = "Full_chromosomes" we are checking whether peaks
    # tend to fall on some chromosomes (i.e. not a the feature level)
    # -------------------------------------------------------------------------

    peak_file = BedTool(peak_file)
    nb_peaks = len(peak_file)

    for i in peak_file:
        if ft_type != "Full_chromosomes":
            my_dict[ft_type][i.chrom]["Nb_trial_or_peak"] += 1
            my_dict[ft_type]["all_chrom"]["Nb_trial_or_peak"] += 1
        else:
            my_dict[ft_type][i.chrom]["Nb_trial_or_peak"] = nb_peaks

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
        if ft_type != "Full_chromosomes":
            my_dict[ft_type][i.chrom]["coverage"] += size
            my_dict[ft_type]["all_chrom"]["coverage"] += size
        else:
            # We will compare "Chromosomes" to the full genome.
            my_dict[ft_type][i.chrom]["coverage"] = chrom_len[i.chrom]

    # -------------------------------------------------------------------------
    # Compute intersections
    # -------------------------------------------------------------------------

    intersections = peak_file.intersect(feature_merge_bo)

    for i in intersections:

        chrom = i.chrom
        my_dict[ft_type][chrom]["Observed"] += 1
        if ft_type != "Full_chromosomes":
            my_dict[ft_type]["all_chrom"]["Observed"] += 1

    file_out_save = make_tmp_file(prefix="peak_anno_intersections_" +
                                         ft_type,
                                  suffix=".bed")
    intersections.saveas(file_out_save.name)

    return my_dict


# -------------------------------------------------------------------------
# The command function
# -------------------------------------------------------------------------


def peak_anno(inputfile=None,
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
              chrom_info=None,
              user_img_file=None,
              page_format=None,
              tmp_dir=None,
              logger_file=None,
              verbosity=True):
    """
    This function is intended to perform statistics on peak intersection. It will compare your peaks to
    classical features (e.g promoter, tts, gene body, UTR,...) and to sets of user provided peaks.
    """

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
    # As one diagram inspect enrichment on a per-chromosome basis
    # we can not accept too much chromosomes
    # -------------------------------------------------------------------------

    chrom_len = chrom_info_as_dict(chrom_info)

    if len(chrom_len.keys()) > 100:
        message("Can't accept more than 100 chromosomes in peak_anno (see --chrom-info).",
                type="ERROR")

    # -------------------------------------------------------------------------
    # The hits variable is a muti-level dict to store the results
    # per chromosome
    #   - e.g my_dict[ft_type][chrom]["Observed"]
    #   - e.g my_dict[ft_type][chrom]["coverage"]
    #   - e.g my_dict[ft_type][chrom]["Nb_trial_or_peaks"]
    # -------------------------------------------------------------------------

    def nested_dict(n, type):
        """"http://stackoverflow.com/questions/29348345"""
        if n == 1:
            return defaultdict(type)
        else:
            return defaultdict(lambda: nested_dict(n - 1, type))

    hits = nested_dict(3, int)

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
                                                "00_peak_anno_diagrams." + page_format,
                                                "00_peak_anno_diagrams_by_chrom." + page_format
                                                ],
                                         force=True)

    data_file, pdf_file, pdf_file_by_chrom = file_out_list

    # -------------------------------------------------------------------------
    # Get the midpoints of the peaks
    # -------------------------------------------------------------------------

    region_mid_point_file = make_tmp_file("peaks_midpoints", ".bed")

    # Loop through peaks
    region_mid_point = BedTool(peak_file).get_midpoints()

    region_mid_point.saveas(region_mid_point_file.name)

    # -------------------------------------------------------------------------
    # Check chromosomes for peaks are defined in the chrom-info file
    # -------------------------------------------------------------------------

    chrom_list = set()
    for i in region_mid_point:
        chrom_list.add(i.chrom)

    for i in chrom_list:
        if i not in chrom_len:
            message("Chromosome " + " i from GTF is undefined in --chrom-info file.",
                    type="ERROR")

    # -------------------------------------------------------------------------
    # Fill the dict with info about basic features include in GTF
    # -------------------------------------------------------------------------

    if not no_basic_feature:
        for feat_type in gtf.get_feature_list(nr=True):
            if feat_type not in ["start_codon", "stop_codon"]:
                gtf_sub = gtf.select_by_key("feature", feat_type, 0)

                gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                                   "gene_id",
                                                   "exon_id"]).sort().merge()  # merging bed file !

                hits = _intersection_results(peak_file=region_mid_point.fn,
                                             feature_bo=gtf_sub_bed,
                                             my_dict=hits,
                                             chrom_len=chrom_len,
                                             ft_type=feat_type)

        # -------------------------------------------------------------------------
        # Get the intergenic regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_intergenic(chrom_info,
                                         0,
                                         0,
                                         chrom_len.keys()).merge()

        hits = _intersection_results(peak_file=region_mid_point.fn,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Intergenic")

        # -------------------------------------------------------------------------
        # Get the intronic regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_introns()

        hits = _intersection_results(peak_file=region_mid_point.fn,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Introns")

        # -------------------------------------------------------------------------
        # Get the promoter regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_tss().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits = _intersection_results(peak_file=region_mid_point.fn,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Promoters")

        # -------------------------------------------------------------------------
        # Get the tts regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_tts().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits = _intersection_results(peak_file=region_mid_point.fn,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Terminator")

        # -------------------------------------------------------------------------
        # Test the whole chromosomes as features
        # -------------------------------------------------------------------------

        gtf_sub_bed = BedTool(
            chrom_info_to_bed_file(
                chrom_info,
                chr_list=chrom_len.keys()))

        hits = _intersection_results(peak_file=region_mid_point.fn,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Full_chromosomes")

    # -------------------------------------------------------------------------
    # if the user request --more-keys (e.g. gene_biotype)
    # -------------------------------------------------------------------------

    if more_keys is not None:
        more_keys_list = more_keys.split(",")
        if len(more_keys_list) > 50:
            message(
                "The selected key in --more-keys should be associated with less than 50 different values.")
        for user_key in more_keys_list:
            user_key_values = set(gtf.extract_data(user_key,
                                                   as_list=True))
            for val in user_key_values:

                gtf_sub = gtf.select_by_key(user_key, val, 0)

                if len(gtf_sub) > 0:
                    gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                                       "gene_id",
                                                       "exon_id"]).sort().merge()  # merging bed file !

                    hits = _intersection_results(peak_file=region_mid_point.fn,
                                                 feature_bo=gtf_sub_bed,
                                                 my_dict=hits,
                                                 chrom_len=chrom_len,
                                                 ft_type=":".join([user_key,
                                                                   val]))
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

            hits = _intersection_results(peak_file=region_mid_point.fn,
                                         feature_bo=BedTool(bed_anno.name),
                                         my_dict=hits,
                                         chrom_len=chrom_len,
                                         ft_type=bed_lab)

    # Store the result into a file
    # before loading it into R
    message("Storing data in: " + data_file.name)

    data_file.write("\t".join(["ft_type",
                               "chrom",
                               "reference_size",
                               "Nb_trial_or_peak",
                               "coverage",

                               "Observed\n"]))

    if len(hits) == 0:
        message("No feature found.", type="ERROR")

    nb_peaks = str(len(region_mid_point))

    for chrom in chrom_len:
        if chrom != "all_chrom":
            hits["Full_chromosomes"][chrom]["Nb_trial_or_peak"] = nb_peaks

    for key1 in hits:
        for key2 in chrom_len:
            if key1 == "Full_chromosomes":
                ref_size = str(chrom_len["all_chrom"])
                nb_trial = nb_peaks
            else:
                ref_size = str(chrom_len[key2])
                nb_trial = str(hits[key1][key2]["Nb_trial_or_peak"])

            if not (key1 == "Full_chromosomes" and key2 == "all_chrom"):
                out_list = [key1,
                            key2,
                            ref_size,
                            nb_trial,
                            str(hits[key1][key2]["coverage"]),
                            str(hits[key1][key2]["Observed"])]

                data_file.write("\t".join(out_list) + "\n")

    close_properly(data_file)

    if user_img_file is not None:
        os.unlink(pdf_file.name)
        os.unlink(pdf_file_by_chrom.name)
        pdf_file = user_img_file
        pdf_file_by_chrom = pdf_file.name.replace("." + page_format, "_by_chrom." + page_format)
        pdf_file_by_chrom = open(pdf_file_by_chrom, 'w')
        if not pdf_file.name.endswith(page_format):
            msg = "Image format: {f}. Please fix.".format(f=page_format)
            message(msg, type="ERROR")

    # -------------------------------------------------------------------------
    # Read the data set
    # -------------------------------------------------------------------------

    d = pd.read_csv(data_file.name, sep="\t", header=0)

    if d.shape[0] == 0:
        message("No lines found in input file.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Get the chromosome order
    # -------------------------------------------------------------------------

    # Compute expected number of intersections
    d['freq'] = d['coverage'] / d['reference_size']
    d['Expected'] = d['freq'] * d['Nb_trial_or_peak']

    # -------------------------------------------------------------------------
    # Compute binomial p.val  (unilateral)
    # -------------------------------------------------------------------------

    for i, _ in d.iterrows():

        if d.loc[i, 'Expected'] > 0:
            if d.loc[i, 'Observed'] > 0:
                log2_ratio = np.log2(d.loc[i, 'Observed'] / d.loc[i, 'Expected'])
            else:
                log2_ratio = np.nan
        else:
            log2_ratio = np.nan

        d.loc[i, 'log2_ratio'] = log2_ratio

        if d.loc[i, 'Nb_trial_or_peak'] > 0:
            pval = binom_test(n=d.loc[i, 'Nb_trial_or_peak'],
                              x=d.loc[i, 'Observed'],
                              p=d.loc[i, 'freq'])

            d.loc[i, 'pval_binom'] = pval
            d.loc[i, 'pval_binom_str'] = "{0:0.3g}".format(pval)

        if d.loc[i, 'Observed'] > d.loc[i, 'Expected']:
            d.loc[i, 'test_type'] = 'Enrichment'
        elif d.loc[i, 'Observed'] < d.loc[i, 'Expected']:
            d.loc[i, 'test_type'] = 'Depletion'
        else:
            d.loc[i, 'test_type'] = 'Unchanged'

        d.loc[i, 'max_y'] = max(d.loc[i, 'Observed'], d.loc[i, 'Expected'])

    # -------------------------------------------------------------------------
    # Save file
    # -------------------------------------------------------------------------

    d.to_csv(open(data_file.name, 'w'), sep="\t", header=True, index=False)

    # -------------------------------------------------------------------------
    # Melt the data frame
    # -------------------------------------------------------------------------

    import pickle
    test_file = open("toto.pick", 'wb')
    pickle.dump(d, test_file)
    test_file.close()
    message('Melting.')
    dm = d.melt(id_vars=[x for x in d.columns if x not in ['Observed',
                                                           'Expected']],
                value_vars=['Observed', 'Expected'])

    dm['variable'] = pd.Categorical(dm['variable'])

    # -------------------------------------------------------------------------
    # Prepare a function for drawing diagram
    # -------------------------------------------------------------------------

    def bar_plot(dm=None,
                 ft_type=None,
                 which_row=None,
                 by_chrom=False,
                 y_axis_trans=None):
        p = ggplot()

        # -------------------------------------------------------------------------
        # Display bars
        # -------------------------------------------------------------------------

        message('Adding bar plot.')

        if not by_chrom:
            aes_plot = aes('ft_type', 'value', fill='variable')
        else:
            aes_plot = aes('chrom', 'value', fill='variable')

        p += geom_bar(data=dm[which_row],
                      mapping=aes_plot,
                      stat='identity',
                      position='dodge',
                      alpha=0.6,
                      show_legend=True)

        p += ylab("Number of overlaps")
        p += theme(legend_title=element_blank(),
                   legend_position="bottom",
                   legend_box_spacing=0.65,
                   legend_key_size=8,
                   legend_text=element_text(size=8),
                   legend_key=element_blank(),
                   axis_title_x=element_blank(),
                   axis_title_y=element_text(colour="#333333",
                                             size=8,
                                             hjust=4,
                                             angle=90,
                                             face="plain"),
                   axis_text_y=element_text(size=8,
                                            margin={'r': 5},
                                            angle=0),
                   axis_text_x=element_text(size=8,
                                            margin={'t': 5, 'r': 5},
                                            ha='right',
                                            angle=45))

        p += scale_y_log10()

        # -------------------------------------------------------------------------
        # Display text (observed vs expected)
        # -------------------------------------------------------------------------

        message('Adding text (observed vs expected)')

        if not by_chrom:
            x_ft_type = 'ft_type'
        else:
            x_ft_type = 'chrom'

        aes_plot = aes(x='ft_type',
                       y='value',
                       label='value',
                       colour='variable')

        dodge_text = position_dodge(width=0.9)

        p += geom_text(data=dm[which_row],
                       mapping=aes_plot,
                       format_string='{0:.2f}',
                       position=dodge_text,
                       angle=0,
                       va='bottom',
                       ha='center',
                       size=4)

        # -------------------------------------------------------------------------
        # Display text (pval)
        # -------------------------------------------------------------------------

        message('Adding text (p-values)')

        aes_plot = aes(x='ft_type',
                       y='max_y',
                       label='pval_binom_str',
                       colour='test_type')

        p += geom_label(data=dm[(which_row) & (dm['variable'] == 'Observed')],
                        mapping=aes_plot,
                        position='identity',
                        angle=0,
                        nudge_y=0.1,
                        va='bottom',
                        ha='center',
                        size=4)

        # -------------------------------------------------------------------------
        # Set color scale
        # -------------------------------------------------------------------------

        col_dict = {'Observed': 'blue',
                    'Expected': '#8A8A8A',
                    'Enrichment': '#990000',
                    'Unchanged': '#000000',
                    'Depletion': '#008800'}

        p += scale_fill_manual(values=col_dict)

        p += scale_color_manual(values=col_dict)

        p += guides(colour=False, fill=guide_legend(ncol=2))

        # p += scale_color_discrete(l=.4)

        # -------------------------------------------------------------------------
        # return
        # -------------------------------------------------------------------------

        return (p)

    # -------------------------------------------------------------------------
    # Compute bar plot by feature
    # -------------------------------------------------------------------------

    p = bar_plot(dm=dm,
                 ft_type=list(d['ft_type'].unique()),
                 which_row=dm['chrom'] == 'all_chrom',
                 by_chrom=False,
                 y_axis_trans=None)

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
               height=pdf_height)
        # dm.to_csv(data_file, sep="\t", header=True, index=False)

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

    test = '''
            
        #peak_anno: chr2 len
        @test "peak_anno_1" {
             result=`rm -Rf peak_annotation; gtftk peak_anno  -i pygtftk/data/simple_02/simple_02.gtf -p pygtftk/data/simple_02/simple_02_peaks.bed -c pygtftk/data/simple_02/simple_02.chromInfo -u 2 -d 2 -K peak_annotation`
          [ "$result" = "" ]
        }
        
        
        #peak_anno: all_chrom len
        @test "peak_anno_2" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep chr2 | cut -f 3 | sort | uniq | perl -npe 's/\\n/,/'`
          [ "$result" = "1700,400," ]
        }
        
        
        #peak_anno: all_chrom len
        @test "peak_anno_3" {
         result=`cat peak_annotation/00_peak_anno_stats_* | grep all_chrom | cut -f 3 | sort | uniq`
          [ "$result" = "1700" ]
        }
        
        #peak_anno: there are 11 guys (midpoints) intersecting CDS
        @test "peak_anno_4" {
         result=`grep -w CDS peak_annotation/00_peak_anno_stats_* | cut -f 6 | perl -npe 's/\\n/,/'`
          [ "$result" = "0,0,0,11,11,0," ]
        }
        
        
        #peak_anno: there is 59 guys (midpoints) intersecting intergenic regions (and 19 on chr1).
        @test "peak_anno_5" {
         result=`grep -i Interg peak_annotation/00_peak_anno_stats_* | cut -f 6 | perl -npe 's/\\n/,/'`
          [ "$result" = "0,0,40,19,59,0," ]
        }
        
        
        
        #peak_anno: size of intronic features is 21 bases
        @test "peak_anno_7" {
         result=`grep -i Introns peak_annotation/00_peak_anno_stats_* | cut -f 5 | sort | uniq | perl -npe 's/\\n/,/'`
          [ "$result" = "0,21," ]
        }
        
        
        #peak_anno: the size of all features (this has been checked)
        @test "peak_anno_8" {
         result=`cat peak_annotation/00_peak_anno_stats_* | cut -f 5 | sort -n | uniq | grep -v cov | perl -npe 's/\\n/,/'`
          [ "$result" = "0,21,48,53,92,100,113,187,200,300,400,700,1587," ]
        }
        
        #peak_anno: 40 peaks tested for chromosome 2
        @test "peak_anno_9" {
         result=`cut -f 1,2,4 peak_annotation/00_peak_anno_*txt | grep Promoter | grep chr2 | cut -f 3`
          [ "$result" -eq 40 ]
        }
        
        #peak_anno: 45 peaks on chromosome 1
        @test "peak_anno_10" {
         result=`cut -f 1,2,4 peak_annotation/00_peak_anno_*txt | grep Promoter | grep chr1 | cut -f 3`
          [ "$result" -eq 45 ]
        }
        
        #peak_anno: Chromosomal coverage of promoter is 48 on chr1
        @test "peak_anno_11" {
         result=`cut -f 1,2,5 peak_annotation/00_peak_anno_*txt | grep Promoter | grep chr1 | cut -f 3`
          [ "$result" -eq 48 ]
        }
        
        
        #peak_anno: Coverage of intergenic regions is 187 nuc
        @test "peak_anno_12" {
         result=`cut -f 1,2,5 peak_annotation/00_peak_anno_*txt | grep Inter | grep chr1 | cut -f 3`
          [ "$result" -eq 187 ]
        }
        
        #peak_anno: Number of peaks midpoints intersecting intergenic is 19
        @test "peak_anno_13" {
         result=`cut -f 1,2,6 peak_annotation/00_peak_anno_*txt | grep Inter | grep chr1 | cut -f 3`
          [ "$result" -eq 19 ]
        }
        
        
        #peak_anno: Coverage of intronic regions is 21 nuc
        @test "peak_anno_14" {
         result=`cut -f 1,2,5 peak_annotation/00_peak_anno_*txt | grep Intro | grep chr1 | cut -f 3`
          [ "$result" -eq 21 ]
        }
        
        #peak_anno: The number of peaks midpoints falling in introns is 1
        @test "peak_anno_15" {
         result=`cut -f 1,2,6 peak_annotation/00_peak_anno_*txt | grep Intro | grep chr1 | cut -f 3`
          [ "$result" -eq 1 ]
        }
        
        # This test does not test peak_anno and currently fail on some systems (issue 184). Its output is already git-controlled so skipping it allows to pursue on following tests.
        ##peak_anno: check more keys
        #@test "peak_anno_15a" {
        # result=`gtftk nb_exons -i pygtftk/data/simple_02/simple_02.gtf -g > pygtftk/data/simple_02/simple_02_nbe.gtf`
        #  [ "$result" = "" ]
        #}
        
        #peak_anno: check more keys
        @test "peak_anno_16" {
         result=`rm -Rf peak_annotation;  gtftk peak_anno  -i pygtftk/data/simple_02/simple_02_nbe.gtf -p pygtftk/data/simple_02/simple_02_peaks.bed  -K peak_annotation --more-keys nb_exons -c pygtftk/data/simple_02/simple_02.chromInfo`
          [ "$result" = "" ]
        }
    
        #peak_anno: check more keys
        @test "peak_anno_17" {
         result=`cat peak_annotation/*txt| grep "nb_exons"| wc -l`
          [ "$result" -eq 18 ]
        }
        
        #peak_anno: check more keys
        @test "peak_anno_18" {
         result=`cat peak_annotation/*txt| grep "nb_exons"| wc -l`
          [ "$result" -eq 18 ]
        }
        
        #peak_anno: check no basic feature
        @test "peak_anno_19" {
         result=`rm -Rf peak_annotation; gtftk peak_anno  -i pygtftk/data/simple_02/simple_02_nbe.gtf -p pygtftk/data/simple_02/simple_02_peaks.bed  -K peak_annotation --more-keys nb_exons -c pygtftk/data/simple_02/simple_02.chromInfo -n`
          [ "$result" = "" ]
        }
        
        #peak_anno: check no basic feature
        @test "peak_anno_20" {
         result=`cat peak_annotation/*txt|  wc -l`
          [ "$result" -eq 24 ]
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
                    test=test,
                    rlib=R_LIBS)
