
# FOR MY JUPYTER KERNEL TESTING ONLY
import os
os.getcwd()
os.chdir('/home/ferre/anaconda3/lib/python3.6/site-packages/pygtftk-0.9.8-py3.6-linux-x86_64.egg/pygtftk')


























#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import math
import os
import re
import sys
import warnings
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd
from mizani.formatters import scientific_format
from plotnine import (ggplot, aes, position_dodge,
                      geom_bar, ggsave, ylab, theme, element_blank, element_text, scale_fill_manual,
                      guides, guide_legend, geom_label, geom_text, geom_point, scale_color_manual,
                      scale_y_continuous)
from scipy.stats import binom_test

from pygtftk.arg_formatter import FileWithExtension, bedFileList
from pygtftk.arg_formatter import bed6_or_bed3
from pygtftk.arg_formatter import checkChromFile
from pygtftk.arg_formatter import int_greater_than_null_or_None, int_greater_than_null
from pygtftk.bedtool_extension import BedTool
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chrom_info_as_dict
from pygtftk.utils import chrom_info_to_bed_file
from pygtftk.utils import close_properly
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message


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




    # TODO : chromfile will return a file, use the arg_formatter module






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

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=int_greater_than_null,
                            default=300,
                            required=False)

    parser_grp.add_argument('-r', '--order-bar',
                            help='The way bar should be ordered.',
                            choices=['log2_ratio','pval_binom'],
                            default='log2_ratio',
                            required=False)

    return parser


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
              dpi=300,
              order_bar=None,
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
    # -------------------------------------------------------------------------

    chrom_len = chrom_info_as_dict(chrom_info)


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
                                                "00_peak_anno_diagrams." + page_format
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
    for i in region_mid_point:
        chrom_list.add(i.chrom)

    for i in chrom_list:
        if i not in chrom_len:
            message("Chromosome " + " i from GTF is undefined in --chrom-info file.",
                    type="ERROR")

    # -------------------------------------------------------------------------
    # Fill the dict with info about basic features include in GTF
    # -------------------------------------------------------------------------


    from functools import partial
    # Import the main function from the stats.intersect module
    from pygtftk.stats.intersect.overlap_stats_shuffling import compute_overlap_stats

    # PARTIAL TODO EXPLAIN
    my_intersection = partial(compute_overlap_stats, chrom_len=chrom_len)


    # OLD CALLER
    # hits = _intersection_results(peak_file=region_mid_point.fn,
    #                              feature_bo=gtf_sub_bed,
    #                              my_dict=hits,
    #                              chrom_len=chrom_len,
    #                              ft_type="Intergenic")

    if not no_basic_feature:
        for feat_type in gtf.get_feature_list(nr=True):
            gtf_sub = gtf.select_by_key("feature", feat_type, 0)

            gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                               "gene_id",
                                               "exon_id"]).sort().merge()  # merging bed file !

            # PERSONAL NOTE : gtf_sub_bed is already a BedTool object, no need to import it (I worked with filepaths previously)


            # REPLACE WITH : hits[feat_type] = pygtftk.stats.intersect.overlap_stats_shuffling # import the function before
            hits[feat_type] = my_intersection(peak_file=region_mid_point.fn, feature_bo=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the intergenic regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_intergenic(chrom_info,
                                         0,
                                         0,
                                         chrom_len.keys()).merge()

        hits["Intergenic"] = my_intersection(peak_file=region_mid_point.fn, feature_bo=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the intronic regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_introns()

        hits["Introns"] = my_intersection(peak_file=region_mid_point.fn, feature_bo=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the promoter regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_tss().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits["Promoters"] = my_intersection(peak_file=region_mid_point.fn, feature_bo=gtf_sub_bed)

        # -------------------------------------------------------------------------
        # Get the tts regions
        # -------------------------------------------------------------------------

        gtf_sub_bed = gtf.get_tts().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits["Terminator"] = my_intersection(peak_file=region_mid_point.fn, feature_bo=gtf_sub_bed)



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
                    ft_type=":".join([user_key,val]) # Key for the dictionary
                    hits[ft_type] = my_intersection(peak_file=region_mid_point.fn,
                                                 feature_bo=gtf_sub_bed)

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

            hits[bed_lab] = my_intersection(peak_file=region_mid_point.fn,
                                         feature_bo=BedTool(bed_anno.name))





    # ------------------ Treating the 'hits' dictionary --------------------- #











    # personal note : bayesian fitting and true intersection can be done here. Just replace this obsolete code

    if len(hits) == 0:
        message("No feature found.", type="ERROR")

    nb_peaks = str(len(region_mid_point))
    # Do all my math here

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

    # -------------------------------------------------------------------------
    # Read the data set and plot it
    # -------------------------------------------------------------------------


    # Personal note : this is simply a pandas version of my output_intersect.txt
    d = pd.read_csv(data_file.name, sep="\t", header=0)

    plot_results(d)



































# NOTE : my code should output something very much like d, so the plotnine
# code below can be re-used as is


def plot_results(d):
    """
    Main plotting function by D.Puthier
    """

    if d.shape[0] == 0:
        message("No lines found in input file.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Compute expected number of intersections
    # -------------------------------------------------------------------------

    d['freq'] = d['coverage'] / d['reference_size']
    d['Expected'] = d['freq'] * d['Nb_trial_or_peak']

    # -------------------------------------------------------------------------
    # Compute binomial p.val  (unilateral)
    # -------------------------------------------------------------------------


    # i are the features
    for i, _ in d.iterrows():

        ## Compute a log ratio
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



        ## Compute the binomial test
        if d.loc[i, 'Nb_trial_or_peak'] > 0:
            pval = binom_test(n=d.loc[i, 'Nb_trial_or_peak'],
                              x=d.loc[i, 'Observed'],
                              p=d.loc[i, 'freq'])

            d.loc[i, 'pval_binom'] = pval
            d.loc[i, 'pval_binom_str'] = "{0:0.3g}".format(pval)

        ## Enriched or depleted ?
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

    def bar_plot(d=None,
                 x_ft_type=None,
                 which_row=None):

        # -------------------------------------------------------------------------
        # Subset the data
        # -------------------------------------------------------------------------

        d_sub = d[which_row].copy()

        # -------------------------------------------------------------------------
        # Compute text position
        # -------------------------------------------------------------------------

        max_y = max(d_sub['Expected'].tolist() + d_sub['Observed'].tolist())

        offset = 10 / 100 * max_y

        for i, _ in d_sub.iterrows():
            ## max_y to be use to plot pval
            max_y_col = max(d_sub.loc[i, 'Observed'],
                            d_sub.loc[i, 'Expected'])

            d_sub.loc[i, 'y_log2_ratio'] = max_y_col + offset * 2.4
            d_sub.loc[i, 'y_pval'] = max_y_col + offset * 1.8
            d_sub.loc[i, 'y_text'] = max_y_col + offset / 10

            ## y_lim to set diagram limits in y coords.
            d_sub.loc[i, 'y_lim'] = max_y + offset * 2.2

        # -------------------------------------------------------------------------
        # Melt the data frame
        # -------------------------------------------------------------------------

        message('Melting.')
        dm = d_sub.melt(id_vars=[x for x in d_sub.columns if x not in ['Expected',
                                                                       'Observed']],
                        value_vars=['Observed', 'Expected'])

        dm['variable'] = pd.Categorical(dm['variable'])

        # -------------------------------------------------------------------------
        # Create a new plot
        # -------------------------------------------------------------------------

        p = ggplot()

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

        aes_plot = aes(x_ft_type, 'value', fill='variable')

        p += geom_bar(data=dm,
                      mapping=aes_plot,
                      stat='identity',
                      alpha=0.6,
                      position='dodge',
                      show_legend=True)

        p += ylab("Number of overlaps")
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

        # -------------------------------------------------------------------------
        # Display text (observed vs expected)
        # -------------------------------------------------------------------------

        message('Adding text (observed vs expected)')

        aes_plot = aes(x=x_ft_type,
                       y='y_text',
                       label='value',
                       colour='variable')

        dodge_text = position_dodge(width=0.9)

        p += geom_text(data=dm,
                       mapping=aes_plot,
                       format_string='{0:.3g}',
                       position=dodge_text,
                       angle=90,
                       va='bottom',
                       ha='center',
                       size=4,
                       show_legend=False)

        # -------------------------------------------------------------------------
        # Display text (pval)
        # -------------------------------------------------------------------------

        message('Adding text (p-values)')

        aes_plot = aes(x=x_ft_type,
                       y='y_pval',
                       label='pval_binom_str',
                       fill='test_type')

        p += geom_label(data=dm[dm['variable'] == 'Observed'],
                        mapping=aes_plot,
                        position='identity',
                        colour="white",
                        angle=0,
                        va='bottom',
                        ha='center',
                        size=4,
                        show_legend=False)

        # -------------------------------------------------------------------------
        # Display text (log2_ratio)
        # -------------------------------------------------------------------------

        message('Adding text (log2_ratio)')

        aes_plot = aes(x=x_ft_type,
                       y='y_log2_ratio',
                       label='log2_ratio_str',
                       fill='test_type')

        p += geom_label(data=dm[dm['variable'] == 'Observed'],
                        mapping=aes_plot,
                        position='identity',
                        colour="white",
                        angle=0,
                        va='bottom',
                        ha='center',
                        size=4,
                        show_legend=False)

        # -------------------------------------------------------------------------
        # Set some constrains on y_lim
        # -------------------------------------------------------------------------

        p += geom_point(data=dm[dm['variable'] == 'Observed'],
                        mapping=aes(x=x_ft_type,
                                    y='y_lim'),
                        alpha=0,
                        colour="white",
                        show_legend=False)

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

        p += guides(colour=False, fill=guide_legend(ncol=2, byrow=True))

        # p += scale_color_discrete(l=.4)

        # -------------------------------------------------------------------------
        # y axis labels in scientific notation
        # -------------------------------------------------------------------------

        p += scale_y_continuous(labels=scientific_format(digits=2))

        # -------------------------------------------------------------------------
        # return
        # -------------------------------------------------------------------------

        return (p)

    # -------------------------------------------------------------------------
    # Compute bar plot by feature
    # -------------------------------------------------------------------------

    p = bar_plot(d=d,
                 x_ft_type='ft_type',
                 which_row=d['chrom'] == 'all_chrom')

    # -------------------------------------------------------------------------
    #
    # Computing page size
    #
    # -------------------------------------------------------------------------

    nb_ft = len(list(d['ft_type'].unique()))

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
    peak_anno(**args)


if __name__ == '__main__':
    main()


else:


    # Do my unitary test here (in bash)
    # Add some doctests in the code itself (I think the module is already loaded, just write the doctests, check with Denis)

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
