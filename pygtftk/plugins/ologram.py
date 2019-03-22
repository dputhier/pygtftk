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
from plotnine import (ggplot, aes, position_dodge, ggtitle,
                      geom_bar, ylab, theme, element_blank,
                      element_text, geom_errorbar, theme_bw,
                      geom_label, save_as_pdf_pages, scale_fill_manual)

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
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

__updated__ = "2019-03-18"
__doc__ = """

 OLOGRAM -- OverLap Of Genomic Regions Analysis using Monte Carlo. Ologram 
 annotates peaks (in bed format) using (i) genomic features extracted 
 from a GTF file (e.g promoter, tts, gene body, UTR...) (ii) genomic regions tagged with 
  particular keys/values in a GTF file (e.g. gene_biotype "protein_coding", 
  gene_biotype "LncRNA"...) or (iii) from a BED file (e.g. user-defined regions).

 Each couple peak file/region is randomly shuffled across the genome (inter-region
 lengths are considered). Then the probability of intersection under the null
 hypothesis (the peaks and this feature are independent) is deduced thanks to
 this Monte Carlo approach.

 The program will return statistics for both the number of intersections and the
 total lengths (in basepairs) of all intersections.

 Authors : Quentin Ferré <quentin.q.ferre@gmail.com>, Guillaume Charbonnier 
 <guillaume.charbonnier@outlook.com> and Denis Puthier <denis.puthier@univ-amu.fr>.
 """

__notes__ = """
 -- Although ologram itself is not RAM-intensive, base pygtftk processing of a full human GTF can require upwards of 8Gb.
 It is recommended you do not run other programs in the meantime on a laptop.

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
    to a Negative Binomial assessed by a Cramer's V score (fit_quality gives 1-V ; as per Cramer (1948) a good fit
    should have a fit quality above (1 - 0.25 = 0.75) if your nb. of shuffles is in the hundreds, but closer to 0.9
    if it is in the thousands or above.

    The p-value of the true intersection under the distribution characterized by the shuffles is also given, under 'p_value'.
    Finally, the log2 fold change between true and shuffles is also given.

 -- If -\-more-keys is used additional region sets will be tested based on the associated key value.
 As an example, if -\-more-keys is set to the 'gene_biotype' (a key generally found in ensembl GTF), the
 region related to 'protein_coding', 'lncRNA' or any other values for that key will be retrieved merged and tested
 for enrichment.

 -- Use -\-no-basic-feature if you want to perform enrichment analysis on focused annotations only (-\-more-bed or -\-more-key).

 -- The goal of the minibatch is to save RAM. Increase the number of minibatches, instead of their size.
 You may need to use very small minibatches if you have large sets of regions.

 -- You can exclude regions from the shuffling. This is done by shuffling across a concatenated "sub-genome" obtained by removing
 the excluded regions, but the same ones will be excluded from the peak_file and the GTF.
 This in Beta for now and will be very time-consuming (hours), especially if you have few CPU cores.
 Try using an exclusion file that is as small (around a thousand elements) as possible.

 -- BETA : About -\-use-markov. This arguments control whether to use Markov model realisations instead of independant shuffles
 for respectively region lengths and inter-region lengths. This is not recommended in the general case and can *very* time-consuming (hours).

 """


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    # --------------------- Main arguments ----------------------------------- #

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Defaults to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')),
                            required=False)

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

    # --------------------- More regions ------------------------------------- #

    parser_grp.add_argument('-b', '--more-bed',
                            help="A list of bed files to be considered as additional genomic annotations.",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'),
                            nargs='*',
                            required=False)

    parser_grp.add_argument('-l', '--more-bed-labels',
                            help="A comma separated list of labels (see --more-bed)",
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-e', '--bed-excl',
                            help='Exclusion file. The chromosomes will be shortened by this much for the shuffles of peaks and features. Can take a long time.'
                                 ' (bed format).',
                            default=None,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'),
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

    parser_grp.add_argument('-m', '--more-keys',
                            help='A comma separated list of key used for labeling the genome. See Notes.',
                            type=str,
                            default=None,
                            required=False)

    parser_grp.add_argument('-n', '--no-basic-feature',
                            help="No statistics for basic features of GTF. Concentrates on --more-bed and --more-keys.",
                            action="store_true",
                            required=False)

    # --------------------- Backend ------------------------------------------ #

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

    parser_grp.add_argument('-ma', '--use-markov',
                            help='Whether to use Markov model realisations instead of independant shuffles. See notes.',
                            action='store_true',
                            required=False)

    # --------------------- Output ------------------------------------------- #

    parser_grp.add_argument('-o', '--outputdir',
                            help='Output directory name.',
                            metavar="DIR",
                            default="ologram_output",
                            type=str)

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

    parser_grp.add_argument('-if', '--user-img-file',
                            help="Provide an alternative path for the main image. ",
                            default=None,
                            nargs=None,
                            type=arg_formatter.FormattedFile(mode='w', file_ext='pdf'),
                            required=False)

    parser_grp.add_argument('-x', '--no-pdf',
                            help="Do not produce any image file. ",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-tp', '--tsv-file-path',
                            help="Provide an alternative path for text output file.",
                            default=None,
                            type=arg_formatter.FormattedFile(mode='w', file_ext='txt'),
                            required=False)

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=arg_formatter.ranged_num(0, None),
                            default=300,
                            required=False)

    # --------------------- Other input arguments----------------------------- #

    parser_grp.add_argument('-z', '--no-gtf',
                            help="No gtf file is provide as input.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-f', '--force-chrom-gtf',
                            help="Discard silently, from GTF, genes outside chromosomes defined in --chrom-info.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-w', '--force-chrom-peak',
                            help="Discard silently, from --peak-file, peaks outside chromosomes defined in --chrom-info.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-q', '--force-chrom-more-bed',
                            help="Discard silently, from --more-bed files, regions outside chromosomes defined in --chrom-info.",
                            action='store_true',
                            required=False)
    return parser


# -------------------------------------------------------------------------
# The command function
# -------------------------------------------------------------------------


def ologram(inputfile=None,
            outputdir=None,
            peak_file=None,
            chrom_info=None,
            tsv_file_path=None,
            more_bed=None,
            more_bed_labels=None,
            no_gtf=False,
            upstream=1000,
            more_keys=None,
            downstream=1000,
            no_basic_feature=False,
            bed_excl=None,
            use_markov=False,
            no_pdf=None,
            pdf_width=5,
            pdf_height=5,
            force_chrom_gtf=False,
            force_chrom_peak=False,
            force_chrom_more_bed=False,
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

    # -------------------------------------------------------------------------
    # Set random seed
    # -------------------------------------------------------------------------

    np.random.seed(seed)

    # -------------------------------------------------------------------------
    # Are we using markov shuffling ?
    # If yes, send a warning to the user.
    # -------------------------------------------------------------------------

    if use_markov:
        message('Using Markov order 2 shuffling.', type='INFO')
        message(
            'Markov-based null is still in beta at the moment and tends to biais the "null" hypothesis towards association.',
            type='WARNING')

    # -------------------------------------------------------------------------
    # If user wants don't provide a GTF
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
    # Load the peak file as pybedtools.BedTool object
    # -------------------------------------------------------------------------

    peak_file = pybedtools.BedTool(peak_file.name)

    # -------------------------------------------------------------------------
    # Check chromosomes for peaks are defined in the chrom-info file
    # Depending on force_chrom_peak, peaks undefined in peak_file may
    # be silently removed.
    # -------------------------------------------------------------------------

    peak_chrom_list = set()

    for i in pybedtools.BedTool(peak_file):
        peak_chrom_list.add(i.chrom)

    if not force_chrom_peak:
        for i in peak_chrom_list:
            if i not in chrom_len:
                msg = "Chromosome " + str(i) + " from peak file is undefined in --chrom-info file. "
                message(msg + 'Please fix --chrom-info file or use --force-chrom-peak.',
                        type="ERROR")
    else:
        peak_file_sub = make_tmp_file(prefix='peaks_x_chrom_info', suffix='.bed')

        for i in peak_file:
            if i.chrom in chrom_len:
                peak_file_sub.write("\t".join(i.fields) + "\n")

        peak_file_sub.close()
        peak_file = BedTool(peak_file_sub.name)

    # -------------------------------------------------------------------------
    # Sort and merge the peaks
    # -------------------------------------------------------------------------
    # Just in case it was not, sort and merge the file.
    # In any case, it should be short compared to the
    # expected total running time.
    peak_file = peak_file.sort().merge()

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
    # rest of the ologram code relie otherwise only on its keys.

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

    if not no_gtf:
        if not no_basic_feature or more_keys:
            gtf = GTF(inputfile)
            gtf_chrom_list = gtf.get_chroms(nr=True)

            # -------------------------------------------------------------------------
            # Check chromosomes from the GTF are defined in the chrom-info file
            # -------------------------------------------------------------------------

            if not force_chrom_gtf:
                for i in gtf_chrom_list:
                    if i not in chrom_len:
                        msg = "Chromosome " + str(i) + " from GTF is undefined in --chrom-info file. "
                        message(msg + "Please check your --chrom-info file or use --force-chrom-gtf",
                                type="ERROR")

            # -------------------------------------------------------------------------
            # Subset the GTF using chromosomes defined in chrom-info file.
            # -------------------------------------------------------------------------

            gtf = gtf.select_by_key("seqid", ",".join(chrom_len.keys()))

            if len(gtf) == 0:
                message("The GTF file does not contain any genomic feature "
                        "falling in chromosomes declared in --chrom-info.",
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
                                         alist=["00_ologram_stats.tsv",
                                                "00_ologram_diagrams.pdf"
                                                ],
                                         force=True)

    data_file, pdf_file = file_out_list

    if no_pdf:
        if user_img_file:
            os.unlink(user_img_file.name)
        os.unlink(pdf_file.name)
        pdf_file = None
    else:
        if user_img_file is not None:

            os.unlink(pdf_file.name)
            pdf_file = user_img_file

            test_path = os.path.abspath(pdf_file.name)
            test_path = os.path.dirname(test_path)

            if not os.path.exists(test_path):
                os.makedirs(test_path)

    if tsv_file_path is not None:

        os.unlink(data_file.name)
        data_file = tsv_file_path

        test_path = os.path.abspath(data_file.name)
        test_path = os.path.dirname(test_path)

        if not os.path.exists(test_path):
            os.makedirs(test_path)

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

    if not no_gtf:
        if not no_basic_feature:
            for feat_type in gtf.get_feature_list(nr=True):
                message("Processing " + str(feat_type), type="INFO")
                gtf_sub = gtf.select_by_key("feature", feat_type, 0)
                gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                                   "gene_id",
                                                   "exon_id"]).sort().merge()  # merging bed file !
                tmp_file = make_tmp_file(prefix=str(feat_type), suffix='.bed')
                gtf_sub_bed.saveas(tmp_file.name)

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

            if not force_chrom_more_bed:
                chrom_list = set()
                for i in BedTool(bed_anno.name):
                    chrom_list.add(i.chrom)

                for i in chrom_list:
                    if i not in chrom_len:
                        message("Chromosome " + str(i) + " is undefined in --more-bed with label " + bed_lab + ".",
                                type="ERROR")
            else:
                bed_anno_sub = make_tmp_file(prefix='more_bed_x_chrom_info' + bed_lab, suffix='.bed')

                n = 0
                for i in BedTool(bed_anno.name):
                    if i.chrom in chrom_len:
                        bed_anno_sub.write("\t".join(i.fields) + "\n")
                        n += 1
                if n == 0:
                    msg = "The --more-bed file " + bed_lab + " is empty after checking for --chrom-info."
                    message(msg + "Please check your --chrom-info file or use --force-chrom-more-bed",
                            type="ERROR")

                bed_anno_sub.close()
                bed_anno = bed_anno_sub

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

    if pdf_file is not None:
        plot_results(d, data_file, pdf_file, pdf_width, pdf_height, dpi)
        close_properly(pdf_file)
    close_properly(data_file)


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

    # -------------------------------------------------------------------------
    # Rename the feature type.
    # When --more-keys is used the key and value are separated by ":".
    # This give rise to long name whose display in the plot is ugly.
    # We can break these names using a "\n".
    # -------------------------------------------------------------------------

    dm["feature_type"] = [x.replace(":", "\n") for x in dm["feature_type"]]

    # -------------------------------------------------------------------------
    # Create a new plot
    # -------------------------------------------------------------------------

    message('Adding bar plot.')

    # This can be used to plot either 'summed_bp_overlaps' or 'nb_intersections'
    def plot_this(statname):

        ## DATA PROCESSING

        # Collect true and shuffled number of the stat being plotted
        data_ni = dm[['feature_type', statname + '_esperance_shuffled', statname + '_true']]
        maximum = data_ni[[statname + '_esperance_shuffled', statname + '_true']].max(axis=1)

        data_ni.columns = ['Feature', 'Shuffled', 'True']  # Rename columns

        # For later purposes (p-value display), collect the fold change.
        fc = data_ni['True'] / (data_ni['Shuffled'] + 1)

        # Now melt the dataframe
        dmm = data_ni.melt(id_vars='Feature')
        dmm.columns = ['Feature', 'Type', statname]

        ## PLOTTING

        # Create plot
        p = ggplot(dmm)
        p += theme_bw()  # Add the black & white theme

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
                           position=position_dodge(.6), size=0.3)

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

        # Compute the colors for the text box : orange if significantly depleted,
        # green if significantly enriched, black otherwise. For display purposes,
        # p<0.05 counts as significant.
        signif_color = pd.Series(['#000000'] * len(text))
        for i in range(len(text)):
            if text[i] < 0.05:  # If significant
                if fc[i] < 1: signif_color[i] = '#f57c00'
                if fc[i] > 1: signif_color[i] = '#43a047'

        text = text.apply(format_pvalue)
        text_pos = (maximum + 0.05 * max(maximum)).append(na_series)
        text_pos.index = range(len(text_pos))
        aes_plot = aes(x='Feature', y=text_pos, label=text)
        p += geom_label(mapping=aes_plot, stat='identity',
                        size=5, boxstyle='round', label_size=0.2,
                        color='white', fill=signif_color)

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

        # Add a nicer set of colors.
        p += scale_fill_manual(values={'Shuffled': '#757575', 'True': '#0288d1'})

        return p

    # Compute the plots for both statistics
    p1 = plot_this('summed_bp_overlaps') + ylab("Nb. of overlapping base pairs") + ggtitle(
        'Total overlap length per region type')
    p2 = plot_this('nb_intersections') + ylab("Number of intersections") + ggtitle(
        'Total nb. of intersections per region type')

    # -------------------------------------------------------------------------
    # Computing page size
    # -------------------------------------------------------------------------

    nb_ft = len(list(d['feature_type'].unique()))

    if pdf_width is None:
        panel_width = 0.6
        pdf_width = panel_width * nb_ft

        if pdf_width > 25:
            pdf_width = 25
            message("Setting --pdf-width to 25 (limit)")

    if pdf_height is None:
        pdf_height = 5

    message("Page width set to " + str(pdf_width))
    message("Page height set to " + str(pdf_height))
    figsize = (pdf_width, pdf_height)

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

        # NOTE : We must manually specify figure size with save_as_pdf_pages
        save_as_pdf_pages(filename=pdf_file.name,
                          plots=[p1 + theme(figure_size=figsize), p2 + theme(figure_size=figsize)],
                          width=pdf_width,
                          height=pdf_height,
                          dpi=dpi)


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    ologram(**args)


if __name__ == '__main__':
    main()


else:

    # 'Bats' tests
    test = '''
        #ologram: get example files
        @test "ologram_0" {
             result=`gtftk get_example -d simple_02 -f '*'`
          [ "$result" = "" ]
        }

        #ologram: run on simple test file
        @test "ologram_1" {
             result=`rm -Rf ologram_output; gtftk ologram -i simple_02.gtf -p simple_02_peaks.bed -c simple_02.chromInfo -u 2 -d 2 -K ologram_output --no-date`
          [ "$result" = "" ]
        }

        #ologram: proper number of true intersections
        @test "ologram_2" {
         result=`cat ologram_output/00_ologram_stats.tsv | grep gene | cut -f 6`
          [ "$result" = "16" ]
        }

        #ologram: proper number of shuffled intersections
        @test "ologram_3" {
         result=`cat ologram_output/00_ologram_stats.tsv | grep gene | cut -f 2`
          [ "$result" = "14.97" ]
        }

        #ologram: overlapping bp
        @test "ologram_4" {
         result=`cat ologram_output/00_ologram_stats.tsv | grep gene | cut -f 12`
          [ "$result" = "75" ]
        }

        #ologram: shuffled overlapping bp
        @test "ologram_5" {
         result=`cat ologram_output/00_ologram_stats.tsv | grep gene | cut -f 8`
          [ "$result" = "61.35" ]
        }

        #ologram: shuffled overlapping bp variance
        @test "ologram_6" {
         result=`cat ologram_output/00_ologram_stats.tsv | grep gene | cut -f 9`
          [ "$result" = "32.94" ]
        }

        #ologram: shuffled overlapping bp fitting
        @test "ologram_7" {
         result=`cat ologram_output/00_ologram_stats.tsv | grep gene | cut -f 10`
          [ "$result" = "0.8227700000000001" ]
        }
        '''

    cmd = CmdObject(name="ologram",
                    message="Statistics on bed file intersections with genomic features.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="annotation",
                    notes=__notes__,
                    updated=__updated__,
                    test=test)
