#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import re
import sys
from collections import defaultdict

from pybedtools import BedTool

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import bedFileList, bedFile
from pygtftk.arg_formatter import checkChromFile
from pygtftk.arg_formatter import int_greater_than_null_or_None
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import check_r_installed
from pygtftk.utils import check_r_packages
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
 
 -- If -\more-keys is used additional region sets will be tested based on the associated key value.
 As an example, if -\more-keys is set to the 'gene_biotype' (a key generally found in ensembl GTF), the
 region related to 'protein_coding', 'lncRNA' or any other value for that key will be retrieved merged and tested.

 -- Use -\no-basic-feature if you want to perform enrichment analysis on focused annotations only (-\more-bed or -\more-key).

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
                            type=bedFile(),
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


# Perform intersection between a bedTools object
# an another file (feature_file) e.g promoter, intron, exon...
# Fill a dictionary of dictionary


def _intersection_results(reg_file=None,
                          feature_bo=None,
                          my_dict=None,
                          chrom_len=None,
                          page_format=None,
                          user_img_file=None,
                          ft_type=None,
                          file_out_list=None):
    """Get feature and peaks and compute intersections. Returns an updated
    dict."""

    # Compute the number of regions of interest falling in
    # each chromosome.

    feature_bo = feature_bo.sort()
    reg_file = BedTool(reg_file)
    nb_peaks = len(reg_file)

    for i in reg_file:
        chrom = i.chrom
        if chrom in chrom_len:
            if ft_type != "Full_chromosomes":
                my_dict[ft_type][chrom]["Nb_trial_or_peak"] += 1
                my_dict[ft_type]["all_chrom"]["Nb_trial_or_peak"] += 1
            else:
                my_dict[ft_type][chrom]["Nb_trial_or_peak"] = nb_peaks

        else:
            for file_out in file_out_list:
                os.remove(file_out.name)
            msg = "Peak file: " + chrom + \
                  " is undefined in chromInfo file. Please fix."
            message(msg, type="ERROR")

    # Compute the nucleotide size of all feature (promoter, exons, introns,...)

    feature_merge_bo = feature_bo.merge()

    # Save features (introns, intergenic,...) for tracability
    peak_anno_tmp_file = make_tmp_file(prefix="peak_anno_" + ft_type + "_merged",
                                       suffix=".bed")
    feature_merge_bo.saveas(peak_anno_tmp_file.name)
    peak_anno_tmp_file.close()

    for i in feature_merge_bo:

        chrom = i.chrom

        if chrom in chrom_len:
            size = i.end - i.start
            if ft_type != "Full_chromosomes":
                my_dict[ft_type][chrom]["coverage"] += size
                my_dict[ft_type]["all_chrom"]["coverage"] += size
            else:
                # We will compare "Chromosomes" to the full genome.
                my_dict[ft_type][chrom][
                    "coverage"] = chrom_len[chrom]
        else:
            for file_out in file_out_list:
                os.remove(file_out.name)
            msg = """The chromosome {chrom} is undefined in: {file}. Please fix."""
            msg = msg.format(chrom=chrom, file=feature_bo.fn)
            message(msg, type="ERROR")

    # Compute intersections
    intersections = reg_file.intersect(feature_merge_bo)

    for i in intersections:

        if i.chrom in chrom_len:

            chrom = i.chrom
            my_dict[ft_type][chrom]["Observed"] += 1
            if ft_type != "Full_chromosomes":
                my_dict[ft_type]["all_chrom"]["Observed"] += 1
        else:
            for file_out in file_out_list:
                os.remove(file_out.name)
            msg = "The chromosome " + chrom + " is undefined."
            message(msg, type="WARNING")

    file_out_save = make_tmp_file(prefix="peak_anno_intersections_" +
                                         ft_type,
                                  suffix=".bed")
    intersections.saveas(file_out_save.name)

    return my_dict


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
    classical features (e.g promoter, tts, gene body, UTR,...) and to sets of user provided pics.
    """

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

    # chrom sizes
    chrom_len = chrom_info_as_dict(chrom_info)

    if len(chrom_len.keys()) > 100:
        message("Can't accept more than 100 chromosomes in peak_anno (see --chrom-info).",
                type="ERROR")

    # A muti-level dict to store the results
    # per chromosome
    # e.g my_dict[ft_type][chrom]["Observed"]
    # e.g my_dict[ft_type][chrom]["coverage"]
    # e.g my_dict[ft_type][chrom]["Nb_trial_or_peaks"]

    def nested_dict(n, type):
        """"http://stackoverflow.com/questions/29348345"""
        if n == 1:
            return defaultdict(type)
        else:
            return defaultdict(lambda: nested_dict(n - 1, type))

    hits = nested_dict(3, int)

    # Read the gtf file
    # Discard any records corresponding to chr not declared in
    # ChromInfo file
    if not no_basic_feature or more_keys:
        gtf = GTF(inputfile).select_by_key("seqid", ",".join(chrom_len.keys()))

        if len(gtf) == 0:
            message("The GTF file does not contain any genomic feature "
                    "falling in chromosomes declared in chromInfo file.",
                    type="ERROR")

    # Check user provided annotations
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

    # Preparing output files
    file_out_list = make_outdir_and_file(out_dir=outputdir,
                                         alist=["00_peak_anno_stats.txt",
                                                "00_peak_anno_diagrams." + page_format,
                                                "00_peak_anno_diagrams_by_chrom." + page_format,
                                                "00_peak_anno_R_code.R"],
                                         force=True)

    data_file, pdf_file, pdf_file_by_chrom, r_code_file = file_out_list

    # Get the midpoints of the peaks
    region_mid_point = make_tmp_file("peaks_midpoints", ".bed")

    # Loop through peaks
    peak_bo = BedTool(peak_file.name)

    for line in peak_bo:

        diff = line.end - line.start

        if diff % 2 != 0:
            # e.g 10-13 (zero based) -> 11-13 one based
            # mipoint is 12 (one-based) -> 11-12 (zero based)
            # e.g 949-1100 (zero based) -> 950-1100 one based
            # mipoint is 1025 (one-based) -> 1024-1025 (zero based)
            # floored division (python 2)...
            line.end = line.start + int(diff / 2) + 1
            line.start = line.end - 1
        else:
            # e.g 10-14 (zero based) -> 11-14 one based
            # mipoint is 12-13 (one-based) -> 11-13 (zero based)
            # e.g 9-5100 (zero based) -> 10-5100 one based
            # mipoint is 2555-2555 (one-based) -> 2554-2555 (zero based)
            # floored division (python 2)...
            # No real center. Take both

            line.start = line.start + int(diff / 2) - 1
            line.end = line.start + 2

        region_mid_point.write(str(line))

    region_mid_point.close()

    # Fill the dict with info about basic features include in GTF
    if not no_basic_feature:
        for feat_type in gtf.get_feature_list(nr=True):
            if feat_type not in ["start_codon", "stop_codon"]:
                gtf_sub = gtf.select_by_key("feature", feat_type, 0)

                gtf_sub_bed = gtf_sub.to_bed(name=["transcript_id",
                                                   "gene_id",
                                                   "exon_id"]).sort().merge()  # merging bed file !

                hits = _intersection_results(reg_file=region_mid_point.name,
                                             feature_bo=gtf_sub_bed,
                                             my_dict=hits,
                                             chrom_len=chrom_len,
                                             ft_type=feat_type,
                                             file_out_list=file_out_list)

        # Get the intergenic regions

        gtf_sub_bed = gtf.get_intergenic(chrom_info,
                                         0,
                                         0,
                                         chrom_len.keys()).merge()

        hits = _intersection_results(reg_file=region_mid_point.name,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Intergenic")

        # Get the intronic regions
        gtf_sub_bed = gtf.get_introns()

        hits = _intersection_results(reg_file=region_mid_point.name,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Introns")

        # Get the promoter regions
        gtf_sub_bed = gtf.get_tss().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits = _intersection_results(reg_file=region_mid_point.name,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Promoters")

        # Get the tts regions
        gtf_sub_bed = gtf.get_tts().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        hits = _intersection_results(reg_file=region_mid_point.name,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="TTS")

        # Test the whole chromosomes as features
        gtf_sub_bed = gtf.get_tts().slop(s=True,
                                         l=upstream,
                                         r=downstream,
                                         g=chrom_info.name).cut([0, 1, 2,
                                                                 3, 4, 5]).sort().merge()

        gtf_sub_bed = BedTool(
            chrom_info_to_bed_file(
                chrom_info,
                chr_list=chrom_len.keys()))

        hits = _intersection_results(reg_file=region_mid_point.name,
                                     feature_bo=gtf_sub_bed,
                                     my_dict=hits,
                                     chrom_len=chrom_len,
                                     ft_type="Full_chromosomes")

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

                    hits = _intersection_results(reg_file=region_mid_point.name,
                                                 feature_bo=gtf_sub_bed,
                                                 my_dict=hits,
                                                 chrom_len=chrom_len,
                                                 ft_type=":".join([user_key,
                                                                   val]))
    # Process user defined annotations
    if more_bed is not None:
        message("Processing user-defined regions (bed format).")
        for bed_anno, bed_lab in zip(more_bed, more_bed_labels):
            hits = _intersection_results(reg_file=region_mid_point.name,
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

    nb_peaks = str(len(BedTool(region_mid_point.name)))

    for chrom in chrom_len:
        if chrom != "all_chrom":
            hits["Full_chromosomes"][chrom]["Nb_trial_or_peak"] = nb_peaks

    for key1 in hits:
        for key2 in chrom_len:
            if key1 == "Full_chromosomes":
                ref_size = str(chrom_len["all_chrom"])
                nb_trial = str(len(BedTool(region_mid_point.name)))
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

    code_body = """

    # Gtftk-like messages
    message <- function(msg){{
      cat(paste("    |--- ",format(Sys.time(),"%R"), "-INFO : ", msg, "\\n", sep=""), file=stderr())
    }}

    None <- 'None'

    ####################################################
    # Load libraries
    ####################################################
    
    library(reshape2)
    library(ggplot2)
    library(scales) # Used for declaration of pseudolog10 axis transformation "trans_new".

    ####################################################
    # Read the dataset
    ####################################################
     
    d <- read.table("{data_file}",
                    head=T,
                    sep="\\t",
                    quote="")

    ####################################################
    # Adjust label names if all features comes from the same key
    ####################################################
    # When using only features from a given key in a gtf,
    # the readability can be increased removing the
    # redundant 'key' used as prefix for all features.

    # Set default xlab for plot:
    xlabel <- 'Features'

    # Check that 'ft_type' either contains ':' the separator between key and value or the 'Full_chromosomes' feature.
    if (length(grep(":|Full_chromosomes",d[,"ft_type"], invert=FALSE)) == length(d[,"ft_type"])){{
        potential_common_key_rows=grep(":",d[,"ft_type"])
        potential_common_key_vector=gsub(":.*$", '', d[potential_common_key_rows,"ft_type"])

        if (length(unique(potential_common_key_vector)) == 1){{
            common_key <- unique(potential_common_key_vector)

            #Remove 'common_key' from 'ft_type'
            # 'as.factor' needed for 'by_chrom' plot to be correctly drawn.
            d[,"ft_type"] <- as.factor(gsub(paste0(common_key,":"), '', d[,"ft_type"]))

            xlabel <- paste("Features in",common_key)

        }}
    }}

    ####################################################
    # Get the chromosome order
    ####################################################

    chr_order <- d$chrom
    chr_order <- sort(chr_order[!duplicated(chr_order)], decreasing = TRUE)
    chr_order_num <- gsub("[^0-9]", "", chr_order)
    chr_order <- chr_order[order(as.numeric(chr_order_num))]
    d$chrom <- factor(d$chrom, levels = chr_order)
    
    # Compute expected number of intersections
    d$freq <- d$coverage / d$reference_size

    d$Expected <- d$freq * d$Nb_trial_or_peak

    ####################################################
    # Compute binomial p.val  (unilateral)
    ####################################################
    
    for(i in 1:nrow(d)){{
      
          if(d$Expected[i] > 0){{
              if(d$Observed[i] > 0){{
                  log2_ratio <- log2(d$Observed[i] / d$Expected[i])
              }}else{{
                  log2_ratio <- NA
              }}
          }}else{{
              log2_ratio <- NA
          }}
    
          d$log2_ratio[i] <- log2_ratio
        
          if(d$Nb_trial_or_peak[i] > 0){{
          pval <- binom.test(alternative="two.sided",
                              n=d$Nb_trial_or_peak[i],
                              x= d$Observed[i],
                              p = d$freq[i])$p.value
          }}else{{
              pval <- NA
          }}
          
          d$pval.binom[i] <- pval
          
          if(d$Observed[i] > d$Expected[i]){{
              d$test_type[i] <- "Enrichment"
          }}else if(d$Observed[i] < d$Expected[i]){{
              d$test_type[i] <- "Depletion"
          }}else{{
                d$test_type[i] <- "Unchanged"
          }}
      
    }}

    ####################################################
    # Melt the data frame
    ####################################################

    dm <- melt(d[ , c("ft_type", "chrom", "Observed", "Expected")],
                 id.vars =c("ft_type", "chrom"))
    
    dm$variable <- factor(dm$variable, levels=c("Expected", "Observed"))

    ####################################################
    # function to compute pseudolog10 transformation of y-axis
    ####################################################
  
    log10p_trans <- function() {{trans_new(
        name="log10p",
        transform=function(x){{log10(x+1)}},
        inverse=function(x){{10^x - 1}},
        breaks = trans_breaks('log10', function(x) 10^x),
        format = trans_format('log10', math_format(10^.x))
        )}}

    ####################################################
    # function to compute diagram
    ####################################################

    bar_plot_with_pval <- function(ft_type=NULL,
                                    which_chrom=NULL,
                                    which_row=NULL,
                                    by_chrom=FALSE,
                                    y_axis_trans=log10p_trans)
                                    {{


          p <- ggplot()
          if(! by_chrom){{
            aes.plot <- aes(ft_type, value, fill=variable)
          }}else{{
            aes.plot <- aes(chrom, value, fill=variable)
          }}
            p <- p + geom_bar(data=dm[which_row, ],
                        aes.plot,
                        stat="identity",
                        position = position_dodge(width = 0.90),
                        alpha = 0.6,
                        show.legend=TRUE)

        p <- p + ylab("Number of overlaps")
        p <- p + xlab(xlab.plot)
        p <- p + theme(legend.title=element_blank())
        p <- p + theme_bw()
        p <- p + theme(legend.position="bottom")
        p <- p + theme(legend.key = element_blank())
        p <- p + theme(legend.title=element_blank())
        p <- p + theme(axis.title.x = element_text(colour="grey20",
                                                    size=10,angle=0,
                                                    face="plain"))
        p <- p + theme(axis.text.x = element_blank())
        p <- p + theme(axis.text.x = element_text(size=10,
                                                  vjust = 0.9,
                                                  hjust = 1,
                                                  angle=45))
        p <- p + theme(axis.text.y = element_blank())
        p <- p + theme(axis.text.y = element_text(size=8,
                                                  angle=0))

        if(! by_chrom){{
          aes.plot <- aes(label = sprintf("%.02f", value),
                          x = ft_type,
                          y = ifelse(value < 1, 1, value),
                          colour=variable)
        }}else{{
          aes.plot <- aes(label = sprintf("%.02f", value),
                          x = chrom,
                          y = ifelse(value < 1, 1, value),
                          colour=variable)
        }}
        
        p <- p + geom_text(data=dm[which_row, ],
                           aes.plot,
                          position = position_dodge(width = 0.8),
                          angle=90,
                          vjust = 0.5,
                          hjust=-0.25,
                          size=2)
    
        p <- p + guides(colour=FALSE, fill=guide_legend(ncol=2)        )

    
        
        ####################################################
        # Add an horizontal segment to indicate p-val
        ####################################################
        
        j <- 1
    
        # keep the highest value (exp/obs)
        # for each genomic feature (prom, UTR,...)
        dm.sub <- dm[which_row, ]
        
        if(! by_chrom){{
          dm.sub <- do.call(rbind,
                        lapply(
                          lapply(split(dm.sub, dm.sub$ft_type),
                                 function(x) x[order(x$value, decreasing = TRUE), ]),
                          head, 1))
        }}else{{
          dm.sub <- do.call(rbind,
                            lapply(
                              lapply(split(dm.sub, dm.sub$chrom),
                                     function(x) x[order(x$value, decreasing = TRUE), ]),
                              head, 1))
        }}
        
        if(! by_chrom){{
          ft_type_list <- levels(as.factor(as.character(dm.sub$ft_type)))
        }}else{{
          ft_type_list <- levels(d$chrom)[levels(d$chrom) %in% levels(as.factor(as.character(dm.sub$chrom)))]
        }}
        
        for(i in 1:length(ft_type_list)){{
 
              y_val <- ifelse(is.na(dm.sub$value[i]), 8, dm.sub$value[i] * 8)
              y_val <- ifelse(y_val < 8 , 8, y_val)
              x <- c(i-0.2, i-0.2, i+0.2, i+0.2)
              y <- c(y_val * 0.95,
                     y_val * 1,
                     y_val * 1,
                     y_val * 0.95)
                        
              z <- ft_type_list[i]
              p <- p + geom_path(data = data.frame(x=x,
                                                   y=y,
                                                   ft_type=z,
                                                   value=NA,
                                                   variable=NA,
                                                   check.names = FALSE),
                                 aes(x = x, y = y)
                                 )
    
        }}
    
    
        ####################################################
        # Add p-val
        ####################################################
            
        for(i in 1:length(ft_type_list)){{
         
             x <- i
             y <- ifelse(is.na(dm.sub$value[i]), 10, dm.sub$value[i] * 12)
             y <- ifelse(y < 10, 10, y)


             if(!by_chrom){{
                    selected <- d$chrom == which_chrom & d$ft_type == ft_type_list[i]
                    lab <- paste(sprintf(" %.03g",
                                         d[selected, ]$pval.binom),
                                         "\\n",
                                         sprintf("%.03g ", d[selected, ]$log2_ratio))
               
                   test_type <- d[selected, ]$test_type
                   d.tmp <- data.frame(x=x,
                                       y=y,
                                       value=lab,
                                       ft_type=ft_type_list[i],
                                       test_type=test_type,
                                       check.names = FALSE)
                                       
                  d.tmp$test_type <- as.factor(d.tmp$test_type)
              }}else{{
                selected <- d$chrom == ft_type_list[i] & d$ft_type == ft_type
                lab <- paste(sprintf(" %.03g",
                                     d[selected, ]$pval.binom),
                                     "\\n",
                                     sprintf("%.03g ", d[selected, ]$log2_ratio))
                
                test_type <- as.character(d[selected, ]$test_type)
                d.tmp <- data.frame(x=x,
                                    y=y,
                                    value=lab,
                                    ft_type=ft_type,
                                    test_type=test_type,
                                    check.names = FALSE)
                
                d.tmp$test_type <- as.factor(d.tmp$test_type)
              }}

             p <- p + geom_text(angle = 90,
                                data=d.tmp,
                                aes(label = value,
                                    x = x,
                                    y = y,
                                    colour=factor(test_type)),
                                    fontface = "bold",
                                position = position_dodge(width = 0.8),
                                vjust = 0.5,
                                hjust=0.5,
                                size=2,
                                show.legend=F)
            #print(d.tmp)
            
        }}
        
    
        ####################################################
        # plot limits
        ####################################################
    
        ylim_max <- max(dm$value) * 25
        ybreaks_max <- as.numeric(format(ylim_max,digits=1))
        p <- p + scale_y_continuous(trans="log10p",breaks=trans_breaks('log10', function(x) 10^x)(c(1,ybreaks_max)))

        ####################################################
        # Colors
        ####################################################
            
        p <- p + scale_fill_manual(values=c("Observed"="blue",
                                            "Expected"="#8A8A8A",
                                            "Enrichment"="#8B0000",
                                            "Unchanged"="#000000",
                                            "Depletion"="#006400"))
    
        p <- p + scale_color_manual(values=c("Observed"="blue",
                                            "Expected"="#8A8A8A",
                                            "Enrichment"="#8B0000",
                                            "Unchanged"="#000000",
                                            "Depletion"="#006400"))
        
        return(p)
        
    }}


    ####################################################
    # Compute bar plot by feature
    ####################################################


                
    xlab.plot <- xlabel
    ft_type <- levels(d$ft_type)
    which_chrom = "all_chrom"
    which_row <- dm$chrom == which_chrom
    pdf_path <- "{pdf_file}"

    pdf_width <- {pdf_width}
    pdf_height <- {pdf_height}
        
    if(pdf_width == 'None')
        pdf_width <- length(ft_type)

    if(pdf_height == 'None')
        pdf_height <- 6
        
    
    ggsave(filename=pdf_path,
         plot=bar_plot_with_pval(ft_type=ft_type,
                               which_chrom=which_chrom,
                               which_row=which_row),
         width=pdf_width,
         height=pdf_height)
             


    
    ####################################################
    # Compute bar plot ~ chromosome
    ####################################################


    pdf_path <- "{pdf_file_by_chrom}"
    {page_format}(pdf_path, width=pdf_width, height=pdf_height)
    by_chrom <- TRUE
    
    plot_list <- list()
    

    for(i in levels(d$ft_type)){{

        message(paste("Preparing profile diagram by chromosome: ", i))
        xlab.plot <- i
        ft_type <- i
        which_row <- dm$chrom != "all_chrom" & dm$ft_type == i
        
        p <- bar_plot_with_pval(ft_type=ft_type,
                                which_chrom=which_chrom,
                                which_row=which_row,
                                by_chrom=by_chrom)

        print(p)

    }}
    
    dev.null <- dev.off()

    
    write.table(d, "{data_file}", sep="\\t", quote=F, col.names=TRUE, row.names=FALSE)
    """.format(data_file=data_file.name,
               pdf_file=pdf_file.name,
               pdf_file_by_chrom=pdf_file_by_chrom.name,
               pdf_width=pdf_width,
               pdf_height=pdf_height,
               page_format=page_format)

    if verbosity:
        message("Printing R code to: " + r_code_file.name)

    print(code_body, file=r_code_file)
    r_code_file.close()

    if verbosity:
        message("Executing R code.")

    # Execute R code.
    check_r_packages(["reshape2", "ggplot2"])
    check_r_installed()

    os.system("cat " + r_code_file.name + "| R --slave")


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
                    fun=peak_anno,
                    desc=__doc__,
                    group="annotation",
                    notes=__notes__,
                    updated=__updated__,
                    test=test,
                    rlib=R_LIBS)
