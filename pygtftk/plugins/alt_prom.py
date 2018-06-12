#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys
from _collections import defaultdict

from pybedtools import BedTool

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import bedFileWithUnambiguousNames
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly, message
from pygtftk.utils import flatten_list
from pygtftk.utils import make_tmp_file

__updated__ = "2018-02-15"
__doc__ = """
Takes a GTF as input and two H3K4me3 bed files. Search for genes with alternative promoters.
 """
__notes__ = """
-- A set of two BED for H3K4me3 should be provided and the program will try to find promoters which are active in one sample
while inactive in the other sample.

-- To call a gene has having an alternative promoter two transcript TSSs and two H3K4me3 bed files from two differents
samples are considered. The two following configurations are called 'alternative promoters'. The first configuration is 
considered as a gain or loss of promoter (depending on the reference condition).The second one is called 
'promoter switch'. These information are contained in two dedicated columns in the output file.

--   
-- |  
-- | Case 1:
-- |              |-\-\->       |-\-\->    
-- |              |             |
-- | -\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
-- |
-- | bed_a     -\-\-\        -\-\-\ 
-- | bed_b     -\-\-\  
-- | 
-- | Case 2:
-- |              |-\-\->       |-\-\->    
-- |              |             |
-- | -\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
-- | bed_a                   -\-\-\ 
-- | bed_b     -\-\-\  
-- |  
-- |  
--

-- Only pair of TSSs whose distance is greater to -\dist-min are tested. By default -\dist-min is set to zero.

-- Peaks should probably be merged before with mergeBed (Bedtools) to force merging of close peaks.

-- The gene_nb_computed_promoters columns contains, for a given gene_id the number of features, merged from both H3K4me3 conditions, that overlap its TSSs. It is supposed to give an estimation of genes that may potentially undergo alternative promoter gain, loss or switch.  

"""


# -------------------------------------------------------------------------
# The main parser
# -------------------------------------------------------------------------


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)
    parser_grp = parser.add_argument_group('Arguments')

    # required but not required...

    parser_grp.add_argument(
        'bed_list',
        help='A list of bed files for H3K4me3.',
        type=bedFileWithUnambiguousNames(),
        nargs='+')

    parser_grp.add_argument('-i', '--inputfile',
                            help="The input GTF file.",
                            default=sys.stdin,
                            metavar="GTF",
                            required=False,
                            type=FileWithExtension('r',
                                                   valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Cc][Oo][Vv]')))

    parser_grp.add_argument('-l', '--labels',
                            help='Bed labels.',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-m', '--dist-min',
                            help='Minimum distance between two TSSs for them to be tested.',
                            default=0,
                            type=int,
                            required=False)

    parser_grp.add_argument('-k', '--keys',
                            help='Transcript-related keys to be exported (comma separated list).',
                            default=None,
                            metavar="DISTANCE",
                            type=str,
                            required=False)

    parser_grp.add_argument('-c', '--cluster-dist',
                            help="Maximum distance between gene TSSs to merged them into a TSS cluster.",
                            default=100,
                            type=int,
                            required=False)

    parser_grp.add_argument('-s', '--select-switch',
                            help="Write only potential promoter switch.",
                            action="store_true")

    parser_grp.add_argument('-g', '--select-gain-or-loss',
                            help="Write only potential promoter gain or loss.",
                            action="store_true")

    return parser


# -------------------------------------------------------------------------
# The main function
# -------------------------------------------------------------------------


def alt_prom(
        inputfile=None,
        outputfile=None,
        bed_list=None,
        keys=None,
        cluster_dist=100,
        labels=None,
        dist_min=1000,
        select_switch=False,
        select_gain_or_loss=False,
        tmp_dir=None,
        logger_file=None,
        verbosity=True):
    """
    Takes a GTF as input to search for genes with alternative promoters.
    """

    # -------------------------------------------------------------------------
    # Create a list of labels.
    # Take user input in account
    # -------------------------------------------------------------------------

    bed_list = [x.name for x in bed_list]

    if len(bed_list) != len(set(bed_list)):
        message("Found the same BED file several times.",
                type="ERROR")

    if len(bed_list) < 2:
        message("At least two bed files are needed.",
                type="ERROR")

    message('Checking labels.')

    if labels is not None:
        labels = labels.split(",")
        # Ensure the number of labels is the same as the number of bed files.
        if len(labels) != len(bed_list):
            message("The number of labels should be the same as the number of"
                    " bed files.", type="ERROR")
        # Ensure labels are non-redondant
        if len(labels) > len(set(labels)):
            message("Labels must be unique.", type="ERROR")
    else:
        labels = []
        for i in range(len(bed_list)):
            labels += [
                os.path.splitext(
                    os.path.basename(
                        bed_list[i]))[0]]

    # -------------------------------------------------------------------------
    # Get tss of transcript
    #
    # -------------------------------------------------------------------------

    gtf = GTF(inputfile)

    # A dict tx -> tss
    tss_dict = gtf.get_tss("transcript_id", as_dict=True, one_based=True)

    # Cluster of tss (force strandness)
    tss_clusters = dict()
    tss_bo = gtf.get_tss().sort().cluster(d=cluster_dist, s=True)
    for line in tss_bo:
        tss_clusters[line[3]] = line[6]

    # Some info for output
    gn_start_end = gtf.select_by_key("feature",
                                     "gene"
                                     ).extract_data("gene_id,chrom,gene_name,start,end,strand",
                                                    no_na=False,
                                                    as_dict_of_lists=True)

    tx_start_end = gtf.select_by_key("feature",
                                     "transcript"
                                     ).extract_data("transcript_id,start,end",
                                                    no_na=False,
                                                    as_dict_of_lists=True)

    if keys is not None:
        tx_info = gtf.select_by_key("feature",
                                    "transcript"
                                    ).extract_data("transcript_id," + keys,
                                                   no_na=False,
                                                   as_dict_of_lists=True)
    # gn to tx mapping
    gn_to_tx = gtf.get_gn_to_tx()

    # gn -> tx1|tx2 -> dist(tx2-tx1)
    gn_tss_dist = defaultdict(lambda: defaultdict(list))

    # tx to gene mapping
    tx_to_gn = gtf.get_tx_to_gn()

    # -------------------------------------------------------------------------
    # Loop over transcripts. Check whether their distance is greater to
    # distmin.
    # -------------------------------------------------------------------------

    for gn, v in gn_to_tx.items():
        for tx1 in v[:-1]:
            for tx2 in v[1:]:
                dist_tx1_tx2_tss = int(tss_dict[tx1]) - int(tss_dict[tx2])
                if abs(dist_tx1_tx2_tss) > dist_min:
                    if (tx2 + "|" + tx1) not in gn_tss_dist[gn]:
                        gn_tss_dist[gn][tx1 + "|" + tx2] = dist_tx1_tx2_tss

    tss_bo = gtf.get_tss(["transcript_id"]).sort()

    tx_to_peak = defaultdict(lambda: defaultdict(list))

    # -------------------------------------------------------------------------
    # Check intersection between peaks from the two samples
    #
    # -------------------------------------------------------------------------

    message("Check intersection between peaks from the two samples")

    # A dict to store bed_file_label -> peak from tissus 1 -> list of intersected peaks in tissues 2
    peak_intersections = defaultdict(lambda: defaultdict(list))

    bf1_bo = BedTool(bed_list[0])
    bf2_bo = BedTool(bed_list[1])

    tmp_intersect = make_tmp_file(prefix="peak_to_peak_intersections",
                                  suffix=".bed")
    tss_bo_x_bf = bf1_bo.intersect(bf2_bo, wb=True)
    tss_bo_x_bf.saveas(tmp_intersect.name)

    for line in tss_bo_x_bf:
        peak_intersections[labels[0]][line[3]] += [line[9]]
        peak_intersections[labels[1]][line[9]] += [line[3]]

    # -------------------------------------------------------------------------
    # Computing the number of promoter per gene
    # This is considered as the number of merged peaks from both conditions
    # intersecting a TSS from that gene
    # -------------------------------------------------------------------------

    merged_peaks = bf1_bo.cat(bf2_bo).merge()

    tss_bo_gid = gtf.get_tss(["gene_id"]).sort()

    tss_x_merged_peaks = tss_bo_gid.intersect(merged_peaks)

    gene_nb_promoters = defaultdict(lambda: 0)

    for line in tss_x_merged_peaks:
        gene_nb_promoters[line.name] += 1

    del gtf

    # -------------------------------------------------------------------------
    # Check whether TSSs intersect a H3K4me3 peak
    #
    # -------------------------------------------------------------------------

    for bf, lab in zip(bed_list, labels):

        message("Reading BED file : " + bf)

        bf_bo = BedTool(bf)

        message("Checking number of columns in BED.")

        if bf_bo.field_count() < 6:
            message("BED files should be in BED6 format",
                    type="ERROR")

        message("Intersecting BED files.")

        tss_bo_x_bf = tss_bo.intersect(bf_bo, wa=True, wb=True)
        tmp_bed = make_tmp_file(prefix="TSS_intersect_" + lab + "_",
                                suffix=".bed")
        tss_bo_x_bf.saveas(tmp_bed.name)

        for line in tss_bo_x_bf:
            tx_to_peak[line[3]][lab] += [line[9]]

    # -------------------------------------------------------------------------
    # Preparing header
    #
    # -------------------------------------------------------------------------

    if keys is not None:
        a = [str(i) + '_tx1' for i in keys.split(",")]
        b = [str(i) + '_tx2' for i in keys.split(",")]
        header_more = flatten_list(zip(a, b))
    else:
        header_more = []

    header_base = ["gene_id",
                   "chrom",
                   "gene_name",
                   "gene_start",
                   "gene_end",
                   "strand",
                   "transcript_id_1",
                   "transcript_id_2",
                   "tx1_start",
                   "tx1_end",
                   "tx2_start",
                   "tx2_end",
                   "tss_tx1",
                   "tss_tx2",
                   "tss_clusters_tx1",
                   "tss_clusters_tx2",
                   "tss_dist",
                   "tss_of_tx1_overlaps_H3K4me3_cond_1",
                   "tss_of_tx2_overlaps_H3K4me3_cond_1",
                   "tss_of_tx1_overlaps_H3K4me3_cond_2",
                   "tss_of_tx2_overlaps_H3K4me3_cond_2",
                   "tss_of_tx1_overlaps_both",
                   "tss_of_tx2_overlaps_both",
                   "potential_gain_or_loss",
                   "potential_switch",
                   "gene_nb_computed_promoters",
                   "gain_or_loss"]

    header = header_base + header_more

    outputfile.write("\t".join(header) + "\n")

    #
    # -------------------------------------------------------------------------
    # Write output
    #
    # -------------------------------------------------------------------------

    message("looping over transcript pair and printing")

    for gn in gn_tss_dist:
        for tup in gn_tss_dist[gn]:
            tx1, tx2 = tup.split("|")

            line = dict(zip(header_base,
                            [gn,
                             gn_start_end[gn][0],
                             gn_start_end[gn][1],
                             gn_start_end[gn][2],
                             gn_start_end[gn][3],
                             gn_start_end[gn][4],
                             tx1,
                             tx2,
                             tx_start_end[tx1][0],
                             tx_start_end[tx1][1],
                             tx_start_end[tx2][0],
                             tx_start_end[tx2][1],
                             str(tss_dict[tx1]),
                             str(tss_dict[tx2]),
                             tss_clusters[tx1],
                             tss_clusters[tx2],
                             str(gn_tss_dist[gn][tup]),
                             ] + ["0"] * 9 + [""]))

            # The number of promoter detected for that gene
            gn_id = tx_to_gn[tx1]
            line["gene_nb_computed_promoters"] = str(gene_nb_promoters[gn_id])

            ## Additional / facultative information about tx1 and tx2
            if keys is not None:
                key_split = keys.split(",")
                for k in range(len(key_split)):

                    if tx1 in tx_info:
                        line[key_split[k] + "_tx1"] = tx_info[tx1][k]
                    else:
                        line[key_split[k] + "_tx1"] = "."

                    if tx2 in tx_info:
                        line[key_split[k] + "_tx2"] = tx_info[tx2][k]
                    else:
                        line[key_split[k] + "_tx2"] = "."

            tx1_peaks = defaultdict(lambda: defaultdict(list))
            tx2_peaks = defaultdict(lambda: defaultdict(list))

            for i in range(len(labels)):

                if tx_to_peak[tx1][labels[i]] != []:
                    line["tss_of_tx1_overlaps_H3K4me3_cond_" + str(i + 1)] = ";".join(tx_to_peak[tx1][labels[i]])
                    for p in tx_to_peak[tx1][labels[i]]:
                        if peak_intersections[labels[i]][p]:
                            line["tss_of_tx1_overlaps_both"] = ";".join(peak_intersections[labels[i]][p])
                else:
                    line["tss_of_tx1_overlaps_H3K4me3_cond_" + str(i + 1)] = "0"

                if tx_to_peak[tx2][labels[i]] != []:
                    line["tss_of_tx2_overlaps_H3K4me3_cond_" + str(i + 1)] = ";".join(tx_to_peak[tx2][labels[i]])
                    for p in tx_to_peak[tx2][labels[i]]:
                        if peak_intersections[labels[i]][p]:
                            line["tss_of_tx2_overlaps_both"] = ";".join(peak_intersections[labels[i]][p])

                else:
                    line["tss_of_tx2_overlaps_H3K4me3_cond_" + str(i + 1)] = "0"

            cond_1 = (line["tss_of_tx1_overlaps_H3K4me3_cond_1"] == "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx1_overlaps_both"] == "0")

            if cond_1:
                line["gain_or_loss"] += "tx1_tss_gained_in_" + labels[1]

            cond_2 = (line["tss_of_tx1_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_2"] == "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx1_overlaps_both"] == "0")

            if cond_2:
                line["gain_or_loss"] += "tx1_tss_lost_in_" + labels[1]

            cond_3 = (line["tss_of_tx1_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_2"] == "0" and
                      line["tss_of_tx2_overlaps_both"] == "0")

            if cond_3:
                line["gain_or_loss"] += "tx2_tss_lost_in_" + labels[1]

            cond_4 = (line["tss_of_tx1_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_1"] == "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx2_overlaps_both"] == "0")

            if cond_4:
                line["gain_or_loss"] += "tx2_tss_gained_in_" + labels[1]

            if cond_1 or cond_2 or cond_3 or cond_4:
                line["potential_gain_or_loss"] = "1"

            cond_1 = (line["tss_of_tx1_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_2"] == "0" and
                      line["tss_of_tx1_overlaps_both"] == "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_1"] == "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx2_overlaps_both"] == "0")

            cond_2 = (line["tss_of_tx2_overlaps_H3K4me3_cond_1"] != "0" and
                      line["tss_of_tx2_overlaps_H3K4me3_cond_2"] == "0" and
                      line["tss_of_tx2_overlaps_both"] == "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_1"] == "0" and
                      line["tss_of_tx1_overlaps_H3K4me3_cond_2"] != "0" and
                      line["tss_of_tx1_overlaps_both"] == "0")

            if cond_1 or cond_2:
                line["potential_switch"] = "1"
                line["gain_or_loss"] = "both"

            full_list = [line[k] for k in header]
            if select_gain_or_loss:
                if line["potential_gain_or_loss"] == "1":
                    outputfile.write("\t".join(full_list) + "\n")
            elif select_switch:
                if line["potential_switch"] == "1":
                    outputfile.write("\t".join(full_list) + "\n")
            else:
                outputfile.write("\t".join(full_list) + "\n")

    close_properly(inputfile, outputfile)


# -------------------------------------------------------------------------
# Call  to main
# -------------------------------------------------------------------------


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    alt_prom(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    

    @test "alt_prom_1" {
     result=`gtftk get_example -d mini_real | gtftk alt_prom H3K4me3_cond_1.bed H3K4me3_cond_2.bed| awk 'BEGIN{FS=OFS="\t"}$3=="CRMP1"||NR==1'| grep ENST00000513911 | grep ENST00000324989 | cut -f 25`
      [ "$result" = "1" ]
    }

    @test "alt_prom_2" {
     result=`gtftk get_example -d mini_real | gtftk alt_prom H3K4me3_cond_1.bed H3K4me3_cond_3.bed| awk 'BEGIN{FS=OFS="\t"}$3=="CRMP1"||NR==1'| grep ENST00000513911 | grep ENST00000324989 | cut -f 24`
      [ "$result" = "1" ]
    }

    @test "alt_prom_3" {
     result=`gtftk get_example -d mini_real | gtftk alt_prom H3K4me3_cond_2.bed H3K4me3_cond_3.bed| awk 'BEGIN{FS=OFS="\t"}$3=="CRMP1"||NR==1'| grep ENST00000513911 | grep ENST00000324989 | cut -f 24`
      [ "$result" = "1" ]
    }
    """

    CmdObject(name='alt_prom',
              message='Search for genes with alternative promoters.',
              parser=make_parser(),
              fun=alt_prom,
              desc=__doc__,
              notes=__notes__,
              updated=__updated__,
              group="annotation",
              test=test)
