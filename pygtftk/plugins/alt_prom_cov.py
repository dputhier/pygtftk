#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys
from _collections import defaultdict

from pybedtools import BedTool
from scipy.stats import fisher_exact

from gtftk.arg_formatter import FileWithExtension
from gtftk.arg_formatter import bedFileWithUnambiguousNames
from gtftk.arg_formatter import int_greater_than_null
from gtftk.bwig.bw_coverage import bw_cov_mp
from gtftk.cmd_object import CmdObject
from gtftk.gtf_interface import GTF
from gtftk.utils import close_properly, message
from gtftk.utils import flatten_list
from gtftk.utils import make_tmp_file

__updated__ = "2018-02-15"
__doc__ = """
Takes a GTF as input, two H3K4me3 bed files and corresponding bigwigs. Search for genes with alternative promoter usage between two conditions.
 """
__notes__ = """
-- This program is very similar to alt_prom. However it will also use bigwig coverage to call alternative promoters.
"""


# -------------------------------------------------------------------------
# The main parser
# -------------------------------------------------------------------------


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)
    parser_grp = parser.add_argument_group('Arguments')

    # required but not required...

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

    parser_grp.add_argument('-b', '--bed-list',
                            help='A list of bed files for H3K4me3 epigenetic mark.',
                            type=bedFileWithUnambiguousNames(),
                            nargs='+',
                            required=True)

    parser_grp.add_argument('-w', '--bigwig_list',
                            help='A list of bigwigs files for H3K4me3 epigenetic mark.',
                            type=argparse.FileType('r'),
                            nargs='+')

    parser_grp.add_argument('-l', '--labels',
                            help='Bed/Bigwig labels.',
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

    parser_grp.add_argument('-p', '--nb-proc',
                            help='Use this many threads to compute coverage.',
                            default=1,
                            type=int_greater_than_null,
                            required=False)

    return parser


# -------------------------------------------------------------------------
# The main function
# -------------------------------------------------------------------------


def alt_prom_cov(
        inputfile=None,
        outputfile=None,
        bed_list=None,
        keys=None,
        cluster_dist=100,
        labels=None,
        dist_min=1000,
        bigwig_list=[],
        nb_proc=1,
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

    if len(bigwig_list) != 2:
        message("The program requires two bigwig files.", type="ERROR")

    if len(bed_list) != 2:
        message("The program requires two bed files.", type="ERROR")

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

        # Ensure the number of labels is the same as the number of bw files.
        if len(labels) != len(bigwig_list):
            message("The number of labels should be the same as the number of"
                    " bigwig files.", type="ERROR")
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

    # -------------------------------------------------------------------------
    # Merge peaks from the two samples
    #
    # -------------------------------------------------------------------------

    message("Merging peaks from the two samples.")

    bf1_bo = BedTool(bed_list[0])
    bf2_bo = BedTool(bed_list[1])

    merged_peaks = bf1_bo.cat(bf2_bo).merge()

    tmp_file = make_tmp_file(prefix="merged_peaks",
                             suffix=".bed")

    n = 1
    for line in merged_peaks:
        out_list = [line.chrom,
                    line.start,
                    line.end,
                    "merged_peak_" + str(n),
                    "0",
                    "."]
        tmp_file.write("\t".join([str(x) for x in out_list]) + "\n")
        n += 1

    tmp_file.close()

    merged_peaks = BedTool(tmp_file.name)

    # -------------------------------------------------------------------------
    # Computing the number of promoter per gene
    # This is considered as the number of merged peaks from both conditions
    # intersecting a TSS from that gene
    # -------------------------------------------------------------------------

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

    message("Checking TSS intersection with merged BED file.")

    tx_to_peak = defaultdict()

    merged_peaks_bo = BedTool(merged_peaks)

    tss_bo_x_bf = tss_bo.intersect(merged_peaks_bo, wa=True, wb=True)
    tmp_bed = make_tmp_file(prefix="TSS_intersect",
                            suffix=".bed")
    tss_bo_x_bf.saveas(tmp_bed.name)

    for line in tss_bo_x_bf:
        tx_to_peak[line[3]] = line[9]

    # -------------------------------------------------------------------------
    # Compute coverage of requested region
    # Each worker will send a file
    # -------------------------------------------------------------------------

    outputfile_list = {}
    message("Starting coverage analysis.")

    bed_cov_file = make_tmp_file(prefix="merged_peaks", suffix=".cov")
    bed_cov = bw_cov_mp(bw_list=[x.name for x in bigwig_list],
                        region_file=open(merged_peaks.fn),
                        labels=labels,
                        bin_nb=1,
                        nb_proc=nb_proc,
                        n_highest=None,
                        zero_to_na=False,
                        pseudo_count=1,
                        verbose=verbosity)

    bed_cov_dict = defaultdict(dict)
    for line in bed_cov:
        field = line.split("\t")
        lab, peak = field[3].split("|")
        bed_cov_dict[lab][peak] = field[4]

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

    header_base = [ "condition_1",
                    "condition_2",
                    "gene_id",
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
                   "coverage_tss_of_tx1_H3K4me3_cond_1",
                   "coverage_tss_of_tx1_H3K4me3_cond_2",
                   "coverage_tss_of_tx2_H3K4me3_cond_1",
                   "coverage_tss_of_tx2_H3K4me3_cond_2",
                   "multiple_promoters",
                   "fisher_exact_test",
                   "tx1_peak",
                   "tx2_peak"]

    header = header_base + header_more

    outputfile.write("\t".join(header) + "\n")

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
                             ] + ["NA"] * 8))

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

            cont_table = [["NA", "NA"], ["NA", "NA"]]

            line["condition_1"] = labels[0]
            line["condition_2"] = labels[1]

            for i in range(len(labels)):

                tx1_peak = None
                tx2_peak = None

                if tx1 in tx_to_peak:
                    tx1_peak = tx_to_peak[tx1]
                    val = bed_cov_dict[labels[i]][tx1_peak]
                    line["coverage_tss_of_tx1_H3K4me3_cond_" + str(i + 1)] = int(float(val))
                    tx1_has_peak = True
                    cont_table[i][0] = int(float(val))
                    line["tx1_peak"] = tx1_peak
                else:
                    line["coverage_tss_of_tx1_H3K4me3_cond_" + str(i + 1)] = "NA"
                    line["tx1_peak"] = "NA"

                if tx2 in tx_to_peak:
                    tx2_peak = tx_to_peak[tx2]
                    val = bed_cov_dict[labels[i]][tx2_peak]
                    line["coverage_tss_of_tx2_H3K4me3_cond_" + str(i + 1)] = int(float(val))
                    tx2_has_peak = True
                    cont_table[i][1] = int(float(val))
                    line["tx2_peak"] = tx2_peak
                else:
                    line["coverage_tss_of_tx2_H3K4me3_cond_" + str(i + 1)] = "NA"
                    line["tx2_peak"] = "NA"

                if tx1_peak is not None and tx2_peak is not None:
                    if tx1_peak != tx2_peak:
                        line["multiple_promoters"] = 1

            flat_list = flatten_list(cont_table)

            if "NA" not in flat_list:
                line["fisher_exact_test"] = fisher_exact(cont_table)[1]
            else:
                line["fisher_exact_test"] = "NA"

            full_list = [str(line[k]) for k in header]
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
    alt_prom_cov(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    

    @test "alt_prom_cov_1" {
     result=`gtftk get_example -d mini_real -f "*" | gtftk alt_prom_cov H3K4me3_cond_1.bed H3K4me3_cond_2.bed| awk 'BEGIN{FS=OFS="\t"}$3=="CRMP1"||NR==1'| grep ENST00000513911 | grep ENST00000324989 | cut -f 25`
      [ "$result" = "1" ]
    }

    @test "alt_prom_cov_2" {
     result=`gtftk get_example -d mini_real | gtftk alt_prom_cov H3K4me3_cond_1.bed H3K4me3_cond_3.bed| awk 'BEGIN{FS=OFS="\t"}$3=="CRMP1"||NR==1'| grep ENST00000513911 | grep ENST00000324989 | cut -f 24`
      [ "$result" = "1" ]
    }

    @test "alt_prom_cov_3" {
     result=`gtftk get_example -d mini_real | gtftk alt_prom_cov H3K4me3_cond_2.bed H3K4me3_cond_3.bed| awk 'BEGIN{FS=OFS="\t"}$3=="CRMP1"||NR==1'| grep ENST00000513911 | grep ENST00000324989 | cut -f 24`
      [ "$result" = "1" ]
    }
    """

    CmdObject(name='alt_prom_cov',
              message='Search for genes with alternative promoters.',
              parser=make_parser(),
              fun=alt_prom_cov,
              desc=__doc__,
              notes=__notes__,
              updated=__updated__,
              group="annotation",
              test=test)
