#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys
from _collections import defaultdict

from pybedtools import BedTool

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import bedFile, checkChromFile
from pygtftk.arg_formatter import int_ge_to_null
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly, make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
Get the genes/transcripts closest to a set of feature. The region to be search around the feature is controled by 
 -\-slop-value. The -\-name argument allows to export also some informations regarding these  genes/transcripts (e.g. biotype...).
"""
__notes__ = """
 -- Features (e.g. peaks, enhancer regions...) are provided in BED6 format.
 -- Output is in BED format with two additional columns: gene list and associated distances.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=FileWithExtension('r',
                                                   valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file (BED).",
                            default=sys.stdout,
                            metavar="BED",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Bb][Ee][Dd]$',
                                                                     '\.[Bb][Ee][Dd]6$')))

    parser_grp.add_argument('-r', '--region-file',
                            help="The input BED file containing regions of interest.",
                            metavar="BED",
                            required=True,
                            action=bedFile)

    parser_grp.add_argument('-t', '--ft-type',
                            help="The feature of interest.",
                            default='gene',
                            choices=['gene', 'transcript'],
                            required=False)

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name in the output BED file.",
                            default="gene_id,gene_name",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-p', '--slop-value',
                            help="Look until this number of nucleotides in 5' and 3' direction",
                            default=0,
                            type=int_ge_to_null,
                            required=False)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-c', '--chrom-info',
                            help='Tabulated file (chr as '
                                 'column 1, sizes as column 2.)',
                            default=None,
                            action=checkChromFile,
                            required=True)

    return parser


def closest_gn_to_feat(inputfile=None,
                       outputfile=None,
                       ft_type="gene",
                       names=None,
                       region_file=None,
                       slop_value=None,
                       chrom_info=None,
                       separator="|",
                       tmp_dir=None,
                       logger_file=None,
                       verbosity=0):
    """
    Get the gene/transcript closest to a set of feature.
    """

    # -------------------------------------------------------------------------
    # Load the GTF an get the regions of interest
    # -------------------------------------------------------------------------

    nms = names.split(",")

    gtf = GTF(inputfile, check_ensembl_format=True)

    ann_bo = gtf.select_by_key("feature",
                               ft_type).to_bed(name=nms)

    # -------------------------------------------------------------------------
    # Load the BED file, slop and intersect with gene/tx
    # -------------------------------------------------------------------------

    region_bo = BedTool(region_file).slop(s=False,
                                          b=slop_value,
                                          g=chrom_info.name)

    tmp_file = make_tmp_file("slopped_regions", ".bed")
    region_bo.saveas(tmp_file.name)

    if region_bo.field_count() < 6:
        message("BED should be in BED6 format.", type="ERROR")

    intersect_bo = region_bo.intersect(ann_bo, wa=True, wb=True)

    feat_to_gen = defaultdict(list)

    for feat in intersect_bo:
        cur_key = [feat[0],
                   int(feat[1]),
                   int(feat[2])] + feat[3:6]
        cur_val = [feat[6],
                   int(feat[7]),
                   int(feat[8])] + feat[9:12]

        feat_to_gen[tuple(cur_key)] += [tuple(cur_val)]

    # -------------------------------------------------------------------------
    # Compute distance to gene/transcript
    # -------------------------------------------------------------------------
    # Compare the positions of two fragments f1,f2
    #
    #    a --------------- b f1
    #
    #         c -------------- d  f2

    feat_to_dist = defaultdict(list)

    for feat, gen_list in feat_to_gen.items():

        for gen in gen_list:

            a = feat[1] + int(slop_value)
            b = feat[2] - int(slop_value)
            c = gen[1]
            d = gen[2]

            if a > d:
                dist_feat = a - d + 1
            elif a < d:
                if b > c or b == c:
                    dist_feat = 0
                else:
                    dist_feat = c - b + 1
            elif a == d:
                dist_feat = 0

            feat_to_dist[feat] += [dist_feat]

    for feat in feat_to_gen:

        output_str = []
        output_str += feat
        gene_list = []
        gene_dist = []

        for gn, dist in zip(feat_to_gen[feat], feat_to_dist[feat]):
            gene_list += [gn[3]]
            gene_dist += [dist]

        gene_list = [x for _, x in sorted(zip(gene_dist, gene_list))]
        gene_dist = sorted(gene_dist)
        gene_dist = [str(x) for x in gene_dist]
        col_right = ",".join(gene_list) + "\t" + ", ".join(gene_dist)

        feat = [str(x) for x in feat]
        outputfile.write("\t".join(feat) + "\t" + col_right + "\n")

    # -------------------------------------------------------------------------
    # Load the BED file, slop and intersect with gene/tx
    # -------------------------------------------------------------------------

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    closest_gn_to_feat(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #closest_gn_to_feat: 
    @test "closest_gn_to_feat_0" {
     result=`gtftk  get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
    
    #closest_gn_to_feat: 
    @test "closest_gn_to_feat_1" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -n gene_id| md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "7883485b27bc9bda4a7c47846c7a0d25" ]
    }
        
    
    """
    CmdObject(name="closest_gn_to_feat",
              message="Get the list of genes/transcripts closest to a set of feature.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              notes=__notes__,
              updated=__updated__,
              desc=__doc__,
              group="annotation",
              test=test)
