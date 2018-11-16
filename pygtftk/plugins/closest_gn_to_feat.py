#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys
from _collections import defaultdict

from pybedtools import BedTool

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import bed6
from pygtftk.arg_formatter import checkChromFile
from pygtftk.arg_formatter import int_ge_to_null
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly, make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
Get the genes/transcripts/tss/tts in the neighborhood of peaks  ordered by distance. The region to be search around the peak is controled by 
 -\-slop-value. The -\-name argument allows to export also some information regarding these  genes/transcripts (e.g. biotype...). You may ask 
 for a gene-centric ouput (print gene/transcript/tss coordinates instead of peak coordinates).
"""

__notes__ = """
 -- Features (e.g. peaks, enhancer regions...) are provided in BED6 format.
 -- If -\-ft_type is set to tss or tts the closest tss or tts of each transcript will be considered.
 -- Output is in BED format with two additional columns: gene list (csv) and associated distance (csv). 
 -- You may ask for -\-gene-centric to print gene/transcript/tss coordinates instead of peak coordinates).
 -- Use -\-uncollapse to wite one peak to gene/distance per line (no csv).
 
>>>>>>> closest_gn_to_feat
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
                            action=bed6)

    parser_grp.add_argument('-t', '--ft-type',
                            help="The feature of interest.",
                            default='gene',
                            choices=['gene', 'transcript', 'tss', 'tts'],
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

    parser_grp.add_argument('-u', '--uncollapse',
                            help='Write one peak to gene/distance per line (no comma separated list).',
                            action="store_true")

    parser_grp.add_argument('-g', '--gene-centric',
                            help='Print --feature-type coordinate and the peak list and associated distance.',
                            action="store_true")

    parser_grp.add_argument('-c', '--chrom-info',
                            help='Tabulated file (chr as column 1, sizes as column 2.)',
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
                       uncollapse=False,
                       tmp_dir=None,
                       logger_file=None,
                       gene_centric=False,
                       verbosity=0):
    """
    Get the gene/transcript closest to a set of feature.
    """

    collapsed = not uncollapse

    # -------------------------------------------------------------------------
    # Load the GTF an get the regions of interest
    # -------------------------------------------------------------------------

    nms = names.split(",")

    gtf = GTF(inputfile, check_ensembl_format=True)

    if ft_type in ['gene', 'transcript']:
        ann_bo = gtf.select_by_key("feature",
                                   ft_type).to_bed(name=nms)
    elif ft_type == 'tss':
        ann_bo = gtf.get_tss(name=nms)
    elif ft_type == 'tts':
        ann_bo = gtf.get_tts(name=nms)
    else:
        message("Unsupported feature (see --ft_type).",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Load the BED file
    # -------------------------------------------------------------------------

    region_bo = BedTool(region_file)

    # -------------------------------------------------------------------------
    # Store the start and ends of peaks. As region_file is a bed6 action
    # names are unambiguous.
    # -------------------------------------------------------------------------

    peak_starts = defaultdict()
    peak_ends = defaultdict()

    for i in region_bo:
        peak_starts[i.name] = i.start
        peak_ends[i.name] = i.end

    # -------------------------------------------------------------------------
    # Slop the regions
    # -------------------------------------------------------------------------

    region_bo = region_bo.slop(s=False,
                               b=slop_value,
                               g=chrom_info.name)

    tmp_file = make_tmp_file("slopped_regions", ".bed")
    region_bo.saveas(tmp_file.name)

    # -------------------------------------------------------------------------
    # Intersect
    # -------------------------------------------------------------------------

    intersect_bo = region_bo.intersect(ann_bo, wa=True, wb=True)

    closest = defaultdict(list)

    for feat in intersect_bo:
        cur_key = [feat[0],
                   int(feat[1]),
                   int(feat[2])] + feat[3:6]
        cur_val = [feat[6],
                   int(feat[7]),
                   int(feat[8])] + feat[9:12]

        closest[tuple(cur_key)] += [tuple(cur_val)]

    # -------------------------------------------------------------------------
    # Compute distance to gene/transcript
    # -------------------------------------------------------------------------
    # Compare the positions of two fragments f1,f2
    #
    #                      a --------------- b f1
    #                     22                25
    #    c -------------- d  f2
    #    21               22

    closest_dist = defaultdict(list)

    for peak, gen_list in closest.items():

        for gen in gen_list:

            a = peak_starts[peak[3]]
            b = peak_ends[peak[3]]
            c = gen[1]
            d = gen[2]

            if a > d:
                dist_peak = a - d + 1
            elif a < d:
                if b >= c:
                    dist_peak = 0
                else:
                    dist_peak = c - b + 1
            elif a == d:
                # zero based and (!) half open
                dist_peak = 1

            closest_dist[peak] += [dist_peak]

    # -------------------------------------------------------------------------
    # Loop over peak or genes
    # -------------------------------------------------------------------------

    if not gene_centric:

        # -------------------------------------------------------------------------
        # Keep only gene names in the list.
        # Order gene lists based on distance.
        # -------------------------------------------------------------------------

        for peak in closest:
            gene_list = [x for _, x in sorted(zip(closest_dist[peak], closest[peak]))]
            gene_dist = sorted(closest_dist[peak])

            gene_dist = [str(x) for x in gene_dist]

            closest[peak] = [x[3] for x in gene_list]
            closest_dist[peak] = gene_dist

        for peak in closest:

            if collapsed:

                # -------------------------------------------------------------------------
                # print csv of genes and dists for each peak
                # -------------------------------------------------------------------------

                col_right = ",".join(closest[peak]) + "\t" + ", ".join(closest_dist[peak])
                peak = [str(x) for x in peak]
                outputfile.write("\t".join(peak) + "\t" + col_right + "\n")

            else:

                # -------------------------------------------------------------------------
                # for each peak print a gene and a dist
                # -------------------------------------------------------------------------

                for gn, dist in zip(closest[peak], closest_dist[peak]):
                    col_right = "\t".join([gn, str(dist)])
                    peak = [str(x) for x in peak]
                    outputfile.write("\t".join(peak) + "\t" + col_right + "\n")


    else:

        # -------------------------------------------------------------------------
        # Put gene as key and peak list as values
        # -------------------------------------------------------------------------

        closest_cp = defaultdict(list)
        closest_dist_cp = defaultdict(list)

        for peak in closest:

            for gn, dist in zip(closest[peak], closest_dist[peak]):
                gn = tuple([str(x) for x in gn])
                closest_cp[gn] += [peak]
                closest_dist_cp[gn] += [dist]

        # -------------------------------------------------------------------------
        # Keep only peak names in the list.
        # Order peak lists based on distance.
        # -------------------------------------------------------------------------

        for gn in closest_cp:
            peak_list = [x for _, x in sorted(zip(closest_dist_cp[gn], closest_cp[gn]))]
            peak_dist = sorted(closest_dist_cp[gn])

            peak_dist = [str(x) for x in peak_dist]

            # get peak name
            closest_cp[gn] = [x[3] for x in peak_list]
            closest_dist_cp[gn] = peak_dist

        if collapsed:

            # -------------------------------------------------------------------------
            # for each gene print csv of peaks and dists
            # -------------------------------------------------------------------------

            for gn in closest_cp:
                col_right = ",".join(closest_cp[gn]) + "\t" + ",".join(closest_dist_cp[gn])
                outputfile.write("\t".join(gn) + "\t" + col_right + "\n")
        else:

            # -------------------------------------------------------------------------
            # for each gene print a peak and a dist
            # -------------------------------------------------------------------------

            for gene in closest_cp:

                for peak, dist in zip(closest_cp[gene], closest_dist_cp[gene]):
                    col_right = "\t".join([peak, str(dist)])
                    outputfile.write("\t".join(gene) + "\t" + col_right + "\n")

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
    
    #closest_gn_to_feat: check gene
    @test "closest_gn_to_feat_1" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -n gene_id| md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "1ea2fbb5a5f10c607836139757c4f6de" ]
    }

    #closest_gn_to_feat: check tx
    @test "closest_gn_to_feat_2" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -t transcript -n transcript_id | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "a90cf09d38d28399ab480df0c258d53f" ]
    }
    
    
    #closest_gn_to_feat: check tss
    @test "closest_gn_to_feat_3" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -t tss -n transcript_id  | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "aa89ecf2a0c86f8ff3790589a0436dcc" ]
    }

    #closest_gn_to_feat: check tts
    @test "closest_gn_to_feat_4" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -t tts -n transcript_id -K toto | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "0d3e7b17314469ecd8b2ffff65acaa29" ]
    }

    #closest_gn_to_feat: check tts, gene centric
    @test "closest_gn_to_feat_5" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -t tss -n transcript_id -K toto -g | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "b76208a2c0a384ba3368af5251090039" ]
    }


    #closest_gn_to_feat: check tts, gene centric, uncollapsed
    @test "closest_gn_to_feat_6" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -t tss -n transcript_id -K toto -gu | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "2e4683b4f7ea8a79c8c7c632bec16e9e" ]
    }

    #closest_gn_to_feat: check gene, peak centric, uncollapsed
    @test "closest_gn_to_feat_7" {
     result=`gtftk closest_gn_to_feat -r simple_peaks.bed6 -i simple.gtf -c simple.chromInfo -p 10 -n gene_id -u | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = "9721e0bf42fea4ca4f3d0738154d69f0" ]
    }

    """
    CmdObject(name="closest_gn_to_feat",
              message="Get the list of genes/transcripts/tss/tts closest to a set of peaks.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              notes=__notes__,
              updated=__updated__,
              desc=__doc__,
              group="annotation",
              test=test)
