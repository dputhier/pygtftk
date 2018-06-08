#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys

from gtftk.arg_formatter import FileWithExtension
from gtftk.arg_formatter import checkChromFile
from gtftk.bedtool_extension import BedTool
from gtftk.bwig.bw_coverage import bw_profile_mp
from gtftk.cmd_object import CmdObject
from gtftk.gtf_interface import GTF
from gtftk.utils import close_properly, message
from gtftk.utils import make_tmp_file

__updated__ = "2018-02-05"
__doc__ = """
Takes a GTF or BED as input. Compute bigwig coverage around the center of the features. This can be interesting to create 
a signal matrix  that capture signal shape.
 """
__notes__ = """
-- Regions were signal can be computed (if GTF file as input): promoter, tts,
introns, intron_by_tx, intergenic regions or any feature available in the GTF file (transcript, exon, gene...).
-- If bed is used as input, each region should have its own name (column 4).
-- The output is a BED file with additional columns 
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
        'bw_list',
        help='A list of Bigwig files (last argument).',
        type=argparse.FileType('r'),
        nargs='+')

    parser_grp.add_argument('-i', '--inputfile',
                            help="The input GTF/BED file. Only GTF file if <stdin> is used.",
                            default=sys.stdin,
                            metavar="GTF/BED",
                            required=False,
                            type=FileWithExtension('r',
                                                   valid_extensions=('\.[Gg][Tt][Ff](\.[Gg][Zz])?$',
                                                                     '\.[Bb][Ee][Dd]$',
                                                                     '\.[Bb][Ee][Dd]3$',
                                                                     '\.[Bb][Ee][Dd]6$')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Tt][Xx][Tt]$',
                                                                     '\.[Cc][Ss][Vv]$',
                                                                     '\.[Tt][Aa][Bb]$',
                                                                     '\.[Tt][Ss][Vv]$',
                                                                     '\.[Cc][Oo][Vv]$',
                                                                     '\.[Bb][Ee][Dd]$')))
    parser_grp.add_argument(
        '-c', '--chrom-info',
        help="Tabulated two-columns file. Chromosomes"
             " as column 1 and sizes as"
             " column 2 ",
        default=None,
        metavar="CHROMINFO",
        action=checkChromFile,
        required=True)

    parser_grp.add_argument(
        '-u', '--upstream',
        help="Number of windows upstream of the midpoint.",
        default=10,
        metavar="UPSTREAM",
        type=int,
        required=False)

    parser_grp.add_argument(
        '-d', '--downstream',
        help="Number of windows downstream of the midpoint.",
        default=10,
        metavar="DOWNSTREAM",
        type=int,
        required=False)

    parser_grp.add_argument(
        '-s', '--window-size',
        type=int,
        default=25,
        help='The size of the window.',
        required=False)

    parser_grp.add_argument(
        '-k', '--nb-proc',
        type=int,
        default=1,
        help='Use this many threads to compute coverage.',
        required=False)

    parser_grp.add_argument(
        '-f', '--ft-type',
        type=str,
        default="promoter",
        help="Region in which coverage is to be computed (promoter, intron, intergenic, tts or any feature defined in the column 3 of the GTF).",
        required=False)

    parser_grp.add_argument(
        '-l', '--labels',
        help='Bigwig labels.',
        default=None,
        type=str,
        required=False)

    parser_grp.add_argument(
        '-m', '--name-column',
        type=str,
        default="transcript_id",
        help="Use this ids (comma separated) to compute the name (4th column in "
             "bed output).",
        required=False)

    parser_grp.add_argument(
        '-p',
        '--pseudo-count',
        type=int,
        default=1,
        help='A pseudo-count to add in case count is equal to 0.')

    parser_grp.add_argument(
        '-nz', '--na-to-zero',
        help='Use zero not NA when region is undefined in bigwig.',
        action='store_true',
        required=False)

    return parser


# -------------------------------------------------------------------------
# The main function
# -------------------------------------------------------------------------

def cov_around_mid(
        inputfile=None,
        outputfile=None,
        bw_list=None,
        labels=None,
        pseudo_count=1,
        ft_type="promoter",
        window_size=10,
        downstream=100,
        na_to_zero=False,
        name_column=None,
        upstream=100,
        chrom_info=None,
        nb_proc=1,
        tmp_dir=None,
        logger_file=None,
        verbosity=True):
    """
    Compute transcript coverage with one or several bigWig.
    """

    # -------------------------------------------------------------------------
    # Create a list of labels.
    # Take user input in account
    # -------------------------------------------------------------------------

    bw_list = [x.name for x in bw_list]

    if len(bw_list) != len(set(bw_list)):
        message("Found the same bigwigs several times.",
                type="ERROR")

    message('Checking labels.')

    if labels is not None:
        labels = labels.split(",")
        # Ensure the number of labels is the same as the number of bw files.
        if len(labels) != len(bw_list):
            message("The number of labels should be the same as the number of"
                    " bigwig files.", type="ERROR")
        # Ensure labels are non-redondant
        if len(labels) > len(set(labels)):
            message("Labels must be unique.", type="ERROR")
    else:
        labels = []
        for i in range(len(bw_list)):
            labels += [
                os.path.splitext(
                    os.path.basename(
                        bw_list[i]))[0]]

    # -------------------------------------------------------------------------
    # Check input file is in bed or GTF format
    #
    # -------------------------------------------------------------------------

    message("Loading input file...")
    if inputfile.name == '<stdin>':
        gtf = GTF(inputfile.name)
        is_gtf = True
    else:
        region_bo = BedTool(inputfile.name)
        if len(region_bo) == 0:
            message("Unable to find requested regions",
                    type="ERROR")

        if region_bo.file_type == 'gff':
            gtf = GTF(inputfile.name)
            is_gtf = True
        else:
            is_gtf = False

    # -------------------------------------------------------------------------
    # Get regions of interest
    #
    # -------------------------------------------------------------------------

    if is_gtf:

        name_column = name_column.split(",")
        if "transcript_id" not in name_column:
            message("'transcript_id' is mandatory in --name-column.")

        message("Getting regions of interest...")

        if ft_type == "intergenic":

            region_bo = gtf.get_intergenic(chrom_info, 0, 0).get_midpoints().slop(s=True,
                                                                                  l=upstream * window_size,
                                                                                  r=downstream * window_size,
                                                                                  g=chrom_info.name).sort()

        elif ft_type == "intron":

            region_bo = gtf.get_introns().get_midpoints().slop(s=True,
                                                               l=upstream * window_size,
                                                               r=downstream * window_size,
                                                               g=chrom_info.name).sort()

        elif ft_type == "intron_by_tx":
            name_column = name_column.split(",")
            region_bo = gtf.get_introns(by_transcript=True,
                                        name_column=name_column).get_midpoints().slop(s=True,
                                                                                      l=upstream * window_size,
                                                                                      r=downstream,
                                                                                      g=chrom_info.name).sort()

        elif ft_type in ["promoter", "tss"]:
            region_bo = gtf.get_tss(name=name_column
                                    ).slop(s=True,
                                           l=upstream * window_size,
                                           r=downstream * window_size,
                                           g=chrom_info.name).sort()

        elif ft_type == "tts":
            region_bo = gtf.get_tts(name=name_column).slop(s=True,
                                                           l=upstream * window_size,
                                                           r=downstream * window_size,
                                                           g=chrom_info.name).sort()

        else:
            region_bo = gtf.select_by_key("feature", ft_type, 0).get_midpoints(name=name_column
                                                                               ).slop(s=True,
                                                                                      l=upstream * window_size,
                                                                                      r=downstream * window_size,
                                                                                      g=chrom_info.name).sort()

        if len(region_bo) == 0:
            message("Unable to find requested regions",
                    type="ERROR")

    else:
        ft_type = "region"
        region_bo = region_bo.slop(s=True,
                                   l=upstream * window_size,
                                   r=downstream * window_size,
                                   g=chrom_info.name).sort()

    region_bed = make_tmp_file(prefix="region", suffix=".bed")

    region_bo.saveas(region_bed.name)

    # -------------------------------------------------------------------------
    # Compute coverage
    #
    # -------------------------------------------------------------------------

    tmp_file = make_tmp_file(prefix="cov_around_mid", suffix=".bed")

    result_bed = bw_profile_mp(in_bed_file=region_bed,
                               nb_proc=nb_proc,
                               big_wig=bw_list,
                               bin_nb=upstream + downstream,
                               pseudo_count=pseudo_count,
                               stranded=False,
                               type=ft_type,
                               labels=labels,
                               outputfile=tmp_file.name,
                               bed_format=True,
                               add_score=True,
                               zero_to_na=not na_to_zero,
                               verbose=verbosity)

    nb_line = 0
    for i in open(result_bed.name):
        outputfile.write(i)
        nb_line += 1

    if nb_line == 0:
        message("No line available in output...",
                type="ERROR")

    close_properly(inputfile, outputfile)


# -------------------------------------------------------------------------
# Call  to main
# -------------------------------------------------------------------------


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    cov_around_mid(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    
    @test "cov_around_mid_1" {
     result=`gtftk get_example | gtftk cov_around_mid -c gtftk/data/simple/simple.chromInfo -i gtftk/data/simple/simple.gtf gtftk/data/simple/simple.bw -u 1 -d 1 -s 1 -f tts -p 0 | grep G0009T001| perl -npe 's/\\t/|/g'`
     [ "$result" = "chr1|1|4|simple|G0009T001|.|-|1.0|1.0" ]
    }
    
    @test "cov_around_mid_2" {
      result=`gtftk get_example | gtftk cov_around_mid -c gtftk/data/simple/simple.chromInfo -i gtftk/data/simple/simple.gtf gtftk/data/simple/simple.bw -u 1 -d 1 -s 1 -f tts -p 0 | grep G0006T001| perl -npe 's/\\t/|/g'`
      [ "$result" = "chr1|20|23|simple|G0006T001|.|-|3.0|3.0" ]
    }
    
    @test "cov_around_mid_3" {
      result=`gtftk get_example | gtftk cov_around_mid -c gtftk/data/simple/simple.chromInfo -i gtftk/data/simple/simple.gtf gtftk/data/simple/simple.bw -u 7 -d 0  -s 1 -f tss -p 0  -K toto | grep G0010T001| perl -npe 's/\\t/|/g'`
        [ "$result" = "chr1|168|176|simple|G0010T001|.|+|3.0|3.0|2.0|2.0|1.0|1.0|0.0" ]
    } 
        
    @test "cov_around_mid_4" {
      result=`gtftk get_example | gtftk cov_around_mid -c gtftk/data/simple/simple.chromInfo -i gtftk/data/simple/simple.gtf gtftk/data/simple/simple.bw -u 5 -d 0  -s 1 -f CDS -p 0  -K toto | perl -npe 's/\\t/|/g' | grep G0001T001`
        [ "$result" = "chr1|125|131|simple|G0001T001|.|+|4.0|3.0|3.0|2.0|2.0" ]
    }     

    @test "cov_around_mid_5" {
      result=`gtftk get_example | gtftk cov_around_mid -c gtftk/data/simple/simple.chromInfo -i gtftk/data/simple/simple.gtf gtftk/data/simple/simple.bw -u 5 -d 0  -s 1 -f exon| grep G0004T001| wc -l`
        [ "$result" -eq 3 ]
    }   
    """

    CmdObject(name='cov_around_mid',
              message='Compute bigwig coverage around midpoints of features.',
              parser=make_parser(),
              fun=cov_around_mid,
              desc=__doc__,
              notes=__notes__,
              updated=__updated__,
              group="coverage",
              test=test)
