#!/usr/bin/env python
"""
 Get the midpoint coordinates for the requested feature.
 Output is bed format.
"""

import argparse
import os
import sys

from pybedtools import BedTool

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF/BED",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('bed', 'gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file (BED).",
                            default=sys.stdout,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='bed'))

    parser_grp.add_argument('-t', '--ft-type',
                            help="The target feature (as found in the 3rd "
                                 "column of the GTF).",
                            default='transcript',
                            type=str,
                            required=False)

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name.",
                            default="transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    return parser


def midpoints(
        inputfile=None,
        outputfile=None,
        ft_type="transcript",
        names="transcript_id",
        separator="|"):
    """
 Get the midpoint coordinates for the requested feature.
    """

    message("Loading input file...")
    if inputfile.name == '<stdin>':
        is_gtf = True
    else:
        region_bo = BedTool(inputfile.name)
        if len(region_bo) == 0:
            message("Unable to find requested regions",
                    type="ERROR")

        if region_bo.file_type == 'gff':
            is_gtf = True
        else:
            is_gtf = False

    if is_gtf:

        gtf = GTF(inputfile.name, check_ensembl_format=False)

        bed_obj = gtf.select_by_key("feature",
                                    ft_type).get_midpoints(name=names.split(","),
                                                           sep=separator)
        for line in bed_obj:
            write_properly(chomp(str(line)), outputfile)

    else:
        for line in region_bo:

            diff = line.end - line.start

            if diff % 2 != 0:
                # e.g 10-13 (zero based) -> 11-13 one based
                # mipoint is 12 (one-based) -> 11-12 (zero based)
                # e.g 949-1100 (zero based) -> 950-1100 one based
                # mipoint is 1025 (one-based) -> 1024-1025 (zero based)
                # floored division (python 2)...
                line.end = line.start + int(diff // 2) + 1
                line.start = line.end - 1
            else:
                # e.g 10-14 (zero based) -> 11-14 one based
                # mipoint is 12-13 (one-based) -> 11-13 (zero based)
                # e.g 9-5100 (zero based) -> 10-5100 one based
                # mipoint is 2555-2555 (one-based) -> 2554-2555 (zero based)
                # floored division (python 2)...
                # No real center. Take both

                line.start = line.start + int(diff // 2) - 1
                line.end = line.start + 2

            outputfile.write(str(line))

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    midpoints(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #midpoints: load dataset
    @test "midpoints_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
        

    #midpoints: test line number
    @test "midpoints_1" {
     result=`gtftk midpoints -i simple.gtf -t exon | wc -l`
      [ "$result" -eq 25 ]
    }
    
    #midpoints: test line transcript
    @test "midpoints_2" {
     result=`gtftk midpoints -i simple.gtf -t transcript | wc -l`
      [ "$result" -eq 15 ]
    }
    
    
    #midpoints: test line gene
    @test "midpoints_3" {
     result=`gtftk midpoints -i simple.gtf -t gene | wc -l`
      [ "$result" -eq 10 ]
    }
    
    #midpoints: test coordinates
    @test "midpoints_4" {
     result=`gtftk midpoints -i simple.gtf -t transcript| cut -f2| perl -npe 's/\\n/,/'`
      [ "$result" = "7,7,27,30,39,54,69,69,110,110,130,130,180,183,215," ]
    }
    
    #midpoints: test coordinates
    @test "midpoints_5" {
     result=`gtftk midpoints -i simple.gtf -t transcript| cut -f3| perl -npe 's/\\n/,/'`
      [ "$result" = "9,9,29,32,40,56,71,71,112,112,132,132,181,185,216," ]
    }
    
    #midpoints: test coordinates
    @test "midpoints_6" {
     result=`gtftk midpoints -i simple.gtf -t exon| cut -f2| perl -npe 's/\\n/,/'`
      [ "$result" = "7,7,22,28,28,33,33,33,43,51,58,65,65,70,70,74,74,110,110,130,130,180,183,211,220," ]
    }
    
    
    #midpoints: test coordinates
    @test "midpoints_7" {
     result=`gtftk midpoints -i simple.gtf -t exon| cut -f3| perl -npe 's/\\n/,/'`
      [ "$result" = "9,9,24,29,29,34,34,34,45,52,59,67,67,71,71,75,75,112,112,132,132,181,185,212,221," ]
    }

    #midpoints: test coordinates
    @test "midpoints_8" {
     result=`gtftk get_example| gtftk midpoints| gtftk bed_to_gtf|  gtftk midpoints | cut -f3| perl -npe 's/\\n/,/'`
      [ "$result" = "9,9,29,32,40,56,71,71,112,112,132,132,181,185,216," ]
    }
    
    #midpoints: test coordinates
    @test "midpoints_9" {
     result=`gtftk midpoints -i simple_peaks.bed | cut -f2| perl -npe 's/\\n/,/g'`
      [ "$result" = "12,18,23,38,68,43," ]
    }
    
    #midpoints: test coordinates
    @test "midpoints_10" {
     result=`gtftk midpoints -i simple_peaks.bed | cut -f3| perl -npe 's/\\n/,/g'`
      [ "$result" = "13,19,24,39,70,44," ]
    }
        
    
    """

    CmdObject(name="midpoints",
              message=" Get the midpoint coordinates for the requested feature.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="coordinates",
              updated=__updated__,
              desc=__doc__,
              test=test)
