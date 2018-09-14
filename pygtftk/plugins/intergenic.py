#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys

from builtins import str

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import checkChromFile
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"
__doc__ = """
 Extract intergenic regions. This command requires a chromInfo file to compute
 the bed file boundaries. The command will print the coordinates of genomic
 regions without any transcript features.
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
                                                   valid_extensions='\.[Bb][Ee][Dd]$'))

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Tabulated two-columns file. Chromosomes"
                                 " as column 1 and their sizes as"
                                 " column 2",
                            default=None,
                            metavar="CHROMINFO",
                            action=checkChromFile,
                            required=True)

    return parser


def intergenic(
        inputfile=None,
        outputfile=None,
        chrom_info=None,
        tmp_dir=None,
        logger_file=None,
        verbosity=0):
    """
 Extract intergenic regions.
    """

    message("Searching for intergenic regions.")

    gtf = GTF(inputfile)

    intergenic_regions = gtf.get_intergenic(chrom_info)

    nb_intergenic_region = 1

    for i in intergenic_regions:
        i.name = "region_" + str(nb_intergenic_region)
        write_properly(chomp(str(i)), outputfile)
        nb_intergenic_region += 1

    close_properly(outputfile, inputfile)


def main():
    """The main program."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    intergenic(**args)


if __name__ == '__main__':
    main()

else:
    test = """
        
    #intergenic: check few coord
    @test "intergenic_1" {
     result=`gtftk intergenic -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo | head -3| tail -1 | cut -f2`
      [ "$result" -eq 47 ]
    }
    
    #intergenic: check few coord
    @test "intergenic_2" {
     result=`gtftk intergenic -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo | wc -l`
      [ "$result" -eq 10 ]
    }
    
    
    #intergenic: check chr
    @test "intergenic_3" {
     result=`gtftk intergenic -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo | cut -f1| sort|uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "chr1,chr2," ]
    }
    """
    CmdObject(name="intergenic",
              message="Extract intergenic regions.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="coordinates",
              desc=__doc__,
              updated=__updated__,
              test=test)
