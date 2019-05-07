#!/usr/bin/env python
"""
 Extract intergenic regions. This command requires a chromInfo file to compute
 the bed file boundaries. The command will print the coordinates of genomic
 regions without any transcript features.
"""

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.arg_formatter import CheckChromFile
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"

__notes__ = '''
 -- -\-chrom-info may also accept 'mm8', 'mm9', 'mm10', 'hg19', 'hg38', 'rn3' or 'rn4'. In this case the 
 corresponding size of conventional chromosomes are used. ChrM is not used.  
'''


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file (BED).",
                            default=sys.stdout,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='bed'))

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Tabulated two-columns file. Chromosomes"
                                 " as column 1 and their sizes as"
                                 " column 2",
                            default=None,
                            metavar="CHROMINFO",
                            action=CheckChromFile,
                            required=True)

    return parser


def intergenic(
        inputfile=None,
        outputfile=None,
        chrom_info=None):
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

    # intergenic: load dataset
    @test "intergenic_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
            
    #intergenic: check few coord
    @test "intergenic_1" {
     result=`gtftk intergenic -i simple.gtf -c simple.chromInfo | head -3| tail -1 | cut -f2`
      [ "$result" -eq 47 ]
    }
    
    #intergenic: check few coord
    @test "intergenic_2" {
     result=`gtftk intergenic -i simple.gtf -c simple.chromInfo | wc -l`
      [ "$result" -eq 10 ]
    }
    
    
    #intergenic: check chr
    @test "intergenic_3" {
     result=`gtftk intergenic -i simple.gtf -c simple.chromInfo | cut -f1| sort|uniq| perl -npe 's/\\n/,/'`
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
              notes=__notes__,
              test=test)
