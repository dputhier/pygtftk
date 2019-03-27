#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
 Select transcripts based on the number of exons.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('-i', '--inputfile',
                        help="Path to the GTF file. Default to STDIN",
                        default=sys.stdin,
                        metavar="GTF",
                        type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser.add_argument('-o', '--outputfile',
                        help="Output file.",
                        default=sys.stdout,
                        metavar="GTF",
                        type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser.add_argument('-m', '--min-exon-number',
                        help="Minimum number of exons.",
                        default=0,
                        type=int)

    parser.add_argument('-M', '--max-exon-number',
                        help="Maximum number of exons.",
                        default=None,
                        type=int)
    return parser


def select_by_nb_exon(inputfile=None,
                      outputfile=None,
                      min_exon_number=None,
                      max_exon_number=None):
    """
    Select transcripts based on the number of exons.
    """

    msg = "Selecting transcript by exon number (range: [{m},{M}])"
    msg = msg.format(m=str(min_exon_number),
                     M=str(max_exon_number))
    message(msg)

    gtf = GTF(inputfile,
              check_ensembl_format=False
              ).select_by_number_of_exons(min_exon_number,
                                          max_exon_number)

    gtf.write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_nb_exon(**args)


if __name__ == '__main__':
    main()


else:

    test = """
    #select_by_nb_exons
    @test "select_by_nb_exons_1" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 1 -M 1 | gtftk select_by_key -k feature -v transcript| wc -l`
      [ "$result" -eq 8 ]
    }
    #select_by_nb_exons
    @test "select_by_nb_exons_2" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 1 -M 1 | gtftk select_by_key -k feature -v exon| wc -l`
      [ "$result" -eq 8 ]
    }

    #select_by_nb_exons
    @test "select_by_nb_exons_3" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 3 -M 3 | gtftk select_by_key -k feature -v exon| wc -l`
      [ "$result" -eq 9 ]
    }

    #select_by_nb_exons
    @test "select_by_nb_exons_4" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 3 -M 3 | gtftk select_by_key -k feature -v transcript| wc -l`
      [ "$result" -eq 3 ]
    }

    #select_by_nb_exons
    @test "select_by_nb_exons_5" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 1 -M 100 | gtftk select_by_key -k feature -v exon| wc -l`
      [ "$result" -eq 25 ]
    }

    #select_by_nb_exons
    @test "select_by_nb_exons_6" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 1 -M 100 | wc -l`
      [ "$result" -eq 60 ]
    }
    
    #select_by_nb_exons
    @test "select_by_nb_exons_7" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 0 -M 0 | wc -l`
      [ "$result" -eq 0 ]
    }

    #select_by_nb_exons
    @test "select_by_nb_exons_8" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 2 -M 2 | gtftk select_by_key -k feature -v exon| wc -l`
      [ "$result" -eq 8 ]
    }

    #select_by_nb_exons
    @test "select_by_nb_exons_9" {
     result=`gtftk get_example | gtftk select_by_nb_exon -m 2 -M 2 | gtftk select_by_key -k feature -v transcript| wc -l`
      [ "$result" -eq 4 ]
    }
    """

    CmdObject(name="select_by_nb_exon",
              message="Select transcripts based on the number of exons.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              updated=__updated__,
              desc=__doc__,
              test=test)
