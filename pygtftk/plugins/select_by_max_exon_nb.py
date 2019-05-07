#!/usr/bin/env python
"""
For each gene select the transcript with the highest number of exons. If ties, select the first encountered.
"""
import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import message

__updated__ = "2018-02-11"


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

    return parser


def select_by_max_exon_nb(inputfile=None,
                          outputfile=None):
    """
    Select transcripts based on the number of exons.
    """

    msg = "Selecting transcript with the highest number of exon for each gene."
    message(msg)

    gtf = GTF(inputfile,
              check_ensembl_format=False
              ).select_by_max_exon_nb()

    gtf.write(outputfile, gc_off=True)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_max_exon_nb(**args)


if __name__ == '__main__':
    main()


else:

    test = """
    #select_by_max_exon_nb
    @test "select_by_max_exon_nb_1" {
     result=`gtftk get_example  -d simple_04 | gtftk select_by_max_exon_nb | grep G0005T001 | wc -l`
      [ "$result" -eq 0 ]
    }

    #select_by_max_exon_nb
    @test "select_by_max_exon_nb_2" {
     result=`gtftk get_example  -d simple_04 | gtftk select_by_max_exon_nb | grep G0004T002 | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #select_by_max_exon_nb
    @test "select_by_max_exon_nb_3" {
     result=`gtftk get_example  -d simple_04 | gtftk select_by_max_exon_nb | grep G0006T002 | wc -l`
      [ "$result" -eq 0 ]
    }
    
    """

    CmdObject(name="select_by_max_exon_nb",
              message="For each gene select the transcript with the highest number of exons.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              updated=__updated__,
              desc=__doc__,
              test=test)
