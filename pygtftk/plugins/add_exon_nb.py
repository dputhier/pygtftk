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
 Add exon number transcript-wise (based on 5' to 3' orientation).
"""


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
                            help="Output file.",
                            default=sys.stdout,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-k', '--exon-numbering-key',
                            help="The name of the key containing the exon numbering.",
                            default="exon_nbr",
                            type=str)

    return parser


def add_exon_nb(inputfile=None,
                outputfile=None,
                exon_numbering_key=None):
    """Add the exon number to each exon (based on 5' to 3' orientation)."""

    message("Calling nb_exons.", type="DEBUG")

    GTF(inputfile.name,
        check_ensembl_format=False
        ).add_exon_number(exon_numbering_key
                          ).write(outputfile, gc_off=True)

    close_properly(inputfile, outputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    add_exon_nb(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #add_exon_nb_1: on sorted gtf. gene from - strand.
    @test "add_exon_nb_1" {
     result=`gtftk get_example| bedtools sort | gtftk add_exon_nb| gtftk select_by_key -k transcript_id -v G0006T001  | gtftk select_by_key -k feature -v exon | gtftk tabulate -H  -k exon_nbr| perl -npe 's/\\n/,/'`
      [ "$result" = "3,2,1," ]
    }
    
    #add_exon_nb_2:  on sorted gtf. gene from + strand.
    @test "add_exon_nb_2" {
     result=`gtftk get_example| bedtools sort | gtftk add_exon_nb| gtftk select_by_key -k transcript_id -v G0004T002  | gtftk select_by_key -k feature -v exon | gtftk tabulate -H  -k exon_nbr| perl -npe 's/\\n/,/'`
      [ "$result" = "1,2,3," ]
    }

    #add_exon_nb_3: on randomized gtf. Gene from - strand.
    @test "add_exon_nb_3" {
     result=`gtftk get_example| perl -MList::Util -e 'print List::Util::shuffle <>'  | gtftk add_exon_nb| gtftk select_by_key -k transcript_id -v G0006T001  | gtftk select_by_key -k feature -v exon | bedtools sort | gtftk tabulate -H  -k exon_nbr| perl -npe 's/\\n/,/'`
      [ "$result" = "3,2,1," ]
    }
    
    #add_exon_nb_4: on randomized gtf. Gene from + strand.
    @test "add_exon_nb_4" {
     result=`gtftk get_example| perl -MList::Util -e 'print List::Util::shuffle <>'  | gtftk add_exon_nb| gtftk select_by_key -k transcript_id -v G0004T002  | gtftk select_by_key -k feature -v exon | bedtools sort | gtftk tabulate -H  -k exon_nbr| perl -npe 's/\\n/,/'`
      [ "$result" = "1,2,3," ]
    }

    #add_exon_nb_5: all
    @test "add_exon_nb_5" {
     result=`gtftk get_example| gtftk add_exon_nb  | bedtools sort| gtftk select_by_key -k feature -v exon| gtftk tabulate -H -k start,exon_nbr| perl -npe 's/\\t/,/g; s/\\n/,/g'`
      [ "$result" = "3,1,3,1,22,3,28,2,28,2,33,1,33,2,33,1,42,1,50,2,57,1,65,1,65,1,71,2,71,2,74,3,74,3,107,1,107,1,125,1,125,1,176,1,180,1,210,2,220,1," ]
    }
    
    #add_exon_nb_5: all
    @test "add_exon_nb_6" {
     result=`gtftk get_example| gtftk add_exon_nb  -o add_exon_nb_6.gtf`
      [ "$result" = "" ]
    }    
    """

    CmdObject(name="add_exon_nb",
              message="Add exon number transcript-wise.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              desc=__doc__,
              group="information",
              test=test)
