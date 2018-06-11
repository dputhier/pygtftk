#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-20"
__doc__ = """
 Select lines from a GTF file based on a boolean test on numeric values.
"""

__notes__ = """
 -- The test should be of the form : "{key} operator value".
 -- e.g. "{transcript_version} > 10"
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('-i', '--inputfile',
                        help="Path to the GTF file. Default to STDIN",
                        default=sys.stdin,
                        metavar="GTF",
                        type=FileWithExtension('r',
                                               valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser.add_argument('-o', '--outputfile',
                        help="Output file.",
                        default=sys.stdout,
                        metavar="GTF",
                        type=FileWithExtension('w',
                                               valid_extensions='\.[Gg][Tt][Ff]$'))

    parser.add_argument('-t', '--test',
                        help='The test to be applied.',
                        default=None,
                        type=str,
                        required=True)

    parser.add_argument('-n', '--na-omit',
                        help='If one of the evaluated values is enclosed in this list (csv), line is skipped.',
                        default=None,
                        type=str,
                        required=False)

    return parser


def select_by_numeric_value(inputfile=None,
                            outputfile=None,
                            test=None,
                            na_omit=None,
                            tmp_dir=None,
                            logger_file=None,
                            verbosity=0):
    """Select lines from a GTF file based on a boolean test on numeric values.
    """

    gtf = GTF(inputfile, check_ensembl_format=False
              ).select_by_numeric_value(test,
                                        na_omit=na_omit,
                                        ).write(outputfile)
    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_numeric_value(**args)


if __name__ == '__main__':
    main()


else:

    test = '''
    #select_by_numeric_value_1
    @test "select_by_numeric_value_1" {
      result=`gtftk join_attr -i pygtftk/data/simple/simple.gtf  -j pygtftk/data/simple/simple.join_mat -k gene_id -m|  gtftk select_by_numeric_value -t 'start < 10 and end > 10 and S1 == 0.5555 and S2 == 0.007e2' -n "."| wc -l`
      [ "$result" -eq 5 ]
    }

    '''

    CmdObject(name="select_by_numeric_value",
              message="Select lines from a GTF file based on a boolean test on numeric values.",
              parser=make_parser(),
              fun=select_by_numeric_value,
              group="selection",
              desc=__doc__,
              updated=__updated__,
              test=test)
