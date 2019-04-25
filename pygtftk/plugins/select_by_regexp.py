#!/usr/bin/env python
"""
 Select lines from a GTF file based on a regexp.
"""

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-20"

__notes__ = """
 -- The default is to try to select feature from conventional human chromosome (chr1..chr22, chrX and chrY) with
 -\-key set to chrom and -\-regexp set to "^chr[0-9XY]+$".
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

    parser.add_argument('-k', '--key',
                        help='The key name',
                        default="chrom",
                        metavar="KEY",
                        type=str)

    parser.add_argument('-r', '--regexp',
                        help='The regular expression.',
                        default="^chr[0-9XY]+$",
                        type=str)

    parser.add_argument('-n', '--invert-match',
                        help='Not/invert match. Selected lines whose requested key '
                             'do not match the regexp.',
                        action="store_true")

    return parser


def select_by_regexp(inputfile=None,
                     outputfile=None,
                     key=None,
                     regexp=None,
                     invert_match=False):
    """Select lines from a GTF file based on attributes and
    associated values.
    """

    GTF(inputfile,
        check_ensembl_format=False).select_by_regexp(key,
                                                     regexp,
                                                     invert_match
                                                     ).write(outputfile,
                                                             gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_regexp(**args)


if __name__ == '__main__':
    main()


else:

    test = '''
    # select_by_regexp
    @test "select_by_regexp_1" {
      result=`gtftk get_example -d mini_real | gtftk select_by_regexp -r "^BCL" -k gene_name | gtftk tabulate -k gene_name -Hun| wc -l`
      [ "$result" -eq 2 ]
    }

    # select_by_regexp
    @test "select_by_regexp_2" {
      result=`gtftk get_example -d mini_real | gtftk select_by_regexp -r "^B.*B$" -k gene_name | gtftk tabulate -k gene_name -Hun| wc -l`
      [ "$result" -eq 1 ]
    }

    # select_by_regexp
    @test "select_by_regexp_3" {
      result=`gtftk get_example -d mini_real | gtftk select_by_regexp -r "^B.*B$" -k gene_name | wc -l`
      [ "$result" -eq 128 ]
    }
    
    # select_by_regexp
    @test "select_by_regexp_4" {
      result=`gtftk get_example -d mini_real | gtftk select_by_regexp -r "^B.*B$" -k gene_name -n | wc -l`
      [ "$result" -eq 137542 ]
    }
          
    '''

    CmdObject(name="select_by_regexp",
              message="Select lines from a GTF file based on a regexp.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              desc=__doc__,
              notes=__notes__,
              updated=__updated__,
              test=test)
