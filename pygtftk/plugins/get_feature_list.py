#!/usr/bin/env python
"""
Get the list of features enclosed in the GTF.
"""

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-02-11"


def make_parser():
    """The parser."""
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
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='txt'))

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating key names.",
                            default="\n",
                            metavar="SEP",
                            type=str)

    return parser


def get_feature_list(
        inputfile=None,
        outputfile=None,
        separator=""):
    """
    Get the list of features enclosed in the GTF.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)
    for i in gtf.get_feature_list(nr=True):
        outputfile.write(str(i) + separator)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_feature_list(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    @test "get_feature_list_1" {
     result=`gtftk get_example -d mini_real | gtftk get_feature_list| wc -l`
      [ "$result" -eq 9 ]
    }

    @test "get_feature_list_2" {
     result=`gtftk get_example -d mini_real | gtftk get_feature_list -s ","`
      [ "$result" = "gene,transcript,exon,five_prime_utr,CDS,start_codon,stop_codon,three_prime_utr,Selenocysteine," ]
    }
    
    """

    CMD = CmdObject(name="get_feature_list",
                    message="Get the list of features enclosed in the GTF.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="information",
                    test=test)
