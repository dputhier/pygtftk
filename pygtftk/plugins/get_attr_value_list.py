#!/usr/bin/env python
"""
Get the list of values observed for an attributes.
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

    parser_grp.add_argument('-k', '--key-name',
                            help="Key name.",
                            default=None,
                            type=str,
                            required=True)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating key names.",
                            default="\n",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-c', '--count',
                            help="Add the counts for each key (separator will be set to '\t' by default).",
                            action="store_true")

    return parser


def get_attr_value_list(
        inputfile=None,
        outputfile=None,
        key_name=None,
        separator="\n",
        count=False):
    """
    Get the list of values observed for an attributes.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)

    if not count:
        for i in gtf.get_attr_value_list(key_name):
            outputfile.write(i + separator)
        close_properly(outputfile, inputfile)
    else:
        if separator == "\n":
            separator = "\t"
        for i in gtf.get_attr_value_list(key_name, count=True):
            outputfile.write(i[0] + separator + i[1] + "\n")
        close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_attr_value_list(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    
    # get_attr_value_list: load dataset
    @test "get_attr_value_list_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
    
    #get_attr_value_list
    @test "get_attr_value_list_1" {
     result=`gtftk get_attr_value_list -i simple.gtf -k feature | perl -npe 's/\\n/,/g'`
      [ "$result" = "CDS,exon,gene,transcript," ]
    }

    #get_attr_value_list
    @test "get_attr_value_list_2" {
     result=`gtftk get_attr_value_list -i simple.gtf -k feature  | wc -l`
      [ "$result" -eq 4 ]
    }

    #get_attr_value_list
    @test "get_attr_value_list_3" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m | gtftk get_attr_value_list -k S1 | wc -l `
      [ "$result" -eq 3 ]
    }
        
    #get_attr_value_list
    @test "get_attr_value_list_4" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m | gtftk get_attr_value_list -k S1 | perl -npe 's/\\n/,/'`
      [ "$result" = "0.2322,0.5555,0.999," ]
    }
        
    """

    CMD = CmdObject(name="get_attr_value_list",
                    message="Get the list of values observed for an attributes.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    group="information",
                    desc=__doc__,
                    test=test)
