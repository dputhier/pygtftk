#!/usr/bin/env python
"""
Get the list of attributes from a GTF file.
"""
import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-20"


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


def get_attr_list(
        inputfile=None,
        outputfile=None,
        separator="\n"):
    """
    Get the list of attributes from a GTF file.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)
    attr_list = gtf.get_attr_list()
    n = 0
    for i in attr_list:
        if n != len(attr_list) - 1:
            outputfile.write(i + separator)
        else:
            outputfile.write(i)
        n += 1
    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_attr_list(**args)


if __name__ == '__main__':
    main()

else:

    test = """


    #get_attr_list: load dataset
    @test "get_attr_list_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
       
       
    #get_attr_list
    @test "get_attr_list_1" {
     result=`gtftk get_attr_list -i simple.gtf| grep exon| cut -f2 `
      [ "$result" = "exon_id" ]
    }

    #get_attr_list
    @test "get_attr_list_2" {
     result=`gtftk get_attr_list -i simple.gtf | wc -l`
      [ "$result" -eq 3 ]
    }

    #get_attr_list
    @test "get_attr_list_3" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m | gtftk get_attr_list| wc -l`
      [ "$result" -eq 5 ]
    }
        
    """

    CMD = CmdObject(name="get_attr_list",
                    message="Get the list of attributes from a GTF file.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    group="information",
                    desc=__doc__,
                    test=test)
