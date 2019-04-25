#!/usr/bin/env python
"""Convert a GTF to various format (still limited)."""

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-20"


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN.",
                            default=sys.stdin,
                            metavar="GTF",
                            required=False,
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="BED/BED3/BED6",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='bed'))

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name.",
                            default="gene_id,transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-m', '--more-names',
                            help="Add this information to the 'name' column of the BED file.",
                            default="",
                            type=str)

    parser_grp.add_argument('-f', '--format',
                            help='Currently one of bed3, bed6',
                            type=str,
                            choices=('bed', 'bed3', 'bed6'),
                            default='bed6',
                            required=False)

    return parser


def convert(inputfile=None,
            outputfile=None,
            format="bed",
            names="gene_id,transcript_id",
            separator="|",
            more_names=''):
    """
 Convert a GTF to various format.
    """

    if format == "bed3":
        gtf = GTF(inputfile, check_ensembl_format=False)

        for i in gtf.extract_data("seqid,start,end", as_list_of_list=True, hide_undef=False, no_na=False):
            i[1] = str(int(i[1]) - 1)
            outputfile.write("\t".join(i) + "\n")

    elif format in ["bed", "bed6"]:
        gtf = GTF(inputfile,
                  check_ensembl_format=False).write_bed(outputfile=outputfile,
                                                        name=names,
                                                        sep=separator,
                                                        more_name=more_names)

    close_properly(outputfile, inputfile)


def main():
    """The main program."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    convert(**args)


if __name__ == '__main__':
    main()

else:

    test = '''

    # convert: load dataset
    @test "convert_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
            
    # Convert:
    @test "convert_1" {
     result=`gtftk convert -i simple.gtf | awk 'BEGIN{FS="\\t"}{print NF}'| sort | uniq`
      [ "$result" -eq 6 ]
    }
    
    
    # Convert: basic.
    @test "convert_2" {
     result=`gtftk convert -f bed3 -i simple.gtf | awk 'BEGIN{FS="\\t"}{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    # Convert: check name.
    @test "convert_3" {
     result=`gtftk convert -i simple.gtf -n gene_id,transcript_id,start | cut -f4| awk 'BEGIN{FS="|"}{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    # Convert: check zero based (bed6)
    @test "convert_4" {
     result=`gtftk convert -i simple.gtf -n gene_id,transcript_id,start | cut -f2| head -n 1`
      [ "$result" -eq 124 ]
    }
    # Convert: check zero based (bed3)
    @test "convert_4" {
     result=`gtftk convert -i simple.gtf -f bed3 | cut -f2| head -n 1`
      [ "$result" -eq 124 ]
    }
    '''

    CmdObject(name="convert",
              message="Convert a GTF to various format including bed.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              desc=__doc__,
              group="conversion",
              test=test)
