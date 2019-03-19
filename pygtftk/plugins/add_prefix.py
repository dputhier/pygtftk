#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-20"
__doc__ = """
  Add a prefix to target values. By default add 'chr' to seqid/chromosome key.
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

    parser_grp.add_argument('-k', '--key',
                            help="The name of the attribute for which a "
                                 "prefix/suffix is to be added to the corresponding"
                                 " values (e.g, gene_id, transcript_id...).",
                            default="chrom",
                            metavar="KEY",
                            type=str)

    parser_grp.add_argument('-t', '--text',
                            help='The character string to add as a prefix to the values.',
                            default="chr",
                            metavar="TEXT",
                            type=str)

    parser_grp.add_argument('-s', '--suffix',
                            help='The character string to add as a prefix to the values.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-f', '--target-feature',
                            help='The name of the target feature.',
                            default="*",
                            type=str,
                            required=False)

    return parser


def add_prefix(inputfile=None,
               outputfile=None,
               key="transcript_id",
               text=None,
               target_feature="*",
               suffix=False):
    """
    Add a prefix to target values.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)

    gtf.add_prefix(target_feature,
                   key,
                   text,
                   suffix).write(outputfile,
                                 gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    add_prefix(**args)


if __name__ == '__main__':
    main()


else:
    test = """
    #add_prefix: test target key
    @test "add_prefix_1" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t toto | awk '$3=="gene"'| grep toto | wc -l`
      [ "$result" -eq 0 ]
    }
    
    
    ##add_prefix: test line number
    @test "add_prefix_2" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t toto | wc -l`
      [ "$result" -eq 70 ]
    }
    
    
    #add_prefix: test target key
    @test "add_prefix_3" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t toto | awk '$3=="gene"'| grep toto | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #add_prefix: test number of line
    @test "add_prefix_4" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t toto | grep toto | wc -l`
      [ "$result" -eq 60 ]
    }
    
    
    #add_prefix: test line number (suffix True)
    @test "add_prefix_5" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t toto -s | wc -l`
      [ "$result" -eq 70 ]
    }
    
    
    #add_prefix: test target key (suffix True)
    @test "add_prefix_6" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t toto -s | awk '$3=="gene"'| grep toto | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #add_prefix: test number of line (suffix True)
    @test "add_prefix_7" {
     result=`gtftk get_example | gtftk add_prefix  -k transcript_id -t toto -s | grep toto | wc -l`
      [ "$result" -eq 60 ]
    }
    
    #add_prefix: test stdin
    @test "add_prefix_8" {
     result=`gtftk get_example | gtftk  add_prefix -k transcript -t toto | wc -l`
      [ "$result" -eq 70 ]
    }

    #add_prefix: test stdin
    @test "add_prefix_9" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t "|blabla" -s -f transcript| gtftk select_by_key -k feature -v transcript,exon| gtftk tabulate -k feature,transcript_id | grep exon | sort | cut -f2 | uniq| perl -ne 's/\\n/,/;print'`
      [ "$result" = "G0001T001,G0001T002,G0002T001,G0003T001,G0004T001,G0004T002,G0005T001,G0006T001,G0006T002,G0007T001,G0007T002,G0008T001,G0009T001,G0009T002,G0010T001," ]
    }

    #add_prefix: test stdin
    @test "add_prefix_10" {
     result=`gtftk get_example | gtftk add_prefix -k transcript_id -t "|blabla" -s -f transcript| gtftk select_by_key -k feature -v transcript,exon| gtftk tabulate -k feature,transcript_id | grep blabla | sort | cut -f1 | uniq`
      [ "$result" = "transcript" ]
    }

    #add_prefix: test large dataset
    @test "add_prefix_11" {
     result=`gtftk get_example -d mini_real | gtftk add_prefix -k transcript_id -t toto | wc -l`
      [ "$result" -eq 137670 ]
    }

    #add_prefix: test large dataset
    @test "add_prefix_12" {
     result=`gtftk get_example -d mini_real | gtftk add_prefix -k transcript_id -t toto| gtftk select_by_key -k feature -v transcript | wc -l`
      [ "$result" -eq 8531 ]
    }

    #add_prefix: test large dataset
    @test "add_prefix_13" {
     result=`gtftk get_example -d mini_real | gtftk add_prefix -k transcript_id -t toto| gtftk select_by_key -k feature -v transcript | gtftk tabulate -H -k transcript_id,gene_id| wc -l`
      [ "$result" -eq 8531 ]
    }
    
    #add_prefix: test large dataset
    @test "add_prefix_14" {
     result=`gtftk get_example -d mini_real | gtftk add_prefix -k transcript_id -t toto| awk '$3=="transcript"'| gtftk tabulate -H -k transcript_id,gene_id| wc -l`
      [ "$result" -eq 8531 ]
    }
    


    """

    CmdObject(name="add_prefix",
              message="Add a prefix or suffix to target values. ",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              desc=__doc__,
              group="editing",
              test=test)
