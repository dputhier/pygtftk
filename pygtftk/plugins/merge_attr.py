#!/usr/bin/env python
"""
 Merge a set of attributes into a destination attribute. Can be
 useful, for instance, to merge gene_name and gene_id values into
 a new key to prepare the GTF for RNA-seq quantification.
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
-- The destination key can be one of the source key, leading to an update of that key.
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

    parser_grp.add_argument('-k', '--src-key',
                            help='comma-separated list of keys to join.',
                            default=None,
                            metavar="KEY",
                            type=str,
                            required=True)

    parser_grp.add_argument('-d', '--dest-key',
                            help='The target key name.',
                            default=None,
                            metavar="KEY",
                            type=str,
                            required=True)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator for the concatenated values.",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-f', '--target-feature',
                            help='The name of the target feature.',
                            default="*",
                            type=str,
                            required=False)
    return parser


def merge_attr(
        inputfile=None,
        outputfile=None,
        src_key="gene_id,transcript_id",
        separator="|",
        target_feature="*",
        dest_key="gene_tx_ids"):
    """
    Merge a set of attributes into a destination attribute.
    """

    GTF(inputfile,
        check_ensembl_format=False
        ).merge_attr(target_feature,
                     src_key,
                     dest_key,
                     separator).write(outputfile,
                                      gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    merge_attr(**args)


if __name__ == '__main__':
    main()


else:
    test = """
   
    #merge_attr: check basic args.
    @test "merge_attr_1" {
     result=`gtftk get_example | gtftk merge_attr -k transcript_id,gene_id -d gene_tx_concat|gtftk tabulate -H -k gene_tx_concat|grep G0009T001| sort | uniq`
      [ "$result" = "G0009T001|G0009" ]
    }

    @test "merge_attr_2" {
     result=`gtftk get_example | gtftk merge_attr -k transcript_id,gene_id -d gene_tx_concat -f transcript| gtftk select_by_key -k feature -v exon| gtftk tabulate -k gene_tx_concat -H | wc -l`
      [ "$result" -eq 0 ]
    }

    @test "merge_attr_3" {
     result=`gtftk get_example | gtftk merge_attr -k transcript_id,gene_id -d gene_tx_concat -f transcript| gtftk select_by_key -k feature -v transcript| gtftk tabulate -k gene_tx_concat -H | sort | uniq | wc -l`
      [ "$result" -eq 15 ]
    }
    
    @test "merge_attr_4" {
     result=`gtftk get_example -d mini_real | gtftk merge_attr -k end,start,transcript_id,gene_id -d gene_tx_concat -f transcript| gtftk select_by_key -k feature -v transcript| wc -l`
      [ "$result" -eq 8531 ]
    }
    
    @test "merge_attr_5" {
     result=`gtftk get_example -d mini_real | gtftk merge_attr -k end,start,transcript_id,gene_id -d gene_tx_concat -f exon| gtftk select_by_key -k feature -v exon| wc -l`
      [ "$result" -eq 64251 ]
    }

    @test "merge_attr_6" {
     result=`gtftk get_example |  gtftk merge_attr -k transcript_id,gene_id -d transcript_id -s "|" -f transcript | gtftk select_by_key -t| gtftk tabulate -k transcript_id -H | perl -ne 'print (length $_,"\\n")'| sort | uniq`
      [ "$result" -eq 16 ]
    }    
      
    @test "merge_attr_7" {
     result=`gtftk get_example |  gtftk merge_attr -k transcript_id,gene_id -d transcript_id -s "|" -f transcript | gtftk select_by_key -e| gtftk tabulate -k transcript_id -H | perl -ne 'print (length $_,"\\n")'| sort | uniq`
      [ "$result" -eq 10 ]
    }    

    @test "merge_attr_8" {
     result=`gtftk get_example |  gtftk merge_attr -k transcript_id,gene_id -d transcript_id -s "|" -f transcript | gtftk select_by_key -g| gtftk tabulate -k transcript_id -H | perl -ne 'print (length $_,"\\n")'| sort | uniq| wc -l`
      [ "$result" -eq 0 ]
    }    
     
    """
    msg = "Merge a set of attributes into a destination attribute."
    CmdObject(name="merge_attr",
              message=msg,
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="editing",
              updated=__updated__,
              notes=__notes__,
              desc=__doc__,
              test=test)
