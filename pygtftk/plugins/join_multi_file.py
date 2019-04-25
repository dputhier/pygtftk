#!/usr/bin/env python
"""
Join attributes from mutiple files.
"""

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-02-05"


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

    parser_grp.add_argument('-k', '--key-to-join',
                            help='The name of the key used to join (e.g transcript_id).',
                            default=None,
                            metavar="KEY",
                            type=str,
                            required=True)

    parser_grp.add_argument('-t', '--target-feature',
                            help='The name(s) of the target feature(s). Comma-separated.',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-m', '--matrix-files',
                            help="A set of matrix files with row names as "
                                 "target keys column names as novel key and each cell as value.",
                            type=argparse.FileType('r'),
                            nargs='+')

    return parser


def join_multi_file(
        inputfile=None,
        outputfile=None,
        target_feature=None,
        key_to_join=None,
        matrix_files=()):
    """
    Join attributes from a set of tabulated files.
    """

    # -----------------------------------------------------------
    #  load the GTF
    # -----------------------------------------------------------

    gtf = GTF(inputfile, check_ensembl_format=False)

    # -----------------------------------------------------------
    #  Check target feature
    # -----------------------------------------------------------

    feat_list = gtf.get_feature_list(nr=True)

    if target_feature is not None:
        target_feature_list = target_feature.split(",")

        for i in target_feature_list:
            if i not in feat_list + ["*"]:
                message("Feature " + i + " not found.",
                        type="ERROR")
    else:
        target_feature = ",".join(feat_list)

    # -----------------------------------------------------------
    #  Do it
    # -----------------------------------------------------------

    for join_file in matrix_files:
        gtf = gtf.add_attr_from_matrix_file(feat=target_feature,
                                            key=key_to_join,
                                            inputfile=join_file.name)
    gtf.write(outputfile,
              gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    join_multi_file(**args)


if __name__ == '__main__':
    main()
else:

    test = """

    #join_multi_file: load dataset
    @test "join_multi_file_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
        
                
    #join_attr: simple test
    @test "join_multi_file_1" {
     result=`gtftk join_multi_file -i simple.gtf -k gene_id -t gene -m simple.join_mat simple.join_mat_2| gtftk select_by_key -g| grep G0003 | gtftk tabulate -k all -s "|"| tail -n 1`
      [ "$result" = "chr1|gtftk|gene|50|61|.|-|.|G0003|0.2322|0.4|A|B" ]
    }
    
    #join_attr: simple test
    @test "join_multi_file_2" {
     result=`gtftk join_multi_file -i simple.gtf -k gene_id  -t gene  -V 2 -m simple.join_mat simple.join_mat_2 simple.join_mat_3 | gtftk select_by_regexp  -k S5 -r "\d+"| gtftk tabulate -Hun -k S5,S6| perl -npe 's/\\t/_/g; s/\\n/;/'`
      [ "$result" = "0.2322_0.4;0.999|0.999_0.6|0.6;0.5555|20_0.7|30;" ]
    }
 
    """

    CmdObject(name="join_multi_file",
              message="Join attributes from mutiple files.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="editing",
              updated=__updated__,
              desc=__doc__,
              test=test)
