#!/usr/bin/env python
"""
 Join attributes from a tabulated file.
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

    parser_grp.add_argument('-j', '--join-file',
                            help="A two columns file with (i) the value for joining "
                                 "(e.g value for transcript_id), (ii) the value for novel "
                                 "key (e.g the coding potential computed value).",
                            default=None,
                            metavar="JOIN_FILE",
                            type=argparse.FileType("r"),
                            required=True)

    parser_grp.add_argument('-H',
                            '--has-header',
                            help="Indicates that the 'join-file' has a header.",
                            action="store_true")

    parser_grp.add_argument('-m',
                            '--matrix',
                            help="'join-file' expect a "
                                 "matrix with row names as target keys column names as novel "
                                 "key and each cell as value.",
                            action="store_true")

    parser_grp.add_argument('-n', '--new-key',
                            help='The name of the novel key. Mutually exclusive with --matrix.',
                            default=None,
                            metavar="NEW_KEY",
                            type=str,
                            required=False)

    parser_grp.add_argument('-t', '--target-feature',
                            help='The name(s) of the target feature(s). comma-separated.',
                            default=None,
                            type=str,
                            required=False)

    return parser


def join_attr(
        inputfile=None,
        outputfile=None,
        join_file=None,
        has_header=False,
        new_key=None,
        target_feature=None,
        key_to_join=None,
        matrix=None):
    """
    Join attributes from a tabulated file.
    """

    # -----------------------------------------------------------
    #  Check argument consistency
    # -----------------------------------------------------------

    if matrix is True:
        if new_key is not None:
            message("--new-key and --matrix are mutually exclusive.",
                    type="ERROR")
    else:
        if new_key is None:
            message("--new-key is required when --matrix is False.",
                    type="ERROR")

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

    if not matrix:

        gtf = gtf.add_attr_from_file(feat=target_feature,
                                     key=key_to_join,
                                     new_key=new_key,
                                     inputfile=join_file.name,
                                     has_header=has_header)
        gtf.write(outputfile,
                  gc_off=True)

    else:

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
    join_attr(**args)


if __name__ == '__main__':
    main()
else:

    test = """

    # join_attr: load dataset
    @test "join_attr_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
       
       
    #join_attr: simple test
    @test "join_attr_1" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join -k gene_id -n bla| gtftk tabulate -H -k gene_id,bla|grep G0009 | sort | uniq | cut -f2`
      [ "$result" = "0.5555" ]
    }
    
    
    #join_attr: simple test
    @test "join_attr_2" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join -k gene_id -n bla| gtftk tabulate -H -k gene_id,bla|grep G0009 | sort | uniq | cut -f2`
      [ "$result" = "0.5555" ]
    }
    
    #join_attr: simple test
    @test "join_attr_3" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join -k gene_id -n bla| wc -l`
      [ "$result" -eq 70 ]
    }
    
    
    #join_attr: test matrix input
    @test "join_attr_4" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m| wc -l`
      [ "$result" -eq 70 ]
    }
    
    #join_attr: simple test
    @test "join_attr_5" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m | grep 0.5555|gtftk  tabulate -H -k gene_id| head -1`
      [ "$result" = "G0009" ]
    }

    #join_attr: simple test
    @test "join_attr_6" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id  -V 2 -m -t exon| grep 0.999| wc -l`
      [ "$result" -eq 6 ]
    }
        
    #join_attr: simple test
    @test "join_attr_7" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id  -V 2 -m -t exon,transcript| grep 0.999| wc -l`
      [ "$result" -eq 8 ]
    }

    #join_attr: simple test
    @test "join_attr_8" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id  -V 2 -m -t exon,transcript,gene| grep 0.999| wc -l`
      [ "$result" -eq 9 ]
    }

    #join_attr: simple test
    @test "join_attr_9" {
     result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id  -V 2 -m -t exon,transcript,gene,CDS | grep 0.999| wc -l`
      [ "$result" -eq 13 ]
    }

    #join_attr: simple test
    @test "join_attr_10" {
     result=`gtftk get_example |  gtftk join_attr -j simple.join -k gene_id -n bla -t gene  -V 2 -K toto| gtftk select_by_regexp -k bla -r "0\..*" | gtftk tabulate -k bla -Hun | perl -npe 's/\\n/|/g'`
      [ "$result" = "0.2322|0.999|0.5555|" ]
    }
    
    #join_attr: simple test
    @test "join_attr_11" {
     result=`gtftk get_example |  gtftk join_attr -j simple.join_with_dup -k gene_id -n bla -t gene  -V 2 | gtftk select_by_regexp -k bla -r "0\."| gtftk tabulate -k bla -Hun | perl -npe 's/\\n/|/'`
      [ "$result" = "0.2322|0.2|0.999|0.5555|0.1|" ]
    }
        

    
    """
    CmdObject(name="join_attr",
              message="Join attributes from a tabulated file.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="editing",
              updated=__updated__,
              desc=__doc__,
              test=test)
