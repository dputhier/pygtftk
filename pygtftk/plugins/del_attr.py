#!/usr/bin/env python

import argparse
import os
import re
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
 Delete one or several attributes from the gtf file.
"""
__notes__ = """
 -- You may also use 'complex' regexp such as : "(^.*_id$|^.*_biotype$)"
 -- Example: gtftk get_example -d mini_real | gtftk del_attr -k "(^.*_id$|^.*_biotype$)" -r -v
"""


def make_parser():
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
                            help='comma-separated list of attribute names or a '
                                 'regular expression (see -r).',
                            default=None,
                            metavar="KEY",
                            type=str,
                            required=True)

    parser_grp.add_argument('-r', '--reg-exp',
                            help='The key name is a regular expression.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-v', '--invert-match',
                            help='Delected keys are those not matching any of'
                                 ' the specified key.',
                            action="store_true",
                            required=False)

    return parser


def del_attr(
        inputfile=None,
        outputfile=None,
        key="transcript_id",
        reg_exp=False,
        invert_match=False):
    """
    Delete extended attributes in the target gtf file. attr_list can be a
    comma-separated list of attributes.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)

    if reg_exp:
        try:
            rgxp = re.compile(key)
        except:
            message("Check the regular expression please.", type="ERROR")
        key_list = [key]
    else:
        key_list = key.split(",")

    for i in gtf:

        feature_keys = i.get_attr_names()

        if not invert_match:
            for k in key_list:
                if not reg_exp:
                    try:
                        del i.attr[k]
                    except:
                        pass
                else:
                    for feat_key in feature_keys:
                        if rgxp.search(feat_key):
                            del i.attr[feat_key]
        else:

            for k in feature_keys:
                if not reg_exp:
                    if k not in key_list:
                        del i.attr[k]
                else:
                    if not rgxp.search(k):
                        del i.attr[k]

        i.write(outputfile)

    close_properly(outputfile, inputfile)


if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    del_attr(**args)

else:

    test = """

    # del_attr: load dataset
    @test "del_attr_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
    
    #del_attr:
    # If you delete almost all extended attributes there are only exon_id left
    @test "del_attr_1" {
     result=`gtftk del_attr -i simple.gtf  -k ccds_id,transcript_id,gene_id| cut -f9| grep -v "^$"| sed 's/ \".*//'| sort | uniq`
      [ "$result" = "exon_id" ]
    }
    
    
    #del_attr: check -v
    @test "del_attr_2" {
     result=`gtftk del_attr -i simple.gtf  -k ccds_id,transcript_id,gene_id -v| grep exon_id| wc -l`
      [ "$result" -eq 0 ]
    }

    #del_attr: check -r
    @test "del_attr_3" {
     result=`gtftk get_example -d mini_real | gtftk del_attr -r -k 'transcript.*'| gtftk tabulate -k "*" -x -V 3| head -1| grep transcript | wc -l`
      [ "$result" -eq 0 ]
    }
        
 
    #del_attr: check -r
    @test "del_attr_4" {
     result=`gtftk get_example -d mini_real | gtftk del_attr -r -k '(trancrip)|(biotype)|(exon_id)'| gtftk tabulate -k "*" -x| head -1| awk 'BEGIN{FS="\\t"}{print NF}'`
      [ "$result" -eq 13 ]
    }
        
           
    #del_attr: check -r
    @test "del_attr_5" {
     result=`gtftk get_example -d mini_real | gtftk del_attr -r -k '(transcript_id)|(gene_id)' -v| awk 'BEGIN{FS="\\t"}{print NF}' | sort | uniq -c | cut -f2 | sed 's/ //g'`
      [ "$result" = "1376709" ]
    }
 
    """

    CmdObject(name="del_attr",
              message="Delete attributes in the target gtf file.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              notes=__notes__,
              desc=__doc__,
              group="editing",
              test=test)
