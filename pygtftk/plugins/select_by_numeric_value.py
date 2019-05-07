#!/usr/bin/env python
"""
 Select lines from a GTF file based on a boolean test on numeric values.
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
 -- The test should be of the form : "{key} operator value".
 -- e.g. "{transcript_version} > 10"
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

    parser.add_argument('-t', '--test',
                        help='The test to be applied.',
                        default=None,
                        type=str,
                        required=True)

    parser.add_argument('-n', '--na-omit',
                        help='If one of the evaluated values is enclosed in this list (csv), line is skipped.',
                        default=None,
                        type=str,
                        required=False)

    return parser


def select_by_numeric_value(inputfile=None,
                            outputfile=None,
                            test=None,
                            na_omit=None):
    """Select lines from a GTF file based on a boolean test on numeric values.
    """

    GTF(inputfile,
        check_ensembl_format=False).eval_numeric(test,
                                                 na_omit=na_omit,
                                                 ).write(outputfile,
                                                         gc_off=True)
    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_numeric_value(**args)


if __name__ == '__main__':
    main()


else:

    test = '''
    
    # select_by_numeric_value: load dataset
    @test "select_by_numeric_value_0" {
     result=`gtftk get_example -f '*' -d simple; gtftk get_example -d mini_real -f "*"`
      [ "$result" = "" ]
    }
    
    #select_by_numeric_value_1
    @test "select_by_numeric_value_1" {
      result=`gunzip -f airway_love.txt.gz`
      [ -s "airway_love.txt" ]
    }
   
    @test "select_by_numeric_value_2" {
      result=`gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m|  gtftk select_by_numeric_value -t 'start < 10 and end > 10 and S1 == 0.5555 and S2 == 0.007e2' -n ".,?"| wc -l`
      [ "$result" -eq 5 ]
    }
    
    
    @test "select_by_numeric_value_3" {
      result=`gtftk join_attr -k gene_id -t gene -j airway_love.txt -i mini_real.gtf.gz -m -V 3 | gtftk select_by_numeric_value -t 'GSM1275862 > 2000 and GSM1275870 == 2527' -n '.,?' | gtftk tabulate -k gene_name,gene_id,gene_biotype,GSM1275862,GSM1275870 | wc -l`
      [ "$result" -eq 2 ]
    }

    @test "select_by_numeric_value_4" {
      result=`gtftk join_attr -H  -k gene_id -t gene -j airway_love.txt -i mini_real.gtf.gz -m -V 3 | gtftk select_by_numeric_value -t 'GSM1275862 > 2000 and GSM1275870 > 100 and GSM1275875 == 2125 ' -n '.,?' | gtftk tabulate -k gene_name,gene_id,gene_biotype,GSM1275862,GSM1275870,GSM1275875 -H | cut -f1`
      [ "$result" = "GDI1"  ]
    }

    @test "select_by_numeric_value_5" {
      result=`gtftk join_attr -H  -k gene_id -t gene -j airway_love.txt -i mini_real.gtf.gz -m -V 3 | gtftk select_by_numeric_value -t 'GSM1275862 > 2000 and GSM1275870 > 100 and GSM1275875 <1000 ' -n '.,?' | gtftk tabulate -k gene_name,gene_id,gene_biotype,GSM1275862,GSM1275870,GSM1275875 -H | cut -f1`
      [ "$result" = "CRABP2"  ]
    }
    '''

    CmdObject(name="select_by_numeric_value",
              message="Select lines from a GTF file based on a boolean test on numeric values.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              desc=__doc__,
              updated=__updated__,
              test=test)
