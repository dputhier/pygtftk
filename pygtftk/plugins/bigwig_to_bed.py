'''
Convert a bigwig to a BED3 format by selecting regions with coverage above --lower-val.
This tool is not part of pygtftk distribution (and thus is not supposed to be maintained).
'''

import argparse
import os
import pyBigWig
import sys

from pygtftk import arg_formatter
from pygtftk.utils import check_boolean_exprs


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True,
                                     description=__doc__)

    parser_grp = parser.add_argument_group('Arguments')

    # --------------------- Main arguments ----------------------------------- #

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the bigwig file",
                            default=None,
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('bigwig')),
                            required=False)

    parser_grp.add_argument('-o', '--outputfile',
                            help='The output BED file.',
                            default=sys.stdout,
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'),
                            required=False)

    parser_grp.add_argument('-c', '--chrom-list',
                            help='The list of chromosomes of interest.',
                            default=None,
                            type=str,
                            nargs='*',
                            required=False)

    parser_grp.add_argument('-e', '--expression',
                            help='A boolean expression where s is the signal (e.g; s > 0 and s < 0.5).',
                            default=None,
                            type=str,
                            required=True)

    return parser


def bigwig_to_bed(inputfile=None,
                  outputfile=None,
                  expression=None,
                  chrom_list=None):
    '''
    Convert a bigwig to a BED3 format.
    '''

    # -------------------------------------------------------------------------
    # Check chromosomes
    # -------------------------------------------------------------------------

    bwig = pyBigWig.open(inputfile.name)

    chr_in_bw = list(bwig.chroms().keys())

    if chrom_list is None:
        chrom_list = chr_in_bw

    else:
        for i in chrom_list:
            if i not in chr_in_bw:
                sys.stderr.write('Chromosome:' + i + " not found. Please check.\n")
                sys.exit()

    # -------------------------------------------------------------------------
    # Check boolean exprs
    # -------------------------------------------------------------------------
    check_boolean_exprs(expression, operand=['s'])

    # -------------------------------------------------------------------------
    # Loop through chroms
    # -------------------------------------------------------------------------

    exp_as_func = eval('lambda s: ' + expression)

    for chrom in chrom_list:
        start = None
        end = None
        sys.stderr.write("Processing chromosome:" + chrom + "\n")
        for interval in bwig.intervals(chrom):
            s = interval[2]
            if start is None:
                if exp_as_func(s):
                    start = interval[0]
                    end = interval[1]
            else:
                if exp_as_func(s):
                    end = interval[1]
                else:
                    a_list = [chrom, str(start), str(end)]
                    a_str = "\t".join(a_list)
                    outputfile.write(a_str + "\n")
                    start = None


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    bigwig_to_bed(**args)


if __name__ == '__main__':
    main()

else:
    test = """
    #bigwig_to_bed
    @test "bigwig_to_bed_0" {
     result=`gtftk get_example -f '*' `
      [ "$result" = "" ]
    }

    #bigwig_to_bed
    @test "bigwig_to_bed_1" {
     result=`gtftk bigwig_to_bed -i simple.2.bw -e 's > 2' | wc -l`
      [ "$result" -eq 7 ]
    }
    
    #bigwig_to_bed
    @test "bigwig_to_bed_2" {
     result=`gtftk bigwig_to_bed -i simple.2.bw -e 's > 2 and  s< 4' | wc -l`
      [ "$result" -eq 10 ]
    }
        
    
    
    """
    from pygtftk.cmd_object import CmdObject

    CmdObject(name="bigwig_to_bed",
              message="Convert a bigwig to a BED3 format.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              group="conversion",
              test=test)
