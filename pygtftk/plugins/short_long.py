# -*- coding: utf-8 -*-


""" Select the shortest mature transcript (i.e without introns) for each gene or the longest if the \
-l arguments is used. """

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF

__doc__ = """ Select the shortest mature transcript (i.e without introns) for each gene or the longest if the -l arguments is used. """

__updated__ = "2018-01-25"

__notes__ = ""


def make_parser():
    """parse"""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Argument')

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

    parser_grp.add_argument('-l', '--longs',
                            help="Take the longest transcript of each gene",
                            action="store_true")

    parser_grp.add_argument('-g', '--keep-gene-lines',
                            help="Add gene lines to the output",
                            action="store_true")

    return parser


def short_long(inputfile=None,
               outputfile=None,
               longs=None,
               keep_gene_lines=False):
    """ Select the shortest transcript for each gene, Or the longuest if the \
-l arguments is used. """

    gtf = GTF(inputfile, check_ensembl_format=False)

    if longs:
        gtf = gtf.select_longuest_transcripts()
    else:
        gtf = gtf.select_shortest_transcripts()

    if not keep_gene_lines:
        gtf = gtf.select_by_key("feature", "gene", 1)

    gtf.write(outputfile,
              gc_off=True)


def main():
    """main"""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    short_long(**args)


if __name__ == '__main__':
    main()

else:
    test = """

    #short_long: load the dataset
    @test "short_long_0" {
     result=`gtftk get_example -f '*' -d simple_03`
      [ "$result" = "" ]
    }
        
    #Test number of output lines
    @test "short_long_1" {
    result=$(gtftk short_long -i simple_short_long.gtf  -g | gtftk select_by_key -t | wc -l)
    [ $result -eq 11 ]
    }
    
    #Checks that the shortest transcripts have been selected (short)
    @test "short_long_2" {
    result=$(gtftk short_long -i simple_short_long.gtf| grep -c -E 'G0001T002|G0002T002|G0003T002|G0006T002|G0008T002|G0011T001')
    [ $result -eq 22 ]
    }
    
    #Checks that the transcripts With the same tss are not present (short)
    @test "short_long_3" {
    result=$(gtftk short_long -i simple_short_long.gtf | perl -ne 'print if (/(G0001T001)|(G0002T001)|(G0003T001)|(G0006T001)|(G0008T001)|(G0011T002)/)')
    [ -z $result ]
    }
    
    #Checks that the longest transcripts have been selected (long)
    @test "short_long_4" {
    result=$(gtftk short_long -i simple_short_long.gtf -l | grep -c -E 'G0001T001|G0002T001|G0003T001|G0006T001|G0008T001|G0011T002')
    [ $result -eq 27 ]
    }
    
    #Checks that the transcripts With the same tss are not present (long)
    @test "short_long_5" {
    result=$(gtftk short_long -i simple_short_long.gtf -l | perl -ne 'print if (/(G0001T002)|(G0002T002)|(G0003T002)|(G0006T002)|(G0008T002)|(G0011T001)/)')
    [ -z $result ]
    }


    #Test number of output lines (genes)
    @test "short_long_6" {
    result=$(gtftk short_long -i simple_short_long.gtf  -g |  gtftk select_by_key -k feature -v gene| wc -l)
    [ $result -eq 11 ]
    }

    #Test number of output lines (genes)
    @test "short_long_7" {
    result=$(gtftk get_example -d mini_real  | gtftk feature_size -t mature_rna | gtftk short_long -l |  gtftk select_by_key -k gene_name -v ISG15 | gtftk select_by_key -t |  gtftk tabulate -H -k feat_size)
    [ $result -eq 788 ]
    }
    
    #Test number of output lines (genes)
    @test "short_long_8" {
    result=$(gtftk get_example -d mini_real  | gtftk feature_size -t mature_rna | gtftk short_long  |  gtftk select_by_key -k gene_name -v ISG15 | gtftk select_by_key -t |  gtftk tabulate -H -k feat_size)
    [ $result -eq 657 ]
    }


    #Test number of output lines (genes)
    @test "short_long_9" {
    result=$(gtftk get_example -d mini_real  | gtftk short_long -l |  gtftk feature_size -t mature_rna |  gtftk select_by_key -k gene_name -v AURKAIP1 | gtftk select_by_key --select-transcripts | gtftk tabulate -Hun -k feat_size)
    [ $result -eq 1072 ]
    }
    
    #Test number of output lines (genes)
    @test "short_long_10" {
    result=$(gtftk get_example -d mini_real  | gtftk short_long  |  gtftk feature_size -t mature_rna |  gtftk select_by_key -k gene_name -v AURKAIP1 | gtftk select_by_key --select-transcripts | gtftk tabulate -Hun -k feat_size)
    [ $result -eq 608 ]
    }
    
    """

    CmdObject(name="short_long",
              message="Get the shortest or longest transcript of each gene",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              desc=__doc__,
              updated=__updated__,
              notes=__notes__,
              test=test)
