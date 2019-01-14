#!/usr/bin/env python


import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import message

__updated__ = "2018-01-20"

__doc__ = """
 Select the most 5' transcript of each gene.
"""

__notes__ = """
  -- If several transcript share the samemost 5' TSS, only one transcript is
  selected.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser()

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

    parser.add_argument('-g', '--keep-gene-lines',
                        help="Add gene lines to the output",
                        action="store_true")
    return parser


def select_most_5p_tx(inputfile=None,
                      outputfile=None,
                      keep_gene_lines=False):
    """
    Select the most 5' transcript of each gene.
    """

    message("Selecting the most 5' transcript of each gene.")

    gtf = GTF(inputfile)

    if keep_gene_lines:
        gtf = gtf.select_5p_transcript()
    else:
        gtf = gtf.select_5p_transcript().select_by_key("feature", "gene", 1)

    gtf.write(outputfile, gc_off=True)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_most_5p_tx(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #select_most_5p_tx
    @test "select_most_5p_tx_1" {
     result=`gtftk get_example -d simple_04 | gtftk select_most_5p_tx| wc -l`
      [ "$result" -eq 41 ]
    }
    
    #select_most_5p_tx
    @test "select_most_5p_tx_2" {
     result=`gtftk get_example -d simple_04 | gtftk select_most_5p_tx -g | wc -l`
      [ "$result" -eq 51 ]
    }

    #select_most_5p_tx
    @test "select_most_5p_tx_3" {
     result=`gtftk get_example -d mini_real  | gtftk select_most_5p_tx |  gtftk select_by_key -k gene_name -v ISG15 | gtftk select_by_key --select-transcripts | gtftk get_5p_3p_coords| cut -f4`
      [ "$result" = "ENSG00000187608|ENST00000624697" ]
    }

    #select_most_5p_tx
    @test "select_most_5p_tx_4" {
     result=`gtftk get_example -d mini_real  | gtftk select_most_5p_tx | gtftk select_by_key -k gene_name -v CRABP2 | gtftk select_by_key --select-transcripts | gtftk get_5p_3p_coords  -t transcript | cut -f4`
      [ "$result" = "ENSG00000143320|ENST00000621784" ]
    }        
    
    #select_most_5p_tx
    @test "select_most_5p_tx_5" {
     result=`gtftk get_example -d mini_real  | gtftk select_most_5p_tx |  gtftk select_by_key -k gene_name -v RAB13 | gtftk select_by_key --select-transcripts | gtftk get_5p_3p_coords  -t transcript  | cut -f 4`
      [ "$result" = "ENSG00000143545|ENST00000495720" ]
    }  


    #select_most_5p_tx
    @test "select_most_5p_tx_6" {
     result=`gtftk get_example -d mini_real | gtftk select_most_5p_tx | gtftk select_by_key -t | gtftk tabulate -k gene_id,gene_name,transcript_id | perl -ne 'print if (/(LRRFIP1)|(MPP1)|(IGSF1)|(SIN3B)|(PITPNM2)|(XIST)/)' | wc -l`
      [ "$result" -eq 6 ]
    } 
    

    #select_most_5p_tx
    @test "select_most_5p_tx_7" {
     result=`gtftk get_example -d mini_real | gtftk select_most_5p_tx | gtftk select_by_key -t | gtftk tabulate -k gene_id,gene_name,transcript_id | perl -ne 'print if (/(LRRFIP1)|(MPP1)|(IGSF1)|(SIN3B)|(PITPNM2)|(XIST)/)' | cut -f 3 | perl -npe 's/\\n/,/'`
      [ "$result" = "ENST00000542210,ENST00000248054,ENST00000308482,ENST00000429829,ENST00000370904,ENST00000393529," ]
    }  

    #select_most_5p_tx
    @test "select_most_5p_tx_8" {
     result=`gtftk get_example -d simple_04 | gtftk select_most_5p_tx | grep "G0005T002"| wc -l`
      [ "$result" -eq 5 ]
    } 

    #select_most_5p_tx
    @test "select_most_5p_tx_9" {
     result=`gtftk get_example -d simple_04 | gtftk select_most_5p_tx | grep "G0004T001"| wc -l`
      [ "$result" -eq 6 ]
    } 
    
        
    
    
    """

    CmdObject(name="select_most_5p_tx",
              message="Select the most 5' transcript of each gene.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              group="selection",
              desc=__doc__,
              notes=__notes__,
              test=test)
