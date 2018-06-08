#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys

from pygtftk.arg_formatter import FileWithExtension
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
                        type=FileWithExtension('r',
                                               valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser.add_argument('-o', '--outputfile',
                        help="Output file.",
                        default=sys.stdout,
                        metavar="GTF",
                        type=FileWithExtension('w',
                                               valid_extensions='\.[Gg][Tt][Ff]$'))

    parser.add_argument('-g', '--keep-gene-lines',
                        help="Add gene lines to the output",
                        action="store_true")
    return parser


def select_most_5p_tx(inputfile=None,
                      outputfile=None,
                      keep_gene_lines=False,
                      tmp_dir=None,
                      logger_file=None,
                      verbosity=0):
    """
    Select the most 5' transcript of each gene.
    """

    message("Selecting the most 5' transcript of each gene.")

    gtf = GTF(inputfile)

    if keep_gene_lines:
        gtf = gtf.select_5p_transcript()
    else:
        gtf = gtf.select_5p_transcript().select_by_key("feature", "gene", 1)

    gtf.write(outputfile)


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
    """

    CmdObject(name="select_most_5p_tx",
              message="Select the most 5' transcript of each gene.",
              parser=make_parser(),
              fun=select_most_5p_tx,
              updated=__updated__,
              group="selection",
              desc=__doc__,
              notes=__notes__,
              test=test)
