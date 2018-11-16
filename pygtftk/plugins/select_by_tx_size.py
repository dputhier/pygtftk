#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys
from builtins import str

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import message

__updated__ = "2018-01-20"

__doc__ = """
 Select transcript based on their size (i.e size of mature/spliced transcript).
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

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

    parser.add_argument('-m', '--min-size',
                        help="Minimum size.",
                        default=0,
                        type=int)

    parser.add_argument('-M', '--max-size',
                        help="Maximum size.",
                        default=1000000000,
                        type=int)

    return parser


def select_by_tx_size(inputfile=None,
                      outputfile=None,
                      min_size=None,
                      max_size=None,
                      tmp_dir=None,
                      logger_file=None,
                      verbosity=0):
    """
    Select features by size.
    """

    msg = "Selecting 'mature/spliced transcript by size (range: [{m},{M}])."
    msg = msg.format(m=str(min_size),
                     M=str(max_size))
    message(msg)

    GTF(inputfile
        ).select_by_transcript_size(min_size,
                                    max_size
                                    ).write(outputfile,
                                            gc_off=True)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_tx_size(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #select_by_tx_size
    @test "select_by_tx_size_1" {
     result=`gtftk get_example | gtftk select_by_tx_size -m 14 | wc -l`
      [ "$result" -eq 6 ]
    }
    #select_by_tx_size
    @test "select_by_tx_size_2" {
     result=`gtftk get_example | gtftk select_by_tx_size -m 15 | wc -l`
      [ "$result" -eq 0 ]
    }
    #select_by_tx_size
    @test "select_by_tx_size_3" {
     result=`gtftk get_example | gtftk select_by_tx_size -m 0 | wc -l`
      [ "$result" -eq 60 ]
    }
    #select_by_tx_size
    @test "select_by_tx_size_4" {
     result=`gtftk get_example | gtftk feature_size -t mature_rna | gtftk select_by_tx_size -m 8 -M 8|wc -l`
      [ "$result" -eq 16 ]
    }
    """

    CmdObject(name="select_by_tx_size",
              message="Select transcript based on their size (i.e size of mature/spliced transcript).",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              desc=__doc__,
              updated=__updated__,
              test=test)
