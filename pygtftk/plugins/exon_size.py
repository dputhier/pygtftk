#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-24"
__doc__ = """
 Add a new key to transcript features containing a comma-separated list of exon sizes.
"""
__notes__ = """
 -- The GTF should be sorted before computation. Use bedtools sortBed.
 -- Sizes are provided in 5'->3' orientation (iff the GTF was previously sorted).
"""


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN.",
                            default=sys.stdin,
                            metavar="GTF",
                            required=False,
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output GTF file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-a',
                            '--key-name',
                            type=str,
                            default="exon_sizes",
                            help="The name of the key.",
                            required=False)

    return parser


def exon_sizes(
        inputfile=None,
        outputfile=None,
        key_name=None):
    """
 Add a new key to transcript features containing a comma-separated list of exon-size.
    """

    gtf = GTF(inputfile)

    all_tx_ids = gtf.get_tx_ids(nr=True)
    tx_to_size_list = dict()
    exons_starts = gtf.select_by_key("feature",
                                     "exon").extract_data("transcript_id,start",
                                                          as_dict_of_merged_list=True,
                                                          no_na=True,
                                                          nr=False)

    if not len(exons_starts):
        message("No exon found.", type="ERROR")

    exons_ends = gtf.select_by_key("feature",
                                   "exon").extract_data("transcript_id,end",
                                                        as_dict_of_merged_list=True,
                                                        no_na=True,
                                                        nr=False)

    strands = gtf.select_by_key("feature",
                                "transcript").extract_data("transcript_id,strand",
                                                           as_dict_of_values=True,
                                                           no_na=True,
                                                           nr=True,
                                                           hide_undef=True)

    for tx_id in all_tx_ids:
        size_list = []
        for s, e in zip(exons_starts[tx_id], exons_ends[tx_id]):
            size = str(int(e) - int(s) + 1)
            size_list += [size]
        if strands[tx_id] == "-":
            size_list = reversed(size_list)
        tx_to_size_list[tx_id] = ",".join(size_list)

    if len(tx_to_size_list):
        gtf = gtf.add_attr_from_dict(feat="transcript",
                                     key="transcript_id",
                                     a_dict=tx_to_size_list,
                                     new_key=key_name)
    gtf.write(outputfile,
              gc_off=True)
    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    exon_sizes(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #exon_sizes
    @test "exon_sizes_1" {
     result=`gtftk get_example  | sortBed | gtftk exon_sizes| gtftk tabulate -k transcript_id,exon_sizes -Hun| grep G0004T001| cut -f2`
      [ "$result" = "4,1,3" ]
    }
    
    #exon_sizes
    @test "exon_sizes_2" {
     result=`gtftk get_example  | sortBed |  gtftk exon_sizes| gtftk tabulate -k transcript_id,exon_sizes -Hun| grep G0006T001| cut -f2`
      [ "$result" = "3,3,4" ]
    }

    #exon_sizes
    @test "exon_sizes_3" {
     result=`gtftk get_example -d mini_real | sortBed | gtftk intron_sizes | gtftk exon_sizes |   gtftk select_by_key -t| grep --color ENST00000475943| gtftk tabulate -k exon_sizes -Hun`
      [ "$result" = "206,197,457" ]
    }
    
    """

    CMD = CmdObject(name="exon_sizes",
                    message=" Add a new key to transcript features containing a comma-separated list of exon sizes. ",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="annotation",
                    test=test)
