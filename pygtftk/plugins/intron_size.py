#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-24"
__doc__ = """
 Add a new key to transcript features containing a comma-separated list of intron-size.
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
                            help="Output file.",
                            default=sys.stdout,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-a',
                            '--key-name',
                            type=str,
                            default="intron_sizes",
                            help="The name of the key.",
                            required=False)

    return parser


def intron_sizes(
        inputfile=None,
        outputfile=None,
        key_name=None):
    """
 Add a new key to transcript features containing a comma-separated list of intron sizes.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)

    all_tx_ids = gtf.get_tx_ids(nr=True)
    intron_bo = gtf.get_introns(by_transcript=True,
                                name=["transcript_id"],
                                intron_nb_in_name=False,
                                feat_name=False)

    strands = gtf.select_by_key("feature",
                                "transcript").extract_data("transcript_id,strand",
                                                           as_dict_of_values=True,
                                                           no_na=True,
                                                           nr=True,
                                                           hide_undef=True)

    intron_size = {tx: [] for tx in all_tx_ids}

    for bed_line in intron_bo:
        intron_size[bed_line.name] += [str(bed_line.end - bed_line.start)]

    for tx_id in intron_size:
        if len(intron_size[tx_id]):
            if strands[tx_id] == "-":
                intron_size[tx_id] = ",".join(reversed(intron_size[tx_id]))
            else:
                intron_size[tx_id] = ",".join(intron_size[tx_id])
        else:
            intron_size[tx_id] = "0"
    if len(intron_size):
        gtf = gtf.add_attr_from_dict(feat="transcript",
                                     key="transcript_id",
                                     a_dict=intron_size,
                                     new_key=key_name)
    gtf.write(outputfile,
              gc_off=True)
    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    intron_sizes(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #intron_sizes
    @test "intron_sizes_1" {
     result=`gtftk get_example -d simple_04 | sortBed | gtftk intron_sizes| gtftk select_by_key -t| gtftk tabulate -k transcript_id,intron_sizes -Hun| grep G0004T001 | cut -f2`
      [ "$result" = "3,2,2" ]
    }
    
    #intron_sizes
    @test "intron_sizes_2" {
     result=`gtftk get_example -d simple_04 | sortBed | gtftk intron_sizes| gtftk select_by_key -t| gtftk tabulate -k transcript_id,intron_sizes -Hun| grep G0005T002 | cut -f2`
      [ "$result" = "2,6" ]
    }
    
    #intron_sizes
    @test "intron_sizes_3" {
     result=`gtftk get_example -d mini_real | sortBed | gtftk intron_sizes | gtftk exon_sizes |   gtftk intron_sizes | gtftk select_by_key -t| grep --color ENST00000475943| gtftk tabulate -k intron_sizes -Hun`
      [ "$result" = "18871,1046" ]
    }
    """

    CMD = CmdObject(name="intron_sizes",
                    message=" Add a new key to transcript features containing a comma-separated list of intron sizes. ",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="annotation",
                    test=test)
