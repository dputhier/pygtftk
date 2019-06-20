#!/usr/bin/env python
"""
 Add the tss number to each transcript (5'->3').
"""

import argparse
import os
import sys
from collections import defaultdict

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly, make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
- For each gene reports as a new key the tss number of each transcript.
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
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-d', '--dest-key',
                            help="The name of the new key.",
                            default='tss_number',
                            type=str,
                            required=False)

    return parser


def tss_numbering(
        inputfile=None,
        outputfile=None,
        dest_key='tss_number'):
    """
    Computes the distance between TSS of gene transcripts.
    """

    gtf = GTF(inputfile, check_ensembl_format=True)

    gn_tss_dist = defaultdict(dict)

    message("Getting TSSs.")
    tss = gtf.get_tss(name=["transcript_id", "gene_id"], as_dict=True)

    for k in tss:
        tx_id, gn_id = k.split("|")
        gn_tss_dist[gn_id][tx_id] = int(tss[k])

    gn_to_tx_to_tss = gtf.get_gn_to_tx(as_dict_of_dict=True)

    message("Numbering TSSs.")

    tss_number_file = make_tmp_file(prefix='tx_to_tss_number', suffix='.txt')

    for gn_id in gn_to_tx_to_tss:
        for tx_id in gn_to_tx_to_tss[gn_id]:
            tss_num = str(gn_to_tx_to_tss[gn_id][tx_id])
            tss_number_file.write(tx_id + "\t" + tss_num + "\n")

    tss_number_file.close()

    gtf = gtf.add_attr_from_file(feat='transcript',
                                 key='transcript_id',
                                 new_key=dest_key,
                                 inputfile=open(tss_number_file.name),
                                 has_header=False)

    gtf.write(outputfile,
              gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    tss_numbering(**args)


if __name__ == '__main__':
    main()

else:

    test = """
        
    #tss_dist 
    @test "tss_numbering_0" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -k gene_id -v ENSG00000175756 | gtftk tss_numbering| gtftk select_by_key -t| gtftk tabulate -x -k transcript_id,start,end,strand,tss_number | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = 'e01c4c1c167c8f2551584a3a6352b48b' ]
    }

    #tss_dist 
    @test "tss_numbering_1" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -k gene_id -v ENSG00000142611 | gtftk tss_numbering| gtftk select_by_key -t| gtftk tabulate -x -k transcript_id,start,end,strand,tss_number | md5sum-lite | perl -npe 's/\\s.*//'`
      [ "$result" = 'c4a9b29bd0385b57ed64cde6325b46fe' ]
    }


    """

    CMD = CmdObject(name="tss_numbering",
                    message=" Add the tss number to each transcript (5'->3').",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    notes=__notes__,
                    desc=__doc__,
                    group="annotation",
                    test=test)
