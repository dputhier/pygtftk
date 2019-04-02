#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message, make_tmp_file

__updated__ = "2018-01-20"
__doc__ = """
Compute the number of transcript per gene.
"""


def make_parser():
    """The program parser."""
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
                            metavar="GTF/TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf', 'txt')))

    parser_grp.add_argument('-f', '--text-format',
                            help="Return a text format.",
                            action="store_true")

    parser_grp.add_argument('-a',
                            '--key-name',
                            type=str,
                            default="nb_tx",
                            help="If gtf format is requested, the name of the key.",
                            required=False)

    return parser


def nb_transcripts(inputfile=None,
                   outputfile=None,
                   text_format=False,
                   key_name=""):
    """
    Compute the number of transcript per gene.
    """

    gtf = GTF(inputfile)

    message("Computing the number of transcript per gene in input GTF file.")

    # Computation of transcript number is performed on exon lines
    # Just in case some transcript lines would be lacking (but they should
    # not...)

    n_tx = gtf.get_gn_to_tx()

    if not text_format:
        tmp_file = make_tmp_file(prefix="nb_tx", suffix=".txt")

    for i in n_tx:
        if not text_format:
            tmp_file.write(i + "\t" + str(len(n_tx[i])) + "\n")
        else:
            outputfile.write(i + "\t" + str(len(n_tx[i])) + "\n")

    if not text_format:
        tmp_file.close()
        gtf.add_attr_from_file(feat="gene",
                               key="gene_id",
                               new_key=key_name,
                               inputfile=tmp_file.name).write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    nb_transcripts(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #nb_transcripts:
    @test "nb_transcripts_1" {
     result=`gtftk get_example| bedtools sort |   gtftk  nb_transcripts  | gtftk select_by_key -k feature -v gene | gtftk tabulate -k gene_id,nb_tx| sort -nk2,2r| perl -npe 's/\\t/,/g; s/\\n/,/g'`
      [ "$result" = "gene_id,nb_tx,G0001,2,G0004,2,G0006,2,G0007,2,G0009,2,G0002,1,G0003,1,G0005,1,G0008,1,G0010,1," ]
    }

    #nb_transcripts:
    @test "nb_transcripts_2" {
     result=`gtftk get_example| bedtools sort |   gtftk nb_transcripts -f | perl -npe 's/\\t/,/g; s/\\n/,/g'`
      [ "$result" = "G0009,2,G0006,2,G0005,1,G0003,1,G0004,2,G0007,2,G0001,2,G0010,1,G0002,1,G0008,1," ]
    }
    
   """

    CmdObject(name="nb_transcripts",
              message="Count the number of transcript per gene.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="information",
              updated=__updated__,
              desc=__doc__,
              test=test)
