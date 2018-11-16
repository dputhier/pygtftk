#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO

__DESC__ = '''If several records share the same ID, print the first encountered, irrespective of the sequence.'''

def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True,
                                     description=__DESC__)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the FASTA file. Default to STDIN",
                            default=sys.stdin,
                            metavar="FASTA",
                            type=argparse.FileType('r'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Path to the output FASTA file. Default to STDOUT",
                            default=sys.stdout,
                            metavar="FASTA",
                            type=argparse.FileType('w'))


    return parser

def get_seq_from_ids(inputfile=None,
               outputfile=None):


    id_printed = set()

    with inputfile as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            if rec.id not in id_printed:
                outputfile.write(">" + rec.id + "\n")
                outputfile.write(str(rec.seq) + "\n")
                id_printed.add(rec.id)
                
            
def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_seq_from_ids(**args)


if __name__ == '__main__':
    main()
