#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.SeqIO import FastaIO

__DESC__ = '''Take a FASTA as input and produced a new FASTA file with record 
ordered as in --id-file.'''

def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True, description=__DESC__)

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

def unfold(inputfile=None,
                     outputfile=None):

    for rec in SeqIO.parse(inputfile.name, "fasta"):
            SeqIO.write(rec, outputfile, "fasta-2line")

def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    unfold(**args)


if __name__ == '__main__':
    main()
