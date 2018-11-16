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

    parser_grp.add_argument('-f', '--id-file',
                            help="The file containing the IDS to be extracted.",
                            default=None,
                            metavar="FASTA",
                            type=argparse.FileType('r'),
                            required=True)

    return parser

def get_seq_from_ids(inputfile=None,
               outputfile=None,
               id_file=None):

    id_list = id_file.readlines()
    id_list = [x.replace(">","").rstrip("\n") for x in id_list]

    record_dict = SeqIO.to_dict(SeqIO.parse(inputfile.name, "fasta"))

    fasta_out = FastaIO.FastaWriter(outputfile, wrap=None)
    
    for i in id_list:
        print(">" + record_dict[i].id)
        print(record_dict[i].seq)

def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_seq_from_ids(**args)


if __name__ == '__main__':
    main()
