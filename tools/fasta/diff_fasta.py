#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.SeqIO import FastaIO

__DESC__ = ''' Compare sequences from two fasta files with records ordered in the same way.'''

def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True, description=__DESC__)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-1', '--inputfile-1',
                            help="Path to the FASTA file.",
                            default=None,
                            metavar="FASTA",
                            type=argparse.FileType('r'),
                            required=True)

    parser_grp.add_argument('-2', '--inputfile-2',
                            help="Path to the FASTA file.",
                            default=None,
                            metavar="FASTA",
                            type=argparse.FileType('r'),
                            required=True)

    parser_grp.add_argument('-o', '--outputfile',
                            help="Path to the output FASTA file. Default to STDOUT",
                            default=sys.stdout,
                            metavar="FASTA",
                            type=argparse.FileType('w'))

    return parser

def get_seq_from_ids(inputfile_1=None,
                     inputfile_2=None,
                     outputfile=None):


    record_dict_1 = SeqIO.to_dict(SeqIO.parse(inputfile_1.name, "fasta"))
    record_dict_2 = SeqIO.to_dict(SeqIO.parse(inputfile_2.name, "fasta"))

    all_key = list(set(list(record_dict_1.keys()) + list(record_dict_2.keys())))
    for i in all_key:
        if record_dict_1[i].seq != record_dict_2[i].seq:
            print("-- " + i + " has not the same sequence in both files...")
            print("->" + record_dict_1[i].seq)
            print("->" + record_dict_2[i].seq)                  
        else:
            print("-- " + i + " has same sequence")

def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_seq_from_ids(**args)


if __name__ == '__main__':
    main()
