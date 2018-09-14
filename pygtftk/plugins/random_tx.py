#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import random
import sys

from builtins import range

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import int_greater_than_null
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import PY2
from pygtftk.utils import PY3
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-30"
__doc__ = """
 Select randomly up to m transcript for each gene.
"""


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=FileWithExtension('r',
                                                   valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="GTF",
                            type=FileWithExtension('w',
                                                   valid_extensions='\.[Gg][Tt][Ff]$'))

    parser_grp.add_argument('-m', '--max-transcript',
                            help="The maximum number of transcripts to select for each gene.",
                            default=1,
                            metavar="MAX",
                            type=int_greater_than_null,
                            required=False)

    parser_grp.add_argument('-s', '--seed-value',
                            help="Seed value for the random number generator.",
                            default=None,
                            metavar="SEED",
                            type=int_greater_than_null,
                            required=False)

    return parser


def random_tx(
        inputfile=None,
        outputfile=None,
        max_transcript=None,
        seed_value=None,
        tmp_dir=None,
        logger_file=None,
        verbosity=0):
    """
    Select randomly up to m transcript for each gene.
    """

    message("loading the GTF.")

    gtf = GTF(inputfile).select_by_key("feature",
                                       "gene",
                                       invert_match=True
                                       )

    message("Getting gene_id and transcript_id")

    gene2tx = gtf.extract_data("gene_id,transcript_id",
                               as_dict_of_merged_list=True,
                               no_na=True,
                               nr=True)

    message("Selecting random transcript")

    if seed_value is not None:
        if PY3:
            # Unfortunatly 'version 1'
            # is not back-compatible...
            random.seed(seed_value, version=1)
        elif PY2:
            random.seed(seed_value)

    tx_to_delete = []

    for gn_id in gene2tx:
        tx_list = gene2tx[gn_id]
        nb_tx = len(tx_list)
        max_cur = min(max_transcript, nb_tx)
        pos_to_keep = random.sample(list(range(len(tx_list))), max_cur)
        tx_list = [j for i, j in enumerate(tx_list) if i not in pos_to_keep]
        tx_to_delete += tx_list

    message("Printing results")

    message("Selecting transcript.")
    gtf = gtf.select_by_key("transcript_id",
                            ",".join(tx_to_delete),
                            invert_match=True
                            ).write(outputfile,
                                    gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    random_tx(**args)


if __name__ == '__main__':
    main()

else:

    test = '''
          
    #random_tx: random_tx should return 1 tx per gene (10 genes)
    @test "random_tx_1" {
     result=`gtftk random_tx -i pygtftk/data/simple/simple.gtf | gtftk select_by_key -k feature -v transcript | gtftk tabulate -k gene_id,transcript_id -H | wc -l`
      [ "$result" -eq 10 ]
    }
    
    #random_tx: this test should return 10
    @test "random_tx_3" {
     result=`gtftk get_example -d mini_real | gtftk random_tx -m 10 -s 111 | gtftk convert_ensembl |  gtftk nb_transcripts | gtftk select_by_key -g | gtftk tabulate -k gene_id,nb_tx -Hun | cut -f 2 | sort -n | tail -n 1`
      [ "$result" -eq 10 ]
    }

    #random_tx: this test should return 20
    @test "random_tx_4" {
     result=`gtftk get_example -d mini_real | gtftk random_tx -m 20  -s 111 | gtftk convert_ensembl |  gtftk nb_transcripts | gtftk select_by_key -g | gtftk tabulate -k gene_id,nb_tx -Hun | cut -f 2 | sort -n | tail -n 1`
      [ "$result" -eq 20 ]
    }

    #random_tx: 
    @test "random_tx_5" {
     result=`gtftk get_example -d mini_real | gtftk random_tx -m 10 -s 111 | gtftk convert_ensembl |  gtftk nb_exons | gtftk select_by_key -t | gtftk tabulate -k transcript_id,nb_exons -Hun -s "," | sort -n -k2,2 -t","| tail -n 1`
      [ "$result" = "ENST00000378016,107" ]
    }
    
    #random_tx: check -m
    @test "random_tx_6" {
     result=`gtftk get_example -d mini_real |  gtftk random_tx  -m 3 -s 111 | gtftk select_by_key -t | gtftk tabulate -k gene_id -H | sort | uniq -c| perl -npe 's/^ +//; s/ /\\t/' | cut -f1 | sort | uniq | sort -n | perl -npe 's/\\n/,/g'`
      [ "$result" = "1,2,3," ]
    }

    #random_tx: check -m
    @test "random_tx_7" {
     result=`gtftk get_example -d mini_real |  gtftk random_tx  -m 4 -s 111 | gtftk select_by_key -t | gtftk tabulate -k gene_id -H | sort | uniq -c| perl -npe 's/^ +//; s/ /\\t/' | cut -f1 | sort | uniq | sort -n | perl -npe 's/\\n/,/g'`
      [ "$result" = "1,2,3,4," ]
    }     
     
    '''

    CMD = CmdObject(name="random_tx",
                    message="Select randomly up to m transcript for each gene.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    group="selection",
                    desc=__doc__,
                    updated=__updated__,
                    test=test)
