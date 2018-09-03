#!/usr/bin/env python
from __future__ import print_function
from builtins import str

import argparse
import random
import sys

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import int_greater_than_null
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
 Select a random list of genes or transcripts. Note that if transcripts
 are requested the 'gene' feature is not returned.
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

    parser_grp.add_argument('-n', '--number',
                            help="The number of transcripts or gene to select.",
                            default=1,
                            metavar="NUMBER",
                            type=int_greater_than_null,
                            required=False)

    parser_grp.add_argument('-t', '--ft-type',
                            help="The type of feature.",
                            default="transcript",
                            choices=["gene", "transcript"],
                            required=False)

    parser_grp.add_argument('-s', '--seed-value',
                            help="Seed value for the random number generator.",
                            default=None,
                            metavar="SEED",
                            type=int_greater_than_null,
                            required=False)

    return parser


def random_list(
        inputfile=None,
        outputfile=None,
        number=None,
        ft_type=None,
        seed_value=None,
        tmp_dir=None,
        logger_file=None,
        verbosity=0):
    """
    Select a random list of genes or transcripts.
    """

    message("loading the GTF.")

    gtf = GTF(inputfile)

    message("Getting ID list.")

    if ft_type == 'gene':
        id_list = set(gtf.extract_data("gene_id", as_list=True))
    else:
        id_list = set(gtf.extract_data("transcript_id", as_list=True))

    if number > len(id_list):
        message("To much feature. Using : " + str(len(id_list)),
                type="WARNING")
        number = len(id_list)

    if seed_value is not None:
        random.seed(seed_value)

    id_list = random.sample(id_list, number)

    message("Printing.")

    my_id = ft_type + "_id"

    gtf.select_by_key(my_id, ",".join(id_list)).write(outputfile)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    random_list(**args)


if __name__ == '__main__':
    main()

else:

    test = '''
       
    #random_list: set seed and return G0010 (4 lines)
    @test "random_list_1" {
     result=`gtftk random_list -i pygtftk/data/simple/simple.gtf -n 1 -t gene -s 123| wc -l`
      [ "$result" -eq 4 ]
    }
    
    #random_list: set seed and return G0008 (4 lines)
    @test "random_list_2" {
     result=`gtftk random_list -i pygtftk/data/simple/simple.gtf -n 1 -t gene -s 111| wc -l`
      [ "$result" -eq 5 ]
    }
    
    #random_list: set seed and return G0007T001 (3 lines)
    @test "random_list_3" {
     result=`gtftk random_list -i pygtftk/data/simple/simple.gtf -n 1  -s 111 -t transcript| wc -l`
      [ "$result" -eq 3 ]
    }
    
    #random_list: set seed and return G0007T001, G0004T002 (10 lines)
    @test "random_list_4" {
     result=`gtftk random_list -i pygtftk/data/simple/simple.gtf -n 2  -s 111 -t transcript | wc -l`
      [ "$result" -eq 10 ]
    }
    
    '''

    CMD = CmdObject(name="random_list",
                    message="Select a random list of genes or transcripts.",
                    parser=make_parser(),
                    fun=random_list,
                    group="selection",
                    updated=__updated__,
                    desc=__doc__,
                    test=test)
