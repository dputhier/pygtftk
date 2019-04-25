#!/usr/bin/env python
"""
 Select a random list of genes or transcripts. Note that if transcripts
 are requested the 'gene' feature is not returned.
"""

import argparse
import os
import random
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"


def make_parser():
    """The parser."""
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
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-n', '--number',
                            help="The number of transcripts or gene to select.",
                            default=1,
                            metavar="NUMBER",
                            type=arg_formatter.ranged_num(lowest=1,
                                                          highest=None,
                                                          linc=True,
                                                          val_type='int'),
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
                            type=arg_formatter.ranged_num(lowest=1,
                                                          highest=None,
                                                          linc=True,
                                                          val_type='int'),
                            required=False)

    return parser


def random_list(
        inputfile=None,
        outputfile=None,
        number=None,
        ft_type=None,
        seed_value=None):
    """
    Select a random list of genes or transcripts.
    """

    message("loading the GTF.")

    gtf = GTF(inputfile)

    message("Getting ID list.")

    if ft_type == 'gene':
        id_list = gtf.extract_data("gene_id", as_list=True, nr=True, hide_undef=True, no_na=True)
    else:
        id_list = gtf.extract_data("transcript_id", as_list=True, nr=True, hide_undef=True, no_na=True)

    if number > len(id_list):
        message("To much feature. Using : " + str(len(id_list)),
                type="WARNING")
        number = len(id_list)

    if seed_value is not None:
        random.seed(seed_value, version=1)

    id_list = random.sample(id_list, number)

    message("Printing.")

    my_id = ft_type + "_id"

    gtf.select_by_key(my_id, ",".join(id_list)).write(outputfile,
                                                      gc_off=True)

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
    #random_list: load dataset
    @test "random_list_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
    
    #random_list: set seed and return G0010 (4 lines)
    @test "random_list_1" {
     result=`gtftk random_list -i simple.gtf -n 1 -t gene -s 123| wc -l`
      [ "$result" -eq 7 ]
    }
    
    #random_list: set seed and return G0007T001 (3 lines)
    @test "random_list_3" {
     result=`gtftk random_list -i simple.gtf -n 1  -s 111 -t transcript| wc -l`
      [ "$result" -eq 3 ]
    }
    
    '''

    CMD = CmdObject(name="random_list",
                    message="Select a random list of genes or transcripts.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    group="selection",
                    updated=__updated__,
                    desc=__doc__,
                    test=test)
