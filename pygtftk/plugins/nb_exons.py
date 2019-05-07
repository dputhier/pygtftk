#!/usr/bin/env python
"""
 Returns the transcript name and number of exons with nb_exons as a novel key for each transcript feature.
"""
import argparse
import os
import sys
from collections import defaultdict

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"


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
                            metavar="TXT/GTF",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf', 'txt')))

    parser_grp.add_argument('-f', '--text-format',
                            help="Return a text format.",
                            action="store_true")

    parser_grp.add_argument('-a',
                            '--key-name',
                            type=str,
                            default="nb_exons",
                            help="The name of the key.",
                            required=False)

    return parser


def nb_exons(inputfile=None,
             outputfile=None,
             key_name=None,
             text_format=False):
    """
    Count the number of exons in the gtf file.
    """

    gtf = GTF(inputfile)
    n_exons = defaultdict(int)

    # -------------------------------------------------------------------------
    # Computing number of  exon for each transcript in input GTF file
    #
    # -------------------------------------------------------------------------

    message("Computing number of exons for each transcript in input GTF file.")

    exon = gtf.select_by_key("feature", "exon")
    fields = exon.extract_data("transcript_id")

    for i in fields:
        tx_id = i[0]
        n_exons[tx_id] += 1

    if text_format:
        for tx_id in n_exons:
            outputfile.write(tx_id + "\t" + str(n_exons[tx_id]) + "\ttranscript\n")
    else:

        if len(n_exons):
            gtf = gtf.add_attr_from_dict(feat="transcript",
                                         key="transcript_id",
                                         a_dict=n_exons,
                                         new_key=key_name)
        gtf.write(outputfile,
                  gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    nb_exons(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    # nb_exons: load dataset
    @test "nb_exons_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }

        #nb_exons: G0004T001 has 3 exons
    @test "nb_exons_1" {
     result=`cat simple.gtf| gtftk nb_exons -f | grep G0004T001| cut -f2`
      [ "$result" -eq 3 ]
    }
    
    #nb_exons: G0005T001 has 2 exons
    @test "nb_exons_2" {
     result=`cat simple.gtf| gtftk nb_exons -f | grep G0005T001| cut -f2`
      [ "$result" -eq 2 ]
    }
    
    #nb_exons: G0001T002 has 1 exons
    @test "nb_exons_3" {
     result=`cat simple.gtf| gtftk nb_exons -f | grep G0001T002| cut -f2`
      [ "$result" -eq 1 ]
    }
    
    #nb_exons: G0006T001 has 3 exons
    # test -g
    @test "nb_exons_4" {
     result=`cat simple.gtf| gtftk nb_exons | gtftk select_by_key -k feature -v transcript | grep G0006T001 | gtftk tabulate -H -k nb_exons`
      [ "$result" -eq 3 ]
    }
    
    #nb_exons: G0005T001 has 2 exons
    # test -g
    @test "nb_exons_5" {
     result=`cat simple.gtf| gtftk nb_exons | gtftk select_by_key -k feature -v transcript | grep G0005T001| gtftk tabulate -H -k nb_exons`
      [ "$result" -eq 2 ]
    }
    
    #nb_exons: no line lost
    # test -g
    @test "nb_exons_6" {
     result=`cat simple.gtf| gtftk nb_exons | wc -l`
      [ "$result" -eq 70 ]
    }
    
    #nb_exons: G0004T001 has 3 exons
    @test "nb_exons_7" {
     result=`cat simple.gtf| gtftk nb_exons | gtftk select_by_key -k feature -v transcript | grep G0004T001| gtftk tabulate -H -k nb_exons`
      [ "$result" -eq 3 ]
    }
    
    #nb_exons: test stdin
    @test "nb_exons_8" {
     result=`cat simple.gtf| gtftk  nb_exons -f | wc -l`
      [ "$result" -eq 15 ]
    }
   
   """

    CmdObject(name="nb_exons",
              message="Count the number of exons by transcript.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="information",
              updated=__updated__,
              desc=__doc__,
              test=test)
