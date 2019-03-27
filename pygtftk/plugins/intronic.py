#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"
__doc__ = """
 Returns a bed file containing the intronic regions. If by_transcript is false
 (default), returns merged genic regions with no exonic overlap ("strict" mode).
 Otherwise, the intronic regions corresponding to each transcript are returned
 (may contain exonic overlap and redundancy).
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
                            help="Output file (BED).",
                            default=sys.stdout,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='bed'))

    parser_grp.add_argument('-b', '--by-transcript',
                            help="The intronic regions are returned for each"
                                 " transcript.",
                            action="store_true")

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name (if -b is used).",
                            default="gene_id,transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (if -b is used).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-w', '--intron-nb-in-name',
                            help="By default intron number is written in 'score' column. Force it to be written in 'name' column.",
                            action="store_true")

    parser_grp.add_argument('-F', '--no-feature-name',
                            help="Don't add the feature name ('intron') in the name column.",
                            action="store_true")

    return parser


def intronic(
        inputfile=None,
        outputfile=None,
        names='transcript_id',
        separator="_",
        intron_nb_in_name=False,
        no_feature_name=False,
        by_transcript=False):
    """
 Extract intronic regions.
    """

    message("Searching for intronic regions.")

    # Need to load if the gtf comes from
    # <stdin>
    gtf = GTF(inputfile, check_ensembl_format=False)

    if not by_transcript:
        introns_bo = gtf.get_introns()

        for i in introns_bo:
            write_properly(chomp(str(i)), outputfile)
    else:

        introns_bo = gtf.get_introns(by_transcript=True,
                                     name=names.split(","),
                                     sep=separator,
                                     intron_nb_in_name=intron_nb_in_name,
                                     feat_name=not no_feature_name)
        for i in introns_bo:
            write_properly(chomp(str(i)), outputfile)

    close_properly(outputfile, inputfile)


def main():
    """The main program."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    intronic(**args)


if __name__ == '__main__':
    main()

else:

    test = '''
    #intronic: load dataset
    @test "intronic_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
    
    #intronic: check coordinates
    @test "intronic_1" {
     result=`gtftk intronic -i simple.gtf | cut -f2,3| perl -npe 's/\\t/|/g; s/\\n/,/g'`
      [ "$result" = "25|27,30|32,35|41,54|56,68|70,71|73,214|219," ]
    }

    #intronic: check number of lines
    @test "intronic_2" {
     result=`gtftk intronic -i simple.gtf | wc -l`
      [ "$result" -eq 7 ]
    }

    #intronic: check number of lines
    @test "intronic_3" {
     result=`gtftk intronic -i simple.gtf -b | cut -f2 | perl -npe  's/\\n/,/'`
      [ "$result" = "54,68,71,68,71,35,25,30,30,214," ]
    }
    
    #intronic: check number of lines
    @test "intronic_4" {
     result=`gtftk intronic -i simple.gtf -b | cut -f3 | perl -npe  's/\\n/,/'`
      [ "$result" = "56,70,73,70,73,41,27,32,32,219," ]
    }
    
    '''

    CmdObject(name="intronic",
              message="Extract intronic regions.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="coordinates",
              desc=__doc__,
              updated=__updated__,
              test=test)
