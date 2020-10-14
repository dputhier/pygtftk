'''
Merge feature from a BED file by strands (i.e produce a BED6 file with strand info).
'''

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.bedtool_extension import BedTool


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True,
                                     description=__doc__)

    parser_grp = parser.add_argument_group('Arguments')

    # --------------------- Main arguments ----------------------------------- #

    parser_grp.add_argument('-i', '--inputfile',
                            help='The input BED file.',
                            default=sys.stdin,
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('bed')),
                            required=False)

    parser_grp.add_argument('-o', '--outputfile',
                            help="The output BED file.",
                            default=sys.stdout,
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('bed')),
                            required=False)

    return parser


def merge_bed_by_strands(inputfile=None,
                         outputfile=None):
    '''
    Merge feature from a BED file by strands (i.e produce a BED6 file with strand info).
    '''

    # -------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------

    in_bed = BedTool(inputfile).sort()
    out_bed = in_bed.merge_by_strand()
    if outputfile.name == "<stdout>":
        for line in out_bed:
            print(line)
    else:
        out_bed.saveas(outputfile.name)


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    merge_bed_by_strands(**args)


if __name__ == '__main__':
    main()

else:
    test = """
    #bigwig_to_bed
    @test "bigwig_to_bed_0" {
     result=`gtftk get_example -f '*' `
      [ "$result" = "" ]
    }

    
    """
    from pygtftk.cmd_object import CmdObject

    CmdObject(name="merge_bed_by_strands",
              message="Merge feature from a BED file by strands.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              group="miscellaneous",
              test=test)
