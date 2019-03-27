# -*- coding: utf-8 -*-
"""Select columns from a tabulated file based on their names."""

import argparse
import os
import re
import sys

from pygtftk import arg_formatter
from pygtftk.utils import chomp
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"
__doc__ = """
  Select columns from a tabulated file based on their names.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="The tabulated file. Default to STDIN",
                            default=sys.stdin,
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='txt'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='txt'))

    parser_grp.add_argument('-c', '--columns',
                            help="The list (csv) of column names.",
                            default=None,
                            type=str,
                            required=True)

    parser_grp.add_argument('-n', '--invert-match',
                            help='Not/invert match.',
                            action="store_true")

    parser_grp.add_argument('-u', '--unique',
                            help='Write non redondant lines.',
                            action="store_true")

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="\t",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-H', '--no-header',
                            help="Don't print the header line.",
                            action="store_true",
                            required=False)

    return parser


def col_from_tab(inputfile=None,
                 outputfile=None,
                 columns=None,
                 invert_match=False,
                 no_header=False,
                 unique=False,
                 separator=None):
    """Select columns from a tabulated file based on their names."""

    line_set = dict()

    if re.search(",", columns):
        columns = columns.split(",")
    else:
        columns = [columns]

    for p, line in enumerate(inputfile):

        line = chomp(line)
        line = line.split(separator)

        if p == 0:

            if not invert_match:

                pos_list = list()

                for i in range(len(columns)):

                    pos = line.index(columns[i]) if columns[i] in line else -1

                    if pos > -1:
                        pos_list.append(pos)
                    else:
                        message("Column " + columns[i] + " not found",
                                type="ERROR")

            else:

                pos_list = list(range(len(line)))

                for i in range(len(columns)):

                    pos = line.index(columns[i]) if columns[i] in line else -1

                    if pos > -1:
                        pos_list.remove(pos)
                    else:
                        message("Column " + columns[i] + " not found",
                                type="ERROR")

            if not no_header:
                out = separator.join([line[k] for k in pos_list])
                write_properly(out, outputfile)
        else:
            out = separator.join([line[k] for k in pos_list])
            if unique:
                if out not in line_set:
                    write_properly(out, outputfile)
                    line_set[out] = 1
            else:
                write_properly(out, outputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    col_from_tab(**args)


if __name__ == '__main__':
    main()


else:
    test = """
    #col_from_tab
    @test "col_from_tab_1" {
     result=`gtftk get_example | gtftk tabulate -k all -x |gtftk col_from_tab -c start,end,seqid| wc -l`
      [ "$result" -eq 71 ]
    }

    #col_from_tab
    @test "col_from_tab_2" {
     result=`gtftk get_example | gtftk tabulate -k all -x |gtftk col_from_tab -H -c start,end,seqid| wc -l`
      [ "$result" -eq 70 ]
    }
    
    #col_from_tab
    @test "col_from_tab_3" {
     result=`gtftk get_example | gtftk tabulate -k all -x |gtftk col_from_tab -c start,end,seqid| awk 'BEGIN{FS=OFS="\\t"}{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    """
    from pygtftk.cmd_object import CmdObject

    CmdObject(name="col_from_tab",
              message="Select columns from a tabulated file based on their names.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              desc=__doc__,
              group="miscellaneous",
              test=test)
