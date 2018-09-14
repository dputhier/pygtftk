#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys
from collections import defaultdict

from builtins import str
from builtins import zip

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-01-20"
__doc__ = """
 Count the number of values/levels for a set of keys.
"""

__notes__ = """
 -- Use -\-uniq to get the count of non-redondant values.
 -- Not available values ('.') are not taken into account.
 -- Use "*" as the key to get the counts fro all keys.
"""


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN.",
                            default=sys.stdin,
                            metavar="GTF",
                            required=False,
                            type=FileWithExtension('r',
                                                   valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Cc][Oo][Vv]')))

    parser_grp.add_argument('-k', '--keys',
                            default="*",
                            help="The set of keys of interest.",
                            type=str,
                            required=False)

    parser_grp.add_argument('-t', '--additional-text',
                            help="A facultative text to be printed in the third "
                                 "column (e.g species name).",
                            default=None,
                            metavar="TEXT",
                            type=str,
                            required=False)

    parser_grp.add_argument('-u', '--uniq',
                            help="Ask for the count of non redondant values.",
                            action="store_true",
                            required=False)
    return parser


def count_key_values(
        inputfile=None,
        outputfile=None,
        keys="gene_id,transcript_id",
        uniq=True,
        tmp_dir=None,
        additional_text=None,
        logger_file=None,
        verbosity=0):
    """
 Count the number values for a set of keys.
    """

    gtf = GTF(inputfile, check_ensembl_format=False)

    if uniq:
        val_list = defaultdict(set)
    else:
        val_list = defaultdict(list)

    if keys == "*":
        key_list = gtf.get_attr_list()
        keys = ",".join(key_list)
    else:
        key_list = keys.split(",")

    for i in gtf.extract_data(keys, as_list_of_list=True):

        for k, v in zip(key_list, i):
            if v == '.':
                continue
            if uniq:
                val_list[k].add(v)
            else:
                val_list[k] += [v]

    for i in key_list:
        if additional_text is None:
            outputfile.write(i + "\t" + str(len(val_list[i])) + "\n")
        else:
            outputfile.write(i + "\t" + str(len(val_list[i])) + "\t" +
                             additional_text + "\n")
    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    count_key_values(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #count
    @test "count_1" {
     result=`gtftk get_example  | gtftk count| grep exon| cut -f2 `
      [ "$result" -eq 25 ]
    }

    """

    CMD = CmdObject(name="count_key_values",
                    message="Count the number values for a set of keys.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="information",
                    test=test)
