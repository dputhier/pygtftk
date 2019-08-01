#!/usr/bin/env python
"""
 Delete one or several attributes from the gtf file.
"""

import argparse
import os
import re
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
 -- You may also use 'complex' regexp such as : "(_id)|(_b.*pe)"
 -- Example: gtftk get_example -d mini_real | gtftk del_attr -k "(^.*_id$|^.*_biotype$)" -r -v
 -- TODO: currently a segfault is thrown when no keys are left after deletion (libgtftk issue #98).
"""


def make_parser():
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

    parser_grp.add_argument('-k', '--key',
                            help='comma-separated list of attribute names or a '
                                 'regular expression (see -r).',
                            default=None,
                            metavar="KEY",
                            type=str,
                            required=True)

    parser_grp.add_argument('-r', '--reg-exp',
                            help='The key name is a regular expression.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-v', '--invert-match',
                            help='Delected keys are those not matching any of'
                                 ' the specified key.',
                            action="store_true",
                            required=False)

    return parser


def del_attr(
        inputfile=None,
        outputfile=None,
        key="transcript_id",
        reg_exp=False,
        invert_match=False):
    """
    Delete extended attributes in the target gtf file. attr_list can be a
    comma-separated list of attributes.
    """

    # ----------------------------------------------------------------------
    # Read the GTF and get the list of attributes
    # ----------------------------------------------------------------------

    gtf = GTF(inputfile, check_ensembl_format=False)

    attr_list = gtf.attr_extended

    # ----------------------------------------------------------------------
    # If regExp, select the corresponding keys
    # ----------------------------------------------------------------------

    if reg_exp:

        key_list = []

        try:
            rgxp = re.compile(key)
        except:
            message("Check the regular expression please.", type="ERROR")

        for attr in attr_list:
            if rgxp.search(attr):
                key_list += [attr]
    else:
        key_list = key.split(",")

    # ----------------------------------------------------------------------
    # If invert-match select all but the selected
    # ----------------------------------------------------------------------

    key_to_del = []
    if invert_match:
        for attr in attr_list:
            if attr not in key_list:
                key_to_del += [attr]
    else:
        key_to_del = key_list

    # ----------------------------------------------------------------------
    # Delete the keys
    # ----------------------------------------------------------------------

    gtf = gtf.del_attr(feat="*",
                       keys=",".join(key_list), force=True).write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    del_attr(**args)

else:

    test = """

    # del_attr: load dataset
    @test "del_attr_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }


    #del_attr:
    # If you delete almost all extended attributes there are only exon_id left
    @test "del_attr_1" {
     result=`gtftk del_attr -i simple.gtf  -k ccds_id,transcript_id| grep ccds_id | wc -l`
      [ "$result" -eq 0 ]
    }

    #del_attr: check -v
    @test "del_attr_2" {
     result=`gtftk del_attr -i simple.gtf  -k ccds_id,transcript_id -v| gtftk tabulate | cut -f 9,10| head -1| perl -npe 's/\t/,/g'`
      [ "$result" = "gene_id,exon_id" ]
    }

    """

    CmdObject(name="del_attr",
              message="Delete attributes in the target gtf file.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              notes=__notes__,
              desc=__doc__,
              group="editing",
              test=test)
