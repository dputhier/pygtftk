#!/usr/bin/env python

import argparse
import os

from pygtftk.cmd_manager import CmdManager
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import message

__updated__ = "2017-09-27"
__doc__ = """
 Search in all command description files those related to a user-defined keyword.
"""


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-k', '--keyword',
                            help="The keyword",
                            default=None,
                            type=str,
                            required=True)

    parser_grp.add_argument('-n', '--notes',
                            help="Look also for the keywords in notes associated to each command.",
                            action='store_true')

    return parser


def apropos(keyword="",
            notes=False):
    """
    Search in all command description files those related to a user-defined keyword.
    """

    out_list = set()

    for i in CmdManager.cmd_obj_list:
        if keyword in CmdManager.cmd_obj_list[i].desc:
            out_list.add(i)
        if notes:
            try:
                if keyword in CmdManager.cmd_obj_list[i].notes:
                    out_list.add(i)
            except:
                pass

    if out_list:
        message(">> Keyword '" + keyword + "' was found in the following command:",
                force=True)
        for i in out_list:
            print("\t- " + i + ".")
    else:
        message(">> Keyword '" + keyword + "' was not found.",
                force=True)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    apropos(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #count
    @test "apropos_1" {
     result=`gtftk apropos -k antisens | wc -l `
      [ "$result" -eq 2 ]
    }
    

    """

    CMD = CmdObject(name="apropos",
                    message="Search in all command description files those related to a user-defined keyword.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="information",
                    test=test)
