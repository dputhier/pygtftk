#!/usr/bin/env python

"""
This is the doc about the command that will appear when gtftk my_command -h
is called...
"""
import os
import sys
import argparse
from pygtftk.cmd_object import CmdObject
from pygtftk import arg_formatter

from pygtftk.utils import message
from pygtftk.utils import make_tmp_file

#-------------------------------------------------------------------------
# Command information
#-------------------------------------------------------------------------

__notes__ = """
-- A note that will appear when 'gtftk my_command -h' will be called.
-- Another note. If you want to refer to long form arguments use '\'. e.g -\-distance.
"""


#-------------------------------------------------------------------------
# First define the function/command arguments.
# Note that the syntax is the same that would be used for a regular program
# implementing an argument parser.
# Make use as possible of argparse.FileType and more complexes types as
# found in gtftk.arg_formatter.
#-------------------------------------------------------------------------

def make_parser():
   parser = argparse.ArgumentParser(add_help=True)

   parser_grp = parser.add_argument_group('Arguments')

   parser_grp.add_argument('-i', '--inputfile',
                           help="A csv file as produced by OLOGRAM.",
                           default=sys.stdin,
                           metavar="CSV",
                           type=arg_formatter.FormattedFile(mode='r', file_ext=('txt')),
                           required=True)


   parser_grp.add_argument('-o', '--outputfile',
                           help="Output file.",
                           default=sys.stdout,
                           metavar="GTF",
                           type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))
   return parser

#-------------------------------------------------------------------------
# Now we declare a main function, as would be done
# for a regular program
#-------------------------------------------------------------------------


# NB: The verbosity, tmp_dir=None and logger_file are mandatory arguments

def html_diagram(inputfile=None,
              outputfile=None):
    """This function will produced dynamic diagrams from OLOGRAM output csv file."""

    message("Reading csv file")



#-------------------------------------------------------------------------
# Now we check if the python interpreter is running this module
# as the main program or whether it is called by the plugin manager.
#-------------------------------------------------------------------------

def main():
    """The main function."""
    args = make_parser().parse_args()
    args = dict(args.__dict__)
    html_diagram(**args)

if __name__ == '__main__':
    message("I am in main")
    main()
else:
    test = """

       # coverage: load dataset
       @test "coverage_0" {
        result=`gtftk get_example -f '*' -d simple`
         [ "$result" = "" ]
       }  
    """

    # Just declare a new command object
    # That will call the command manager.
    # With the user-passed arguments.
    # Available groups are: editing, information, selection, conversion,
    # coordinates, annotation, sequences, coverage,
    # and miscellaneous.

    cmd = CmdObject(name="html_diagram",
                    message="This command produces dynamic diagrams from OLOGRAM output csv files.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    group="ologram",
                    desc=__doc__,
                    notes=__notes__,
                    test=test)
