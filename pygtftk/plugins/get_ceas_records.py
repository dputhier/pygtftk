#!/usr/bin/env python
"""
 Convert a CEAS sqlite file back into a GTF.
"""
import argparse
import gzip
import os
import shutil
import sqlite3
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import make_tmp_file, message

__updated__ = "2018-01-20"


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the CEAS file.",
                            default=None,
                            type=argparse.FileType('r'),
                            required=True)

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            type=arg_formatter.FormattedFile(mode='w', file_ext='gtf'))

    parser_grp.add_argument('-t', '--target-table',
                            help="The target table.",
                            default="GeneTable",
                            type=str,
                            required=False)

    parser_grp.add_argument('-s', '--show-tables',
                            help="Only list tables and exit.",
                            action="store_true",
                            required=False)

    return parser


def get_ceas_records(
        inputfile=None,
        outputfile=None,
        show_tables=False,
        target_table='GeneTable'):
    """
    Convert a CEAS sqlite file back into a flat file.
    """

    # ----------------------------------------------------------------------
    # load the CEAS file
    # ----------------------------------------------------------------------

    if inputfile.name.endswith('gz'):
        tmp_file = make_tmp_file(prefix='ceas_gunzip', suffix='.txt')
        with gzip.open(inputfile.name, 'rb') as f_in:
            with open(tmp_file.name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        inputfile = open(tmp_file.name)

    conn = sqlite3.connect(inputfile.name)
    cursor = conn.cursor()

    # ----------------------------------------------------------------------
    # A func to get the list of tables
    # ----------------------------------------------------------------------

    def get_tables(cursor):

        out_list = list()
        cursor.execute('SELECT name from sqlite_master where type= "table"')

        for rec in cursor.fetchall():
            out_list += [rec[0]]

        return out_list

    # ----------------------------------------------------------------------
    # Get table list
    # ----------------------------------------------------------------------

    tables = get_tables(cursor)

    # ----------------------------------------------------------------------
    # To show table
    # ----------------------------------------------------------------------

    if show_tables:

        for tab in tables:
            outputfile.write(tab + "\n")
        sys.exit()
    # ----------------------------------------------------------------------
    # loop through records
    # Each line contains:
    #   chrom,name,strand,txStart,txEnd,cdsStart,
    #   cdsEnd,exonCount,exonStarts,exonEnds,name
    # ----------------------------------------------------------------------

    # Check tables exists

    if target_table not in tables:
        message('Table is undefined', type="ERROR")

    for rec in cursor.execute('SELECT * FROM % s' % target_table):
        for rec in cursor.fetchall():
            out_list = []
            for elemnt in rec:
                out_list += [str(elemnt)]
            outputfile.write("\t".join(out_list) + "\n")


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_ceas_records(**args)


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

    CMD = CmdObject(name="get_ceas_records",
                    message="Convert a CEAS sqlite file back into a flat file.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="miscellaneous",
                    test=test)
