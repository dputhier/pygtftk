#!/usr/bin/env python
"""

"""
from __future__ import print_function

import argparse
import sys

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"
__doc__ = """
 Convert a GTF to tabulated format.
"""
__notes__ = """
 -- To refer to default keys use: seqid,source,feature,start,end,frame,gene_id...
 -- Note that 'all' or '*' are special keys that can be used to convert the whole GTF into a tabulated file. Thanks @fafa13.
"""


def make_parser():
    """The program CLI."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_mut = parser_grp.add_mutually_exclusive_group(required=False)

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
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
                                                                     '\.[Tt][Ss][Vv]')))
    parser_grp.add_argument('-s', '--separator',
                            help="The output field separator.",
                            default="\t",
                            metavar="SEPARATOR",
                            type=str)

    parser_grp.add_argument('-k', '--key',
                            help='A comma separated list of key names.',
                            default="*",
                            metavar="KEY,KEY,...",
                            type=str,
                            required=False)

    parser_grp.add_argument('-u', '--unique',
                            help='Print a non redondant list of lines.',
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-H', '--no-header',
                            help="Don't print the header line.",
                            default=False,
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-n', '--no-undef',
                            help="Don't print lines containing '.' (undefined values)",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-b', '--no-basic',
                            help="In case key is set to 'all' or '*', don't write basic attributes.",
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-t', '--select-transcript-ids',
                            help='A shortcuts for "-k transcript_id".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-g', '--select-gene_ids',
                            help='A shortcuts for "-k gene_id".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-a', '--select-gene-names',
                            help='A shortcuts for "-k gene_name".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-e', '--select-exon-ids',
                            help='A shortcuts for "-k exon_ids".',
                            action="store_true",
                            required=False)

    return parser


def tabulate(inputfile=None,
             outputfile=None,
             key=None,
             tmp_dir=None,
             no_undef=False,
             unique=False,
             no_basic=False,
             select_gene_ids=False,
             select_gene_names=False,
             select_transcript_ids=False,
             select_exon_ids=False,
             separator="\t",
             logger_file=None,
             no_header=False,
             verbosity=0):
    """Convert a GTF to tabulated format.
    """

    # ----------------------------------------------------------------------
    # Check mode
    # ----------------------------------------------------------------------

    if select_transcript_ids:
        key = "transcript_id"

    elif select_gene_ids:
        key = "gene_id"

    elif select_gene_names:
        key = "gene_id"

    elif select_exon_ids:
        key = "exon_id"

    # ----------------------------------------------------------------------
    # REad GTF and process
    # ----------------------------------------------------------------------

    gtf = GTF(inputfile, check_ensembl_format=False)

    if key in ["all", "*"]:
        if no_basic:
            attr_list = gtf.get_attr_list(add_basic=False)
        else:
            attr_list = gtf.get_attr_list(add_basic=True)
        tab = gtf.extract_data(attr_list)
    else:
        tab = gtf.extract_data(key)

    if not no_header:
        write_properly(separator.join(tab.colnames),
                       outputfile)

    if not unique:
        if no_undef:
            for i in tab:
                if any([True for x in i.fields if x == "."]):
                    continue
                i.write(outputfile, separator)
        else:
            for i in tab:
                i.write(outputfile, separator)

    else:
        printed = {}
        if no_undef:
            for i in tab:
                t = tuple(i)
                if t not in printed:
                    if any([True for x in i.fields if x == "."]):
                        continue
                    i.write(outputfile, separator)
                printed[t] = 1
        else:
            for i in tab:
                t = tuple(i)
                if t not in printed:
                    i.write(outputfile, separator)
                printed[t] = 1

    close_properly(outputfile, inputfile)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    tabulate(**args)


if __name__ == '__main__':
    main()

else:

    test = '''
    # tabulate: check column number
    @test "tabulate_1" {
     result=`gtftk tabulate -H -i  pygtftk/data/simple/simple.gtf -k transcript_id,gene_id,start| awk -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    # tabulate: check column number
    @test "tabulate_2" {
     result=`gtftk tabulate -i  pygtftk/data/simple/simple.gtf -k transcript_id,gene_id,gene_id,start| awk -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 4 ]
    }
    
    
    # tabulate: check separator
    @test "tabulate_3" {
     result=`gtftk tabulate -H -i  pygtftk/data/simple/simple.gtf -k transcript_id,gene_id,gene_id,start -s ":" | awk -F ":" '{print NF}'| sort | uniq`
      [ "$result" -eq 4 ]
    }

    # tabulate: check -u
    @test "tabulate_4" {
     result=`gtftk tabulate  -i  pygtftk/data/simple/simple.gtf -k transcript_id,gene_id,start -uH| wc -l`
      [ "$result" -eq 44 ]
    }
    
    # tabulate: check -n
    @test "tabulate_5" {
     result=`gtftk tabulate  -i  pygtftk/data/simple/simple.gtf -k transcript_id,gene_id,start -unH| wc -l`
      [ "$result" -eq 34 ]
    }

    # tabulate: check -un
    @test "tabulate_6" {
     result=`gtftk tabulate  -i  pygtftk/data/simple/simple.gtf -k transcript_id,gene_id,start -un| wc -l`
      [ "$result" -eq 35 ]
    }

    # tabulate: check all as key
    @test "tabulate_7" {
     result=`gtftk get_example | gtftk tabulate -k all | awk  -F "\\t" '{print NF}'| wc -l `
      [ "$result" -eq 71 ]
    }

    # tabulate: check all as key
    @test "tabulate_8" {
     result=`gtftk get_example | gtftk tabulate -k all | awk  -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 12 ]
    }
    
    # tabulate: check "*" as key and -b
    @test "tabulate_9" {
     result=`gtftk get_example | gtftk tabulate -k "*" -b  | awk  -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 4 ]
    }
    
    '''

    cmd = CmdObject(name="tabulate",
                    message="Convert a GTF to tabulated format.",
                    parser=make_parser(),
                    fun=tabulate,
                    desc=__doc__,
                    group="conversion",
                    updated=__updated__,
                    notes=__notes__,
                    test=test)
