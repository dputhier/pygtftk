#!/usr/bin/env python
"""
 Convert a GTF to tabulated format.
"""
import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"

__notes__ = """
 -- Warning: by default tabulate will discard any line for which one of the selected key is not defined. Use -x (-\-accept-undef) to print them.
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
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='txt'))

    parser_grp.add_argument('-s', '--separator',
                            help="The output field separator.",
                            default="\t",
                            metavar="SEPARATOR",
                            type=str)

    parser_grp.add_argument('-k', '--key',
                            help='A comma-separated list of key names.',
                            default="*",
                            metavar="KEY,KEY...",
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

    parser_grp.add_argument('-n', '--no-unset',
                            help="Don't print lines containing '.' (unset values)",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-x', '--accept-undef',
                            help="Print line for which the key is undefined (i.e, '?', does not exists).",
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
             no_unset=False,
             unique=False,
             no_basic=False,
             accept_undef=False,
             select_gene_ids=False,
             select_gene_names=False,
             select_transcript_ids=False,
             select_exon_ids=False,
             separator="\t",
             no_header=False):
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

    no_undef = False
    if not accept_undef:
        no_undef = True
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
        message("Writing header")
        write_properly(separator.join(tab.colnames),
                       outputfile)

    message("Writing")

    try:
        if not unique:
            if no_unset:
                if no_undef:
                    for i in tab:
                        if any([True for x in i.fields if x in [".", "?"]]):
                            continue
                        i.write(outputfile, separator)
                else:
                    for i in tab:
                        if any([True for x in i.fields if x in ["."]]):
                            continue
                        i.write(outputfile, separator)

            else:
                if no_undef:
                    for i in tab:
                        if any([True for x in i.fields if x in ["?"]]):
                            continue
                        i.write(outputfile, separator)
                else:
                    for i in tab:
                        i.write(outputfile, separator)

        else:
            printed = {}
            if no_unset:
                if no_undef:
                    for i in tab:
                        t = tuple(i)
                        if t not in printed:
                            if any([True for x in i.fields if x in [".", "?"]]):
                                continue
                            i.write(outputfile, separator)
                        printed[t] = 1
                else:
                    for i in tab:
                        t = tuple(i)
                        if t not in printed:
                            if any([True for x in i.fields if x in ["."]]):
                                continue
                            i.write(outputfile, separator)
                        printed[t] = 1
            else:
                if no_undef:
                    for i in tab:
                        t = tuple(i)
                        if t not in printed:
                            if any([True for x in i.fields if x in ["?"]]):
                                continue
                            i.write(outputfile, separator)
                        printed[t] = 1
                else:
                    for i in tab:
                        t = tuple(i)
                        if t not in printed:
                            i.write(outputfile, separator)
                        printed[t] = 1

    except (BrokenPipeError, IOError):
        def _void_f(*args, **kwargs):
            pass

        message("Received a boken pipe signal", type="WARNING")
        sys.stdout.write = _void_f
        sys.stdout.flush = _void_f

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

    # tabulate: load dataset
    @test "tabulate_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
       
    
    # tabulate: check column number
    @test "tabulate_1" {
     result=`gtftk tabulate -H -i  simple.gtf -k transcript_id,gene_id,start| awk -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    # tabulate: check column number
    @test "tabulate_2" {
     result=`gtftk tabulate -i  simple.gtf -k transcript_id,gene_id,gene_id,start| awk -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 4 ]
    }
    
    
    # tabulate: check separator
    @test "tabulate_3" {
     result=`gtftk tabulate -H -i  simple.gtf -k transcript_id,gene_id,gene_id,start -s ":" | awk -F ":" '{print NF}'| sort | uniq`
      [ "$result" -eq 4 ]
    }

    # tabulate: check -u
    @test "tabulate_4" {
     result=`gtftk tabulate  -i  simple.gtf -k transcript_id,gene_id,start -uHx| wc -l`
      [ "$result" -eq 44 ]
    }
    
    # tabulate: check -n
    @test "tabulate_5" {
     result=`gtftk tabulate  -i  simple.gtf -k transcript_id,gene_id,start -unH| wc -l`
      [ "$result" -eq 34 ]
    }

    # tabulate: check -un
    @test "tabulate_6" {
     result=`gtftk tabulate  -i  simple.gtf -k transcript_id,gene_id,start -un| wc -l`
      [ "$result" -eq 35 ]
    }

    # tabulate: check all as key
    @test "tabulate_7" {
     result=`gtftk tabulate -k all -x -i  simple.gtf| wc -l `
      [ "$result" -eq 71 ]
    }

    # tabulate: check all as key
    @test "tabulate_8" {
     result=`gtftk tabulate -k all -x -i simple.gtf| awk  -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 12 ]
    }
    
    # tabulate: check "*" as key and -b
    @test "tabulate_9" {
     result=`gtftk tabulate -k "*" -b  -x -i  simple.gtf | awk  -F "\\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 4 ]
    }
    
    # tabulate: check "*" as key and -b
    @test "tabulate_10" {
     result=`gtftk tabulate -i  simple.gtf -k all -x | awk  -F "\t" '{print NF}'| sort | uniq`
      [ "$result" -eq 12 ]
    }
    
    
    
    '''

    cmd = CmdObject(name="tabulate",
                    message="Convert a GTF to tabulated format.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="conversion",
                    updated=__updated__,
                    notes=__notes__,
                    test=test)
