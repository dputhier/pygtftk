#!/usr/bin/env python
"""
 Convert a bed file to a gtf. This will make the poor bed feel as if it was a
 big/fat gtf (but with lots of empty fields...sniff). May be helpful sometimes...
"""

import argparse
import os
import sys

from pybedtools import BedTool

import pygtftk.utils
from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-02-11"


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the poor BED file to would like to "
                                 "behave as if it was a GTF.",
                            default=sys.stdin,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bed'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-t', '--ft-type',
                            help="The type of features you are trying to "
                                 "mimic...",
                            default='transcript',
                            type=str,
                            required=False)

    parser_grp.add_argument('-s', '--source',
                            help="The source of annotation.",
                            default='Unknown',
                            type=str,
                            required=False)
    return parser


def bed_to_gtf(
        inputfile=None,
        outputfile=None,
        ft_type="transcript",
        source="Unknown"):
    """
 Convert a bed file to a gtf. This will make the poor bed feel as if it was a
 nice gtf (but with lots of empty fields...). May be helpful sometimes...
    """

    message("Converting the bed file into GTF file.")

    if inputfile.name == '<stdin>':
        tmp_file = make_tmp_file(prefix="input_bed", suffix=".bed")
        for i in inputfile:
            write_properly(chomp(str(i)), tmp_file)

        tmp_file.close()
        inputfile.close()

        bed_obj = BedTool(tmp_file.name)
    else:
        bed_obj = BedTool(inputfile.name)

    n = 1
    for i in bed_obj:

        if i.strand == "":
            i.strand = "."
        if i.name == "":
            i.name = str("feature_" + str(n))
        if i.score == "":
            i.score = "0"

        if ft_type == "exon":
            key_value = "gene_id \"" + i.name + "\"; " + \
                        "transcript_id \"" + i.name + "\"; " + \
                        "exon_id \"" + i.name + "\";"
        elif ft_type == "gene":
            key_value = "gene_id \"" + i.name + "\";"
        else:
            key_value = "gene_id \"" + i.name + "\"; " + \
                        "transcript_id \"" + i.name + "\";"

        if pygtftk.utils.ADD_CHR == 1:
            chrom_out = "chr" + i.chrom
        else:
            chrom_out = i.chrom

        list_out = [chrom_out,
                    source,
                    ft_type,
                    str(i.start + 1),
                    str(i.end),
                    str(i.score),
                    i.strand,
                    ".",
                    key_value]

        write_properly("\t".join(list_out), outputfile)

        n += 1

    close_properly(outputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    bed_to_gtf(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #bed_to_gtf: load dataset
    @test "bed_to_gtf_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }

    #bed_to_gtf: test stdin and ine number
    @test "bed_to_gtf_1" {
     result=`gtftk select_by_key -i simple.gtf -k feature -v transcript| gtftk convert -f bed | gtftk bed_to_gtf | wc -l`
      [ "$result" -eq 15 ]
    }
    
    #bed_to_gtf: test column number
    @test "bed_to_gtf_2" {
     result=`gtftk select_by_key -i simple.gtf -k feature -v transcript| gtftk convert -f bed | gtftk bed_to_gtf | awk 'BEGIN{FS="\\t"}{print NF}' | sort | uniq`
      [ "$result" -eq 9 ]
    }

    """

    CMD = CmdObject(name="bed_to_gtf",
                    message="Convert a bed file to a gtf but with lots of empty fields...",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="conversion",
                    test=test)
