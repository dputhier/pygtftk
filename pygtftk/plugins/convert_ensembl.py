#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly

__updated__ = "2018-02-06"

__doc__ = """
 Convert the GTF file to ensembl format. It will essentially add a 'transcript' feature and 'gene' feature when required.
 This command can be viewed as a 'groomer' command for those starting with a non ensembl GTF.
"""
__notes__ = """
    -- The gtftk program is designed to handle files in ensembl GTF format. This means that the GTF file provided to
    gtftk must contain transcript and gene feature/lines. They will be used to get access to transcript and gene
    coordinates whenever needed. This solution was chosen to define a reference GTF file format for gtftk (since Ensembl format is probably the most widely used).
    
    -- Almost all commands of gtftk use transcript_id or gene_id as keys to perform operation on genomic coordinates.
    One of the most common issue when working with  gene coordinates is the lack  of non ambiguous gene or transcript names
    For instance, a refSeq sequence ID used as transcript_id can be associated to  several chromosomal locations as a
    sequence may be duplicated. These identifiers are ambiguous and thus should be avoid. Use UCSC or ensembl IDs instead.
    
"""


def make_parser():
    """The parser."""
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

    parser_grp.add_argument('-n', '--no-check-gene-chr',
                            help="By default the command raise an error if several chromosomes are associated with the same gene_id. Disable this behaviour (but you should better think about what it means...).",
                            action="store_true")

    return parser


def convert_ensembl(inputfile=None,
                    outputfile=None,
                    no_check_gene_chr=False):
    """
    Convert the GTF file to ensembl format.
    """

    GTF(inputfile,
        check_ensembl_format=False
        ).convert_to_ensembl(check_gene_chr=not no_check_gene_chr,
                             ).write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    convert_ensembl(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    #count: test additional args
    @test "convert_ensembl_1" {
     result=`gtftk get_example | awk '$3 != "transcript"' | awk '$3 != "gene"' | gtftk convert_ensembl | gtftk select_by_key -k feature -v gene | cut -f4 | perl -npe 's/\\n/,/'`
      [ "$result" = "125,180,50,65,33,22,107,210,3,176," ]
    }

    @test "convert_ensembl_2" {
     result=`gtftk get_example | awk '$3 != "transcript"' | awk '$3 != "gene"' | gtftk convert_ensembl | gtftk select_by_key -k feature -v gene | cut -f5 | perl -npe 's/\\n/,/'`
      [ "$result" = "138,189,61,76,47,35,116,222,14,186," ]
    }

    @test "convert_ensembl_3" {
     result=`gtftk get_example | awk '$3 != "transcript"' | awk '$3 != "gene"' | gtftk convert_ensembl | gtftk select_by_key -k feature -v transcript | cut -f4 | perl -npe 's/\\n/,/'`
      [ "$result" = "125,125,180,50,65,65,33,22,28,107,107,210,3,3,176," ]
    }

    @test "convert_ensembl_4" {
     result=`gtftk get_example | awk '$3 != "transcript"' | awk '$3 != "gene"' | gtftk convert_ensembl | gtftk select_by_key -k feature -v transcript | cut -f5 | perl -npe 's/\\n/,/'`
      [ "$result" = "138,138,189,61,76,76,47,35,35,116,116,222,14,14,186," ]
    }

    @test "convert_ensembl_5" {
     result=`gtftk get_example | awk '$3 != "transcript"' | awk '$3 != "gene"' | gtftk convert_ensembl | wc -l`
      [ "$result" -eq 70 ]
    }
    
    # regenerate delete a transcript line
    @test "convert_ensembl_6" {
     result=`gtftk get_example | grep -v "transcript.*gene_id.*G0001T002"| wc -l `
      [ "$result" -eq 69 ]
    }

    # regenerate delete the transcript line
    @test "convert_ensembl_7" {
     result=`gtftk get_example | grep -v "transcript.*gene_id.*G0001T002"| gtftk convert_ensembl | wc -l `
      [ "$result" -eq 70 ]
    }


    # regenerate delete a gene line
    @test "convert_ensembl_8" {
     result=`gtftk get_example | grep -v "gene.*gene_id.*G0010"| wc -l `
      [ "$result" -eq 69 ]
    }

    # regenerate delete the gene line
    @test "convert_ensembl_9" {
     result=`gtftk get_example | grep -v "gene.*gene_id.*G0010"| gtftk convert_ensembl | wc -l `
      [ "$result" -eq 70 ]
    }

    # Check the md5sum-lite signature of get_example
    @test "convert_ensembl_10" {
     result=`gtftk get_example | md5sum-lite | sed 's/ .*//' `
      [ "$result" = "679aa6be7ee8d8402f4d05e05d2b49d5" ]
    }
       
    # Check that the md5sum-lite signature is the same after regenerating...
    @test "convert_ensembl_11" {
     result=`gtftk get_example | grep -v "gene.*gene_id.*G0010"| gtftk convert_ensembl | md5sum-lite | sed 's/ .*//'`
      [ "$result" = "679aa6be7ee8d8402f4d05e05d2b49d5" ]
    }
           
    # Delete all genes and transcripts, regenerate, check md5sum-lite...
    @test "convert_ensembl_12" {
     result=`gtftk get_example | awk '$3 != "transcript"' | awk '$3 != "gene"' | gtftk convert_ensembl  | md5sum-lite | sed 's/ .*//'`
      [ "$result" = "679aa6be7ee8d8402f4d05e05d2b49d5" ]
    }
        
            
    """

    CMD = CmdObject(name="convert_ensembl",
                    message="Convert the GTF file to ensembl format. Essentially add 'transcript'/'gene' features.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    notes=__notes__,
                    desc=__doc__,
                    group="conversion",
                    test=test)
