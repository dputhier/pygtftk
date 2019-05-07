#!/usr/bin/env python
"""
 Select lines from a GTF file based on attributes and
 associated values.
"""

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-31"

__notes__ = '''
-- select_by_key only returns lines for which the key is defined (i.e. exists) even with -\-invert-match.
'''


def make_parser():
    """The program parser."""

    parser = argparse.ArgumentParser(add_help=True)
    parser_grp = parser.add_argument_group('Arguments')
    parser_mut = parser.add_mutually_exclusive_group(required=False)

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
                            help='The key name.',
                            default=None,
                            metavar="KEY",
                            type=str,
                            required=False)

    parser_mut.add_argument('-v', '--value',
                            help='A comma-separated list of values.',
                            default=None,
                            metavar="VALUE",
                            type=str)

    parser_mut.add_argument('-f', '--file-with-values',
                            help='A file containing values as a single column.',
                            default=None,
                            metavar="FILE",
                            type=argparse.FileType("r"))

    parser_grp.add_argument('-c', '--col',
                            help='The column number (one-based) that contains the values in the file. File is tab-delimited.',
                            default=1,
                            metavar="COL",
                            type=arg_formatter.ranged_num(lowest=1,
                                                          highest=None,
                                                          linc=True,
                                                          val_type='int'),
                            required=False)

    parser_grp.add_argument('-n', '--invert-match',
                            help='Not/invert match. Selected lines whose requested ke'
                                 'y is not associated with the requested value.',
                            action="store_true")

    parser_grp.add_argument('-b', '--bed-format',
                            help='Ask for bed format output.',
                            action="store_true")

    parser_grp.add_argument('-m', '--names',
                            help="If Bed output. The key(s) that should be used as name.",
                            default="gene_id,transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="If Bed output. The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-l', '--log',
                            help='Print some statistics about selected features. To be used in conjunction with -V 1/2.',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-t', '--select-transcripts',
                            help='A shortcuts for "-k feature -v transcript".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-g', '--select-genes',
                            help='A shortcuts for "-k feature -v gene".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-e', '--select-exons',
                            help='A shortcuts for "-k feature -v exon".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-d', '--select-cds',
                            help='A shortcuts for "-k feature -v CDS".',
                            action="store_true",
                            required=False)

    parser_mut.add_argument('-a', '--select-start-codon',
                            help='A shortcuts for "-k feature -v start_codon".',
                            action="store_true",
                            required=False)

    return parser


def select_by_key(inputfile=None,
                  outputfile=None,
                  key=None,
                  value=None,
                  invert_match=False,
                  file_with_values=None,
                  col=0,
                  select_transcripts=False,
                  select_genes=False,
                  select_exons=False,
                  select_cds=False,
                  select_start_codon=False,
                  bed_format=False,
                  log=False,
                  separator="|",
                  names="transcript_id"):
    """Select lines from a GTF file based on attributes and
    associated values.
    """

    # ----------------------------------------------------------------------
    # Check mode
    # ----------------------------------------------------------------------

    if select_transcripts:
        key = "feature"
        value = "transcript"

    elif select_cds:
        key = "feature"
        value = "CDS"

    elif select_start_codon:
        key = "feature"
        value = "start_codon"

    elif select_genes:
        key = "feature"
        value = "gene"

    elif select_exons:
        key = "feature"
        value = "exon"

    elif file_with_values is None:
        if key is None or value is None:
            message("Key and value are mandatory. Alternatively use -e/t/g/f or -f with -k.",
                    type="ERROR")

    elif file_with_values is not None:
        if key is None:
            message("Please set -k.", type="ERROR")
        if value is not None:
            message("The -f and -v arguments are mutually exclusive.", type="ERROR")

    # ----------------------------------------------------------------------
    # Load file with value
    # ----------------------------------------------------------------------

    gtf = GTF(inputfile, check_ensembl_format=False)
    all_values = gtf.extract_data(key, as_list=True, no_na=True, nr=True)

    if log:
        feat_before = len(gtf)

    if not file_with_values:
        value_list = value.split(",")
        gtf = gtf.select_by_key(key, value, invert_match)
    else:
        value_list = []

        for line in file_with_values:
            cols = line.split("\t")
            value_list += [cols[col - 1]]
        file_with_values.close()
        file_with_values = open(file_with_values.name)

        gtf = gtf.select_by_key(key=key,
                                invert_match=invert_match,
                                file_with_values=file_with_values,
                                col=col)

    if log:

        not_found = list(set(value_list) - set(all_values))
        feat_after = len(gtf)
        pct = feat_after / feat_before * 100

        message("Number of features before selection: %d" % feat_before)
        message("Fraction of feature selected: %.2f%%" % pct)

        if len(not_found):
            nfj = ",".join(not_found)
            max_letter = min(len(nfj), 50)
            if len(nfj) > 50:
                etc = "..."
            else:
                etc = ""
            message("Values not found: [" +
                    ",".join(not_found)[:max_letter] + etc + "].")
        else:
            message("Values not found: [].")

    # ----------------------------------------------------------------------
    # Write GTF file
    # ----------------------------------------------------------------------

    if not bed_format:

        gtf.write(outputfile,
                  gc_off=True)

    else:
        nb_tokens = len(names.split(","))
        keys = "seqid,start,end," + names + ",score,strand"
        nb_fields = len(keys.split(","))

        for i in gtf.extract_data_iter_list(keys, zero_based=True):
            outputfile.write("\t".join([i[0],
                                        i[1],
                                        i[2],
                                        separator.join(i[3:(3 + nb_tokens)]),
                                        i[nb_fields - 2],
                                        i[nb_fields - 1],
                                        ]) + "\n")

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_key(**args)


if __name__ == '__main__':
    main()


else:

    test = '''

    # select_by_key: load dataset
    @test "select_by_key_0" {
     result=`gtftk get_example -f '*' -d simple; gtftk get_example -f '*' -d simple_02`
      [ "$result" = "" ]
    }
  
    # Select_by_key: gene_id selection. Nb Lines.
    @test "select_by_key_1" {
      result=`gtftk select_by_key  -k gene_id -v G0003 -i simple.gtf | wc -l`
      [ "$result" -eq 5 ]
    }
    
    # Select_by_key: gene_id selection. Nb columns.
    @test "select_by_key_2" {
      result=`gtftk select_by_key  -k gene_id -v G0003 -i simple.gtf | awk 'BEGIN{FS="\\t"}{print NF}' | sort | uniq`
      [ "$result" -eq 9 ]
    }
     
    # Select_by_key: Feature selection. Nb Lines.
    @test "select_by_key_3" {
     result=`gtftk select_by_key  -k feature -v gene -i simple.gtf | wc -l`
      [ "$result" -eq 10 ]
    }
    
    # Select_by_key: Feature selection. Nb columns.
    @test "select_by_key_4" {
     result=`gtftk select_by_key  -k feature -v gene -i simple.gtf | awk 'BEGIN{FS="\\t"}{print NF}' | sort | uniq`
      [ "$result" -eq 9 ]
    }
    
    # Select_by_key: Feature selection. --output
    @test "select_by_key_5" {
     result=`gtftk select_by_key  -k gene_id -v G0003 -i simple.gtf  -o /tmp/test_select_by_key.gtf; cat /tmp/test_select_by_key.gtf | wc -l`
      [ "$result" -eq 5 ]
    }
    
    # Select_by_key: Feature selection. --file-with-values -col
    @test "select_by_key_6" {
     result=`rm -f /tmp/test_select_by_key.gtf; gtftk select_by_key  -k gene_id  -i simple.gtf  -f simple.geneList -c 3| wc -l`
      [ "$result" -eq 18 ]
    }
    
    
    # Select_by_key: not operator
    @test "select_by_key_7" {
     result=`gtftk select_by_key -i simple.gtf -k feature -v gene -n| wc -l`
      [ "$result" -eq 60 ]
    }
    
    # Select_by_key: not operator
    @test "select_by_key_8" {
     result=`gtftk select_by_key -i simple.gtf -k feature -v gene -n| awk '$3=="gene"' | wc -l`
      [ "$result" -eq 0 ]
    }
    
    # Select_by_key: not operator
    @test "select_by_key_9" {
     result=`gtftk select_by_key -i simple.gtf -k feature -v transcript -n | wc -l`
      [ "$result" -eq 55 ]
    }
    
    # Select_by_key: not operator
    @test "select_by_key_10" {
     result=`gtftk select_by_key -i simple.gtf -k feature -v transcript -n | awk '$3=="transcript"' | wc -l`
      [ "$result" -eq 0 ]
    }
    
    # Select_by_key: bed conversion
    @test "select_by_key_11" {
     result=`gtftk select_by_key -k feature -v exon -i simple_02.gtf -m transcript_id,gene_id,exon_id  -s "||" -b | wc -l`
      [ "$result" -eq 25 ]
    }
    
    # Select_by_key: not operator with CSV lists.... Expecting 45...
    @test "select_by_key_12" {
     result=`gtftk get_example -f join > simple_join.txt ; gtftk get_example| gtftk select_by_key -f simple_join.txt -c 1 -k gene_id -n | wc -l`
      [ "$result" -eq 45 ]
    }
    
    # Select_by_key: bed conversion
    @test "select_by_key_13" {
     result=`gtftk select_by_key -k feature -v transcript -i simple_02.gtf -m transcript_id,gene_id,exon_id -b  -s "||"| wc -l`
      [ "$result" -eq 15 ]
    }
    
    # Select_by_key: bed conversion
    @test "select_by_key_14" {
     result=`gtftk select_by_key -k feature -v transcript -i simple_02.gtf -m transcript_id,gene_id,exon_id  -b -s "||"|  cut -f2 | sort -n| perl -npe 's/\\n/,/'`
      [ "$result" = "2,2,21,27,32,49,64,64,106,106,124,124,175,179,209," ]
    }
    
    # Select_by_key: test with a 'large' dataset
    @test "select_by_key_15" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -k feature -v transcript| wc -l`
      [ "$result" -eq 8531 ]
    }
    
    # Select_by_key: test with a 'large' dataset
    @test "select_by_key_16" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -k feature -v exon | wc -l`
      [ "$result" -eq 64251 ]
    }
        
    '''

    CmdObject(name="select_by_key",
              message="Select lines from a GTF based on attributes and "
                      "values.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              desc=__doc__,
              updated=__updated__,
              test=test)
