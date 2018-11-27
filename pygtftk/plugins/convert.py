#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys
from builtins import str

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import is_comment

__updated__ = "2018-01-20"
__doc__ = """Convert a GTF to various format (still limited)."""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN.",
                            default=sys.stdin,
                            metavar="GTF",
                            required=False,
                            type=arg_formatter.gtf_rwb('r'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="BED/BED3/BED6",
                            type=arg_formatter.bed_rw('w'))

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name.",
                            default="gene_id,transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-m', '--more-names',
                            help="Add this information to the 'name' column of the BED file.",
                            default="",
                            type=str)

    parser_grp.add_argument('-f', '--format',
                            help='Currently one of bed3, bed6',
                            type=str,
                            choices=('bed', 'bed3', 'bed6'),
                            default='bed6',
                            required=False)

    return parser


def convert(inputfile=None,
            outputfile=None,
            format="bed",
            names="gene_id,transcript_id",
            separator="|",
            more_names=''):
    """
 Convert a GTF to various format.
    """

    if format == "bed3":
        gtf = GTF(inputfile, check_ensembl_format=False)

        for i in gtf.extract_data("seqid,start,end"):
            # First line is the header (#...)
            # Should be skipped
            if not is_comment(i[0]):
                # Should be zero based
                i[1] = str(int(i[1]) - 1)
                i.write(outputfile)

    elif format in ["bed", "bed6"]:
        gtf = GTF(inputfile, check_ensembl_format=False)

        nb_tokens = len(names.split(",")) + len(more_names.split(','))

        keys = "seqid,start,end," + names + ",score,strand"

        data = gtf.extract_data(keys)

        for i in data:
            i[1] = str(int(i[1]) - 1)

            if len(more_names):
                i[3:(2 + nb_tokens)] = [separator.join(more_names.split(',') +
                                                       i[3:(2 + nb_tokens)])]
            else:
                i[3:(2 + nb_tokens)] = [separator.join(i[3:(2 + nb_tokens)])]

            i.write(outputfile)

    """
    if format == "bed12":
        
        if inputfile.name == "<stdin>":
            tmp_file = make_tmp_file(suffix=".bed")
            with tmp_file as tf:
                for i in inputfile:
                    tf.write(i)
            tmp_file.close()
            inputfile = tmp_file

        gtf = GTF(inputfile.name)
        all_tx_id = gtf.get_tx_ids(nr=True)
        
        a_bedtool = BedTool(inputfile.name).sort()

        gtf = GTF(a_bedtool.fn)

        exons = gtf.select_by_key("feature", "exon")
        
        starts = exons.extract_data("transcript_id,start",
                                    zero_based=True,
                                    as_dict_of_merged_list=True)

        ends = exons.extract_data("transcript_id,end",
                                  as_dict_of_merged_list=True)

        info = gtf.select_by_key("feature",
                                 "transcript").extract_data("transcript_id,seqid,start,end,score,strand",
                                                            zero_based=True,
                                                            as_dict_of_lists=True)


        cds_start_dict = gtf.select_by_key("feature",
                                           "CDS").extract_data("transcript_id,start",
                                                         zero_based=True,
                                                         as_dict_of_merged_list=True)
        cds_end_dict = gtf.select_by_key("feature",
                                           "CDS").extract_data("transcript_id,end",
                                                         zero_based=True,
                                                         as_dict_of_merged_list=True)
                                           
        for tx_id in cds_start_dict:
            cds_start_dict[tx_id] = sorted(cds_start_dict[tx_id])[0]
            

        for tx_id in cds_end_dict:
            cds_end_dict[tx_id] = sorted(cds_end_dict[tx_id], reverse=True)[0]
            
            
        for tx_id in all_tx_id:
            
            exon_size = []
            for start, end in zip(starts[tx_id], ends[tx_id]):
                exon_size += [str(int(end) - int(start) + 1)]

            if tx_id not in cds_start_dict:
                cds_start_dict[tx_id] = "."
                cds_end_dict[tx_id] = "."

            token = [info[tx_id][0],
                     str(info[tx_id][1]),
                     str(info[tx_id][2]),
                     tx_id,
                     str(info[tx_id][3]),
                     str(info[tx_id][4]),
                     str(cds_start_dict[tx_id]),
                     str(cds_end_dict[tx_id]),
                     '255,0,0',
                     str(len(exon_size)),
                     ",".join(starts[tx_id]) + ",",
                     ",".join(ends[tx_id]) + ","]
            write_properly("\t".join(token), outputfile)
    """
    close_properly(outputfile, inputfile)


def main():
    """The main program."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    convert(**args)


if __name__ == '__main__':
    main()

else:

    test = '''

    # convert: load dataset
    @test "convert_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
            
    # Convert:
    @test "convert_1" {
     result=`gtftk convert -i simple.gtf | awk 'BEGIN{FS="\\t"}{print NF}'| sort | uniq`
      [ "$result" -eq 6 ]
    }
    
    
    # Convert: basic.
    @test "convert_2" {
     result=`gtftk convert -f bed3 -i simple.gtf | awk 'BEGIN{FS="\\t"}{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    # Convert: check name.
    @test "convert_3" {
     result=`gtftk convert -i simple.gtf -n gene_id,transcript_id,start | cut -f4| awk 'BEGIN{FS="|"}{print NF}'| sort | uniq`
      [ "$result" -eq 3 ]
    }
    
    # Convert: check zero based
    @test "convert_4" {
     result=`gtftk convert -i simple.gtf -n gene_id,transcript_id,start | cut -f2| head -n 1`
      [ "$result" -eq 124 ]
    }
    '''

    CmdObject(name="convert",
              message="Convert a GTF to various format including bed.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              desc=__doc__,
              group="conversion",
              test=test)
