#!/usr/bin/env python
"""
 Find transcripts whose body/TSS/TTS region extended in 5' and 3'
 (-u/-d) overlaps with any transcript from another gene. Strandness is not
 considered by default. Used --invert-match to find those that do not overlap. If
 --annotate-gtf is used, all lines of the input GTF file will be printed and a new
 key containing the list of overlapping transcripts will be added to the transcript
 features/lines (key will be 'overlapping_*' with * one of body/TSS/TTS). The --annotate-gtf and
 --invert-match arguments are mutually exclusive.
"""
import argparse
import os
import sys
from collections import defaultdict

from pygtftk import arg_formatter
from pygtftk.arg_formatter import CheckChromFile
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-24"

__notes__ = '''
 -- -\-chrom-info may also accept 'mm8', 'mm9', 'mm10', 'hg19', 'hg38', 'rn3' or 'rn4'. In this case the 
 corresponding size of conventional chromosomes are used. ChrM is not used.  
'''


def make_parser():
    """The program parser."""
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

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Chromosome information. A tabulated two-columns"
                                 " file with chromosomes as column 1 and sizes as"
                                 " column 2",
                            default=None,
                            metavar="CHROMINFO",
                            action=CheckChromFile,
                            required=True)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the region in 5' by a given value (int)."
                                 " Used to define the region around the TSS/TTS.",
                            default=1500,
                            metavar="UPSTREAM",
                            type=int,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the region in 3' by a given value (int)."
                                 " Used to define the region around the TSS/TTS.",
                            default=1500,
                            metavar="DOWNSTREAM",
                            type=int,
                            required=False)

    parser_grp.add_argument('-t', '--feature-type',
                            help="The feature of interest.",
                            choices=['transcript', 'promoter', 'tts'],
                            default='transcript',
                            type=str,
                            required=False)

    parser_grp.add_argument('-s', '--same-strandedness',
                            help="Require same strandedness",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-S', '--diff-strandedness',
                            help="Require different strandedness",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-n', '--invert-match',
                            help="Not/Invert match.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-a', '--annotate-gtf',
                            help="All lines of the original GTF will be printed.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-k',
                            '--key-name',
                            type=str,
                            default=None,
                            help="The name of the key.",
                            required=False)

    parser_grp.add_argument('-b',
                            '--bool',
                            action="store_true",
                            help="When --annotate-gtf is used use 0 or 1 as key values (instead of overlapping transcripts id).",
                            required=False)

    parser_grp.add_argument('-@',
                            '--annotate-all',
                            action="store_true",
                            help="When --annotate-gtf annotate all transcripts (default value would be '0').",
                            required=False)

    return parser


def overlapping(
        inputfile=None,
        outputfile=None,
        key_name=None,
        upstream=1500,
        downstream=1500,
        chrom_info=None,
        feature_type='transcript',
        same_strandedness=False,
        diff_strandedness=False,
        annotate_gtf=False,
        bool=False,
        annotate_all=False,
        invert_match=False):
    """
Description: Find transcripts whose body/TSS/TTS do or do not overlap with any
transcript from another gene.
    """

    # ----------------------------------------------------------------------
    # Prepare key names
    # ----------------------------------------------------------------------

    if annotate_gtf:
        if key_name is None:
            key_info = ["overlap",
                        feature_type,
                        "u" + str(upstream / 1000) + "k",
                        "d" + str(downstream / 1000) + "k"
                        ]
            key_name = "_".join(key_info)

        if invert_match:
            message("--annotate-gtf and --invert-match are "
                    "mutually exclusive.",
                    type="ERROR")

    if same_strandedness and diff_strandedness:
        message("--same-strandedness and --diff-strandedness are "
                "mutually exclusive.",
                type="ERROR")

    message("Using -u " + str(upstream))
    message("Using -d " + str(downstream))

    overlapping_tx = defaultdict(list)

    # Load the GTF so that it won't be lost
    # if GTF stream comes from stdin
    gtf = GTF(inputfile)

    message("Getting transcript in bed format")

    tx_feat = gtf.select_by_key("feature",
                                "transcript")

    if annotate_all:
        overlapping_tx = gtf.extract_data(keys=["transcript_id"], as_dict=True, default_val="0")
        for i in overlapping_tx:
            overlapping_tx[i] = []

    # ----------------------------------------------------------------------
    # Get transcript limits
    # ----------------------------------------------------------------------

    tx_bed = tx_feat.to_bed(name=["transcript_id", "gene_id"], sep="||")

    message("Getting " + feature_type + " and 'slopping'.")

    if feature_type == "transcript":

        bed_obj = tx_bed.slop(s=True,
                              l=upstream,
                              r=downstream,
                              g=chrom_info.name).cut([0, 1, 2, 3, 4, 5])

    elif feature_type == "promoter":

        bed_obj = tx_feat.get_tss(name=["transcript_id", "gene_id"],
                                  sep="||").slop(s=True,
                                                 l=upstream,
                                                 r=downstream,
                                                 g=chrom_info.name).cut([0, 1,
                                                                         2, 3,
                                                                         4, 5])

    elif feature_type == "tts":

        bed_obj = tx_feat.get_tts(name=["transcript_id", "gene_id"],
                                  sep="||").slop(s=True,
                                                 l=upstream,
                                                 r=downstream,
                                                 g=chrom_info.name).cut([0, 1,
                                                                         2, 3,
                                                                         4, 5])
    else:
        message("Not implemented yet", type="ERROR")

    tmp_file = make_tmp_file(feature_type + "_slopped_region", ".bed")
    bed_obj.saveas(tmp_file.name)

    overlap_regions = bed_obj.intersect(tx_bed,
                                        wb=True,
                                        s=same_strandedness,
                                        S=diff_strandedness)

    tmp_file = make_tmp_file(feature_type + "_overlapping_regions", ".bed")
    overlap_regions.saveas(tmp_file.name)

    for i in overlap_regions:

        tx_other, gn_other = i.fields[9].split("||")
        tx_id, gene_id = i.fields[3].split("||")
        if gene_id != gn_other:
            overlapping_tx[tx_id] += [tx_other]

    if bool:
        for k, _ in overlapping_tx.items():
            if not len(overlapping_tx[k]):
                overlapping_tx[k] = "0"
            else:
                overlapping_tx[k] = "1"

    if not invert_match:

        if not annotate_gtf:
            value = ",".join(set(overlapping_tx.keys()))
            gtf.select_by_key("transcript_id",
                              value).write(outputfile,
                                           gc_off=True)
        else:

            if len(overlapping_tx):
                gtf = gtf.add_attr_from_dict(feat="transcript",
                                             key="transcript_id",
                                             a_dict=overlapping_tx,
                                             new_key=key_name)
            gtf.write(outputfile,
                      gc_off=True)

    else:
        values = ",".join(set(overlapping_tx.keys()))
        gtf.select_by_key("transcript_id",
                          values,
                          invert_match).write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    overlapping(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    # overlapping: load dataset
    @test "overlapping_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
        
    #overlapping: tts of G0005T001,G0010T001 are overlapping other tx.
    @test "overlapping_1" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t tts -u 2 -d 2 | gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0005T001,G0010T001," ]
    }
    
    
    #overlapping: tss/promoter of G0002T001,G0006T001,G0006T002, are overlapping other tx.
    @test "overlapping_2" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t promoter -u 2 -d 2 | gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0002T001,G0006T001,G0006T002," ]
    }
    
    #overlapping: there is no tss/promoter same overlap with a tx from anothe gene on the other strand.
    @test "overlapping_3" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t promoter -u 2 -d 2 -S | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #overlapping: tss/promoter of G0002T001,G0006T001,G0006T002, are overlapping other tx on the same strand.
    @test "overlapping_4" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t promoter -u 2 -d 2 -s | gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0002T001,G0006T001,G0006T002," ]
    }
    
    #overlapping: tts of G0005T001,G0010T001 are overlapping other tx on the same strand.
    @test "overlapping_5" {
     result=`gtftk overlapping -i simple.gtf  -c  simple.chromInfo -V 1 -t tts -u 2 -d 2 -s | gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0005T001,G0010T001," ]
    }
    
    
    #overlapping: there is no tts overlapping other tx on the opposite same strand.
    @test "overlapping_6" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t tts -u 2 -d 2 -S | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #overlapping: ensure the right number of exons are provided as output.
    @test "overlapping_7" {
     result=`gtftk overlapping -i simple.gtf  -c  simple.chromInfo -V 1 -t tts -u 2 -d 2 | gtftk nb_exons -f | grep G0005T001| cut -f 2`
      [ "$result" -eq 2 ]
    }
    
    #overlapping: ensure -a is working. Checking nb lines
    @test "overlapping_8" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo  -t tts -u 2 -d 2 -a| wc -l`
      [ "$result" -eq 70 ]
    }
    

    #overlapping: tts of G0005T001,G0010T001 are overlapping other tx on the same strand.
    @test "overlapping_9" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo  -t tts -u 2 -d 2  --annotate-gtf | grep overlap_tts | gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id| grep -v transcript_id| perl -npe 's/\\n/,/'`
      [ "$result" = "G0005T001,G0010T001," ]
    }

    #Check -b -@
    @test "overlapping_10" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t tts -u 2 -d 2 -k over -a -@ -b | gtftk tabulate -H -k transcript_id,over| awk '$2==1'| cut -f 1 | sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0005T001,G0010T001," ]
    }
        

    #Check -b -@
    @test "overlapping_11" {
     result=`gtftk overlapping -i simple.gtf -c  simple.chromInfo -V 1 -t tts -u 2 -d 2 -k over -a -@ -b | gtftk tabulate -H -k transcript_id,over| awk '$2==0'| cut -f 1 | sort | uniq| wc -l`
      [ "$result" -eq 13     ]
    }
        
    """

    CmdObject(name="overlapping",
              message="Find (non)overlapping transcripts.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="annotation",
              updated=__updated__,
              desc=__doc__,
              notes=__notes__,
              test=test)
