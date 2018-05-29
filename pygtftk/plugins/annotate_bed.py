#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys
from collections import OrderedDict

from pybedtools import BedTool

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import bedFile
from pygtftk.arg_formatter import checkChromFile
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-03-03"

__doc__ = """
 Annotate a list of BED files. The output is a text file with peak names and overlapped features in the GTF. 
"""

__notes__ = """
 -- Overlaps are reported without respect to strand.
 -- The program checks overlapping with included features (e.g. transcript, exon, CDS, gene...). Use convert_ensembl first if required. 
 -- You can ask more information about each feature property by using -m (e.g. gene_biotype or any additional key produced by gftk sub commands).

"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('bed_list',
                            help='A list of BED files (last argument).',
                            default=None,
                            metavar="BED",
                            type=bedFile(),
                            nargs='+')

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
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Cc][Oo][Vv]',
                                                                     '\.[Bb][Ee][Dd]')))

    parser_grp.add_argument('-l', '--labels',
                            help='BED file labels.',
                            default=None,
                            type=str,
                            metavar="LABELS",
                            required=False)

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Chromosome information. A tabulated two-columns"
                                 " file with chromosomes as column 1 and sizes as"
                                 " column 2",
                            default=None,
                            metavar="CHROMINFO",
                            action=checkChromFile,
                            required=True)

    parser_grp.add_argument('-m', '--name-column',
                            type=str,
                            default="transcript_id",
                            help="Use this ids (comma separated) to compute the name (4th column in "
                                 "bed output).",
                            required=False)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the region in 5' by a given value (int)."
                                 "Used to define the promoter regions and regions around the tts.",
                            default=1500,
                            metavar="UPSTREAM",
                            type=int,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the region in 3' by a given value (int)."
                                 "Used to define the promoter regions and regions around the tts.",
                            default=1500,
                            metavar="DOWNSTREAM",
                            type=int,
                            required=False)

    return parser


def annotate_bed(
        inputfile=None,
        outputfile=None,
        bed_list=None,
        upstream=None,
        labels=None,
        downstream=None,
        chrom_info=None,
        name_column=None,
        tmp_dir=None,
        logger_file=None,
        verbosity=0):
    """

    """

    name_column = name_column.split(",")

    # -------------------------------------------------------------------------
    # Create a list of labels.
    # Take user input into account
    # -------------------------------------------------------------------------

    bed_list = [x.name for x in bed_list]

    if len(bed_list) != len(set(bed_list)):
        message("Found the same bigwigs several times.",
                type="ERROR")

    message('Checking labels.')

    if labels is not None:
        labels = labels.split(",")
        # Ensure the number of labels is the same as the number of bw files.
        if len(labels) != len(bed_list):
            message("The number of labels should be the same as the number of"
                    " BED files.", type="ERROR")
        # Ensure labels are non-redondant
        if len(labels) > len(set(labels)):
            message("Labels must be unique.", type="ERROR")
    else:
        labels = []
        for i in range(len(bed_list)):
            labels += [
                os.path.splitext(
                    os.path.basename(
                        bed_list[i]))[0]]

    # ----------------------------------------------------------------------
    # Load the GTF
    # ----------------------------------------------------------------------

    gtf = GTF(inputfile)

    # ----------------------------------------------------------------------
    # Check intersect with all features
    # ----------------------------------------------------------------------

    target_bed_list = OrderedDict()

    message("Getting regions of interest...")

    target_bed_list['intergenic'] = gtf.get_intergenic(chrom_file=chrom_info,
                                                       upstream=0,
                                                       downstream=0,
                                                       chr_list=None,
                                                       feature_name="intergenic").sort()

    target_bed_list['intron_by_tx'] = gtf.get_introns(by_transcript=True,
                                                      name=name_column,
                                                      feat_name=True,
                                                      feat_name_last=True,
                                                      ).sort()

    target_bed_list['tss'] = gtf.get_tss(name=name_column,
                                         feature_name='tss').sort()

    target_bed_list['promoter'] = gtf.get_tss(name=name_column,
                                              feature_name='promoter').slop(s=True,
                                                                            l=upstream,
                                                                            r=downstream,
                                                                            g=chrom_info.name).sort()

    target_bed_list['around_tts'] = gtf.get_tts(name=name_column,
                                                feature_name='around_tts').slop(s=True,
                                                                                l=upstream,
                                                                                r=downstream,
                                                                                g=chrom_info.name).sort()

    target_bed_list['tts'] = gtf.get_tts(name=name_column,
                                         feature_name='tts').sort()

    for i in gtf.get_feature_list(nr=True):
        target_bed_list[i] = gtf.select_by_key("feature", i, 0).to_bed(name=name_column,
                                                                       add_feature_type=True)

    for bed_file, label in zip(bed_list, labels):
        for target in target_bed_list:
            bed_obj = BedTool(bed_file)
            overlap_regions = bed_obj.intersect(BedTool(target_bed_list[target]), wb=True)
            tmp_file = make_tmp_file(target + "_overlaps", ".bed")
            overlap_regions.saveas(tmp_file.name)

            for i in overlap_regions:
                i = chomp(label + "\t" + str(i))
                outputfile.write(i + "\n")

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    annotate_bed(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    @test "annotate_bed_1" {
     result=`gtftk get_example | gtftk annotate_bed -c pygtftk/data/simple/simple.chromInfo pygtftk/data/simple/simple_peaks.bed -u 2 -d 2 | grep -w peak_2| head -1| cut -f 9`
      [ "$result" = "intergenic" ]
    }



    @test "annotate_bed_2" {
     result=`gtftk get_example | gtftk annotate_bed -c pygtftk/data/simple/simple.chromInfo pygtftk/data/simple/simple_peaks.bed -u 2 -d 2 | grep -w peak_2| wc -l `
      [ "$result" -eq 2 ]
    }

    @test "annotate_bed_2" {
     result=`gtftk get_example | gtftk annotate_bed -c pygtftk/data/simple/simple.chromInfo pygtftk/data/simple/simple_peaks.bed -u 2 -d 2 | grep -w peak_1 | grep tss| wc -l`
      [ "$result" -eq 2 ]
    }


    """

    CmdObject(name="annotate_bed",
              message="Annotate a list of BED files.",
              parser=make_parser(),
              fun=annotate_bed,
              group="annotation",
              updated=__updated__,
              notes=__notes__,
              desc=__doc__,
              test=test)
