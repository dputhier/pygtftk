#!/usr/bin/env python
"""
 Select transcripts/gene overlapping a given locations.
"""

import argparse
import os
import re
import sys

from pybedtools.bedtool import BedTool

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
 -- A transcript is defined here as the genomic region from TSS to TTS including introns.
 -- This function will return the transcript and all its associated elements (exons, utr...)
 even if only a fraction (e.g intron) of the transcript is overlapping the feature.
 -- If -/-ft-type is set to 'gene' returns the gene and all its associated elements.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('-i', '--inputfile',
                        help="Path to the GTF file. Default to STDIN",
                        default=sys.stdin,
                        metavar="GTF",
                        type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser.add_argument('-o', '--outputfile',
                        help="Output file.",
                        default=sys.stdout,
                        metavar="GTF",
                        type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_mut = parser.add_mutually_exclusive_group(required=True)

    parser_mut.add_argument('-l', '--location',
                            help="List of chromosomal locations (chr:start-end[,chr:start-end]). 0-based",
                            default=None,
                            metavar="LOC",
                            type=str)

    parser_mut.add_argument('-f', '--location-file',
                            help="Bed file with chromosomal location.",
                            default=None,
                            metavar="BEDFILE",
                            type=argparse.FileType('r'))

    parser.add_argument('-t', '--ft-type',
                        help="The feature of interest.",
                        default="transcript",
                        choices=["transcript", "gene"],
                        type=str)

    parser.add_argument('-n', '--invert-match',
                        help='Not/invert match. Select transcript not overlapping.',
                        action="store_true")

    return parser


def select_by_loc(inputfile=None,
                  outputfile=None,
                  location=None,
                  ft_type=None,
                  invert_match=False,
                  location_file=None):
    """
 Select transcripts overlapping a given locations.
    """

    if invert_match:
        invert_match = 1
    else:
        invert_match = 0

    chroms = list()
    starts = list()
    ends = list()

    loc_n = 0

    if location_file is None:

        token = location.split(",")

        for i in range(len(token)):
            if not re.search(r".*:\d+\-\d+", token[i]):
                message("Check --location format. Should be chrom:start-end.",
                        type="ERROR")
            loc_n += 1
            chrom = [re.search("(.*):", token[i]).group(1)][0]
            start = [re.search(":(.*)-", token[i]).group(1)][0]
            end = [re.search("-(.*)$", token[i]).group(1)][0]

            message("Location #" + str(loc_n) +
                    " at :" +
                    chrom + ":" + start + "-" + end)

            chroms += [chrom]
            starts += [str(int(start) + 1)]
            ends += [end]
    else:
        for i in BedTool(location_file.name):
            chroms += [i.chrom]
            starts += [str(i.start + 1)]
            ends += [str(i.end)]

    chrom_str = ",".join(chroms)
    start_str = ",".join(starts)
    end_str = ",".join(ends)

    gtf = GTF(inputfile)

    id_list = gtf.select_by_key("feature",
                                ft_type).select_by_loc(chrom_str,
                                                       start_str,
                                                       end_str
                                                       ).extract_data(keys=ft_type + '_id',
                                                                      nr=True,
                                                                      as_list=True)

    if not invert_match:

        gtf.select_by_key(ft_type + '_id',
                          ",".join(id_list),
                          invert_match=invert_match).write(outputfile,
                                                           gc_off=True)
    else:
        if ft_type == 'transcript':

            gtf.select_by_key('transcript_id',
                              ",".join(id_list),
                              invert_match=invert_match,
                              no_na=True
                              ).write(outputfile,
                                      gc_off=True)
        else:
            gtf.select_by_key('gene_id',
                              ",".join(id_list),
                              invert_match=invert_match
                              ).write(outputfile,
                                      gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_loc(**args)


if __name__ == '__main__':

    main()

else:

    test = """

    #select_by_loc: copying datasets
    @test "select_by_loc_0" {
     result=`gtftk get_example -d simple -f "*"`
      [ "$result" = "" ]
    }
    
    #select_by_loc: various coordinate tests
    @test "select_by_loc_1" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:176-186 | gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0002T001,G0010T001," ]
    }
    
    #select_by_loc: various coordinate tests
    @test "select_by_loc_2" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:176-177 | gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0010T001," ]
    }
    
    #select_by_loc: just outside G0008T001
    @test "select_by_loc_3" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:223-262 | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #select_by_loc: just inside G0008T001
    @test "select_by_loc_4" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:222-261 |gtftk tabulate -FH -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "" ]
    }
    
    #select_by_loc: this region contains 4 transcripts
    @test "select_by_loc_5" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:103-142 |gtftk tabulate -H -k transcript_id| sort | uniq| perl -npe 's/\\n/,/'`
      [ "$result" = "G0001T001,G0001T002,G0007T001,G0007T002," ]
    }
    
    #select_by_loc: should support multiple regions
    @test "select_by_loc_6" {
     result=`gtftk get_example | gtftk select_by_key -k feature -v transcript | gtftk  select_by_loc -l chr1:10-15,chr1:40-50  | wc -l `
      [ "$result" -eq 4 ]
    }
    
    #select_by_loc: should support multiple regions from a file
    @test "select_by_loc_7" {
     result=`gtftk get_example | gtftk select_by_key -k feature -v transcript | gtftk  select_by_loc -f simple.loc_bed  | wc -l `
      [ "$result" -eq 4 ]
    }

    #select_by_loc: should support multiple regions from a file
    @test "select_by_loc_8" {
     result=`gtftk get_example | gtftk  select_by_loc -f simple.loc_bed| wc -l`
      [ "$result" -eq 14 ]
    }
    
    
    #select_by_loc:
    @test "select_by_loc_9" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:210-220 | wc -l`
      [ "$result" -eq 4 ]
    }

    #select_by_loc: should support multiple regions from a file
    @test "select_by_loc_10" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:24-25 | wc -l`
      [ "$result" -eq 7 ]
    }

    #select_by_loc: with -t gene
    @test "select_by_loc_11" {
     result=`gtftk select_by_loc -i simple.gtf -l chr1:24-25 -t gene| wc -l`
      [ "$result" -eq 13 ]
    }

    #select_by_loc:
    @test "select_by_loc_12" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:175-176   |   gtftk select_by_key -k feature -v transcript | wc -l`
      [ "$result" -eq 1 ]
    }
    
    #select_by_loc:
    @test "select_by_loc_13" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:174-175   |   gtftk select_by_key -k feature -v transcript | wc -l`
      [ "$result" -eq 0 ]
    }
    
    #select_by_loc: check -n is working

    @test "select_by_loc_14" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:175-176   -n | gtftk select_by_key -k transcript_id -v G0010T001  |  wc -l`
      [ "$result" -eq 0 ]
    }
    
    @test "select_by_loc_15" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:175-176   | gtftk select_by_key -k transcript_id -v G0010T001  |  wc -l`
      [ "$result" -eq 3 ]
    }
    
    @test "select_by_loc_16" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:175-176    -t gene  -n | wc -l`
      [ "$result" -eq 66 ]
    }

    @test "select_by_loc_17" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:100-120 | wc -l`
      [ "$result" -eq 6 ]
    }
                
    @test "select_by_loc_18" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:100-120 -n | wc -l`
      [ "$result" -eq 54 ]
    }
    

    @test "select_by_loc_19" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:100-120 -t gene | wc -l`
      [ "$result" -eq 7 ]
    }
    
        
    @test "select_by_loc_20" {
     result=`gtftk get_example | gtftk select_by_loc -l chr1:100-120 -t gene -n | wc -l`
      [ "$result" -eq 63 ]
    }
    
    """

    CmdObject(name="select_by_loc",
              message="Select transcript/gene overlapping a genomic feature.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              group="selection",
              notes=__notes__,
              desc=__doc__,
              test=test)
