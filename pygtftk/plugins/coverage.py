#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys

import pandas as pd
from builtins import range
from pybedtools import BedTool

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import checkChromFile
from pygtftk.bwig.bw_coverage import bw_cov_mp
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly, message
from pygtftk.utils import make_tmp_file

__updated__ = "2018-02-05"
__doc__ = """
Takes a GTF as input to compute bigwig coverage in regions of interest (promoter,
transcript body, intron, intron_by_tx, tts...) or a BED6 to focus on user-defined regions. If --n-highest
is used the program will compute the coverage of each bigwig based on the average value of
the n windows (--nb-window) with the highest coverage values.
 """
__notes__ = """
-- Regions were signal can be computed (if GTF file as input): promoter/tss, tts,
introns, intron_by_tx, intergenic regions or any feature available in the GTF file (transcript, exon, gene...).
-- If -\matrix-out is selected, the signal for each bigwig will be provided in a dedicated column. Otherwise, signal for each bigwig is provided through a dedicated line.
-- If bed is used as input, each region should have its own name (column 4).
"""


# -------------------------------------------------------------------------
# The main parser
# -------------------------------------------------------------------------


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)
    parser_grp = parser.add_argument_group('Arguments')

    # required but not required...

    parser_grp.add_argument(
        'bw_list',
        help='A list of Bigwig file (last argument).',
        type=FileWithExtension('r', valid_extensions=('\.[Bb][Ww]$',
                                                      '\.[Bb][Ii][Gg][Ww][Ii][Gg]$')),
        nargs='+')

    parser_grp.add_argument('-i', '--inputfile',
                            help="The input GTF/BED file. Only GTF file if <stdin> is used.",
                            default=sys.stdin,
                            metavar="GTF/BED",
                            required=False,
                            type=FileWithExtension('r',
                                                   valid_extensions=('\.[Gg][Tt][Ff](\.[Gg][Zz])?$',
                                                                     '\.[Bb][Ee][Dd]$',
                                                                     '\.[Bb][Ee][Dd]3$',
                                                                     '\.[Bb][Ee][Dd]6$')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Cc][Oo][Vv]')))
    parser_grp.add_argument(
        '-c', '--chrom-info',
        help="Tabulated two-columns file. Chromosomes"
             " as column 1 and sizes as"
             " column 2 ",
        default=None,
        metavar="CHROMINFO",
        action=checkChromFile,
        required=True)

    parser_grp.add_argument(
        '-u', '--upstream',
        help="Extend the regions in 5' by a given value (int).",
        default=0,
        metavar="UPSTREAM",
        type=int,
        required=False)

    parser_grp.add_argument(
        '-d', '--downstream',
        help="Extend the regions in 3' by a given value (int).",
        default=0,
        metavar="DOWNSTREAM",
        type=int,
        required=False)

    parser_grp.add_argument(
        '-w', '--nb-window',
        type=int,
        default=1,
        help='Split the region into w bins (see -n).',
        required=False)

    parser_grp.add_argument(
        '-k', '--nb-proc',
        type=int,
        default=1,
        help='Use this many threads to compute coverage.',
        required=False)

    parser_grp.add_argument(
        '-f', '--ft-type',
        type=str,
        default="promoter",
        help="Region in which coverage is to be computed (promoter, intron, intergenic, tts or any feature defined in the column 3 of the GTF).",
        required=False)

    parser_grp.add_argument(
        '-l', '--labels',
        help='Bigwig labels.',
        default=None,
        type=str,
        required=False)

    parser_grp.add_argument(
        '-m', '--name-column',
        type=str,
        default="transcript_id",
        help="Use this ids to compute the name (4th column in "
             "bed output).",
        required=False)

    parser_grp.add_argument(
        '-p',
        '--pseudo-count',
        type=int,
        default=1,
        help='A pseudo-count to add in case count is equal to 0.')

    parser_grp.add_argument(
        '-n', '--n-highest',
        type=int,
        default=None,
        help='For each bigwig, use the n windows with higher values to compute coverage.',
        required=False)

    parser_grp.add_argument(
        '-x',
        '--matrix-out',
        action="store_true",
        help='Matrix output format. Bigwigs as column names features as rows.',
        required=False)

    parser_grp.add_argument(
        '-zn', '--zero-to-na',
        help='Use NA not zero when region is undefined in bigwig or below window size.',
        action='store_true',
        required=False)

    parser_grp.add_argument(
        '-a',
        '--key-name',
        type=str,
        default="cov",
        help="If gtf format is requested, the name of the key.",
        required=False)

    parser_grp.add_argument(
        '-s', '--stat',
        type=str,
        choices=['mean', 'sum'],
        default='mean',
        help='The statistics to be computed for each region.',
        required=False)

    return parser


# -------------------------------------------------------------------------
# The main function
# -------------------------------------------------------------------------

def coverage(
        inputfile=None,
        outputfile=None,
        bw_list=None,
        labels=None,
        pseudo_count=1,
        nb_window=1,
        ft_type="promoter",
        n_highest=None,
        downstream=1000,
        key_name="cov",
        zero_to_na=False,
        name_column=None,
        upstream=1000,
        chrom_info=None,
        nb_proc=1,
        matrix_out=False,
        tmp_dir=None,
        logger_file=None,
        stat='mean',
        verbosity=True):
    """
    Compute transcript coverage with one or several bigWig.
    """

    # -------------------------------------------------------------------------
    # Create a list of labels.
    # Take user input in account
    # -------------------------------------------------------------------------

    bw_list = [x.name for x in bw_list]

    if len(bw_list) != len(set(bw_list)):
        message("Found the same bigwigs several times.",
                type="ERROR")

    message('Checking labels.')

    if labels is not None:
        labels = labels.split(",")
        # Ensure the number of labels is the same as the number of bw files.
        if len(labels) != len(bw_list):
            message("The number of labels should be the same as the number of"
                    " bigwig files.", type="ERROR")
        # Ensure labels are non-redondant
        if len(labels) > len(set(labels)):
            message("Labels must be unique.", type="ERROR")
    else:
        labels = []
        for i in range(len(bw_list)):
            labels += [
                os.path.splitext(
                    os.path.basename(
                        bw_list[i]))[0]]

    # -------------------------------------------------------------------------
    # Check the number of windows
    #
    # -------------------------------------------------------------------------

    if n_highest is None:
        n_highest = nb_window

    message('Number of bins: %d' % nb_window)
    message('N highest values: %d' % n_highest)

    if n_highest > nb_window:
        message('The number of window used for computing the score'
                ' (-n) can not be greater than the number of'
                ' windows (-w)', type="ERROR")
        sys.exit()

    # -------------------------------------------------------------------------
    # Check input file is in bed or GTF format
    #
    # -------------------------------------------------------------------------

    message("Loading input file...")
    if inputfile.name == '<stdin>':
        gtf = GTF(inputfile.name)
        is_gtf = True
    else:
        region_bo = BedTool(inputfile.name)
        if len(region_bo) == 0:
            message("Unable to find requested regions",
                    type="ERROR")

        if region_bo.file_type == 'gff':
            gtf = GTF(inputfile.name)
            is_gtf = True
        else:
            is_gtf = False

    # -------------------------------------------------------------------------
    # Get regions of interest
    #
    # -------------------------------------------------------------------------

    if is_gtf:

        message("Getting regions of interest...")

        if ft_type == "intergenic":

            region_bo = gtf.get_intergenic(chrom_info, 0, 0).slop(s=True,
                                                                  l=upstream,
                                                                  r=downstream,
                                                                  g=chrom_info.name).sort()

        elif ft_type == "intron":

            region_bo = gtf.get_introns().slop(s=True,
                                               l=upstream,
                                               r=downstream,
                                               g=chrom_info.name).sort()

        elif ft_type == "intron_by_tx":

            region_bo = gtf.get_introns(by_transcript=True
                                        ).slop(s=True,
                                               l=upstream,
                                               r=downstream,
                                               g=chrom_info.name).sort()

        elif ft_type == "promoter":

            region_bo = gtf.get_tss(name=["transcript_id"]).slop(s=True,
                                                                 l=upstream,
                                                                 r=downstream,
                                                                 g=chrom_info.name).sort()

        elif ft_type == "tts":

            region_bo = gtf.get_tts(name=["transcript_id"]).slop(s=True,
                                                                 l=upstream,
                                                                 r=downstream,
                                                                 g=chrom_info.name).sort()

        else:
            name_column = name_column.split(",")
            region_bo = gtf.select_by_key(
                "feature",
                ft_type, 0
            ).to_bed(name=name_column).slop(s=True,
                                            l=upstream,
                                            r=downstream,
                                            g=chrom_info.name).sort()

        if len(region_bo) == 0:
            message("Unable to find requested regions",
                    type="ERROR")

    else:
        region_bo = region_bo.slop(s=True,
                                   l=upstream,
                                   r=downstream,
                                   g=chrom_info.name).sort()

    region_bed = make_tmp_file(prefix="region", suffix=".bed")

    region_bo.saveas(region_bed.name)

    # -------------------------------------------------------------------------
    # Compute coverage
    #
    # -------------------------------------------------------------------------

    result_bed = bw_cov_mp(bw_list=bw_list,
                           region_file=open(region_bed.name),
                           labels=labels,
                           bin_nb=nb_window,
                           pseudo_count=pseudo_count,
                           zero_to_na=zero_to_na,
                           nb_proc=nb_proc,
                           n_highest=n_highest,
                           stat=stat,
                           verbose=verbosity)

    if matrix_out:
        result_bed.close()

        df_first = pd.read_csv(result_bed.name, sep="\t", header=None)

        df_first = df_first.ix[:, [0, 1, 2, 3, 5, 4]]
        """
        bwig = [z.split("|")[0] for z in df_first[3]]
        bwig_nr = list(set(bwig))
        """

        df_list = []
        """
        for i in range(len(bwig_nr)):
        """
        for i in range(len(labels)):
            # create a sub data frame containing the coverage values of the
            # current bwig
            str_to_find = r"^" + labels[i] + r"\|"
            tmp_df = df_first[df_first[3].str.match(str_to_find)].copy()
            to_replace = r"^" + labels[i] + r"\|"
            tmp_df.iloc[:, 3] = tmp_df.iloc[:, 3].replace(to_replace,
                                                          r"", regex=True)

            df_list += [tmp_df]

        df_final = df_list.pop(0)

        for i in df_list:
            # Add columns to df final by joining on
            # chrom, start, end, transcript_id, strand
            df_final = df_final.merge(i.iloc[:,
                                      list(range(6))], on=[0, 1,
                                                           2, 3, 5])

        df_final.columns = ["chrom",
                            "start",
                            "end",
                            "name",
                            "strand"] + labels

        df_final.to_csv(outputfile, sep="\t", index=False)

    else:
        nb_line = 0
        for i in result_bed:
            outputfile.write(i)
            nb_line += 1

        if nb_line == 0:
            message("No line available in output...",
                    type="ERROR")

    close_properly(inputfile, outputfile)


# -------------------------------------------------------------------------
# Call  to main
# -------------------------------------------------------------------------


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    coverage(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    
    #coverage: test coverage of tts
    @test "coverage_1" {
     result=`gtftk coverage -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -f tts  -p 0  pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "1.0,1.0,3.0,2.0,3.0,2.0,0.0,0.0,2.0,2.0,1.0,1.0,1.0,2.0,1.0," ]
    }
    
    #coverage: test coverage of tss
    @test "coverage_2" {
     result=`gtftk coverage -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0 pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.0,0.0,2.0,2.0,2.0,0.0,0.0,0.0,1.0,1.0,4.0,4.0,0.0,0.0,1.0," ]
    }
    
    
    #coverage: test coverage of tss (-u, -d)
    @test "coverage_3" {
     result=`gtftk coverage -u 3 -d 1 -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0  pygtftk/data/simple/simple.bw| cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.8,0.8,2.2,2.2,2.0,0.0,0.0,0.0,1.6,1.6,4.0,4.0,0.4,0.0,0.8," ]
    }
    
    #coverage: test coverage of intron
    @test "coverage_4" {
     result=`gtftk coverage -f intron -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0 pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "3.0,1.0,2.666667,0.0,0.0,0.0,2.0," ]
    }
    
    
    #coverage: test coverage of intergenic
    @test "coverage_5" {
     result=`gtftk coverage -f intergenic -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0 pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "1.0,1.857143,2.0,0.0,0.866667,3.5,1.72973,2.15,1.230769,0.0," ]
    }
    
    #coverage: transcript
    @test "coverage_6" {
     result=`gtftk coverage -f transcript -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0 pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.666667,0.666667,2.5,2.0,2.666667,0.333333,0.0,0.0,1.0,1.0,2.142857,2.142857,0.181818,0.8,1.615385," ]
    }
    
    #coverage: test coverage of promoter
    @test "coverage_7" {
     result=`gtftk coverage -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0  pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.0,0.0,2.0,2.0,2.0,0.0,0.0,0.0,1.0,1.0,4.0,4.0,0.0,0.0,1.0," ]
    }
    
    #coverage: test coverage of peaks
    @test "coverage_8" {
     result=`gtftk coverage -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo  -f tts  -p 0 -K toto -f transcript -m transcript_id,gene_id,exon_id pygtftk/data/simple/simple.bw | cut -f5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.0,2.0,3.333333,2.666667,2.666667,0.0," ]
    }
    
    #coverage: Same results should be obatined with <stdin>
    @test "coverage_9" {
     result=`gtftk coverage -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -p 0  pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.0,0.0,2.0,2.0,2.0,0.0,0.0,0.0,1.0,1.0,4.0,4.0,0.0,0.0,1.0," ]
    }
    
    #coverage: Same results should be obatined with <stdin>
    @test "coverage_10" {
     result=`cat pygtftk/data/simple/simple.gtf | gtftk coverage -c pygtftk/data/simple/simple.chromInfo  -p 0  pygtftk/data/simple/simple.bw | cut -f 5| perl -npe 's/\\n/,/'`
      [ "$result" = "0.0,0.0,2.0,2.0,2.0,0.0,0.0,0.0,1.0,1.0,4.0,4.0,0.0,0.0,1.0," ]
    }
    
    #coverage: Using peaks we expect 6 lines
    @test "coverage_11" {
     result=`gtftk coverage -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo  -f tts  -p 0 -K toto  pygtftk/data/simple/simple.bw| wc -l`
      [ "$result" -eq 6 ]
    }
    
    #coverage: If two bigwig are provided we expect twice lines (bed)
    @test "coverage_12" {
     result=`gtftk get_example -f 2.bw 2>/dev/null; gtftk coverage -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo  -f tts  -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw| wc -l`
      [ "$result" -eq 12 ]
    }
    
    #coverage: If two bigwig are provided we expect twice lines (gtf)
    @test "coverage_13" {
     result=`gtftk coverage -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo  -f tts  -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw | wc -l`
      [ "$result" -eq 30 ]
    }
    
    
    #coverage: check label is working (1)
    @test "coverage_14" {
     result=`gtftk coverage -i pygtftk/data/simple/simple.gtf -c pygtftk/data/simple/simple.chromInfo    -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw -l s1,s2 | cut -f4 | sed 's/|.*//'| sort | uniq -c| perl -npe 's/ +//g; s/\\n/,/'`
      [ "$result" = "15s1,15s2," ]
    }
    
    #coverage: check label is working (2)
    @test "coverage_15" {
     result=`gtftk coverage -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo    -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw -l s1,s2 | cut -f4 | sed 's/|.*//'| sort | uniq -c| perl -npe 's/ +//g; s/\\n/,/'`
      [ "$result" = "6s1,6s2," ]
    }
    
    #coverage: Using peaks (as bed3) we expect 6 lines (here * 2)
    @test "coverage_16" {
     result=`gtftk coverage -i pygtftk/data/simple/simple_peaks.bed3 -c pygtftk/data/simple/simple.chromInfo    -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw -l s1,s2| wc -l`
      [ "$result" -eq 12 ]
    }
    
    #coverage: here we expect 18 lines
    @test "coverage_17" {
     result=`cp simple.2.bw simple.3.bw; gtftk coverage -i pygtftk/data/simple/simple_peaks.bed3 -c pygtftk/data/simple/simple.chromInfo    -p 0 -K toto  pygtftk/data/simple/simple.bw simple.3.bw  pygtftk/data/simple/simple.2.bw -l s1,s2,s3 | wc -l`
      [ "$result" -eq 18 ]
    }
    
    #coverage: 7 lines for introns
    @test "coverage_18" {
     result=`gtftk get_example|  gtftk coverage -c pygtftk/data/simple/simple.chromInfo  -f intron -u 0 -d 0  pygtftk/data/simple/simple.bw| wc -l`
      [ "$result" -eq 7 ]
    }
    
    #coverage: 10 lines for intergenic
    @test "coverage_19" {
     result=`gtftk get_example|  gtftk coverage -c pygtftk/data/simple/simple.chromInfo  -f intergenic -u 0 -d 0  pygtftk/data/simple/simple.bw| wc -l`
      [ "$result" -eq 10 ]
    }
    
    #coverage: 20 lines for CDS
    @test "coverage_20" {
     result=`gtftk get_example|  gtftk coverage -c pygtftk/data/simple/simple.chromInfo  -f CDS -u 0 -d 0  pygtftk/data/simple/simple.bw| wc -l`
      [ "$result" -eq 20 ]
    }
    
    #coverage: test -s 
    @test "coverage_21" {
     result=`gtftk coverage -s sum -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo  -f tts  -p 0 -K toto  pygtftk/data/simple/simple.bw | sortBed |  cut -f5 | perl -npe 's/\\n/,/'`
      [ "$result" = "0.0,6.0,10.0,8.0,8.0,0.0," ]
    }

    #coverage: test -s 
    @test "coverage_22" {
     result=`gtftk coverage -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo    -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw  -l s1,s2 -x | cut -f7 | perl -npe 's/\\n/,/'`
      [ "$result" = "s2,0.0,2.0,3.333333,2.666667,2.666667,0.0," ]
    }
      
            
    #coverage: test -s 
    @test "coverage_23" {
     result=`gtftk coverage -i pygtftk/data/simple/simple_peaks.bed -c pygtftk/data/simple/simple.chromInfo    -p 0 -K toto  pygtftk/data/simple/simple.bw pygtftk/data/simple/simple.2.bw  -l s1,s2 -x | cut -f6 | perl -npe 's/\\n/,/'`
      [ "$result" = "s1,0.0,2.0,3.333333,2.666667,2.666667,0.0," ]
    }
          
    """

    CmdObject(name='coverage',
              message='Compute bigwig coverage in body, promoter, tts...',
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              notes=__notes__,
              updated=__updated__,
              group="coverage",
              test=test)
