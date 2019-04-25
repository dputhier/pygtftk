#!/usr/bin/env python
"""
 Create a matrix storing the bigwig coverage computed from binned regions.
"""

import argparse
import os
import sys
import zipfile

import pandas as pd
import pyBigWig
from pybedtools import BedTool

import pygtftk
from pygtftk import arg_formatter
from pygtftk.arg_formatter import CheckChromFile
from pygtftk.bwig.bw_coverage import bw_profile_mp
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
 -- -\-chrom-info may also accept 'mm8', 'mm9', 'mm10', 'hg19', 'hg38', 'rn3' or 'rn4'. In this case the 
 corresponding size of conventional chromosomes are used. ChrM is not used.  
"""


# NEED TO CHECK wether region names are uniq. Or force them to be ...

# -------------------------------------------------------------------------
# First define the function/command arguments.
# -------------------------------------------------------------------------


def make_parser():
    """The main parser."""

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="A GTF file or bed file. A GTF if <stdin>.",
                            default=sys.stdin,
                            metavar="GTF/BED",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('bed', 'gtf', 'gtf.gz')))

    parser_grp.add_argument('-y', '--bigwiglist',
                            help='A list of Bigwig files (last argument).',
                            type=arg_formatter.FormattedFile(mode='r', file_ext='bigwig'),
                            nargs='+',
                            required=True)

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file name (.zip extension will be added).",
                            default=sys.stdout,
                            metavar="GTF/TXT",
                            required=True,
                            type=argparse.FileType('w'))

    parser_grp.add_argument('-l', '--labels',
                            help='Bigwig labels (i.e short name version for '
                                 'plotting).',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-t', '--ft-type',
                            help='If input is a GTF, the region to analyse.',
                            default='promoter',
                            choices=['promoter',
                                     'tts',
                                     'transcript',
                                     'user_regions',
                                     'single_nuc'],
                            type=str,
                            required=False)

    parser_grp.add_argument('-p', '--pseudo-count',
                            help='Pseudo-count to add to all values.',
                            default=0,
                            type=arg_formatter.ranged_num(lowest=0, highest=None,
                                                          val_type="float", linc=True),
                            required=False)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the region of interest in 5' by a given value.",
                            default=1000,
                            type=arg_formatter.ranged_num(lowest=0, highest=None,
                                                          val_type="int", linc=True),
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the region of interest in 3' by a given value.",
                            default=1000,
                            type=arg_formatter.ranged_num(lowest=0, highest=None,
                                                          val_type="int", linc=True),
                            required=False)

    parser_grp.add_argument('-c', '--chrom-info',
                            help='Tabulated file (chr as '
                                 'column 1, sizes as column 2.)',
                            default=None,
                            action=CheckChromFile,
                            required=True)

    parser_grp.add_argument('-w', '--bin-nb',
                            help='Split the region into w bins.',
                            default=100,
                            type=arg_formatter.ranged_num(lowest=1, highest=None,
                                                          val_type="int", linc=True),
                            required=False)

    parser_grp.add_argument('-k', '--nb-proc',
                            help='Use this many threads to compute coverage.',
                            default=1,
                            type=arg_formatter.ranged_num(lowest=1, highest=None,
                                                          val_type="int", linc=True),
                            required=False)

    parser_grp.add_argument('-b', '--bin-around-frac',
                            help="Fraction of bins used in 5' and 3' regions.",
                            default=0.1,
                            type=float,
                            required=False)

    parser_grp.add_argument('-zn', '--zero-to-na',
                            help='Use NA not zero when region is undefined in bigwig.',
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-nst', '--no-stranded',
                            help='The bins should not be oriented relative to strand.',
                            action='store_true',
                            required=False)

    return parser


# -------------------------------------------------------------------------
# Now we declare a main function
# -------------------------------------------------------------------------

def mk_matrix(
        inputfile=None,
        outputfile=None,
        bigwiglist=None,
        ft_type=None,
        pseudo_count=0,
        upstream=1000,
        downstream=1000,
        bin_around_frac=0.1,
        chrom_info=None,
        bin_nb=100,
        nb_proc=None,
        labels=None,
        no_stranded=False,
        zero_to_na=False):
    """
 Description: Create a matrix to be used by 'profile' and 'heatmap' commands.
    """

    # -------------------------------------------------------------------------
    # Check argument consistency
    #
    # -------------------------------------------------------------------------

    if ft_type in ['single_nuc', 'promoter', 'tts']:
        region_size = upstream + downstream + 1
        if region_size < bin_nb:
            message("The region (-u/-d) needs to be extended given the number "
                    "of bins (--bin-nb)",
                    type="ERROR")

    # -------------------------------------------------------------------------
    # Check output file name does not ends with .zip
    #
    # -------------------------------------------------------------------------

    if outputfile.name.endswith(".zip"):
        outfn = outputfile.name.replace(".zip", "")
        outputfile = open(outfn, "w")

    # -------------------------------------------------------------------------
    # Check input file is in bed or GTF format
    #
    # -------------------------------------------------------------------------

    message("Loading input file...")
    if inputfile.name == '<stdin>':
        gtf = GTF(inputfile.name)
        is_gtf = True
        if ft_type == 'user_regions':
            message("--ft-type can not be set to user_regions"
                    " when a gtf is provided.", type="ERROR")
    else:
        try:

            region_bo = BedTool(inputfile.name)
            len(region_bo)
        except IndexError:
            message("Unable to read the input file. Check format",
                    type="ERROR")
        if len(region_bo) == 0:
            message("Unable to find requested regions",
                    type="ERROR")

        if region_bo.file_type == 'gff':
            message('Loading the GTF file.')
            gtf = GTF(inputfile.name)
            is_gtf = True
        else:
            is_gtf = False

            if ft_type != 'user_regions' and ft_type != 'single_nuc':
                message("Set --ft-type to 'user_regions' or 'single_nuc'"
                        " when using input bed file.",
                        type="ERROR")
            # Check that the strand is provided and
            # check it is located in the right column
            # (not checked by BedTool...).
            if region_bo.field_count() < 6:
                if not no_stranded:
                    message("Strand is undefined. Use -nst.", type="ERROR")
            else:
                region_name = dict()
                for i in region_bo:
                    if region_name.get(i.name, None) is None:
                        region_name[i.name] = 1
                    else:
                        message("Regions in bed file should have "
                                "unique identifier (col 4).",
                                type="ERROR")
                    if i.strand[0] not in ['.', '+', '-']:
                        message("Strand should be one of '+','-' or '.'.",
                                type="ERROR")
                    if ft_type == 'single_nuc':
                        if i.end - i.start != 1:
                            message("Region length should be 1 nucleotide "
                                    "long when 'single_nuc' is set. Use 'user_regions'.",
                                    type="ERROR")
                    elif ft_type == 'user_regions':
                        if i.end - i.start == 1:
                            message("Region length should not be 1 nucleotide "
                                    "long when 'user_regions' is set. Use 'single_nuc'.",
                                    type="ERROR")

    # -------------------------------------------------------------------------
    # Create a list of labels for the diagrams.
    # Take user input in account
    # -------------------------------------------------------------------------
    message('Checking labels.')

    if labels is not None:
        labels = labels.split(",")
        # Ensure the number of labels is the same as the number of bw files.
        if len(labels) != len(bigwiglist):
            message("The number of labels should be the same as the number of"
                    " bigwig files.", type="ERROR")
        # Ensure labels are non-redondant
        if len(labels) > len(set(labels)):
            message("Labels must be unique.", type="ERROR")
    else:
        labels = []
        for i in range(len(bigwiglist)):
            labels += [
                os.path.splitext(
                    os.path.basename(
                        bigwiglist[i].name))[0]]

    # -------------------------------------------------------------------------
    #
    # Get the requested transcrit lines in bed format
    # Tx are restricted to those found on chromosome
    # declared in the bigwig file.
    # -------------------------------------------------------------------------
    message('Getting the list of chromosomes declared in bigwig files.')
    bw_chrom = list()
    for i in bigwiglist:
        bw_chrom += list(pyBigWig.open(i.name).chroms().keys())

    bed_col = [0, 1, 2, 3, 4, 5]

    if is_gtf:

        message('Selecting chromosomes declared in bigwig from gtf.')
        tmp = gtf.select_by_key("feature",
                                "transcript"
                                ).select_by_key("seqid",
                                                ",".join(bw_chrom))

        tmp = gtf.select_by_key("feature",
                                "transcript")
        tmp_tx_name = tmp.extract_data("transcript_id", as_list=True)

        # If several trancript records are associated to
        # the same transcript_id, raise an error.
        if len(tmp_tx_name) > len(set(tmp_tx_name)):
            message('Transcripts should have a unique identifier.',
                    type="ERROR")

        message('Selecting requested regions.')

        # ----------------------------------------------------------------------
        #
        # Slop tss and promoters.
        # No need if transcript was requested (it will be flanked by upstream
        # and doswnstream regions later on).
        # ----------------------------------------------------------------------

        if ft_type == 'transcript':
            message("Getting transcript boundaries (input gtf).")

            main_region_bo = tmp.to_bed(name=["transcript_id"])

        elif ft_type == 'promoter':

            message("Getting promoter regions [-%d,+%d]." % (upstream,
                                                             downstream))

            main_region_bo = tmp.get_tss(name=["transcript_id"]).slop(s=True,
                                                                      l=upstream,
                                                                      r=downstream,
                                                                      g=chrom_info.name)

        elif ft_type == 'tts':

            main_region_bo = tmp.get_tts(name=["transcript_id"]).slop(s=True,
                                                                      l=upstream,
                                                                      r=downstream,
                                                                      g=chrom_info.name)

    else:
        message("Loading regions")

        if ft_type == 'user_regions':
            main_region_bo = BedTool(inputfile.name).cut(bed_col)
        elif ft_type == 'single_nuc':
            main_region_bo = BedTool(inputfile.name).cut(bed_col).slop(s=True,
                                                                       l=upstream,
                                                                       r=downstream,
                                                                       g=chrom_info.name)
        else:
            message("Unknown method.")

    # Save for tracability
    main_region_bed = make_tmp_file(prefix="region" + ft_type, suffix=".bed")
    main_region_bo.saveas(main_region_bed.name)

    # -------------------------------------------------------------------------
    #
    # Print a header in the output file
    #
    # -------------------------------------------------------------------------
    message("Preparing comments")

    comments = "#"
    comments += "ft_type:" + ft_type + ";"
    comments += "from:" + str(upstream) + ";"
    comments += "to:" + str(downstream) + ";"
    comments += "labels:" + ",".join(labels) + ";"

    # -------------------------------------------------------------------------
    # Compute coverage of requested region
    # Each worker will send a file
    # -------------------------------------------------------------------------

    outputfile_list = {}
    message("Using %d bins for main region." % bin_nb)

    tmp_file = bw_profile_mp(in_bed_file=main_region_bed.name,
                             nb_proc=nb_proc,
                             big_wig=[x.name for x in bigwiglist],
                             bin_nb=bin_nb,
                             pseudo_count=pseudo_count,
                             stranded=not no_stranded,
                             type="main",
                             labels=labels,
                             outputfile=outputfile.name,
                             zero_to_na=zero_to_na,
                             verbose=pygtftk.utils.VERBOSITY)

    outputfile_list["main"] = tmp_file

    # -------------------------------------------------------------------------
    # If transcript was requested
    # we must process flanking regions
    # We need to retrieve coverage of promoter [-upstream, 0]
    # as transcript coverage window size will depend on transcript length.
    # For promoter the length of windows will be fixed.
    # -------------------------------------------------------------------------

    if ft_type in ['transcript', 'user_regions']:

        # Number of bins for TTS and TSS
        around_bin_nb = int(round(bin_nb * bin_around_frac))
        if around_bin_nb < 1:
            around_bin_nb = 1

        if upstream > 0:

            if ft_type == 'transcript':
                message("Getting promoter (using %d bins)." % around_bin_nb)
                ups_region_bo = tmp.get_tss(name=["transcript_id"]
                                            ).slop(s=True,
                                                   l=upstream,
                                                   r=-1,
                                                   g=chrom_info.name).cut(bed_col)

            else:
                message("Getting upstream regions (%d bins)." % around_bin_nb)
                ups_region_bo = main_region_bo.flank(s=True,
                                                     l=upstream,
                                                     r=0,
                                                     g=chrom_info.name)

            upstream_bed_file = make_tmp_file(
                prefix="upstream_region" + ft_type,
                suffix=".bed")

            ups_region_bo.saveas(upstream_bed_file.name)

            tmp_file = bw_profile_mp(in_bed_file=upstream_bed_file.name,
                                     nb_proc=nb_proc,
                                     big_wig=[
                                         x.name for x in bigwiglist],
                                     bin_nb=around_bin_nb,
                                     pseudo_count=pseudo_count,
                                     stranded=not no_stranded,
                                     type="upstream",
                                     labels=labels,
                                     outputfile=outputfile.name,
                                     zero_to_na=zero_to_na,
                                     verbose=pygtftk.utils.VERBOSITY)

            outputfile_list["upstream"] = tmp_file

        if downstream > 0:

            if ft_type == 'transcript':
                message("Getting TTS (using %d bins)." % around_bin_nb)
                dws_region_bo = tmp.get_tts(name=["transcript_id"]
                                            ).slop(s=True,
                                                   l=-1,
                                                   r=downstream,
                                                   g=chrom_info.name).cut(bed_col)
            else:
                message(
                    "Getting downstream regions (%d bins)." %
                    around_bin_nb)

                dws_region_bo = main_region_bo.flank(s=True,
                                                     l=0,
                                                     r=downstream,
                                                     g=chrom_info.name)
            dws_bed_file = make_tmp_file(prefix="dowstream_region" + ft_type,
                                         suffix=".bed")

            dws_region_bo.saveas(dws_bed_file.name)

            tmp_file = bw_profile_mp(in_bed_file=dws_bed_file.name,
                                     nb_proc=nb_proc,
                                     big_wig=[
                                         x.name for x in bigwiglist],
                                     bin_nb=around_bin_nb,
                                     pseudo_count=pseudo_count,
                                     stranded=not no_stranded,
                                     type="downstream",
                                     labels=labels,
                                     outputfile=outputfile.name,
                                     zero_to_na=zero_to_na,
                                     verbose=pygtftk.utils.VERBOSITY)

            outputfile_list["downstream"] = tmp_file

    # -------------------------------------------------------------------------
    #
    # Merge file using pandas
    #
    # -------------------------------------------------------------------------

    message("Reading (pandas): " + outputfile_list["main"].name, type="DEBUG")
    df_main = pd.read_csv(outputfile_list["main"].name, sep="\t")
    # save strand and end
    # They will re-joined added later
    df_copy = df_main[['bwig', 'chrom', 'gene',
                       'strand', 'start', 'end']]

    df_start = df_main.pop('start')
    df_end = df_main.pop('end')

    if "upstream" in outputfile_list:
        message("Merging upstream file")
        message(
            "Reading (pandas): " +
            outputfile_list["upstream"].name,
            type="DEBUG")
        df_up = pd.read_csv(outputfile_list["upstream"].name, sep="\t")
        df_up = df_up.drop(['start', 'end'], 1)
        df_main = df_up.merge(df_main.loc[:, df_main.columns], on=['bwig',
                                                                   'chrom',
                                                                   'gene',
                                                                   'strand'])

    if "downstream" in outputfile_list:
        message("Merging downstream file")
        message(
            "Reading (pandas): " +
            outputfile_list["downstream"].name,
            type="DEBUG")
        df_dws = pd.read_csv(outputfile_list["downstream"].name, sep="\t")
        df_dws = df_dws.drop(['start', 'end'], 1)
        df_main = df_main.merge(df_dws.loc[:, df_dws.columns], on=['bwig',
                                                                   'chrom',
                                                                   'gene',
                                                                   'strand'])

    # join start and end.
    df_main = df_main.merge(df_copy.loc[:, df_copy.columns], on=['bwig',
                                                                 'chrom',
                                                                 'gene',
                                                                 'strand'])
    df_start = df_main.pop('start')
    df_end = df_main.pop('end')
    df_main.insert(2, 'start', df_start)
    df_main.insert(3, 'end', df_end)

    message("Writing to file")
    outputfile.close()

    with open(outputfile.name, 'a') as f:
        f.write(comments + "\n")
        df_main.to_csv(f, sep="\t", index=False,
                       mode='a',
                       columns=df_main.columns, na_rep='NA')

    # -------------------------------------------------------------------------
    #
    # Compress
    #
    # -------------------------------------------------------------------------

    message("Compressing")
    path = os.path.abspath(outputfile.name)
    filename = os.path.basename(path)
    message("filename: " + filename, type="DEBUG")
    zip_filename = filename + '.zip'
    message("zip_filename: " + zip_filename, type="DEBUG")
    zip_path = os.path.join(os.path.dirname(path), zip_filename)
    message("zip_path: " + zip_path, type="DEBUG")

    with zipfile.ZipFile(zip_path, 'w', allowZip64=True) as zf:
        zf.write(filename=path, arcname=filename)

    for i in outputfile_list:
        message("deleting " + outputfile_list[i].name)
        os.remove(outputfile_list[i].name)
    os.remove(outputfile.name)

    close_properly(inputfile, outputfile)


def main():
    """Main function."""
    my_parser = make_parser()
    args = my_parser.parse_args()
    args = dict(args.__dict__)
    mk_matrix(**args)


if __name__ == '__main__':
    main()


else:

    test = '''

        #mk_matrix: load dataset
        @test "mk_matrix_0" {
         result=`gtftk get_example -d simple -f "*"`
          [ "$result" = "" ]
        }
        
        #mk_matrix: test  unstranded
        @test "mk_matrix_1" {
         result=`gtftk mk_matrix -nst -p 0 -w 4 -u 2 -d 1 -c simple.chromInfo -i simple.gtf -y simple.bw -o simple_mat; unzip -u  simple_mat.zip &> /dev/null ; cat simple_mat | grep -v "#" | grep -v "main"| sort -k4,4n |cut -f7-10| perl -npe 's/\\t/|/g; s/\\n/,/'| sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "0.0|0.0|1.0|1.0,0.0|0.0|1.0|1.0,3.0|2.0|2.0|2.0,3.0|2.0|2.0|2.0,2.0|2.0|2.0|2.0,0.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,2.0|2.0|1.0|1.0,2.0|2.0|1.0|1.0,4.0|4.0|4.0|4.0,4.0|4.0|4.0|4.0,1.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,1.0|1.0|1.0|1.0," ]
        }
        
        
        #mk_matrix: test stranded
        @test "mk_matrix_2" {
         result=` rm -f simple_mat*; gtftk mk_matrix -p 0 -w 4 -u 2 -d 1 -c simple.chromInfo -i simple.gtf -y simple.bw -o simple_mat; unzip -u  simple_mat.zip &> /dev/null ;  cat simple_mat | grep -v "#" | grep -v "main"| sort -k4,4n | sort -k4,4n |cut -f7-10| perl -npe 's/\\t/|/g; s/\\n/,/'| sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "1.0|1.0|0.0|0.0,1.0|1.0|0.0|0.0,2.0|2.0|2.0|3.0,2.0|2.0|2.0|3.0,2.0|2.0|2.0|2.0,0.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,2.0|2.0|1.0|1.0,2.0|2.0|1.0|1.0,4.0|4.0|4.0|4.0,4.0|4.0|4.0|4.0,1.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,1.0|1.0|1.0|1.0," ]
        }
        
        #mk_matrix: test pseudo-count
        @test "mk_matrix_3" {
         result=`gtftk mk_matrix -p 1 -w 4 -u 2 -d 1 -c simple.chromInfo -i simple.gtf -y simple.bw -o simple_mat; unzip -u  simple_mat.zip &> /dev/null ;  cat simple_mat | grep -v "#" | grep -v "main" |  sort -k4,4n | sort -k4,4n |cut -f7-10| perl -npe 's/\\t/|/g; s/\\n/,/'| sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "2.0|2.0|1.0|1.0,2.0|2.0|1.0|1.0,3.0|3.0|3.0|4.0,3.0|3.0|3.0|4.0,3.0|3.0|3.0|3.0,1.0|1.0|1.0|1.0,1.0|1.0|1.0|1.0,1.0|1.0|1.0|1.0,3.0|3.0|2.0|2.0,3.0|3.0|2.0|2.0,5.0|5.0|5.0|5.0,5.0|5.0|5.0|5.0,2.0|1.0|1.0|1.0,1.0|1.0|1.0|1.0,2.0|2.0|2.0|2.0," ]
        }
        
        #mk_matrix: test transcript
        @test "mk_matrix_4" {
         result=`gtftk mk_matrix -p 1 -w 4 -d 0  -u 0 -t transcript -c simple.chromInfo -i simple.gtf -y simple.bw -o simple_mat; unzip -u  simple_mat.zip &> /dev/null ;  cat simple_mat | grep -v "#" | grep -v "main"|  sort -k4,4n | sort -k4,4n |cut -f7-10| perl -npe 's/\\t/|/g; s/\\n/,/' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "1.0|1.666667|2.0|2.0,1.0|1.666667|2.0|2.0,3.0|3.0|4.333333|4.0,3.5|3.0|2.5|3.0,3.666667|4.333333|3.0|3.666667,1.0|1.0|1.0|2.333333,1.0|1.0|1.0|1.0,1.0|1.0|1.0|1.0,2.0|2.0|2.0|2.0,2.0|2.0|2.0|2.0,4.666667|3.333333|2.666667|2.4,4.666667|3.333333|2.666667|2.4,1.0|1.0|1.0|1.4,1.0|1.0|1.5|2.75,2.5|3.0|3.0|2.0," ]
        }
        
        
        #mk_matrix: test transcript unstranded
        @test "mk_matrix_5" {
         result=`gtftk mk_matrix -nst -p 1 -V 1 -w 4 -K toto  -d 0  -u 0 -t transcript -c simple.chromInfo -i simple.gtf -y simple.bw -o simple_mat; unzip -u  simple_mat.zip &> /dev/null ;  cat simple_mat | grep -v "#" | grep -v "main"|  sort -k4,4n | sort -k4,4n |cut -f8,9,10,11| perl -npe 's/\\t/|/g; s/\\n/,/' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "2.0|1.666667|1.0,2.0|1.666667|1.0,4.333333|3.0|3.0,2.5|3.0|3.5,3.0|4.333333|3.666667,1.0|1.0|1.0,1.0|1.0|1.0,1.0|1.0|1.0,2.0|2.0|2.0,2.0|2.0|2.0,3.333333|2.666667|2.4,3.333333|2.666667|2.4,1.0|1.0|1.4,1.0|1.5|2.75,3.0|3.0|2.5," ]
        }
        
        #mk_matrix: every file contain 15 transcripts
        @test "mk_matrix_6" {
         result=`rm -Rf simple_mat*; rm -Rf mk_matrix_6; gtftk mk_matrix -i simple.gtf -u 2 -d 2 -t transcript -w 4 -y simple.bw -o simple_mat -c simple.chromInfo -K mk_matrix_6 -k 1; wc -l mk_matrix_6/* | grep 15 | wc -l`
          [ "$result" -eq 9 ]
        }
        
        
        #mk_matrix: test downstream + upstream (transcript)
        @test "mk_matrix_7" {
         result=`rm -Rf  simple_mat*; rm -Rf mk_matrix_7; gtftk mk_matrix -i simple.gtf -u 2 -d 2 -t transcript -w 4 -y simple.bw -o simple_mat -c simple.chromInfo ; unzip -u  simple_mat.zip &>/dev/null; cut -f8,9,10,11 simple_mat| tail -14| perl -npe 's/\\t/|/g; s/\\n/,/g' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "3.666667|2.333333|1.666667|1.4,0.0|0.0|0.5|1.75,0.0|0.0|0.0|1.333333,0.0|0.0|0.0|0.0,0.0|0.0|0.0|0.0,2.666667|3.333333|2.0|2.666667,2.0|2.0|3.333333|3.0,2.5|2.0|1.5|2.0,1.0|1.0|1.0|1.0,1.0|1.0|1.0|1.0,1.5|2.0|2.0|1.0,0.0|0.666667|1.0|1.0,0.0|0.666667|1.0|1.0,0.0|0.0|0.0|0.4," ]
        }
        
        #mk_matrix: test upstream.
        @test "mk_matrix_8" {
         result=`rm -Rf  simple_mat*; rm -Rf toto; gtftk mk_matrix -i simple.gtf -u 2 -d 2 -t transcript -w 4 -y simple.bw -o simple_mat -c simple.chromInfo ; unzip -u  simple_mat.zip &>/dev/null; cat simple_mat| cut -f7| grep -v "#" | grep -v "upstream" | perl -npe 's/\\t/|/g; s/\\n/,/g' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "4.0,4.0,0.0,0.0,0.0,0.0,2.0,2.0,2.0,2.0,2.0,1.0,1.0,1.0,0.5," ]
        }
        
        #mk_matrix: test downstream.
        @test "mk_matrix_9" {
         result=`rm -Rf  simple_mat*; rm -Rf toto; gtftk mk_matrix -i simple.gtf -u 2 -d 2 -t transcript -w 4 -y simple.bw -o simple_mat -c simple.chromInfo ; unzip -u  simple_mat.zip &>/dev/null; cat simple_mat| cut -f12| grep -v "#" | perl -npe 's/\\t/|/g; s/\\n/,/g' |  sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "downstream_1,2.0,2.0,2.0,2.0,0.0,0.0,1.0,2.5,3.0,3.0,3.0,1.5,1.0,1.0,2.0," ]
        }
        
        
        #mk_matrix: test coordinates.
        @test "mk_matrix_10" {
         result=`cut -f3 simple_mat | tail -14| sort |uniq| perl -npe 's/\\n/,/g' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "106,124,175,179,2,209,21,27,32,49,64," ]
        }
        
        #mk_matrix: test NA values
        @test "mk_matrix_11" {
         result=`cat simple.gtf | sed 's/chr1/chr2/' |gtftk mk_matrix -u 2 -d 2 -t transcript -w 4  -o simple_mat -c simple.chromInfo -V 1 -zn   -y simple.bw  ; unzip -u  simple_mat.zip &>/dev/null; cut -f7- simple_mat | tail -14| perl -npe 's/\\t/\\n/g'| sort | uniq`
          [ "$result" = "NA" ]
        }
        
        #mk_matrix: test header
        @test "mk_matrix_12" {
         result=`cat simple.gtf | gtftk mk_matrix -u 2 -d 2 -t transcript -w 4  -o simple_mat -c simple.chromInfo -V 1 -zn   -y simple.bw  ; unzip -u  simple_mat.zip &>/dev/null; cut -f7- simple_mat| head -2| tail -1| perl -npe 's/\\t/|/g'`
          [ "$result" = "upstream_1|main_1|main_2|main_3|main_4|downstream_1" ]
        }
        
        #mk_matrix: test --bin-frac
        @test "mk_matrix_13" {
         result=`rm -Rf simple_mat*; cat simple.gtf | sed 's/chr1/chr2/' |gtftk mk_matrix -u 5 -d 5 -b 0.5 -t transcript -w 10  -o simple_mat -c simple.chromInfo -V 1 -zn  -y simple.bw  ; unzip -u  simple_mat.zip &>/dev/null; cut -f7- simple_mat| head -2| tail -1| perl -npe 's/\\t/|/g'`
          [ "$result" = "upstream_1|upstream_2|upstream_3|upstream_4|upstream_5|main_1|main_2|main_3|main_4|main_5|main_6|main_7|main_8|main_9|main_10|downstream_1|downstream_2|downstream_3|downstream_4|downstream_5" ]
        }
        
        
        #mk_matrix: test orientation
        @test "mk_matrix_14" {
         result=`rm -Rf simple_mat*; cat simple.gtf  |gtftk mk_matrix -u 5 -d 5 -t transcript -w 5  -o simple_mat -c simple.chromInfo -V 1 -zn   -y simple.bw  ; unzip -u  simple_mat.zip &>/dev/null; cat simple_mat| grep G0005T001 simple_mat| cut -f8-12| perl -npe 's/\\t/|/g' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "2.0|3.333333|3.333333|2.0|2.666667" ]
        }
        
        #mk_matrix: test NO orientation (-nst)
        @test "mk_matrix_15" {
         result=`rm -Rf simple_mat*; cat simple.gtf  |gtftk mk_matrix -nst -u 5 -d 5 -t transcript -w 5  -o simple_mat -c simple.chromInfo -V 1 -zn   -y simple.bw  ; unzip -u  simple_mat.zip &>/dev/null; cat simple_mat| grep G0005T001 simple_mat| cut -f8-12| perl -npe 's/\\t/|/g' | sed 's/0000000005//g' | sed 's/69999999995/7/g'`
          [ "$result" = "2.666667|2.0|3.333333|3.333333|2.0" ]
        }
      
    '''

    CmdObject(name="mk_matrix",
              message="Compute a coverage matrix (see profile).",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="coverage",
              updated=__updated__,
              desc=__doc__,
              notes=__notes__,
              test=test)
