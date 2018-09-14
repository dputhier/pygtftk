#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import shutil
import sys

from builtins import str
from builtins import zip

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import globbedFileList
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message, make_tmp_file

__updated__ = "2018-01-20"
__doc__ = """
 Get feature sequences in a flexible fasta format from a GTF file.
"""
__notes__ = """
 -- The sequences are returned in 5' to 3' orientation.
 -- If you want to use wildcards, use quotes :e.g. 'foo/bar*.fa'.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=FileWithExtension('r',
                                                   valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$'))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output FASTA file.",
                            default=sys.stdout,
                            metavar="FASTA",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Ff][Aa][Ss][Tt][Aa]$',
                                                                     '\.[Ff][Aa]$')))

    parser_grp.add_argument('-g', '--genome',
                            help="The genome in fasta format. Accept path with wildcards (e.g. *.fa).",
                            default=None,
                            action=globbedFileList,
                            required=True)

    parser_grp.add_argument('-s', '--separator',
                            help="To separate info in header.",
                            default="|",
                            type=str,
                            required=False)

    parser_grp.add_argument('-l', '--label',
                            help="A set of key for the header that will be extracted from the transcript line.",
                            default="feature,transcript_id,gene_id,seqid,start,end",
                            type=str,
                            required=False)

    parser_grp.add_argument('-t', '--feature-type',
                            help="The feature type (one defined in column 3).",
                            type=str,
                            default="exon",
                            required=False)

    parser_grp.add_argument('-n', '--no-rev-comp',
                            help="Don't reverse complement sequence "
                                 "corresponding to gene on minus strand.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-e', '--explicit',
                            help="Write explicitly the name of the keys in the header.",
                            action="store_true",
                            required=False)

    return parser


def get_feat_seq(inputfile=None,
                 outputfile=None,
                 genome=None,
                 feature_type="exon",
                 separator="",
                 no_rev_comp=False,
                 explicit=False,
                 label="",
                 tmp_dir=None,
                 logger_file=None,
                 verbosity=0):
    """
    Description: Get transcripts sequences in fasta format from a GTF file.
    """

    # -------------------------------------------------------------------------
    # Check args consistancy.
    #
    # -------------------------------------------------------------------------

    if feature_type in ["transcript", "gene"]:
        message("Please use get_tx_seq for transcript sequence.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Check chrom to avoid segfault
    # https://github.com/dputhier/libgtftk/issues/27
    # -------------------------------------------------------------------------

    genome_chr_list = []

    message("%d fasta files found." % len(genome))

    if len(genome) == 1:
        message("Checking fasta file chromosome list")
        genome = genome[0]
        with genome as genome_file:
            for i in genome_file:
                if i.startswith(">"):
                    i = i.rstrip("\n")
                    genome_chr_list += [i[1:]]
    else:
        message("Merging fasta files")
        tmp_genome = make_tmp_file(prefix="genome", suffix=".fa")
        with tmp_genome as tg:
            for curr_file in genome:
                message("Merging %s" % curr_file.name)
                with curr_file as cf:
                    shutil.copyfileobj(cf, tg, 1024 * 1024 * 100)

        message("Checking fasta file chromosome list")
        genome = open(tmp_genome.name, "r")
        with genome as genome_file:
            for i in genome_file:
                if i.startswith(">"):
                    i = i.rstrip("\n")
                    genome_chr_list += [i[1:]]

    rev_comp = not no_rev_comp

    gtf = GTF(inputfile)

    gtf_chr_list = gtf.get_chroms(nr=True)

    # Check chrom to avoid segfault
    # https://github.com/dputhier/libgtftk/issues/27
    message("Comparing chromosomes from GTF and Fasta files.")
    gtf_chr_list_found = [x for x in gtf_chr_list if x in genome_chr_list]

    if len(gtf_chr_list_found) == 0:
        message("Chromosome from GTF were not found in fasta file",
                type="ERROR")

    if len(gtf_chr_list_found) != len(gtf_chr_list):
        not_found = [x for x in gtf_chr_list if x not in gtf_chr_list_found]
        message("Some chromosomes were not found in the fasta file: %s" % ",".join(not_found),
                type="ERROR")

    # -------------------------------------------------------------------------
    # Retrieving fasta sequences
    #
    # -------------------------------------------------------------------------

    message("Retrieving fasta sequences.")
    fasta_seq = gtf.get_sequences(genome=genome.name,
                                  intron=True,
                                  rev_comp=rev_comp)

    tx_gtf = gtf.select_by_key("feature", "transcript")

    label = ",".join([x for x in label.split(',') if x != 'feature'])

    if label != '':
        tx_info = tx_gtf.extract_data("transcript_id," +
                                      label,
                                      as_dict_of_lists=True)

    else:
        tx_info = dict()

    for i in fasta_seq.iter_features(feat=feature_type):

        if not explicit:
            feat_info = [feature_type] + [str(i.start)] + [str(i.end)]
            header = separator.join(tx_info[i.transcript_id] + feat_info)
        else:
            feat_info = ["feature_type=" + feature_type] + \
                        ["feature_start=" + str(i.start)] + \
                        ["feature_end=" + str(i.end)]
            header = [str(x[0]) + "=" + x[1]
                      for x in zip(label.split(","), tx_info[i.transcript_id])]
            header = separator.join(header + feat_info)
        outputfile.write(">" + header + "\n")
        outputfile.write(i.sequence + "\n")

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_feat_seq(**args)


if __name__ == '__main__':
    main()


else:

    test = """

    #get_feat_seq:
    @test "get_feat_seq_0.0" {
     result=`echo 1`
      [ "$result" = "1" ]
    }
    
    #get_feat_seq:
    @test "get_feat_seq_0.1" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  exon  | grep "G0003T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "aatta,gcttg," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_1" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  exon -n | grep "G0003T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "caagc,taatt," ]
    }
    
    #get_feat_seq:
    @test "get_feat_seq_2" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  exon  | grep "G0003T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "aatta,gcttg," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_3" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  CDS  | grep "G0004T002" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "tct,g,gc," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_4" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  CDS  | grep "G0006T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "ct,att,acat," ]
    }
    
    #get_feat_seq:
    @test "get_feat_seq_5" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  CDS -n | grep "G0006T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "atgt,aat,ag," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_6" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa  -l feature,transcript_id,start -t  CDS -n | grep "G0006T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "atgt,aat,ag," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_7" {
     result=`gtftk get_feat_seq -i pygtftk/data/simple/simple.gtf -g pygtftk/data/simple/simple.fa   -t  CDS  | grep G0006T001 | perl -npe 's/\\n/,/g'`
      [ "$result" = ">G0006T001|G0006|chr1|22|35|CDS|33|34,>G0006T001|G0006|chr1|22|35|CDS|28|30,>G0006T001|G0006|chr1|22|35|CDS|22|25," ]
    }
    
        
    """

    CmdObject(name="get_feat_seq",
              message="Get feature sequence (e.g exon, UTR...).",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              group="sequences",
              notes=__notes__,
              test=test)
