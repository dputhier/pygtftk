#!/usr/bin/env python

import argparse
import os
import re
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
 Get feature sequences (i.e. column 3) in a flexible fasta format from a GTF file. 
"""
__notes__ = """
 -- The sequences are returned in 5' to 3' orientation.
 -- If you want to use wildcards, use quotes: e.g. 'foo/bar*.fa'.
 -- See get_tx_seq for mature RNA sequence.
 -- If -\-unique is used if a header was already encountered the record won't be print. 
 Take care to use unambiguous identifiers in the header.
"""


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
                            help="Output FASTA file.",
                            default=sys.stdout,
                            metavar="FASTA",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='fasta'))

    parser_grp.add_argument('-g', '--genome',
                            help="The genome in fasta format.",
                            default=None,
                            metavar="FASTA",
                            required=True,
                            type=arg_formatter.FormattedFile(mode='r', file_ext='fasta'))

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

    parser_grp.add_argument('-r', '--rev-comp-to-header',
                            help="Indicate in the header whether sequence was rev-complemented.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-u', '--unique',
                            help="Don't write redondant IDS.",
                            action="store_true",
                            required=False)

    return parser


def get_feat_seq(inputfile=None,
                 outputfile=None,
                 genome=None,
                 feature_type="exon",
                 separator="",
                 no_rev_comp=False,
                 label="",
                 rev_comp_to_header=False,
                 unique=False):
    """
    Description: Get transcripts sequences in fasta format from a GTF file.
    """

    # -------------------------------------------------------------------------
    # Should sequences be reverse-complemented
    # -------------------------------------------------------------------------

    force_strandedness = not no_rev_comp

    # -------------------------------------------------------------------------
    # Check chrom to avoid segfault
    # https://github.com/dputhier/libgtftk/issues/27
    # -------------------------------------------------------------------------

    if genome.name.endswith(".gz"):
        message("Genome in gz format is not currently supported.", type="ERROR")

    genome_chr_list = []

    message("Fasta files found: %s" % genome.name)

    message("Checking fasta file chromosome list")

    with genome as geno:
        for i in geno:
            if i.startswith(">"):
                i = i.rstrip("\n")
                genome_chr_list += [i[1:]]

    gtf = GTF(inputfile, check_ensembl_format=False)

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

    feat_seq = gtf.select_by_key("feature",
                                 feature_type
                                 ).to_bed(name=label.split(","),
                                          sep=separator
                                          ).sequence(fi=genome.name,
                                                     name=True,
                                                     s=force_strandedness)

    id_printed = set()

    to_print = True

    for nb_line, line in enumerate(open(feat_seq.seqfn)):
        if line.startswith(">"):

            # This (+/-) may be added by pybedtool
            # but can be accessed though --label
            line = re.sub("\(\+\)$", "", line)
            line = re.sub("\(\-\)$", "", line)

            if rev_comp_to_header:
                if force_strandedness:
                    line = line + separator + "rev_comp"
                else:
                    line = line + separator + "no_rev_comp"

            if unique:
                if line in id_printed:
                    to_print = False
            if to_print:
                outputfile.write(line)
                id_printed.add(line)



        else:
            if not to_print:
                to_print = True
            else:
                outputfile.write(line)

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

    #get_feat_seq: load dataset
    @test "get_feat_seq_0" {
     result=`gtftk get_example -f '*' -d simple; gtftk get_example -f '*' -d mini_real_10M; if [ ! -f chr1_hg38_10M.fa ]; then gunzip -f chr1_hg38_10M.fa.gz; fi `
      [ "$result" = "" ]
    }
       
    #get_feat_seq:
    @test "get_feat_seq_1" {
     result=`gtftk get_feat_seq -i simple.gtf -g simple.fa  -l feature,transcript_id,start -t  exon  | grep "G0003T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "gcttg,aatta," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_2" {
     result=`gtftk get_feat_seq -i simple.gtf -g simple.fa  -l feature,transcript_id,start -t  exon -n | grep "G0003T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "caagc,taatt," ]
    }
    


    #get_feat_seq:
    @test "get_feat_seq_3" {
     result=`gtftk get_feat_seq -i simple.gtf -g simple.fa  -l feature,transcript_id,start -t  CDS  | grep "G0004T002" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "tct,g,gc," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_4" {
     result=`gtftk get_feat_seq -i simple.gtf -g simple.fa  -l feature,transcript_id,start -t  CDS  | grep "G0006T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "acat,att,ct," ]
    }
    
    #get_feat_seq:
    @test "get_feat_seq_5" {
     result=`gtftk get_feat_seq -i simple.gtf -g simple.fa  -l feature,transcript_id,start -t  CDS -n | grep "G0006T001" -A 1| perl -ne  'chomp, print $_,"," if(/^[AaTtCcGg]+$/)'`
      [ "$result" = "atgt,aat,ag," ]
    }

    #get_feat_seq:
    @test "get_feat_seq_7" {
     result=`gtftk get_feat_seq -i simple.gtf -g simple.fa   -t  CDS  | grep G0006T001 | perl -npe 's/\\n/,/g'`
      [ "$result" = ">CDS|G0006T001|G0006|chr1|21|25,>CDS|G0006T001|G0006|chr1|27|30,>CDS|G0006T001|G0006|chr1|32|34," ]
    }
 
    #check with a real dataset (no rev-comp) minus strand
    @test "get_feat_seq_8" {
     result=`gtftk select_by_key -f ids_minus_exon.txt -k exon_id -i mini_real_10M.gtf.gz | gtftk get_feat_seq -g chr1_hg38_10M.fa -l exon_id -u | perl -ne 'print uc $_' > observed_sequence_minus_exon.fa`
      [ -f  observed_sequence_minus_exon.fa ]
    }

    #get_feat_seq:
    @test "get_feat_seq_9" {
     result=`diff observed_sequence_minus_exon.fa expected_sequence_minus_exon.fa`
      [ "$result" = "" ]
    }    
    
     #check with a real dataset (no rev-comp) plus strand
    @test "get_feat_seq_10" {
     result=`gtftk select_by_key -f ids_plus_exon.txt -k exon_id -i mini_real_10M.gtf.gz | gtftk get_feat_seq -g chr1_hg38_10M.fa -l exon_id -u | perl -ne 'print uc $_' > observed_sequence_plus_exon.fa`
      [ -f  observed_sequence_plus_exon.fa ]
    }

    #get_feat_seq:
    @test "get_feat_seq_11" {
     result=`diff observed_sequence_plus_exon.fa expected_sequence_plus_exon.fa`
      [ "$result" = "" ]
    }       

    #check with a real dataset (no rev-comp) minus strand
    @test "get_feat_seq_12" {
     result=`gtftk select_by_key -f ids_minus_exon.txt -k exon_id -i mini_real_10M.gtf.gz | gtftk get_feat_seq -g chr1_hg38_10M.fa -l exon_id -u -n | perl -ne 'print uc $_' > observed_sequence_minus_exon_no_rv.fa`
      [ -f  observed_sequence_minus_exon_no_rv.fa ]
    }

    #get_feat_seq:
    @test "get_feat_seq_13" {
     result=`diff observed_sequence_minus_exon_no_rv.fa expected_sequence_minus_exon_no_rv.fa`
      [ "$result" = "" ]
    }    

     #check with a real dataset (with rev-comp) plus strand
    @test "get_feat_seq_14" {
     result=`gtftk select_by_key -f ids_plus_exon.txt -k exon_id -i mini_real_10M.gtf.gz | gtftk get_feat_seq -g chr1_hg38_10M.fa -l exon_id -u -n | perl -ne 'print uc $_' > observed_sequence_plus_exon.fa`
      [ -f  observed_sequence_plus_exon.fa ]
    }

    #get_feat_seq:
    @test "get_feat_seq_15" {
     result=`diff observed_sequence_plus_exon.fa expected_sequence_plus_exon.fa`
      [ "$result" = "" ]
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
