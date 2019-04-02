#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.arg_formatter import CheckChromFile
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-24"
__doc__ = """
 Find transcripts with divergent promoters. These transcripts will be defined here
 as those whose promoter region (defined by -u/-d) overlaps with the tss of
 another gene in reverse/antisens orientation. This may be useful to select
 coding genes in head-to-head orientation or LUAT as described in "Divergent
 transcription is associated with promoters of transcriptional regulators"
 (Lepoivre C, BMC Genomics, 2013). The output is a GTF with an additional key
 ('divergent') whose value is set to '.' if the gene has no antisens transcript
 in its promoter region. If the gene has an antisens transcript in its promoter
 region the 'divergent' key is set to the identifier of the transcript whose tss
 is the closest relative to the considered promoter. The tss to tss distance is
 also provided as an additional key (dist_to_divergent).
"""

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
                            help="Tabulated two-columns file. Chromosomes"
                                 " as column 1 and their sizes as"
                                 " column 2",
                            default=None,
                            metavar="CHROMINFO",
                            action=CheckChromFile,
                            required=True)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the promoter in 5' by a given value (int)."
                                 " Defines the region around the tss.",
                            default=1500,
                            metavar="UPSTREAM",
                            type=int,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the region in 3' by a given value (int)."
                                 " Defines the region around the tss.",
                            default=1500,
                            metavar="DOWNSTREAM",
                            type=int,
                            required=False)

    parser_grp.add_argument('-n', '--no-annotation',
                            help="Do not annotate the GTF. Just select the divergent transcripts.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-S', '--no-strandness',
                            help="Do not consider strandness (only look whether the promoter from a transcript overlaps with the promoter from another gene).",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-a',
                            '--key-name',
                            type=str,
                            default=None,
                            help="The name of the key.",
                            required=False)

    return parser


def divergent(
        inputfile=None,
        outputfile=None,
        key_name=None,
        upstream=1500,
        downstream=1500,
        chrom_info=None,
        no_strandness=False,
        no_annotation=False):
    """
Find transcript with divergent promoters.
    """

    message("Using -u " + str(upstream) + ".")
    message("Using -d " + str(downstream) + ".")

    tx_with_divergent = dict()
    dist_to_divergent = dict()
    tss_pos = dict()

    message("Loading GTF.")

    gtf = GTF(inputfile)

    message("Getting transcript coordinates.")

    tx_feat = gtf.select_by_key("feature",
                                "transcript")
    message("Getting tss coordinates.")

    tss_bo = tx_feat.get_tss(name=["transcript_id", "gene_id"],
                             sep="||")

    # get tss position
    for i in tss_bo:
        tx_id_tss, gn_id_tss = i.name.split("||")
        tss_pos[tx_id_tss] = int(i.start)

    message("Getting promoter coordinates.")

    promoter_bo = tss_bo.slop(s=True,
                              l=upstream,
                              r=downstream,
                              g=chrom_info.name).cut([0, 1,
                                                      2, 3,
                                                      4, 5])
    message("Intersecting...")

    if no_strandness:
        prom_with_tss_bo = promoter_bo.intersect(tss_bo,
                                                 wb=True,
                                                 s=False,
                                                 S=False)
    else:
        prom_with_tss_bo = promoter_bo.intersect(tss_bo,
                                                 wb=True,
                                                 s=False,
                                                 S=True)

    tmp_file = make_tmp_file("promoter_slop", ".bed")
    promoter_bo.saveas(tmp_file.name)
    tmp_file = make_tmp_file("promoter_intersection_with_tss_as_", ".bed")
    prom_with_tss_bo.saveas(tmp_file.name)

    for i in prom_with_tss_bo:

        tx_id_tss, gn_id_tss = i.fields[9].split("||")
        tx_id_prom, gene_id_prom = i.fields[3].split("||")

        if gene_id_prom != gn_id_tss:
            if tx_id_prom in tx_with_divergent:
                dist = abs(tss_pos[tx_id_prom] - tss_pos[tx_id_tss])
                if dist < dist_to_divergent[tx_id_prom]:
                    dist_to_divergent[tx_id_prom] = dist
                    tx_with_divergent[tx_id_prom] = tx_id_tss
            else:

                dist = abs(tss_pos[tx_id_prom] - tss_pos[tx_id_tss])
                dist_to_divergent[tx_id_prom] = dist
                tx_with_divergent[tx_id_prom] = tx_id_tss

    if not no_annotation:

        if key_name is None:
            key_name = "divergent"
            key_name_dist = "dist_to_divergent"
        else:
            key_name_dist = "dist_" + key_name

        if len(tx_with_divergent):
            gtf = gtf.add_attr_from_dict(feat="transcript",
                                         key="transcript_id",
                                         a_dict=tx_with_divergent,
                                         new_key=key_name)

            gtf = gtf.add_attr_from_dict(feat="transcript",
                                         key="transcript_id",
                                         a_dict=dist_to_divergent,
                                         new_key=key_name_dist)

        gtf.write(outputfile,
                  gc_off=True)

    else:
        gtf.select_by_key("transcript_id",
                          ",".join(list(tx_with_divergent.keys()))).write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    divergent(**args)


if __name__ == '__main__':
    main()

else:
    test = """

    #divergent: load dataset
    @test "divergent_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
        
    #divergent: the region as 2 divergent transcripts
    @test "divergent_1" {
     result=`gtftk divergent -i simple.gtf -c simple.chromInfo -u 4 -d 4 | gtftk select_by_key -k feature -v transcript | grep 'dist_to_divergent..4' | gtftk tabulate -H -k transcript_id| sort| perl -npe 's/\\n/,/'`
      [ "$result" = "G0003T001,G0004T001,G0004T002," ]
    }
    
    #divergent: the region as 2 divergent transcripts
    @test "divergent_2" {
     result=`gtftk divergent -i simple.gtf -c simple.chromInfo -u 4 -d 4 | grep 'dist_to_divergent..4' |gtftk select_by_key -k feature -v transcript| gtftk tabulate -H -k transcript_id| sort| perl -npe 's/\\n/,/'`
      [ "$result" = "G0003T001,G0004T001,G0004T002," ]
    }
    
    #divergent: the number of features is as expected.
    @test "divergent_3" {
     result=`gtftk divergent -i simple.gtf -c simple.chromInfo -u 4 -d 4 | wc -l`
      [ "$result" -eq 70 ]
    }
    
    #divergent: the number of exons is as expected.
    @test "divergent_4" {
     result=`gtftk divergent -i simple.gtf -c simple.chromInfo -u 4 -d 4 | awk '$3=="exon"'| wc -l`
      [ "$result" -eq 25 ]
    }
    
    #divergent: this region contains 4 divergent tx
    @test "divergent_5" {
     result=`gtftk divergent -u 18 -d 18 -c simple.chromInfo -i simple.gtf |  gtftk select_by_key -k feature -v  transcript| grep "dist_to_divergent \\"[0-9]"| gtftk tabulate -H -k transcript_id,dist_to_divergent,divergent| wc -l`
      [ "$result" -eq 4 ]
    }
    
    """
    CmdObject(name="divergent",
              message="Find transcripts with divergent promoters.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              updated=__updated__,
              group="annotation",
              notes=__notes__,
              test=test)
