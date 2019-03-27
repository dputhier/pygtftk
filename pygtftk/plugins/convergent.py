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

__updated__ = "2018-01-20"

__doc__ = """
 Find transcripts with convergent tts. These transcripts will be defined here
 as those whose tts region (defined by -u/-d) overlaps with the tts of
 another gene in reverse/antisens orientation. The output is a GTF with an
 additional key ('convergent') whose value is set to '.' if the gene has no
 convergent transcript in its tts region. If the gene has an antisens transcript
 in its tts region the 'convergent' key is set to the identifier of the
 transcript whose tts is the closest relative to the considered tts.
 The tts to tts distance is also provided as an additional key (dist_to_convergent).
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
                                 " as column 1 and sizes as"
                                 " column 2 ",
                            default=None,
                            metavar="CHROMINFO",
                            action=CheckChromFile,
                            required=True)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extends the tts in 5' by a given value (int)."
                                 " Defines the region around the tts.",
                            default=1500,
                            metavar="UPSTREAM",
                            type=int,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extends the region in 3' by a given value (int)."
                                 " Defines the region around the tts.",
                            default=1500,
                            metavar="DOWNSTREAM",
                            type=int,
                            required=False)

    return parser


def convergent(
        inputfile=None,
        outputfile=None,
        upstream=1500,
        downstream=1500,
        chrom_info=None):
    """
    Find transcript with convergent tts.
    """

    message("Using -u " + str(upstream) + ".")
    message("Using -d " + str(downstream) + ".")

    tx_to_convergent_nm = dict()
    dist_to_convergent = dict()
    tts_pos = dict()

    message("Loading GTF.")

    gtf = GTF(inputfile)

    message("Getting transcript coordinates.")

    tx_feat = gtf.select_by_key("feature", "transcript")

    message("Getting tts coordinates.")

    tts_bo = tx_feat.get_tts(name=["transcript_id", "gene_id"],
                             sep="||")

    # get tts position
    for i in tts_bo:
        tx_id_ov, gn_id_ov = i.name.split("||")
        tts_pos[tx_id_ov] = int(i.start)

    message("Getting tts coordinates.")

    tts_region_bo = tts_bo.slop(s=True,
                                l=upstream,
                                r=downstream,
                                g=chrom_info.name).cut([0, 1,
                                                        2, 3,
                                                        4, 5])

    message("Intersecting...")
    tts_intersect_bo = tts_region_bo.intersect(tts_bo,
                                               wb=True,
                                               s=False,
                                               S=True)

    tmp_file = make_tmp_file("tts_slop", ".bed")
    tts_region_bo.saveas(tmp_file.name)
    tmp_file = make_tmp_file("tts_slop_intersection_with_tts_as_", ".bed")
    tts_intersect_bo.saveas(tmp_file.name)

    for i in tts_intersect_bo:

        tx_id_main, gene_id_main = i.fields[3].split("||")
        tx_id_ov, gn_id_ov = i.fields[9].split("||")

        if gene_id_main != gn_id_ov:
            if tx_id_main in tx_to_convergent_nm:
                dist = abs(tts_pos[tx_id_main] - tts_pos[tx_id_ov])
                if dist < dist_to_convergent[tx_id_main]:
                    dist_to_convergent[tx_id_main] = dist
                    tx_to_convergent_nm[tx_id_main] = tx_id_ov
            else:
                dist = abs(tts_pos[tx_id_main] - tts_pos[tx_id_ov])
                dist_to_convergent[tx_id_main] = dist
                tx_to_convergent_nm[tx_id_main] = tx_id_ov

    if len(tx_to_convergent_nm):
        gtf = gtf.add_attr_from_dict(feat="transcript",
                                     key="transcript_id",
                                     a_dict=tx_to_convergent_nm,
                                     new_key="convergent")

        gtf = gtf.add_attr_from_dict(feat="transcript",
                                     key="transcript_id",
                                     a_dict=dist_to_convergent,
                                     new_key="dist_to_convergent")

    gtf.write(outputfile,
              gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    convergent(**args)


if __name__ == '__main__':
    main()


else:

    test = """
    
    #convergent: load dataset
    @test "convergent_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }

    #convergent: this region contains 3 convergent tx
    @test "convergent_1" {
     result=`gtftk convergent -K toto -i simple.gtf  -c simple.chromInfo -u 24 -d 24| gtftk select_by_key -t |gtftk tabulate -H -k transcript_id,dist_to_convergent -s ","| grep ",[0-9]"| wc -l`
     [ "$result" -eq 3 ]
    }

    #convergent: 70 lines in output
    @test "convergent_2" {
     result=`gtftk convergent -K toto -i simple.gtf  -c simple.chromInfo -u 12 -d 12| wc -l`
      [ "$result" -eq 70 ]
    }
    
    #convergent: 25 exons
    @test "convergent_3" {
     result=`gtftk convergent -K toto -i simple.gtf  -c simple.chromInfo -u 12 -d 12| gtftk select_by_key -k feature -v exon| wc -l`
      [ "$result" -eq 25 ]
    }
    """

    CmdObject(name="convergent",
              message="Find transcripts with convergent tts.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              desc=__doc__,
              group="annotation",
              notes=__notes__,
              test=test)
