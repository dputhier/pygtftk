#!/usr/bin/env python
"""
 Add the tss number to each transcript (5'->3').
"""

import argparse
import os
import sys
from collections import defaultdict

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly, make_tmp_file
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
- For each transcript feature add a key containing its tss number relative to the most 5'.
- See -\-add-nb-tss-to-gene to add the number of different tss to each gene feature.
"""


def make_parser():
    """The parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN.",
                            default=sys.stdin,
                            metavar="GTF",
                            required=False,
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-k', '--key-name',
                            help="The name of the new key.",
                            default='tss_number',
                            type=str,
                            required=False)

    parser_grp.add_argument('-c', '--compute-dist',
                            help="Add a key indicating the distance to the first tss (the most 5').",
                            action='store_true')

    parser_grp.add_argument('-d', '--key-name-dist',
                            help="If --compute-dist is true a name for that key.",
                            default='dist_to_first_tss',
                            type=str,
                            required=False)

    parser_grp.add_argument('-g', '--add-nb-tss-to-gene',
                            help="Add the number of different tss to each gene",
                            action='store_true')

    parser_grp.add_argument('-l', '--gene-key',
                            help="The name of the key if --add-nb-tss-to-gene is selected.",
                            default='nb_tss',
                            type=str,
                            required=False)

    return parser


def tss_numbering(
        inputfile=None,
        outputfile=None,
        compute_dist=False,
        key_name='tss_number',
        key_name_dist='dist_to_first_tss',
        add_nb_tss_to_gene=False,
        gene_key='nb_tss'):
    """
    Computes the distance between TSS of gene transcripts.
    """

    gtf = GTF(inputfile, check_ensembl_format=True)

    gn_tss_dist = defaultdict(dict)

    message("Getting TSSs.")
    tss = gtf.get_tss(name=["transcript_id"], as_dict=True)
    tx_to_gn = gtf.get_tx_to_gn()

    for k in tss:
        gn_id = tx_to_gn[k]
        gn_tss_dist[gn_id][k] = int(tss[k])

    # if_dict_of_dict is true, get_gn_to_tx() returns a dict of dict
    # that maps gene_id to transcript_id and transcript_id to TSS
    # numbering (1 for most 5', then 2...). For transcripts having
    # the same TSSs, the tss number will be the same.
    gn_to_tx_to_tss = gtf.get_gn_to_tx(as_dict_of_dict=True)

    message("Numbering TSSs.")

    tss_number_file = make_tmp_file(prefix='tx_to_tss_number', suffix='.txt')

    gn_how_many_tss = dict()

    for gn_id in gn_to_tx_to_tss:
        for tx_id in gn_to_tx_to_tss[gn_id]:
            tss_num = str(gn_to_tx_to_tss[gn_id][tx_id])
            tss_number_file.write(tx_id + "\t" + tss_num + "\n")
            if gn_id not in gn_how_many_tss:
                gn_how_many_tss[gn_id] = tss_num
            else:
                if int(tss_num) > int(gn_how_many_tss[gn_id]):
                    gn_how_many_tss[gn_id] = tss_num

    tss_number_file.close()

    gtf = gtf.add_attr_from_file(feat='transcript',
                                 key='transcript_id',
                                 new_key=key_name,
                                 inputfile=open(tss_number_file.name),
                                 has_header=False)

    if add_nb_tss_to_gene:

        gn_how_many_tss_file = make_tmp_file(prefix='gn_how_many_tss', suffix='.txt')

        for a_key, a_val in gn_how_many_tss.items():
            gn_how_many_tss_file.write(a_key + "\t" + a_val + "\n")

        gn_how_many_tss_file.close()

        gtf = gtf.add_attr_from_file(feat='gene',
                                     key='gene_id',
                                     new_key=gene_key,
                                     inputfile=open(gn_how_many_tss_file.name),
                                     has_header=False)

    if compute_dist:
        gn_to_tx_ordered_by_tss = gtf.get_gn_to_tx(ordered_5p=True)
        tss_dist_file = make_tmp_file(prefix='tx_tss_dist_to_first_tss', suffix='.txt')

        for gn_id in gn_to_tx_to_tss:
            tx_list = gn_to_tx_ordered_by_tss[gn_id]
            tx_first = tx_list.pop(0)
            # The first tss as distance 0 to the
            # first tss...
            tss_dist_file.write(tx_first + "\t0\n")
            for tx_id in tx_list:
                dist_to_first = abs(int(tss[tx_first]) - int(tss[tx_id]))
                tss_dist_file.write(tx_id + "\t" + str(dist_to_first) + "\n")

        tss_dist_file.close()

        gtf = gtf.add_attr_from_file(feat='transcript',
                                     key='transcript_id',
                                     new_key=key_name_dist,
                                     inputfile=open(tss_dist_file.name),
                                     has_header=False)

    gtf.write(outputfile,
              gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    tss_numbering(**args)


if __name__ == '__main__':
    main()

else:

    test = """
        
    #tss_numbering 
    @test "tss_numbering_0" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -k gene_id -v ENSG00000175756 | gtftk tss_numbering| gtftk select_by_key -t| gtftk tabulate -x -k transcript_id,start,end,strand,tss_number | md5 -r | perl -npe 's/\\s.*//'`
      [ "$result" = 'e01c4c1c167c8f2551584a3a6352b48b' ]
    }

    #tss_numbering 
    @test "tss_numbering_1" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -k gene_id -v ENSG00000142611 | gtftk tss_numbering| gtftk select_by_key -t| gtftk tabulate -x -k transcript_id,start,end,strand,tss_number | md5 -r | perl -npe 's/\\s.*//'`
      [ "$result" = 'c4a9b29bd0385b57ed64cde6325b46fe' ]
    }

    # G0011T001	250	267	-	2	13
    # G0011T002	250	280	-	1	0
    # tss_numbering 
    @test "tss_numbering_2" {
     result=`gtftk get_example -d simple_03 | gtftk tss_numbering -c | gtftk select_by_key -t | gtftk tabulate -x -k transcript_id,start,end,strand,tss_number,dist_to_first_tss | grep G0011 |md5 -r | perl -npe 's/\\s.*//'`
      [ "$result" = '7e950a2fbe9c832736f1f6def27ea6de' ]
    }


    # G0002T001	180	189	+	1	0
    # G0002T002	185	189	+	2	5
    # tss_numbering 
    @test "tss_numbering_2" {
     result=`gtftk get_example -d simple_03 | gtftk tss_numbering -c | gtftk select_by_key -t | gtftk tabulate -x -k transcript_id,start,end,strand,tss_number,dist_to_first_tss | grep G0002 |md5 -r | perl -npe 's/\\s.*//'`
      [ "$result" = '636ea1db16f5a44c3c6d83b46f372655' ]
    }

   
    
    """

    CMD = CmdObject(name="tss_numbering",
                    message=" Add the tss number to each transcript (5'->3').",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    notes=__notes__,
                    desc=__doc__,
                    group="annotation",
                    test=test)
