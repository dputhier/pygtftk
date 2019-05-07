# -*- coding: utf-8 -*-
"""
If several transcripts of a gene share the same TSS, select one transcript per TSS.
"""

import os

__updated__ = "2018-01-20"

import argparse
import operator
import sys
from collections import defaultdict

from pybedtools import BedTool

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

__notes__ = """
-- The alphanumeric order of transcript_id is used to select the representative of a TSS.
"""


def make_parser():
    """ Parser """
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Argument')

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

    return parser


def rm_dup_tss(inputfile=None,
               outputfile=None):
    """If several transcripts of a gene share the same tss, select only one."""

    # ----------------------------------------------------------------------
    # Get the TSS
    # ----------------------------------------------------------------------

    gtf = GTF(inputfile)
    tss_bo = gtf.get_tss(["gene_id", "transcript_id"])

    # ----------------------------------------------------------------------
    # Sort the file by name (4th col) to ensure reproducibility between calls.
    # ----------------------------------------------------------------------

    with open(tss_bo.fn) as f:
        lines = [line.split('\t') for line in f]

    tmp_file = make_tmp_file(prefix="tss_sorted_by_tx_id", suffix=".bed")

    for line in sorted(lines, key=operator.itemgetter(3)):
        tmp_file.write('\t'.join(line))

    tmp_file.close()

    tss_bo = BedTool(tmp_file.name)

    # ----------------------------------------------------------------------
    # Get the list of non redundant TSSs
    # ----------------------------------------------------------------------

    gene_dict = defaultdict(dict)
    to_delete = []

    message("Looking for redundant TSS (gene-wise).")

    for line in tss_bo:

        tss = line.start
        name = line.name
        gene_id, tx_id = name.split("|")

        if gene_id in gene_dict:
            if tss not in gene_dict[gene_id]:
                gene_dict[gene_id][tss] = tx_id
            else:
                to_delete += [tx_id]
        else:
            gene_dict[gene_id][tss] = tx_id

    message("Deleted transcripts: " + ",".join(to_delete[1:min(10,
                                                               len(to_delete)
                                                               )]) + "...",
            type="DEBUG")

    # ----------------------------------------------------------------------
    # Write
    # ----------------------------------------------------------------------

    gtf.select_by_key("feature",
                      "gene",
                      invert_match=True
                      ).select_by_key("transcript_id",
                                      ",".join(to_delete),
                                      invert_match=True).write(outputfile,
                                                               gc_off=True)


def main():
    """main"""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    rm_dup_tss(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    
    #Number of output line
    @test "rm_dup_tss_1" {
    result=$(gtftk get_example -d simple_05 |  gtftk rm_dup_tss | wc -l)
    [ $result -eq 45 ]
    }
    
    #Number of output line
    @test "rm_dup_tss_2" {
    result=$(gtftk get_example -d simple_05 | perl -MList::Util -e 'print List::Util::shuffle <>' | gtftk rm_dup_tss | wc -l)
    [ $result -eq 45 ]
    }
    
    
    #Check no duplicate TSS exists
    @test "rm_dup_tss_3" {
    result=$(gtftk get_example -d simple_05 |  gtftk rm_dup_tss  | perl -ne 'print if (/(G0001T002)|(G0003T002)|(G0004T002)|(G0006T002)|(G0007T002)|(G0008T002)/)')
    [ -z $result ]
    }

    #Check mini_real example.
    @test "rm_dup_tss_4" {
    result=$(gtftk get_example -d mini_real -f gtf | gtftk rm_dup_tss | gtftk select_by_key  -t |  gtftk get_5p_3p_coords -n gene_id | grep "\+" | cut -f2,4 | sort | uniq -c | sort -n| perl -npe 's/^ +//; s/ +/\\t/' | cut -f 1| sort | uniq)
    [ $result -eq 1 ]
    }

    #Check without rmdup (27 tx with the same TSS)
    @test "rm_dup_tss_5" {
    result=$(gtftk get_example -d mini_real  | gtftk select_by_key -t | gtftk get_5p_3p_coords -n gene_name | cut -f2,4 | awk 'BEGIN{n=0};{ if($2=="PCDH15" && $1=="54801290"){n++}}END{print n}')
    [ $result -eq 27 ]
    }
    
    
    #Check with rmdup (now 1 tx with the same TSS)
    @test "rm_dup_tss_6" {
    result=$(   gtftk get_example -d mini_real  | gtftk rm_dup_tss | gtftk select_by_key -t | gtftk get_5p_3p_coords -n gene_name | cut -f2,4 | awk 'BEGIN{n=0};{ if($2=="PCDH15" && $1=="54801290"){n++}}END{print n}')
    [ $result -eq 1 ]
    }
        
    """

    CmdObject(name="rm_dup_tss",
              message="If several transcripts of a gene share the same \
              TSS, select only one representative.",
              parser=make_parser(),
              group="selection",
              fun=os.path.abspath(__file__),
              notes=__notes__,
              updated=__updated__,
              desc=__doc__,
              test=test)
