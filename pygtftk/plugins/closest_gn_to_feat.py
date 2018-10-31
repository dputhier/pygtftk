#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import sys

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import message
from pygtftk.utils import write_properly
from pygtftk.arg_formatter import int_ge_to_null
from pygtftk.arg_formatter import bedFile,checkChromFile
from pybedtools import BedTool


__updated__ = "2018-01-20"
__doc__ = """
Get the gene/transcript closest to a set of feature. 
"""
__notes__ = """
 -- Features are provided in BED format.
 -- Output is in BED format.
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
                            help="Output file (BED).",
                            default=sys.stdout,
                            metavar="BED",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Bb][Ee][Dd]$',
                                                                     '\.[Bb][Ee][Dd]6$')))

    parser_grp.add_argument('-r', '--region-file',
                            help="The input BED file containing regions of interest.",
                            default=None,
                            metavar="BED",
                            required=False,
                            type=bedFile('r'))

    parser_grp.add_argument('-t', '--ft-type',
                            help="The feature of interest.",
                            default='gene',
                            choices=['gene', 'transcript'],
                            required=False)

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name in the output BED file.",
                            default="gene_id,gene_name",
                            metavar="NAME",
                            type=str)


    parser_grp.add_argument('-p', '--slop-value',
                            help="Look until this number of nucleotides in 5' and 3' direction",
                            default=0,
                            type=int_ge_to_null,
                            required=False)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-e', '--explicit',
                            help="Write explicitly the name of the keys in the header.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-c', '--chrom-info',
                            help='Tabulated file (chr as '
                                 'column 1, sizes as column 2.)',
                            default=None,
                            action=checkChromFile,
                            required=True)


    return parser



def closest_gn_to_feat(inputfile=None,
                     outputfile=None,
                     ft_type="gene",
                     names=None,
                     region_file=None,
                     slop_value=None,
                     chrom_info=None,
                     separator="|",
                     explicit=False,
                     tmp_dir=None,
                     logger_file=None,
                     verbosity=0):
    """
    Get the gene/transcript closest to a set of feature.
    """

    # -------------------------------------------------------------------------
    # Load the GTF an get the regions of interest
    # -------------------------------------------------------------------------

    nms = names.split(",")

    gtf = GTF(inputfile, check_ensembl_format=False)

    ann_bo = gtf.select_by_key("feature",
                                ft_type).to_bed(name=names)

    if not len(ann_bo):
        message("Requested feature (" + ft_type + ") could not be found. Use convert_ensembl maybe.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Load the BED file, slop and intersect with gene/tx
    # -------------------------------------------------------------------------


    region_bo = BedTool(region_file.name).slop(s=False,
                                            b=slop_value,
                                            g=chrom_info.name).cut([0, 1,
                                                                    2, 3,
                                                                    4, 5])


    for line in region_bo.intersect(ann_bo):
        print(line)

    close_properly(outputfile, inputfile)



def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    closest_gn_to_feat(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    

    #get_5p_3p_coord: test transpose
    @test "closest_gn_to_feat_11" {
     result=`gtftk  closest_gn_to_feat -p 10 -e -m bla -i simple.gtf -n transcript_id,gene_id,gene_name| head -1 | cut -f4`
      [ "$result" = "transcript_id=G0001T002|gene_id=G0001|gene_name=.|more_name=bla" ]
    }
    
    
    
    """
    CmdObject(name="closest_gn_to_feat",
              message="Get the gene/transcript closest to a set of feature.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              notes=__notes__,
              updated=__updated__,
              desc=__doc__,
              group="coordinates",
              test=test)
