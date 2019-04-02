#!/usr/bin/env python

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import chomp
from pygtftk.utils import close_properly
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message
from pygtftk.utils import write_properly

__updated__ = "2018-01-20"
__doc__ = """
 Get the size and limits (start/end) of features enclosed in the GTF.
 The feature can be of any type (as found in the 3rd column of the GTF)
 or 'mature_rna' to get transcript size (i.e without introns).
 If bed format is requested returns the limits in bed format and the size as a score.
 Otherwise output GTF file with 'feat_size' as a new key and size as value.
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
                            help="Output file (BED).",
                            default=sys.stdout,
                            metavar="GTF/BED",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('bed', 'gtf')))

    parser_grp.add_argument('-t', '--ft-type',
                            help="A target feature (as found in the 3rd "
                                 "column of the GTF) or 'mature_rna' to get transcript size (without introns).",
                            default='transcript',
                            type=str,
                            required=False)

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name if bed is requested.",
                            default="transcript_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-k', '--key-name',
                            help="Name for the new key if GTF format is requested.",
                            default='feat_size',
                            metavar="KEY_NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)

    parser_grp.add_argument('-b', '--bed',
                            help="Output bed format.",
                            action="store_true",
                            required=False)

    return parser


def feature_size(
        inputfile=None,
        outputfile=None,
        ft_type="transcript",
        names="transcript_id",
        key_name='feature_size',
        separator="|",
        bed=False):
    """
 Get the size and limits (start/end) of features enclosed in the GTF. If bed
 format is requested returns the limits zero-based half open and the size as a score.
 Otherwise output GTF file with 'feat_size' as a new key and size as value.
    """

    message("Computing feature sizes.")

    gtf = GTF(inputfile)

    feat_list = gtf.get_feature_list(nr=True) + ['mature_rna']

    if ft_type not in feat_list + ["*"]:
        message("Unable to find requested feature.", type="ERROR")

    names = names.split(",")

    if ft_type != 'mature_rna':

        if bed:
            bed_obj = gtf.select_by_key("feature", ft_type).to_bed(name=names,
                                                                   sep=separator,
                                                                   add_feature_type=True)

            for i in bed_obj:
                i.score = str(i.end - i.start)
                write_properly(chomp(str(i)), outputfile)
        else:

            tmp_file = make_tmp_file(prefix="feature_size", suffix=".txt")

            elmt = gtf.extract_data("feature,start,end", as_list_of_list=True, no_na=False, hide_undef=False)

            for i in elmt:
                if i[0] != ft_type and ft_type != "*":
                    tmp_file.write("?\n")
                else:
                    tmp_file.write(str(int(i[2]) - int(i[1]) + 1) + "\n")

            tmp_file.close()

            gtf.add_attr_column(tmp_file, key_name).write(outputfile, gc_off=True)



    else:

        tx_size = gtf.get_transcript_size()

        if bed:
            bed_obj = gtf.select_by_key("feature",
                                        'transcript').to_bed(['transcript_id'] + names,
                                                             add_feature_type=False,
                                                             sep=separator,
                                                             more_name=['mature_rna'])

            for i in bed_obj:
                names = i.name.split(separator)
                tx_id = names.pop(0)
                i.score = tx_size[tx_id]
                i.name = separator.join(names)
                write_properly(chomp(str(i)), outputfile)
        else:

            if len(tx_size):
                gtf = gtf.add_attr_from_dict(feat="transcript",
                                             key="transcript_id",
                                             a_dict=tx_size,
                                             new_key=key_name)

            gtf.write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    feature_size(**args)


if __name__ == '__main__':
    main()

else:

    test = """

    # feature_size: load dataset
    @test "feature_size_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }
    
     
    #feat_size: check line number
    @test "feature_size_1" {
     result=`gtftk feature_size -i simple.gtf| wc  -l`
      [ "$result" -eq 70 ]
    }
    
    #feat_size: check number of columns
    @test "feature_size_2" {
     result=`gtftk feature_size -i simple.gtf| awk 'BEGIN{FS="\\t"}{print NF}'| sort| uniq`
      [ "$result" -eq 9 ]
    }
    
    #feat_size: check the right elements are flagged
    @test "feature_size_3" {
     result=`gtftk feature_size -i simple.gtf | grep feat_size| wc -l`
      [ "$result" -eq 15 ]
    }
    
    #feat_size: check the right elements are flagged
    @test "feature_size_4" {
     result=`gtftk feature_size -i simple.gtf -t gene | grep feat_size | wc -l`
      [ "$result" -eq 10 ]
    }
    
    #feat_size: check size is ok.
    @test "feature_size_5" {
     result=`gtftk feature_size -i simple.gtf -t gene| grep feat_size| grep G0001 | sed 's/.*feat_size "//'| sed 's/".*//'`
      [ "$result" -eq 14 ]
    }
    
    #feat_size: check line number with -b.
    @test "feature_size_6" {
     result=`gtftk feature_size -i simple.gtf -t gene -b| wc -l`
      [ "$result" -eq 10 ]
    }
    
    #feat_size: check number of columns
    @test "feature_size_7" {
     result=`gtftk feature_size -i simple.gtf -t gene -b| awk 'BEGIN{FS="\\t"}{print NF}'| sort | uniq`
      [ "$result" -eq 6 ]
    }
    
    
    #feat_size: check size
    @test "feature_size_8" {
     result=`gtftk feature_size -i simple.gtf -t gene -b | cut -f 5| sort -n | perl -npe 's/\\n/,/'`
      [ "$result" = "10,10,11,12,12,12,13,14,14,15," ]
    }
    
    #feat_size: check size of mature rna (bed)
    @test "feature_size_9" {
     result=`gtftk feature_size -i simple.gtf  -t mature_rna -b -n transcript_id,gene_id |cut -f 5| sort | uniq | perl -npe 's/\\n/,/'`
      [ "$result" = "10,11,12,14,6,8,9," ]
    }
    
    #feat_size: check size of mature rna (gtf)
    @test "feature_size_10" {
     result=`gtftk feature_size -i simple.gtf  -t mature_rna  -n transcript_id,gene_id | gtftk select_by_key -k feature -v transcript| gtftk tabulate -k feat_size -H | sort | uniq  | perl -npe 's/\\n/,/'`
      [ "$result" = "10,11,12,14,6,8,9," ]
    }

    #feat_size: check size of mature rna (gtf)
    @test "feature_size_11" {
     result=`gtftk feature_size -i simple.gtf  -t mature_rna  -n transcript_id,gene_id | wc -l`
      [ "$result" -eq 70 ]
    }
    
    
    """
    CMD = CmdObject(name="feature_size",
                    message="Compute the size of features enclosed in the GTF.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    updated=__updated__,
                    desc=__doc__,
                    group="information",
                    test=test)
