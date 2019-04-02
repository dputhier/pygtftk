#!/usr/bin/env python


import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"

__doc__ = """
  Compute the locations of donor and acceptor splice sites.
"""

__notes__ = """
 This will return a single position, which corresponds to the most 5' and/or the most 3' intronic region.
 If the gtf file does not contain exon numbering you can compute it using the add_exon_nb command.
 The score column of the bed file contains the number of the closest exon relative to the splice site.
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
                            help="Output file.",
                            default=sys.stdout,
                            metavar="BED",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='bed')
                            )

    parser_grp.add_argument('-k', '--exon-numbering-key',
                            help="The name of the key containing the exon numbering (exon_number in ensembl)",
                            default="exon_number",
                            type=str)

    parser_grp.add_argument('-n', '--names',
                            help="The key(s) that should be used as name.",
                            default="exon_id,transcript_id,gene_id",
                            metavar="NAME",
                            type=str)

    parser_grp.add_argument('-s', '--separator',
                            help="The separator to be used for separating name elements (see -n).",
                            default="|",
                            metavar="SEP",
                            type=str)
    return parser


def splicing_site(inputfile=None,
                  outputfile=None,
                  exon_numbering_key=False,
                  names="exon_id,transcript_id,gene_id",
                  separator="\t"):
    """
    Compute the locations of splice donor are acceptor  sites. You may extend them in 3' and 5' depending on your needs.
    """

    gtf = GTF(inputfile)

    nb_exons = gtf.nb_exons()

    info = "feature,seqid,start,end,transcript_id," + exon_numbering_key
    info += ",strand," + names

    exon_info = gtf.extract_data_iter_list(info)

    for i in exon_info:

        if i[0] == "exon":
            if i[5] == ".":
                message("Some exon lines do not contain any numbering. "
                        "Use add_exon_nb or set --exon-numbering-key to the proper key.",
                        type="ERROR")

            if i[6] == "+":
                if int(i[5]) < nb_exons[i[4]]:
                    out = [i[1],
                           i[3],
                           str(int(i[3]) + 1),
                           separator.join(["donor"] + i[7:]),
                           i[5],
                           i[6]]
                    outputfile.write("\t".join(out) + "\n")

                if int(i[5]) > 1:
                    out = [i[1],
                           str(int(i[2]) - 2),
                           str(int(i[2]) - 1),
                           separator.join(["acceptor"] + i[7:]),
                           i[5],
                           i[6]]
                    outputfile.write("\t".join(out) + "\n")

            elif i[6] == "-":

                if int(i[5]) > 1:
                    out = [i[1],
                           i[3],
                           str(int(i[3]) + 1),
                           separator.join(["acceptor"] + i[7:]),
                           i[5],
                           i[6]]
                    outputfile.write("\t".join(out) + "\n")

                if int(i[5]) < nb_exons[i[4]]:
                    out = [i[1],
                           str(int(i[2]) - 2),
                           str(int(i[2]) - 1),
                           separator.join(["donor"] + i[7:]),
                           i[5],
                           i[6]]

                    outputfile.write("\t".join(out) + "\n")

    close_properly(outputfile, inputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    splicing_site(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #splicing_site: all
    @test "splicing_site_1" {
     result=`gtftk get_example | bedtools sort | gtftk add_exon_nb |  gtftk splicing_site -k exon_nbr| cut -f5| perl -npe 's/\\n/,/g'`
      [ "$result" = "3,2,2,2,1,2,1,1,2,1,1,1,2,2,2,2,3,3,2,1," ]
    }

    #splicing_site: all
    @test "splicing_site_2" {
     result=`gtftk get_example | bedtools sort | gtftk add_exon_nb |  gtftk splicing_site -k exon_nbr| cut -f2 | sort -n| perl -npe 's/\\n/,/g'`
      [ "$result" = "25,26,30,30,31,31,35,40,54,55,68,68,69,69,71,71,72,72,214,218," ]
    }

    #splicing_site: all
    @test "splicing_site_3" {
     result=`gtftk get_example | bedtools sort | gtftk add_exon_nb |  gtftk splicing_site -k exon_nbr| cut -f3 | sort -n| perl -npe 's/\\n/,/g'`
      [ "$result" = "26,27,31,31,32,32,36,41,55,56,69,69,70,70,72,72,73,73,215,219," ]
    }
   """

    CmdObject(name="splicing_site",
              message="Compute the locations of donor and acceptor splice sites.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              updated=__updated__,
              group="coordinates",
              desc=__doc__,
              test=test,
              notes=__notes__)
