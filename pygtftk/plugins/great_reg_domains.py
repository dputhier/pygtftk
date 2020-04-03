# -*- coding: utf-8 -*-
""" Given a GTF and a GO term computes labeled regions using GREAT 'association rule'. """

import argparse
import os
import sys

from pygtftk import arg_formatter
from pygtftk.biomart import Biomart
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import message, chrom_info_as_dict, make_tmp_file

__updated__ = "2018-01-25"

__notes__ = """

-- This tool represents an attempt to process genomic annotations in GTF format in order to produced 
a set of 'labeled' regions with the same rules as those described in GREAT (Genomic Regions Enrichment of Annotations
Tool) documentation. We can not warrant that the procedure is exactly the same.
-- The tool only currently supports 'basal_plus_extension' association rule.
-- Operations are performed on a transcript TSS basis rather than a gene TSS basis. To our knowledge, the way 
representative TSS are selected for each gene is not described in GREAT paper nor documentation. As a consequence 
the output contains several regulatory domain for a given gene (one per TSS).
-- The tool does not include 'curated regulatory domains' as proposed by GREAT.
-- The tool does not allow background regions to be imported.
-- The resulting BED file can be used for instance as OLOGRAM input (using -z/--no-gtf and -b/--more-bed) to check the GREAT
results by assessing whether the enrichment is also significant regarding the number of overlapping nucleotides. 
"""


def make_parser():
    """parse"""
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
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('bed')))

    parser_grp.add_argument('-g', '--go-id',
                            help='The GO ID (e.g GO:0003700). Default is to return all genes.',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-s', '--species',
                            help='The dataset/species. Use select_by_go to get the list of species',
                            default="hsapiens",
                            type=str,
                            required=False)

    parser_grp.add_argument('-u', '--upstream',
                            help="Extend the TSS in 5' by a given value.",
                            default=5000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-d', '--downstream',
                            help="Extend the TSS 3' by a given value. ",
                            default=1000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-t', '--distal',
                            help="Maximum extension in one direction",
                            default=1000000,
                            type=int,
                            required=False)

    parser_grp.add_argument('-c', '--chrom-info',
                            help="Tabulated two-columns file. "
                                 "Chromosomes as column 1, sizes as column 2",
                            default=None,
                            metavar="TXT",
                            action=arg_formatter.CheckChromFile,
                            required=True)

    parser_grp.add_argument('-m', '--mode',
                            help="The type of 'association rule'. ",
                            default="basal_plus_extension",
                            choices=["basal_plus_extension"],
                            metavar="TXT",
                            action=arg_formatter.CheckChromFile,
                            required=False)

    parser_grp.add_argument('-p1', '--http-proxy',
                            help='Use this http proxy (not tested/experimental).',
                            default='',
                            type=str,
                            required=False)

    parser_grp.add_argument('-p2', '--https-proxy',
                            help='Use this https proxy (not tested/experimental).',
                            default='',
                            type=str,
                            required=False)

    return parser


# -------------------------------------------------------------------------
# The XML query
# -------------------------------------------------------------------------

XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
            
    <Dataset name = "{species}_gene_ensembl" interface = "default" >
        <Filter name = "go_parent_term" value = "{go}"/>
        <Attribute name = "ensembl_gene_id" />
    </Dataset>
</Query>"""


# -------------------------------------------------------------------------
# Main function definition
# -------------------------------------------------------------------------


def great_reg_domains(inputfile=None,
                      outputfile=None,
                      go_id="GO:0003700",
                      species="hsapiens",
                      upstream=1000,
                      downstream=1000,
                      chrom_info=None,
                      distal=1000000,
                      mode='basal_plus_extension',
                      http_proxy=None,
                      https_proxy=None):
    """ Given a GTF and a GO term, attempt compute labeled regions using GREAT 'association rule'. """

    # -------------------------------------------------------------------------
    # chrom_len will store the chromosome sizes.
    # -------------------------------------------------------------------------

    chrom_len = chrom_info_as_dict(chrom_info)

    # -------------------------------------------------------------------------
    # Read the GTF
    # -------------------------------------------------------------------------

    gtf = GTF(inputfile, check_ensembl_format=False)

    # -------------------------------------------------------------------------
    # Get the TSSs -- Extend them by upstream/dowstream
    # -------------------------------------------------------------------------

    message("Defining basal regulatory domains.", type="INFO")
    basal_reg_bed = gtf.get_tss(name=['gene_id', 'gene_name']).slop(s=True,
                                                                    l=upstream,
                                                                    r=downstream,
                                                                    g=chrom_info.name).sort()

    basal_reg_bed_file = make_tmp_file(prefix='basal_reg', suffix='.bed')
    basal_reg_bed.saveas(basal_reg_bed_file.name)

    if mode == 'basal_plus_extension':
        # -------------------------------------------------------------------------
        # Search for upstream limits of each basal_reg_bed
        # Here we ignore overlapping  basal_reg_bed as the way they
        # are proceded is not documented in GREAT to our knowledge
        # -------------------------------------------------------------------------
        message("Defining regulatory domains upstream regions.", type="INFO")

        regulatory_region_start = dict()
        regulatory_region_end = dict()
        chroms = dict()
        strands = dict()

        basal_reg_bed_upstream = basal_reg_bed.closest(basal_reg_bed,
                                                       # Ignore features in B that overlap A
                                                       io=True,
                                                       # In addition to the closest feature in B report distance
                                                       # use negative distances to report upstream features.
                                                       # Report distance with respect to A.
                                                       # When A is on the - strand, "upstream" means B has a
                                                       # higher(start, stop).
                                                       D="a",
                                                       # Ignore features in B that are downstream of features in A
                                                       id=True,
                                                       # How ties are handled. "first"  Report the first tie
                                                       t="first",
                                                       # Require that the query and the closest hit have different names/gene_ids.
                                                       N=True)

        basal_reg_bed_upstream_file = make_tmp_file(prefix='basal_reg_bed_upstream', suffix='.bed')
        basal_reg_bed_upstream.saveas(basal_reg_bed_upstream_file.name)

        for line in basal_reg_bed_upstream:

            gene_id = line.name
            strand = line.strand
            end = line.end
            start = line.start
            gene_id = "|".join([gene_id, str(start), str(end), strand])
            chroms[gene_id] = line.chrom
            strands[gene_id] = strand

            if strand == '+':

                # if the feature chromosome in B is
                # '.' we have reached the start of the chr
                if line.fields[6] == '.':
                    regulatory_region_start[gene_id] = max(0, line.start - distal)
                else:
                    padding = min(distal, abs(int(line.fields[12])) - 1)
                    regulatory_region_start[gene_id] = line.start - padding

            elif strand == '-':
                # if the feature chromosome in B is
                # '.' we have reached the end of the chr
                if line.fields[6] == '.':
                    regulatory_region_end[gene_id] = min(int(chrom_len[line.chrom]), line.end + distal)
                else:
                    padding = min(distal, abs(int(line.fields[12])) - 1)
                    regulatory_region_end[gene_id] = line.end + padding
            else:
                message("Cannot process genes without strand", type="WARNING")
                message("Please check:" + gene_id, type="ERROR")

        # -------------------------------------------------------------------------
        # Search for downstream limits of each basal_reg_bed
        # Here we ignore overlapping  basal_reg_bed as the way they
        # are proceded is not documented in GREAT to our knowledge
        # -------------------------------------------------------------------------
        message("Defining regulatory domains downstream regions.", type="INFO")

        basal_reg_bed_downstream = basal_reg_bed.closest(basal_reg_bed,
                                                         # Ignore features in B that overlap A
                                                         io=True,
                                                         # In addition to the closest feature in B report distance
                                                         # use negative distances to report upstream features.
                                                         # Report distance with respect to A.
                                                         # When A is on the - strand, "upstream" means B has a
                                                         # higher(start, stop).
                                                         D="a",
                                                         # Ignore features in B that are upstream of features in A
                                                         iu=True,
                                                         # How ties are handled. "first"  Report the first tie
                                                         t="first",
                                                         # Require that the query and the closest hit have different names/gene_ids.
                                                         N=True)

        basal_reg_bed_downstream_file = make_tmp_file(prefix='basal_reg_bed_upstream', suffix='.bed')
        basal_reg_bed_downstream.saveas(basal_reg_bed_downstream_file.name)

        for line in basal_reg_bed_downstream:

            gene_id = line.name
            strand = line.strand
            end = line.end
            start = line.start
            gene_id = "|".join([gene_id, str(start), str(end), strand])
            chroms[gene_id] = line.chrom
            strands[gene_id] = strand

            if strand == '+':
                # if the feature chromosome in B is
                # '.' we have reached the start of the chr
                if line.fields[6] == '.':
                    regulatory_region_end[gene_id] = min(int(chrom_len[line.chrom]), line.end + distal)
                else:
                    padding = min(distal, abs(int(line.fields[12])) - 1)
                    regulatory_region_end[gene_id] = line.end + padding
            elif strand == '-':
                if line.fields[6] == '.':
                    # sys.stderr.write(str(line.start - distal + 1) + "\n")
                    # sys.stderr.write(gene_id + "\n")
                    regulatory_region_start[gene_id] = max(0, line.start - distal)
                else:
                    padding = min(distal, abs(int(line.fields[12])) - 1)
                    regulatory_region_start[gene_id] = max(0, line.start - padding)
            else:
                message("Cannot process genes without strand", type="WARNING")
                message("Please check:" + gene_id, type="ERROR")
            # print(regulatory_region_start)

    else:
        message("Only 'basal_plus_extension' association rule is currently supported.", type='ERROR')

    # -------------------------------------------------------------------------
    # Print the regulatory regions of all genes
    # By default print all genes
    # -------------------------------------------------------------------------

    if go_id is None:
        for gene_id in regulatory_region_start:
            outlist = [chroms[gene_id],
                       str(regulatory_region_start[gene_id]),
                       str(regulatory_region_end[gene_id]),
                       gene_id.split("|")[0],
                       "0",
                       strands[gene_id]]

            outputfile.write("\t".join(outlist) + "\n")
    else:

        # -------------------------------------------------------------------------
        # Get the list of gene/transcript associated with a particular GO term
        # -------------------------------------------------------------------------

        message("Getting Gene Ontology annotations.")

        if not go_id.startswith("GO:"):
            go_id = "GO:" + go_id

        is_associated = set()

        bm = Biomart(http_proxy=http_proxy,
                     https_proxy=https_proxy)

        bm.get_datasets('ENSEMBL_MART_ENSEMBL')

        if species + "_gene_ensembl" not in bm.datasets:
            message("Unknow dataset/species.", type="ERROR")

        bm.query({'query': XML.format(species=species, go=go_id)})

        for i in bm.response.content.decode().split("\n"):
            i = i.rstrip("\n")
            if i != '':
                is_associated.add(i)

        for gene_id in regulatory_region_start:
            gene_id_short = gene_id.split("|")[0]
            if gene_id_short in is_associated:
                outlist = [chroms[gene_id],
                           str(regulatory_region_start[gene_id]),
                           str(regulatory_region_end[gene_id]),
                           gene_id.split("|")[0],
                           "0",
                           strands[gene_id]]

                outputfile.write("\t".join(outlist) + "\n")


def main():
    """main"""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    great_reg_domains(**args)


if __name__ == '__main__':
    main()

else:
    test = """

    #great_reg_domains: load the dataset
    @test "great_reg_domains_0" {
     result=`gtftk get_example -f '*' -d simple`
      [ "$result" = "" ]
    }

    #great_reg_domains: check that this dataset always provide the same result
    @test "great_reg_domains_1" {
     result=`gtftk great_reg_domains -i simple.gtf -c simple -u 1 -d 1 -t 20 | md5 -r | sed 's/ .*//'`
      [ "$result" = "7c91eb8597871f0c9d79a51566be0a38" ]
    }

    #great_reg_domains: check simply the number of lines
    @test "great_reg_domains_2" {
     result=`gtftk great_reg_domains -i simple.gtf -c simple -u 1 -d 1 -t 20 | wc -l | perl -npe  's/\\s+//'`
      [ "$result" = "10" ]
    }    
    """

    CmdObject(name="great_reg_domains",
              message="Attempt to compute labeled regions using GREAT 'association rule'",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="miscellaneous",
              desc=__doc__,
              updated=__updated__,
              notes=__notes__,
              test=test)
