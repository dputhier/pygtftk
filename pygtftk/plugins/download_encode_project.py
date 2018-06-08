# -*- coding: utf-8 -*-
"""Select columns from a tabulated file based on their names."""

import argparse
import json
import multiprocessing
import os
import re
import sys
from collections import defaultdict
from itertools import repeat

reload(sys)
sys.setdefaultencoding('utf8')

import requests
import wget
from json2html import *

from gtftk.arg_formatter import FileWithExtension
from gtftk.arg_formatter import int_greater_than_null
from gtftk.utils import chomp
from gtftk.utils import intervals
from gtftk.utils import make_outdir_and_file
from gtftk.utils import make_tmp_file
from gtftk.utils import message
from gtftk.utils import to_alphanum, call_nested_dict_from_list

__updated__ = "2018-01-20"

__doc__ = """
  Download files from encode project (http://www.encodeproject.org/).
"""

__notes__ = """
-- You need to to provide an input file obtained from gtftk search_encode_project.   
-- Currently output file name (see -m) can be constructed by concatening the following information: 'biosample_summary', 
'assay_title', 'target_label', 'accession', 'target_gene_name','technical_replicates' and 'biological_replicates', 'organism'.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="A list of IDs (encodeproject.org).",
                            default=sys.stdin,
                            metavar="TXT",
                            type=FileWithExtension('r',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Jj][Ss][Oo][Nn]')),
                            required=False)

    parser_grp.add_argument('-o', '--outputdir',
                            help='Output directory name.',
                            metavar="DIR",
                            default="encode_dataset",
                            type=str)

    parser_grp.add_argument('-f', '--output-type',
                            help='The output type to be downloaded.',
                            type=str,
                            choices=['alignments', 'reads', 'unfiltered alignments', 'peaks',
                                     'signal p-value', 'stable peaks',
                                     'fold change over control'],
                            required=False,
                            default="bam")

    parser_grp.add_argument('-k', '--nb-proc',
                            help='The number of parallel threads.',
                            type=int_greater_than_null,
                            required=False,
                            default=1)

    parser_grp.add_argument('-j', '--to-json',
                            help='Dump sample info into a json file (use -K to indicate a folder).',
                            action="store_true")

    parser_grp.add_argument('-a', '--assembly',
                            help='The assembly (e.g hg19, GRCh38...).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-n', '--dry-run',
                            help="Don't actually download file. Just print URLs.",
                            action="store_true")

    parser_grp.add_argument('-m', '--file-name-info',
                            help='Defines the info used to construct file name. See notes',
                            type=str,
                            default="biosample_summary,target_label,accession,technical_replicates,biological_replicates,run_type",
                            required=False)

    parser_grp.add_argument('-s', '--file-name-sep',
                            help='Use to separate info in file name.',
                            type=str,
                            default="~",
                            required=False)

    return parser


def _dld_encode_worker(input_values):
    """
    A function to download a set of encode files.
    """
    (span, output_type, smp_list, to_json, assembly, dry_run, outputdir, file_name_info,
     fn_sep, organism, verbose) = input_values

    smp_list = smp_list.split(",")
    file_name_info = file_name_info.split(",")

    for id_file in smp_list[span[0]:span[1]]:

        mesg = "Processing sample %s (format %s) ."
        mesg = mesg % (id_file, output_type)
        message(mesg, type="INFO")

        # Force return from the server in JSON format
        HEADERS = {'accept': 'application/json'}

        # Retrieve all samples
        URL = "https://www.encodeproject.org/experiments/{0}/?format=json".format(id_file)

        go_on = True

        try:
            response = requests.get(URL, headers=HEADERS)
        except requests.exceptions.Timeout as err:
            message("Encountered requests.exceptions.Timeout")
            go_on = False
        except requests.exceptions.TooManyRedirects:
            message("Encountered requests.exceptions.TooManyRedirects")
            go_on = False
        except requests.exceptions.RequestException as e:
            message("Encountered requests.exceptions.RequestException")
            go_on = False

        if not go_on:
            continue

        response_json_dict = response.json()

        if to_json:
            html_file = make_tmp_file(prefix=id_file, suffix=".html")
            html_file.write(json2html.convert(json=response_json_dict))
            html_file.close()
            json_file = make_tmp_file(prefix=id_file, suffix=".json")
            json_file.write(json.dumps(response_json_dict, indent=4, separators=(',', ': ')))
            json_file.close()

        run_list = defaultdict(dict)

        for rec in response_json_dict['files']:

            if rec['output_type'] == output_type:
                do_it = False
                if assembly is not None:
                    if rec['assembly'] == assembly:
                        do_it = True
                else:
                    do_it = True

                    if do_it:

                        try:
                            acc = rec[u'accession']
                        except:
                            acc = rec[u'external_accession']
                        run_list[acc]['href'] = call_nested_dict_from_list(rec, ['href'])

                        run_list[acc]['url'] = 'https://www.encodeproject.org' + run_list[acc]['href']

                        run_list[acc]['technical_replicates'] = call_nested_dict_from_list(rec,
                                                                                           ['technical_replicates'])
                        run_list[acc]['biological_replicates'] = call_nested_dict_from_list(rec,
                                                                                            ['biological_replicates'])
                        run_list[acc]['assay_title'] = call_nested_dict_from_list(response_json_dict,
                                                                                  ['assay_title'])
                        run_list[acc]['target_label'] = call_nested_dict_from_list(response_json_dict,
                                                                                   ['target', 'label'])
                        run_list[acc]['target_gene_name'] = call_nested_dict_from_list(response_json_dict,
                                                                                       ['target', 'gene_name'])
                        run_list[acc]['biosample_summary'] = call_nested_dict_from_list(response_json_dict,
                                                                                        ['biosample_summary'])
                        run_list[acc]['run_type'] = call_nested_dict_from_list(response_json_dict, ['run_type'])
                        run_list[acc]['accession'] = acc
                        run_list[acc]['run_type'] = call_nested_dict_from_list(response_json_dict, ['run_type'])
                        run_list[acc]['organism'] = response_json_dict['replicates'][0]['library'][0]['biosample'][0][
                            'organism']

                        if dry_run:
                            for acc in run_list:
                                message("FOUND : " + run_list[acc]['url'], force=True)


                        else:
                            url_to_name = dict()
                            for acc in run_list:
                                re_search = re.search('(\.[^\.]*?)(\.gz)?$', run_list[acc]['url'])
                                if re_search:
                                    name_file = fn_sep.join([to_alphanum(run_list[acc][k]) for k in file_name_info])
                                    name_file += re_search.group(0)
                                    url_to_name[run_list[acc]['url']] = name_file
                                else:
                                    print("not found" + run_list[acc]['url'])

                            for u, n in url_to_name.items():
                                message("Downloading " + u)
                                try:
                                    wget.download(url=u, out=os.path.join(outputdir, n), bar=None)
                                    message("Completed: " + u)
                                except IOError as err:
                                    message("Problem encountered while downloading : " + u)
                                    print(err)

                    if dry_run:
                        return False
                    else:
                        return True


def download_encode_project(inputfile=None,
                            output_type=None,
                            tmp_dir=None,
                            outputdir=None,
                            to_json=False,
                            organism="9606",
                            logger_file=None,
                            nb_proc=1,
                            assembly="hg19",
                            dry_run=False,
                            file_name_info=None,
                            file_name_sep=None,
                            verbosity=None):
    """Download files from encode project (http://www.encodeproject.org/)."""

    # -------------------------------------------------------------------------
    # Preparing output file
    # -------------------------------------------------------------------------

    file_out_list = make_outdir_and_file(out_dir=outputdir,
                                         alist=["sample_info.txt"],
                                         force=True)

    sample_info_txt = file_out_list

    smp_list = []
    for p, line in enumerate(inputfile):

        line = chomp(line)
        fields = line.split("\t")

        if p == 0:
            if len(fields) > 1:
                multiple_columns = True

                for i in range(len(fields)):

                    pos = fields.index(fields[i]) if fields[i] == "accession" else -1

                    if pos > -1:
                        break

                if pos == -1:
                    message("The matrix contains multiple columns. Column 'accession' not found",
                            type="ERROR")
            else:
                if line != 'accession':
                    smp_list += [line]

        else:

            if multiple_columns:
                smp_list += [fields[pos]]
            else:
                smp_list += [line]

    n_smp_to_proceed = len(smp_list)

    if nb_proc > n_smp_to_proceed:
        nb_proc = n_smp_to_proceed
        message("Using only {0} processors (one per sample).".format(nb_proc), type="WARNING")

    message("Received " + str(n_smp_to_proceed) +
            " sample to be proceeded.")

    tokens = intervals(range(n_smp_to_proceed), nb_proc)

    pool = multiprocessing.Pool(nb_proc)
    ret = pool.map_async(_dld_encode_worker,
                         zip(tokens,
                             repeat(output_type),
                             repeat(",".join(smp_list)),
                             repeat(to_json),
                             repeat(assembly),
                             repeat(dry_run),
                             repeat(outputdir),
                             repeat(file_name_info),
                             repeat(file_name_sep),
                             repeat(organism),
                             repeat(verbosity))).get(9999999)

    if not ret:
        message("Exiting")
        sys.exit()


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    download_encode_project(**args)


if __name__ == '__main__':
    main()


else:
    test = """
    #download_encode_project
    @test "download_encode_project_1" {
     result=`gtftk get_example | wc -l`
      [ "$result" -eq 75 ]
    }
    
    """

    from gtftk.cmd_object import CmdObject

    CmdObject(name="download_encode_project",
              message="Download data from https://www.encodeproject.org.",
              parser=make_parser(),
              fun=download_encode_project,
              updated=__updated__,
              desc=__doc__,
              notes=__notes__,
              group="miscellaneous",
              test=test)
