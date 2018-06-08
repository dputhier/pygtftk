# -*- coding: utf-8 -*-
"""Select columns from a tabulated file based on their names."""

import argparse
import os
import pickle

import sys

reload(sys)
sys.setdefaultencoding('utf8')

import cloudpickle
import requests
from json2html import *

from gtftk.arg_formatter import FileWithExtension
from gtftk.utils import message

__updated__ = "2018-01-20"

__doc__ = """
  Search encode project (http://www.encodeproject.org/) for experiments.
"""

__notes__ = """
  -- Each contrains applied on biosample-type, biosample-summary (...) are applied using the 'and' logical operator.
  -- The column names of the output matrix are the following: 'id', 'accession', 'organism', 'target_gene_name', 
  'target_label', 'project', 'level_name', 'description', 'biosample_summary', 'biosample_type', 'biosample_term_name',
  'aliases', 'laboratory', 'audit_warning', 'audit_internal_action', 'audit_error', 'replicates'.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-u', '--update',
                            help='Update the ENCODE dataset (i.e download new version)',
                            action="store_true")

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="TXT",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Jj][Ss][Oo][Nn]')))
    parser_grp.add_argument('-a', '--assay-term-name',
                            help='Assay type (e.g. ChIP).',
                            type=str,
                            required=False,
                            choices=['ChIP-seq'],
                            default="ChIP-seq",
                            nargs='+')

    parser_grp.add_argument('-b', '--biosample-type',
                            help='Only experiments with biosample type containing one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-y', '--biosample-summary',
                            help='Only experiments with biosample summary containing one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-d', '--description',
                            help='Only experiments with description containing one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-p', '--project',
                            help='Only experiments whose project name contains one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-l', '--laboratory',
                            help='Only experiments whose lab name contains one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-x', '--aliases',
                            help='Only if experiment alias contains one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-c', '--accession',
                            help='Only if experiment accession contains one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-t', '--target-label',
                            help='Only if experiment "target label" entry one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-g', '--target-gene-name',
                            help='Only if experiment "target gene name" contains one of these strings (csv).',
                            type=str,
                            required=False,
                            default=None)

    parser_grp.add_argument('-j', '--json-file-to-html',
                            help="Print the json file (restricted to experiments) into this html file.",
                            default=None,
                            metavar="GTF/BED",
                            type=FileWithExtension('w',
                                                   valid_extensions=('\.[Hh][Tt][Mm]([Ll])?$')))

    parser_grp.add_argument('-m', '--organism',
                            help='The organism scientific name (e.g "Homo sapiens").',
                            type=str,
                            default=None,
                            required=False)

    parser_grp.add_argument('-n', '--accession-only',
                            help='Print only accession number of experiments.',
                            action="store_true")
    return parser


def search_encode_project(outputfile=None,
                          update=False,
                          assay_term_name="",
                          biosample_type=None,
                          biosample_summary=None,
                          description=None,
                          laboratory=None,
                          project=None,
                          accession_only=False,
                          organism=None,
                          target_label=None,
                          target_gene_name=None,
                          accession=None,
                          aliases=None,
                          json_file_to_html=None,
                          tmp_dir=None,
                          logger_file=None,
                          verbosity=None

                          ):
    """Search and retrieve data from the encode project."""

    # -------------------------------------------------------------------------
    # Get experiment list from encode if not previously dlded
    # -------------------------------------------------------------------------

    config_dir_user = os.path.join(os.path.expanduser("~"), ".gtftk")
    pick_path = os.path.join(config_dir_user, "search_encode_project.pick")

    reload_data = True

    if not os.path.exists(pick_path) or update:
        reload_data = False
        message("No json file found. Downloading json file from ENCODE.", type="WARNING")

        # Force return from the server in JSON format
        HEADERS = {'accept': 'application/json'}

        # Retrieve all samples
        URL = "https://www.encodeproject.org/search/?format=json&limit=all"

        go_on = True

        try:
            response = requests.get(URL, headers=HEADERS)
        except requests.exceptions.Timeout as err:
            message("Encountered requests.exceptions.Timeout", type="ERROR")

        except requests.exceptions.TooManyRedirects:
            message("Encountered requests.exceptions.TooManyRedirects", type="ERROR")

        except requests.exceptions.RequestException as e:
            message("Encountered requests.exceptions.RequestException", type="ERROR")

        # Extract the JSON response as a python dict
        response_json_dict = response.json()

        f_handler = open(pick_path, "w")
        pick = cloudpickle.CloudPickler(f_handler)
        pick.dump((response_json_dict))
        message("Json file stored in :" + pick_path)
        f_handler.close()

    # -------------------------------------------------------------------------
    # Load experiment from file if previously dumped
    # -------------------------------------------------------------------------

    message("assay_term_name = " + assay_term_name)

    if reload_data:
        message("Loading ENCODE data.")
        message("Opening file " + pick_path)
        f_handler = open(pick_path, "r")
        try:
            response_json_dict = pickle.load(f_handler)
        except EOFError as error:
            message(error)
            message("Encountered EOFError error", type="ERROR")

    # -------------------------------------------------------------------------
    # Subsetting to experiment
    # -------------------------------------------------------------------------

    response_json_dict["@graph"] = [k for k in response_json_dict["@graph"] if 'Experiment' in k['@type']]

    # -------------------------------------------------------------------------
    # what we want to extract
    # -------------------------------------------------------------------------

    info = ['@id', 'accession', 'organism', 'target_gene_name', 'target_label', 'project', 'level_name',
            'description', 'biosample_summary', 'biosample_type', 'biosample_term_name',
            'aliases', 'laboratory', 'audit_warning', 'audit_internal_action', 'audit_error', 'replicates', 'organism']

    # -------------------------------------------------------------------------
    # Searching relevant info in json DB
    # -------------------------------------------------------------------------

    message("Searching json DB.")

    if accession_only:
        outputfile.write("accession\n")
    else:
        outputfile.write("\t".join([x.replace('@', '') for x in info]) + "\n")

    for pos, rec in enumerate(response_json_dict["@graph"]):

        if assay_term_name == rec["assay_term_name"]:

            out_list = dict()

            for k in info:
                if k in ['audit_internal_action', 'audit_warning', 'audit_error']:
                    if k == 'audit_internal_action':
                        key_audit = 'INTERNAL_ACTION'
                    elif k == 'audit_warning':
                        key_audit = 'WARNING'
                    elif k == 'audit_error':
                        key_audit = 'ERROR'
                    try:
                        audit_list = []
                        for a in rec['audit'][key_audit]:
                            string = "category=" + a['category'].decode('utf-8') + "|level=" + str(a['level'])
                            audit_list += [string]
                        out_list[k] = ";".join(audit_list)
                    except:
                        out_list[k] = 'NA'

                elif k == 'replicates':
                    try:
                        rep_list = []
                        for a in rec['replicates']:
                            string = "id=" + a['@id'].decode('utf-8') + "|stage=" + a['library']['biosample'][
                                'life_stage']
                            rep_list += [string]
                        out_list[k] = ";".join(rep_list)
                    except:
                        out_list[k] = 'NA'

                    try:
                        rep_list = []
                        for a in rec['replicates']:
                            string = a['library']['biosample']['organism']['scientific_name']
                            rep_list += [string]
                        out_list['organism'] = ";".join(rep_list)
                    except:
                        out_list[k] = 'NA'


                elif k == 'project':
                    try:
                        out_list['project'] = rec['award']['project'].decode('utf - 8')
                    except:
                        out_list['project'] = 'NA'

                elif k == 'lab':
                    try:
                        out_list['laboratory'] = rec['lab']['title'].decode('utf - 8')
                    except:
                        out_list['laboratory'] = 'NA'

                elif k == 'aliases':
                    try:
                        out_list['aliases'] = rec['lab']['title'].decode('utf - 8')
                    except:
                        out_list['aliases'] = 'NA'

                elif k == 'target_gene_name':
                    try:
                        out_list['target_gene_name'] = rec['target']['gene_name'].decode('utf - 8')
                    except:
                        out_list['target_gene_name'] = 'NA'

                elif k == 'target_label':
                    try:
                        out_list['target_label'] = rec['target']['label'].decode('utf - 8')
                    except:
                        out_list['target_label'] = 'NA'

                elif k == 'organism':
                    # see 'replicates'
                    pass
                else:
                    try:
                        out_list[k] = rec[k].decode('utf-8')
                    except:
                        out_list[k] = 'NA'

            cond = ((biosample_type is None or
                     [x for x in biosample_type.split(',') if x.lower() in out_list['biosample_type'].lower()]) and
                    (biosample_summary is None or
                     [x for x in biosample_summary.split(',') if
                      x.lower() in out_list['biosample_summary'].lower()]) and
                    (description is None or
                     [x for x in description.split(',') if x.lower() in out_list['description'].lower()]) and
                    (project is None or
                     [x for x in project.split(',') if x.lower() in out_list['project'].lower()]) and
                    (aliases is None or
                     [x for x in aliases.split(',') if x.lower() in out_list['aliases'].lower()]) and
                    (accession is None or
                     [x for x in accession.split(',') if x.lower() in out_list['accession'].lower()]) and
                    (laboratory is None or
                     [x for x in laboratory.split(',') if x.lower() in out_list['laboratory'].lower()]) and
                    (target_label is None or
                     [x for x in target_label.split(',') if x.lower() in out_list['target_label'].lower()]) and
                    (organism is None or
                     [x for x in organism.split(',') if x.lower() in out_list['organism'].lower()]) and
                    (target_gene_name is None or
                     [x for x in target_gene_name.split(',') if
                      x.lower() in out_list['target_gene_name'].lower()])

                    )

            if cond:
                if accession_only:
                    outputfile.write(out_list['accession'] + "\n")
                else:
                    full_list = [str(out_list[k]) for k in info]
                    outputfile.write("\t".join(full_list) + "\n")

    # -------------------------------------------------------------------------
    # Convert json file to html if requested
    # -------------------------------------------------------------------------

    if json_file_to_html is not None:
        json_file_to_html.write(json2html.convert(json=response_json_dict))


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    search_encode_project(**args)


if __name__ == '__main__':
    main()


# TODO
else:
    test = """
    #search_encode_project
    @test "search_encode_project_1" {
     result=`gtftk get_example | wc -l`
      [ "$result" -eq 75 ]
    }

    """
    from gtftk.cmd_object import CmdObject

    CmdObject(name="search_encode_project",
              message="Search data from https://www.encodeproject.org.",
              parser=make_parser(),
              fun=search_encode_project,
              updated=__updated__,
              desc=__doc__,
              notes=__notes__,
              group="miscellaneous",
              test=test)
