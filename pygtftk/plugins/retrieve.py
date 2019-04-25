#!/usr/bin/env python
"""
 Retrieve a GTF file from ensembl.
"""

import argparse
import os
import re
import sys

import ftputil
from ftputil.error import FTPOSError

import pygtftk
from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import message

__updated__ = "2018-01-31"


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-s', '--species-name',
                            help="The species name.",
                            type=str,
                            metavar="SPECIES",
                            default="homo_sapiens",
                            required=False)

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file (gtf.gz).",
                            metavar="GTF.GZ",
                            type=arg_formatter.FormattedFile(mode='w', file_ext='gtf.gz'))

    parser_grp.add_argument('-e', '--ensembl-collection',
                            help="Which ensembl collection to interrogate"
                                 "('vertebrate', 'protists', 'fungi', 'plants', 'metazoa').",
                            default="vertebrate",
                            choices=[
                                'vertebrate',
                                'protists',
                                'fungi',
                                'plants',
                                'metazoa'],
                            type=str)

    parser_grp.add_argument('-r', '--release',
                            help="Release number. By default, the latest.",
                            default=None,
                            metavar="RELEASE",
                            type=str)

    parser_grp.add_argument('-l', '--list-only',
                            help="If selected, list available species.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-hs', '--hide-species-name',
                            help="If selected, hide species names upon listing.",
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-c', '--to-stdout',
                            help='This option specifies that output will go to '
                                 'the standard output stream, leaving downloaded'
                                 ' files intact (or not, see -d).',
                            action='store_true',
                            required=False)

    parser_grp.add_argument('-d', '--delete',
                            help='Delete the gtf file after processing (e.g if -c is used).',
                            action='store_true',
                            required=False)

    return parser


def retrieve(species_name='homo_sapiens',
             outputfile=None,
             release=None,
             to_stdout=False,
             list_only=False,
             delete=False,
             hide_species_name=None,
             ensembl_collection='vertebrate'):
    """Retrieve a GTF file from ensembl.

    :Example:
    >>> # retrieve("Xenopus_tropicalis")
    """

    if outputfile is None:
        outputdir = os.getcwd()
    else:
        outputdir = os.path.dirname(os.path.abspath(outputfile.name))

    if species_name is None and not list_only:
        message("Choose --species-name or --list-only.", type='ERROR')

    if outputfile is not None:
        if not os.path.exists(
                os.path.dirname(os.path.abspath(outputfile.name))):
            message("Output directory does not exists. Exiting.", type='ERROR')
        else:
            if os.path.isdir(outputfile.name):
                message("Output file is a directory !.", type='ERROR')

    # Will contain the url pointing to the
    # requested gtf.
    target_gtf = None

    # -------------------------------------------------------------------------
    # Check ensembl repository
    # -------------------------------------------------------------------------

    if ensembl_collection == 'vertebrate':

        host = "ftp.ensembl.org"
        user = "anonymous"  # votre identifiant
        password = "anonymous@gtftk.fr"
    elif ensembl_collection in ['protists', 'fungi', 'plants', 'metazoa']:

        host = "ftp.ensemblgenomes.org"
        user = "anonymous"
        password = "anonymous@gtftk.fr"

    try:
        ftp = ftputil.FTPHost(host, user, password)
        if pygtftk.utils.VERBOSITY:
            message("Connected to ensembl FTP website.")
    except FTPOSError as err:
        message(str(err))
        message("Unable to connect (FTPOSError).", type="ERROR")

    try:
        ftp.chdir('/pub')
        message("Successfully change directory to pub")
    except:
        message("Unable to change directory to 'pub'.",
                type="ERROR")

    if ensembl_collection in ['protists', 'fungi', 'plants', 'metazoa']:
        try:
            ftp.chdir(ensembl_collection)
            message("Successfully change directory to " + ensembl_collection)
        except:
            message("Unable to change directory to '%s'." % ensembl_collection,
                    type="ERROR")

    try:
        all_releases = ftp.listdir(ftp.curdir)
    except Exception as e:
        print(str(e))
        message("Unable to list directory.",
                type="ERROR")

    if release is not None:
        release_dir = "release-" + release
        if release_dir not in all_releases:
            message("This release number could not be found. Aborting",
                    type="ERROR")
    else:

        version_list = []

        for ver in all_releases:
            regexp = re.compile("release-(\d+)")
            hit = regexp.search(ver)
            if hit:
                version_list += [int(hit.group(1))]
        release = max(version_list)
        release_dir = "release-" + str(release)
        message("Latest version is %d." % release)

    try:
        ftp.chdir(release_dir)
        message("Changed release directory: %s" % release_dir,
                type="DEBUG")
    except:
        message("Unable to change directory to '%s'." % release_dir,
                type="ERROR")

    ftp.chdir('gtf')

    try:
        all_species = ftp.listdir(ftp.curdir)
        all_species = [x for x in all_species if ftp.path.isdir(x)]
    except:
        message("Unable to list directory.",
                type="ERROR")

    if list_only:

        species_list = []
        url_list = []

        for sp in all_species:
            gtfs = [x for x in ftp.listdir(sp) if x.endswith('.gtf.gz')]

            for gtf in gtfs:
                species_list += [sp]
                current_url = 'ftp://' + host + ftp.getcwd() + '/'
                url_list += [current_url + sp + "/" + gtf]

        for sp, url in zip(species_list, url_list):
            if hide_species_name:
                print(url)
            else:
                print(sp.ljust(50) + url)

        sys.exit()
    else:

        if species_name not in all_species:
            message("Species could not be found for release: %s" % str(release))
            message("Trying species name in lower case.")
            species_name = species_name.lower()
            if species_name not in all_species:
                message("Species could not be found for release: %s" % str(release),
                        type="ERROR")

        ftp.chdir(species_name)
        gtf_list = ftp.listdir(ftp.curdir)

        # choice 1 (only regular chromosome)
        gtf_sub = [x for x in gtf_list if x.endswith("chr.gtf.gz")]

        # choice 2 should be ! choice 1 and ! 'ab_initio'.
        # Should be default gtf
        gtf_sub_2 = [x for x in gtf_list if "abinitio.gtf.gz" not in x]
        gtf_sub_2 = [x for x in gtf_sub_2 if x.endswith(".gtf.gz")]
        if gtf_sub:
            gtf_sub_2.remove(gtf_sub[0])

        # Choice 3 abinitio
        gtf_sub_3 = [x for x in gtf_list if x.endswith("abinitio.gtf.gz")]
        # Choice 4:
        # Any gtf

        if len(gtf_sub) > 0:
            target_gtf = gtf_sub[0]
        elif len(gtf_sub_2) > 0:
            target_gtf = gtf_sub_2[0]
        elif len(gtf_sub_3) > 0:
            target_gtf = gtf_sub_3[0]
        else:
            gtf_sub = [x for x in gtf_list if x.endswith(".gtf.gz")]
            target_gtf = gtf_sub[0]

    # -------------------------------------------------------------------------
    # Download if requested
    # -------------------------------------------------------------------------

    if target_gtf is not None:
        if not list_only:
            message("Downloading GTF file : " + target_gtf)

            ftp.download(target_gtf,
                         target_gtf)

            os.rename(target_gtf, os.path.join(outputdir, target_gtf))

            if to_stdout:
                gtf = GTF(os.path.join(outputdir, target_gtf))

                gtf.write("-", gc_off=True)

            if delete:
                os.remove(os.path.join(outputdir, target_gtf))
            else:
                if outputfile is not None:
                    message("Renaming.")
                    os.rename(os.path.join(outputdir, target_gtf),
                              outputfile.name)

    else:
        message("Species could not be found for release: " + release,
                type='ERROR')


def main():
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    retrieve(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #retrieve
    @test "retrieve_1" {
     result=`gtftk retrieve -V 3 -s trypanosoma_brucei -e protists -r 34 -cd  | wc -l`
      [ "$result" -eq 54607 ]
    }

    #retrieve
    @test "retrieve_2" {
     result=`gtftk retrieve -V 3  -l -r 87 | wc -l`
      [ "$result" -eq 218 ]
    }

    #retrieve
    @test "retrieve_3" {
     result=`gtftk retrieve -V 3  -l -r 33 -e protists | wc -l`
      [ "$result" -eq 74 ]
    }

    #retrieve
    @test "retrieve_4" {
     result=`gtftk retrieve -V 3  -e fungi -s trichoderma_reesei -r 34 -cd | wc -l`
      [ "$result" -eq 93776 ]
    }

    #retrieve
    @test "retrieve_5" {
     result=`gtftk retrieve -V 3  -s sorex_araneus -r 87 -cd| wc -l`
      [ "$result" -eq 401171 ]
    }

    """

    CmdObject(name="retrieve",
              message="Retrieve a GTF file from ensembl.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              group="information",
              updated=__updated__,
              test=test)
