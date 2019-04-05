#!/usr/bin/env python

import argparse
import glob
import os
import shutil
import sys

import pygtftk
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import get_example_file
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
 Print example files including GTF.
"""

__notes__ = '''

-- Use format '*' to get all files from a dataset.
'''


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-d', '--dataset',
                            help="Select a dataset.",
                            type=str,
                            choices=[
                                "simple",
                                "mini_real",
                                "mini_real_noov_rnd_tx",
                                "tiny_real",
                                "hg38_chr1",
                                "simple_02",
                                "simple_03",
                                "simple_04",
                                "simple_05",
                                "simple_06",
                                "mini_real_10M",
                                "control_list"],
                            default="simple",
                            required=False)

    parser_grp.add_argument('-f', '--format',
                            help="The dataset format.",
                            type=str,
                            choices=["*", "gtf", "bed", "bw", "bam", "join",
                                     "join_mat", "chromInfo",
                                     "tsv",
                                     "fa", "fa.idx", "genes", "geneList",
                                     "2.bw", "genome"],
                            default="gtf",
                            required=False)

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="OUTPUT",
                            type=argparse.FileType('w'))

    parser_grp.add_argument('-l', '--list',
                            help="Only list files of a dataset.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-a', '--all-dataset',
                            help="Get the list of all datasets.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-q', '--quiet',
                            help="Don't write any message when copying files.",
                            action="store_true",
                            required=False)
    return parser


def get_example(outputfile=None,
                dataset=None,
                format="gtf",
                quiet=False,
                list=False,
                all_dataset=False):
    """
    Print example gtf files.
    """

    message("Printing example...")

    if all_dataset:
        message("The following datasets were found:")
        path_dataset = os.path.join(pygtftk.__path__[0],
                                    "data")
        dataset_list = glob.glob(os.path.join(path_dataset, "*"))
        dataset_list = [x for x in dataset_list if "__init__" not in x]
        for i in dataset_list:
            print("\t- " + os.path.basename(i))

    elif list:

        path_dataset = os.path.join(pygtftk.__path__[0],
                                    "data",
                                    dataset)

        file_dataset = os.listdir(path_dataset)

        for i in file_dataset:
            print("\t" + i)
    else:
        if format == "gtf":
            try:
                gtf = GTF(get_example_file(datasetname=dataset,
                                           ext=format)[0], check_ensembl_format=False)
            except:
                try:
                    gtf = GTF(get_example_file(datasetname=dataset,
                                               ext=format + ".gz")[0], check_ensembl_format=False)
                except:
                    message("No GTF file found for this dataset.",
                            type="ERROR")

            gtf.write(outputfile, gc_off=True)

        elif format in ["fa", "join", "join_mat", "genome", "chromInfo", "genes", "geneList"]:
            try:
                infile = open(get_example_file(datasetname=dataset,
                                               ext=format)[0], "r")
            except:
                message("Unable to find example file.", type="ERROR")

            for line in infile:
                outputfile.write(line)

        elif format == "*":

            file_path = glob.glob(os.path.join(pygtftk.__path__[0],
                                               'data',
                                               dataset,
                                               "*"))
            file_path = [x for x in file_path if "__" not in x]
            target_path = os.path.join(pygtftk.__path__[0], 'data', dataset)
            if not quiet:
                message("Copying from :" + target_path)

            for i in file_path:
                if not os.path.exists(os.path.join(os.getcwd(), os.path.basename(i))):
                    if not quiet:
                        message("Copying file : " + os.path.basename(i), force=True)
                    shutil.copy(i, ".")
                else:
                    if not quiet:
                        message("Copy canceled, file already exist: " + os.path.basename(i), force=True)

        else:

            if not quiet:
                message("Copying ", force=True)
            file_path = glob.glob(os.path.join(pygtftk.__path__[0],
                                               'data',
                                               dataset,
                                               "*" + '.' + format))

            file_path = [x for x in file_path if "__" not in x]
            target_path = os.path.join(pygtftk.__path__[0], 'data', dataset)
            message("Copying from :" + target_path)

            for i in file_path:
                if not os.path.exists(os.path.join(os.getcwd(), os.path.basename(i))):
                    if not quiet:
                        message("Copying file : " + os.path.basename(i), force=True)
                    shutil.copy(i, ".")
                else:
                    if not quiet:
                        message("Copy canceled, file already exist:" + os.path.basename(i), force=True)

    close_properly(outputfile)


def main():
    """The main function."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    get_example(**args)


if __name__ == '__main__':
    main()

else:

    test = """
    #get_example
    @test "get_example_1" {
     result=`gtftk get_example| wc -l`
      [ "$result" -eq 70 ]
    }
    """

    CmdObject(name="get_example",
              message="Get example files including GTF.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              desc=__doc__,
              notes=__notes__,
              group="information",
              updated=__updated__,
              test=test)
