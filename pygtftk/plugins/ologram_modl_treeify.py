#!/usr/bin/env python
"""
Turns a result of OLOGRAM-MODL multiple overlap (tsv file) in a tree for easier visualisation

See the /pygtftk/plugins/ologram.py file, as well as the documentation, for more information about OLOGRAM.
"""

import argparse

import numpy as np
import os
import pandas as pd
import warnings
from distutils.spawn import find_executable

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.stats.intersect.modl import tree
from pygtftk.utils import message

__updated__ = ''' 2020-06-17 '''

__notes__ = """
 -- Turns a result of OLOGRAM-MODL multiple overlap (tsv file) in a tree for easier visualisation.

 -- This is the preferred representation for OLOGRAM-MODL results. Each node represents 
 a combination, with its number of overlapping basepairs in true data (S) and the corresponding 
 fold change and p-value compared to the shuffles.

 -- Result tsv files can be manually edited (ie. removing combinations) before passing them to this plugin

 -- For a quick filtering, it is possible to show only the top T combinations sorted by total basepairs in real data.
"""


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Complete path to the OLOGRAM output file",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('tsv')),
                            required=True)

    parser_grp.add_argument('-o', '--output',
                            help="Pdf file name",
                            default=None,
                            nargs=None,
                            type=arg_formatter.FormattedFile(mode='w', file_ext='pdf'),
                            required=True)

    parser_grp.add_argument('-l', '--query-label',
                            help="Name of the query for display",
                            default="Query",
                            type=str,
                            required=False)

    parser_grp.add_argument('-t', '--top-s',
                            help="Optional. Only the top t combinations sorted by total basepairs in real data will be displayed.",
                            default=-1,
                            type=int,
                            required=False)

    parser_grp.add_argument('-mh', '--min_inheritance',
                            help="Optional. A combination will be added to a shorter parent in the tree only if it represents at least a proportion MH of its true base pairs (between 0 and 1).",
                            default=-1,
                            type=float,
                            required=False)


    return parser


def ologram_modl_treeify(inputfile=None, output=None,
    query_label="Query", top_s = -1, min_inheritance = 0):
    # -------------------------------------------------------------------------
    # Check graphviz is installed
    # to avoid ugly error messages.
    # -------------------------------------------------------------------------
    if not find_executable("dot"):
        message("gtftk ologram_modl_treeify needs graphviz to be installed on your system.",
                type="ERROR")

    # -------------------------------------------------------------------------
    # Reading the input dataframe
    # -------------------------------------------------------------------------

    message("Reading OLOGRAM-MODL results from : " + inputfile.name)

    # Read dataframe to create the found_combis dictionary
    df_res = pd.read_csv(inputfile.name, sep='\t', header=0, index_col=None)
    # Pval set to 0 or -1 are changed to 1e-320 and NaN respectively
    df_res.loc[df_res['summed_bp_overlaps_pvalue'] == 0, 'summed_bp_overlaps_pvalue'] = 1e-320
    df_res.loc[df_res['summed_bp_overlaps_pvalue'] == -1, 'summed_bp_overlaps_pvalue'] = np.nan

    # Optional : use only the top T combinations, sorted by total basepairs in real data
    if top_s != -1:
        df_res = df_res.nlargest(top_s, 'summed_bp_overlaps_true').reset_index(drop = True)


    # -------------------------------------------------------------------------
    # Produce and save the visualization
    # -------------------------------------------------------------------------

    # Create a Library tree with the combinations, like in the MODL algorithm itself    
    L = tree.Library()
    L.build_nodes_for_words_from_ologram_result_df(df_res, query_label)
    L.assign_nodes(min_inheritance = min_inheritance)

    # Now call the function to produce the visualisation
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree.output_visualize(L, output.name)
        message("Saving diagrams to files with root : " + output.name)


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    ologram_modl_treeify(**args)


if __name__ == '__main__':
    main()

else:

    # 'Bats' tests
    test = '''
    #ologram: get example files
    @test "ologram_modl_treeify_0" {
         result=`gtftk get_example -d ologram_2 -f '*'`
      [ "$result" = "" ]
    }

    #ologram: produce tree
    @test "ologram_modl_treeify_1" {
         result=`gtftk ologram_modl_treeify -i multiple_overlap_trivial_ologram_stats.tsv -o tree_multi.pdf -l QueryLabel`
      [ "$result" = "" ]
    }
    '''

    cmd = CmdObject(name="ologram_modl_treeify",
                    message="Build a tree representation from an OLOGRAM-MODL multiple combinations result files (tsv).",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="ologram",
                    notes=__notes__,
                    updated=__updated__,
                    test=test)
