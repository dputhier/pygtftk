#!/usr/bin/env python
"""
Merge a set of OLOGRAM outputs into a single output. Build a heatmap from the results.

See the /pygtftk/plugins/ologram.py file, as well as the documentation, for more information about OLOGRAM.
"""

import argparse
import os
import re
import warnings

import numpy as np
import pandas as pd
from plotnine import (ggplot, aes, geom_tile, theme_bw,
                      element_blank, theme, geom_point,
                      element_text, save_as_pdf_pages,
                      scale_color_gradientn, labs,
                      scale_fill_gradient2)

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import message

__updated__ = ''' 2019-05-27 '''

__notes__ = """
-- By default, labels in the diagram are derived from the name of the enclosing folder. E.g. if file is a/b/c/00_ologram_stats.tsv, 'c' will be used as label.
-- Otherwise use -\-labels to set the labels.

-- Squares without a diamond mean the p-value was NaN due to poor fitting. This is mostly the case for higher-order combis in multiple overlaps that were so rare that they are not encountered in the shuffles.
"""


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('inputfiles',
                            help="Complete paths to the OLOGRAM output files",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('txt')),
                            nargs='+')

    parser_grp.add_argument('-pw', '--pdf-width',
                            help='Output pdf file width (inches).',
                            type=arg_formatter.ranged_num(0, None),
                            default=None,
                            required=False)

    parser_grp.add_argument('-ph', '--pdf-height',
                            help='Output pdf file height (inches).',
                            type=arg_formatter.ranged_num(0, None),
                            default=None,
                            required=False)

    parser_grp.add_argument('-o', '--output',
                            help="Pdf file name.",
                            default=None,
                            nargs=None,
                            type=arg_formatter.FormattedFile(mode='w', file_ext='pdf'),
                            required=True)

    parser_grp.add_argument('-l', '--labels',
                            help="A comma separated list of labels.",
                            default=None,
                            type=str,
                            required=False)

    return parser



def ologram_merge_stats(inputfiles=None,
                        pdf_width=None,
                        pdf_height=None,
                        output=None,
                        labels=None):
    # -------------------------------------------------------------------------
    # Check user provided labels
    # -------------------------------------------------------------------------

    if labels is not None:

        labels = labels.split(",")

        for elmt in labels:
            if not re.search("^[A-Za-z0-9_]+$", elmt):
                message(
                    "Only alphanumeric characters and '_' allowed for --more-bed-labels",
                    type="ERROR")
        if len(labels) != len(inputfiles):
            message("--labels: the number of labels should be"
                    " the same as the number of input files ", type="ERROR")

        if len(labels) != len(set(labels)):
            message("Redundant labels not allowed.", type="ERROR")

    # -------------------------------------------------------------------------
    # Loop over input files
    # -------------------------------------------------------------------------

    df_list = list()
    df_label = list()

    for pos, infile in enumerate(inputfiles):
        message("Reading file : " + infile.name)
        # Read the dataset into a temporay dataframe
        df_tmp = pd.read_csv(infile, sep='\t', header=0, index_col=None)
        # Change name of 'feature_type' column.
        df_tmp = df_tmp.rename(index=str, columns={"feature_type": "Feature"})
        # Assign the name of the dataset to a new column

        if labels is None:
            file_short_name = os.path.basename(os.path.normpath(os.path.dirname(infile.name)))
            df_label += [file_short_name]
        else:
            file_short_name = labels[pos]
            df_label += [labels[pos]]

        df_tmp = df_tmp.assign(**{"dataset": [file_short_name] * df_tmp.shape[0]})
        # Pval set to 0 or -1 are changed to 1e-320 and NaN respectively
        df_tmp.loc[df_tmp['summed_bp_overlaps_pvalue'] == 0, 'summed_bp_overlaps_pvalue'] = 1e-320
        df_tmp.loc[df_tmp['summed_bp_overlaps_pvalue'] == -1, 'summed_bp_overlaps_pvalue'] = np.nan
        # Compute -log10(pval)
        df_tmp = df_tmp.assign(**{"-log_10(pval)": -np.log10(df_tmp.summed_bp_overlaps_pvalue)})

        # Which p-values are signifcant ?
        # TODO: For now, draws all p-values. Add Benjamini-Hochberg correction, and distinguish between NaN and 0.
        df_tmp = df_tmp.assign(**{"pval_signif": df_tmp.summed_bp_overlaps_pvalue > 0})

        # Add the df to the list to be subsequently merged
        df_list += [df_tmp]



    if len(set(df_label)) < len(df_label):
        message('Enclosing directories are ambiguous and cannot be used as labels. You may use "--labels".',
                type="ERROR")

    # -------------------------------------------------------------------------
    # Concatenate dataframes (row bind)
    # -------------------------------------------------------------------------

    message("Merging dataframes.")
    df_merged = pd.concat(df_list, axis=0)

    # -------------------------------------------------------------------------
    # Plotting
    # -------------------------------------------------------------------------

    message("Plotting")
    my_plot = ggplot(data=df_merged,
                     mapping=aes(y='Feature', x='dataset'))
    my_plot += geom_tile(aes(fill = 'summed_bp_overlaps_log2_fold_change'))
    my_plot += scale_fill_gradient2()
    my_plot += labs(fill = "log2(fold change) for summed bp overlaps")

    # Points for p-val. Must be after geom_tile()
    my_plot += geom_point(data = df_merged.loc[df_merged['pval_signif']],
        mapping = aes(x='dataset',y='Feature',color = '-log_10(pval)'), size=4, shape ='D', inherit_aes = False)
    my_plot += scale_color_gradientn(colors = ["#160E00","#FFB025","#FFE7BD"])
    my_plot += labs(color = "-log10(p-value)")

    # Theming
    my_plot += theme_bw()
    my_plot += theme(panel_grid_major=element_blank(),
                     axis_text_x=element_text(rotation=90),
                     panel_border=element_blank(),
                     axis_ticks=element_blank())

    # -------------------------------------------------------------------------
    # Saving
    # -------------------------------------------------------------------------

    message("Saving")
    nb_ft = len(list(df_merged['Feature'].unique()))
    nb_datasets = len(list(df_merged['dataset'].unique()))

    if pdf_width is None:
        panel_width = 0.6
        pdf_width = panel_width * nb_datasets

        if pdf_width > 100:
            pdf_width = 100
            message("Setting --pdf-width to 100 (limit)")

    if pdf_height is None:
        panel_height = 0.6
        pdf_height = panel_height * nb_ft

        if pdf_height > 500:
            pdf_height = 500
            message("Setting --pdf-height to 500 (limit)")

    message("Page width set to " + str(pdf_width))
    message("Page height set to " + str(pdf_height))
    figsize = (pdf_width, pdf_height)

    # -------------------------------------------------------------------------
    # Turn warning off. Both pandas and plotnine use warnings for deprecated
    # functions. I need to turn they off although I'm not really satisfied with
    # this solution...
    # -------------------------------------------------------------------------

    def fxn():
        warnings.warn("deprecated", DeprecationWarning)

    # -------------------------------------------------------------------------
    # Saving
    # -------------------------------------------------------------------------

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()

        message("Saving diagram to file : " + output.name)
        message("Be patient. This may be long for large datasets.")

        # NOTE : We must manually specify figure size with save_as_pdf_pages
        save_as_pdf_pages(filename=output.name,
                          plots=[my_plot + theme(figure_size=figsize)],
                          width=pdf_width,
                          height=pdf_height)


def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    ologram_merge_stats(**args)


if __name__ == '__main__':
    main()

else:

    # 'Bats' tests
    test = '''
    #ologram: get example files
    @test "ologram_merge_stats_0" {
         result=`gtftk get_example -d ologram_1 -f '*'`
      [ "$result" = "" ]
    }

    #ologram: produce heatmap
    @test "ologram_merge_stats_1" {
         result=`gtftk ologram_merge_stats H3K4me3_ologram_stats.tsv H3K36me3_ologram_stats.tsv H3K79me2_ologram_stats.tsv -o merged_ologram.pdf --labels H3K4me3,H3K36me3,H3K79me2`
      [ "$result" = "" ]
    }
    '''

    cmd = CmdObject(name="ologram_merge_stats",
                    message="Build a heatmap from several ologram output files (tsv).",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="ologram",
                    notes=__notes__,
                    updated=__updated__,
                    test=test)
