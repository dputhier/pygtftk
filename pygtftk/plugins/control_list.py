"""
    The ``control_list`` plugin
    ============================

    Based on a reference gene list, returns a list of genes matched for
    signal. The expression/signal values should be provided for all genes through
    the in_file argument.
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from pandas import Categorical
from plotnine import (aes, xlab,
                      ylab, geom_jitter,
                      geom_rug, facet_wrap,
                      theme, element_blank,
                      theme_bw, scale_fill_manual, geom_violin)
from plotnine import ggplot
from plotnine.exceptions import PlotnineError

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import chomp
from pygtftk.utils import is_hex_color
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import message

__updated__ = "2018-01-20"
__doc__ = """
 Based on a reference gene list (or more generally IDs) this command tries to extract a set of
 other genes/IDs matched for signal/expression. The --reference-gene-file contains
 the list of reference IDs while the -\inputfile contains a tuple gene/signal for all genes.
"""

__notes__ = """
    -- -\-infile is a two columns tabulated file. The first column contains the list of ids (including reference IDs)
    and the second column contains the expression/signal values. This file should contain no header.
    -- Think about discarding any unwanted IDs from -\-infile before calling control_list.
"""


def make_parser():
    """The program parser."""

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('--in-file', '-i',
                            metavar='TXT',
                            help='A two columns tab-file. See notes.',
                            default=None,
                            type=arg_formatter.FormattedFile(mode='r', file_ext='txt'),
                            required=True)

    parser_grp.add_argument('--reference-gene-file', '-r',
                            metavar='TXT',
                            help='The file containing the reference gene list (1 column,'
                                 ' transcript ids).'
                                 ' No header.',
                            default=None,
                            type=arg_formatter.FormattedFile(mode='r', file_ext='txt'),
                            required=True)

    parser_grp.add_argument('--out-dir', '-o',
                            help='Name of the output directory.',
                            type=str,
                            metavar='DIR',
                            default="control_list",
                            required=False)

    parser_grp.add_argument('--log2', '-l',
                            help='If selected, data will be log transformed.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('--pseudo-count', '-p',
                            help='The value for a pseudo-count to be added.',
                            type=float,
                            default=0,
                            required=False)

    parser_grp.add_argument('-pw', '--page-width',
                            help='Output pdf file width (e.g. 7 inches).',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=None,
                                                          linc=False,
                                                          val_type='float'),
                            default=None,
                            required=False)

    parser_grp.add_argument('-ph', '--page-height',
                            help='Output  file height (e.g. 5 inches).',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=None,
                                                          linc=False,
                                                          val_type='float'),
                            default=None,
                            required=False)

    parser_grp.add_argument('-pf', '--page-format',
                            help='Output file format.',
                            choices=['pdf', 'png'],
                            default='pdf',
                            required=False)

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=None,
                                                          linc=False,
                                                          val_type='int'),
                            default=300,
                            required=False)

    parser_grp.add_argument('--skip-first', '-s',
                            help='Indicates that infile hase a header.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('--rug', '-u',
                            help='Add rugs to the diagram.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('--jitter', '-j',
                            help='Add jittered points.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-if', '--user-img-file',
                            help="Provide an alternative path for the image.",
                            default=None,
                            type=argparse.FileType("w"),
                            required=False)

    parser_grp.add_argument('-c', '--set-colors',
                            help='Colors for the two sets (comma-separated).',
                            default="#b2df8a,#6a3d9a",
                            type=str,
                            required=False)

    return parser_grp


def control_list(in_file=None,
                 out_dir=None,
                 reference_gene_file=None,
                 log2=False,
                 page_width=None,
                 page_height=None,
                 user_img_file=None,
                 page_format=None,
                 pseudo_count=1,
                 set_colors=None,
                 dpi=300,
                 rug=False,
                 jitter=False,
                 skip_first=False):
    # -------------------------------------------------------------------------
    #
    # Check in_file content
    #
    # -------------------------------------------------------------------------

    for p, line in enumerate(in_file):

        line = chomp(line)
        line = line.split("\t")

        if len(line) > 2:
            message("Need a two columns file.",
                    type="ERROR")
        if skip_first:
            if p == 0:
                continue
        try:
            fl = float(line[1])
        except:
            msg = "It seems that column 2 of input file"
            msg += " contains non numeric values. "
            msg += "Check that no header is present and that "
            msg += "columns are ordered properly. "
            msg += "Or use '--skip-first'. "
            message(msg, type="ERROR")

        if log2:
            fl = fl + pseudo_count
            if fl <= 0:
                message("Can not log transform negative/zero values. Add a pseudo-count.",
                        type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Check colors
    #
    # -------------------------------------------------------------------------

    set_colors = set_colors.split(",")

    if len(set_colors) != 2:
        message("Need two colors. Please fix.", type="ERROR")

    mcolors_name = mcolors.cnames

    for i in set_colors:
        if i not in mcolors_name:
            if not is_hex_color(i):
                message(i + " is not a valid color. Please fix.", type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Preparing output files
    #
    # -------------------------------------------------------------------------

    # Preparing pdf file name
    file_out_list = make_outdir_and_file(out_dir, ["control_list.txt",
                                                   "reference_list.txt",
                                                   "diagnostic_diagrams." + page_format],
                                         force=True)

    control_file, reference_file_out, img_file = file_out_list

    if user_img_file is not None:

        os.unlink(img_file.name)
        img_file = user_img_file

        if not img_file.name.endswith(page_format):
            msg = "Image format should be: {f}. Please fix.".format(f=page_format)
            message(msg, type="ERROR")

        test_path = os.path.abspath(img_file.name)
        test_path = os.path.dirname(test_path)

        if not os.path.exists(test_path):
            os.makedirs(test_path)

    # -------------------------------------------------------------------------
    #
    # Read the reference list
    #
    # -------------------------------------------------------------------------

    try:
        reference_genes = pd.read_csv(reference_gene_file.name, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        message("No genes in --reference-gene-file.", type="ERROR")

    reference_genes.rename(columns={reference_genes.columns.values[0]: 'gene'}, inplace=True)

    # -------------------------------------------------------------------------
    #
    # Delete duplicates
    #
    # -------------------------------------------------------------------------

    before = len(reference_genes)
    reference_genes = reference_genes.drop_duplicates(['gene'])
    after = len(reference_genes)

    msg = "%d duplicate lines have been deleted in reference file."
    message(msg % (before - after))

    # -------------------------------------------------------------------------
    #
    # Read expression data and add the pseudo_count
    #
    # -------------------------------------------------------------------------

    if skip_first:
        exp_data = pd.read_csv(in_file.name, sep="\t",
                               header=None, index_col=None,
                               skiprows=[0], names=['exprs'])
    else:

        exp_data = pd.read_csv(in_file.name, sep="\t", names=['exprs'], index_col=0)

    exp_data.exprs = exp_data.exprs.values + pseudo_count

    # -------------------------------------------------------------------------
    #
    # log transformation
    #
    # -------------------------------------------------------------------------

    ylabel = 'Expression'

    if log2:
        if len(exp_data.exprs.values[exp_data.exprs.values == 0]):
            message("Can't use log transformation on zero or negative values. Use -p.",
                    type="ERROR")
        else:
            exp_data.exprs = np.log2(exp_data.exprs.values)
            ylabel = 'log2(Expression)'

    # -------------------------------------------------------------------------
    #
    # Are reference gene found in control list
    #
    # -------------------------------------------------------------------------

    # Sort in increasing order
    exp_data = exp_data.sort_values('exprs')

    #  Vector with positions indicating which in the
    # expression data list are found in reference_gene

    reference_genes_found = [x for x in reference_genes['gene'] if x in exp_data.index]

    msg = "Found %d genes of the reference in the provided signal file" % len(reference_genes_found)
    message(msg)

    not_found = [x for x in reference_genes['gene'] if x not in exp_data.index]

    if len(not_found):
        if len(not_found) == len(reference_genes):
            message("Genes from reference file where not found in signal file (n=%d)." % len(not_found), type="ERROR")
        else:
            message("List of reference genes not found :%s" % not_found)
    else:
        message("All reference genes were found.")

    # -------------------------------------------------------------------------
    #
    # Search for genes with matched signal
    #
    # -------------------------------------------------------------------------

    exp_data_save = exp_data.copy()

    control_list = list()

    nb_candidate_left = exp_data.shape[0] - len(reference_genes_found)

    message("Searching for genes with matched signal.")

    if nb_candidate_left < len(reference_genes_found):
        message("Not enough element to perform selection. Exiting", type="ERROR")

    for i in reference_genes_found:
        not_candidates = reference_genes_found + control_list
        not_candidates = list(set(not_candidates))

        diff = abs(exp_data.loc[i] - exp_data)
        control_list.extend(diff.loc[np.setdiff1d(diff.index, not_candidates)].idxmin(axis=0, skipna=True).tolist())

    # -------------------------------------------------------------------------
    #
    # Prepare a dataframe for plotting
    #
    # -------------------------------------------------------------------------

    message("Preparing a dataframe for plotting.")

    reference = exp_data_save.loc[reference_genes_found].sort_values('exprs')
    reference = reference.assign(genesets=['Reference'] * reference.shape[0])

    control = exp_data_save.loc[control_list].sort_values('exprs')
    control = control.assign(genesets=['Control'] * control.shape[0])

    data = pd.concat([reference, control])
    data['sets'] = pd.Series(['sets' for x in data.index.tolist()], index=data.index)
    data['genesets'] = Categorical(data['genesets'])

    # -------------------------------------------------------------------------
    #
    # Diagnostic plots
    #
    # -------------------------------------------------------------------------

    p = ggplot(data, aes(x='sets', y='exprs', fill='genesets'))

    p += scale_fill_manual(values=dict(zip(['Reference', 'Control'], set_colors)))

    p += geom_violin(color=None)

    p += xlab('Gene sets') + ylab(ylabel)

    p += facet_wrap('~genesets')

    if rug:
        p += geom_rug()

    if jitter:
        p += geom_jitter()

    p += theme_bw()
    p += theme(axis_text_x=element_blank())

    # -------------------------------------------------------------------------
    # Turn warning off. Both pandas and plotnine use warnings for deprecated
    # functions. I need to turn they off although I'm not really satisfied with
    # this solution...
    # -------------------------------------------------------------------------

    def fxn():
        warnings.warn("deprecated", DeprecationWarning)

    # -------------------------------------------------------------------------
    #
    # Saving
    #
    # -------------------------------------------------------------------------

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()
        message("Saving diagram to file : " + img_file.name)
        message("Be patient. This may be long for large datasets.")

        try:
            p.save(filename=img_file.name, width=page_width, height=page_height, dpi=dpi, limitsize=False)
        except PlotnineError as err:
            message("Plotnine message: " + err.message)
            message("Plotnine encountered an error.", type="ERROR")

    # -------------------------------------------------------------------------
    #
    # write results
    #
    # -------------------------------------------------------------------------

    exp_data_save.loc[reference_genes_found].sort_values('exprs').to_csv(reference_file_out.name, sep="\t")
    exp_data_save.loc[control_list].sort_values('exprs').to_csv(control_file.name, sep="\t")


if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    control_list(**args)

else:

    test = '''
    #control_list: load dataset
    @test "control_list_0" {
     result=`gtftk get_example -f '*' -d control_list`
      [ "$result" = "" ]
    }
    
    #control_list
    @test "control_list_1" {
      result=`gtftk control_list -i control_list_data.txt -r control_list_reference.txt -D ; cat control_list/control_list.txt | cut -f2| perl -npe 's/\\n/,/'`
      [ "$result" = "exprs,2.02,4.04,6.06," ]
    }
        
    '''

    cmd = CmdObject(name="control_list",
                    message="Returns a list of gene matched for expression based on reference values.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    updated=__updated__,
                    notes=__notes__,
                    group="miscellaneous",
                    test=test)
