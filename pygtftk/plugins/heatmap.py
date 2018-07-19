#!/usr/bin/env python
from __future__ import division

import argparse
import os
import re
import warnings
import zipfile
from collections import OrderedDict

import numpy as np
import pandas as pd
import plotnine
from numpy import random
from pandas import Categorical
from plotnine import element_blank
from plotnine import element_line
from plotnine import element_rect
from plotnine import element_text
from plotnine import facet_grid
from plotnine import ggplot, aes, geom_tile
from plotnine import scale_fill_gradientn
from plotnine import scale_x_continuous
from plotnine import theme
from plotnine import theme_bw
from plotnine.labels import ggtitle
from scipy.cluster.vq import vq, kmeans
from scipy.stats import iqr

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import float_greater_than_null
from pygtftk.arg_formatter import float_grt_than_null_and_lwr_than_one
from pygtftk.arg_formatter import int_greater_than_null
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import chomp
from pygtftk.utils import mad
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import message
from pygtftk.utils import pos_max_val

__updated__ = "2018-01-20"
__doc__ = """
 Create a heatmap from mk_matrix result. This can be used to visualize coverage around or along a set of genomic features.
"""

__notes__ = """
    -- The program makes call to ggplot (geom_raster). It should be used with a limited number of features (e.g by selecting one transcript per gene). Otherwise it will run in memory issues.
    -- The program proposes various ways to split the dataset (kmeans, global kmeans, equally sized classes, chromosomes, user-defined transcript sets...).
    -- The program also proposes various way to organize the rows inside each sub-panel (mean, median, variance, standard deviation, median absolute deviation, inter-quartile range, maximum, minimum, position with maximum value, user_defined...).
    -- Kmeans and ordering are performed based on the first appearing bigwig. The leading bigwig can be changed using -bo.
"""


def make_parser():
    """The main parser."""

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help='A zipped of matrix file as produced by mk_matrix.',
                            default=None,
                            metavar='MATRIX',
                            type=FileWithExtension('r', '\.[Zz][Ii][Pp]'),
                            required=True)

    parser_grp.add_argument('-o', '--out-dir',
                            help='Output directory name.',
                            default="heatmap_gtftk",
                            type=str)

    parser_grp.add_argument('-t', '--transcript-file',
                            help="A two columns file with the transcripts"
                                 " of interest and their classes.",
                            default=None,
                            type=argparse.FileType("r"),
                            required=False)

    parser_grp.add_argument('-s', '--order-fun',
                            help="The statistics used for row ordering.",
                            default="mean",
                            choices=["mean", "median", "var",
                                     "sd", "mad", "IQR", "max", "min",
                                     "l_r", "r_l", "user_defined"],
                            type=str,
                            required=False)

    parser_grp.add_argument('-tl', '--to-log',
                            action="store_true",
                            help="Control whether the data should be log2-transform before plotting.",
                            required=False)

    parser_grp.add_argument('-rn', '--show-row-names',
                            action='store_true',
                            help='Show row names (need a limited set of genes).')

    """
    parser_grp.add_argument('-di', '--distance',
                            choices=("pearson", "euclidean", "maximum",
                                     "manhattan", "canberra",
                                     "binary", "abspearson",
                                     "abscorrelation", "correlation", "spearman", "kendall"),
                            default="pearson",
                            help='The distance to be used for k-means clustering.',
                            required=False)
    """

    parser_grp.add_argument('-bo', '--bwig-order-user',
                            help='A comma-separated list indicating bwig ordering.',
                            default=None,
                            required=False)

    parser_grp.add_argument('-y', '--y-factor',
                            help='The factor to use for y/second dimension.',
                            default='eq_sizes',
                            type=str,
                            choices=['kmeans',
                                     'gkmeans',
                                     'eq_sizes',
                                     'chrom',
                                     'tx_classes',
                                     'signal',
                                     'eq_sizes'],
                            required=False)

    parser_grp.add_argument('-c', '--color-palette',
                            type=str,
                            # default="#FFF7FB,#ECE2F0,#D0D1E6,#A6BDDB,#67A9CF,#3690C0,#02818A,#016450",
                            # default='#1A1835,#15464E,#2B6F39,#757B33,#C17A70,#D490C6,#C3C1F2,#CFEBEF',
                            # default="#0000FF,#00FFFF,#80FF80,#FFFF00,#FF0000",
                            # default="#0000AA,#0000FF,#00FFFF,#80FF80,#FFFF00,#FF0000,#AA0000",
                            # default='#4575b4,#74add1,#abd9e9,#e0f3f8,#fee090,#fdae61,#f46d43,#d73027',
                            # default='#67001f,#b2182b,#d6604d,#f4a582,#fddbc7,#f7f7f7,#d1e5f0,#92c5de,#4393c3,#2166ac,#053061',
                            # default="#2b83ba,#abdda4,#fdae61,#d7191c",
                            # default="#bababa,#f4a582,darkviolet",
                            # default="#0000BF,#0000FF,#0080FF,#00FFFF,#40FFBF,#80FF80,#BFFF40,#FFFF00,#FF8000,#FF0000,#BF0000", #matlab.like2
                            # default="#D58C52,#BF5D4E,#A92E4A,#930047", #
                            # jaime
                            default="#d73027,#fc8d59,#fee090,#e0f3f8,#91bfdb,#253494",
                            help='A set of colors to create an interpolated color palette.',
                            required=False)

    parser_grp.add_argument('-n', '--nb-class',
                            type=int_greater_than_null,
                            default=1,
                            help='Split the dataset into nb class based on mean expression level (exprs) or kmeans.',
                            required=False)

    parser_grp.add_argument('-pw', '--page-width',
                            help='Output pdf file width (inches).',
                            type=float_greater_than_null,
                            default=7,
                            required=False)

    parser_grp.add_argument('-ph', '--page-height',
                            help='Output odf file height (inches).',
                            type=int_greater_than_null,
                            default=7,
                            required=False)

    parser_grp.add_argument('-pf', '--page-format',
                            help='Output file format.',
                            choices=['pdf', 'png'],
                            default='pdf',
                            required=False)

    parser_grp.add_argument('-ti', '--title',
                            help='A title for the diagram.',
                            default="",
                            type=str,
                            required=False)

    parser_grp.add_argument('-xl', '--xlabel',
                            help='X axis label.',
                            default="Selected genomic regions",
                            type=str,
                            required=False)

    parser_grp.add_argument('-ml', '--max-line',
                            help='Add a line that underline the maximum values across rows (to be used with --order-fun 5p-3p).',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-fo', '--force-tx-class',
                            help='Force even if some transcripts from --transcript-file were not found.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-ms', '--min-signal',
                            help='All lines without a sum of bin values equal or greater to --min-signal will be deleted.',
                            type=float,
                            default=5,
                            required=False)

    parser_grp.add_argument('-ul',
                            '--upper-limit',
                            type=float_grt_than_null_and_lwr_than_one,
                            default=0.95,
                            help='Upper limit based on quantile computed from unique values.',
                            required=False)

    parser_grp.add_argument('-nm',
                            '--normalization-method',
                            choices=['none', 'pct'],
                            default='none',
                            help='The normalization method : pct = (x_i - min(x))/(max(x) - min(x)).',
                            required=False)

    parser_grp.add_argument('-if', '--user-img-file',
                            help="Provide an alternative path for the image.",
                            default=None,
                            type=argparse.FileType("w"),
                            required=False)

    parser_grp.add_argument('-ry', '--rotate-y-label',
                            help="Rotate the y label",
                            type=int,
                            default=0,
                            required=False)

    parser_grp.add_argument('-rx', '--rotate-x-label',
                            help="Rotate the x label",
                            type=int,
                            default=0,
                            required=False)

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=int_greater_than_null,
                            default=300,
                            required=False)

    return parser


def heatmap(inputfile=None,
            out_dir=None,
            transcript_file=None,
            order_fun=None,
            to_log=True,
            show_row_names=None,
            upper_limit=1,
            rotate_y_label=False,
            rotate_x_label=False,
            bwig_order_user=None,
            y_factor=None,
            nb_class=None,
            normalization_method='',
            color_palette=None,
            title='',
            page_width=7,
            min_signal=1,
            page_height=7,
            page_format='pdf',
            user_img_file=None,
            dpi=300,
            max_line=False,
            tmp_dir=None,
            xlabel="Selected genomic regions",
            logger_file=None,
            verbosity=False,
            force_tx_class=None
            ):
    # -------------------------------------------------------------------------
    #
    # Check args
    #
    # -------------------------------------------------------------------------

    if order_fun == "user_defined":
        if y_factor != 'tx_classes':
            message("if --order-fun is set to 'user_defined',"
                    " --y-factor should be set to 'tx_classes'.",
                    type="ERROR")

    if transcript_file is not None:
        if y_factor != 'tx_classes':
            message("if providing a transcript file"
                    " --y-factor should be set to 'tx_classes'.",
                    type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Input and output should not be the same (see yasmina issue)
    #
    # -------------------------------------------------------------------------

    if not inputfile.name.endswith('.zip'):
        message("Not a valid zip file (*.zip).", type="ERROR")

    base_input = os.path.split(os.path.abspath(inputfile.name))[1]
    base_input = base_input.replace(".zip", "")
    base_output = os.path.split(os.path.abspath(out_dir))[1]

    if base_output in [base_input, base_input]:
        message("The input file and output directory should have different names.",
                type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Unzipping input file
    #
    # -------------------------------------------------------------------------

    dir_name = os.path.dirname(os.path.abspath(inputfile.name))

    message("Uncompressing :" + dir_name,
            type="DEBUG")
    # input_zip = zipfile.ZipFile(inputfile.name, "r")
    # input_zip = input_zip.extractall()

    try:
        with zipfile.ZipFile(inputfile.name) as zf:
            zf.extractall(dir_name)
    except:
        message("Problem encountered when unzipping...",
                type="ERROR")
    inputfile_main = open(os.path.join(dir_name, zf.namelist()[0]), "r")
    message("Reading :" + inputfile_main.name,
            type="DEBUG")

    # -------------------------------------------------------------------------
    #
    # Retrieving info from the matrix file
    #
    # -------------------------------------------------------------------------

    message("Getting configuration info from input file.")

    input_file_tx = set()
    infile_chrom = set()
    infile_bwig = set()
    header = ""

    for line_number, line in enumerate(inputfile_main):

        # comment (line 0)
        if line_number == 0:
            header = chomp(line.lstrip("#"))
            header = header.rstrip(";")
            continue
        # skip header (line 1)
        elif line_number > 1:
            line = chomp(line)
            field = line.split("\t")
            tx_id = field[4]
            chrom = field[1]
            input_file_tx.add(tx_id)
            infile_chrom.add(chrom)
            infile_bwig.add(field[0])

    message("BigWigs found : " + ",".join(list(infile_bwig)))

    # -------------------------------------------------------------------------
    #
    # Parse the header
    #
    # -------------------------------------------------------------------------
    header = [x.split(":") for x in header.split(";")]

    config = dict()
    for x in header:
        config[x[0]] = x[1]

    # -------------------------------------------------------------------------
    #
    # Check arguments: --transcript-file
    #
    # -------------------------------------------------------------------------

    tx_to_class = OrderedDict()
    tx_class_list = OrderedDict()

    if y_factor == 'tx_classes':
        if transcript_file is not None:

            message("Reading transcript class file (--transcript_file).")

            for line in transcript_file:
                if line in ('\n', '\r\n'):
                    continue

                line = chomp(line)
                fields = line.split("\t")
                tx_name = fields[0].strip()

                try:
                    tx_to_class[tx_name] = chomp(fields[1])
                    tx_class_list[fields[1]] = 1
                except:
                    message("The file provided to --target-tx-file"
                            " should contain two columns (transcript and "
                            "class).", type="ERROR")

            nb_class = len(list(set(tx_class_list.keys())))

            if nb_class == 0:
                message("No transcript found in file provided through "
                        "--target-tx-file.", type="ERROR")

            else:
                message("Found : " + str(nb_class) + " classes.")

            # ------------------------------------------------------------------
            # Check the transcripts are found in the GTF...
            # ------------------------------------------------------------------

            intersect = set(tx_to_class.keys()).intersection(input_file_tx)
            nb_intersect = len(intersect)

            if nb_intersect != len(tx_to_class.keys()):
                not_found = list(set(tx_to_class.keys()) - input_file_tx)
                message(not_found[0] + "...", type="WARNING")

                if force_tx_class:
                    msg_type = "WARNING"
                else:
                    msg_type = "ERROR"

                message("Some transcripts (n={n}) from --transcript-class where not"
                        " found in input file (overriden with --force-tx-class).".format(
                    n=len(not_found)),
                    type=msg_type)

                if force_tx_class:
                    tx_to_class = {
                        k: tx_to_class[k] for k in tx_to_class if k in intersect}
                    nb_class = len(list(set(tx_to_class.values())))

                    if nb_class == 0:
                        message("No transcript found in file provided through "
                                "--target-tx-file.", type="ERROR")

            else:
                message("Found %d transcripts of interest in "
                        "input file." % nb_intersect)
        else:
            message(
                "Please provide --transcript-file if --x-factor or "
                "--y-factor is set to tx_classes.",
                type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Prepare output files
    #
    # -------------------------------------------------------------------------

    img_file = config['ft_type'] + "_u%s_d%s." + page_format
    img_file = img_file % (config['from'], config['to'])

    file_out_list = make_outdir_and_file(out_dir,
                                         [img_file,
                                          "data_long_format.txt",
                                          "transcript_order_and_class.txt"],
                                         force=True)

    img_file, data_file, tx_order_file_out = file_out_list

    if user_img_file is not None:
        os.unlink(img_file.name)
        img_file = user_img_file
        if not img_file.name.endswith(page_format):
            msg = "Image format: {f}. Please fix.".format(f=page_format)
            message(msg, type="ERROR")

        test_path = os.path.abspath(img_file.name)
        test_path = os.path.dirname(test_path)

        if not os.path.exists(test_path):
            os.makedirs(test_path)

    # -------------------------------------------------------------------------
    #
    # Colors for profiles
    #
    # -------------------------------------------------------------------------

    # Colors for the heatmap
    color_palette_list = color_palette.split(",")
    if len(color_palette_list) < 2:
        message("Need more than 2 colors for heatmap color palette.",
                type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Check bwig ordering
    #
    # -------------------------------------------------------------------------

    # bwig_order_user
    if bwig_order_user is None:
        bwig_order = list(infile_bwig)
    else:
        bwig_order = bwig_order_user.split(",")
        cond = [True for x in bwig_order if x not in infile_bwig]
        if any(cond):
            message("Fix --bwig-order. Unknown bwig. Should be one of: " + ",".join(infile_bwig),
                    type="ERROR")
        if len(set(bwig_order)) != len(bwig_order):
            message("Fix --bwig-order. Duplicates not allowed.",
                    type="ERROR")
        if any([x not in infile_bwig for x in bwig_order]):
            message("Fix --bwig-order. Some bigwigs were not found.",
                    type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Prepare diagram
    #
    # -------------------------------------------------------------------------

    first_bigwig = bwig_order[0]

    # -------------------------------------------------------------------------
    #
    # Read tab-delimited file
    #
    # -------------------------------------------------------------------------

    data = pd.read_csv(inputfile_main.name, sep="\t", header=1)

    # -------------------------------------------------------------------------
    #
    # Read transcript file
    #
    # -------------------------------------------------------------------------

    # all tx/gene names of the dataframe
    all_tx = list(OrderedDict.fromkeys(data['gene'].tolist()))

    if y_factor == 'tx_classes':

        # -------------------------------------------------------------------------
        # Get the transcript classes
        # -------------------------------------------------------------------------

        message("Reading transcript file.")
        df_classes = pd.read_csv(transcript_file, sep='\t', header=None)
        message("Deleting duplicates in transcript-file.")
        df_classes = df_classes.drop_duplicates(subset=[0])
        tx_ordering = df_classes[0].tolist()
        tx_classes = OrderedDict(zip(df_classes[0], df_classes[1]))

        # -------------------------------------------------------------------------
        # Select the transcript of interest and add classes info to the data.frame
        # -------------------------------------------------------------------------

        message("Checking how many genes where found in the transcript list.")

        nb_retained = len([x for x in all_tx if x in tx_ordering])

        msg = "Keeping {a} transcript out of {b}.".format(a=nb_retained, b=len(all_tx))
        message(msg)

        # subsetting
        data = data[[True if x in tx_ordering else False for x in data['gene'].tolist()]]
        data = data.assign(tx_classes=[tx_classes[x] for x in data['gene'].tolist()])

    else:
        data = data.assign(tx_classes=["1" for x in data['gene'].tolist()])

    # -------------------------------------------------------------------------
    #
    # Select the bwig of interest
    #
    # -------------------------------------------------------------------------

    data = data[data['bwig'].isin(bwig_order)]

    # -------------------------------------------------------------------------
    #
    # Find coverage columns
    #
    # -------------------------------------------------------------------------

    pos_order = []

    for i in data.columns:
        if re.search('(main_\d+)|(main_\d+)|(main_\d+)', i):
            pos_order += [i]

    bin_nb_main = len([x for x in data.columns if "main" in x])
    bin_nb_ups = len([x for x in data.columns if "upstream" in x])
    bin_nb_dws = len([x for x in data.columns if "downstream" in x])
    bin_nb_total = bin_nb_ups + bin_nb_main + bin_nb_dws

    # -------------------------------------------------------------------------
    #
    # Delete row without enough signal
    #
    # -------------------------------------------------------------------------

    # binding bwig data by column for testing the sum over rows
    bwig_nr = data['bwig'].unique()
    df_cbind = data.loc[data.bwig == bwig_nr[0], ['gene'] + pos_order]

    for k in bwig_nr[1:]:
        tmp = data.loc[data.bwig == k, ['gene'] + pos_order]
        df_cbind = df_cbind.merge(tmp, on='gene', how='left', suffixes=(k, 'previous'))

    df_cbind = df_cbind.set_index('gene')

    nb_gene_before = len(all_tx)
    msg = "Number of lines before filtering (--min-signal): {a}.".format(a=nb_gene_before)
    message(msg)

    df_cbind = df_cbind.assign(row_sum=df_cbind.sum(axis=1))
    df_cbind = df_cbind[df_cbind.row_sum >= min_signal].drop('row_sum', axis=1)
    all_tx = list(OrderedDict.fromkeys(df_cbind.index))

    msg = "Deleted lines (not enough signal, see -ms) : {a}.".format(a=nb_gene_before - len(all_tx))
    message(msg)
    msg = "Line kept: {a}.".format(a=len(all_tx))
    message(msg)

    # subset the main dataframe
    data = data.loc[[True if x in df_cbind.index else False for x in data.gene]]

    # -------------------------------------------------------------------------
    #
    # Melting (data -> dm)
    #
    # -------------------------------------------------------------------------

    dm = data.melt(id_vars=['tx_classes', 'bwig', 'chrom', 'start', 'end', 'gene', 'strand'], value_vars=pos_order)
    dm = dm.rename(columns={'variable': 'pos', 'value': 'exprs'})
    dm['bwig'] = Categorical(dm['bwig'], categories=bwig_order, ordered=True)

    # -------------------------------------------------------------------------
    #
    # ceiling
    #
    # -------------------------------------------------------------------------

    if upper_limit < 1:
        message('Ceiling')
        for k in dm['bwig'].unique():
            tmp = dm[dm.bwig == k]['exprs'].tolist()
            qu = np.percentile(list(set(tmp)), upper_limit * 100)
            tmp = [qu if x > qu else x for x in tmp]
            dm.loc[dm.bwig == k, 'exprs'] = tmp

    # -------------------------------------------------------------------------
    #
    # Normalize/transform
    #
    # -------------------------------------------------------------------------

    if to_log:
        if dm[dm.exprs == 0].shape[0] > 0:
            message("Zero value detected. Adding a pseudocount (+1) before log transformation.")
            dm['exprs'] = np.array(dm['exprs']) + 1

        message("Converting to log2.")
        dm['exprs'] = np.log2(dm['exprs'])
        ylab = "log2(Signal)"

    else:
        ylab = "Signal"

    if normalization_method == 'pct':
        message('Normalizing (ranging)')
        for k in dm['bwig'].unique():
            tmp = np.array(dm[dm.bwig == k]['exprs'].tolist())
            tmp_norm = (tmp - min(tmp)) / (max(tmp) - min(tmp)) * 100
            dm.loc[dm.bwig == k, 'exprs'] = tmp_norm

        ylab = "scaled(" + ylab + ", %)"

    # -------------------------------------------------------------------------
    #
    # Defining functions for row ordering
    #
    # -------------------------------------------------------------------------

    if order_fun == 'mean':
        fun_rows = np.nanmean
    elif order_fun == 'median':
        fun_rows = np.nanmedian
    elif order_fun == 'var':
        fun_rows = np.nanvar
    elif order_fun == 'sd':
        fun_rows = np.nanstd
    elif order_fun == 'mad':
        fun_rows = mad
    elif order_fun == 'IQR':
        fun_rows = iqr
    elif order_fun == 'mad':
        fun_rows = mad
    elif order_fun == 'min':
        fun_rows = np.nanmin
    elif order_fun == 'max':
        fun_rows = np.nanmax
    elif order_fun == 'l_r':
        fun_rows = pos_max_val
    elif order_fun == 'r_l':
        fun_rows = np.nanargmin

    # -------------------------------------------------------------------------
    #
    # Compute gene ordering
    #
    # -------------------------------------------------------------------------

    message("Computing gene ordering.")
    if order_fun == 'user_defined':
        message("Ordering based on --transcript-file.")
        dm['gene'] = Categorical(dm['gene'],
                                 categories=tx_ordering,
                                 ordered=True)
    else:

        message("Ordering based on : " + order_fun + ".")
        tmp = data.loc[data.bwig == first_bigwig, ['gene'] + pos_order]
        tmp = tmp.assign(fun_order=tmp[pos_order].apply(fun_rows, axis=1))
        tmp = tmp.sort_values('fun_order', ascending=False)
        dm['gene'] = Categorical(dm['gene'].tolist(),
                                 categories=tmp['gene'].tolist(),
                                 ordered=True)

    # -------------------------------------------------------------------------
    #
    # y-factor (tx/gene classes)
    #
    # -------------------------------------------------------------------------

    msg = "Computing gene classes n={a}.".format(a=nb_class)
    message(msg)

    random.seed((1000, 2000))

    if y_factor == 'gkmeans':
        msg = "Performing global kmeans  with : {a} classes.".format(a=nb_class)
        message(msg)

        # set row names and convert to np_array
        matrix = df_cbind.values

        # replace nans for 0 to avoid errors
        if np.any(np.isnan(matrix)):
            matrix[np.isnan(matrix)] = 0

        centers = kmeans(matrix, nb_class)[0]

        # map objects to the centroids
        cluster_labels = vq(matrix, centers)[0]
        km_dict = OrderedDict(zip(df_cbind.index, cluster_labels))
        dm = dm.assign(y_factor=[km_dict[x] for x in dm.gene])

        message("Nb elements per class : " + str(Categorical(cluster_labels).value_counts().tolist()))


    elif y_factor == 'kmeans':

        msg = "Performing kmeans (based on  {a} signal) with : {b} classes."
        msg = msg.format(a=first_bigwig, b=nb_class)
        message(msg)

        # Kmeans is performed on the first bigwig
        df_first_bw = data.loc[data.bwig == first_bigwig, ['gene'] + pos_order]

        # set row names and convert to np_array
        df_first_bw = df_first_bw.set_index('gene')
        matrix = df_first_bw.values

        # replace nans for 0 to avoid errors
        if np.any(np.isnan(matrix)):
            matrix[np.isnan(matrix)] = 0

        centers = kmeans(matrix, nb_class)[0]

        # map objects to the centroids
        cluster_labels = vq(matrix, centers)[0]
        km_dict = OrderedDict(zip(df_first_bw.index, cluster_labels))
        dm = dm.assign(y_factor=[km_dict[x] for x in dm.gene])

        message("Nb elements per class : " + str(Categorical(cluster_labels).value_counts().tolist()))

    elif y_factor == 'signal':

        if nb_class > 1:
            tmp = data.loc[data.bwig == first_bigwig, ['gene'] + pos_order]
            signal = tmp[pos_order].apply(fun_rows, axis=1)
            row_classes = pd.cut(signal, nb_class)
            gene2classes = OrderedDict(zip(tmp.index, row_classes))
            dm = dm.assign(y_factor=[gene2classes[x] for x in dm.gene])

        message("Nb elements per class : " + str(row_classes.value_counts().tolist()))

    # by eq_sizes

    elif y_factor == 'eq_sizes':

        if nb_class > 1:

            tmp = data.loc[data.bwig == first_bigwig, ['gene'] + pos_order]
            signal = tmp[pos_order].apply(fun_rows, axis=1).tolist()

            classes = np.array_split(np.argsort(signal), nb_class)

            k = dict()

            for i in range(nb_class):
                for j in classes[i]:
                    k[tmp.iloc[j].loc['gene']] = i

            dm = dm.assign(y_factor=[k[x] for x in dm.gene])

        else:
            dm = dm.assign(y_factor=['1' for x in dm.gene])

        message("Nb elements per class : " + str(Categorical([x[1] for x in k.items()]).value_counts().tolist()))

    dm.y_factor = Categorical(dm.y_factor.tolist())

    # -------------------------------------------------------------------------
    #
    # Column ordering
    #
    # -------------------------------------------------------------------------

    message("Computing column ordering.")
    pos_order = sorted(sorted(dm.pos.unique(), reverse=True), key=lambda x: int(x.split("_")[1]))
    dm.pos = Categorical(dm.pos.tolist(), categories=pos_order, ordered=True)

    # -------------------------------------------------------------------------
    #
    # Prepare x axis. Turning x axis into continuous scale if needed
    #
    # -------------------------------------------------------------------------

    fr = int(config['from'])
    to = int(config['to'])

    if config['ft_type'] in ["transcript", "user_regions"]:

        dm = dm.assign(extra=dm.pos.cat.codes)
        seq = np.linspace(0, 100, len(dm.extra.unique()))
        dm.extra = seq[dm.extra]
        dm.extra = [float(x) for x in dm.extra]
        dm.pos = dm.extra


    else:
        dm = dm.assign(extra=dm.pos.cat.codes)
        seq = np.linspace(-fr, to, len(dm.pos.unique()))
        dm.extra = seq[dm.extra]
        dm.extra = [float(x) for x in dm.extra]
        dm.pos = dm.extra

    dm = dm.drop('extra', axis=1)

    # -------------------------------------------------------------------------
    #
    # Preparing diagram
    #
    # -------------------------------------------------------------------------

    message("Preparing diagram")

    p = ggplot(data=dm, mapping=aes('pos',
                                    'gene')) + geom_tile(aes(fill='exprs'))

    p += theme_bw()
    p += theme(legend_text=element_text(size=6),
               panel_grid_major=element_blank(),
               panel_grid_minor=element_blank(),
               panel_border=element_rect(colour="black", size=1),
               legend_key_size=2,
               legend_position="top",
               legend_key=element_rect(colour="white"),
               axis_text_y=element_text(colour="#333333",
                                        size=4),
               axis_line=element_line(size=0.1, colour=""),
               axis_text_x=element_text(colour="#333333",
                                        size=6,
                                        angle=65,
                                        face="plain"),
               axis_title_x=element_text(colour="#333333",
                                         size=8,
                                         angle=0,
                                         face="plain"),
               axis_title_y=element_text(colour="#333333",
                                         size=8, angle=90,
                                         face="plain"),
               strip_text_y=element_text(size=7,
                                         colour='white',
                                         angle=rotate_y_label),
               strip_text_x=element_text(colour='white',
                                         size=7,
                                         angle=rotate_x_label),
               strip_background=element_rect(colour="#000000",
                                             fill="#000000")
               )

    p += scale_fill_gradientn(colors=color_palette_list,
                              name="Signal", na_value="#222222")

    p += plotnine.labels.xlab(xlabel)
    p += plotnine.labels.ylab("Genes")

    if not show_row_names:
        p += theme(axis_text_y=element_blank(),
                   axis_ticks_major_y=element_blank(),
                   axis_ticks_minor_y=element_blank())

    if config['ft_type'] in ["transcript", "user_regions"]:

        if config['from']:

            if config['to']:
                ticks = [0, bin_nb_ups / 2] + list(np.linspace(bin_nb_ups, bin_nb_main + bin_nb_ups, 11)) + [
                    bin_nb_total - bin_nb_dws / 2, bin_nb_total]
                ticks = [x / bin_nb_total * 100 for x in ticks]
                labels = [-fr,
                          round(-fr / 2, 0)
                          ] + [str(x) + "%" for x in np.linspace(0, 100, 11)] + [round(to / 2, 0), to]


            else:
                ticks = [0, bin_nb_ups / 2] + list(np.linspace(bin_nb_ups, bin_nb_total, 11))
                ticks = [x / bin_nb_total * 100 for x in ticks]

                labels = [- fr,
                          round(-fr / 2, 0)
                          ] + [str(x) + "%" for x in np.linspace(0, 100, 6)]
        else:
            if config['to']:

                ticks = list(np.linspace(0, bin_nb_main, 11)) + [bin_nb_total - bin_nb_dws / 2,
                                                                 bin_nb_total,
                                                                 bin_nb_total - bin_nb_dws / 2,
                                                                 bin_nb_total]
                ticks = [x / bin_nb_total * 100 for x in ticks]

                labels = [str(x) + "%" for x in np.linspace(0, 100, 6)] + [to / 2, to]


            else:
                ticks = list(np.linspace(0, bin_nb_total, 6))
                ticks = [x / bin_nb_total * 100 for x in ticks]
                labels = [str(x) + "%" for x in np.linspace(0, 100, 6)]

        p += scale_x_continuous(expand=[0, 0], breaks=ticks, labels=labels)

    else:
        p += scale_x_continuous(expand=[0, 0])

    p += facet_grid("y_factor ~ bwig ", scales="free_y", space="free_y")

    p += ggtitle(title)

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
        p.save(filename=img_file.name, width=page_width, height=page_height, dpi=dpi)
        dm.to_csv(data_file, sep="\t", header=True, index=False)


if __name__ == '__main__':

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    heatmap(**args)

else:

    test = '''

    #heatmap: prepare dataset
    @test "heatmap_1" {
     result=`gtftk get_example -d mini_real -f '*'; gtftk overlapping -i mini_real.gtf.gz -c hg38.genome  -n > mini_real_noov.gtf; gtftk random_tx -i mini_real_noov.gtf  -m 1 -s 123 > mini_real_noov_rnd_tx.gtf`
      [ -s "hg38.genome" ]
    }

    #heatmap: make mini_real_promoter
    @test "heatmap_2" {
     result=`gtftk mk_matrix -i mini_real_noov_rnd_tx.gtf -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_promoter`
      [ -s "mini_real_promoter.zip" ]
    }

    #heatmap: make mini_real_promoter
    @test "heatmap_3" {
     result=`describe()`
      [ -s "example_10.png" ]
    }

    #heatmap: make mini_real_promoter
    @test "heatmap_4" {
     result=`gtftk heatmap -D -i mini_real_promoter.zip -o heatmap_prom_2  -tl  -n 5 --y-factor kmeans -ul 0.9 -nm pct -pf png -if example_11.png -c "#e66101,#fdb863,#f7f7f7,#b2abd2,#5e3c99"`
      [ -s "example_11.png" ]
    }
        
    #heatmap: make mini_real_promoter
    @test "heatmap_5" {
     result=`gtftk heatmap -D -i mini_real_promoter.zip -o heatmap_prom_3  -tl  -n 5 --bwig-order-user  H3K36me3,H3K4me3,H3K79me --y-factor kmeans -ul 0.75 -nm pct -pf png -if example_12.png`
      [ -s "example_12.png" ]
    }
    
    #heatmap: make mini_real_promoter
    @test "heatmap_6" {
     result=`gtftk heatmap -D -i mini_real_promoter.zip -o heatmap_prom_4 -n 5 -tl --bwig-order-user  H3K36me3,H3K4me3,H3K79me --y-factor gkmeans -s max -ul 0.75 -nm pct -pf png -if example_13.png`
      [ -s "example_13.png" ]
    }
        
    
    #heatmap: make mini_real_promoter
    @test "heatmap_7" {
     result=`gtftk heatmap -D -i mini_real_promoter.zip -o heatmap_prom_5 -n 5 -tl --y-factor eq_sizes  -s mean -ul 0.75 -nm pct -pf png -if example_14.png -c "#0000AA,#0055FF,#00AAFF,#40FFFF,#80FFBF,#BFFF80,#FFFF40,#FFAA00,#FF5500,#AA0000"`
      [ -s "example_14.png" ]
    }

    #heatmap: make mini_real_promoter
    @test "heatmap_8" {
     result=`gtftk heatmap -D -i mini_real_promoter.zip -o heatmap_prom_6 -n 5 -tl --y-factor eq_sizes  -s l_r -ul 0.75 -nm pct -pf png -if example_15.png`
      [ -s "example_15.png" ]
    }
    

    #heatmap: make mini_real_promoter
    @test "heatmap_9" {
     result=`gtftk tabulate -i mini_real.gtf.gz -k transcript_id,transcript_biotype -Hun | perl -ne  'print if(/(protein_coding)|(lincRNA)/)'  > tx_classes.txt; head -n 100 tx_classes.txt >  tx_classes_100.txt`
      [ -s "tx_classes_100.txt" ]
    }
    
        
    #heatmap: make mini_real_promoter
    @test "heatmap_10" {
     result=` gtftk heatmap  -i mini_real_promoter.zip  -t tx_classes.txt -y tx_classes  -tl  -fo -o hh -s user_defined -c "#F9DA6B,#f03b20"  -pf png -if example_16.png`
      [ -s "example_16.png" ]
    }

    #heatmap: make mini_real_promoter
    @test "heatmap_11" {
     result=`gtftk heatmap  -i mini_real_promoter.zip  -t tx_classes_100.txt -y tx_classes -tl -fo -o hh -s user_defined --show-row-names -c "#c51b7d,#e9a3c9,#fde0ef,#f7f7f7,#e6f5d0,#a1d76a,#4d9221"  -pf png -if example_17.png`
      [ -s "example_17.png" ]
    }
        
    '''

    cmd = CmdObject(name="heatmap",
                    message="Create a heatmap from mk_matrix result.",
                    parser=make_parser(),
                    fun=heatmap,
                    desc=__doc__,
                    updated=__updated__,
                    notes=__notes__,
                    group="coverage",
                    test=test)
