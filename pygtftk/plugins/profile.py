#!/usr/bin/env python
from __future__ import division

import argparse
import os
import re
import shutil
import tempfile
import warnings
import zipfile
from collections import OrderedDict

import matplotlib as mpl
import numpy as np
import pandas as pd
import plotnine
from matplotlib import cm
from pandas import Categorical
from plotnine import aes
from plotnine import element_blank
from plotnine import element_rect
from plotnine import element_text
from plotnine import facet_wrap
from plotnine import geom_line
from plotnine import geom_rect
from plotnine import ggplot
from plotnine import ggtitle
from plotnine import guide_legend
from plotnine import guides
from plotnine import scale_x_continuous
from plotnine import theme
from plotnine import ylab

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import float_greater_than_null
from pygtftk.arg_formatter import float_grt_than_null_and_lwr_than_one
from pygtftk.arg_formatter import int_greater_than_null
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import chomp
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

R_LIB = 'ggplot2,reshape2,grid,data.table,plyr'

__updated__ = "2018-01-20"
__doc__ = """
 Produces bigWig coverage profiles using calls to ggplot2.
"""

__notes__ = """
 -- The ranging normalization method [1] implies the following transformation: 
 -- -  (x_i - min(x))/(max(x) - min(x)).
 -- Think about using normalized bigWig files as input to mk_matrix. This
 will limit the requirement for an additional normalization step (see
 Deeptools for a set of useful methods implemented in bamCoverage/bamCompare).
"""

__references__ = """
 -- [1] Numerical Ecology - second Edition - P. Legendre, L. Legendre (1998) Elsevier.
"""


def make_parser():
    """The main parser."""

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help='A zip file containing a matrix as produced by mk_matrix.',
                            default=None,
                            metavar='MATRIX',
                            type=FileWithExtension('r', '\.[Zz][Ii][Pp]'),
                            required=True)

    parser_grp.add_argument('-o', '--out-dir',
                            help='Output directory name.',
                            default="draw_profile",
                            metavar="DIR",
                            type=str)

    parser_grp.add_argument('-t', '--transcript-file',
                            help="A two columns file with the transcripts"
                                 " of interest and their classes.",
                            default=None,
                            type=argparse.FileType("r"),
                            required=False)

    parser_grp.add_argument('-s', '--stat',
                            help="The statistics to be computed.",
                            default="mean",
                            choices=["mean", "median", "sum", "min", "max"],
                            type=str,
                            required=False)

    parser_grp.add_argument('-c', '--profile-colors',
                            help='Colors.',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-d', '--color-order',
                            help='Factor ordering. Comma separated bwig labels or tx classes.',
                            default=None,
                            type=str)

    parser_grp.add_argument('-g', '--group-by',
                            help='The variable used for grouping.',
                            default='bwig',
                            type=str,
                            choices=['bwig', 'tx_classes', 'chrom'],
                            required=False)

    parser_grp.add_argument('-f', '--facet-var',
                            help='The variable to be used for splitting into facets.',
                            default=None,
                            type=str,
                            choices=['bwig', 'tx_classes', 'chrom'],
                            required=False)

    parser_grp.add_argument('-pw', '--page-width',
                            help='Output pdf file width (inches).',
                            type=int_greater_than_null,
                            default=7,
                            required=False)

    parser_grp.add_argument('-ph', '--page-height',
                            help='Output odf file height (inches).',
                            type=int_greater_than_null,
                            default=5,
                            required=False)

    parser_grp.add_argument('-pf', '--page-format',
                            help='Output file format.',
                            choices=['pdf', 'png'],
                            default='pdf',
                            required=False)

    parser_grp.add_argument('-lw', '--line-width',
                            help='Line width.',
                            type=float_greater_than_null,
                            default=1.25,
                            required=False)

    parser_grp.add_argument('-sc', '--strip-color',
                            help='Strip background color.',
                            default="#989898",
                            type=str,
                            required=False)

    parser_grp.add_argument('-x', '--xlab',
                            help='X axis label.',
                            default="\nSelected genomic regions",
                            type=str,
                            required=False)

    parser_grp.add_argument('-at', '--axis-text',
                            help='Size of axis text.',
                            default=8,
                            type=int,
                            required=False)

    parser_grp.add_argument('-st', '--strip-text',
                            help='Size of strip text.',
                            default=8,
                            type=int,
                            required=False)

    parser_grp.add_argument('-fc', '--facet-col',
                            help='Number of facet columns.',
                            default=None,
                            type=int,
                            required=False)

    parser_grp.add_argument('-fo', '--force-tx-class',
                            help='Force even if some transcripts from --transcript-file were not found.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-if', '--user-img-file',
                            help="Provide an alternative path for the image.",
                            default=None,
                            type=argparse.FileType("w"),
                            required=False)

    parser_grp.add_argument('-ul',
                            '--upper-limit',
                            type=float_grt_than_null_and_lwr_than_one,
                            default=0.95,
                            help='Upper limit based on quantile computed from unique values.',
                            required=False)

    parser_grp.add_argument('-nm',
                            '--normalization-method',
                            choices=['none', 'ranging'],
                            default='none',
                            help='The normalization method performed on a per bigwig basis.',
                            required=False)

    parser_grp.add_argument('-tl', '--to-log',
                            action="store_true",
                            help="Control whether the data should be log2-transform before plotting.",
                            required=False)

    parser_grp.add_argument('-ti', '--title',
                            help='A title for the diagram.',
                            default="",
                            type=str,
                            required=False)

    parser_grp.add_argument('-dpi', '--dpi',
                            help='Dpi to use.',
                            type=int_greater_than_null,
                            default=300,
                            required=False)

    parser_grp.add_argument('-th', '--theme-plotnine',
                            choices=[
                                '538',
                                'bw',
                                'grey',
                                'gray',
                                'linedraw',
                                'light',
                                'dark',
                                'minimal',
                                'classic',
                                'void',
                                'test',
                                'matplotlib',
                                'seaborn',
                                'xkcd'],
                            default='bw',
                            help="The theme for ggplot2 diagram.",
                            required=False)

    return parser


def draw_profile(inputfile=None,
                 out_dir=None,
                 group_by='bwig',
                 color_order=None,
                 transcript_file=None,
                 transform=None,
                 normalization_method=None,
                 to_log=False,
                 upper_limit=0.95,
                 quantiles=False,
                 profile_colors=None,
                 page_width=None,
                 title=None,
                 page_height=None,
                 page_format='pdf',
                 user_img_file=None,
                 tmp_dir=None,
                 facet_col=None,
                 strip_color="#707070",
                 force_tx_class=False,
                 stat="mean",
                 facet_var=None,
                 xlab="Selected genomic regions",
                 axis_text=8,
                 strip_text=8,
                 line_width=1,
                 theme_plotnine='bw',
                 dpi=300,
                 logger_file=None,
                 verbosity=False
                 ):
    # -------------------------------------------------------------------------
    #
    # The selected theme
    #
    # -------------------------------------------------------------------------

    theme_plotnine = 'theme_' + theme_plotnine

    # -------------------------------------------------------------------------
    #
    # The selected stat
    #
    # -------------------------------------------------------------------------

    if stat == "mean":
        stat_fun = np.mean
    elif stat == "median":
        stat_fun = np.median
    elif stat == "min":
        stat_fun = np.min
    elif stat == "max":
        stat_fun = np.max
    elif stat == "sum":
        stat_fun = np.sum

    # -------------------------------------------------------------------------
    #
    # facet and group_by should be different
    #
    # -------------------------------------------------------------------------

    if facet_var == group_by:
        message("--facet-var and --group-by should be different.",
                type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Check argument consistency
    #
    # -------------------------------------------------------------------------

    if facet_var == 'tx_classes' or group_by == 'tx_classes':
        if transcript_file is None:
            message("Please provide --transcript-file",
                    type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Input and output should not be the same
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

    # Use a temp file to avoid concurrency issues
    dir_name = tempfile.mkdtemp(prefix='GTFtk_matrix_')
    message("Uncompressing : " + dir_name,
            type="DEBUG")

    try:
        with zipfile.ZipFile(inputfile.name) as zf:
            zf.extractall(dir_name)
    except:
        message("Problem encountered when unzipping...",
                type="ERROR")

    inputfile_main = open(os.path.join(dir_name, zf.namelist()[0]), "r")
    message("Reading : " + inputfile_main.name,
            type="DEBUG")

    # -------------------------------------------------------------------------
    #
    # Retrieving info from the matrix file
    #
    # -------------------------------------------------------------------------

    message("Getting configuration info from input file.")

    input_file_tx = set()
    input_file_chrom = set()
    input_file_bwig = set()
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
            input_file_chrom.add(chrom)
            input_file_bwig.add(field[0])

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
    # Check arguments: --facet, --group-by according to the number of bigwig
    #
    # -------------------------------------------------------------------------
    # If one is analyzing more than one bigwig
    # the "bwig" factor should appear in  --facet or --group-by
    if len(input_file_bwig) > 1:
        if facet_var != "bwig":
            if group_by != "bwig":
                message("If more than on bigWig is analyzed, --facet or --group-by should be set to 'bwig'.",
                        type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Check arguments: --transcript-class
    #
    # -------------------------------------------------------------------------

    tx_to_class = dict()
    class_list = set()

    if transcript_file is not None:

        message("Reading transcript class file (--transcript_file).")

        for line in transcript_file:
            if line in ('\n', '\r\n'):
                continue

            line = chomp(line)
            fields = line.split("\t")
            tx_name = fields[0]

            try:
                tx_to_class[tx_name] = chomp(fields[1])
                class_list.add(fields[1])
            except:
                message("The file provided to --target-tx-file"
                        " should contain two columns (transcript and "
                        "class).", type="ERROR")

        nb_tx_classes = len(list(set(tx_to_class.values())))

        if nb_tx_classes == 0:
            message("No transcript found in file provided through "
                    "--target-tx-file.", type="ERROR")

        else:
            message("Found : " + str(nb_tx_classes) + " classes.")

        # ----------------------------------------------------------------------
        # Check the transcripts are found in the GTF...
        # ----------------------------------------------------------------------

        intersect = set(tx_to_class.keys()).intersection(input_file_tx)
        nb_intersect = len(intersect)

        if nb_intersect != len(tx_to_class.keys()):
            not_found = list(set(tx_to_class.keys()) - input_file_tx)
            message(not_found[0] + "...", type="WARNING")

            if force_tx_class:
                msg_type = "WARNING"
            else:
                msg_type = "ERROR"

            message("Some transcripts (n={n}) provided with --transcript-class where not"
                    " found in input file (overriden with --force-tx-class).".format(
                n=len(not_found)),
                type=msg_type)

            if force_tx_class:
                tx_to_class = {
                    k: tx_to_class[k] for k in tx_to_class if k in intersect}
                nb_tx_classes = len(list(set(tx_to_class.values())))

                if nb_tx_classes == 0:
                    message("No transcript found in file provided through "
                            "--target-tx-file.", type="ERROR")

        else:
            message("Found %d transcripts of interest in "
                    "input file." % nb_intersect)

    else:

        transcript_file_tmp = make_tmp_file(prefix="transcript_classes")
        class_list.add("all_regions")
        nb_tx_classes = 1

        for tx_id in input_file_tx:
            transcript_file_tmp.write(tx_id + "\tall_regions\n")
            tx_to_class[tx_id] = "all_regions"

        transcript_file_tmp.close()
        transcript_file = open(transcript_file_tmp.name, "r")

    # -------------------------------------------------------------------------
    #
    # Colors for profiles
    #
    # -------------------------------------------------------------------------

    def get_list_of_colors_mpl(number, pal='nipy_spectral'):

        colormap = cm.get_cmap(pal, lut=number)
        colors = [mpl.colors.rgb2hex(colormap(i)) for i in np.linspace(0., 1., number)]
        return colors

    if profile_colors is None:

        if group_by == 'bwig':
            profile_colors = get_list_of_colors_mpl(len(input_file_bwig))
        elif group_by == 'tx_classes':
            profile_colors = get_list_of_colors_mpl(len(class_list))
        elif group_by == 'chrom':
            profile_colors = get_list_of_colors_mpl(len(input_file_chrom))


    else:
        profile_colors = profile_colors.split(",")

    # -------------------------------------------------------------------------
    #
    # Colors orders
    #
    # -------------------------------------------------------------------------

    if color_order is None:
        if group_by == 'bwig':
            color_order = ",".join(input_file_bwig)
        elif group_by == 'tx_classes':
            color_order = ",".join(class_list)
        elif group_by == 'chrom':
            color_order = ",".join(list(input_file_chrom))
        else:
            message("color_order is undefined.", type="ERROR")
        color_order_list = color_order.split(",")

    else:
        color_order_list = color_order.split(",")
        color_order_pb = False

        if group_by == 'bwig':
            if len(color_order_list) != len(input_file_bwig):
                color_order_pb = True
            if len(set(color_order_list)) != len(set(input_file_bwig)):
                color_order_pb = True
            for co in color_order_list:
                if co not in input_file_bwig:
                    color_order_pb = True

        elif group_by == 'tx_classes':
            if len(color_order_list) != len(class_list):
                color_order_pb = True
            if len(set(color_order_list)) != len(set(class_list)):
                color_order_pb = True
            for co in color_order_list:
                if co not in class_list:
                    color_order_pb = True

        elif group_by == 'chrom':
            if len(color_order_list) != len(list(input_file_chrom)):
                color_order_pb = True
            if len(set(color_order_list)) != len(set(list(input_file_chrom))):
                color_order_pb = True
            for co in color_order_list:
                if co not in list(input_file_chrom):
                    color_order_pb = True

        else:
            color_order_pb = True
        if color_order_pb:
            message("Please, check --color-order.", type="ERROR")

    if group_by == 'bwig':
        if len(input_file_bwig) > len(profile_colors):
            msg = "Need more colors for displaying bigwig files (n=%d)"
            message(msg % len(input_file_bwig),
                    type="ERROR")
        for curr_item in color_order_list:
            if curr_item not in input_file_bwig:
                message("Color order: Found undefined bigwig labels (" + curr_item + ")... Please Check.",
                        type="WARNING")
                message("Use one of :" + ",".join(input_file_bwig) + ".",
                        type="ERROR")

    elif group_by == 'tx_classes':
        if nb_tx_classes > len(profile_colors):
            msg = "Need more colors for displaying transcript classes (n=%d)"
            message(msg % nb_tx_classes,
                    type="ERROR")
        for curr_item in color_order_list:
            if curr_item not in class_list:
                message("{a} not found in {b}... Please check.".format(a=str(curr_item),
                                                                       b=str(class_list)),
                        type="WARNING")
                message("The class name may contain special characters...",
                        type="WARNING")

                message("Found undefined class labels... Please check.",
                        type="ERROR")
    elif group_by == 'chrom':
        if len(input_file_chrom) > len(profile_colors):
            msg = "Need more colors for displaying chromosome classes (n=%d)"
            message(msg % len(input_file_chrom),
                    type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Prepare output files
    #
    # -------------------------------------------------------------------------

    if config['ft_type'] == 'promoter':

        img_file = "promoter_u%s_d%s." + page_format
        img_file = img_file % (config['from'], config['to'])

    elif config['ft_type'] == 'tts':

        img_file = "tts_u%s_d%s." + page_format
        img_file = img_file % (config['from'], config['to'])

    elif config['ft_type'] == 'transcript':

        img_file = "transcript_u%s_d%s." + page_format
        img_file = img_file % (config['from'], config['to'])

    elif config['ft_type'] == 'user_regions':
        img_file = "user_regions_u%s_d%s." + page_format
        img_file = img_file % (config['from'], config['to'])

    elif config['ft_type'] == 'single_nuc':
        img_file = "user_positions_u%s_d%s." + page_format
        img_file = img_file % (config['from'], config['to'])

    file_out_list = make_outdir_and_file(out_dir,
                                         ["profile_stats.txt",
                                          img_file],
                                         force=True)

    data_file, img_file = file_out_list

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
    # Read tab-delimited file
    #
    # -------------------------------------------------------------------------

    data = pd.read_csv(inputfile_main.name, sep="\t", header=1)

    # -------------------------------------------------------------------------
    #
    # Find coverage columns
    #
    # -------------------------------------------------------------------------

    pos_order = []

    for i in data.columns:
        if re.search('(main_\d+)|(upstream_\d+)|(downstream_\d+)', i):
            pos_order += [i]

    bin_nb_main = len([x for x in data.columns if "main" in x])
    bin_nb_ups = len([x for x in data.columns if "upstream" in x])
    bin_nb_dws = len([x for x in data.columns if "downstream" in x])
    bin_nb_total = bin_nb_ups + bin_nb_main + bin_nb_dws

    # -------------------------------------------------------------------------
    #
    # Read transcript file
    #
    # -------------------------------------------------------------------------

    # all tx/gene names of the dataframe
    all_tx = list(OrderedDict.fromkeys(data['gene'].tolist()))

    if group_by == 'tx_classes' or facet_var == 'tx_classes':

        # -------------------------------------------------------------------------
        # Get the transcript classes
        # -------------------------------------------------------------------------

        message("Reading transcript file.")
        df_classes = pd.read_csv(transcript_file.name, sep='\t', header=None)
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
    # Melting (data -> dm)
    #
    # -------------------------------------------------------------------------

    dm = data.melt(id_vars=['tx_classes', 'bwig', 'chrom', 'start', 'end', 'gene', 'strand'], value_vars=pos_order)
    dm = dm.rename(columns={'variable': 'pos', 'value': 'exprs'})
    dm['bwig'] = Categorical(dm['bwig'])

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
        y_lab = "log2(Signal)"

    else:
        y_lab = "Signal"

    if normalization_method == 'ranging':
        message('Normalizing (ranging)')
        for k in dm['bwig'].unique():
            tmp = np.array(dm[dm.bwig == k]['exprs'].tolist())
            tmp_norm = (tmp - min(tmp)) / (max(tmp) - min(tmp)) * 100
            dm.loc[dm.bwig == k, 'exprs'] = tmp_norm

        y_lab = "scaled(" + y_lab + ", %)"

    # -------------------------------------------------------------------------
    #
    # Compute the statistics
    #
    # -------------------------------------------------------------------------

    if facet_var is None:
        df_groups = [group_by] + ['pos']
    else:
        df_groups = [group_by] + [facet_var] + ['pos']

    dm = dm.groupby(df_groups, as_index=False).agg({'exprs': [stat_fun, np.std]})

    # -------------------------------------------------------------------------
    #
    # Multi-indexed dataframes (various column name levels) are not easy to
    # handle
    # -------------------------------------------------------------------------

    dm.columns = ["_".join(x) if len(x) > 1 and x[1] != '' else x[0] for x in dm.columns.ravel()]

    # fixed a weird behavior of pandas
    dm.columns = [x if x != "exprs_amin" else 'exprs_min' for x in dm.columns]
    dm.columns = [x if x != "exprs_amax" else 'exprs_max' for x in dm.columns]

    # -------------------------------------------------------------------------
    #
    # Column ordering
    #
    # -------------------------------------------------------------------------

    message("Computing column ordering.")

    pos_order = []

    tmp = [x for x in dm.pos.unique() if 'upstream' in x]
    pos_order += sorted(sorted(tmp, reverse=True), key=lambda x: int(x.split("_")[1]))

    tmp = [x for x in dm.pos.unique() if 'main' in x]
    pos_order += sorted(sorted(tmp, reverse=True), key=lambda x: int(x.split("_")[1]))

    tmp = [x for x in dm.pos.unique() if 'downstream' in x]
    pos_order += sorted(sorted(tmp, reverse=True), key=lambda x: int(x.split("_")[1]))

    dm.pos = Categorical(dm.pos.tolist(), categories=pos_order, ordered=True)

    # -------------------------------------------------------------------------
    #
    # Turning x axis into continuous scale if needed
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

    p = ggplot(data=dm,
               mapping=aes(x='pos',
                           y="exprs_" + stat,
                           color=group_by))

    p += geom_line(size=line_width)

    p += ylab(y_lab)
    # -------------------------------------------------------------------------
    #
    # Theming
    #
    # -------------------------------------------------------------------------

    message("Theming and ordering. Please be patient...")
    p += theme(legend_title=element_blank())
    theme_plotnine = getattr(plotnine, theme_plotnine)
    p += theme_plotnine()

    p += theme(legend_position="top",
               legend_title=element_blank(),
               legend_key=element_rect(colour="white", fill="white"),
               legend_text=element_text(size=8),
               axis_text=element_text(size=axis_text, angle=40, hjust=1),
               strip_text_x=element_text(size=strip_text, colour='white'),
               strip_background=element_rect(fill=strip_color)
               )

    p += guides(col=guide_legend(ncol=5))

    p += ggtitle(title)

    # -------------------------------------------------------------------------
    #
    # Preparing x axis
    #
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    #
    # Adding rectangle overlays to highlight 5' and 3' regions
    #
    # -------------------------------------------------------------------------

    if config['ft_type'] in ["transcript", "user_regions"]:
        if config['from']:
            message("Highlighting upstream regions")

            rectangles = {'xmin': [0],
                          'xmax': [bin_nb_ups / bin_nb_total * 100],
                          'ymin': dm["exprs_" + stat].min(),  # np.Inf is not accepted.
                          'ymax': dm["exprs_" + stat].max()  # np.Inf is not accepted.
                          }

            rectangles = pd.DataFrame(data=rectangles)

            p += geom_rect(data=rectangles,
                           mapping=aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax'),
                           fill='lightslategray',
                           alpha=0.3,
                           inherit_aes=False,
                           show_legend=False)

        if config['to']:
            message("Highlighting downstream regions")
            rectangles = {'xmin': [(bin_nb_total - bin_nb_dws) / bin_nb_total * 100],
                          'xmax': [100],
                          'ymin': dm["exprs_" + stat].min(),  # np.Inf is not accepted.
                          'ymax': dm["exprs_" + stat].max()  # np.Inf is not accepted.
                          }

            rectangles = pd.DataFrame(data=rectangles)

            p += geom_rect(data=rectangles,
                           mapping=aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax'),
                           fill='lightslategray',
                           alpha=0.3,
                           inherit_aes=False,
                           show_legend=False)

    # --------------------------------------------------------------------------
    #
    #
    #
    # --------------------------------------------------------------------------

    if facet_var is not None:
        # , scales="free_y", space="free_y"
        p += facet_wrap("~ " + facet_var)

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

    # Delete temporary dir
    shutil.rmtree(dir_name)


if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    draw_profile(**args)

else:

    test = '''
            
            #profile: prepare dataset
            @test "profile_1" {
             result=`gtftk get_example -d mini_real -f '*'; gtftk overlapping -i mini_real.gtf.gz -c hg38.genome  -n > mini_real_noov.gtf; gtftk random_tx -i mini_real_noov.gtf  -m 1 -s 123 > mini_real_noov_rnd_tx.gtf`
              [ -s "hg38.genome" ]
            }
        
            #profile: prepare dataset
            @test "profile_2" {
             result=`gtftk mk_matrix -i mini_real_noov_rnd_tx.gtf -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_promoter_pr`
              [ -s "mini_real_promoter_pr.zip" ]
            }
        
            #profile: test
            @test "profile_3" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -o profile_prom_1 -pf png -if example_01.png`
              [ -s "example_01.png" ]
            }
        
            #profile: make mini_real_promoter
            @test "profile_4" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -o profile_prom_1 -pf png -if example_01.png`
              [ -s "example_01.png" ]
            }
        
            #profile: make tss.bed
            @test "profile_5" {
             result=`gtftk select_by_key -i mini_real_noov_rnd_tx.gtf -k feature -v transcript |  gtftk 5p_3p_coord > tss.bed`
              [ -s "tss.bed" ]
            }
        
            #profile: make single_nuc
            @test "profile_6" {
             result=`gtftk mk_matrix -u 5000 -d 5000 -i tss.bed -w 200 -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_single_nuc_pr -c hg38.genome -t single_nuc`
              [ -s "mini_real_single_nuc_pr.zip" ]
            }
        
            #profile: test single_nuc
            @test "profile_7" {
             result=`gtftk profile -i mini_real_single_nuc_pr.zip -o profile_prom_1a -pf png -if example_01a.png`
              [ -s "example_01a.png" ]
            }
        
            #profile: make mini_real_tx
            @test "profile_8" {
             result=`gtftk mk_matrix -i mini_real_noov_rnd_tx.gtf -t transcript  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tx_pr`
              [ -s "mini_real_tx_pr.zip" ]
            }
        
            #profile: make mini_real_tx
            @test "profile_9" {
             result=`gtftk select_by_key -i mini_real_noov_rnd_tx.gtf -k feature -v transcript | gtftk convert -f bed6 > mini_real_rnd_tx.bed`
              [ -s "mini_real_rnd_tx.bed" ]
            }
        
            #profile: test mini_real_tx
            @test "profile_10" {
             result=`gtftk profile -D -i mini_real_tx_pr.zip -o profile_tx_1 -pf png -if example_02.png`
              [ -s "example_02.png" ]
            }
        
            #profile: make mini_real_user_def
            @test "profile_11" {
             result=`gtftk mk_matrix --bin-around-frac 0.5 -i mini_real_rnd_tx.bed -t user_regions  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_user_def`
              [ -s "mini_real_user_def.zip" ]
            }
        
            #profile: test mini_real_user_def
            @test "profile_12" {
             result=`gtftk profile -D -i mini_real_user_def.zip -o profile_udef_4  -pf png -if example_04.png`
              [ -s "example_04.png" ]
            }
        
            #profile: test mini_real_user_def
            @test "profile_13" {
             result=`gtftk profile -D -nm ranging -i mini_real_user_def.zip -o profile_udef_5  -pf png -if example_04b.png`
              [ -s "example_04b.png" ]
            }
                
            #profile: test mini_real_user_def
            @test "profile_14" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -g tx_classes -f bwig -o profile_prom_2  -ph 5 -c "#23AF36" -pf png -if example_05.png`
              [ -s "example_05.png" ]
            }
        
            #profile: create dataset
            @test "profile_15" {
             result=`gtftk tabulate -k transcript_id,gene_biotype -i mini_real_noov_rnd_tx.gtf -H | sort | uniq | perl -ne 'print if (/(protein_coding)|(lincRNA)|(antisense)|(processed_transcript)/)'> tx_classes.txt`
              [ -s "tx_classes.txt" ]
            }
        
            #profile: create dataset
            @test "profile_16" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -g tx_classes -f bwig -o profile_prom_2  -ph 5 -c "#23AF36" -pf png -if example_05.png`
              [ -s "example_05.png" ]
            }
        
            #profile: create dataset
            @test "profile_17" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -g bwig -f tx_classes  -o profile_prom_3  -ph 4 -c "#66C2A5,#FC8D62,#8DA0CB" -t tx_classes.txt  -pf png -if example_06.png`
              [ -s "example_06.png" ]
            }
        
            #profile: create dataset
            @test "profile_18" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -g tx_classes -f bwig  -o profile_prom_4  -ph 4 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -t tx_classes.txt  -pf png -if example_07.png`
              [ -s "example_07.png" ]
            }
        
            #profile: create dataset
            @test "profile_19" {
             result=`gtftk mk_matrix --bin-around-frac 0.5 -i mini_real_noov_rnd_tx.gtf -t transcript  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tx_pr_2`
              [ -s "mini_real_tx_pr_2.zip" ]
            }
            
            #profile: create dataset
            @test "profile_20" {
             result=`gtftk profile -D -i mini_real_tx_pr_2.zip -g tx_classes -f bwig  -o profile_tx_3 -pw 12  -ph 7 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -t tx_classes.txt  -pf png -if example_08.png`
              [ -s "example_08.png" ]
            }
        
            #profile: create dataset
            @test "profile_21" {
             result=`gtftk profile -D -i mini_real_promoter_pr.zip -g bwig -f chrom  -o profile_prom_5  -ph 15 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF"   -pf png -if example_09.png`
              [ -s "example_09.png" ]
            }
             
            #profile: create dataset
            @test "profile_22" {
             result=`gtftk profile -th classic -D -i mini_real_promoter_pr.zip -g bwig -f chrom  -o profile_prom_5  -ph 15 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF"   -pf png -if example_09b.png`
              [ -s "example_09b.png" ]
            }
            
        
            '''

    cmd = CmdObject(name="profile",
                    message="Create coverage profile using a bigWig as input.",
                    parser=make_parser(),
                    fun=draw_profile,
                    desc=__doc__,
                    updated=__updated__,
                    notes=__notes__,
                    references=__references__,
                    group="coverage",
                    test=test,
                    rlib=R_LIB)
