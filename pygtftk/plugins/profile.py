#!/usr/bin/env python
"""
 Produces bigWig coverage profiles using calls to plotnine graphic package.
"""

import argparse
import os
import re
import shutil
import sys
import warnings
import zipfile
from collections import OrderedDict
from zipfile import BadZipFile

import matplotlib as mpl
import numpy as np
import pandas as pd
import plotnine
from matplotlib import cm
from matplotlib import colors as mcolors
from pandas import Categorical
from plotnine import aes, geom_text, scale_x_continuous, scale_color_manual, guide_legend, guides, ggtitle, \
    element_rect, element_blank, element_text, element_line, theme, facet_wrap, geom_rect
from plotnine import geom_line
from plotnine import ggplot
from plotnine import xlab
from plotnine import ylab
from plotnine.exceptions import PlotnineError

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import ALL_MPL_PALETTES
from pygtftk.utils import GTFtkError
from pygtftk.utils import chomp
from pygtftk.utils import is_hex_color
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import make_tmp_dir
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
 -- Think about using normalized bigWig files as input to mk_matrix. This
 will limit the requirement for an additional normalization step (see
 Deeptools for a set of useful methods implemented in bamCoverage/bamCompare).
"""

'''
 -- The ranging normalization method [1] implies the following transformation:
 -- -  (x_i - min(x))/(max(x) - min(x)).
'''

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
                            type=arg_formatter.FormattedFile(mode='r', file_ext='zip'),
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

    parser_grp.add_argument('-e', '--confidence-interval',
                            help="Add a confidence interval to estimate standard error of the mean.",
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-c', '--profile-colors',
                            help='Colors.',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-d', '--color-order',
                            help='Factor ordering. comma-separated bwig labels or tx classes.',
                            default=None,
                            type=str)

    parser_grp.add_argument('-g', '--group-by',
                            help='The variable used for grouping.',
                            default=None,
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
                            help='Output pdf file width (e.g. 7 inches).',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=None,
                                                          val_type="float",
                                                          linc=False),
                            default=None,
                            required=False)

    parser_grp.add_argument('-ph', '--page-height',
                            help='Output  file height (e.g. 5 inches).',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=None,
                                                          val_type="float",
                                                          linc=False),
                            default=None,
                            required=False)

    parser_grp.add_argument('-pf', '--page-format',
                            help='Output file format.',
                            choices=['pdf', 'png'],
                            default='pdf',
                            required=False)

    parser_grp.add_argument('-lw', '--line-width',
                            help='Line width.',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=None,
                                                          val_type="float",
                                                          linc=False),
                            default=1.25,
                            required=False)

    parser_grp.add_argument('-bc', '--border-color',
                            help='Border color for the plot.',
                            default="#777777",
                            type=str,
                            required=False)

    parser_grp.add_argument('-x', '--x-lab',
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

    parser_grp.add_argument('-u', '--subset-bwig',
                            help='Use only a subset of the bigwigs for plotting',
                            default=None,
                            type=str,
                            required=False)

    parser_grp.add_argument('-fc', '--facet-col',
                            help='Number of facet columns.',
                            default=4,
                            type=arg_formatter.ranged_num(lowest=1,
                                                          highest=None,
                                                          val_type="int",
                                                          linc=True),
                            required=False)

    parser_grp.add_argument('-w', '--show-group-number',
                            help='Show the number of element per group.',
                            action="store_true",
                            required=False)

    parser_grp.add_argument('-if', '--user-img-file',
                            help="Provide an alternative path for the image.",
                            default=None,
                            type=argparse.FileType("w"),
                            required=False)

    parser_grp.add_argument('-ul',
                            '--upper-limit',
                            type=arg_formatter.ranged_num(lowest=0,
                                                          highest=1,
                                                          val_type="float",
                                                          linc=False),
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
                            type=arg_formatter.ranged_num(lowest=50,
                                                          highest=None,
                                                          val_type="int",
                                                          linc=True),
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
                                'seaborn'],
                            default='bw',
                            help="The theme for plotnine diagram.",
                            required=False)

    parser_grp.add_argument('-m', '--palette',
                            help='A color palette (see: https://tinyurl.com/ydacyfxx).',
                            default="nipy_spectral",
                            required=False)

    parser_grp.add_argument('-l', '--list-bwig',
                            help='List the bigwig files in the matrix file..',
                            action="store_true")

    return parser


def profile(inputfile=None,
            out_dir=None,
            group_by='bwig',
            color_order=None,
            transcript_file=None,
            normalization_method=None,
            to_log=False,
            upper_limit=0.95,
            profile_colors=None,
            palette='nipy_spectral',
            page_width=None,
            title=None,
            page_height=None,
            page_format='pdf',
            user_img_file=None,
            facet_col=None,
            border_color="#BBBBBB",
            stat="mean",
            facet_var=None,
            x_lab="Selected genomic regions",
            axis_text=8,
            strip_text=8,
            subset_bwig=None,
            show_group_number=False,
            line_width=1,
            theme_plotnine='bw',
            list_bwig=False,
            confidence_interval=False,
            dpi=300):
    # -------------------------------------------------------------------------
    #
    # Pandas version is sometimes problematic
    #
    # -------------------------------------------------------------------------

    message("Using pandas version " + pd.__version__, type="DEBUG")
    message("Pandas location " + pd.__file__, type="DEBUG")

    message("Using numpy version " + np.__version__, type="DEBUG")
    message("Pandas numpy " + np.__file__, type="DEBUG")

    message("Using plotnine version " + plotnine.__version__, type="DEBUG")
    message("Pandas plotnine " + plotnine.__file__, type="DEBUG")

    # -------------------------------------------------------------------------
    #
    # The selected theme
    #
    # -------------------------------------------------------------------------

    theme_plotnine = 'theme_' + theme_plotnine

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
    dir_name = make_tmp_dir(prefix='profile_matrix_')

    # Uncompressing the file
    message("Uncompressing : " + dir_name,
            type="DEBUG")

    try:
        with zipfile.ZipFile(inputfile.name) as zf:
            zf.extractall(dir_name)

    except BadZipFile:
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

    if list_bwig:
        message("Bigwig list: " + ",".join(input_file_bwig), force=True)
        sys.exit(0)

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
    # Check whether we are using  a subset
    #
    # -------------------------------------------------------------------------

    if subset_bwig is not None:
        subset_bwig_list = subset_bwig.split(",")
        if len(subset_bwig_list) == 0:
            message("--subset-bwig is empty", type="ERROR")
        for i in subset_bwig_list:
            if i not in input_file_bwig:
                message("Bwig " + i + " not found. Check --subset-bwig", type="ERROR")
        input_file_bwig = [x for x in input_file_bwig if x in subset_bwig_list]

    # -------------------------------------------------------------------------
    #
    # Check arguments: --facet, --group-by according to the number of bigwig
    #
    # -------------------------------------------------------------------------
    # If one is analyzing more than one bigwig
    # the "bwig" factor should appear in  --facet or --group-by

    if group_by is None:
        group_by = "bwig"
        message("--group-by not set. Choosing 'bwig'.", type="WARNING")

    if group_by == facet_var:
        message("--facet-var and --group-by should be different.",
                type="ERROR")

    if len(input_file_bwig) > 1:
        if facet_var != "bwig":
            if group_by != "bwig":
                message("If more than one bigWig is analyzed, --facet or --group-by should be set to 'bwig'.",
                        type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Read coverage file
    #
    # -------------------------------------------------------------------------

    data = pd.read_csv(inputfile_main.name, sep="\t", header=1)
    # import pickle
    # data = pickle.load(file=open('dm_pickle.pic'))

    if data.shape[0] == 0:
        message("No lines found in input file.",
                type="ERROR")

    # -------------------------------------------------------------------------
    #
    # Read transcript file
    #
    # -------------------------------------------------------------------------

    # all tx/gene names of the dataframe
    all_tx = list(OrderedDict.fromkeys(data['gene'].tolist()))

    if group_by == 'tx_classes' or facet_var == 'tx_classes':
        if transcript_file is not None:
            # -------------------------------------------------------------------------
            # Get the transcript classes
            # -------------------------------------------------------------------------

            message("Reading transcript file.")
            df_classes = pd.read_csv(transcript_file.name, sep='\t', header=None)

            if df_classes.shape[0] == 0:
                message("No lines found in transcript file.",
                        type="ERROR")

            if len(df_classes.columns) < 2:
                message("The transcript file should contain at least two columns.",
                        type="ERROR")

            message("Deleting duplicates in transcript-file.")
            df_classes = df_classes.drop_duplicates(subset=[0])
            tx_ordering = df_classes[0].tolist()
            tx_classes = OrderedDict(list(zip(df_classes[0], df_classes[1])))
            class_list = set(df_classes[1])

            # -------------------------------------------------------------------------
            # Select the transcript of interest and add classes info to the data.frame
            # -------------------------------------------------------------------------

            message("Checking how many transcript where found in the transcript list.")

            nb_retained = len([x for x in all_tx if x in tx_ordering])

            msg = "Keeping {a} transcript out of {b} in input transcript list.".format(a=nb_retained, b=len(all_tx))
            message(msg)

            # subsetting
            data = data[[True if x in tx_ordering else False for x in data['gene'].tolist()]]
            data = data.assign(tx_classes=[tx_classes[x] for x in data['gene'].tolist()])
        else:
            data = data.assign(tx_classes=["All transcripts" for x in data['gene'].tolist()])
            class_list = ["All transcripts"]
    else:
        data = data.assign(tx_classes=["All transcripts" for x in data['gene'].tolist()])
        class_list = ["All transcripts"]

    # -------------------------------------------------------------------------
    #
    # Colors for profiles
    #
    # -------------------------------------------------------------------------

    if palette not in ALL_MPL_PALETTES:
        message("Sorry but the palette is unknown.")

    def get_list_of_colors_mpl(number, pal=palette):

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
            raise GTFtkError('--group-by is unknown')


    else:
        profile_colors = profile_colors.split(",")
        profile_colors = [x.strip() for x in profile_colors]

        if subset_bwig is not None:
            profile_colors = profile_colors[:len(subset_bwig_list)]

        mcolors_name = mcolors.cnames

        for i in profile_colors:
            if i not in mcolors_name:
                if not is_hex_color(i):
                    message(i + " is not a valid color. Please fix.", type="ERROR")

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
        color_order = color_order.split(",")

    else:
        color_order = color_order.split(",")
        color_order_pb = False

        if group_by == 'bwig':
            if len(color_order) != len(input_file_bwig):
                color_order_pb = True
            if len(set(color_order)) != len(set(input_file_bwig)):
                color_order_pb = True
            for co in color_order:
                if co not in input_file_bwig:
                    color_order_pb = True

        elif group_by == 'tx_classes':
            if len(color_order) != len(class_list):
                color_order_pb = True
            if len(set(color_order)) != len(set(class_list)):
                color_order_pb = True
            for co in color_order:
                if co not in class_list:
                    color_order_pb = True

        elif group_by == 'chrom':
            if len(color_order) != len(list(input_file_chrom)):
                color_order_pb = True
            if len(set(color_order)) != len(set(list(input_file_chrom))):
                color_order_pb = True
            for co in color_order:
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
        for curr_item in color_order:
            if curr_item not in input_file_bwig:
                message("Color order: Found undefined bigwig labels (" + curr_item + ")... Please Check.",
                        type="WARNING")
                message("Use one of : " + ",".join(input_file_bwig) + ".",
                        type="ERROR")

    elif group_by == 'tx_classes':
        if len(class_list) > len(profile_colors):
            msg = "Need more colors for displaying transcript classes (n=%d)"
            message(msg % len(class_list),
                    type="ERROR")
        for curr_item in color_order:
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
    else:
        raise GTFtkError("--group-by is unknown.")

    if len(color_order) < len(profile_colors):
        profile_colors = profile_colors[:len(color_order)]

    message("Color order : " + str(color_order), type="DEBUG")
    message("Profile color : " + str(profile_colors), type="DEBUG")

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
    else:
        raise GTFtkError("Unknown feature type.")

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
    # Find coverage columns
    #
    # -------------------------------------------------------------------------

    message("Searching coverage columns.")

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
    # Melting (data -> dm)
    #
    # -------------------------------------------------------------------------

    message("Melting.")
    dm = data.melt(id_vars=['tx_classes', 'bwig', 'chrom', 'start', 'end', 'gene', 'strand'], value_vars=pos_order)
    dm = dm.rename(columns={'variable': 'pos', 'value': 'exprs'})
    dm['bwig'] = Categorical(dm['bwig'])

    # -------------------------------------------------------------------------
    #
    # Subsetting dm if required
    #
    # -------------------------------------------------------------------------

    dm = dm[dm['bwig'].isin(input_file_bwig)]

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
    # Compute the statistics to be plotted
    #
    # -------------------------------------------------------------------------

    if facet_var is None:
        df_aggr_var = [group_by] + ['pos']
    else:
        df_aggr_var = [group_by] + [facet_var] + ['pos']

    dm = dm.groupby(df_aggr_var, as_index=False).agg({'exprs': ['mean', 'median', 'mad',
                                                                'std', 'min', 'max',
                                                                'sum', len]})

    # -------------------------------------------------------------------------
    #
    # Multi-indexed dataframes (various column name levels) are not easy to
    # handle
    # -------------------------------------------------------------------------

    dm.columns = [x[1] if len(x) > 1 and x[1] != '' else x[0] for x in dm.columns.ravel()]

    # -------------------------------------------------------------------------
    #
    # Compute confidence interval
    #
    # -------------------------------------------------------------------------

    if stat == "mean":
        std_err = dm['std'] / np.sqrt(dm['len'])
        dm['ci_low'] = dm['mean'] - 1.96 * std_err
        dm['ci_high'] = dm['mean'] + 1.96 * std_err
    elif stat == "median":
        dm['ci_low_robust'] = dm['median'] - 1.96 * dm['mad'] / np.sqrt(dm['len'])
        dm['ci_high_robust'] = dm['median'] + 1.96 * dm['mad'] / np.sqrt(dm['len'])

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

    fr = round(int(config['from']), 0)
    to = round(int(config['to']), 0)

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
    # Compute the number of observations per panel/facet
    #
    # -------------------------------------------------------------------------

    if show_group_number:
        df_aggr_var.remove('pos')
        df_aggr_var += ['len']
        dm_nb = dm
        dm_nb = dm_nb.groupby(df_aggr_var, as_index=False).agg({stat: ['max', 'min']})
        dm_nb.columns = [x[1] if len(x) > 1 and x[1] != '' else x[0] for x in dm_nb.columns.ravel()]
        max_val = dm_nb['max'].iloc[dm_nb['max'].idxmax(),]
        min_val = dm_nb['min'].iloc[dm_nb['min'].idxmax(),]
        y_pos = np.linspace(start=(max_val - min_val) / 1.5, stop=max_val, num=len(dm[group_by].unique()))
        for i, j in zip(y_pos, dm[group_by].unique()):
            tmp = [i] * dm_nb['max'][dm_nb[group_by] == j].shape[0]
            dm_nb.loc[dm_nb[group_by] == j, 'max'] = tmp
        if config['ft_type'] in ["transcript", "user_regions"]:
            dm_nb = dm_nb.assign(x=pd.Series(np.repeat(85, dm_nb.shape[0])).values)
        else:
            dm_nb = dm_nb.assign(x=pd.Series(np.repeat(int(config["to"]) - int(config["to"]) / 3,
                                                       dm_nb.shape[0])).values)
        dm_nb = dm_nb.assign(hjust=pd.Series(np.repeat(0, dm_nb.shape[0])).values)
        dm_nb.len = dm_nb.len.astype(int)
        dm_nb['nb_obs'] = 'n=' + dm_nb['len'].astype(str)

        # transfert this to the main dataframe to avoid problems...

        if facet_var is None:
            dm = dm.join(dm_nb.set_index([group_by]), on=[group_by], rsuffix="_grp")
        else:
            dm = dm.join(dm_nb.set_index([group_by, facet_var]), on=[group_by, facet_var], rsuffix="_grp")

    # -------------------------------------------------------------------------
    #
    # add colors to matrix
    #
    # -------------------------------------------------------------------------
    group2cols = dict(list(zip(color_order, profile_colors)))
    dm['color_palette'] = [group2cols[x] for x in dm[group_by]]

    # -------------------------------------------------------------------------
    #
    # Preparing diagram
    #
    # -------------------------------------------------------------------------
    dm[group_by] = Categorical(dm[group_by], ordered=True, categories=color_order)
    dm['color_palette'] = Categorical(dm['color_palette'], ordered=True, categories=profile_colors)

    message("Preparing diagram")

    p = ggplot(data=dm,
               mapping=aes(x='pos',
                           y=stat,
                           color=group_by))

    p += ylab(y_lab)
    p += xlab(x_lab)

    # -------------------------------------------------------------------------
    #
    # Adding ribbon to delineate sd
    #
    # -------------------------------------------------------------------------

    # TODO: improve this. Particularly buggy (I guess) with the current version of
    # plotnine

    if confidence_interval:
        if facet_var is not None:
            if stat == "mean":
                for i in range(len(color_order)):
                    for j in dm[facet_var].unique():
                        dm_sub = dm[(dm[group_by] == color_order[i]) & (dm[facet_var] == j)].sort_values('pos')
                        p += plotnine.geom_ribbon(data=dm_sub,
                                                  mapping=aes(ymin='ci_low',
                                                              ymax='ci_high'),
                                                  show_legend=False,
                                                  fill=list(dm_sub.color_palette.unique())[0],
                                                  color=None,
                                                  alpha=0.3)
            elif stat == "median":
                for i in range(len(color_order)):
                    for j in dm[facet_var].unique():
                        dm_sub = dm[(dm[group_by] == color_order[i]) & (dm[facet_var] == j)].sort_values('pos')
                        p += plotnine.geom_ribbon(data=dm_sub,
                                                  mapping=aes(ymin='ci_low_robust',
                                                              ymax='ci_high_robust'),
                                                  show_legend=False,
                                                  fill=list(dm_sub.color_palette.unique())[0],
                                                  color=None,
                                                  alpha=0.3)
        else:
            if stat == "mean":
                for i in range(len(color_order)):
                    dm_sub = dm[dm[group_by] == color_order[i]].sort_values('pos')
                    p += plotnine.geom_ribbon(data=dm_sub,
                                              mapping=aes(ymin='ci_low',
                                                          ymax='ci_high'),
                                              show_legend=False,
                                              fill=list(dm_sub.color_palette.unique())[0],
                                              color=None,
                                              alpha=0.3)
            elif stat == "median":
                for i in range(len(color_order)):
                    dm_sub = dm[dm[group_by] == color_order[i]].sort_values('pos')
                    p += plotnine.geom_ribbon(data=dm_sub,
                                              mapping=aes(ymin='ci_low_robust',
                                                          ymax='ci_high_robust'),
                                              show_legend=False,
                                              fill=list(dm_sub.color_palette.unique())[0],
                                              color=None,
                                              alpha=0.3)

    # -------------------------------------------------------------------------
    #
    # Theming
    #
    # -------------------------------------------------------------------------

    message("Theming and ordering. Please be patient...")
    p += theme(legend_title=element_blank())
    theme_plotnine_fun = getattr(plotnine, theme_plotnine)
    p += theme_plotnine_fun()

    # remove major/minor grid due to
    # weird placements by default
    # in this plotnine version

    p += theme(legend_position="top",
               legend_title=element_blank(),
               legend_key=element_rect(colour="white", fill="white"),
               legend_text=element_text(size=8),
               axis_text_x=element_text(size=axis_text, angle=40, hjust=1.5),
               axis_text_y=element_text(size=axis_text),
               axis_ticks=element_line(colour=border_color),
               axis_line=element_line(colour=border_color, size=1),
               axis_line_y=element_line(colour=border_color, size=1),
               panel_spacing_x=0.3,
               panel_spacing_y=0.3,
               strip_text_x=element_text(size=strip_text, colour='white'),
               strip_background=element_rect(fill=border_color, colour=border_color),
               panel_border=element_rect(colour=border_color, size=1),
               panel_grid_minor=element_blank()
               )

    p += ggtitle(title)
    p += guides(col=guide_legend(ncol=5))

    # -------------------------------------------------------------------------
    #
    # Preparing x axis
    #
    # -------------------------------------------------------------------------

    message("Preparing x axis")

    if config['ft_type'] in ["transcript", "user_regions"]:

        if config['from']:

            if config['to']:
                ticks = [0] + list(np.linspace(bin_nb_ups, bin_nb_main + bin_nb_ups, 11)) + [bin_nb_total]
                ticks = [x / bin_nb_total * 100 for x in ticks]
                labels = [int(-fr)] + [str(int(round(x, 0))) + "%" for x in np.linspace(0, 100, 11)] + [int(to)]


            else:
                ticks = [0] + list(np.linspace(bin_nb_ups, bin_nb_total, 11))
                ticks = [x / bin_nb_total * 100 for x in ticks]

                labels = [int(- fr)
                          ] + [str(int(round(x, 0))) + "%" for x in np.linspace(0, 100, 6)]
        else:
            if config['to']:

                ticks = list(np.linspace(0, bin_nb_main, 11)) + [bin_nb_total]
                ticks = [x / bin_nb_total * 100 for x in ticks]

                labels = [str(int(round(x, 0))) + "%" for x in np.linspace(0, 100, 6)] + [int(to)]


            else:
                ticks = list(np.linspace(0, bin_nb_total, 6))
                ticks = [x / bin_nb_total * 100 for x in ticks]
                labels = [str(int(round(x, 0))) + "%" for x in np.linspace(0, 100, 6)]

        p += scale_x_continuous(expand=[0, 0], breaks=ticks, labels=labels)


    else:
        p += scale_x_continuous(expand=[0, 0])

    # -------------------------------------------------------------------------
    #
    # Adding geom_line
    #
    # -------------------------------------------------------------------------

    p += geom_line(size=line_width)

    # --------------------------------------------------------------------------
    #
    # Add facets
    #
    # --------------------------------------------------------------------------

    if facet_var is not None:
        # update facet_col
        facet_col = min(facet_col, len(dm[facet_var].unique()))
        p += facet_wrap("~ " + facet_var, ncol=facet_col)
    else:
        facet_col = 1

    message("facet_col " + str(facet_col))

    # -------------------------------------------------------------------------
    #
    # Adding rectangle overlays to highlight 5' and 3' regions
    #
    # -------------------------------------------------------------------------
    # Using the same dataframe as plotnine is still
    # very fragile.
    if config['ft_type'] in ["transcript", "user_regions"]:
        message("Highlighting upstream regions")

        dm['xmin_fr'] = [0 for x in range(dm.shape[0])]
        dm['xmax_fr'] = [bin_nb_ups / bin_nb_total * 100 for x in range(dm.shape[0])]
        dm['ymin_fr'] = [-np.inf for x in range(dm.shape[0])]
        dm['ymax_fr'] = [np.inf for x in range(dm.shape[0])]

        dm['xmin_to'] = [(bin_nb_total - bin_nb_dws) / bin_nb_total * 100 for x in range(dm.shape[0])]
        dm['xmax_to'] = [100 for x in range(dm.shape[0])]
        dm['ymin_to'] = [-np.inf for x in range(dm.shape[0])]
        dm['ymax_to'] = [np.inf for x in range(dm.shape[0])]

        if facet_var is None:
            dm_sub = dm.head(1)
        else:
            dm_sub = dm.drop_duplicates(subset=facet_var)

        if config['from']:
            p += geom_rect(data=dm_sub,
                           mapping=aes(xmin='xmin_fr',
                                       xmax='xmax_fr',
                                       ymin='ymin_fr',
                                       ymax='ymax_fr'),
                           fill='lightslategray',
                           color=None,
                           alpha=0.3,
                           show_legend=False)

        if config['to']:
            p += geom_rect(data=dm_sub,
                           mapping=aes(xmin='xmin_to',
                                       xmax='xmax_to',
                                       ymin='ymin_to',
                                       ymax='ymax_to'),
                           fill='lightslategray',
                           color=None,
                           alpha=0.3,
                           show_legend=False)

    # -------------------------------------------------------------------------
    #
    # Add number of observations per group
    #
    # -------------------------------------------------------------------------

    if show_group_number:
        if facet_var is None:
            p += geom_text(data=dm.drop_duplicates(subset=[group_by]),
                           mapping=aes(x='x',
                                       y='max_grp',
                                       label='nb_obs',
                                       hjust='hjust',
                                       color=group_by),
                           show_legend=False,
                           size=6,
                           ha='left')
        else:
            p += geom_text(data=dm.drop_duplicates(subset=[group_by, facet_var]),
                           mapping=aes(x='x',
                                       y='max_grp',
                                       label='nb_obs',
                                       hjust='hjust',
                                       color=group_by),
                           show_legend=False,
                           size=6,
                           ha='left')

    # --------------------------------------------------------------------------
    #
    # Apply colors
    #
    # --------------------------------------------------------------------------

    p += scale_color_manual(values=dict(list(zip(color_order, profile_colors))), name='Groups')

    # -------------------------------------------------------------------------
    # Turn warning off. Both pandas and plotnine use warnings for deprecated
    # functions. I need to turn they off although I'm not really satisfied with
    # this solution...
    # -------------------------------------------------------------------------

    def fxn():
        warnings.warn("deprecated", DeprecationWarning)

    # -------------------------------------------------------------------------
    #
    # Computing page size
    #
    # -------------------------------------------------------------------------

    if page_width is None:
        if config['ft_type'] in ["transcript", "user_regions"]:
            panel_width = 4
        else:
            panel_width = 3

        page_width = panel_width * facet_col

    if page_height is None:
        panel_height = 2
        if facet_var is None:
            page_height = panel_height
        else:
            panel_nb = len(dm[facet_var].unique()) / facet_col
            modulo = len(dm[facet_var].unique()) % facet_col
            if modulo:
                panel_nb += 1
            page_height = panel_nb * panel_height

    message("Page width set to " + str(page_width))
    message("Page height set to " + str(page_height))

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

        dm.to_csv(data_file, sep="\t", header=True, index=False)

    # Delete temporary dir
    shutil.rmtree(dir_name)


if __name__ == '__main__':
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    profile(**args)

else:

    test = '''
    #profile: prepare dataset
    @test "profile_1" {
     result=`gtftk get_example -d mini_real -f '*'; gtftk overlapping -i mini_real.gtf.gz -c hg38  -n > mini_real_noov.gtf; gtftk random_tx -i mini_real_noov.gtf  -m 1 -s 123 > mini_real_noov_rnd_tx.gtf`
      [ -s "hg38.genome" ]
    }

    @test "profile_2" {
     result=`gtftk get_example -d mini_real | gtftk select_by_key -t |  gtftk tabulate | gtftk col_from_tab -c transcript_id,gene_biotype | awk '$2=="protein_coding" || $2=="lincRNA" || $2=="transcribed_unprocessed_pseudogene"' > tx_classes.txt`
      [ -s "hg38.genome" ]
    }
                            
    #profile: prepare dataset
    @test "profile_3" {
     result=`gtftk mk_matrix -k 5 -i mini_real_noov_rnd_tx.gtf -d 5000 -u 5000 -w 200 -c hg38  -l  H3K4me3,H3K79me,H3K36me3 -y ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_promoter_pr`
      [ -s "mini_real_promoter_pr.zip" ]
    }

    #profile: test
    @test "profile_4" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -o profile_prom_1 -pf png -if example_01.png`
      [ -s "example_01.png" ]
    }


    #profile: make mini_real_promoter
    @test "profile_5" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -o profile_prom_1 -pf png -if example_01.png`
      [ -s "example_01.png" ]
    }

    #profile: make tss.bed
    @test "profile_6" {
     result=`gtftk select_by_key -i mini_real_noov_rnd_tx.gtf -k feature -v transcript |  gtftk get_5p_3p_coords > tss.bed`
      [ -s "tss.bed" ]
    }

    #profile: make single_nuc
    @test "profile_7" {
     result=`gtftk mk_matrix -u 5000 -d 5000 -i tss.bed -w 200 -l  H3K4me3,H3K79me,H3K36me3 -y ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_single_nuc_pr -c hg38 -t single_nuc`
      [ -s "mini_real_single_nuc_pr.zip" ]
    }

    #profile: test single_nuc
    @test "profile_8" {
     result=`gtftk profile -i mini_real_single_nuc_pr.zip -o profile_prom_1a -pf png -if example_01a.png`
      [ -s "example_01a.png" ]
    }

    #profile: make mini_real_tx
    @test "profile_9" {
     result=`gtftk mk_matrix  -k 5  -i mini_real_noov_rnd_tx.gtf -t transcript  -d 5000 -u 5000 -w 200 -c hg38  -l  H3K4me3,H3K79me,H3K36me3 -y ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tx_pr`
      [ -s "mini_real_tx_pr.zip" ]
    }

    #profile: make mini_real_tx
    @test "profile_10" {
     result=`gtftk select_by_key -i mini_real_noov_rnd_tx.gtf -k feature -v transcript | gtftk convert -f bed6 > mini_real_rnd_tx.bed`
      [ -s "mini_real_rnd_tx.bed" ]
    }

    #profile: test mini_real_tx
    @test "profile_11" {
     result=`gtftk profile -D -i mini_real_tx_pr.zip -o profile_tx_1 -pf png -if example_02.png`
      [ -s "example_02.png" ]
    }

    #profile: make mini_real_user_def
    @test "profile_12" {
     result=`gtftk mk_matrix  -k 5  --bin-around-frac 0.5 -i mini_real_rnd_tx.bed -t user_regions  -d 5000 -u 5000 -w 200 -c hg38  -l  H3K4me3,H3K79me,H3K36me3 -y ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_user_def`
      [ -s "mini_real_user_def.zip" ]
    }

    #profile: test mini_real_user_def
    @test "profile_13" {
     result=`gtftk profile -D -i mini_real_user_def.zip -o profile_udef_4  -pf png -if example_04.png`
      [ -s "example_04.png" ]
    }

    #profile: test mini_real_user_def
    @test "profile_14" {
     result=`gtftk profile -D -nm ranging -i mini_real_user_def.zip -o profile_udef_5  -pf png -if example_04b.png`
      [ -s "example_04b.png" ]
    }
        
    #profile: test mini_real_user_def
    @test "profile_15" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -g tx_classes -f bwig -o profile_prom_2  -ph 5  -pf png -if example_05.png`
      [ -s "example_05.png" ]
    }

    #profile: create dataset
    @test "profile_16" {
     result=`gtftk tabulate -k transcript_id,gene_biotype -i mini_real_noov_rnd_tx.gtf -H | sort | uniq | perl -ne 'print if (/(protein_coding)|(lincRNA)|(antisense)|(processed_transcript)/)'> tx_classes.txt`
      [ -s "tx_classes.txt" ]
    }

    #profile: create dataset
    @test "profile_17" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -g tx_classes -f bwig -o profile_prom_2  -ph 5  -pf png -if example_05.png`
      [ -s "example_05.png" ]
    }

    #profile: create dataset
    @test "profile_18" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -g bwig -f tx_classes  -o profile_prom_3  -ph 4 -c "#66C2A5,#FC8D62,#8DA0CB" -t tx_classes.txt  -pf png -if example_06.png`
      [ -s "example_06.png" ]
    }

    #profile: create dataset
    @test "profile_19" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -g tx_classes -f bwig  -o profile_prom_4  -ph 4 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -t tx_classes.txt  -pf png -if example_07.png`
      [ -s "example_07.png" ]
    }

    #profile: create dataset
    @test "profile_20" {
     result=`gtftk mk_matrix  -k 5  --bin-around-frac 0.5 -i mini_real_noov_rnd_tx.gtf -t transcript  -d 5000 -u 5000 -w 200 -c hg38  -l  H3K4me3,H3K79me,H3K36me3 -y ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tx_pr_2`
      [ -s "mini_real_tx_pr_2.zip" ]
    }
    
    #profile: create dataset
    @test "profile_21" {
     result=`gtftk profile -D -i mini_real_tx_pr_2.zip -g tx_classes -f bwig  -o profile_tx_3 -pw 12  -ph 7 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -t tx_classes.txt  -pf png -if example_08.png`
      [ -s "example_08.png" ]
    }

    #profile: create dataset
    @test "profile_22" {
     result=`gtftk profile -D -i mini_real_promoter_pr.zip -g bwig -f chrom  -o profile_prom_5  -ph 15 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF"   -pf png -if example_09.png`
      [ -s "example_09.png" ]
    }
     
    #profile: create dataset
    @test "profile_23" {
     result=`gtftk profile -th classic -D -i mini_real_promoter_pr.zip -g bwig -f chrom  -o profile_prom_5  -ph 15 -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF"   -pf png -if example_09b.png`
      [ -s "example_09b.png" ]
    }
    
        
            '''

    cmd = CmdObject(name="profile",
                    message="Create coverage profile using a bigWig as input.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    updated=__updated__,
                    notes=__notes__,
                    references=__references__,
                    group="coverage",
                    test=test)
