#!/usr/bin/env python
import argparse
import os
import shutil
import tempfile
import zipfile

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.arg_formatter import float_greater_than_null
from pygtftk.arg_formatter import float_grt_than_null_and_lwr_than_one
from pygtftk.arg_formatter import int_greater_than_null
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import check_r_packages, make_tmp_file
from pygtftk.utils import chomp
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import message

R_LIB = 'ggplot2,reshape2,grid,data.table,plyr'

__updated__ = "2018-01-20"
__doc__ = """
 Produces bigWig coverage profiles using calls to ggplot2.
"""

__notes__ = """
 -- The ranging nomalization method [1] implies the following transformation: 
 -- -  (x_i - min(x))/(max(x) - min(x)).
 -- Think about using normalized bigWig files as input to mk_matrix. This
 will limit the requirement for an additional normalization step (see
 Deeptools for a set of useful methods implemented in bamCoverage/bamCompare).
""".format(r=R_LIB)

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
                            choices=["mean", "median", "sd", "mad", "IQR"],
                            type=str,
                            required=False)

    parser_grp.add_argument('-c', '--profile-colors',
                            help='Set of colors for profiles.',
                            default="#1b9e77,#d95f02,#7570b3,#e7298a,#66a61e,#e6ab02,#a6761d,#666666",
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

    parser_grp.add_argument('-f', '--facet',
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

    parser_grp.add_argument('-th', '--theme',
                            choices=[
                                'bw',
                                'grey',
                                'gray',
                                'linedraw',
                                'light',
                                'dark',
                                'minimal',
                                'classic',
                                'void',
                                'test'],
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
                 page_height=None,
                 page_format='pdf',
                 user_img_file=None,
                 tmp_dir=None,
                 facet_col=None,
                 strip_color="#707070",
                 force_tx_class=False,
                 stat="mean",
                 facet=None,
                 xlab="Selected genomic regions",
                 axis_text=8,
                 strip_text=8,
                 line_width=1,
                 theme='bw',
                 logger_file=None,
                 verbosity=False
                 ):
    # -------------------------------------------------------------------------
    #
    # The selected theme
    #
    # -------------------------------------------------------------------------

    theme = 'theme_' + theme

    # -------------------------------------------------------------------------
    #
    # facet and group_by should be different
    #
    # -------------------------------------------------------------------------

    if facet == group_by:
        message("--facet and --group-by should be different.",
                type="ERROR")

    # -------------------------------------------------------------------------
    #
    # We need some R packages
    #
    # -------------------------------------------------------------------------

    message("Checking R package.")
    check_r_packages(R_LIB.split(","))

    # -------------------------------------------------------------------------
    #
    # Convert some args
    #
    # -------------------------------------------------------------------------

    if to_log:
        to_log = "T"
    else:
        to_log = "F"

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

    # Use a temp file to avoid concurrency issues
    dir_name = tempfile.mkdtemp(prefix='GTFtk_matrix_')
    message("Uncompressing in directory: " + dir_name,
            type="DEBUG")

    try:
        with zipfile.ZipFile(inputfile.name) as zf:
            zf.extractall(dir_name)
    except:
        message("Problem encountered when unzipping...",
                type="ERROR")

    inputfile_main = open(os.path.join(dir_name, zf.namelist()[0]), "r")
    message("Reading from file: " + inputfile_main.name,
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
    # Convert some args
    #
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    #
    # Check arguments: --facet, --group-by according to the number of bigwig
    #
    # -------------------------------------------------------------------------
    # If one is analyzing more than one bigwig
    # the "bwig" factor should appear in  --facet or --group-by
    if len(input_file_bwig) > 1:
        if facet != "bwig":
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

    profile_colors = profile_colors.split(",")

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

    if user_img_file is None:
        file_out_list = make_outdir_and_file(out_dir,
                                             [img_file, "R_diagram_code.R"],
                                             force=True)

        img_file, r_code_file = file_out_list

    else:
        file_out_list = make_outdir_and_file(out_dir,
                                             ["R_diagram_code.R"],
                                             force=True)

        r_code_file = file_out_list[0]
        img_file = user_img_file
        if not img_file.name.endswith(page_format):
            msg = "Image format: {f}. Please fix.".format(f=page_format)
            message(msg, type="ERROR")
        test_path = os.path.abspath(img_file.name)
        test_path = os.path.dirname(test_path)

        if not os.path.exists(test_path):
            os.makedirs(test_path)

    # ------------------------------------------------------------------
    # Graphics with a call to R
    # ------------------------------------------------------------------

    r_code = """
    # Input variables
    # ==============================
    
    None <- 'None'
    ft_type <- '{ft_type}'
    from <- {fr}
    to <- {to}
    stat <- {stat}
    groups <- '{group_by}'
    facet <- '{facet}'
    group_order <- '{color_order}'
    group_order <- strsplit(group_order, ",")[[1]]
    inputfile <- '{inputfile}'
    profile_color <- unlist(strsplit("{profile_color}", ","))[1:length(group_order)]
    axis_text <- {axis_text}
    normalization.method <- '{normalization_method}'
    upper.limit <- {upper_limit}
    to.log <- {to_log}
    transcript_file <- '{transcript_file}'
    img_file <- '{img_file}'
    page_width <- {page_width}
    page_height <- {page_height}
    xlab<- '{xlab}'
    facet_col <- {facet_col}


    
    # Load libraries and declare functions
    # =====================================
    
    suppressWarnings(suppressMessages(library(ggplot2)))
    suppressWarnings(suppressMessages(library(reshape2)))
    suppressWarnings(suppressMessages(library(data.table)))
    suppressWarnings(suppressMessages(library(grid)))
    suppressWarnings(suppressMessages(library(plyr)))
    
    # Gtftk-like messages
    message <- function(msg){{
      cat(paste("    |--- ",format(Sys.time(),"%R"), "-INFO : ", msg, "\n", sep=""))
    }}
    

    # Return a themed ggplot object
    #-------------------------------
    
    get_ggplot <- function(p,
                           axis.text=axis_text,
                           df=NULL,
                           title=""){{
    
      message("Theming and ordering. Please be patient...")
    
      p <- p + theme(legend.title=element_blank())
      p <- p + {theme}()

      p <- p + theme(legend.text=element_text(size=8))
      p <- p + theme(legend.title=element_blank())
      p <- p + theme(legend.position = "bottom",
                     legend.key = element_rect(colour = "white"))
      p <- p + theme(axis.text=element_text(size=axis.text))
      p <- p + guides(col = guide_legend(ncol=5))
      p <- p + theme(axis.text.x = element_blank())
      p <- p + theme(axis.text.x = element_text(angle = 40, hjust = 1))
      p <- p + theme( legend.text=element_text(size=8),
                      legend.position = "bottom",
                      legend.key = element_rect(colour = "white"),
                      strip.text.x = element_text(size = {strip_text}, colour = 'white'),
                      strip.background = element_rect(fill="{strip_color}"))
      
      #p <- p + theme(axis.title.x = element_text(margin = margin(t = 20)))

      
      if(ft_type %in% c("transcript","user_regions")){{

         if(from){{
          
          if(to){{
            breaks <- c(0, bin_nb_ups/2, seq(bin_nb_ups,
                        bin_nb_main + bin_nb_ups, length.out=11),
                        bin_nb_total - bin_nb_dws/2, bin_nb_total) / bin_nb_total * 100
            labels <- c(- from,  round(- from/2,0),
                        paste(seq(0,100, length.out=11), "%", sep=""), round(to/2,0), to)
          }}else{{
            breaks <- c(0, bin_nb_ups/2,
                        seq(bin_nb_ups,  bin_nb_total, length.out=11)) / bin_nb_total * 100
            labels <- c(- from,
                        round(- from/2,0), paste(seq(0,100, length.out=6), "%", sep=""))
          }}
          
        }}else{{
          
          if(to){{
            breaks <- c(seq(0, bin_nb_main, length.out=11),
                        bin_nb_total - bin_nb_dws/2, bin_nb_total) / bin_nb_total * 100
            labels <- c(paste(seq(0,100, length.out=6), "%", sep=""), to/2, to)
          }}else{{
            breaks <- seq(from=0, to=bin_nb_total, length.out=6) / bin_nb_total * 100
            labels <- paste(seq(0,100, length.out=6), "%", sep="")
          }}
          
        }}
        
        p <- p + scale_x_continuous(expand=c(0,0), breaks=breaks, labels=labels)

        if(from){{

          rectangles <- data.frame(xmin = c(0),
                                   xmax = c((bin_nb_ups)/ bin_nb_total * 100 ),
                                   ymin = -Inf,
                                   ymax = Inf)

          p <- p + geom_rect(data=rectangles,
                           aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                           fill='lightslategray', alpha=0.3, inherit.aes = FALSE,
                           show.legend=F)
        }}

        if(to){{
          rectangles <- data.frame(xmin = c(bin_nb_total - bin_nb_dws) / bin_nb_total * 100,
                                   xmax = c(100),
                                   ymin = -Inf,
                                   ymax = Inf)
          p <- p + geom_rect(data=rectangles,
                             aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                             fill='lightslategray', alpha=0.3, inherit.aes = FALSE,
                             show.legend=F)
        }}
          
      }}else{{
        p <- p + scale_x_continuous(expand=c(0,0))
      }}
      
      p <- p + scale_color_manual(values=profile_color)
      
      p <- p + guides(col = guide_legend(ncol=5))
      p <- p + theme(panel.spacing.x = unit(1, "lines")) + ggtitle(title)
      return(p)
    }}
    
    # Compute simple overlayed profile
    # ---------------------------------
    
    # mean computation won't take NA into account
    compute_stat <- function(x) apply(x[,5:ncol(x)], 2, stat, na.rm = TRUE)
    
    
    simple_profile <- function(groups='bwig',
                               pos_order=NULL,
                               group_order=NULL,
                               title='',
                               ylab='Signal',
                               page_width=page_width,
                               page_height=page_height){{

      d_split <- split(d_merge, as.factor(as.character(d_merge[,groups])))
      d_split_mean <- data.frame(lapply(d_split, compute_stat), check.names = FALSE)
      d_split_mean_melt <- suppressWarnings(suppressMessages(melt(data.frame(d_split_mean,pos=rownames(d_split_mean), check.names = FALSE))))
      
      
      if(ft_type %in% c("transcript","user_regions")){{
        d_split_mean_melt$pos <- factor(d_split_mean_melt$pos,
                                        levels=pos_order,
                                        ordered=T)
      }}else{{
        d_split_mean_melt$pos <- factor(d_split_mean_melt$pos,
                                        levels=pos_order,
                                        ordered=T)
        levels(d_split_mean_melt$pos) <- seq(from= - from, to=to, length.out = length(pos_order))
        d_split_mean_melt$pos <- as.numeric(as.character(d_split_mean_melt$pos))
      }}
      
      d_split_mean_melt$variable <- factor(d_split_mean_melt$variable,
                                           levels=group_order,
                                           ordered=T)
                                           
      if(ft_type %in% c("transcript","user_regions")){{
        #continuous scale
        levels(d_split_mean_melt$pos) <- seq(0, 100, length.out=bin_nb_total )
        d_split_mean_melt$pos <- as.numeric(as.character(d_split_mean_melt$pos))
      }}

      p <- ggplot(d_split_mean_melt, aes(x=pos, y=value, colour = variable)) + xlab(xlab)
      p <- ggplot(d_split_mean_melt, aes(x=pos, y=value, colour = variable)) + xlab(xlab)
      p <- p + geom_line(aes(group = variable), size={line_width})
      p <- p + ylab(ylab)

      message(paste("Page height: ",  page_height, collapse=""))
      message(paste("Page width: ",  page_width, collapse=""))
      
      ggsave(filename=img_file,
             plot=get_ggplot(p, df=d_split_mean_melt, title=title),
             width=page_width,
             height=page_height)
      
      return(d_split_mean_melt)
      
    }}
    
    
    # Compute facetted overlayed profile
    # -----------------------------------
    
    facet_profile <- function(group_by="bwig",
                              split_by_1="tx_classes",
                              group_by_order=NULL,
                              split_by_1_order=group_order,
                              x_pos_order=pos_order,
                              title="",
                              ylab="Signal",
                              page_width=page_width,
                              page_height=page_height){{
      
        d_split <- split(d_merge, as.factor(as.character(d_merge[,group_by])))
          
        for(i in names(d_split)){{
            d_split[[i]] <- split(d_split[[i]], as.factor(as.character(d_split[[i]][,split_by_1])))
        }}
        
        new_names <- paste(names(d_split[[1]]), " (#", lapply(d_split[[1]],nrow), ")", sep="")
        names(new_names) <- names(d_split[[1]])
        
        
        for(i in names(d_split)){{
            d_split[[i]] <-  data.frame(lapply(d_split[[i]], compute_stat), check.names = FALSE)
        }}
        
        for(i in names(d_split)){{
            d_split[[i]] <-  suppressWarnings(suppressMessages(melt(data.frame(d_split[[i]], pos=rownames(d_split[[i]]), group_by=i, check.names = FALSE))))
        }}
        
        d_split_mean_melt <- do.call("rbind", d_split)
        
        if(ft_type %in% c("transcript","user_regions")){{
            d_split_mean_melt$pos <- factor(d_split_mean_melt$pos,
                                            levels=pos_order,
                                            ordered=T)
        }}else{{
            d_split_mean_melt$pos <- factor(d_split_mean_melt$pos,
                                            levels=pos_order,
                                            ordered=T)
            levels(d_split_mean_melt$pos) <- seq(from= - from, to=to, length.out = length(pos_order))
            d_split_mean_melt$pos <- as.numeric(as.character(d_split_mean_melt$pos))
        }}
        
        
        d_split_mean_melt$group_by <- factor(d_split_mean_melt$group_by,
                                             levels=group_by_order,
                                             ordered=T)
        
        index <- grep("variable", colnames(d_split_mean_melt))
        colnames(d_split_mean_melt)[index] <- "split_by_1"
          
        d_split_mean_melt$split_by_1 <- factor(new_names[d_split_mean_melt$split_by_1],
                                                 levels=new_names[split_by_1_order],
                                                 ordered=T)
        if(ft_type %in% c("transcript","user_regions")){{
             levels(d_split_mean_melt$pos) <- seq(0, 100, length.out=bin_nb_total )
             d_split_mean_melt$pos <- as.numeric(as.character(d_split_mean_melt$pos))
        }}
        p <- ggplot(d_split_mean_melt, aes(x=pos, y=value, colour = group_by)) + xlab(xlab)
        p <- p + ylab(ylab)
        p <- p + geom_line(aes(group = group_by), size={line_width}) + facet_wrap(~ split_by_1, ncol = facet_col)
        
      message(paste("Page height: ",  page_height, collapse=""))
      message(paste("Page width: ",  page_width, collapse=""))

      ggsave(filename=img_file,
             plot=get_ggplot(p, df=d_split_mean_melt, title=title),
             width=page_width,
             height=page_height)

        return(d_split_mean_melt)
      
    }}
    
    # Load dataset
    #================
    
    message("Preparing profile diagram.")
    message("Reading input files.")
    
    
    d <- as.data.frame(fread(inputfile,
                            sep='\\t',
                            header=T,
                            skip=1,
                            showProgress=FALSE))
    # delete start and end columns
    pos_start_end <- grep("(start)|(end)", colnames(d), perl=T)
    d <- d[, - pos_start_end]
    
    # Merge upstream and downstream bins
    # -----------------------------------
    # (deprecated)
    d_merge <- d
    
    
    # Get the transcript classes
    # ---------------------------
    df_tx_class <- read.table(transcript_file, sep='\\t', head=F, colClasses = "character")
    tx <- as.character(df_tx_class[,1])
    tx_classes <- as.character(df_tx_class[,2])
    names(tx_classes)  <- tx
    
    # Select the transcript of interest and add class info to the data.frame
    # -----------------------------------------------------------------------
    all_tx <- as.character(d_merge$gene)
    d_merge <- d_merge[all_tx %in% tx,]
    all_tx <- as.character(d_merge$gene)
    d_merge <- cbind(tx_classes[all_tx], d_merge)
    colnames(d_merge)[1] <- "tx_classes"
    
    # Save some columns
    # ------------------
    gene.info <- d_merge[,c(5)]
    d_merge <- d_merge[, -c(5)]
    
    # Store level order (pos)
    # ------------------------
    pos_order <- colnames(d_merge)[5:ncol(d_merge)]
    
    # compute bin number upstream, main, downstream
    # -----------------------------------------------
    bin_nb_main <- length(grep("main", colnames(d_merge)))
    bin_nb_ups <- length(grep("upstream", colnames(d_merge)))
    bin_nb_dws <- length(grep("downstream", colnames(d_merge)))
    bin_nb_total <- bin_nb_ups + bin_nb_main + bin_nb_dws


    # ceiling
    # -----------------------------------------------

    num_cols <- 5:ncol(d_merge)
    
    if(upper.limit < 1){{
        message('Ceiling')
        for(k in unique(d_merge$bwig)){{
            tmp <- d_merge[d_merge$bwig == k, num_cols]
            qu <- quantile(unique(as.vector(as.matrix(tmp))), upper.limit)
            tmp[tmp > qu] <- qu
            d_merge[d_merge$bwig == k, num_cols] <- tmp
        }}
    }}


    # Normalize/transform
    # -----------------------------------------------

    if(to.log){{
      tmp <- d_merge[,num_cols]
      tmp <- tmp[tmp == 0]
      if(length(tmp) > 0){{
        message("Zero value detected. Adding a pseudocount (+1) before log transformation.")
        d_merge[,num_cols] <- d_merge[,num_cols] + 1
      }}
      
      d_merge[,num_cols] <- log2(d_merge[,num_cols])
      ylab <- "log2(Signal)"
      
    }} else{{
       ylab <- "Signal"
    }}

    
    if(normalization.method == 'ranging'){{
        message('Normalizing (percentage)')
        for(k in unique(d_merge$bwig)){{
            tmp <- d_merge[d_merge$bwig == k, num_cols]
            tmp.norm <- (tmp - min(tmp))/(max(tmp) - min(tmp)) * 100
            d_merge[d_merge$bwig == k, num_cols] <- tmp.norm
        }}
        
        ylab <- paste("scaled(", ylab, ", %)", sep="")
    }}
    

    # Main calls
    # =========================
    

    # Output
    # -------------------
    
    if(facet == 'None'){{
        
        if(page_width == 'None'){{

            page_width <- 2
        }}

        if(facet_col == 'None'){{
            facet_col <- 4
        }}


        if(page_height == 'None'){{
            page_height <- 4
        }}
        tmp <- simple_profile(groups=groups,
                          pos_order=pos_order,
                          group_order=group_order,
                          title="stat: '{stat}', group: '{group_by}'\n ceiling: {upper_limit}, log: {to_log}, norm: '{normalization_method}'",
                          ylab=ylab,
                          page_width=page_width,
                          page_height=page_height)
    }}else{{
    

        if(page_width == 'None'){{

            if (facet_col != 'None'){{
                message("Inside")
                page_width <- 2 * facet_col
            }}else{{
                
                n <- length(sort(unique(d_merge[,"{facet}"])))
                page_width <- 2 * n
                message(paste("Page height: ",  page_height, collapse=""))
                message(paste("Page width: ",  page_width, collapse=""))
            }}
                
        }}

        if(page_height == 'None'){{
            if (facet_col != 'None'){{
                n <- length(sort(unique(d_merge[,"{facet}"])))
                page_height <- ceiling(n / facet_col) * 4
            }}else{{
                page_height <- 4
            }}
        }}

        if(facet_col == 'None'){{
            facet_col <- 4
        }}
                  
        tmp <- facet_profile(group_by=groups,
                             split_by_1="{facet}",
                             group_by_order=group_order,
                             split_by_1_order=sort(unique(d_merge[,"{facet}"])),
                             x_pos_order=pos_order,
                             title="stat: {stat}, group: {group_by}, facets: {facet}.",
                             ylab=ylab,
                             page_width=page_width,
                             page_height=page_height)
              
    }}
    
    message("Plots saved..")
    
    """.format(ft_type=config['ft_type'],
               fr=config['from'],
               to=config['to'],
               inputfile=inputfile_main.name,
               img_file=img_file.name,
               transcript_file=transcript_file.name,
               group_by=group_by,
               color_order=color_order,
               profile_color=",".join(profile_colors),
               strip_color=strip_color,
               stat=stat,
               facet=facet,
               xlab=xlab,
               page_width=page_width,
               strip_text=strip_text,
               axis_text=axis_text,
               facet_col=facet_col,
               line_width=line_width,
               upper_limit=upper_limit,
               normalization_method=normalization_method,
               to_log=to_log,
               theme=theme,
               page_height=page_height)

    message("Printing R code to: " + r_code_file.name)

    r_code_file.write(r_code)
    r_code_file.close()

    message("Executing R code.")

    # Execute R code.
    os.system("cat " + r_code_file.name + "| R --slave")

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
