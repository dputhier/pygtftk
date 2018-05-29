"""
    The ``control_list`` plugin
    ============================

    Based on a reference gene list, returns a list of genes matched for
    signal. The expression/signal values should be provided for all genes through
    the in_file argument.
"""

import argparse
import os

from pygtftk.arg_formatter import FileWithExtension
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import chomp
from pygtftk.utils import make_outdir_and_file
from pygtftk.utils import message

R_LIB = 'beanplot,ggplot2,reshape2'

__updated__ = "2018-01-20"
__doc__ = """
 Based on a reference gene list (or more generally IDs) this command tries to extract a set of
 other genes/IDs matched for signal/expression. The --reference-gene-file contains
 the list of reference IDs while the -\inputfile contains a tuple gene/signal for all genes.
"""

__notes__ = """
    -- -\infile is a two columns tabulated file. The first column contains the list of ids (including reference IDs)
    and the second column contains the expression/signal values. This file should contain no header.
    -- Think about discarding any unwanted IDs from -\infile before calling control_list.
"""


def make_parser():
    """The program parser."""

    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('--in-file', '-i',
                            metavar='TXT',
                            help='A two columns tab-file. See notes.',
                            default=None,
                            type=FileWithExtension('r',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]')),
                            required=True)

    parser_grp.add_argument('--referenceGeneFile', '-r',
                            metavar='TXT',
                            help='The file containing the reference gene list (1 column,'
                                 ' transcript ids).'
                                 ' No header.',
                            default=None,
                            type=FileWithExtension('r',
                                                   valid_extensions=('\.[Tt][Xx][Tt]',
                                                                     '\.[Cc][Ss][Vv]',
                                                                     '\.[Tt][Aa][Bb]',
                                                                     '\.[Tt][Ss][Vv]',
                                                                     '\.[Ii][Dd][Ss]')),
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

    return parser_grp


def control_list(in_file=None,
                 out_dir=None,
                 referenceGeneFile=None,
                 log2=False,
                 pseudo_count=1,
                 tmp_dir=None,
                 logger_file=None,
                 verbosity=None):
    for p, line in enumerate(in_file):

        line = chomp(line)
        line = line.split("\t")

        try:
            fl = float(line[1])

            if log2:
                fl = fl + pseudo_count
                if fl <= 0:
                    message("Can not log transform negative/zero values.",
                            type="ERROR")

        except:
            msg = "It seems that column 2 of input file"
            msg += " contains non numeric values."
            msg += "Check that no header is present and that "
            msg += "columns are ordered properly."
            message(msg, type="ERROR")

    message("Selecting genes (R call).")

    # Preparing pdf file name
    file_out_list = make_outdir_and_file(out_dir, ["control_list.txt",
                                                   "reference_list.txt",
                                                   "diagnostic_diagrams.pdf",
                                                   "R_code_control_list.R",
                                                   "command.log"],
                                         force=True)

    control_file, reference_file_out, pdf_file, r_code_file, log_file = file_out_list

    code_body = """

            ########################################
            # Load library
            ########################################

            suppressWarnings(suppressMessages(library(beanplot)))
            suppressWarnings(suppressMessages(library(ggplot2)))
            suppressWarnings(suppressMessages(library(reshape2)))

            ########################################
            # Function declaration
            ########################################

            message <- function(msg){{
                cat(paste("    |--- ", msg, "\\n", sep=""), file=stderr())
            }}

            ########################################
            ## Get the list of reference genes
            ########################################
            
            reference_genes <- as.character(read.table('{reference_file}',
                                            header=F)[,1])

            ########################################
            ## Get expression data
            ########################################
            
            exp_data <- read.table('{signal_file}', sep="\\t", head =F)
            
            exp_data_vec <- exp_data[,2] + {pseudo_count}
            
            ########################################
            ## Log transformation
            ########################################
            
            to.log <- '{log2}'
            
            if(to.log == 'True'){{
                if(length(exp_data_vec[exp_data_vec == 0])){{
                    message("Can't use log transformation on zero or negative values. Use -p. Exiting.")
                    q("no", 1, FALSE)
                }}
                exp_data_vec <- log2(exp_data_vec)
            }}
            
            names(exp_data_vec) <- exp_data[,1]

            # Now we have sorted expression data with gene/tx names.
            exp_data_vec <- sort(exp_data_vec)

            # T/F vector indicating which in the
            # expression data list are found in reference_gene
            which_reference_genes <- names(exp_data_vec) %in% reference_genes

            # convert the T/F vector to positions indicating wich position in
            # expression data is a reference gene/tx

            which_reference_genes <- which(which_reference_genes)

            message(paste("Found ", length(which_reference_genes),
                    " genes of the reference in the provided signal file", sep=""))
            
            not_found <- !(reference_genes %in% names(exp_data_vec))
            if(length(reference_genes[not_found]) > 0){{
                message(paste("List of reference genes not found :", reference_genes[not_found]))
            }}else{{
                message("All reference genes were found.")
            }}

            ########################################
            ## Search for gene with matched signal
            ########################################
            
            control_list<- c()

            nb.candidate.left <- length(exp_data_vec) -  length(which_reference_genes)
            
            if(nb.candidate.left < length(which_reference_genes) ){{
                message("Not enough element to perform selection. Exiting")
                q("no", 1, FALSE)
            }}

            cpt <- 1
            
            candidates <- exp_data_vec

            for(i in which_reference_genes){{
                p <- i
                not_candidate_pos <- unique(c(which_reference_genes, control_list))
                candidates[not_candidate_pos] <- NA
                diff <- abs(exp_data_vec[p] - candidates)
                control_list[cpt] <- which.min(diff)
                cpt <- cpt + 1
                
            }}


            write.table(cbind(exp_data_vec[which_reference_genes]),
                    "{reference_file_out}",
                    sep="\t",
                    quote=F,
                    col.names=NA)

            write.table(cbind(exp_data_vec[control_list]),
                    "{control_file}",
                    sep="\t",
                    quote=F,
                    col.names=NA)
                    

            message("Preparing diagnostic plots.")
            m <- as.data.frame(cbind(
                                exp_data_vec[control_list],
                                exp_data_vec[which_reference_genes]
                )
            )
                            
            colnames(m) <- c("Control", "Reference")
            
            
            pdf("{pdf_file}")
            
            # Preparing beanplot (side=both)
            beanplot(m,
                    col = list("blue", "darkgrey"),
                    border="white",
                    log="",
                    ll=0.13,
                    what=c(0,1,0,1),
                    side="both"
            )
            grid(col="grey", nx=0, ny=NULL)
            beanplot(m,
                    col = list("blue", "darkgrey"),
                    border="white",
                    log="", main="Beanplots of Control and Reference values",
                    ll=0.13,
                    what=c(0,1,0,1),
                    add=T,
                    side="both"
            )
            
            # Preparing beanplot (side=default)
            beanplot(m,
                    col = list("blue", "darkgrey"),
                    border="white",
                    log="",
                    ll=0.13,
                    what=c(0,1,0,1)
            )
            grid(col="grey", nx=0, ny=NULL)
            beanplot(    m,
                    col = list("blue", "darkgrey"),
                    border="white",
                    log="", main="Beanplots of Control and Reference values",
                    ll=0.13,
                    what=c(0,1,0,1),
                    add=T
            )
            
            # Preparing boxplot
            boxplot(m, col=c("blue","darkgrey"))
            grid(col="grey", nx=0, ny=NULL)

            # Preparing qqplot
            plot(    sort(m$Control),
                    sort(m$Reference),
                    pch=16, xlab="Control",
                    ylab="Reference",
                    main="QQplot of control and reference values",
                    panel.first=grid()
            )
            
            # Preparing histograms
            
            m.m <- melt(m, id.vars = c(NULL))
            col <- m.m$variable
            levels(col) <- c("blue","darkgrey")
            p <- ggplot(m.m,
                        aes(    x=value,
                                color=variable,
                                fill=variable))
            p <- p + geom_bar(position="dodge")
            p <- p + scale_fill_manual(values=levels(col))
            p <- p + scale_color_manual(values=levels(col))
            p <- p + theme(legend.title=element_blank())
            
            suppressMessages(print(p))
            
            out <- dev.off()
    """.format(signal_file=in_file.name,
               reference_file=referenceGeneFile.name,
               log2=log2,
               pseudo_count=str(pseudo_count),
               pdf_file=pdf_file.name,
               control_file=control_file.name,
               reference_file_out=reference_file_out.name)

    message("Printing R code to: " + r_code_file.name,
            type="DEBUG")

    r_code_file.write(code_body)
    r_code_file.close()

    message("Executing R code.", type="DEBUG")

    # Execute R code.
    os.system("cat " + r_code_file.name + "| R --slave")


if __name__ == '__main__':

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    control_list(**args)

else:

    test = '''
    #control_list
    @test "control_list_1" {
      result=`gtftk control_list -i pygtftk/data/control_list/control_list_data.txt -r pygtftk/data/control_list/control_list_reference.txt -D ; cat control_list/control_list.txt | cut -f2| perl -npe 's/\\n/,/'`
      [ "$result" = "V1,2.02,4.04,6.06," ]
    }
    
    #control_list
    @test "control_list_2" {
      result=` rm -Rf control_list`
      [ "$result" = "" ]
    }
            
    '''

    cmd = CmdObject(name="control_list",
                    message="Returns a list of gene matched for expression based on reference values.",
                    parser=make_parser(),
                    fun=control_list,
                    desc=__doc__,
                    updated=__updated__,
                    notes=__notes__,
                    group="miscellaneous",
                    test=test,
                    rlib=R_LIB)
