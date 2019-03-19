"""
A module to compute bigwig coverage over a set of regions (bed).
"""

import multiprocessing
import os
import pyBigWig
import sys
from builtins import range
from builtins import str
from builtins import zip
from itertools import repeat
from tempfile import NamedTemporaryFile

import numpy as np
from pybedtools import BedTool

import pygtftk
from pygtftk.utils import GTFtkError
from pygtftk.utils import add_prefix_to_file
from pygtftk.utils import close_properly
from pygtftk.utils import flatten_list
from pygtftk.utils import intervals
from pygtftk.utils import make_tmp_file
from pygtftk.utils import message

# -------------------------------------------------------------------------
# TMP_FILE_POOL_MANAGER stores temporary file name
# make_tmp_file_pool is function that add temporary files to TMP_FILE_POOL_MANAGER
# TMP_FILE_POOL_MANAGER will be updated by workers (in contrast to a global
# variable)
# -------------------------------------------------------------------------

TMP_FILE_POOL_MANAGER = multiprocessing.Manager().list()


def make_tmp_file_pool(prefix='tmp',
                       suffix='',
                       store=True,
                       dir=None):
    """
    This

    :Example:

    >>> from pygtftk.utils import make_tmp_file_pool
    >>> tmp_file = make_tmp_file_pool()
    >>> assert os.path.exists(tmp_file.name)
    >>> tmp_file = make_tmp_file_pool(prefix="pref")
    >>> assert os.path.exists(tmp_file.name)
    >>> tmp_file = make_tmp_file_pool(suffix="suf")
    >>> assert os.path.exists(tmp_file.name)

    """

    dir_target = None

    if dir is None:
        if pygtftk.utils.TMP_DIR is not None:
            if not os.path.exists(pygtftk.utils.TMP_DIR):
                msg = "Creating directory {d}."
                message(msg.format(d=pygtftk.utils.TMP_DIR), type="INFO")
                os.mkdir(pygtftk.utils.TMP_DIR)
                dir_target = pygtftk.utils.TMP_DIR

    else:
        if not os.path.exists(dir):
            msg = "Creating directory {d}."
            message(msg.format(d=dir), type="INFO")
            os.mkdir(dir)
            dir_target = dir

    tmp_file = NamedTemporaryFile(delete=False,
                                  mode='w',
                                  prefix=prefix + "_pygtftk_",
                                  suffix=suffix,
                                  dir=dir_target)

    if store:
        TMP_FILE_POOL_MANAGER.append(tmp_file.name)

    return tmp_file


# -------------------------------------------------------------------------
# Define the job of a single worker whose job is to compute mean coverage
# or coverage profile using a bed and a bigwig as input.
# -------------------------------------------------------------------------


def _big_wig_coverage_worker(input_values):
    """
    This function compute bigwig coverage. The input_values arguments is a
    tuple that contains various input parameters. 'span' is a tuple that
    correspond to a fraction (from, to) of the bedfile to be processed. Each
    worker will process all bigwig filesbut it will only process a fraction
    (span) of the bed file regions


    :param span: the fraction (lines) of the bed file [from, to] to be processed.
    :param bw_list: the list of bigWig files to be processed.
    :param region_bed_file_name: the bed file containing the region for which coverage is to be computed.
    :param bin_nb: The number of bin into which the region should be splitted.
    If the number of nucleotides is < nbBin a warning is printed.
    :param pseudo_count: A value for the pseudo_count.
    :param n_highest: compute the score based on the n highest values in the bins.
    :param profile: compute coverage profile not a single coverage value (mean).
    :param stranded: controls whether the profile should be ordered based on
    strand.
    :param type: This string will be added to the output to indicate the type
    of region (e.g tss, promoter...).
    :param label: Bigwig labels (i.e short name version)
    :param zero_to_na: Use NA not zero when region is undefined in bigwig.
    :param stat: mean (default) or sum.
    :param verbose: run in verbose mode.

    """

    (span, bw_list,
     region_bed_file_name,
     bin_nb, pseudo_count,
     n_highest, profile,
     stranded, type, label, zero_to_na, stat, verbose) = input_values

    pc = pseudo_count

    if not profile:
        if n_highest is None:
            n_highest = bin_nb
        results = list()
    else:
        if bin_nb < 1:
            bin_nb = 1
        matrix_file = make_tmp_file_pool(prefix="worker_coverage_", suffix=".txt")

    for cpt, big_wig in enumerate(bw_list):

        try:
            bigwig = pyBigWig.open(big_wig)
            if not bigwig.isBigWig():
                message("Not a bigwig file :" + big_wig, type="ERROR")
        except:
            message("Not a bigwig file :" + big_wig, type="ERROR")

        mesg = "Computing coverage for %s (chunks : #%s , type : %s, lab : %s)."
        mesg = mesg % (os.path.basename(big_wig), str(span[1] - span[0]), type, label[cpt])
        message(mesg, type="INFO")

        # Load the regions for which the coverage is to be processed.

        tx_bed = BedTool(region_bed_file_name)

        # The fraction of bed file
        # to be processed
        (from_here, to_here) = span

        nb = 0
        nb_to_do = to_here - from_here

        for i in tx_bed[slice(from_here, to_here)]:

            nb += 1

            if nb == nb_to_do:
                p_name = str(multiprocessing.current_process().name)
                message(p_name + " has processed " + str(nb) + " regions")

            if (i.end - i.start) < bin_nb:

                if pygtftk.utils.WARN_REGION_SIZE:
                    pygtftk.utils.WARN_REGION_SIZE = False
                    message("Encountered regions shorter than bin number.",
                            type="WARNING")
                    message(i.name +
                            " has length : " +
                            str(i.end -
                                i.start), type="WARNING")
                    message("They will be set to NA or --pseudo-count depending on --zero-to-na.",
                            type="WARNING")
                    message("Filter them out please.",
                            type="WARNING")

                if zero_to_na:
                    out = ['NA'] * bin_nb
                else:
                    out = [pc] * bin_nb

            else:

                try:
                    """
                    bw_cov = bigwig.stats(i.chrom,
                                          i.start,
                                          i.end,
                                          nBins=bin_nb)
                    """

                    bw_cov = bigwig.values(i.chrom,
                                           i.start,
                                           i.end)

                    out = []
                    size = i.end - i.start

                    for range_curr in intervals(list(range(size)), bin_nb, silent=True):

                        interval_cur = bw_cov[range_curr[0]:range_curr[1]]

                        if not zero_to_na:
                            interval_cur = [k if not np.isnan(k) else 0 for k in interval_cur]

                        if stat == 'mean':
                            out += [round(sum(interval_cur) / (range_curr[1] - range_curr[0]), 6)]
                        elif stat == 'sum':
                            out += [round(sum(interval_cur), 6)]
                        else:
                            raise GTFtkError("Stat should be 'sum' or 'mean'.")

                    if zero_to_na:
                        out = ['NA' if np.isnan(k) else k + pc for k in out]

                    else:
                        out = [pc if np.isnan(k) else k + pc for k in out]

                except:
                    if pygtftk.utils.WARN_UNDEF:
                        pygtftk.utils.WARN_UNDEF = False

                        mesg = "Encountered regions undefined in bigWig file."
                        message(mesg, type="WARNING")
                        mesg = '%s:%s-%s' % (i.chrom, str(i.start), str(i.end))
                        message(mesg)

                    if zero_to_na:
                        out = ['NA'] * bin_nb
                    else:
                        out = [pc] * bin_nb

            # Prepare output
            if i.name in ["", "."]:
                name = "|".join([i.chrom,
                                 str(i.start),
                                 str(i.end)])
            else:
                name = i.name

            if i.strand == "":
                strand = "."
            else:
                strand = i.strand

            # Print profiles
            if profile:

                # Data should be oriented in 5' -> 3'
                if stranded:

                    if i.strand == '-':
                        out = out[::-1]

                out = [str(x) for x in out]

                out_text = [label[cpt],
                            i.chrom,
                            str(i.start),
                            str(i.end),
                            str(i.name),
                            i.strand]
                out_text = out_text + out
                out_text = "\t".join(out_text)
                matrix_file.write(out_text + "\n")

            else:

                out = sorted(out, reverse=True)
                out = out[0:n_highest]

                if 'NA' not in out:
                    out = sum(out) / len(out)
                else:
                    out = 'NA'

                results.append("\t".join([i.chrom,
                                          str(i.start),
                                          str(i.end),
                                          label[cpt] + "|" + name,
                                          str(out),
                                          strand]) + "\n")

    if profile:
        matrix_file.close()
        return matrix_file.name

    else:
        return results


def bw_cov_mp(bw_list=None,
              region_file=None,
              labels=None,
              bin_nb=None,
              nb_proc=None,
              n_highest=None,
              zero_to_na=False,
              pseudo_count=None,
              stat='mean',
              verbose=False):
    """
    Compute bigwig coverage (multi-processed) for a set of regions.

    :param bw_list: the list of bigWig files to be processed.
    :param region_file: the bed file containing the region for which coverage is to be computed.
    :param labels: shortname for bigwigs.
    :param bin_nb: The number of bin into which the region should be splitted.
    :param nb_proc: Number of threads to be used.
    :param n_highest: compute the mean coverage based on the n highest values in the bins.
    :param pseudo_count: The value for a pseudo-count.
    :param verbose: run in verbose mode.
    :param stat: mean (default) or sum.
    :param zero_to_na: Convert missing values to NA, not zero.


    Returns a file.

    """

    n_region_to_proceed = len(BedTool(region_file.name))

    message("Received " +
            str(n_region_to_proceed) +
            " regions to proceed for each bigwig")

    tokens = intervals(list(range(n_region_to_proceed)), nb_proc)

    pool = multiprocessing.Pool(nb_proc)
    coverage_list = pool.map_async(_big_wig_coverage_worker,
                                   list(zip(tokens,
                                            repeat(bw_list),
                                            repeat(region_file.name),
                                            repeat(bin_nb),
                                            repeat(pseudo_count),
                                            repeat(n_highest),
                                            repeat(False),
                                            repeat(False),
                                            repeat(None),
                                            repeat(labels),
                                            repeat(zero_to_na),
                                            repeat(stat),
                                            repeat(verbose)))).get(9999999)

    if False in coverage_list:
        sys.stderr.write("Aborting...")
        sys.exit()

    # Unlist the list of list

    coverage_list = [item for sublist in coverage_list for item in sublist]

    tmp_file = make_tmp_file(prefix="region_coverage",
                             suffix=".bed")
    for i in coverage_list:
        tmp_file.write(i)

    tmp_file.close()

    return open(tmp_file.name)


def bw_profile_mp(in_bed_file=None,
                  nb_proc=None,
                  big_wig=None,
                  bin_nb=None,
                  pseudo_count=0,
                  stranded=True,
                  type=None,
                  labels=None,
                  outputfile=None,
                  zero_to_na=False,
                  bed_format=False,
                  add_score=False,
                  stat='mean',
                  verbose=False):
    """
    Compute bigwig profile for a set of regions.

    :param in_bed_file: the bed file containing the region for which coverage is to be computed.
    :param nb_proc: Number of threads to be used.
    :param big_wig: The bigWig files to be processed.
    :param bin_nb: The number of bin into which the region should be splitted.
    :param pseudo_count: The value for a pseudo-count.
    :param stranded: controls whether the profile should be ordered based on strand.
    :param type: This string will be added to the output to indicate the type of region (e.g tss, promoter...).
    :param labels: shortname for bigwigs.
    :param outputfile: output file name.
    :param zero_to_na: Convert missing values to NA, not zero.
    :param bed_format: Force Bed format. Default is to write columns in the following way: bwig, chrom, start, end, gene/feature, strand...
    :param add_score: add a 'score' column ("."). Just for downstream compatibility).
    :param stat: mean (default) or sum.
    :param verbose: run in verbose mode.

    Returns a file.

    """

    outputfile = add_prefix_to_file(infile=outputfile,
                                    prefix=type + "_")

    outputfile = open(outputfile, "w")

    n_region_to_proceed = len(BedTool(in_bed_file))

    message("Received " +
            str(n_region_to_proceed) +
            " regions to proceed for each bigwig")

    # 'Split' the file into multiple fragment
    tokens = intervals(list(range(n_region_to_proceed)), nb_proc)

    # Computing coverage of features.
    # Each worker will return a file
    pool = multiprocessing.Pool(processes=nb_proc)

    # Write a header
    if bed_format:
        prefix = []
    else:
        prefix = ["bwig"]

    suffix = [type + "_" + str(x + 1) for x in range(bin_nb)]

    if add_score:
        score_h = ["score"]
        score = ["."]
    else:
        score_h = []
        score = []

    outputfile.write("\t".join(prefix + ["chrom",
                                         "start",
                                         "end",
                                         "gene"] + score_h + ["strand"] + suffix) + "\n")

    if nb_proc > 1:
        argss = list(zip(tokens,
                         repeat(big_wig),
                         repeat(in_bed_file),
                         repeat(bin_nb),
                         repeat(pseudo_count),
                         repeat(None),
                         repeat(True),
                         repeat(stranded),
                         repeat(type),
                         repeat(labels),
                         repeat(zero_to_na),
                         repeat(stat),
                         repeat(verbose)))

        for res_file_list in pool.map_async(_big_wig_coverage_worker,
                                            argss).get(999999):

            for cur_file in flatten_list([res_file_list], outlist=[]):

                with open(cur_file) as infile:
                    for i in infile:
                        if bed_format:
                            i = i.split("\t")
                            outputfile.write("\t".join(i[1:4] + [i[0] + "|" + i[4]] + score + i[5:]))
                        else:
                            outputfile.write(i)
    # Don't use pool.
    else:
        for tok in tokens:
            res_file = _big_wig_coverage_worker((tok,
                                                 big_wig,
                                                 in_bed_file,
                                                 bin_nb,
                                                 pseudo_count,
                                                 None,
                                                 True,
                                                 stranded,
                                                 type,
                                                 labels,
                                                 zero_to_na,
                                                 stat,
                                                 verbose))

            with open(res_file) as infile:
                for i in infile:
                    if bed_format:
                        i = i.split("\t")
                        outputfile.write("\t".join(i[1:4] + [i[0] + "|" + i[4]] + ["."] + i[5:]))
                    else:
                        outputfile.write(i)

    close_properly(outputfile)

    return outputfile
