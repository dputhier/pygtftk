"""
This is the main function to compute random shuffles to assess statistical
significance overlap of two sets of genomic regions, provided as BED files.
"""

import functools as ft
import time
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np
import pybedtools

from pygtftk.stats.intersect import create_shuffles as cs
from pygtftk.stats.intersect import negbin_fit as nf
from pygtftk.stats.intersect.overlap import overlap_regions as oc
from pygtftk.stats.intersect.read_bed import read_bed_as_list as read_bed
from pygtftk.utils import message


################################################################################
# -------------------------------- MINIBATCH --------------------------------- #
################################################################################


def compute_all_intersections_minibatch(Lr1, Li1, Lr2, Li2,
                                        all_chrom1, all_chrom2,
                                        minibatch_size,
                                        use_markov_shuffling,
                                        nb_threads):
    """
    Main processing function. Computes a minibatch of shuffles for the given parameters.

    This function will be called by the hub function `compute_overlap_stats` to
    create a batch of shuffled "fake" BED files.
    Lr1,Li1,Lr2,Li2,all_chrom1,all_chrom2 are all outputs from the
    bed_Li1to_lists_of_intervals() function calls : those are the lists of region
    lengths, inter-region lengths, and chromosomes for each of the two input files.

    :param Lr1: An output from the bed_to_lists_of_intervals() function calls.
    :param Li1: An output from the bed_to_lists_of_intervals() function calls.
    :param Lr2: An output from the bed_to_lists_of_intervals() function calls.
    :param Li2: An output from the bed_to_lists_of_intervals() function calls.
    :param all_chrom1:  An output from the bed_to_lists_of_intervals() function calls.
    :param all_chrom2: An output from the bed_to_lists_of_intervals() function calls.
    :param minibatch_size: The size of the batchs for shuffling.
    :param use_markov_shuffling: Use a classical or a order-2 Markov shuffling.
    :param nb_threads: number of threads.


    """

    # --------------------- Generate and shuffle batches  -------------------- #
    # We generate a matrix with the batches and shuffle them independantly
    # for both bed files

    # Get chromosomes that are common between both BEDs
    all_chroms = np.intersect1d(all_chrom1, all_chrom2)
    all_chroms = [str(x) for x in all_chroms]  # Revert from numpy list to traditional python list

    shuffled_Lr1_batches, shuffled_Li1_batches = dict(), dict()
    shuffled_Lr2_batches, shuffled_Li2_batches = dict(), dict()

    ## Wrapper to make the code cleaner
    # Tile the list of length (repeat) as many times as we want shuffles, then
    # shuffle the rows independantly.

    # We can use a classical or a order-2 Markov shuffling. This results in different wrappers here :
    if not use_markov_shuffling:
        def batch_and_shuffle_list(l): return cs.shuffle(np.tile(l, (minibatch_size, 1)))
    if use_markov_shuffling:
        def batch_and_shuffle_list(l): return cs.markov_shuffle(np.tile(l, (minibatch_size, 1)), nb_threads=nb_threads)

    # NOTE for improvement : if new types of shuffles are added, the corresponding
    # wrappers should be added here.

    # Produce the shuffles on a chromosome basis
    start = time.time()
    for chrom in all_chroms:
        shuffled_Lr1_batches[chrom] = batch_and_shuffle_list(Lr1[chrom])
        shuffled_Li1_batches[chrom] = batch_and_shuffle_list(Li1[chrom])
        shuffled_Lr2_batches[chrom] = batch_and_shuffle_list(Lr2[chrom])
        shuffled_Li2_batches[chrom] = batch_and_shuffle_list(Li2[chrom])
    stop = time.time()
    message('Batch generated and shuffled in ' + str(stop - start) + ' s.', type='DEBUG')

    # -------------------- Convert batches into BED files -------------------- #
    start = time.time()
    batch_to_bedlist_with_params = ft.partial(cs.batch_to_bedlist, all_chroms=all_chroms, minibatch_size=minibatch_size,
                                              nb_threads=nb_threads)
    bedsA, bedsB = batch_to_bedlist_with_params(shuffled_Lr1_batches, shuffled_Li1_batches, shuffled_Lr2_batches,
                                                shuffled_Li2_batches)
    stop = time.time()
    message('Batch converted to fake beds in : ' + str(stop - start) + ' s.', type='DEBUG')

    # -------------------- Processing intersections -------------------------- #
    # Using our custom cython intersect, process intersection between each pair
    # of 'fake bed files'
    start = time.time()
    all_intersections = oc.compute_intersections_cython(bedsA, bedsB, all_chroms, nb_threads)
    stop = time.time()
    message('All intersections computed by custom code in : ' + str(stop - start) + ' s.', type='DEBUG')

    return all_intersections


################################################################################
# ---------------------------------- CORE ------------------------------------ #
################################################################################

def compute_overlap_stats(bedA, bedB,
                          chrom_len,
                          minibatch_size, minibatch_nb,
                          bed_excl,
                          use_markov_shuffling,
                          nb_threads,
                          ft_type):
    """
    This is the hub function to compute overlap statistics through Monte Carlo
    shuffling with integration of the inter-region lengths.

    The function will generate shuffled BEDs from bedA and bedB independantly,
    and compute intersections between those shuffles. As such, it gives an
    estimation of the intersections under the null hypothesis (the sets of
    regions given in A and B are independant).

    Author : Quentin FERRE <quentin.q.ferre@gmail.com>

    :param bedA: The first bed file.
    :param bedB: The second bed file.
    :param chrom_len: the dictionary of chromosome lengths
    :param minibatch_size: the size of the minibatch for shuffling.
    :param minibatch_nb: The number of minibatchs.
    :param bed_excl: The regions to be excluded.
    :param use_markov_shuffling: Use a classical or a order-2 Markov shuffling.
    :param nb_threads: Number of threads.
    :param ft_type: The name of the feature.

    """

    message('Beginning shuffling for ' + ft_type)
    message('BedA: ' + bedA.fn, type='DEBUG')
    message('Nb. features in BedA: ' + str(len(bedA)), type='DEBUG')
    message('BedB: ' + bedB.fn, type='DEBUG')
    message('Nb. features in BedB: ' + str(len(bedB)), type='DEBUG')
    message('BATCHES : ' + str(minibatch_nb) + ' batches of ' + str(minibatch_size) + ' shuffles.', type='DEBUG')
    message('Total number of shuffles : ' + str(minibatch_nb * minibatch_size) + '.', type='DEBUG')
    message('NB_THREADS = ' + str(nb_threads) + '.', type='DEBUG')

    # --------------------- Read list of intervals --------------------------- #
    start = time.time()

    # Just in case, force type and merge bedA ; same for bedB
    bed_A_as_pybedtool = pybedtools.BedTool(bedA).sort().merge()
    bed_B_as_pybedtool = pybedtools.BedTool(bedB).sort().merge()

    # If there is an exclusion to be done, do it.

    # NOTE : exclusion on the peak file (bedA) has been moved to ologram itself to avoid repetition. Same thing for the chromsizes.
    # This means that when there is an exclusion to be done, this function must do in on bedB only.
    if bed_excl is not None:
        exclstart = time.time()
        message('Performing exclusion on the second file, proceeding. This may take a few minutes.', type='INFO')

        bed_B_as_pybedtool = read_bed.exclude_concatenate(bed_B_as_pybedtool, bed_excl, nb_threads)

        exclstop = time.time()
        message('Exclusion completed in ' + str(exclstop - exclstart) + ' s.', type='DEBUG')

    # Abort if there are less than 2 remaining regions in bedA and bedB
    if (len(bed_A_as_pybedtool) < 2) | (len(bed_B_as_pybedtool) < 2):
        message(
            'Less than 2 remaining regions in one of the BED files. This is likely due to either : one of the considered features has very few peaks falling inside of it, or all the regions are in areas marked in the exclusion file. ologram will discard this particular pair.',
            type='WARNING')

        # Return a result dict full of -1
        result_abort = OrderedDict()
        result_abort['nb_intersections_expectation_shuffled'] = 0;
        result_abort['nb_intersections_variance_shuffled'] = 0
        result_abort['nb_intersections_negbinom_fit_quality'] = -1;
        result_abort['nb_intersections_log2_fold_change'] = 0
        result_abort['nb_intersections_true'] = 0;
        result_abort['nb_intersections_pvalue'] = -1;
        result_abort['summed_bp_overlaps_expectation_shuffled'] = 0;
        result_abort['summed_bp_overlaps_variance_shuffled'] = 0
        result_abort['summed_bp_overlaps_negbinom_fit_quality'] = -1;
        result_abort['summed_bp_overlaps_log2_fold_change'] = 0
        result_abort['summed_bp_overlaps_true'] = 0;
        result_abort['summed_bp_overlaps_pvalue'] = -1
        return result_abort

    # Proper reading of the bed file as a list of intervals
    Lr1, Li1, all_chrom1 = read_bed.bed_to_lists_of_intervals(bed_A_as_pybedtool, chrom_len)
    Lr2, Li2, all_chrom2 = read_bed.bed_to_lists_of_intervals(bed_B_as_pybedtool, chrom_len)
    stop = time.time()
    message('BED files read as lists of intervals in ' + str(stop - start) + ' s', type='DEBUG')

    grand_start = time.time()

    ################################ MINIBATCH  ################################
    # Generate all intersections for a shuffled batch of size n

    minibatches = [minibatch_size for i in range(minibatch_nb)]
    all_intersections = list()
    for k in range(len(minibatches)):
        # Display of current progress
        message("--- Minibatch nb. : " + str(k + 1) + " / " + str(minibatch_nb), type='DEBUG')

        all_intersections = all_intersections + compute_all_intersections_minibatch(Lr1, Li1, Lr2, Li2, all_chrom1,
                                                                                    all_chrom2,
                                                                                    minibatches[k],
                                                                                    use_markov_shuffling, nb_threads)
    message('All intersections have been generated.', type='DEBUG')

    # --------------- Compute statistics on the intersections ---------------- #

    # NOTE For future improvement, since the shuffling itself is done chromosome
    # by chromosome, making some statistics 'by chromosome' should be possible.

    start = time.time()
    with Pool(nb_threads) as p:
        stats = p.map(read_bed.compute_stats_for_intersection, all_intersections)
    stop = time.time()
    message('Statistics on overlaps computed in : ' + str(stop - start) + ' s.', type='DEBUG')

    # Unpack the stats
    bp_overlaps = [s[0] for s in stats]
    # unlisted_bp_overlaps = [item for sublist in bp_overlaps for item in sublist] # TODO : only useful if we plot them individually
    summed_bp_overlaps = [sum(x) for x in bp_overlaps]
    intersect_nbs = [s[1] for s in stats]

    # NOTE FOR IMPROVEMENT : it would be interesting to return the average size
    # of an overlap as well, per shuffle. Since our intersection algorithm returns
    # details about the intersections like `bedtools intersect` would, this could
    # be computed without much hassle.

    # ------ Fitting of a Negative Binomial distribution on the shuffles ----- #
    # Only relevant for classical shuffle, not Markov
    start = time.time()

    # We saw experimentally that Markov shuffles do not fit the Neg Binom model.
    if use_markov_shuffling:
        ps = pn = -1

    else:
        # Renaming expectations and variances
        expectation_fitted_summed_bp_overlaps, variance_fitted_summed_bp_overlaps = np.mean(summed_bp_overlaps), np.var(
            summed_bp_overlaps)
        expectation_fitted_intersect_nbs, variance_fitted_intersect_nbs = np.mean(intersect_nbs), np.var(intersect_nbs)

        ## Check that there is a good adjustment.
        # This is done using 1 minus Cramer's V score ; a good adjustment should return a value close to 1
        # See doc of check_negbin_adjustment() for details
        # NOTE Checking adjustment is meaningless if the expectation is zero
        if expectation_fitted_summed_bp_overlaps == 0:
            ps = -1
        else:
            ps = nf.check_negbin_adjustment(summed_bp_overlaps, expectation_fitted_summed_bp_overlaps,
                                            variance_fitted_summed_bp_overlaps)

        if expectation_fitted_intersect_nbs == 0:
            pn = -1
        else:
            pn = nf.check_negbin_adjustment(intersect_nbs, expectation_fitted_intersect_nbs,
                                            variance_fitted_intersect_nbs)

    # Send warnings when there is a poor fit
    if (ps < 0.75) | (pn < 0.75):
        message(ft_type + ': there may be a poor fit for this feature. '
                          'Check fit quality in the results. This is likely due '
                          'to there being too few regions.',
                type='WARNING')

    # -------------------------- True intersections -------------------------- #
    # Now, calculating the actual p-value for the number of intersections and the
    # total number of overlapping base pairs

    ## True intersection
    # true_intersection = bedA.intersect(bedB)
    true_intersection = bed_A_as_pybedtool.intersect(
        bed_B_as_pybedtool)  # Perform intersection with the exclusion regions removed !

    true_intersect_nb = len(true_intersection)
    true_bp_overlaps = sum([x.length for x in true_intersection])

    # Compute the p-values using the distribution fitted on the shuffles
    # Do not do this for the Markov shuffling, as it is likely a multi-variable fit (see notes)

    # We can only use a Neg Binom p-val if we can fit it, and that is not the case for
    # the Markov shuffle or if the expectation is too small : we must use an empirical p-value

    if (ps == -1) | (pn == -1):
        # NOTE : maybe re-use the empirical p-value later. For now return -1
        # pval_intersect_nb = nf.empirical_p_val(true_intersect_nb, intersect_nbs)
        # pval_bp_overlaps = nf.empirical_p_val(true_bp_overlaps, summed_bp_overlaps)
        pval_intersect_nb = -1
        pval_bp_overlaps = -1

    else:
        pval_intersect_nb = nf.negbin_pval(true_intersect_nb,
                                           expectation_fitted_intersect_nbs,
                                           variance_fitted_intersect_nbs,
                                           ft_type=ft_type)
        pval_bp_overlaps = nf.negbin_pval(true_bp_overlaps,
                                          expectation_fitted_summed_bp_overlaps,
                                          variance_fitted_summed_bp_overlaps,
                                          ft_type=ft_type)

    stop = time.time()
    message('Negative Binomial distributions fitted in : ' + str(stop - start) + ' s (' + ft_type + ').', type='DEBUG')

    grand_stop = time.time()

    message('--- Total time : ' + str(grand_stop - grand_start) + ' s for feature : ' + ft_type + ' ---')
    message('Total time does not include BED reading, as it does not scale with batch size.', type='DEBUG')

    # ------------------------------------------------------------------------
    # Draft code for diagnostic plots of the distribution of each statistic in
    # the shuffles. Kept for potential future improvement.
    # ------------------------------------------------------------------------
    # import matplotlib.pyplot as plt
    #
    # ## Number of overlapping base pairs
    # # Sum by batch
    # plt.figure() ; plt.hist(summed_bp_overlaps, bins=50)
    # plt.savefig(outputdir+'/'+name+'_nb_overlapping_bp_sum_by_batch.png')
    # # All individual lines across all batches
    # plt.figure() ; plt.hist(np.array(unlisted_bp_overlaps), bins=300)
    # plt.savefig(outputdir+'/'+name+'_nb_overlapping_bp_individual.png')
    #
    # # Length of overlapping REFERENCE regions (from BED_B) with which the query (from BED_A) intersected
    # if lengths_wb is not None:
    #     plt.figure() ; plt.hist(np.array(lengths_wb), bins=300)
    #     plt.savefig(outputdir+'/'+name+'_length_reference_regions_intersect.png')
    # # For comparison, length of all REFERENCE regions (from BED_B)
    # lengths = [r.length for r in pybedtools.BedTool(BED_B)]
    # plt.figure() ; plt.hist(lengths, bins=300)
    # plt.savefig(outputdir+'/'+name+'_length_reference_regions_all.png')
    #
    # # Number of intersections
    # plt.figure() ; plt.hist(intersect_nbs, bins=50)
    # plt.savefig(outputdir+'/'+name+'_nb_intersections.png')

    ### Result as a dictionary of statistics
    # WARNING Be careful to use the same order as result_abort, above !
    result = OrderedDict()

    # Number of intersections
    result['nb_intersections_expectation_shuffled'] = '{:.2f}'.format(np.mean(intersect_nbs))
    result['nb_intersections_variance_shuffled'] = '{:.2f}'.format(np.var(intersect_nbs))

    result['nb_intersections_negbinom_fit_quality'] = '{:.5f}'.format(pn)

    if np.mean(intersect_nbs) == 0:
        ni_fc = 0  # Do not divide by zero !
    else:
        ni_fc = true_intersect_nb / np.mean(intersect_nbs)
    if ni_fc != 0: ni_fc = np.log2(ni_fc)  # Apply log transformation
    result['nb_intersections_log2_fold_change'] = '{:.5f}'.format(ni_fc)

    result['nb_intersections_true'] = true_intersect_nb
    result['nb_intersections_pvalue'] = '{0:.4g}'.format(pval_intersect_nb)

    # Summed number of overlapping basepairs
    result['summed_bp_overlaps_expectation_shuffled'] = '{:.2f}'.format(np.mean(summed_bp_overlaps))
    result['summed_bp_overlaps_variance_shuffled'] = '{:.2f}'.format(np.var(summed_bp_overlaps))

    result['summed_bp_overlaps_negbinom_fit_quality'] = '{:.5f}'.format(ps)

    if np.mean(summed_bp_overlaps) == 0:
        sbp_fc = 0  # Do not divide by zero !
    else:
        sbp_fc = true_bp_overlaps / np.mean(summed_bp_overlaps)
    if sbp_fc != 0: sbp_fc = np.log2(sbp_fc)  # Apply log transformation
    result['summed_bp_overlaps_log2_fold_change'] = '{:.5f}'.format(sbp_fc)

    result['summed_bp_overlaps_true'] = true_bp_overlaps
    result['summed_bp_overlaps_pvalue'] = '{0:.4g}'.format(pval_bp_overlaps)

    return result
