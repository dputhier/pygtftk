"""
This is the main function to compute random shuffles to assess statistical
significance overlap of two sets of genomic regions, provided as BED files.

Called directly by the ologram.py plugin. All other functions calls descend from
this one.

Author : Quentin Ferré <quentin.q.ferre@gmail.com>
"""

import functools as ft
import gc
import multiprocessing
import time
from collections import OrderedDict
import concurrent.futures as cf
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pybedtools

from pygtftk.stats.intersect import create_shuffles as cs
from pygtftk.stats.intersect import overlap_stats_compute as osc
from pygtftk.stats.intersect.overlap import overlap_regions as oc
from pygtftk.stats.intersect.read_bed import read_bed_as_list as read_bed
from pygtftk.utils import message


################################################################################
# -------------------------------- MINIBATCH --------------------------------- #
################################################################################


def compute_all_intersections_minibatch(Lr1, Li1, Lrs, Lis,
                                        all_chrom1, all_chrom2,
                                        minibatch_size,
                                        use_markov_shuffling,
                                        keep_intact_in_shuffling,
                                        nb_threads, seed=42):
    """
    Main processing function. Computes a minibatch of shuffles for the given parameters.

    This function will be called by the hub function `compute_overlap_stats` to
    create a batch of shuffled "fake" BED files.
    The Lrs and Li and all_chroms are all outputs from the
    bed_to_lists_of_intervals() function calls : those are the lists of region
    lengths, inter-region lengths, and chromosomes for each of the two input files.

    Lrs and Lis are lists containing [Lr2] and [Li2] if there is only one query
    set, but can contain [Lr2, Lr3, ...] if there is more than one.

    :param Lr1: An output from the bed_to_lists_of_intervals() function calls.
    :param Li1: An output from the bed_to_lists_of_intervals() function calls.
    :param Lrs: A list of outputs from the bed_to_lists_of_intervals() function calls.
    :param Lis: A list of outputs from the bed_to_lists_of_intervals() function calls.
    :param all_chrom1:  An output from the bed_to_lists_of_intervals() function calls.
    :param all_chrom2: An output from the bed_to_lists_of_intervals() function calls.
    :param minibatch_size: The size of the batchs for shuffling.
    :param use_markov_shuffling: Use a classical or a order-2 Markov shuffling.
    :param keep_intact_in_shuffling: Among the Lrs (and Lis), those whose number/order is here will be kept intact during the shuffling (won't be shuffled)
    :param nb_threads: number of threads.

    """

    # Re-seed RNG
    np.random.seed(seed)

    # --------------------- Generate and shuffle batches  -------------------- #
    # We generate a matrix with the batches and shuffle them independantly
    # for both bed files

    # Get chromosomes that are common between both BEDs
    all_chroms = np.intersect1d(all_chrom1, all_chrom2)
    all_chroms = [str(x) for x in all_chroms]  # Revert from numpy list to traditional python list

    shuffled_Lr1_batches, shuffled_Li1_batches = dict(), dict()

    shuffled_Lrs_batches = list()
    shuffled_Lis_batches = list()
    for set in range(len(Lis)):
        shuffled_Lrs_batches += [dict()]
        shuffled_Lis_batches += [dict()]

    ## Wrapper to make the code cleaner
    # Tile the list of length (repeat) as many times as we want shuffles, then
    # shuffle the rows independantly.

    # We can use a classical or a order-2 Markov shuffling. This results in different wrappers here :
    if not use_markov_shuffling:
        def batch_and_shuffle_list(l): return cs.shuffle(np.tile(l, (minibatch_size, 1)))
    if use_markov_shuffling:
        def batch_and_shuffle_list(l): return cs.markov_shuffle(np.tile(l, (minibatch_size, 1)), nb_threads=nb_threads)

    # --------------------------------------------------------------------------
    # NOTE for improvement : 
    #   - If new types of shuffles are added, the corresponding wrappers should 
    # be added here. All that matters is that these shufflers return a shuffled 
    # matrix of distances. 
    #   - We could also implement weighted shuffles, or shuffles that rely
    # on an external seed (for example, ensuring that long regions fall in the
    # middle of the chromosomes for whatever reason). All that matters is
    # to plug them here as lambdas.
    #   - Relatedly, to have a common seed or remember information between the
    # shuffles (ie. shuffle sets A and B the same way), we just need to make
    # batch_and_shuffle_list take two arguments : (l, parameter) and change the
    # function call below.
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # NOTE for improvement : We may wish to shuffle across all chromosomes,
    # instead of chromosome-by-chromosome.
    # A workaround to do so is implemented here by exchanging the base Lrs and
    # Lis (the true ones) before producing the N shuffles
    # --------------------------------------------------------------------------
    # TODO: This is not used yet (the code is commented and no arguments call it)
    # because we need to agree on implementation parameters : namely, when
    # shuffling across all chromosomes, should we put the same number of regions
    # on each chromosome ? Follow a Poisson law ? etc.
    # --------------------------------------------------------------------------
    """
    def shuffle_L_across_chrom(L):
        '''
        Takes a dictionary of lists and returns a new dictionary where the elements
        of all the lists have been shuffled across all lists.

        This simply consists in flattening and rebuilding the dictionary.
        '''

        size_per_chrom = dict()
        flattened_list = list()

        # Flatten the dictionary
        for chrom in L.keys():
            L_for_this_chrom = L[chrom] # This is a list !

            # TODO: This is what needs to be fixed. For now it rebuilds the
            # dictionary with the same number or arguments as before
            size_per_chrom[chrom] = len(L_for_this_chrom)

            flattened_list += L_for_this_chrom

        # Rebuild the dictionary
        for chrom in L.keys():
            # Pick as many elements from the flattened list as there were before
            # and remove them from the flattened list
            L[chrom] = [flattened_list.pop(np.random.randint(len(flattened_list))) for _ in range(size_per_chrom[chrom])]

        return L

    # Shuffle the *true* Lrs and Lis, for all.
    if shuffle_across_all_chrom:

        Lr1 = shuffle_L_across_chrom(Lr1)
        Li1 = shuffle_L_across_chrom(Li1)

        for set in range(len(Lis)):
            Lrs[set] = shuffle_L_across_chrom(Lrs[set])
            Lis[set] = shuffle_L_across_chrom(Lis[set])
    """

    """
    TODO: We could also use this to output stats per chromosome : In results,
    instead of having "exons", "intergenic", etc., have "exons_chr1", 
    "exons_chr2", "intergenic_chr1", etc.
    """

    # Translating the comma-separated string of "keep_intact_in_shuffling" into a list of ids
    if keep_intact_in_shuffling is None: keep_intact_in_shuffling = []
    else: keep_intact_in_shuffling = keep_intact_in_shuffling.split(',')

    # Produce the shuffles on a chromosome basis
    start = time.time()
    for chrom in all_chroms:
        shuffled_Lr1_batches[chrom] = batch_and_shuffle_list(Lr1[chrom])
        shuffled_Li1_batches[chrom] = batch_and_shuffle_list(Li1[chrom])

        for set_id in range(len(Lis)):

            # If this set id is supposed to be kept intact in the shuffling, do not shuffle it and keep its regions at the same locations
            if set_id in keep_intact_in_shuffling:
                shuffled_Lrs_batches[set_id][chrom] = Lrs[set_id][chrom]

            else:
                # Some BEDs may have no peaks on certain chromosomes.
                # Watch for KeyError exception for this case.
                try:
                    shuffled_Lrs_batches[set_id][chrom] = batch_and_shuffle_list(Lrs[set_id][chrom])
                except KeyError:
                    shuffled_Lrs_batches[set_id][chrom] = np.tile([0], (minibatch_size, 1))

                try:
                    shuffled_Lis_batches[set_id][chrom] = batch_and_shuffle_list(Lis[set_id][chrom])
                except KeyError:
                    shuffled_Lis_batches[set_id][chrom] = np.tile([0, 0], (minibatch_size, 1))

    stop = time.time()
    message('Batch generated and shuffled in ' + str(stop - start) + ' s.', type='DEBUG')

    # -------------------- Convert batches into BED files -------------------- #
    start = time.time()
    batch_to_bedlist_with_params = ft.partial(cs.batch_to_bedlist, all_chroms=all_chroms, minibatch_size=minibatch_size,
                                              nb_threads=nb_threads)

    bedsA = batch_to_bedlist_with_params(shuffled_Lr1_batches, shuffled_Li1_batches)

    bedsB = list()

    for k in range(len(Lis)):
        bedsB += [batch_to_bedlist_with_params(shuffled_Lrs_batches[k], shuffled_Lis_batches[k])]

    stop = time.time()
    message('Batch converted to fake beds in : ' + str(stop - start) + ' s.', type='DEBUG')

    # -------------------- Processing intersections -------------------------- #
    # Using our custom cython intersect, process intersection between each pair
    # of 'fake bed files'
    start = time.time()
    # WARNING : bedsA is a list of beds, but bedB is a list of list of beds !
    all_intersections = oc.compute_intersections_cython(bedsA, bedsB, all_chroms, nb_threads)
    stop = time.time()
    message('All intersections computed by custom code in : ' + str(stop - start) + ' s.', type='DEBUG')

    return all_intersections


################################################################################
# ---------------------------------- CORE ------------------------------------ #
################################################################################

class ComputingIntersectionPartial(object):
    """
    This is a replacer for a wrapper for compute_all_intersections_minibatch, needed for multiprocessing below
    This was needed because the argument that changes in multiproessing in minibatch_len, and it's not the leftmost argument, 
    so functools.partial cannot be used
    """

    # Remember the parameters
    def __init__(self, Lr1, Li1, Lrs, Lis, all_chrom1, all_chrom2, use_markov_shuffling, keep_intact_in_shuffling, nb_threads,):
        # Parameters for compute_all_intersections_minibatch
        self.Lr1 = Lr1
        self.Li1 = Li1
        self.Lrs = Lrs
        self.Lis = Lis
        self.all_chrom1 = all_chrom1
        self.all_chrom2 = all_chrom2
        self.use_markov_shuffling = use_markov_shuffling
        self.keep_intact_in_shuffling = keep_intact_in_shuffling
        self.nb_threads = nb_threads

        
    # Callable
    def __call__(self, minibatch_len, seed, id):
        my_result = compute_all_intersections_minibatch(self.Lr1, self.Li1, self.Lrs, self.Lis, self.all_chrom1,
                                                        self.all_chrom2, minibatch_len, self.use_markov_shuffling,
                                                        self.keep_intact_in_shuffling, self.nb_threads, seed=seed)

        message("--- Minibatch nb. : " + str(id) + " is complete.")

        return my_result


def compute_overlap_stats(bedA, bedsB,
                          chrom_len,
                          minibatch_size, minibatch_nb,
                          bed_excl,
                          use_markov_shuffling,
                          keep_intact_in_shuffling,
                          nb_threads,
                          ft_type,
                          multiple_overlap_target_combi_size=None,
                          multiple_overlap_max_number_of_combinations=None,
                          multiple_overlap_custom_combis=None,
                          draw_histogram=False,
                          exact = False):
    """
    This is the hub function to compute overlap statistics through Monte Carlo
    shuffling with integration of the inter-region lengths.

    The function will generate shuffled BEDs from bedA and bedB independantly,
    and compute intersections between those shuffles. As such, it gives an
    estimation of the intersections under the null hypothesis (the sets of
    regions given in A and B are independant).

    Author : Quentin FERRE <quentin.q.ferre@gmail.com>

    :param bedA: The first bed file in pybedtools.Bedtool format.
    :param bedsB: The second bed file. Can also be a list of bed files for multiple overlaps.
    :param chrom_len: the dictionary of chromosome lengths
    :param minibatch_size: the size of the minibatch for shuffling.
    :param minibatch_nb: The number of minibatchs.
    :param bed_excl: The regions to be excluded.
    :param use_markov_shuffling: Use a classical or a order-2 Markov shuffling.
    :param keep_intact_in_shuffling: those numbers in bedsB will be kept intact/fixed during the shuffliing
    :param nb_threads: Number of threads.
    :param ft_type: The name of the feature.
    :param multiple_overlap_target_combi_size: For multiple overlaps, maximum number of sets in the output combinations.
    :param multiple_overlap_max_number_of_combinations: For multiple overlaps, maximum number of combinations to consider. This will use the MOLD mining algorithm. Do not ask for too many.
    :param multiple_overlap_custom_combis: For multiple overlaps, skips combination mining and computes stats directly for these combinations. Path to a file to be read as NumPy matrix.
    :param draw_histogram: if True, draws a temp file histogram for each combi
    :param exact: if True, when performing the statistics, an observation of A+B+C counts as an observation of A+B

    """

    # ------------------------------------------------------------------------ #
    # CAPITAL - If bedB is a singleton, make it a list in bedsB (as in, plural)
    multiple_overlaps_were_originally_requested = isinstance(bedsB, list)  # Rememeber if we queried multiple overlaps
    bedsB = [bedsB] if not isinstance(bedsB, list) else list(bedsB)
    # ------------------------------------------------------------------------ #

    message('Beginning the computation of overlap stats for ' + str(ft_type))
    message('BedA: ' + bedA.fn, type='DEBUG')
    message('Nb. features in BedA: ' + str(len(bedA)), type='DEBUG')
    for bedB in bedsB:
        message('BedB: ' + bedB.fn, type='DEBUG')
        message('Nb. features in BedB: ' + str(len(bedB)), type='DEBUG')
    message('BATCHES : ' + str(minibatch_nb) + ' batches of ' + str(minibatch_size) + ' shuffles.', type='DEBUG')
    message('Total number of shuffles : ' + str(minibatch_nb * minibatch_size) + '.', type='DEBUG')
    message('NB_THREADS = ' + str(nb_threads) + '.', type='DEBUG')

    # --------------------- Read list of intervals --------------------------- #
    start = time.time()

    # Just in case, force type and merge bedA ; same for bedB
    bedA = pybedtools.BedTool(bedA).sort().merge()
    bedsB = [pybedtools.BedTool(bedB).sort().merge() for bedB in bedsB]

    stop = time.time()
    message('BED files merged and sorted via PyBedtools in ' + str(stop - start) + ' s', type='DEBUG')

    # If there is an exclusion to be done, do it.

    # NOTE : exclusion on the peak file (bedA) has been moved to ologram itself to avoid repetition. Same thing for the chromsizes.
    # This means that when there is an exclusion to be done, this function must do in on bedB only.
    if bed_excl is not None:
        exclstart = time.time()
        message('Performing exclusion on the second files. This may take a moment.', type='INFO')

        bedsB = [read_bed.exclude_concatenate(bedB, bed_excl, nb_threads) for bedB in bedsB]

        exclstop = time.time()
        message('Exclusion completed in ' + str(exclstop - exclstart) + ' s.', type='DEBUG')

    was_more_than_one_bedB = (len(bedsB) > 1)  # Remember if there were more than 1 BED in bedsB

    # Remove any bed in bedsB with less than 2 region
    bedsB_final = []
    for bedB in bedsB:
        if (len(bedB) >= 2):
            bedsB_final += [bedB]
        else:
            bedsB_final += [
                pybedtools.BedTool("chr1\t1\t1\n" + str(bedB), from_string=True)
            ]
            message(
                'Less than 2 remaining regions in one of the second BED files. This is likely due to either : one of the considered features has very few peaks falling inside of it, or all the regions are in areas marked in the exclusion file. ologram will adding a ghost, null-length region to proceed.',
                type='WARNING')
    bedsB = bedsB_final

    # If more_bed_multiple_overlap is True, ensure there are at least 2 sets
    if multiple_overlaps_were_originally_requested and (len(bedsB) < 2):
        raise ValueError(
            "--more-bed-multiple-overlap was used, but (after potential exclusion) less than 2 sets of regions are to be compared against the query.")

    # Abort if there are less than 2 remaining regions in bedA or if bedsB is empty
    if (len(bedA) < 2) | (len(bedsB) == 0):
        message(
            'Less than 2 remaining regions in the query bed file, or no file remaining in the set of regions to be compared against. Aborting.',
            type='WARNING')

        # Return a result dict full of -1
        result_abort = OrderedDict()
        result_abort['nb_intersections_expectation_shuffled'] = 0
        result_abort['nb_intersections_variance_shuffled'] = 0
        result_abort['nb_intersections_negbinom_fit_quality'] = -1
        result_abort['nb_intersections_log2_fold_change'] = 0
        result_abort['nb_intersections_true'] = 0
        result_abort['nb_intersections_pvalue'] = -1
        result_abort['summed_bp_overlaps_expectation_shuffled'] = 0
        result_abort['summed_bp_overlaps_variance_shuffled'] = 0
        result_abort['summed_bp_overlaps_negbinom_fit_quality'] = -1
        result_abort['summed_bp_overlaps_log2_fold_change'] = 0
        result_abort['summed_bp_overlaps_true'] = 0
        result_abort['summed_bp_overlaps_pvalue'] = -1
        result_abort['combination_order'] = 0
        result_abort['nb_intersections_empirical_pvalue'] = -1
        result_abort['summed_bp_overlaps_empirical_pvalue'] = -1
        result_abort['beta_summed_bp_overlaps_pvalue_ad_hoc_for_deep_sampling_only'] = -1

        # If it was a multiple overlap : return a nested dict, otherwise return a classical dict
        if was_more_than_one_bedB:
            return {"multiple_beds": result_abort}
        else:
            return result_abort

    start = time.time()

    ## Proper reading of the bed files as a list of intervals, for bedA and also
    # all files in bedsB
    Lr1, Li1, all_chrom1 = read_bed.bed_to_lists_of_intervals(bedA, chrom_len)
    all_chrom1 = all_chrom1.astype(str)

    Lrs = list()
    Lis = list()

    all_chrom2 = list()
    for k in range(len(bedsB)):
        Lrs_toappend, Lis_toappend, all_chrom2_toappend = read_bed.bed_to_lists_of_intervals(bedsB[k], chrom_len)
        Lrs += [Lrs_toappend]
        Lis += [Lis_toappend]

        # WARNING note that if there are multiple files in bedsB, all_chrom2 is overwritten every time ! NO, NOT ANYMORE !
        all_chrom2 = list(all_chrom2) + list(all_chrom2_toappend)
        all_chrom2 = np.unique(all_chrom2)

    # Type safety : force cast to string to prevent cases where chromosomes names are sometimes read as integers
    all_chrom1 = all_chrom1.astype(str)
    all_chrom2 = all_chrom2.astype(str)

    stop = time.time()
    message('BED files read as lists of intervals in ' + str(stop - start) + ' s', type='DEBUG')

    grand_start = time.time()

    ################################ MINIBATCH  ################################
    # Generate all intersections for a shuffled batch of size n

    minibatches = [minibatch_size for i in range(minibatch_nb)]

    ## Compute intersections for each minibatch, multiprocessed

    # Result queue
    all_intersections = list()  # Final result list, to be filled at the end
    pool = ProcessPoolExecutor(nb_threads)  # Process pool

    # Prepare one random seed per batch. This is done to prevent a problem in
    # multiprocessing : since this function will be mutithreaded, we must ensure
    # each thread has been given a different random seed.
    seeds = [np.random.randint(2 ** 32) for _ in range(len(minibatches))]

    # Create a sort-of partial call
    compute_intersection_partial = ComputingIntersectionPartial(Lr1, Li1, Lrs, Lis, all_chrom1, all_chrom2,
                                                                use_markov_shuffling, keep_intact_in_shuffling, nb_threads)

    # Submit to the pool of processes
    futures = list()
    message("We will perform a total of " + str(len(minibatches)) + " batches of " + str(minibatch_size) + " shufflings.")
    for i in range(len(minibatches)):
        futures += [pool.submit(compute_intersection_partial, minibatch_len=minibatches[i], seed=seeds[i], id = i)]
    pool.shutdown() # Release the resources as soon as you are done with those, we won't submit any more jobs to you



    # Transpose the results into all_intersections
    for future in futures:
        try:
            all_intersections += future.result()
        except cf.process.BrokenProcessPool:
            message('A process in the process pool was terminated abruptly while the future was running or pending. This likely means you ran out of RAM. Try restarting the command with fewer threads or smaller minibatches', type="ERROR")



    ## Cleanup
    del pool
    time.sleep(0.5); gc.collect()

    message("Total number of shuffles, reminder : " + str(len(all_intersections)), type='DEBUG')
    message("Number of intersections in the first shuffle, for comparison : " + str(len(all_intersections[0])),
            type='DEBUG')

    message('All intersections have been generated.', type='INFO')

    # The `all_intersections` objects contains all the computed overlaps, 
    # one per shuffle. All shuffles are concatenated.


    # --------------- Compute statistics on the intersections ---------------- #

    # NOTE For future improvement, since the shuffling itself is done chromosome
    # by chromosome, making some statistics 'by chromosome' should be possible.

    # Fitting of a Negative Binomial distribution on the shuffles is only relevant for classical shuffle, not Markov. We saw experimentally that Markov shuffles do not fit the Neg Binom model.
    nofit = False
    if use_markov_shuffling: nofit = True

    # If there was only a single file in bedsB, just pass `all_intersections` to osc.stats_single
    # as the `all_intersections_for_this_combi` object
    if len(bedsB) == 1:
        # Precompute the true intesections between bedA and bedsB once and for all
        true_intersection = osc.compute_true_intersection(bedA, bedsB)

        result = osc.stats_single(all_intersections_for_this_combi=all_intersections,
                                  true_intersection=true_intersection, ft_type=ft_type,
                                  nofit=nofit, draw_histogram=draw_histogram)

    # Otherwise we must call osc.stats_multiple_overlap() which will split `all_intersections` into
    # one object per relevant combination (see function documentation and source for more details)
    if len(bedsB) > 1:
        message('More than one set of regions was provided. Performing statistics on multiple overlaps.')

        result = osc.stats_multiple_overlap(all_overlaps=all_intersections,
                                            bedA=bedA, bedsB=bedsB, all_feature_labels=ft_type,
                                            nb_threads=nb_threads, nofit=nofit,
                                            multiple_overlap_target_combi_size=multiple_overlap_target_combi_size,
                                            multiple_overlap_max_number_of_combinations=multiple_overlap_max_number_of_combinations,
                                            multiple_overlap_custom_combis=multiple_overlap_custom_combis,
                                            draw_histogram=draw_histogram, exact=exact)
        # ft_type, in this case, should be a list of the respective names of all files in bedsB

        # NOTE : in this case, `result` is a dictionary of results giving one 'result'
        # object per interesting combination. This will be separated into the relevant
        # results objects in the main ologram.py code

    # Just in case, explicitly free memory.
    # Theoretically it should have been emptied (pop) by the multiple overlap function above
    del all_intersections
    gc.collect()

    grand_stop = time.time()

    message('--- Total time : ' + str(grand_stop - grand_start) + ' s for feature(s) : ' + str(ft_type) + ' ---')
    message('Total time does not include BED reading, as it does not scale with batch size.', type='DEBUG')

    # This object may be one result or several. Will be treated later (see ologram.py).
    return result
