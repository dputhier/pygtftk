"""
Compute overlap statistics on shuffled sets.
Called by overap_stats_shuffling.compute_overlap_stats(), hence the name.

Those functions tend to take as input lists of shuffles and output statistics.

Author : Quentin Ferré <quentin.q.ferre@gmail.com>
"""

import functools
import gc
import multiprocessing
import time
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import concurrent.futures as cf
from multiprocessing import Pool
import bisect
import copy
import hashlib
import ctypes
import random

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import nbinom

from pygtftk.stats import negbin_fit as nf
from pygtftk.stats import beta as pbeta
from pygtftk.stats.intersect.modl import dict_learning as dl
from pygtftk.stats.intersect.overlap import overlap_regions as oc
from pygtftk.stats.multiprocessing import multiproc as mpc

from pygtftk.utils import make_tmp_file
from pygtftk.utils import message


# ----------------------------- Useful wrappers ------------------------------ #

def compute_true_intersection(bedA, bedsB):
    """
    Returns the custom-computed true intersection between bedA and all in bedsB combined, where
    bedA is a pybedtools.BedTool object and bedsB is a list of such objects.

    Returns also the intersection flags.
    """

    # bedA must be a single file, but bedsB must be a list
    if not isinstance(bedsB, list): raise ValueError(
        "compute_true_intersection was passed a bedsB which is not a list.")

    # Convert bedA and bedsB from pyBedtools files to lists. This is slow, but only called once per analysis.
    def pybedtool_to_pythonlist(bed): return bed.to_dataframe().iloc[:, :3].values.tolist()

    bedA = pybedtool_to_pythonlist(bedA)
    bedsB = [pybedtool_to_pythonlist(bedB) for bedB in bedsB]

    # Pass them to our custom algorithm to get the true intersections
    beds = tuple([bedA] + bedsB)  # Process all beds at once
    all_chrom = set([line[0] for line in bedA])  # Work only on all chromosomes of bedA


    # Sanity check: if all_chrom is an empty list, throw an error
    if not all_chrom:
        msg = "Trying to call `find_intersection` with an empty `all_chrom`. This likely means there are no common chromosomes between your query and reference sets, or that one of the query or references is empty."
        message(msg, type="ERROR")

    true_intersection = oc.find_intersection(beds, all_chrom, return_flags=True)

    return true_intersection


def compute_stats_for_intersection(myintersect):
    """
    Wrapper to compute all stats we want on a single intersect result object.
    The argument (myintersect) is a single bedfile, either as a pybedtools intersect
    result, or a list of tuples/lists like [['chr1',10,20],['chr2',30,40]].
    """
    bp_overlap = [x[2] - x[1] for x in myintersect]

    # Override : if there were no overlaps, return [0] instead of [] to prevent NaNs cropping up later
    if bp_overlap == []: bp_overlap = [0]

    intersect_nb = len(myintersect)
    stats = (bp_overlap, intersect_nb)

    return stats

    # ------------------------ Notes for improvement ------------------------- #
    # TODO: the true way to be fast would be to write this in Cython
    # NOTE: careful, chromomosomes in the intersections are now encoded strings, 
    # not raw strings. So be careful if stats are to be made on them.


################################################################################
# ----------------------- Process sets of intersections ---------------------- #
################################################################################

def merge_consecutive_intersections_in_intersections_list(inter_list):
    """
    Merges the consecutive intersections in a list-of-intersections object.

    Note that this will discard the flags (4th element, np array) and any other element in the list of lists.

    This is NOT an in-place operation, and returns the new list.
    """

    nb_intersections = len(inter_list)

    # This will replace the old shuffles
    new_list = []

    was_written = True # Begin as if we had just written something

    current_chrom = ""
    current_start = -1
    current_end = -1

    pending_chrom = ""
    pending_start = -1
    pending_end = -1

    for i in range(nb_intersections):

        current_chrom = inter_list[i][0]
        current_start = inter_list[i][1]
        current_end = inter_list[i][2]

        # If the end of the list has been reached...
        try:
            next_chrom = inter_list[i+1][0]
            next_start = inter_list[i+1][1]
            next_end = inter_list[i+1][2]
        except:
            next_chrom, next_start, next_end = None, None, None


        # Resume recording the current intersections if one was just written
        if was_written:
            pending_chrom = current_chrom
            pending_start = current_start
            pending_end = current_end
            was_written = False

        # If same chrom and my end == next beginning
        if (current_chrom == next_chrom) and (current_end == next_start):
            pending_end = next_end
        # Otherwise just add the pending
        else:
            new_list.append((pending_chrom, pending_start, pending_end))
            was_written = True

    return new_list


def merge_consecutive_intersections_in_all_overlaps_lists(aiqc):
    r"""
    Merges the consecutive intersections in a list-of-intersection-batches object.

    This is an in-place operation: it will not return a new list, but change the original !

    Example:

    >>> all_intersections_queried_for_this_combi = [ [(b'chr1',1,100),(b'chr1',100,200),(b'chr1',200,300)], [(b'chr1',100,200),(b'chr1',200,300),(b'chr1',600,700)], [(b'chr1',100,200),(b'chr2',200,300)] ]
    >>> expected = [ [(b'chr1',1,300)], [(b'chr1',100,300),(b'chr1',600,700)], [(b'chr1',100,200),(b'chr2',200,300)] ]
    >>> merge_consecutive_intersections_in_all_overlaps_lists(all_intersections_queried_for_this_combi)
    >>> assert all_intersections_queried_for_this_combi == expected

    """
  
    # For each shuffle...
    nb_shuffles = len(aiqc)

    for s in range(nb_shuffles):

        # Record this new list. Overwrite the previous one to save memory.
        aiqc[s] = merge_consecutive_intersections_in_intersections_list(aiqc[s])



def stats_single(all_intersections_for_this_combi, true_intersection,
                 ft_type='some feature', nofit=False, 
                 this_combi_only=None, draw_histogram=False):
                 #was_directly_passed_stats = False):
    """
    Compute statistics such as total number of overlapping base pairs for a given feature.

    :param intersections_for_this_combi: an object returned by compute_all_intersections_minibatch giving all intersections between the shuffled bedA and bedsB.
    :param true_intersection: the result of overlap_stats_compute.compute_true_intersection(bedA, bedsB) where bedA is the query and bedB is the list of the bed files between whose shuffles the aforementioned intersections have been computed. Used here to calculate the true intersections between them and calculate a Neg Binom p-value.
    :param ft_type: for debug messages, which feature/combi are we currently processing ?
    :param nofit: if True, does not do Negative Binomial fitting
    :param this_combi_only: a list of flags (e.g. [1,0,0,1]) corresponding to expected flags in the intersections, one per file (see find_intersection() source and documentation). If not None, we will consider only intersections that have this flag for the number of true intersections and true overlapping basepairs
    :param draw_histogram: if True, draws a temp file histogram for each combi
    """
 
    message('Processing overlaps for ' + ft_type, type='DEBUG')

    start = time.time()

    # Merge consecutive intersections to prevent double counts due to flags changing in the middle of an intersection (when a new set begins or ceases to overlap)
    merge_consecutive_intersections_in_all_overlaps_lists(all_intersections_for_this_combi)

    # Compute the statistics
    stats = [compute_stats_for_intersection(myintersect) for myintersect in all_intersections_for_this_combi]

    # Unpack the stats.
    bp_overlaps = [s[0] for s in stats]  # Those are the individual overlap lengths, a list of lists
    summed_bp_overlaps = [sum(x) for x in bp_overlaps]  # Sum by shuffle
    intersect_nbs = [s[1] for s in stats]


    # ------------------------------------------------------------------------ #
    # NOTE To speed up later, the most efficient way would be to compute
    # statistics directly in the Cython code.
    # TODO: Use this if Cython returned directly the stats.
    # if was_directly_passed_stats:
    #     bp_overlaps, summed_bp_overlaps, intersect_nbs = all_intersections_for_this_combi
    # ------------------------------------------------------------------------ #

    # Merge consecutive intersections in the true one as well
    true_intersection = merge_consecutive_intersections_in_intersections_list(true_intersection)

    ## True intersection statistics
    true_intersect_nb = len(true_intersection)
    true_bp_overlaps = sum([x[2] - x[1] for x in true_intersection])


    stop = time.time()

    message(ft_type + '- Statistics on overlaps computed in : ' + str(stop - start) + ' s.', type='DEBUG')

    # NOTE FOR IMPROVEMENT : it would be interesting to return the average size
    # of an overlap as well, per shuffle. Since our intersection algorithm returns
    # details about the intersections like `bedtools intersect` would, this could
    # be computed without much hassle.



    # ------ Fitting of a Negative Binomial distribution on the shuffles ----- #
    start = time.time()

    # Fitting can be disabled from the main function (for now, mainly relevant if we used a Markov model instead of a classical one.)
    if nofit:
        ps = pn = -1
        expectation_fitted_summed_bp_overlaps, variance_fitted_summed_bp_overlaps = -1
        expectation_fitted_intersect_nbs, variance_fitted_intersect_nbs = -1

    else:
        # Renaming expectations and variances
        expectation_fitted_summed_bp_overlaps, variance_fitted_summed_bp_overlaps = np.mean(summed_bp_overlaps), np.var(summed_bp_overlaps)
        expectation_fitted_intersect_nbs, variance_fitted_intersect_nbs = np.mean(intersect_nbs), np.var(intersect_nbs)

        # If we were passed an empty `all_intersections_for_this_combi` list, the expectations will be NaN.
        # Manually set them to zero to prevent trying to fit a Negative Binomial on them.
        if len(all_intersections_for_this_combi) == 0:
            expectation_fitted_summed_bp_overlaps = 0
            expectation_fitted_intersect_nbs = 0
            ps = pn = -1

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

    # Compute the p-values using the distribution fitted on the shuffles.
    # Do not do this for the Markov shuffling, as it is likely a multi-variable fit (see notes)
    # We can only use a Neg Binom p-val if we can fit it, and that is not the case for
    # the Markov shuffle or if the expectation is too small.


    # To be used if the fitting is bad: add the empirical p-value (proportion of shuffles with values as extreme)
    empirical_pval_intersect_nb = nf.empirical_p_val(true_intersect_nb, intersect_nbs)
    empirical_pval_bp_overlaps = nf.empirical_p_val(true_bp_overlaps, summed_bp_overlaps)



    # Return -1 for the p-value if the fitting was bad.
    if (ps == -1) | (pn == -1):
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
    message(ft_type + '- Negative Binomial distributions fitted in: ' + str(stop - start) + ' s.',
            type='DEBUG')


    
    # --------------------------------------------------------------------------
    # Draft code for diagnostic plots of the distribution of each statistic in
    # the shuffles. Kept for potential future improvement.
    # --------------------------------------------------------------------------

    # Drawing is computationally expensive, make it optional
    if draw_histogram:

        # This may return an error for very long ft_type strings, so wrap it in a try-except block
        # TODO: Fix this properly
        # TODO: currently, those files remain in /tmp because this function is subprocessed.
        #I must prepare a temp file manager like make_tmp_file_pool() to keep them in the directory specified by -K
        start = time.time()
        try:
            hist_S = make_tmp_file(prefix='histogram_' + ft_type + '_S_sum_by_shuffle', suffix='.png')

            ## Number of overlapping base pairs
            # Sum by shuffle
            plt.figure()
            mean = expectation_fitted_summed_bp_overlaps
            var = variance_fitted_summed_bp_overlaps
            r = mean ** 2 / (var - mean)
            p = 1 / (mean / r + 1)

            # Plot the histogram.
            BINS = 200  #TODO: put it higher conditionally 
            de = plt.hist(summed_bp_overlaps, bins=BINS)[1]
            # Plot the PDF.
            xmin, xmax = min(de), max(de)
            x = np.linspace(xmin, xmax, BINS)

            # Steps should always be the same between bins, except for rounding errors
            steps = [0] + [x[i] - x[i - 1] for i in range(1, len(x))]

            try:
                d = [nbinom.cdf(x[i], r, p) - nbinom.cdf(x[i - 1], r, p) for i in range(1, len(x))]

                # Position the 'd' correctly : move one pace right, and
                # put it in te middle of each bar
                x = [x[i] - steps[i]/2 for i in range(len(x))]

                d = np.array([0] + d) * len(summed_bp_overlaps)
            
            except:
                d = [0] * BINS

            plt.plot(x, d, 'k', linewidth=1.2)
            plt.savefig(hist_S.name)
            plt.close()

        except:
            pass

        stop = time.time()
        message(ft_type + '- Drew histogram in : ' + str(stop - start) + ' s.',
                type='DEBUG')
  

    

    ## Now return the result as a dictionary of statistics for this calculation.
    # NOTE Be careful to use the same order as result_abort, in overlap_stats_shuffling.py !
    result = OrderedDict()

    # Number of intersections
    result['nb_intersections_expectation_shuffled'] = '{:.2f}'.format(expectation_fitted_intersect_nbs)
    result['nb_intersections_variance_shuffled'] = '{:.2f}'.format(variance_fitted_intersect_nbs)

    result['nb_intersections_negbinom_fit_quality'] = '{:.5f}'.format(pn)

    if expectation_fitted_intersect_nbs == 0:
        ni_fc = true_intersect_nb  # Do not divide by zero ! Use the true value as fold change
    else:
        ni_fc = true_intersect_nb / expectation_fitted_intersect_nbs

    if ni_fc != 0: ni_fc = np.log2(ni_fc)  # Apply log transformation
    result['nb_intersections_log2_fold_change'] = '{:.5f}'.format(ni_fc)

    result['nb_intersections_true'] = true_intersect_nb
    result['nb_intersections_pvalue'] = '{0:.4g}'.format(pval_intersect_nb)

    # Summed number of overlapping basepairs
    result['summed_bp_overlaps_expectation_shuffled'] = '{:.2f}'.format(expectation_fitted_summed_bp_overlaps)
    result['summed_bp_overlaps_variance_shuffled'] = '{:.2f}'.format(variance_fitted_summed_bp_overlaps)

    result['summed_bp_overlaps_negbinom_fit_quality'] = '{:.5f}'.format(ps)

    if expectation_fitted_summed_bp_overlaps == 0:
        sbp_fc = true_bp_overlaps  # Do not divide by zero ! Use the true value as fold change
    else:
        sbp_fc = true_bp_overlaps / expectation_fitted_summed_bp_overlaps

    if sbp_fc != 0: sbp_fc = np.log2(sbp_fc)  # Apply log transformation
    result['summed_bp_overlaps_log2_fold_change'] = '{:.5f}'.format(sbp_fc)

    result['summed_bp_overlaps_true'] = true_bp_overlaps
    result['summed_bp_overlaps_pvalue'] = '{0:.4g}'.format(pval_bp_overlaps)


    # Remember the order of the combinations (meaning the number of open sets)
    # NOTE May be different when intra-set overlaps are introduced
    my_order = np.sum(this_combi_only)
    if my_order is None: my_order = 1
    result['combination_order'] = str(my_order)

    # Also put the empirical p-value (proportion of shuffles where a value as extreme is found)
    # NOTE Added as last columns so it does not bother the existing column tests
    result['nb_intersections_empirical_pvalue']   = '{0:.4g}'.format(empirical_pval_intersect_nb)
    result['summed_bp_overlaps_empirical_pvalue'] = '{0:.4g}'.format(empirical_pval_bp_overlaps)

    # Use an ad-hoc beta approximation for the beta-binomial
    try:
        beta_pval_bp_overlaps = pbeta.beta_pval(true_bp_overlaps, summed_bp_overlaps)
    except:
        # Return -1 if any problem
        # TODO Remove this band-aid
        beta_pval_bp_overlaps = -1
    result['beta_summed_bp_overlaps_pvalue_ad_hoc_for_deep_sampling_only'] = '{0:.4g}'.format(beta_pval_bp_overlaps)

    #message(ft_type + '- Result dump : ' + str(result), type='DEBUG')

    return result


# -------------------------- Multiple overlap sets --------------------------- #


## ------------ Helper objects

class Counter(object):
    """
    Multiprocessing-compatible counter that can be incremented.
    """
    def __init__(self, initval=0):
        self.val = multiprocessing.Value('i', initval)
        self.lock = multiprocessing.Lock()

    def increment(self):
        with self.lock: self.val.value += 1

    def value(self):
        with self.lock: return self.val.value


def get_items_by_indices_in_list(indices, mylist):
    r"""
    Quick way to perform multiple indexing on a list.
    Unlike the operator.itemgetter method, will always return a list

    Example:
    >>> assert get_items_by_indices_in_list(indices = [0], mylist = [1,2,3]) == [1]

    """
    accessed_mapping = map(mylist.__getitem__, indices)
    return list(accessed_mapping)


class HashableArray(np.ndarray):
    """
    A subclass of NumPy's ndarray that can be used as dictionary key.
    A dictionary of such arrays will consume less RAM than a tuple, especially when pickled.

    This is however immutable, and should not be used on arrays that you intend to modify.

    NOTE It is hashed at creation, so creation can take a bit of time.
    Use a tuple if that is unacceptable.
    """
    def __new__(cls, input_array): 
        if isinstance(input_array, HashableArray):
            return input_array
        else:
            return np.asarray(input_array).view(cls)

    def __array_finalize__(self, obj):
        # Set to read-only to prevent modifications
        self.flags.writeable = False

        # Hash the array
        myhash = int(hashlib.sha1(bytes(self)).hexdigest(), 16)

        # Set myhash to a hash of the byte representation if not already present
        if obj is None: return
        self.myhash = getattr(obj, 'myhash', myhash)  

    def __hash__(self):
        # Custom attributes, such as myhash, are lost during multiprocessing.
        # If that happens, recalculate the hash as needed and store it again
        try:
            return self.myhash
        except AttributeError:
            self.myhash = int(hashlib.sha1(bytes(self)).hexdigest(), 16)
        return self.myhash

    def __eq__(self, other):
        # Force other to cast to a HashableArray
        other = HashableArray(other) 

        # Compare hashes
        return self.__hash__() == other.__hash__()


def which_combis_to_get_from(combi, all_possible_combis, exact):
    r"""
    To build an index later, determine the supplementary combis : meaning, for each
    combination, which other combinations needs to be queried.

    Relevant if we desire inexact or exact combis;
    Reminded : a combi A is a "inexact" match for combi B if all flags of combi B
    are also present in combi A, along with potentially others.

    Returns an index giving the supplementary combinations.

    Example :

    >>> all_combis = [(1,1,0),(1,1,1),(0,0,1)]
    >>> c,r = which_combis_to_get_from((1,1,0), all_combis, exact = False)
    >>> assert r.tolist() == [1]

    """

    combi = tuple(combi) # Force conversion to tuple

    # For all combinations, in the same order, get a list of the IDs of the matching ones
    matching_vector = list()

    i = 0
    for c in all_possible_combis:
        # IMPORTANT : get all the combis that match the original combi to create
        # the index, but do not add the original combi itself to the index.
        # Indeed, the DictionaryWithIndex object always adds the elements of
        # the original combi itself, so they should not be in the index.
        # We pass the current exact flag : this way, it works whether exact is False or True,
        # returning an empty index in the latter case

        # NOTE the combinations in all_possible_combis are converted to tuples
        # before this function call

        if oc.does_combi_match_query(c, combi, exact=exact):
            if c != combi: 
                matching_vector += [i]
        i = i+1

    # The matching vector is a Python list, convert it to numpy array to save some RAM
    matching_vector = np.asarray(matching_vector, dtype= np.uint64)
    message("Computed exactitude index for combi "+str(combi), type="DEBUG")

    return combi, matching_vector




from pygtftk.stats.intersect import dummy
# Store all the global variables into a dummy Python module.
def initProcessMultiproc(shared_array_all_combis_base, combi_nb, combi_size, nb_combis_done, tot_number_interesting_combis):
    dummy.shared_array_all_combis_base = shared_array_all_combis_base
    dummy.combi_nb = combi_nb
    dummy.combi_size = combi_size
    dummy.nb_combis_done = nb_combis_done
    dummy.tot_number_interesting_combis = tot_number_interesting_combis


def index_all_these(combis_to_index, exact):
    """
    Helper function to run the global numpy analogue of 
    which_combis_to_get_from on a list of combis.

    This is not a pure function and relies on several external parameters
    that will be named as globals later, before this function is called.
    Those parameters are stored in the "dummy" module.

    NOTE : the oc.NPARRAY_which_combis_match_with already integrates a check
    and will not return the same index as the query
    """
    mappings = list()

    # Recreate the all_combis arrays from the shared buffer
    # This is a GLOBAL and was not passed to the worker
    all_combis = np.frombuffer(dummy.shared_array_all_combis_base, dtype=ctypes.c_uint)
    all_combis = all_combis.reshape(dummy.combi_nb, dummy.combi_size)

    for combi in combis_to_index:
      
        # Call the exactitude computation
        matching_list = oc.NPARRAY_which_combis_match_with(all_combis, combi, exact)
        matching_vector = np.asarray(matching_list, dtype= np.uint64)

        mappings += [(combi, matching_vector)]


        # Debug prints
        message("Index combination "+str(combi), type = "DEBUG")
        dummy.nb_combis_done.increment()
        message("Combis done: "+str(dummy.nb_combis_done.value())+'/'+str(dummy.tot_number_interesting_combis), type = "DEBUG")

    return mappings



class CombinationExactMapping():
    r"""
    A wrapper containing a list of all combinations and a vector giving, for each
    combination, the ID of all other combinations that are a match for it.

    For example, if working with inexact combinations, [1,1,1] is a match when querying [1,1,0]

    This is done to save RAM when working with very long combinations.

    Example :

    >>> all_combis = [(1,1,0),(1,1,1)]
    >>> mapping = {(1,1,0):[1]}
    >>> cm = CombinationExactMapping(all_combis, mapping)
    >>> assert cm.get_all((1,1,0)) == [(1,1,1)]
    >>> assert cm.get_all((1,1,1)) == []

    """

    def __init__(self, all_combis, mapping):
        # The mapping should be a list of the IDs of combinations to return, like in the example
        self.all_combis = all_combis
        self.mapping = defaultdict(list, mapping)

    def get_all(self, combi):
        ids = self.mapping[combi]
        return [self.all_combis[i] for i in ids]

    def __repr__(self):
        return "all_combis="+self.all_combis.__repr__()+";"+"mapping="+self.mapping.__repr__()





class DictionaryWithIndex():
    r"""
    A wrapper allowing to query a dictionary so that an item can also return values from several other items.
    The dictionary will return concatenated values from several keys when asked for.

    It will return values from itself PLUS everything in the index

    Specifically designed for lists of overlaps in OLOGRAM-MODL, meaning
    the final result will be a list with the same number elements, whose i-th element for i in 1..n
    will be the concatenated i-th elements of all candidate lists


    Example :

    >>> all_combis = ['A','B','C']
    >>> mapping = {'A':[1,2]}
    >>> index = CombinationExactMapping(all_combis, mapping)
    >>> d = { 'A':[[1,1],[1],[1,1]], 'B':[[2],[2,2],[2]], 'C':[[3,3],[3],[3]] }
    >>> di = DictionaryWithIndex(d, index)
    >>> assert di.get_all('A') == [[1,1,2,3,3], [1,2,2,3], [1,1,2,3]]
    >>> assert di.get_all('B') == [[2],[2,2],[2]]

    """

    def __init__(self, data, index, data_default_factory=lambda: [],
        will_store_an_all_overlaps_object = False):
        """
        :param data: The dictionary to be wrapped
        :param index: a CombinationExactMapping object giving for each key in data the *other* keys to be also used when calling each key
        :param data_default_factory: A function producing the default value of data after wrapping. Defaults to `lambda:[]`
        """      

        # Index and data should both be defaultdict to handle missing cases
        self.index = index

        ## If we are storing all_overlaps, use a special dedicated structure
        if will_store_an_all_overlaps_object:

            self.using_cython_all_overlaps = True

            # Store the data in a dedicated structure
            # This is a Cython structure with an underlying NumPy array that can be
            # accessed by several processes
            self.data = mpc.PYWRAPPER_DictionaryOfOverlapsWithSharedNparrayStorage(
                data, data_default_factory)
        else:
            self.using_cython_all_overlaps = False
            self.data = defaultdict(data_default_factory, data)

        gc.collect()  # Force garbage collecting in case of large dictionaries

    def get_all(self, key):
        # Retrieve all combinations from the index
        all_keys_to_get = self.index.get_all(key)

        # Start with the value for the key itself
        # Important to use copy here to not modify the original
        # Rq : DictionaryOfOverlapsWithSharedNparrayStorage always returns a copy anyways, copy.deepcopy is not necessary here
        if self.using_cython_all_overlaps: val_concat = self.data[key]
        else: val_concat = copy.deepcopy(self.data[key])

        # For each other key to be added...
        for k in all_keys_to_get:

            buffer = self.data[k]

            # Merge all elements of the supplementary key to the original key
            for i in range(len(self.data[key])):
                val_concat[i] += buffer[i]

        return val_concat

    def get_simple_concatenation(self, key):
        """
        For the true overlaps. Do not merge individual shuffles, there are no shuffles to be merged.
        Simply concatenate the lists.
        """
        all_keys_to_get = self.index.get_all(key)
        if self.using_cython_all_overlaps: val_concat = self.data[key]
        else: val_concat = copy.deepcopy(self.data[key])
        for k in all_keys_to_get: val_concat += self.data[k]
        return val_concat

    def __repr__(self):
        return "index="+self.index.__repr__()+";"+"data="+self.data.__repr__()


def get_index_if_present(mylist, x):
    r"""
    Locate the leftmost value exactly equal to x in a sorted list, otherwise returns None

    Example:

    >>> L = [1,2,5,6,8]
    >>> assert get_index_if_present(L, 2) == 1
    >>> assert get_index_if_present(L, 4) == None

    """
    i = bisect.bisect_left(mylist, x)
    if i != len(mylist) and mylist[i] == x:
        return i
    return None



class SparseListOfLists:
    r"""
    Container for a list of this type :
        [
            [elements], [], [], [], [other_elements], [], [], ...
        ]

    Meaning, a list of lists where many of the lists inside will be empty.
    This helps save RAM.

    Example:

    >>> l = SparseListOfLists()
    >>> l.put([])
    >>> l.put(['Hello'])
    >>> l.put([])
    >>> assert [i for i in l] == [[],['Hello'],[]]
    >>> assert l[0] == []
    >>> assert l[1] == ['Hello']
    >>> l[0] = ['Ha']
    >>> l[1] = ['Ho']
    >>> assert [i for i in l] == [['Ha'],['Ho'],[]]
    >>> assert list(l) == [i for i in l]

    """

    def __init__(self, starting_index = 0):
        self.elements = []
        self.full_slots = [] # Which slots do these overlaps correspond to ?
        self.current_index = starting_index

    def put(self, element):
        """
        Add an element at the tail of the SparseListOfLists
        """
        # If `element` is an empty list, False or None, do not add it and simply
        # increment the index.
        # In effect, we have added `[]` to the list
        if not element: 
            pass
        else: 
            self.elements.append(element)
            self.full_slots.append(self.current_index)
        self.current_index += 1


    def __getitem__(self, key):
        if key >= self.current_index:
            raise IndexError

        pos = get_index_if_present(self.full_slots, key)
        if pos is not None: 
            return self.elements[pos]
        else:
            return list()   # Return an empty list if the slot was not occupied


    def __setitem__(self, key, item):
        """
        Modify an element that is already in the SparseListOfLists
        """

        pos = get_index_if_present(self.full_slots, key)

        # If modifying a slot that contains an element
        if pos is not None:
            self.elements[pos] = item

        # If modifying a slot that contains an empty list
        else:
            # If adding a non-empty list:
            if item: 
                # Get the position where a given element should be inserted in a sorted list to keep the list sorted
                p = bisect.bisect_left(self.full_slots, key)
                # Add the elements there
                self.full_slots.insert(p, key)
                self.elements.insert(p, item)
            # If adding an empty list, do nothing
            else:
                pass

        # If adding to a slot that has not been reached
        if key >= self.current_index:
            raise IndexError


    def __len__(self):
        return max(self.current_index, 0)

    def __iter__(self):
        self.reading_index = 0
        return self

    def __next__(self):
        if self.reading_index < self.current_index:
            result = self.__getitem__(self.reading_index)
            self.reading_index +=1
            return result 
        else:
            raise StopIteration

    def __repr__(self):
        return str([i.__repr__() for i in self])




## ------ Main function

def stats_multiple_overlap(all_overlaps, bedA, bedsB, all_feature_labels, nb_threads=8, nofit=False,
                           # Parameters for the finding of interesting combis
                           multiple_overlap_target_combi_size=None,
                           multiple_overlap_max_number_of_combinations=None,
                           multiple_overlap_custom_combis=None,
                           draw_histogram=False,
                           exact=False):
    """
    Instead of returning one set of overlap stats per query type (ie. exons, gens, bedfile1, bedfile2, etc...)
    it will return one per multiple overlap (ie. here there was a peak for bedfile1+gene, or bedfile1+bedfile2+gene, etc.)

    Only do so for "common/interesting combis", found by dict learning on original data.

    :param all_overlaps: The list of all overlaps computed for all shuffles
    """

    ## ------  Parameters for the finding of interesting combis

    # Exactitude: should an intersection of A+B+C count towards looking for A+B?
    # It depends on a user-specified parameter called `exact`.
    # By default, yes, meaning we look for "inexact" combis.
    # We will look instead for exclusive combis if the user manually specifies 
    # multiple_overlap_target_combi_size equal to the number of sets.


    # ------------- Override combinations
    # If custom_combis are set, skip all the combination mining above.
    # Directly read a text file.
    if multiple_overlap_custom_combis is not None:
        message('Working on custom combinations for multiple overlap.')

        # Read NumPy matrix and cast to regular Python list
        interesting_combis_matrix = np.loadtxt(multiple_overlap_custom_combis.name, dtype=int)
        interesting_combis_matrix = np.atleast_2d(interesting_combis_matrix)    # Ensure it is always 2D, even if only one element

        interesting_combis = [tuple(combi) for combi in interesting_combis_matrix]





    # Only mine for combinations if custom combinations were NOT specified
    else:

        # --------------- Mining for interesting combinations ---------------- #

        true_intersection = compute_true_intersection(bedA, bedsB)
        flags_matrix = np.array([i[3] for i in true_intersection])

        # Keep only combis WITH query (first element is not 0)
        flags_matrix = flags_matrix[flags_matrix[:, 0] != 0]

        import time  # Needed to re-import here for some reason ?
        start = time.time()

        """
        # TODO Add the abundances as a debug message ?
        u, count = np.unique(flags_matrix, return_counts = True, axis=0)
        count_sort_ind = np.argsort(-count)
        df = pd.DataFrame(u[count_sort_ind])
        df['COUNT'] = count[count_sort_ind] 
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):   print(df)
        """


        # If the combi size (`multiple_overlap_target_combi_size`) and the max 
        # number of combinations (`multiple_overlap_max_number_of_combinations`) 
        # were not set manually, they default to -1.

        # Default multiple_overlap_max_number_of_combinations is -1, meaning
        # MODL should not be applied and we should return all combinations
        if multiple_overlap_max_number_of_combinations == -1:
            message('Multiple-overlap combinations number was not restricted, skipping itemset mining step.')

            # In this case the "interesting combis" are all combis encountered in the true data
            interesting_combis, _ = np.unique(flags_matrix, axis=0, return_counts=True)

            # If multiple_overlap_target_combi_size however is not -1, still apply the restriction
            if multiple_overlap_target_combi_size != -1:
                interesting_combis = np.array(
                    [combi for combi in interesting_combis if sum(combi) <= multiple_overlap_target_combi_size])
                message('Restricted to combinations of maximum size' + str(multiple_overlap_target_combi_size) + '.')


        else:
            message('Mining for multiple-overlap interesting combinations using dictionary learning.')

            # Use the MODL algorithm
            combi_miner = dl.Modl(flags_matrix,
                                  multiple_overlap_target_combi_size,
                                  multiple_overlap_max_number_of_combinations,
                                  nb_threads)  # Critical nb_threads for this

            """
            # NOTE for improvement : any combination mining algorihm could be slotted here instead !
            # TODO: Perhaps offer the option to use our draft of implementation of apriori
            combi_miner = Apriori(flags_matrix, min_support)
            """

            interesting_combis = combi_miner.find_interesting_combinations()

            # TODO: Other possibility, simply take the most common combis
            # This can also be done by the ologram_modl_treeify plugin
            """
            all_combis, counts_per_combi = np.unique(flags_matrix, axis=0, return_counts = True)
            most_frequent_combis_idx = (-counts_per_combi).argsort()[:multiple_overlap_max_number_of_combinations]
            interesting_combis = all_combis[most_frequent_combis_idx]
            """

            stop = time.time()
            message(
                'Interesting combinations found via dictionary learning in : ' + str(stop - start) + ' s among ' + str(
                    all_feature_labels) + '.', type='DEBUG')

    # TODO: Make this debug print prettier
    message("Interesting combinations were " + str(interesting_combis), type="DEBUG")
    message("There were "+str(len(interesting_combis))+" interesting combinations.")

    ## Interesting combis sanity checks
    interesting_combis_final = []
    for combi in interesting_combis:
        combi_as_list = list(combi)

        # Enforce the presence of the query in all combinations of interest
        combi_as_list[0] = 1
        interesting_combis_final += [tuple(combi_as_list)]

    interesting_combis = interesting_combis_final

    # Should never happen, but just in case.
    if len(interesting_combis) == 0:
        raise ValueError("No combination of interest found. Try disabling --multiple-overlap-max-number-of-combinations or increasing the number of shuffles.")

    ## Precompute the true intesections between bedA and bedsB once and for all
    message("Computing all true intersections...")
    true_intersection = compute_true_intersection(bedA, bedsB)

    # stats_single() requires a list of all true intersections be passed for each combi.
    # Split the true intersections and prepare those

    # Read intersection flag as tuples for the interesting combis only

    true_intersections_per_combi = defaultdict(list)
    for inter in true_intersection:
        # NOTE Enforcing np.uint32 type everywhere
        combi = HashableArray(
            np.array(inter[3], dtype = np.uint32)
        ) 
        
        true_intersections_per_combi[combi] += [inter]






    # ----------- Splitting all the intersections for all shuffles ----------- #
    
    global overlaps_per_combi   # Reserve global

    ## Split the list of overlap regions per combis

    import time  # Needed to re-import here for some reason ?
    time.sleep(0.1); gc.collect() # Garbage collect

    start = time.time()

    message("Splitting all overlaps computed for all shuffles by combination...")


    ## ---- Split all overlaps computed for all shuffles by combination

    # NOTE I use pop here, so we begin with the last element, since pop() is O(1).
    # It is necessary to do so to save RAM.
    # WARNING : this will empty the original all_intersections list, since it was passed
    # by reference (?) as all_overlaps !
    total_to_split = len(all_overlaps)

    # By default, missing keys have an empty list
    overlaps_per_combi = defaultdict(lambda: SparseListOfLists(0))

    starting_index = 0


    while all_overlaps:
        # Get a minibatch by popping
        intersections_for_this_shuffle = all_overlaps.pop()

        # For each minibatch, make 'filtered' minibatches for each possible overlap flag, containing only overlaps with these
        # flags. Then, add each filtered minibatch one at a time to the overlaps_per_combi dictionary.
        filtered_minibatches = defaultdict(list)

        while intersections_for_this_shuffle:
            overlap = intersections_for_this_shuffle.pop()

            # What are the flags (ie. sets present) for this overlap ? For example, "A+C but not B" is (1,0,1)
            flags_for_this_overlap = HashableArray(overlap[3]) # Make a Hashable array as lists cannot be dict keys
            # NOTE This should have np.uint32 type

            ## Add exact matches
            # Don't need to remember the flags since I group them by flag anyways 
            filtered_minibatches[flags_for_this_overlap] += [overlap[0:3]]

        ## Now add the filtered minibatches as elements of a list
        # overlaps_per_combi then contains lists of lists : one list per original
        # minibatch which contains a list of all overlaps WITH THIS FLAG for this minibatch

        # For the combinations already encountered and also encountered here, add the filtered batch
        for combi, filtered_batch in filtered_minibatches.items():
            overlaps_per_combi[combi].put(filtered_batch)

        # For the combinations already encountered but NOT encountered here, add an empty batch
        combis_already_encountered_but_not_here = [
            c for c in overlaps_per_combi if c not in filtered_minibatches
        ]
        for c in combis_already_encountered_but_not_here:
            overlaps_per_combi[c].put([])


        # For the combinations not yet encountered, update the default factory to include an additional empty batch     
        starting_index += 1
        new_factory_expression = 'lambda: SparseListOfLists('+str(starting_index)+')'
        overlaps_per_combi.default_factory = eval(new_factory_expression)

        message("Splitting done for "+str(starting_index)+'/'+str(total_to_split), type = "DEBUG")


    # Cleanup
    del intersections_for_this_shuffle
    del filtered_minibatches
    del combis_already_encountered_but_not_here
    del all_overlaps

    ## Partial matches
    # We have registered all exact combis, now add partial matches. Partial matches are defined as "including all flags",
    # meaning (1,1,1,0) will be a match if we query (1,1,0,0) since it contains all its flags, but will not be a match for (1,0,0,1)
    # Iterate over all combinations, and if they match add contents of combi B to combi A

    ## -------------- Compute the index of combis to be fetched ------------- ##
    # Relevant for partial matches.

    time.sleep(0.1); gc.collect() # Garbage collect

    message("Computing index of exact/inexact combinations...")

    # Compute the index against all combis found, but also against the interesting combis and all true combis. Relevant mostly for inexact combis.
    # Also enforce type when creating the list
    all_combis = set(
        list(overlaps_per_combi.keys()) + interesting_combis + list(true_intersections_per_combi.keys())
    )
    all_combis = tuple(set(
        [
            HashableArray(
                np.asarray(c, dtype=np.uint32) 
            ) for c in all_combis
        ]
    ))

    
    # Reserve some globals
    global shared_array_all_combis_base  # It must be a global variable to be shared by the processes

    global combi_nb
    global combi_size

    global nb_combis_done
    global tot_number_interesting_combis

    tot_number_interesting_combis = len(interesting_combis)
    combi_nb = len(all_combis)
    combi_size = len(all_combis[0])

    message("We will index " + str(tot_number_interesting_combis)+"*"+str(len(all_combis)) + " combinations. This can be very long (minutes, hour) for longer combinations.")


    # NOTE We do not need to index all combis : the only ones that will ever be queried are the `interesting_combis`. Those need to 
    # be fully indexed against `all_combis.
    # However but not every combination in `all_combis` : we do not care what we would need to get if we were to query a combination C
    # that is in `all_combis`, but not in `interesting_combis`


    ## Create a shared array containing all combis, to be passed during multiprocessing
    # NOTE This must NOT be passed as an argument to the processes ! Instead it must
    # be a global variable that they can call upon !
    shared_array_all_combis_base = multiprocessing.Array(ctypes.c_uint, combi_nb*combi_size,
        lock = False) # Must disable the lock to permit shared access. Fine since it is read-only.

    nb_combis_done = Counter(0)

    # Temporary reference to populate it, this WILL NOT BE PASSED
    # Now populate it with the combinations
    shared_array_all_combis = np.frombuffer(shared_array_all_combis_base, dtype=ctypes.c_uint)
    shared_array_all_combis = shared_array_all_combis.reshape(combi_nb, combi_size)
    for i in range(len(all_combis)):
        shared_array_all_combis[i,:] = all_combis[i]

    message("Shared array with all combinations populated...", type = "DEBUG")

    del shared_array_all_combis # Now delete this array so the reference is free, just in case






    ## Initialize the multiprocessing
    pool = ProcessPoolExecutor(nb_threads,
    
        # For the shared variables
        initializer=initProcessMultiproc,
        initargs=(shared_array_all_combis_base, combi_nb, combi_size, nb_combis_done, tot_number_interesting_combis)
    
    )
    futures = list()

    # Divide the interesting combis into as many batches as threads, and remove empty batches
    combi_id = list(range(len(interesting_combis)))
    random.shuffle(combi_id)    # Shuffle so that combis that are longer to process should be more uniformly spread
    batches_of_combis_to_index_id = np.array_split(combi_id, nb_threads)
    batches_of_combis_to_index_id = [b for b in batches_of_combis_to_index_id if len(b)]
    # NOTE Using batches turned out to be critical to performance, otherwise too
    # much time is lost

    ## Now submit the jobs
    while batches_of_combis_to_index_id :

        # Get the corresponding combis with the IDs
        b = batches_of_combis_to_index_id.pop()
        combis_to_index = get_items_by_indices_in_list(b, interesting_combis)      

        # Enforce type
        combis_to_index = [
            HashableArray(
                np.array(combi, dtype = np.uint32)
                ) for combi in combis_to_index
        ]

        # Submit
        futures += [pool.submit(index_all_these, 
                combis_to_index = combis_to_index, 
                exact = exact)]


    # Release the resources as soon as you are done with those, we won't submit any more jobs to you
    pool.shutdown() 
    message("Exactitude computation jobs submitted.")

    
    ## Transpose the results into `mappings`
    mappings = list()
    for future in futures:

        try:
            mappings.append(future.result())
        except cf.process.BrokenProcessPool:
            message('A process in the process pool was terminated abruptly while the future was running or pending. This likely means you ran out of RAM. Try restarting the command with fewer threads or smaller minibatches', type="ERROR")
 
    # Unlist mappings
    mappings = [mapping for sublist in mappings for mapping in sublist]


    # Convert the mappings into a sparser object for storage    
    mappings = CombinationExactMapping(all_combis, dict(mappings))


    ## Cleanup
    del shared_array_all_combis_base

    # Need to reinitialize heap to clear memory
    # TODO: supposedly no longer needed in Python 3.8, double check that
    multiprocessing.heap.BufferWrapper._heap = multiprocessing.heap.Heap()
    gc.collect()

    message("Index unpacked. Repartition of overlaps...")





    ## Finally, create a DictionaryWithIndex object to hold all intersections
    # If exact = False, when querying this dictionary using get_all(c), all 
    # combis that are an inexact match for c will also count.
    overlaps_per_combi = DictionaryWithIndex(overlaps_per_combi, mappings,
                                             data_default_factory=overlaps_per_combi.default_factory,
                                             will_store_an_all_overlaps_object=True) # Use a special Numpy-array backed structure that permits multiprocessing
    # NOTE: Ensure that we overwrite the original object to save memory !

    ## Do the same for the true intersections
    true_intersections_per_combi = DictionaryWithIndex(true_intersections_per_combi, mappings,
                                                       data_default_factory=true_intersections_per_combi.default_factory)
    # We will pass those true_intersections to stats_single and resume as usual

    stop = time.time()
    message('All computed overlaps for the shuffles split by combination in : ' + str(stop - start) + ' s among ' + str(
        all_feature_labels) + '.', type='DEBUG')




    ## ------------------- Enrichment for each combination ------------------ ##
    # Now call stats_single on each.


    ## Result queue
    mana = multiprocessing.Manager()
    result_queue = mana.Queue()
    all_results = dict()  # Final result dict, to be filled when emptying the queue


    ## General idea : spawn a process per batch that will get the overlaps, process them, and move on to the next combi

    def compute_those_stats(combis, result_queue):

        # Quick way to do so only one combination at a time
        for c in combis:

            # ----- For each combination, get its full form, and process it
            combi = interesting_combis[c]
            
            # Convert the combi into a user-friendly string
            indices = [i for i in range(1, len(combi)) if combi[i] != 0]  # 0-based !
            combi_list = [all_feature_labels[i - 1] for i in indices]

            if combi[0] != 0: combi_list = ['Query'] + combi_list

            # Add '...' to the combi if not exact, to show the user that others TFs
            # could also be present in these intersections
            if exact: combi_human_readable = '[' + ' + '.join(combi_list) + ']'
            if not exact: combi_human_readable = '[' + ' + '.join(combi_list) + ' + ... ]'

            message('Will compute statistics for the combination : ' + str(combi_human_readable), type='DEBUG')

            # Convert to a HashableArray because lists cannot be dict keys
            # NOTE Enforcing np.uint32 type everywhere, otherwise the hashes will be
            # different to those in the overlaps_per_combi
            combi_key = HashableArray(
                np.array(combi, dtype = np.uint32)
            )


            ## Collect all shuffles for this combination, taking exactitude into account

            # I have redone this using a new Cython object with a parallel-accessible array
            # NOTE this is a bit slow, but should stillbe faster than pickling and will be distributed across the processes             
            s = time.time()      
            list_overlaps_shuffled_for_this_combi = overlaps_per_combi.get_all(combi_key)     
            e = time.time()
            message("Collected shuffled intersections for " + combi_human_readable + " in " + str(e-s) + "seconds.", type = "DEBUG")

            # TODO: In the future, have Cython compute the statistics directly
            # from the underlying NumPy array
            # list_overlaps_shuffled_for_this_combi = overlaps_per_combi.return_directly_stats_for_all(combi_key)

            true_intersections_for_this_combi = true_intersections_per_combi.get_simple_concatenation(combi_key)

            # Compute the result
            myresult = stats_single(list_overlaps_shuffled_for_this_combi,
                true_intersections_for_this_combi,
                combi_human_readable,
                nofit, combi, draw_histogram)

            # Recording the result
            result_queue.put(
                (combi_human_readable, myresult)
            )

            # Cleanup
            del list_overlaps_shuffled_for_this_combi
            del true_intersections_for_this_combi

    


    ## Divide the combis into as many threads
    # NOTE : we divide the IDs, not the combis themselves, to save RAM
    combi_id = list(range(len(interesting_combis)))
    random.shuffle(combi_id)    # Shuffle so that combis that are longer to process should be more uniformly spread
    multiproc_batches_of_combis = np.array_split(combi_id, nb_threads)

    combis_done = []
    processes = [None] * nb_threads

    if nb_threads > 1:
      
        ## Submit the batches of computations
        # For each thread...
        for i in range(nb_threads):

            # Get the corresponding combi IDs
            try: b = multiproc_batches_of_combis.pop()
            except: b = []
     
            # If the batch is not empty
            if len(b):

                processes[i] = multiprocessing.Process(
                    target = functools.partial(
                        compute_those_stats, combis = b, result_queue = result_queue
                        )
                )
                processes[i].daemon = True # Prevent zombie processes
                processes[i].start()
                
                message("Submitted a new batch of statistics computations.", type = 'DEBUG')

                

        ## Results collecting: empty the queue whenever possible until all combinations have been processed 
        while len(combis_done) < len(interesting_combis):

            # If the queue is not empty, get all results inside
            while not result_queue.empty():  

                combi_human_readable, result = result_queue.get()
                combis_done += [combi_human_readable]

                # Add the results to the final result dict
                all_results[combi_human_readable] = result

                message("Finished statistics for combi: " + str(combi_human_readable), type='DEBUG')
                message("Combination " + str(len(combis_done)) + "/" + str(len(interesting_combis)) + " done.")

            time.sleep(1) # Don't saturate the CPU by flooding with requests

            #message("Combinations remaining: "+str([c for c in interesting_combis_human_readable if c not in combis_done]), type = 'DEBUG')
            # Careful, `interesting_combis_human_readable` is not exposed in the current version of the code




    # OVERRIDE : if single-threaded, don't use multiprocessing to save RAM
    # NOTE: for now, keeping it for comparison
    else:

        # Makebatches of 2-3 combis instead       
        multiproc_batches_of_combis = np.array_split( range(len(interesting_combis)), int(0.3*len(interesting_combis))+1 )

        while len(combis_done) < len(interesting_combis):
            try: b = multiproc_batches_of_combis.pop()
            except: b = []
            if len(b):
 
                compute_those_stats(combis = b, result_queue = result_queue)


            while (not result_queue.empty()):  
                combi_human_readable, result = result_queue.get()
                combis_done += [combi_human_readable]
                all_results[combi_human_readable] = result
                message("Finished statistics for combi: " + str(combi_human_readable), type='DEBUG')
                message("Combination " + str(len(combis_done)) + "/" + str(len(interesting_combis)) + " done.")




    ## Cleanup
    message("Garbage collection and cleanup...", type = 'DEBUG')
    time.sleep(1)
    # NOTE: less necessary, since this is the last chunk of code executed, but I do it out of precaution nevertheless.

    for p in range(nb_threads): 
        if processes[p] is not None:   # 'None' are potential excedent processes, if nb_threads was too large
            processes[p].join()
            processes[p].close()
    processes.clear()   

    del result_queue
    del mana

    del overlaps_per_combi

    gc.collect()


    # Et voilà !
    return all_results