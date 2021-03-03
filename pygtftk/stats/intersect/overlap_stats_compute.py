"""
Compute overlap statistics on shuffled sets.

Called by overap_stats_shuffling.compute_overlap_stats(), hence the name.

Those functions tend to take as input lists of shuffles and output statistics
"""

import functools
import gc
import multiprocessing
import time
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
import bisect
import copy

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import nbinom

from pygtftk.stats import negbin_fit as nf
from pygtftk.stats.intersect.modl import dict_learning as dl
from pygtftk.stats.intersect.overlap import overlap_regions as oc
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


################################################################################
# ----------------------- Process sets of intersections ---------------------- #
################################################################################

def stats_single(all_intersections_for_this_combi, true_intersection,
                 ft_type='some feature', nofit=False, this_combi_only=None):
    """
    Compute statistics such as total number of overlapping base pairs for a given feature.

    :param intersections_for_this_combi: an object returned by compute_all_intersections_minibatch giving all intersections between the shuffled bedA and bedsB.
    :param true_intersection: the result of overlap_stats_compute.compute_true_intersection(bedA, bedsB) where bedA is the query and bedB is the list of the bed files between whose shuffles the aforementioned intersections have been computed. Used here to calculate the true intersections between them and calculate a Neg Binom p-value.
    :param ft_type: for debug messages, which feature/combi are we currently processing ?
    :param nofit: if True, does not do Negative Binomial fitting
    :param this_combi_only: a list of flags (e.g. [1,0,0,1]) corresponding to expected flags in the interescetions, one per file (see find_intersection() source and documentation). If not None, we will consider only intersections that have this flag for the number of true intersections and true overlapping basepairs
    """

    message('Processing overlaps for ' + ft_type, type='DEBUG')

    start = time.time()
    stats = [compute_stats_for_intersection(myintersect) for myintersect in all_intersections_for_this_combi]

    stop = time.time()
    message(ft_type + '- Statistics on overlaps computed in : ' + str(stop - start) + ' s.', type='DEBUG')

    # Unpack the stats.
    bp_overlaps = [s[0] for s in stats]  # Those are the individual overlap lengths, a list of lists
    summed_bp_overlaps = [sum(x) for x in bp_overlaps]  # Sum by shuffle
    intersect_nbs = [s[1] for s in stats]

    # NOTE FOR IMPROVEMENT : it would be interesting to return the average size
    # of an overlap as well, per shuffle. Since our intersection algorithm returns
    # details about the intersections like `bedtools intersect` would, this could
    # be computed without much hassle.

    # ------ Fitting of a Negative Binomial distribution on the shuffles ----- #
    start = time.time()

    # Fitting can be disabled from the main function (for now, mainly relevant if we used a Markov model instead of a classical one.)
    if nofit:
        ps = pn = -1


    else:
        # Renaming expectations and variances
        expectation_fitted_summed_bp_overlaps, variance_fitted_summed_bp_overlaps = np.mean(summed_bp_overlaps), np.var(
            summed_bp_overlaps)
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
    true_intersect_nb = len(true_intersection)
    true_bp_overlaps = sum([x[2] - x[1] for x in true_intersection])

    # Compute the p-values using the distribution fitted on the shuffles.
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
    message(ft_type + '- Negative Binomial distributions fitted in : ' + str(stop - start) + ' s.',
            type='DEBUG')

    # --------------------------------------------------------------------------
    # Draft code for diagnostic plots of the distribution of each statistic in
    # the shuffles. Kept for potential future improvement.
    # --------------------------------------------------------------------------

    # TODO: Drawing is computationally expensive, make it optional


    # This may return an error for very long ft_type strings, so wrap it in a try-except block
    # TODO: Fix this properly

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
        BINS = 100
        de = plt.hist(summed_bp_overlaps, bins=BINS)[1]
        # Plot the PDF.
        xmin, xmax = min(de), max(de)
        x = np.linspace(xmin, xmax, BINS)
        try:
            d = [nbinom.cdf(x[i], r, p) - nbinom.cdf(x[i - 1], r, p) for i in range(1, len(x))]
            d = np.array([0] + d) * len(summed_bp_overlaps)
        except:
            d = [0] * BINS

        plt.plot(x, d, 'k', linewidth=2)
        plt.savefig(hist_S.name)
        plt.close()
    except:
        pass

    """
    TODO: currently, those files remain in /tmp because this function is subprocessed.
    I must prepare a temp file manager like make_tmp_file_pool() to keep them in the directory specified by -K
    """

    # Now return the result as a dictionary of statistics for this calculation.

    # WARNING Be careful to use the same order as result_abort, in overlap_stats_shuffling.py !
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

    #message(ft_type + '- Result dump : ' + str(result), type='DEBUG')

    return result


# -------------------------- Multiple overlap sets --------------------------- #


## ------------ Helper objects

class ComputingStatsCombiPartial(object):
    """
    This is a wrapper to compute statistics for one combination
    """

    # Remember the parameters
    def __init__(self,
                 all_intersections_for_this_combi, ft_type, true_intersection, this_combi_only, nofit,
                 #result_queue,  # The result queue
                 combi_human_readable):
        # Parameters for compute_stats_single
        self.all_intersections_for_this_combi = all_intersections_for_this_combi
        self.ft_type = ft_type
        self.true_intersection = true_intersection
        self.this_combi_only = this_combi_only
        self.nofit = nofit

        #self.result_queue = result_queue  # Result queue

        self.combi_human_readable = combi_human_readable

    # Callable
    def __call__(self):

        try:
            my_result = stats_single(self.all_intersections_for_this_combi,
                                    self.true_intersection, self.ft_type, 
                                    self.nofit, self.this_combi_only)
        except Exception as e:
            message("Exception during stats_single call"+e, type = "DEBUG")

        # Put as tuple so we may extract it later and put it in dict, we know which combi was processed
        #self.result_queue.put((self.combi_human_readable, my_result))
        #del my_result
        return (self.combi_human_readable, my_result)



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

        c = tuple(c) # Force conversion to tuple

        if oc.does_combi_match_query(c, combi, exact=exact):
            if c != combi: 
                matching_vector += [i]
        i = i+1

    # The matching vector is a Python list, convert it to numpy array to save some RAM
    matching_vector = np.asarray(matching_vector, dtype= np.uint64)

    message("Computed exactitude for combi "+str(combi)+"...", type="DEBUG")

    return combi, matching_vector


def index_all_these(combis_to_index, all_combis, exact, my_result_queue):
    """
    Helper function to run which_combis_to_get_from on a list of combis and put the results in a queue
    """
    for combi in combis_to_index:
        res = which_combis_to_get_from(combi = combi, 
                all_possible_combis = all_combis, exact = exact)

        my_result_queue.put(res)



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
    >>> d = {'A':[[1,1],[1],[1,1]],'B':[[2],[2,2],[2]],'C':[[3,3],[3],[3]]}
    >>> di = DictionaryWithIndex(d, index)
    >>> assert di.get_all('A') == [[1,1,2,3,3], [1,2,2,3], [1,1,2,3]]
    >>> assert di.get_all('B') == [[2],[2,2],[2]]

    """

    def __init__(self, data, index, data_default_factory=lambda: []):
        """
        :param data: The dictionary to be wrapped
        :param index: a CombinationExactMapping object giving for each key in data the *other* keys to be also used when calling each key
        :param data_default_factory: A function producing the default value of data after wrapping. Defaults to `lambda:[]`
        """
        # Index and data should both be defaultdict to handle missing cases
        self.data = defaultdict(data_default_factory, data)
        self.index = index

        gc.collect()  # Force garbage collecting in case of large dictionaries

    def get_all(self, key):
        # Retrieve all combinations from the index
        all_keys_to_get = self.index.get_all(key)

        # Start with the value for the key itself
        # Important to use copy here to not modify the original
        val_concat = copy.deepcopy(self.data[key])

        # For each other key to be added...
        for k in all_keys_to_get:

            # Merge all elements of the supplementary key to the original key
            for i in range(len(self.data[key])):
                val_concat[i] += self.data[k][i]

        return val_concat

    def get_simple_concatenation(self, key):
        """
        For the true overlaps. Do not merge individual shuffles, there are no shuffles to be merged.
        Simply concatenate the lists.
        """
        all_keys_to_get = self.index.get_all(key)
        val_concat = copy.deepcopy(self.data[key])
        for k in all_keys_to_get: val_concat += self.data[k]
        return val_concat




def do_all_calls(my_calls, my_result_queue):
    """
    Helper function that simply runs all functools.partial calls in my_calls and 
    puts the results in my_result_queue.
    """

    for call in my_calls:

        try:
            my_result = call()
            my_result_queue.put(my_result)
        except Exception as e:
            message("Exception raised by do_all_calls(): "+str(e), type = "DEBUG")



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
        return max(self.current_index - 1, 0)

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



## ------ Main function

def stats_multiple_overlap(all_overlaps, bedA, bedsB, all_feature_labels, nb_threads=8, nofit=False,
                           # Parameters for the finding of interesting combis
                           multiple_overlap_target_combi_size=None,
                           multiple_overlap_max_number_of_combinations=None,
                           multiple_overlap_custom_combis=None,
                           ):
    """
    Instead of returning one set of overlap stats per query type (ie. exons, gens, bedfile1, bedfile2, etc...)
    it will return one per multiple overlap (ie. here there was a peak for bedfile1+gene, or bedfile1+bedfile2+gene, etc.)

    Only do so for "common/interesting combis", found by dict learning on original data.

    :param all_overlaps: The list of all overlaps computed for all shuffles
    """

    ## ------  Parameters for the finding of interesting combis

    # If the combi size (`multiple_overlap_target_combi_size`) and the max number of combinations (`multiple_overlap_max_number_of_combinations`) were not set manually, they default to -1 meaning no restrictions.

    # Exactitude : should an intersection of A+B+C count towards looking for A+B ?
    # By default, yes, meaning we look for "inexact" combis.
    # We will look instead for exclusive combis if the user manually specifies 
    # multiple_overlap_target_combi_size equal to the number of sets

    # Rk number of sets is len(bedsB) +1, let's not forget query !

    if multiple_overlap_target_combi_size == (len(bedsB) + 1):
        exact = True
    else:
        exact = False

    # ------------- Override combinations
    # If custom_combis are set, skip all the combination mining above.
    # Directly read a text file.
    if multiple_overlap_custom_combis is not None:
        message('Working on custom combinations for multiple overlap.')

        # Read NumPy matrix and cast to regular Python list
        interesting_combis_matrix = np.loadtxt(multiple_overlap_custom_combis.name, dtype=int)
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
        raise ValueError(
            "No combination of interest found. Try disabling --multiple-overlap-max-number-of-combinations or increasing the number of shuffles.")

    ## Precompute the true intesections between bedA and bedsB once and for all
    message("Computing all true intersections...")
    true_intersection = compute_true_intersection(bedA, bedsB)

    # stats_single() requires a list of all true intersections be passed for each combi.
    # Split the true intersections and prepare those

    # Read intersection flag as tuples for the interesting combis only

    true_intersections_per_combi = defaultdict(list)
    for inter in true_intersection:
        combi = tuple(inter[3])
        true_intersections_per_combi[combi] += [inter]






    # ------------- Splitting all the intersections for all shuffles
    # Split the list of overlap regions per combis

    import time  # Needed to re-import here for some reason ?

    message("Pausing for 5 seconds to let RAM garbage collection run.")
    time.sleep(5)
    gc.collect()

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

        for overlap in intersections_for_this_shuffle:
            # What are the flags (ie. sets present) for this overlap ? For example, "A+C but not B" is (1,0,1)
            flags_for_this_overlap = tuple(overlap[3])  # Convert to a tuple because lists cannot be dict keys

            ## Add exact matches
            # Don't need to remember the flags since I group them by intersection 
            filtered_minibatches[flags_for_this_overlap] += [overlap[0:3]]

        ## Now add the filtered minibatches as elements of a list
        # overlaps_per_combi then contains lists of lists : one list per original
        # minibatch which contains a list of all overlaps WITH THIS FLAG for this minibatch

        # For the combinations already encountered and also encountered here, add the filtered batch
        for combi, filtered_batch in filtered_minibatches.items():
            overlaps_per_combi[combi].put(filtered_batch)

        # For the combinations already encountered but NOT encountered here, add an empty batch
        oc_keys = set(overlaps_per_combi.keys())
        fm_keys = set(filtered_minibatches.keys())
        combis_already_encountered_but_not_here = [
            c for c in oc_keys if c not in fm_keys
        ]
        for c in combis_already_encountered_but_not_here:
            overlaps_per_combi[c].put([])


        # For the combinations not yet encountered, update the default factory to include an additional empty batch     
        starting_index += 1
        new_factory_expression = 'lambda: SparseListOfLists('+str(starting_index)+')'
        overlaps_per_combi.default_factory = eval(new_factory_expression)

        message("Splitting done for "+str(starting_index)+'/'+str(total_to_split), type = "DEBUG")


    ## Partial matches
    # We have registered all exact combis, now add partial matches. Partial matches are defined as "including all flags",
    # meaning (1,1,1,0) will be a match if we query (1,1,0,0) since it contains all its flags, but will not be a match for (1,0,0,1)
    # Iterate over all combinations, and if they match add contents of combi B to combi A

    ## Compute the index of combis to be fetched
    # Relevant for partial matches.

    message("Pausing for 5 seconds to let RAM garbage collection run.")
    time.sleep(5)
    gc.collect()

    message("Computing index of exact/inexact combinations...")

    # Compute the index for all combis found, but also for the interesting combis and all true combis. Relevant mostly for inexact combis.
    all_combis = list(
        set(list(overlaps_per_combi.keys()) + interesting_combis + list(true_intersections_per_combi.keys()))
    )

    message("We will index " + str(len(interesting_combis))+"*"+str(len(all_combis)) + " combinations. This can be very long (minutes) for longer combinations.")


    # NOTE We do not need to index all combis : the only ones that will ever be queried are the `interesting_combis`. Those need to 
    # be fully indexed against `all_combis.
    # However but not every combination in `all_combis` : we do not care what we would need to get if we were to query a combination C
    # that is in `all_combis`, but not in `interesting_combis`

    
    # Multiprocessing objects
    mana = multiprocessing.Manager()
    result_queue = mana.Queue()
    pool = ProcessPoolExecutor(nb_threads)

    # Divide the interesting combis into as amny batches as threads, and remove empty batches
    multiproc_batches_of_combis_exactitude = np.array_split(np.array(interesting_combis), nb_threads)
    multiproc_batches_of_combis_exactitude = [batch_array for batch_array in multiproc_batches_of_combis_exactitude if (batch_array.ndim and batch_array.size)]
  

    mappings = [] # Final results list, to be filled when emptying the queue

    jobs_in_progress = 0

    # Submit to, and empty the queue whenever possible until all combinations have been processed   
    while len(mappings) < len(interesting_combis):

        if jobs_in_progress < nb_threads:

            # Try to get a new batch
            try: current_batch = multiproc_batches_of_combis_exactitude.pop()
            except: current_batch = []

            if current_batch != []:

                pool.submit(index_all_these, combis_to_index = current_batch,
                    all_combis = all_combis, exact = exact,
                    my_result_queue = result_queue) 

                message("Submitted a batch of exactitude computations to the queue...")
                jobs_in_progress += 1

            time.sleep(0.01)

        # Monitor the queue and empty it whenever possible
        if (not result_queue.empty()):  

            result = result_queue.get()
            mappings += [result]

            jobs_in_progress -= 1
            time.sleep(0.01)

        else:
            message("Waiting for indexing of exact/inexact combinations...", type = 'DEBUG')
            time.sleep(0.5)

    pool.shutdown()
    

       
        
    # Convert the mappings into a sparser object for storage    
    final_mapping = CombinationExactMapping(all_combis, dict(mappings))
    message("Index computed. Repartition of overlaps...")



    ## Finally, create a DictionayWithIndex object to hold all intersections
    # If exact = False, when querying this dictionary using get_all(c), all 
    # combis that are an inexact match for c will also count.
    overlaps_per_combi = DictionaryWithIndex(overlaps_per_combi, final_mapping,
                                             data_default_factory=overlaps_per_combi.default_factory)
    # Overwrite original object to save memory !

    ## Do the same for the true intersections
    true_intersections_per_combi = DictionaryWithIndex(true_intersections_per_combi, final_mapping,
                                                       data_default_factory=true_intersections_per_combi.default_factory)
    # We will pass those true_intersections to stats_single and resume as usual

    stop = time.time()
    message('All computed overlaps for the shuffles split by combination in : ' + str(stop - start) + ' s among ' + str(
        all_feature_labels) + '.', type='DEBUG')

    # ------------- Enrichment for each combination
    # Now call stats_single on each.

    ## Result queue
    mana = multiprocessing.Manager()
    result_queue = mana.Queue()
    all_results = dict()  # Final result dict, to be filled when emptying the queue
    pool = ProcessPoolExecutor(nb_threads)  # Process pool

    ## ---- Now prepare the jobs for submission

    def combis_to_partials(combis):
        """
        This is NOT a pure function and manipulates external objects. It's simply a helper to make the code more legible
        """

        jobs = list()

        for combi in combis:

            # Convert the combi into a user-friendly string
            indices = [i for i in range(1, len(combi)) if combi[i] != 0]  # 0-based !
            combi_list = [all_feature_labels[i - 1] for i in indices]

            if combi[0] != 0: combi_list = ['Query'] + combi_list

            # Add '...' to the combi if not exact, to show the user that others TFs
            # could also be present in these intersections
            if exact: combi_human_readable = '[' + ' + '.join(combi_list) + ']'
            if not exact: combi_human_readable = '[' + ' + '.join(combi_list) + ' + ... ]'

            message('Will compute statistics for the combination : ' + str(combi_human_readable), type='DEBUG')

            combi_key = tuple(combi)  # Convert to a tuple because lists cannot be dict keys

            # Collect all shuffles for this combination, taking exactitude into account
            list_overlaps_shuffled_for_this_combi = overlaps_per_combi.get_all(combi_key)

            # Create a sort-of partial call
            compute_stats_combi_partial = ComputingStatsCombiPartial(list_overlaps_shuffled_for_this_combi,
                                                                    combi_human_readable,
                                                                    true_intersections_per_combi.get_simple_concatenation(combi_key),
                                                                    combi, nofit,
                                                                    combi_human_readable)

            # Add to the jobs to be executed, that will be returned
            jobs += [compute_stats_combi_partial]
        
        return jobs


    # Split the interesting_combis into batches of combis, one per process 
    # (maybe 20 times that just to be safe and not send too many combis at once
    # to the pool, with their accompanying all_overlaps)
    multiproc_batches_of_combis = np.array_split(np.array(interesting_combis), 20*nb_threads)

    # Remove empty batches
    multiproc_batches_of_combis = [batch_array for batch_array in multiproc_batches_of_combis if (batch_array.ndim and batch_array.size)]


    ## Submit to, and empty the queue whenever possible until all combinations have been processed   
    combis_done = []
    jobs_in_progress = 0

    if nb_threads > 1:

        while len(combis_done) < len(interesting_combis):

            # One cannot directly monitor this through `pool` object, so instead I set an external `jobs_in_progress` variable
            # If there are less pools being used than there are available, we can submit a new batch of jobs
            if jobs_in_progress < nb_threads:

                # Only if there are combinations left to be processed
                try: current_batch = multiproc_batches_of_combis.pop()
                except: current_batch = []

                # If the batch to be submitted is not empty...
                # Now submit an entire batch of combis : prepare a list of partials and submit it to the queue
                if current_batch != []:

                    #message("Will send this batch of combinations: "+str(current_batch), type = "DEBUG")
                    current_calls = combis_to_partials(current_batch)

                    try:
                        pool.submit(do_all_calls, my_calls = current_calls, my_result_queue = result_queue) 
                    pool.submit(do_all_calls, my_calls = current_calls, my_result_queue = result_queue) 
                        pool.submit(do_all_calls, my_calls = current_calls, my_result_queue = result_queue) 
                        message("Submitted a new batch of statistics computations.", type = 'DEBUG')
                        # Rk : the submit() function returns a Future object. It is not kept here, but could be.
                    except Exception as e:
                        message("Exception when submitting a batch of stats computations: "+str(e), type = 'DEBUG')

                    jobs_in_progress += 1 # A job was submitted, increment jobs_in_progress

                time.sleep(0.01)


            # If the queue is empty, wait a bit and try again, to not saturate the CPU with requests
            if (not result_queue.empty()):  
        if (not result_queue.empty()):  
            if (not result_queue.empty()):  

                combi_human_readable, result = result_queue.get()
                combis_done += [combi_human_readable]

                # Add the results to the final result dict
                all_results[combi_human_readable] = result

                message("Finished statistics for combi: " + str(combi_human_readable), type='DEBUG')
                message("Combination " + str(len(combis_done)) + "/" + str(len(interesting_combis)) + "done.")

                jobs_in_progress -= 1 # A job was completed, decrement jobs_in_progress

                time.sleep(0.01)

            else:
                time.sleep(1)
                message("Waiting for results in the result queue...", type = 'DEBUG')
                #message("Combinations remaining: "+str([c for c in interesting_combis_human_readable if c not in combis_done]), type = 'DEBUG')
                # Careful, `interesting_combis_human_readable` is not exposed in the current version of the code


    # OVERRIDE : if single-threaded, don't use multiprocessing to save RAM
    else:

        # Make batches of 2-3 combis instead
        multiproc_batches_of_combis = np.array_split(np.array(interesting_combis), int(0.3*len(interesting_combis)))

        while len(combis_done) < len(interesting_combis):

            # Only if there are combinations left to be processed
            try: current_batch = multiproc_batches_of_combis.pop()
            except: current_batch = []

            if current_batch != []:
                do_all_calls(my_calls = current_calls, my_result_queue = result_queue) 

            # If the queue is empty, wait a bit and try again, to not saturate the CPU with requests
            while (not result_queue.empty()):  

                combi_human_readable, result = result_queue.get()
                combis_done += [combi_human_readable]

                # Add the results to the final result dict
                all_results[combi_human_readable] = result

                message("Finished statistics for combi: " + str(combi_human_readable), type='DEBUG')
                message("Combination " + str(len(combis_done)) + "/" + str(len(interesting_combis)) + "done.")

    # Cleanup
    del result_queue

    pool.shutdown()
    del pool
    del mana
    time.sleep(1)
    message("Pause for 1 second to let garbage collection run...", type = 'DEBUG')
    gc.collect()

    return all_results
