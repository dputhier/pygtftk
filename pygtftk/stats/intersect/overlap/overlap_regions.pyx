# ## DEBUG NOTE Add the three lines below to be able to run this file independantly (ie. in a Jupyter Kernel)
# import pyximport; pyximport.install()    # Allow to compile Cython code at runtime
# # Load Jupyter Cython extension. Only do it once.
# %load_ext cython
# # Declare that the following "cell" is Cython code
# %%cython

"""
Algorithm to compute the overlap between two sets of genomic regions (as BED files).
Meant to achieve the same functionality as `bedtools intersect`, except with multiple overlaps.
"""

import copy
from multiprocessing import Pool
import gc

import numpy as np
cimport numpy as np
from libc.math cimport isnan
from libc.stdlib cimport malloc, free

import pandas as pd
import pybedtools

from functools import partial


# Import multiprocessing function and structures
cimport pygtftk.stats.multiprocessing.multiproc as multiproc
from pygtftk.stats.multiprocessing.multiproc_structs cimport *




cdef multiproc.FUNC_2DARRAY_RESULT find_overlap(long long* melted, long long* shape, long long number_of_sets) nogil:
    """
    Sweep line algorithm to find overlaps.

    Takes as input a melted array, where one column is a position and the other is a status :
    0 means 'Beginning of a region in A' ; 1 means 'Ending of a region in A'
    2 means 'Beginning of a region in B' ; 3 means 'Ending of a region in B'
    and etc. for C, D, E, ...
    The array MUST be sorted by its first column (position).

    This is an algorithm of the "sweep line" family, which entails using
    "event points" (see code). The problem itself is akin to the "matching
    brackets" problem.

    This function benefits massively from Cython implementation, as it is a
    comparatively simple algorithm where Python overhead is most damning and
    Cython translation is very straightforward. This is a pure Cython version for multiprocessing.

    Remark : as overlaps will only be computed once per shuffle and as we
    compute for whole sets at once, the use of Interval Trees is not justified.


    The function generates a 2D array, with each line having concatenated the start, end
    and then successive flags.


    Author : Quentin FERRE <quentin.q.ferre@gmail.com>
    """

    cdef long i
    cdef int whichset
    cdef long long current_position
    cdef long long status
    cdef bint overlapping
    cdef long long overlap_begin
    cdef long long overlap_end


    cdef int om
    cdef long long s, currently_open


    # Create a result object
    cdef multiproc.FUNC_2DARRAY_RESULT result
    cdef long long xmax = shape[0]
    result.result_shape[0] = xmax
    result.result_shape[1] = 2+number_of_sets # start, end, and the list of flags
    result.result_array = <long long *> malloc(result.result_shape[0]*result.result_shape[1]*sizeof(long long));

    # result.result_array is the overlaps array

    # Allocate result array. Keep track in C of the number of lines, then update
    # the shape accordingly.
    cdef long long current_overlap_number = 0



    # MAJOR TODO: REGIONS ID
    # Use a set, as there will be few elements in them as it is computationally
    # inexpensive to query and change them, to change which regions are currently open.
    # Elements of the set will be (set_id,region_id)
    # cdef set open_regions
    # This is unoptimized currently as it is a "hotfix" to potentially keep region ID, which would be important to in the future do some stats that are specific to region type !



    # OPEN REGIONS
    # Use an numpy array of fixed size, giving [number of regions open in track A, in track B, etc.]
    cdef int* open_regions_flags = <int*> malloc(number_of_sets*sizeof(int))
    for k in xrange(number_of_sets): open_regions_flags[k] = 0
    cdef int* previous_region_flags = <int*> malloc(number_of_sets*sizeof(int))


    # Start at 0.
    current_position = 0
    overlapping = False


    # Continue for each event point : start or end of anyone in A or B
    for i in xrange(shape[0]):

        # We must copy the region flags at the previous position, these are the flags that will get added when we record an overlap
        for k in xrange(number_of_sets): previous_region_flags[k] = open_regions_flags[k]

        # Reminded : arr[i,j] is arr[(i*arr.shape[1])+j]
        current_position = melted[(i*shape[1])+1]
        status = melted[(i*shape[1])+0]

        # "Open" statuses : even numbers
        if status % 2 == 0:
            whichset = <int>(status/2)
            open_regions_flags[whichset] += 1
            # TODO: get the current line number (region_id) and add the `(whichset, region_id) to open_regions

        # "Close" statuses : odd numbers
        if status % 2 != 0:
            whichset = <int>((status-1)/2)
            open_regions_flags[whichset] -= 1
            # TODO: get the current line number (region_id) and REMOVE the correct element FROM open_regions

        # TODO: It would also be very possible to put a score dependant on the region (if we have region_id) instead of '1' in the open_regions_flags, for future improvements.

        # Get number of CURRENTLY open flags

        currently_open = 0
        for s in xrange(number_of_sets): currently_open += open_regions_flags[s]

        # At any time, if at least 2 regions are open : overlap begins
        # & also if we were not overlapping at the previous position
        if (currently_open > 1) & (overlapping == False) :
            overlapping = True
            overlap_begin = current_position

        # As we are browsing by event point, at each iteration something happens
        # (a region opens or closes) so the status of the overlap will change;
        # so we must record the previous overlap
        if (overlapping == True):
            overlap_end = current_position
            #print(current_position, ':', open_regions_flags) # Debug print

            # Only write an overlap if we have moved since setting the overlap == True flag :
            if overlap_begin != overlap_end:

                # NOTE : we must use copy(open_regions_flags) instead of the
                # array itself, otherwise it appends the original full zeroes
                # for some reason.

                ## Write on the correct line
                # Overlap begin
                result.result_array[(current_overlap_number*result.result_shape[1])+0] = overlap_begin
                # Overlap end
                result.result_array[(current_overlap_number*result.result_shape[1])+1] = overlap_end
                # Flags
                for om in xrange(number_of_sets):
                    result.result_array[(current_overlap_number*result.result_shape[1])+2+om] = previous_region_flags[om]

                # Increment overlap number
                current_overlap_number += 1

                # TODO: When we will keep regions_id do something like this but in C : overlaps.append((overlap_begin,overlap_end, open_regions_flags, open_regions))

            # If there is only 1 or 0 regions open, set overlapping to False.
            # Otherwise, it means there are still two regions open and a new "overlap" portion with the new regions currently open
            if currently_open < 2:
              overlapping = False
            else:
              overlap_begin = current_position


    # Finally, truncate the result array to have only the correct number of overla
    result.result_shape[0] = current_overlap_number

    try:
      return result

    finally:
      free(open_regions_flags)
      free(previous_region_flags)

    # This algorithm supports N-fold intersection.
    # The overlaps returned will be consecutive : eg from 15 to 20 there are 3
    # overlapping regions, then from 21 to 43 there are 5, etc.

    # WARNING : this can play havoc with the N statistic.
    # As we note each time how many regions in each set are open, we can record intra-set and inter-set overlap
    # Intra-set overlaps can also be supported, but will likely also play havoc with the N statistic. Just remove all the merge() calls on
    # pybedtools.BedTool objects to be able to use it.





cdef melt_numpy_array(np.ndarray[np.int64_t, ndim = 2] m):
    """
    More efficient function to melt a numpy array.
    Will output an array with two columns giving respectively the value and its
    column of origin.

    WARNING : Values of '-1' are assumed to be NaNs and are removed.

    Author : Quentin FERRE <quentin.q.ferre@gmail.com>
    """

    cdef long imax
    cdef int jmax

    # Get shape of origin array for dynamic memory allocation
    imax = m.shape[0]
    jmax = m.shape[1]

    # Allocate memory in C for the result arrays...
    cdef long long *result_val = <long long *> malloc(imax * jmax * sizeof(long long))
    cdef long long *result_j   = <long long *> malloc(imax * jmax * sizeof(long long))

    # ... and cast to memoryview, refer to the underlying buffer without copy
    # This will be used to cast them back to NumPy arrays efficiently
    cdef long long[:] result_val_view = <long long[:(imax * jmax)]>result_val
    cdef long long[:] result_j_view = <long long[:(imax * jmax)]>result_j

    cdef long i
    cdef int j
    cdef long k

    cdef np.ndarray result

    try:
        # Do columns first, then lines, so the resulting list will be partially
        # sorted - this speeds up the later mergesort.
        k = 0
        for j in range(jmax):
            for i in range(imax):
                val = m[i,j]
                result_val[k] = val
                result_j[k] = j
                k += 1

        # Now for the result : convert them into numpy arrays, concatenate and return
        result = np.vstack((np.asarray(result_j_view), np.asarray(result_val_view)))
        # Remove rows with values of -1 (previously NaN)
        result = result[:,result[1] != -1].transpose()

        return result

    finally:
        # Return the allocated memory
        free(result_val)
        free(result_j)





cpdef find_intersection(tuple beds, all_chrom, return_flags = True, debug = False, nb_threads = 1):
    r"""
    When given a fake bed (as a *list of tuples*, not a text file to avoid text
    overhead) return the intersection computed using our implemented algorithm.

    WARNING : does not work if beds is a list of pybedtools bedfile. It must be a
    list of tuples.

    It will be faster if the BED files were sorted beforehand.

    :Example:

    >>> import numpy as np
    >>> import numpy.testing as npt
    >>> import pyximport; pyximport.install(reload_support=True) # doctest: +ELLIPSIS
    (None, ...
    >>> from pygtftk.stats.intersect.overlap import overlap_regions as oc
    >>> bedA = [('chr1',100,200), ('chr1',300,400), ('chr1',1000,1200)]
    >>> bedB = [('chr1',150,250), ('chr1',350,380)]
    >>> bedC = [('chr1',180,200), ('chr2',100,200)]
    >>> all_chrom = ['chr1']
    >>> expected_result = [('chr1', 150, 180, np.array([1,1,0])), ('chr1', 180, 200, np.array([1,1,1])), ('chr1', 350, 380, np.array([1,1,0]))]
    >>> beds = (bedA, bedB, bedC)
    >>> inter = oc.find_intersection(beds, all_chrom)
    >>> assert all([[npt.assert_array_equal(xa, xb) for xa, xb in zip(a,b)] for a, b in zip(inter, expected_result)])

    """


    # Ensure all_chrom is in a list format, sorted in a predictable order
    all_chrom = sorted(list(all_chrom))

    # Read those fake beds as pandas dataframes
    cdef list dfs = list()
    # This loop below also takes 25% of this step's time !
    for set in range(len(beds)):
        dfs += [pd.DataFrame(beds[set])]


    cdef list overlaps = list()

    cdef list chrlines_per_df = list()

    # This loop below also takes 25% of this step's time !
    for df in dfs:
        chrlines_per_df += [{c:np.where(df[0].values == c)[0] for c in all_chrom}]

    # Iterate over all chromosomes
    cdef list melted_for_all_chroms = list()
    cdef list dfs_for_this_chrom
    cdef np.ndarray mval
    cdef np.ndarray res
    for chrom in all_chrom:

        # Take only the lines whose chromosome match one of those in all_chroms
        dfs_for_this_chrom = list()
        for i in range(len(dfs)):
            lines = chrlines_per_df[i][chrom]
            dfs_for_this_chrom += [dfs[i].iloc[lines,1:3]]


        # Concatenate the dataframe and melt it into a sorted list of
        # event point along with their status (start for A, start for B, ...)
        # As they are concatenated in order, the column 0 is a start for A, 1 is end for A, 2 is start for B, 4 is start for C, etc.
        m = pd.concat(dfs_for_this_chrom, axis=1)
        m.columns = range(len(m.columns))

        # Drop NAs that arise from from not having the same number of lines in A and B
        m.fillna(-1, inplace = True) # Replace NaN with -1 so they are not lost in the typecasting


        # Pass the values to a custom Cython function for fast melting
        mval = m.values.astype(np.longlong) # Sanity check : the positions must be read as np.longlong ; redundant unless strings were accidentally passed (it happened).
        res = melt_numpy_array(mval)
        melted = pd.DataFrame(res, columns = ['variable','value'])


        number_of_sets = int(m.shape[1] / 2) # Two columns per set
        # TODO: keep the line number in the original dataframe to be able to
        # retrace the regions which were shuffled, and pass it as
        # `region_ids` to the cython overlap algorithm

        # We use mergesort because it will be faster since the 4 columns were sorted themselves originally
        melted = melted.sort_values(by='value', kind='mergesort')
        melted.columns = ['status','position']
        #melted['position'] = melted['position'].astype(np.longlong) # TODO is this line really necessary ? I already do that above.
        melted.index = range(len(melted.index))


        # TODO: ADD DETAILS ABOUT THE SORTING. IT IS NOT TRUE THAT THE 4 COLUMNS
        # WERE SORTED (endpoints weren't) BUT, if coming from our shuffles :
        # - All start columns were sorted
        # - All end columns were sorted if we have no overlapping intervals, and
        # nearly sorted (K-sorted) otherwise.
        # These properties mean that mergesort will not have its worst-case 
        # complexity of O(n log n)

        # Add the corresponding numpy dataframe to the list of all melted
        melted_for_all_chroms += [melted.values]



    # Now multiprocess across chromosomes
    """
    NOTE : this part represents only 1% of this step's time. I leave the multiprocessing because it is an interesting example,
    but the rewriting in pure C needed for it to be applied has made multiprocessing counter productive : it's slower than serial execution.
    I may remove it later.
    """
    numbers_of_sets = [number_of_sets] * len(melted_for_all_chroms)

    ov_allchroms = multiproc.apply_func_multiproc_cython(melted_for_all_chroms, numbers_of_sets, find_overlap, nb_threads)

    # And concatenate the results, by chromosome
    chrom_id = 0

    for ov in ov_allchroms:
        # TODO: Keep more information here from the returned overlaps, such as 
        # which regions are involved in each overlap; the code supports it !

        # Classical overlaps only give the coordinates without the flags.
        # If we want the flags, keep them : they are the overlap[2] term
        if not return_flags: 
            overlaps = overlaps + [(all_chrom[chrom_id], overlap[0], overlap[1]) for overlap in ov]
        else: 

            overlaps = overlaps + [
                (
                    all_chrom[chrom_id], overlap[0], overlap[1],
                    np.array(overlap[2:2+number_of_sets], dtype = np.uint32)
                ) for overlap in ov
            ]

        chrom_id = chrom_id + 1


    # Ensure memory is freed
    del ov_allchroms 
    #gc.collect()   # Do not garbage collect here, too slow. Move it to some place else in the future.

    return overlaps




################################################################################
# ---------------------- Run and process intersections ----------------------- #
################################################################################

"""
This section contains only wrappers for other code.
"""

def run_intersect_cython(bed_tuple, all_chroms, nb_threads):
    intersect = find_intersection(bed_tuple, all_chroms, nb_threads = nb_threads)
    return intersect

def compute_intersections_cython(beds_A, beds_B, all_chroms, nb_threads):
    """
    Given two lists of bed files in "list of sets" format (output of batch_to_bedlist)
    and a number of threads, will return their respective intersections (ie.
    intersection of the first of list A with the first of list B (or all of list B is multiple overlaps), second of list
    A with second of list B, etc.)
    """

    # beds_A is a list of beds (has batch_size elements), but beds_B is a list
    # of lists of beds (has nb_sets elements which have batch_size elements)
    # Merge each element of beds_A with the correspondig list of beds in beds_B,
    # since we compute overlaps as a block.

    data = list()
    for i in range(len(beds_A)):

        beds_for_A = [beds_A[i]]

        # For each set in bedsB, get the i-th shuffle and append
        beds_for_B = []
        for j in range(len(beds_B)):
            beds_for_B += [beds_B[j][i]]

        data += [tuple(beds_for_A + beds_for_B)]


    # NOTE for improvement : we could return the intersections separately for each chromosme
    run_intersect_partial = partial(run_intersect_cython, all_chroms = all_chroms, nb_threads = nb_threads)

    all_intersections = [run_intersect_partial(datum) for datum in data]

    return all_intersections





# ---------------------------- Utility functions ----------------------------- #

cpdef bint does_combi_match_query(tuple combi, tuple query, bint exact = False):
    r"""
    Wrapper function to determine if a vector of open flags has the same open
    flags ('contains') as a query.
    If the 'exact' flag is False, it only requires that all flags of combi be open (higher than 0) in query. If True, the same flags must be open and closed.

    Made into a function to ensure consistent behavior across all code by setting the 'exact' flag once.

    >>> from pygtftk.stats.intersect.overlap.overlap_regions import does_combi_match_query
    >>> does_combi_match_query((1,1,0),(1,0,0), exact = False)
    True
    >>> does_combi_match_query((1,1,0),(1,1,1), exact = False)
    False
    >>> does_combi_match_query((1,1,0),(1,0,0), exact = True)
    False
    >>> does_combi_match_query((1,2,1),(1,1,1), exact = True)
    True
    
    """

    cdef int niter = len(combi)

    if exact:
        for i in range(niter):

            if not combi[i]:
                if query[i]: return False
            else:
                if not query[i]: return False

    if not exact:
        for i in range(niter):
            if query[i] and not combi[i]: return False

    return True





## Multiprocessed exactitude computation

# The Python and NumPy types should match, but better safe than sorry
ctypedef fused short_integer_any:
    int
    unsigned int
    np.npy_int32
    np.npy_uint32


cpdef list NPARRAY_which_combis_match_with(short_integer_any[:,:] combis_array, 
                        short_integer_any[:] query,
                        bint exact = False):
    r"""
    This is another version of which_combis_to_get_from, which instead takes
    an array of combinations (combis_arrays) and returns the line numbers
    of the combinations in combis_array that matched the query.

    Example:

    >>> from pygtftk.stats.intersect.overlap.overlap_regions import NPARRAY_which_combis_match_with
    >>> import numpy as np
    >>> query = np.array((1,1,0,0), dtype = np.uint32)
    >>> combis_array = np.array([(1,0,1,1),(1,1,1,0),(1,1,2,1),(1,1,0,0)], dtype = np.uint32)
    >>> indexes = NPARRAY_which_combis_match_with(combis_array, query, exact = False)
    >>> assert indexes == [1,2]
    >>> indexes_2 = NPARRAY_which_combis_match_with(combis_array, query, exact = True)
    >>> assert indexes_2 == []

    """

    cdef Py_ssize_t niter = query.shape[0]
    cdef Py_ssize_t ncombis = combis_array.shape[0]

    cdef int i
    cdef int j

    cdef list results = []

    # This signals when to stop each iteration and move on to the next combi
    cdef bint badsignal = False

    # Do not add the combination if it is completely equal to the query
    # We must register at least one difference
    cdef perfect_equality = True
   
    if exact:
        for j in range(ncombis):
            badsignal = False
            perfect_equality = True

            for i in range(niter):

                #TODO: optimize by integrating in main loop
                if combis_array[j,i] != query[i]: 
                    perfect_equality = False

                    if not combis_array[j,i]:
                        if query[i]: 
                            badsignal = True; break  #return False
                    else:
                        if not query[i]: 
                            badsignal = True; break # return False

            if (not badsignal) and (not perfect_equality):
                results.append(j)

    if not exact:
        for j in range(ncombis):
            badsignal = False
            perfect_equality = True
    
            for i in range(niter):

                if query[i] and not combis_array[j,i]: 
                    badsignal = True; break # return False

                #TODO: optimize by integrating in main loop
                else:
                    if combis_array[j,i] != query[i]: 
                        perfect_equality = False
            
            if (not badsignal) and (not perfect_equality):
                results.append(j)
            
    return results