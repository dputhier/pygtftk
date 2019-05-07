#%%cython
#import pyximport; pyximport.install()
#%load_ext cython

"""
Algorithm to compute the overlap between two sets of genomic regions (as BED files).
Meant to achieve the same functionality as `bedtools intersect --sorted`.
Rudimentary support for multiple overlaps is present.
"""

from multiprocessing import Pool

import numpy as np
cimport numpy as np

import pandas as pd

import pybedtools


cdef find_overlap(np.ndarray melted, int number_of_sets, np.ndarray region_ids = None):
    """
    Sweep line algorithm to find overlaps.

    Takes as input a melted array, where one column is a position and the other is a status :
    0 means 'Beginning of a region in A' ; 1 means 'Ending of a region in A'
    2 means 'Beginning of a region in B' ; 3 means 'Ending of a region in B'
    and etc. for C, D, E, ...
    The array MUST be sorted by its first column (position).

    This is an algorithm of the "sweep line" family, which entails using
    "event points" (see code). The problem itself is akin to the "matching brackets"
    problem.

    This function benefits massively from Cython implementation, as it is a
    comparatively simple algorithm where Python overhead is most damning and
    Cython translation is very straightforward.

    Remark : as overlaps will only be computed once per shuffle, the use of
    Interval Trees (whose cost would be offset by faster interval query) is
    not justified.
    """
    cdef long current_position
    cdef bint overlapping
    cdef list overlaps
    cdef int status
    cdef long overlap_begin
    cdef long overlap_end
    cdef np.ndarray rangei
    cdef np.ndarray open_regions_flags
    cdef int whichset
    cdef list all_flags = list()


    # OPEN REGIONS
    # Use an numpy array of fixed size, giving [number of regions open in track A, in track B, etc.]
    open_regions_flags = np.array([0] * number_of_sets)


    # REGIONS ID
    # Use a set, as there will be few elements in them as it is computationally
    # inexpensive to query and change them, to change which regions are currently open.
    # Elements of the set will be (set_id,region_id)
    # cdef set open_regions
    # This is unoptimized currently as it is a "hotfix" to potentially keep region ID.


    # Start at 0.
    current_position = 0
    overlapping = False
    overlaps = list()

    # Continue for each event point : start or end of anyone in A or B
    rangei = np.arange(melted.shape[0])
    for i in rangei:

        current_position = melted[i,1]
        status = melted[i,0]

        # "Open" statuses : even numbers
        if status % 2 == 0:
            whichset = int(status/2)
            open_regions_flags[whichset] += 1
            # TODO get the current line number (region_id) and add the `(whichset, region_id) to open_regions

        # "Close" statuses : odd numbers
        if status % 2 != 0:
            whichset = int((status-1)/2)
            open_regions_flags[whichset] -= 1
            # TODO get the current line number (region_id) and REMOVE the correct element FROM open_regions


        # At any time, if at least 2 regions are open : overlap begins
        # & also if we were not overlapping at the previous position
        if (np.sum(open_regions_flags) > 1) & (overlapping == False) :
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
                overlaps.append((overlap_begin,overlap_end, open_regions_flags))
                # TODO When we will keep regions_id do this : overlaps.append((overlap_begin,overlap_end, open_regions_flags, open_regions))


            # If there is only 1 or 0 regions open, set overlapping to False.
            # Otherwise, it means there are still two regions open and a new "overlap" portion with the new regions currently open
            if np.sum(open_regions_flags) < 2:
              overlapping = False
            else:
              overlap_begin = current_position

    return overlaps

    # This algorithm supports N-fold intersection.
    # The overlaps returned will be consecutive : eg from 15 to 20 there are 3 overlapping regions, then from 21 to 43 there are 5, etc.
    # WARNING : this can play havoc with the N statistic. Statistic theory for such a thing is not developed yet.
    # As we note each time how many regions in each set are open, we can record intra-set and iter-set overlap







def find_intersection(beds, all_chrom):
    """
    When given a fake bed (as a *list of tuples*, not a text file to avoid text
    overhead) return the intersection computed using our implemented algorithm.

    It will be faster if the BED files were sorted beforehand.

    >>> import pyximport; pyximport.install(reload_support=True) # doctest: +ELLIPSIS
    (None, ...
    >>> from pygtftk.stats.intersect.overlap import overlap_regions as oc
    >>> bedA = [('chr1',100,200),('chr1',300,400),('chr1',1000,1200)]
    >>> bedB = [('chr1',150,250),('chr1',350,380)]
    >>> bedC = [('chr1',180,200)]
    >>> all_chrom = ['chr1']
    >>> expected_result = [('chr1', 150, 180), ('chr1', 180, 200), ('chr1', 350, 380)]
    >>> beds = (bedA, bedB, bedC)
    >>> inter = oc.find_intersection(beds, all_chrom)
    >>> assert expected_result == inter

    """


    # Read those fake beds as pandas dataframes
    dfs = list()
    for set in range(len(beds)):
        dfs += [pd.DataFrame(beds[set])]

    overlaps = list()

    chrlines_per_df = list()
    for df in dfs:
        chrlines_per_df += [{c:np.where(np.array(df[0]) == c)[0] for c in all_chrom}]



    for chrom in all_chrom:
        # Take only the lines whose chromosome match one of those in all_chroms
        dfs_for_this_chrom = list()
        for i in range(len(dfs)):
            current_df = dfs[i]
            lines = chrlines_per_df[i][chrom]
            dfs_for_this_chrom += [current_df.loc[lines,1:3]]



        # Concatenate the dataframe and melt it into a sorted list of
        # event point along with their status (start for A, start for B, ...)
        # As they are concatenated in order, the column 0 is a start for A, 1 is end for A, 2 is start for B, 4 is start for C, etc.
        m = pd.concat(dfs_for_this_chrom, axis=1)
        m.columns = range(len(m.columns))
        melted = m.melt().dropna()
        # Drop NAs that arise from from not having the same number of lines in A and B


        number_of_sets = int(m.shape[1] / 2) # Two columns per set
        # TODO keep the line number in the original dataframe to be able to
        # retrace the regions which were shuffled, and pass it as
        # `region_ids` to the cython overlap algorithm

        # We use mergesort because it will be faster since the 4 columns were sorted themselves
        melted = melted.sort_values(by='value',kind='mergesort')
        melted.columns = ['status','position']
        melted['position'] = melted['position'].astype(np.long)
        melted.index = range(len(melted.index))

        ov = find_overlap(melted.values, number_of_sets)
        overlaps = overlaps + [(chrom, overlap[0], overlap[1]) for overlap in ov]
        # TODO Keep more information here from the returned overlaps, such as which regions are involved in each overlap; the code supports it

    return overlaps







################################################################################
# ---------------------- Run and process intersections ----------------------- #
################################################################################

"""
This section contains only wrappers for other code.
"""

from functools import partial

def run_intersect_cython(bed_tuple, all_chroms):
    intersect = find_intersection(bed_tuple, all_chroms)
    return intersect

def compute_intersections_cython(beds_A, beds_B, all_chroms, nb_threads = 8):
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
    run_intersect_partial = partial(run_intersect_cython, all_chroms = all_chroms)

    with Pool(nb_threads) as p:
        all_intersections = p.map(run_intersect_partial, data)

    return all_intersections
