"""
Algorithm to compute the overlap between two sets of genomic regions (as BED files).
Meant to achieve the same functionality as `bedtools intersect --sorted`.
"""

from multiprocessing import Pool

import numpy as np
cimport numpy as np

import pandas as pd

import pybedtools


cdef find_overlap(np.ndarray melted):
    """
    Sweep line algorithm mixed with the matching brackets problem to find overlaps.

    Takes as input a melted array, where one column is a position and the other is a status :
    0 means 'Beginning of a region in A' ; 1 means 'Ending of a region in A'
    2 means 'Beginning of a region in B' ; 3 means 'Ending of a region in B'
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
    cdef int current_position
    cdef bint openA
    cdef bint openB
    cdef bint overlapping
    cdef list overlaps
    cdef int status
    cdef int overlap_begin
    cdef int overlap_end
    cdef np.ndarray rangei

    # Start at 0.
    current_position = 0
    openA = False
    openB = False
    overlapping = False
    overlaps = list()

    # Continue until we hit an event point : start or end of anyone in A or B
    rangei = np.arange(melted.shape[0])
    for i in rangei:

        current_position = melted[i,1]
        status = melted[i,0]
        # if it's A's start : openA = True
        if status == 0 : openA = True
        # if it's A's end : openA = False
        elif status == 1 : openA = False
        # same for B
        elif status == 2 : openB = True
        elif status == 3 : openB = False

        # At any time, if openA & openB : overlap begins
        if (openA & openB) :
            overlapping = True
            overlap_begin = current_position
        # If overlapping = True & (openA set to false or openB set to false) : stop and record overlap
        if (overlapping == True) & ((openA == False) | (openB == False)):
            overlapping = False
            overlap_end = current_position
            overlaps.append((overlap_begin,overlap_end))

    return overlaps


    ### NOTES for improvement

    # This algorithm could be adapted for N-fold intersection by replacing openA
    # and openB with a vector of open flags (1 flag per set) and replacing conditions :
    #   'openA & openB' becomes 'at least one open'
    #   when any `open` flag changes (any region in any set opens or closes), record overlap

    # If the sets are not merged and can have overlapping regions, the algorithm
    # could use a stack : incrementing a 'currently_open' variable whenever a
    # 'start' position is found and reducing it if a 'stop' position is found.




def find_intersection(fake_bed_A,fake_bed_B, all_chrom):
    """
    When given a fake bed (as a *list of tuples*, not a text file to avoid text
    overhead) return the intersection computed using our implemented algorithm.

    It will be faster if the BED files were sorted beforehand.

    >>> import pyximport; pyximport.install(reload_support=True) # doctest: +ELLIPSIS
    (None, ...
    >>> from pygtftk.stats.intersect.overlap import overlap_regions as oc
    >>> bedA = [('chr1',100,200),('chr1',300,400),('chr1',1000,1200)]
    >>> bedB = [('chr1',150,250),('chr1',350,380)]
    >>> all_chrom = ['chr1']
    >>> expected_result = [('chr1', 150, 200), ('chr1', 350, 380)]
    >>> inter = oc.find_intersection(bedA, bedB, all_chrom)
    >>> assert expected_result == inter

    """
    # Read those fake beds as pandas dataframes
    a_df = pd.DataFrame(fake_bed_A)
    b_df = pd.DataFrame(fake_bed_B)

    overlaps = list()

    chrlines_a = {c:np.where(np.array(a_df[0]) == c)[0] for c in all_chrom}
    chrlines_b = {c:np.where(np.array(b_df[0]) == c)[0] for c in all_chrom}

    for chrom in all_chrom:
        # Take only the lines whose chromosome match one of those in all_chroms
        a = a_df.loc[chrlines_a[chrom],1:3]
        b = b_df.loc[chrlines_b[chrom],1:3]

        # Concatenate the dataframe and melt it into a sorted list of
        # event point along with their status (start for A, start for B, ...)
        m = pd.concat([a,b], axis=1)
        m.columns = range(len(m.columns))
        melted = m.melt().dropna()
        # Drop NAs that arise from from not having the same number of lines in A and B

        # We use mergesort because it will be faster since the 4 columns were sorted themselves
        melted = melted.sort_values(by='value',kind='mergesort')
        melted.columns = ['status','position']
        melted.index = range(len(melted.index))

        ov = find_overlap(melted.values)
        overlaps = overlaps + [(chrom, start, end) for (start, end) in ov]

    return overlaps













################################################################################
# ---------------------- Run and process intersections ----------------------- #
################################################################################

"""
This section contains only wrappers for other code.
"""

from functools import partial

def run_intersect_cython(bed_tuple, all_chroms):
    bedA, bedB = bed_tuple
    intersect = find_intersection(bedA,bedB, all_chroms)
    return intersect

def compute_intersections_cython(beds_A, beds_B,all_chroms, nb_threads = 8):
    """
    Given two lists of bed files in "list of sets" format (output of batch_to_bedlist)
    and a number of threads, will return their respective intersections (ie.
    intersection of the first of list A with the first of list B, second of list
    A with second of list B, etc.)
    """
    data = [(beds_A[i], beds_B[i]) for i in range(len(beds_A))]

    # NOTE for improvement : we could return the intersections separately for each chromosme
    run_intersect_partial = partial(run_intersect_cython, all_chroms = all_chroms)

    with Pool(nb_threads) as p:
        all_intersections = p.map(run_intersect_partial, data)

    return all_intersections
