# distutils: language = c++
# distutils: sources = exclude.cpp


"""
A set of functions to turn a BED file into a list of intervals, and exclude
certain regions to create concatenated sub-chromosomes.
"""

from multiprocessing import Pool
from functools import partial
from collections import Counter

import pybedtools
import cython
import numpy as np
cimport numpy as np
import pandas as pd

from pygtftk.utils import message




################################################################################
# -------------------------- Reading bed files ------------------------------- #
################################################################################

def bed_to_lists_of_intervals(bed, chromsizes):
    r"""
    Reads a bed file (as a pybedtools.BedTool object) and returns respectively
    two dictionaries, with the list of region lengths and interregions lengths
    (resp. Lr and Li), as well as a list of all chromosomes.

    The dictionaries group the lists of lengths by chromosome.

    A dictionary of chromosome sizes must also be provided to compute the
    distance between the last feature and the chromosome end.

    For example, Li['chr1'] is the list of distances between regions in chr1.

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.stats.intersect.read_bed.read_bed_as_list import bed_to_lists_of_intervals
    >>> import pybedtools
    >>> import numpy as np
    >>> import warnings
    >>> warnings.simplefilter(action='ignore', category=FutureWarning)
    >>> import numpy.testing as npt
    >>> f = pybedtools.BedTool(get_example_file("simple","bed")[0])
    >>> c = get_example_file(ext="chromInfo")[0]
    >>> from pygtftk.utils import chrom_info_as_dict
    >>> cl = chrom_info_as_dict(open(c, "r"))
    >>> result = bed_to_lists_of_intervals(f,cl)
    >>> npt.assert_array_equal(result[0]['chr1'], np.array([ 5, 10]))
    >>> npt.assert_array_equal(result[1]['chr1'], np.array([ 10, 25, 250]))
    >>> npt.assert_array_equal(result[2], np.array(['chr1']))

    """

    # Convert bedfile to pandas dataframe
    bed = bed.to_dataframe()

    Lr = dict()
    Li = dict()
    
    bed.chrom = bed.chrom.apply(str)
    all_chrom = np.unique(bed.chrom)

    for chrom in all_chrom:

        # Select only the features on this chromosome
        features_on_this_chrom = bed[bed.chrom == chrom].index

        lr = list()
        li = list()

        previous_feature_stop = 0 # This way, the distance between chromosome beginning and first feature is covered
        for f in features_on_this_chrom:
            lr.append(bed.at[f, 'end'] - bed.at[f, 'start'])
            li.append(bed.at[f, 'start'] - previous_feature_stop)
            previous_feature_stop = bed.at[f, 'end']

        # Add the missing inter-region distance between the last feature and the chromosome end
        last_li = chromsizes[str(chrom)] - previous_feature_stop # Query str(chrom) because chrom may be numeric, but the keys of chromsizes passed by pytftk are always strings.
        if last_li < 0:
            last_li = 0
            message('Warning - You have a bed file with features after the end of chromosome "'+str(chrom)+'" !', type='INFO')
        li.append(last_li)

        # Remember to convert chrom to a string (it may have been a numeric, but the rest of the program requires it to be a string)
        Lr[str(chrom)] = np.array(lr)
        Li[str(chrom)] = np.array(li)

    return Lr, Li, all_chrom






################################################################################
# ----------------------- Exclusion and concatenation ------------------------ #
#################################################################################

def exclude_chromsizes(exclusion, chromsizes):
    r"""
    Shortens the chromsome sizes (given as a dictionary) by the total length of
    each excluded region (given as a BedTool file).

    >>> from pygtftk.utils import get_example_file
    >>> from collections import OrderedDict
    >>> from pygtftk.stats.intersect.read_bed.read_bed_as_list import exclude_chromsizes
    >>> import pybedtools
    >>> import numpy.testing as npt
    >>> c = get_example_file(ext="chromInfo")[0]
    >>> from pygtftk.utils import chrom_info_as_dict
    >>> cl = chrom_info_as_dict(open(c, "r"))
    >>> e_string = 'chr1\t0\t100\nchr2\t0\t300'
    >>> e = pybedtools.BedTool(e_string,from_string=True).sort().merge()
    >>> result = exclude_chromsizes(e,cl)
    >>> assert result == OrderedDict([('chr1', 200), ('chr2', 300), ('all_chrom', 900)])

    """

    exclusion = exclusion.to_dataframe()
    for _, excl in exclusion.iterrows():
        excl_length = abs(excl['end'] - excl['start'])
        chromsizes[excl['chrom']] = chromsizes[excl['chrom']] - excl_length
    return chromsizes




# ----------------------- Exclusion for each bed file ------------------------ #
# The function to perform it is C++ code, which we interface through Cython.

# Declare the interface to C++ code
cdef extern from "exclude.hpp" namespace "exclusion":
  void cpp_excludeConcatenateForThisChrom(long long* bedfile_starts, long long* bedfile_ends,
                                        long long* exclusion_starts, long long* exclusion_ends,
                                        long long* result_starts, long long* result_ends,
                                        long bed_size, long excl_size)


cdef cpp_wrapper_for_exclusion(np.ndarray[longlong, ndim=1, mode="c"] bedfile_start_nparray,
                np.ndarray[longlong, ndim=1, mode="c"] bedfile_end_nparray,
                np.ndarray[longlong, ndim=1, mode="c"] exclusion_start_nparray,
                np.ndarray[longlong, ndim=1, mode="c"] exclusion_end_nparray,
                np.ndarray[longlong, ndim=1, mode="c"] result_start_nparray,
                np.ndarray[longlong, ndim=1, mode="c"] result_end_nparray,
                long bed_size,
                long excl_size):
    """
    Create a cdef-defined  wrapper to call the C++ code. Needed because if we call
    cpp_excludeConcatenateForThisChrom() directly in a Python function (defined
    with `def`) it will try to pass Python objects.
    """

    cpp_excludeConcatenateForThisChrom(&bedfile_start_nparray[0],    &bedfile_end_nparray[0],
                                        &exclusion_start_nparray[0], &exclusion_end_nparray[0],
                                        &result_start_nparray[0],    &result_end_nparray[0],
                                        bed_size, excl_size)
    # Remember that in C++ you must actually pass a pointer to the first value, not an array.
    # Also note that the operations are performed in-place (we are passing references)
    # so we now work directly on result_start_nparray and result_end_nparray

    return None


def exclude_concatenate_for_this_chrom(chrom, exclusion, bedfile):
    """
    Subfunction of exclude_concatenate, for one chromosome only. Used for
    multiprocessing.

    Please see the documentation and code comments of exclude_concatenate
    for more information.
    """

    message('Exclusion in progress for '+str(chrom),type='DEBUG')

    ### Take PARTIAL bedfiles and exclusion : only for the current chromosome
    bedfile = bedfile[bedfile.chrom == chrom]
    exclusion = exclusion[exclusion.chrom == chrom]

    # If there is no region in either bedfile or exclusion for this chromosome,
    # abort and return an empty dataframe ; do not call C++
    if (len(bedfile) == 0) | (len(exclusion) == 0):
      return pd.DataFrame(columns = ['chrom','start','end'])

    # Convert bedfile and exclusion to 4 numpy arrays, containing only start
    # and end (since we work on a single chromosome)
    bedfile_start_nparray = bedfile[['start']].values.ravel()
    bedfile_end_nparray = bedfile[['end']].values.ravel()
    exclusion_start_nparray = exclusion[['start']].values.ravel()
    exclusion_end_nparray = exclusion[['end']].values.ravel()

    # Create partial_result arrays for the modifications to be applied to
    # (otherwise, as soon as one exclusion for one region) has been done, we
    # would be comparing values from two different coordinates systems
    result_start_nparray = np.copy(bedfile_start_nparray)
    result_end_nparray = np.copy(bedfile_end_nparray)

    # Convert those arrays to ensure that they are np.ndarray[long long, ndim=2, mode="c"]
    # NOTE : np.longlong is int64
    allarrays = [bedfile_start_nparray, bedfile_end_nparray, exclusion_start_nparray, exclusion_end_nparray, result_start_nparray, result_end_nparray]
    bedfile_start_nparray, bedfile_end_nparray, exclusion_start_nparray, exclusion_end_nparray, result_start_nparray, result_end_nparray = (np.array(array, dtype = np.longlong, order = 'C') for array in allarrays)

    # Force bed_size and excl_size into np.longs (int32)
    bed_size  = np.long(bedfile_start_nparray.shape[0])
    excl_size = np.long(exclusion_start_nparray.shape[0])

    try:
      # C++ call here (through the wrapper)
      cpp_wrapper_for_exclusion(bedfile_start_nparray, bedfile_end_nparray,
                exclusion_start_nparray, exclusion_end_nparray,
                result_start_nparray, result_end_nparray,
                bed_size, excl_size)
    except Exception as e:
      message("The C++ code for exclusion computation raised an exception. Let there be panic.", type = 'ERROR')


    # Reformat partial result back by adding the chrom again and re-making it a pandas df with proper colnames
    partial_result = pd.DataFrame(np.transpose([[chrom] * bed_size, result_start_nparray, result_end_nparray]), columns = ['chrom','start','end'])

    # Convert to numeric, just in case
    partial_result[["start", "end"]] = partial_result[["start", "end"]].apply(pd.to_numeric)


    ## Sanity checks
    # In C++, we marked regions to drop by making both coordinates equal to zero.
    # Some other cases may also result in both coordinates being equal to zero.
    # Drop those regions, as pybedtools does not like this specific case.
    to_drop = (partial_result.start == 0) & (partial_result.end == 0)
    partial_result.drop(partial_result.index[to_drop], inplace = True)

    return partial_result







def exclude_concatenate(bedfile, exclusion, nb_threads = 8):
    r"""
    When given a bedfile (in pybedtools BedFile format) and an exclusion bed file
    (in pybedtools BedFile format), will shorten the original bedfile by concatenation.
    Those two arguments must be BedTool objects from pybedtools.

    This means the regions defined in `exclusion` will be considered removed
    and the chromosme stitched back in a shorter version of itself, with the
    coordinates shifted backwards to represent that.
    Example :
        chr1 100 200
        chr1 300 400
    If we exclude 'chr1 150 300' the file becomes :
        chr1 100 150
        chr1 150 250

    Remarks :
        - The multiprocessing will pass a copy of the BED file to each process,
        which can consume a lot of RAM.


    >>> from pygtftk.stats.intersect.read_bed.read_bed_as_list import exclude_concatenate
    >>> import pybedtools
    >>> import numpy.testing as npt
    >>> part_1 = 'chr1\t100\t200\nchr2\t100\t200\n'
    >>> part_2 = 'chr3\t100\t200\nchr4\t100\t200\n'
    >>> part_3 = 'chr5\t100\t200\nchr6\t100\t200\nchr7\t100\t200\n'
    >>> e = pybedtools.BedTool(part_1 + part_2 + part_3, from_string=True)
    >>> part_1 = 'chr1\t50\t80\nchr1\t10000\t10100\nchr2\t50\t150\n'
    >>> part_2 = 'chr3\t120\t180\nchr4\t150\t250\nchr5\t250\t350\nchr6\t1\t300'
    >>> f = pybedtools.BedTool(part_1 + part_2, from_string=True)
    >>> result = exclude_concatenate(f,e)
    >>> assert str(result) == 'chr1\t50\t80\nchr1\t9900\t10000\nchr2\t50\t100\nchr4\t100\t150\nchr5\t150\t250\nchr6\t1\t200\n'

    """

    # Raw edition does not work in pybedtools, so need to use pandas dataframe
    # instead. Also, merge and sort the files before, just in case they were not.
    bedfile = bedfile.sort().merge()
    bedfile = bedfile.to_dataframe()
    exclusion = exclusion.sort().merge()
    exclusion = exclusion.to_dataframe()

    ### Exclude regions chromosome by chromosome, with multiprocessing
    all_chroms = list(exclusion.chrom) # All chromosomes in exclusion

    # To avoid wasted time in the multiprocessing, sort the chromosomes by
    # number of peaks.
    # Python's Pool map() function will split the list of arguments into chunks,
    # which can be a problem since it can result in one thread having, for
    # instance, only short chromosomes resulting in a waste of the first thread's
    # potential. To correct this, chunksize is set to 1. This will be sligtly
    # less efficient but saves time here because not all tasks are equally
    # computationally expensive (chromosomes of different lengths)
    occ = dict(Counter(all_chroms))
    all_chroms = sorted(occ.keys(), key = lambda k: occ[k])
    all_chroms.reverse()

    # TODO: if RAM turns out to be critical, do not pass the entire
    # 'exclusion' and 'bedfile' dataframes but subset by chromosome before.
    # In most use cases however it should be sufficient.

    with Pool(nb_threads) as p:
        partial_exclusion = partial(exclude_concatenate_for_this_chrom, exclusion=exclusion, bedfile=bedfile)
        list_of_partial_results = p.map(partial_exclusion, all_chroms, chunksize=1)

    result = pd.concat(list_of_partial_results, ignore_index = True)

    # Convert the dataframe back into a bedfile and return it
    result_bedfile = pybedtools.BedTool.from_dataframe(result)
    result_bedfile = result_bedfile.sort().merge() # Needed due to multiprocessing
    return result_bedfile
