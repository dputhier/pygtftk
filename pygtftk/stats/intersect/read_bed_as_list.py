"""
A set of functions to turn a BED file into a list of intervals, and exclude
certain regions to create concatenated sub-chromosomes.
"""

from multiprocessing import Pool
from functools import partial
from collections import Counter

import pybedtools
import numpy as np
import pandas as pd

from pygtftk.utils import message




################################################################################
# -------------------------- Reading bed files ------------------------------- #
################################################################################

def bed_to_lists_of_intervals(bed, chromsizes):
    """
    Reads a bed file (as a pybedtools.BedTool object) and returns respectively
    two dictionaries, with the list of region lengths and interregions lengths
    (resp. Lr and Li), as well as a list of all chromosomes.

    The dictionaries group the lists of lengths by chromosome.

    A dictionary of chromosome sizes must also be provided to compute the
    distance between the last feature and the chromosome end.

    For example, Li['chr1'] is the list of distances between regions in chr1.

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.stats.intersect.read_bed_as_list import bed_to_lists_of_intervals
    >>> import pybedtools
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

    # Convert bedfile to pandas array
    bed = bed.to_dataframe()

    Lr = dict()
    Li = dict()

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
        last_li = chromsizes[chrom] - previous_feature_stop
        if last_li < 0:
            last_li = 0
            message('Warning - You have a bed file with features after the end of chromosome "'+str(chrom)+'" !', type='INFO')
        li.append(last_li)


        Lr[chrom] = np.array(lr)
        Li[chrom] = np.array(li)

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
    >>> from pygtftk.stats.intersect.read_bed_as_list import exclude_chromsizes
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




def exclude_concatenate_for_this_chrom(chrom,exclusion,bedfile):
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


    # WARNING Must use a copy and not remove elements one by one, because that
    # would shift the position and now you are comparing positions in two
    # different coordinates systems (between 'exclusion' with the original ones
    # and the shifted 'bedfile')
    partial_result = bedfile.copy()

    # For each region in 'exclusion' :
    for _, excl in exclusion.iterrows():

        excl_length = abs(excl['end'] - excl['start'])

        # Gain some time by selecting only the rows on which we will operate
        filtered_rows = bedfile[(bedfile.chrom == excl.chrom) & (bedfile.end >= excl['start'])].index

        ### TREATING BEDFILE
        for i in filtered_rows:

            # Rq : I use '<' and '>=' so do not use '<=' or '>' if you modify this, else not all conditions will be covered

            # WARNING Since this is an iterative algorithm, we must always
            # compute the conditions and deltas from the old values in bedfile,
            # but modify (ie. apply deltas) the values from result by always
            # writing the new value of result as a function of the previous
            # value of result, otherwise you are comparing positions from two
            # different coordinates sets.

            # For sanity check : do not check a line if it has been removed
            check_for_zero = True


            # all regions where region_start is under exclu_start but region_end is higher than exclu_start BUT lower than excl_end: truncate by setting region_end to exclu_start
            if (bedfile.at[i, 'start'] < excl['start']) & (excl['end'] >= bedfile.at[i, 'end'] >= excl['start']):
                truncate_by = bedfile.at[i, 'end'] - excl['start']
                partial_result.at[i, 'end'] = partial_result.at[i, 'end'] - truncate_by

            # all which contain the excluded region (start before and end after) : shorten the end by the region length
            elif (bedfile.at[i, 'start'] < excl['start']) & (bedfile.at[i, 'end'] >= excl['end']):
                partial_result.at[i, 'end'] = partial_result.at[i, 'end'] - excl_length

            # all regions where region_start > excl_start but region_end < excl_end (so are included) : eliminate those
            # Warning : must be before the 5th test
            elif (bedfile.at[i, 'start'] >= excl['start']) & (bedfile.at[i, 'end'] < excl['end']):
                check_for_zero = False
                partial_result.drop(i, inplace=True)

            # all regions where region_start is higher than excl_start but lower than excl_end and region_end is higher than excl_end : truncate by setting region_start to excl_end and also region_end = region_end - nb_of_nt_of_region_that_are_in_excl
            elif (bedfile.at[i, 'start'] >= excl['start']) & (bedfile.at[i, 'start'] < excl['end']) & (bedfile.at[i, 'end'] >= excl['end']):

                # Compute some utils
                region_length_before_truncating = partial_result.at[i, 'end'] - partial_result.at[i, 'start']
                nb_of_bp_of_region_that_are_in_excl = (excl['end'] - bedfile.at[i, 'start'])

                # Move start point
                forward_by = bedfile.at[i, 'start'] - excl['start']
                partial_result.at[i, 'start'] = partial_result.at[i, 'start'] - forward_by

                # Move end point to 'new start point + new length'
                new_length = region_length_before_truncating - nb_of_bp_of_region_that_are_in_excl
                partial_result.at[i, 'end'] = partial_result.at[i, 'start'] + new_length

            # all regions where region_start and region_end are both higher than excl_end : move by setting region_start = region_start - excl_length and region_end = region_end - excl_length
            elif (bedfile.at[i, 'start'] >= excl['end']):
                partial_result.at[i, 'start'] = partial_result.at[i, 'start'] - excl_length
                partial_result.at[i, 'end'] = partial_result.at[i, 'end'] - excl_length


            # WARNING pybedtools does not like when both start and end are equal to zero. Set them to 1 if that's the case.
            if check_for_zero : # Do not check if the line has been dropped
                if (partial_result.at[i, 'start'] == 0) & (partial_result.at[i, 'end'] == 0):
                    partial_result.at[i, 'start'] = 1
                    partial_result.at[i, 'end'] = 1


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
        - This version is highly inefficient (1 second per excluded feature)
        but is only run once per analysis, so it will be improved later.
        - The multiprocessing will pass a copy of the BED file to each process,
        which can consume a lot of RAM.

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.stats.intersect.read_bed_as_list import exclude_concatenate
    >>> import pybedtools
    >>> import numpy.testing as npt
    >>> f = pybedtools.BedTool(get_example_file("simple","bed")[0])
    >>> e = pybedtools.BedTool('chr1\t12\t45',from_string=True)
    >>> result = exclude_concatenate(f,e)
    >>> assert str(result[0]) == 'chr1\t10\t17\n'
    """

    # Raw edition does not work in pybedtools, so need to use pandas dataframe instead.
    # Also, merge and sort the files before, just in case they were not.
    bedfile = bedfile.sort().merge()
    bedfile = bedfile.to_dataframe()
    exclusion = exclusion.sort().merge()
    exclusion = exclusion.to_dataframe()

    ### Exclude regions chromosome by chromosome, with multiprocessing
    all_chroms = list(exclusion.chrom) # All chromosomes in exclusion

    # To avoid wasted time in the multiprocessing, sort the chromosomes by number of peaks
    # Furthermore, python Pool map() function will split the list of arguments into chunks which can be a problem since it can result in one thread having only short chromosomes
    # and one only long chromosome, resulting in wasting the first thred's potential. To correct this, chunksize is set to 1. This will be sligtly less efficient but saves time
    # here because not all tasks are as computationally expensive.
    occ = dict(Counter(all_chroms))
    all_chroms = sorted(occ.keys(), key = lambda k: occ[k])
    all_chroms.reverse()

    # TODO for later : if RAM turns out to be critical, do not pass the entire
    # 'exclusion' and 'bedfile' dataframes but subset by chromosome before.
    # In most use cases however it should be sufficient.

    with Pool(nb_threads) as p:
        partial_exclusion = partial(exclude_concatenate_for_this_chrom,exclusion=exclusion,bedfile=bedfile)
        list_of_partial_results = p.map(partial_exclusion, all_chroms, chunksize=1)

    result = pd.concat(list_of_partial_results, ignore_index = True)

    # Convert the dataframe back into a bedfile and return it
    result_bedfile = pybedtools.BedTool.from_dataframe(result)
    result_bedfile = result_bedfile.sort().merge() # Needed due to multiprocessing
    return result_bedfile





################################################################################
# ----------------- Compute statistics on the intersections ------------------ #
################################################################################

def compute_stats_for_intersection(myintersect):
    """
    Wrapper to compute all stats we could want on a single intersect result object.
    The argument (myintersect) is a single bedfile, either as a pybedtools intersect
    result, or a list of tuples.
    """
    bp_overlap = [x[2] - x[1] for x in myintersect]
    intersect_nb = len(myintersect)
    stats = (bp_overlap, intersect_nb)
    return stats
