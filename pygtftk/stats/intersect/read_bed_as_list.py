"""
Describe the module
"""

import pybedtools
import numpy as np
import pandas as pd


# -------------------------- Reading bed files ------------------------------- #

def bed_to_lists_of_intervals(bed,chromsizes):
    """
    Reads a bed file and returns two dictionaries, with the list of region lengths
    and interregions lengths (resp. Lr and Li).
    The bed file must be a pybedtools object.
    The dictionaries group the lists of lengths by chromosome.

    A dictionary of chromosome sizes must also be provided to compute the
    distance between the last feature and the chromosome end

    For example, Li['chr1'] is the list of distances between regions in chr1.
    """
    # bed = pybedtools.BedTool(bed_path)

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
            lr.append(bed.at[f,'end'] - bed.at[f,'start'])
            li.append(bed.at[f,'start'] - previous_feature_stop)
            previous_feature_stop = bed.at[f,'end']

        # Add the missing inter-region distance between the last feature and the chromosome end
        last_li = chromsizes[chrom] - previous_feature_stop
        if last_li < 0:
            last_li = 0
            print('Warning - You have a bed file with features after the end of chromosome "'+str(chrom)+'" !')
        li.append(last_li)


        Lr[chrom] = np.array(lr)
        Li[chrom] = np.array(li)

    return Lr, Li, all_chrom



# -------------------------- Exclusion and concatenation ------------------------------- #

def exclude_chromsizes(exclusion,chromsizes):
    """
    Shortens the chromsome sizes (given as a dictionary) by the total length of
    each excluded region (given as a BedTool file).
    """
    exclusion = exclusion.to_dataframe()
    for ie, excl in exclusion.iterrows() :
        excl_length = abs(excl['end'] - excl['start'])
        chromsizes[excl['chrom']] = chromsizes[excl['chrom']] - excl_length
    return chromsizes




def exclude_concatenate(bedfile, exclusion, chromsizes):
    """
    When given a bedfile (in pybedtools BedFile format) and an exclusion bed file
    (in pybedtools BedFile format), will shorten the original bedfile by concatenation.
    The chromsizes (must be a dictionary) will also be appropriately shortened.
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

    Remark : This version is highly inefficient (1 second per excluded feature)
    but is only run once per analysis, so it will be improved later.
    """


    # Raw edition does not work in pybedtools, so need to use pandas dataframe instead.
    bedfile = bedfile.to_dataframe()
    exclusion = exclusion.to_dataframe()

    # WARNING Must use a copy and not remove elements one by one, because that
    # would shift the position and now you are comparing positions in two
    # different coodinates systems (between 'exclusion' with the original ones
    # and the shifted 'bedfile')
    result = bedfile.copy()

    # For each region in 'exclusion' :
    for ie, excl in exclusion.iterrows() :

        excl_length = abs(excl['end'] - excl['start'])

        # Gain some time by selecting only the rows on which we will operate
        filtered_rows = bedfile[(bedfile.chrom == excl.chrom) & (bedfile.end >= excl['start'])]
        filtered_rows = filtered_rows.index

        ### TREATING BEDFILE
        for i in filtered_rows:

            # Rq : I use '<' and '>=' so do not use '<=' or '>' if you modify this, else not all conditions will be covered

            # WARNING Since this is an iterative algorithm, you must always modify the values from result and always write the new value of result
            # as a function of the previous value of result, otherwise you are comparing positions from two different coordinates sets.

            # all regions where region_start is under exclu_start but region_end is higher thean exclu_start : truncate by setting region_end to exclu_start
            if (bedfile.at[i,'start'] < excl['start']) & (bedfile.at[i,'end'] >= excl['start']) :
                truncate_by = bedfile.at[i,'end'] - excl['start']
                result.at[i,'end'] = result.at['end'] - truncate_by

            # all which contain the excluded region (start before and end after) : shorten the end by the region length
            elif (bedfile.at[i,'start'] < excl['start']) & (bedfile.at[i,'end'] >= excl['end']):
                result.at[i,'end'] = result.at[i,'end'] - excl_length

            # TODO it seems hromsomes are shortened twice instead of once !!

            # all regions where region_start > excl_start but region_end < excl_end (so are included) : eliminate those
            elif (bedfile.at[i,'start'] >= excl['start']) & (bedfile.at[i,'end'] < excl['end']) :
                result.drop(i, inplace=True)


            # all regions where region_start is higher than excl_start but lower than excl_end and region_end is higher than excl_end : truncate by setting region_start to excl_end and also region_end = region_end - nb_of_nt_of_region_that_are_in_excl
            elif (bedfile.at[i,'start'] >= excl['start']) & (bedfile.at[i,'start'] < excl['end']) & (bedfile.at[i,'end'] >= excl['end']) :
                forward_by = bedfile.at[i,'start'] - excl['start']
                result.at[i,'start'] = result.at[i,'start'] - forward_by

                nb_of_bp_of_region_that_are_in_excl = (excl['end'] - bedfile.at[i,'start'])
                result.at[i,'end'] = result.at[i,'end'] - nb_of_bp_of_region_that_are_in_excl

            # all regions where region_start and region_end are both higher than excl_end : move by setting region_start = region_start - excl_length and region_end = region_end - excl_length
            elif (bedfile.at[i,'start'] >= excl['end']) :
                result.at[i,'start'] = result.at[i,'start'] - excl_length
                result.at[i,'end'] = result.at[i,'end'] - excl_length

    # Convert the dataframe back into a bedfile and return it
    result_bedfile = pybedtools.BedTool.from_dataframe(result)
    return result_bedfile




# ----------------- Compute statistics on the intersections ------------------ #


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
