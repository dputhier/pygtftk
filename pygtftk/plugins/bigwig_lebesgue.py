"""
For use with OLOGRAM in intra-set overlaps.

This will produce a bed file by performing a Lebesgue integration (along the
y - or signal - axis) of a bigWig file. The intent is to emulate a bed file
containing one region per read, where instead of a single read each region here
represents a pack of reads.

This "resampled" bed file has much less total regions/reads, meaning it can be
reasonably used in OLOGRAM to account for the signal via intra-set overlap.


Author : Quentin Ferré <quentin.q.ferre@gmail.com>

"""

"""
This is a WORK-IN-PROGRESS. To be used once intra-set overlaps are implemented in OLOGRAM.

TODO Here are the steps that would be required :

1) Remove the merging of BEDs regions (calls to bedfile.merge()), and ensure negative inter-regions distances are kept.

2) Remember that overlap_regions.find_intersection() can understand several region open at once and can for example return flags of [2,1] if A has 2 open regions.

3) Dictionary learning is non binary so it can learn proportions. But remember that sum(V**2) = 1 always for the words, so for proportion have a set that contains always 1, so it that set has a value of 0.2 in the word you know you need to multiply by 5 for example. Another idea is to look closely at the coefficients in the encoding U

4) In the display of the combination (combi_human_readable), add the factor, ie. "[2*A + B]"

"""



# # ------------------------------- PARAMETERS --------------------------------- #
# bwf = "./ENCFF431HAA_H3K36me3_K562_sub.bw"

# # Lebesgue integration : how many bins to divide the range into ?
# BINS = 10
# # Result intervals longer than this will be partitioned
# MAX_INTERVAL_LENGTH = 250

# result_filepath = bwf + ".lebesgue_"+str(BINS)+"_bins_"+str(MAX_INTERVAL_LENGTH)+"bp_max.bed"
# # ---------------------------------------------------------------------------- #


# import pyBigWig
# import numpy as np
# import time


# def one_runs(a):
#     """
#     Find successive '1' on a numpy vector.
#     """
#     # Create an array that is 1 where a is 1, and pad each end with an extra 0
#     # which is standing in for "false"
#     isone = np.concatenate(([0], np.equal(a, 1).view(np.int8), [0]))
#     absdiff = np.abs(np.diff(isone))
#     # Runs start and end where absdiff is 1.
#     ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
#     return ranges


# def inclusive_slice(start, end, step):
#     """
#     Like range, but keeps the end.
#     """
#     points = list()
#     pos = start
#     while pos < end:    # use strict inferior because we manually add the end
#         points += [pos]
#         pos += step
#     points += [end]
#     return points


# start = time.time()
# print("Processing...")

# # Read the bigWig file
# bw = pyBigWig.open(bwf)

# # Open result file
# result_file = open(result_filepath,'w+')

# # Get all chromosomes
# all_chroms = list(bw.chroms().keys())

# # Retrieve min and max values to compute the bins
# h = bw.header()
# minval = h["minVal"]
# maxval = h["maxVal"]

# ## Lebesgue integral
# # Divide the range (remove zero and add back 1 as a minimum)
# delta = int((maxval-minval)/BINS)
# bins = np.array([1]+ list(range(minval,maxval,delta))[1:])

# # TODO multiprocess

# # Do the following for each chromosome...
# for current_chromosome in all_chroms :

#     # Get all intervals for this chromosome
#     bw_intervals = bw.intervals(current_chromosome)

#     # Get all remarkable points (starts and stops)
#     remarkable_points = [r[0] for r in bw_intervals] # Get all the starts...
#     remarkable_points += [bw_intervals[-1][1] - 1] # ... and the last stop.
#     # Substract 1 from the last stop to get the last known nucleotide

#     # Convert to numpy array
#     remarkable_points = np.array(remarkable_points)

#     # Now use a sweep line
#     flags = np.array([0] * len(bins))

#     # Result numpy array
#     result = np.zeros(shape=(len(remarkable_points),len(bins)))

#     for i in range(len(remarkable_points)):
#         point = remarkable_points[i]

#         signal = bw.values(current_chromosome,point,point+1)[0]

#         # Get the indices of all bins smaller or equal to the signal
#         signal_flags = np.nonzero(bins <= signal)[0].astype(int)
#         unsignal_flags = np.nonzero(bins > signal)[0].astype(int)

#         flags[signal_flags] = 1
#         flags[unsignal_flags] = 0

#         # Record this
#         result[i,:] = flags


#     # Finally do the lebesgue integration
#     # For each column, get the ranges where there are consecutive nonzeros
#     all_intervals = list()

#     for j in range(result.shape[1]):
#         column_bin = result[:,j]
#         intervals = one_runs(column_bin)

#         # Translate the intervals back into genomic coordinates
#         intervals_genomic = remarkable_points[intervals]

#         # Record them
#         all_intervals += [tuple(i) for i in intervals_genomic]


#     # Option to truncate : intervals longer than MAX_INTERVAL_LENGTH bp should
#     # be cut every X base pairs.
#     sliced_intervals = []
#     for interval in all_intervals:
#         start, end = interval
#         sliced_points = r = inclusive_slice(start, end, MAX_INTERVAL_LENGTH)
#         cut_intervals = [(r[i],r[i+1]) for i in range(len(r)-1)]
#         sliced_intervals += cut_intervals

#     # Finally, write each (sliced) interval for all bins into a bed file
#     for interval in sliced_intervals:
#         start = interval[0]
#         end = interval[1]
#         to_write = str(current_chromosome)+'\t'+str(start)+'\t'+str(end)+'\n'
#         result_file.write(to_write)


# # TODO : to help choose the bins, return the compression loss in decibels, like for MP3 compression

# # Close result file
# result_file.close()


# end = time.time()
# print("Lebesgue integration performed in "+str(end-start)+" s.")
