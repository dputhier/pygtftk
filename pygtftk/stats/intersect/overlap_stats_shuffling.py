"""
A module to
"""



#
# # FOR MY JUPYTER KERNEL TESTING ONLY
# import os
# os.getcwd()
# os.chdir('/home/ferre/anaconda3/lib/python3.6/site-packages/pygtftk-0.9.8-py3.6-linux-x86_64.egg/pygtftk')
#
#



import numpy as np
import pandas as pd
import pybedtools
import scipy
import scipy.stats

from multiprocessing import Pool
import functools as ft
import time, sys

from pygtftk.utils import message




#
#
# # WIP : this is used during development to import Cython code
# import pyximport; pyximport.install(reload_support=True)
#
#




from pygtftk.stats.intersect import read_bed_as_list as read_bed
from pygtftk.stats.intersect import overlap_regions as oc
from pygtftk.stats.intersect import create_shuffles as cs
from pygtftk.stats.intersect import negbin_fit as nf

################################## Functions ###################################


# -------------------------------- MINIBATCH --------------------------------- #



def compute_all_intersections_minibatch(Lr1,Li1,Lr2,Li2,all_chrom1,all_chrom2,
                                        minibatch_size,
                                        use_markov_shuffling,
                                        nb_threads):
    """
    Main function. Computes a minibatch of shuffles for the given parameters.

    Lr1,Li1,Lr2,Li2,all_chrom1,all_chrom2 are all outputs from the
    bed_to_lists_of_intervals function calls right above in the code.

    'minibatch_size' (int) and 'use_markov_shuffling' (bool) are the other parameters.

    'nb_threads' is the nb of threads for the multiprocessing
    """

    # --------------------- Generate and shuffle batches  ---------------- #
    # We generate a matrix with the batches and shuffle them independantly
    # for both bed files

    # Get chromosomes that are common between both BEDs
    all_chroms = np.intersect1d(all_chrom1, all_chrom2)
    all_chroms = [str(x) for x in all_chroms] # Revert from numpy list to traditional python list

    shuffled_Lr1_batches, shuffled_Li1_batches = dict(), dict()
    shuffled_Lr2_batches, shuffled_Li2_batches = dict(), dict()

    ## Wrapper to make the code cleaner
    # Tile the list of length (repeat) as many times as we want shuffles, then
    # shuffle the rows independantly.

    # We can use a classical or a order-2 Markov shuffling. This results in different wrappers here :
    if not use_markov_shuffling :
        def batch_and_shuffle_list(l): return cs.shuffle(np.tile(l, (minibatch_size,1)))
    if use_markov_shuffling :
        def batch_and_shuffle_list(l): return cs.markov_shuffle(np.tile(l, (minibatch_size,1)), nb_threads = nb_threads)

    # Produce the shuffles on a chromosome basis
    start = time.time()
    for chrom in all_chroms :
        shuffled_Lr1_batches[chrom] = batch_and_shuffle_list(Lr1[chrom])
        shuffled_Li1_batches[chrom] = batch_and_shuffle_list(Li1[chrom])
        shuffled_Lr2_batches[chrom] = batch_and_shuffle_list(Lr2[chrom])
        shuffled_Li2_batches[chrom] = batch_and_shuffle_list(Li2[chrom])
    stop = time.time()
    message('Batch generated and shuffled in '+str(stop-start)+' s', type='DEBUG')



    # --------------------- Convert batches into BED files ----------------------- #
    start = time.time()
    batch_to_bedlist_with_params = ft.partial(cs.batch_to_bedlist,all_chroms=all_chroms, minibatch_size=minibatch_size, nb_threads = nb_threads)
    bedsA, bedsB = batch_to_bedlist_with_params(shuffled_Lr1_batches, shuffled_Li1_batches, shuffled_Lr2_batches, shuffled_Li2_batches)
    stop = time.time()
    message('Batch converted to fake beds in : '+str(stop-start)+' s', type='DEBUG')


    # --------------------- Processing intersections ----------------------------- #
    # Using our custom cython intersect, process intersection between each pair of
    # 'fake bed files'
    start = time.time()
    all_intersections = oc.compute_intersections_cython(bedsA, bedsB, all_chroms, nb_threads)
    stop = time.time()
    message('All intersections computed by our custom Cython in : '+str(stop-start)+' s', type='DEBUG')




    return all_intersections




























################################################################################

#
# Of course inmport hem one by one, like 'from pygtftk.stats.intersect import nb_fit as nf'
#
#
# PUT ALL OTHER RELEVANT SNIPPERS IN THIS DIRECTORY AND IMPORT THEM !!!

#
# # Placeholders for debugging
# from pygtftk.utils import chrom_info_as_dict
# import pybedtools
#
# chrom_len = chrom_info_as_dict(open('/home/ferre/Git/region_overlap/data/hg38_chrom_size.txt'))
# bedA = pybedtools.BedTool('/home/ferre/Git/region_overlap/data/homo_sapiens_hg38_r92_promoter.bed.gz')
# bedB = pybedtools.BedTool('/home/ferre/Git/region_overlap/data/ENCFF112BHN_H3K4me3_K562.bed.gz')
#
#
# minibatch_size, minibatch_nb =25,8
# bed_excl = pybedtools.BedTool('/home/ferre/Git/region_overlap/data/poub_excl.bed')
#
# use_markov_shuffling = False
#
# nb_threads = 8











# MUST EXPLAIN IN COMMENTS THAT THIS FILE IS THE ROOT OF ALL THE stats.intersect MODULE !!


from pygtftk.utils import message

# NO DEFAULT ARGS HERE, DEFAULTS COME FROM
def compute_overlap_stats(bedA, # corresponds to the old argument 'peak_file=region_mid_point.fn'
                        bedB, # corresponds to the old argument 'feature_bo=gtf_sub_bed'

                        chrom_len,
                        minibatch_size, minibatch_nb,
                        bed_excl,
                        use_markov_shuffling,
                        nb_threads):
    """
    Doc. Add link to the paper.

    This is the hub function to compute overlap statsitics based on Monte Carlo shuffling with integration of the inter-region lengths.


    Author : Quentin Ferré <quentin.q.ferre@gmail.com>
    """


    message('Beginning shuffling for a given set of features...')
    message('BATCHES : '+str(minibatch_nb)+' batches of '+str(minibatch_size)+' shuffles', type='DEBUG')
    message('Total number of shuffles : '+ str(minibatch_nb*minibatch_size), type='DEBUG')
    message('NB_THREADS = ' + str(nb_threads), type='DEBUG')


    # --------------------- Read list of intervals --------------------------- #

    start = time.time()

    # Just in case, force type and merge bedA among itself and bedB same.
    bed_A_as_pybedtool = pybedtools.BedTool(bedA).merge()
    bed_B_as_pybedtool = pybedtools.BedTool(bedB).merge()


    # If there is an exclusion to be done, do it
    if bed_excl != None:
        exclusion = pybedtools.BedTool(bed_excl) # Just in case

        chrom_len = read_bed.exclude_chromsizes(exclusion,chrom_len) # Shorten the chrom_len only once, and separately

        bed_A_as_pybedtool = read_bed.exclude_concatenate(bed_A_as_pybedtool, exclusion, chrom_len)
        bed_B_as_pybedtool = read_bed.exclude_concatenate(bed_B_as_pybedtool, exclusion, chrom_len)




    # Raise exception if there are less than 2 remainin regions in bedA and bedB (or rather raise an exception with gtftk.error)
    if (len(bed_A_as_pybedtool) < 2) | (len(bed_B_as_pybedtool) < 2):
        raise ValueError('Less than 2 remaining regions in either bed.')



    # Proper reading of the bed file as a list of intervals
    Lr1, Li1, all_chrom1 = read_bed.bed_to_lists_of_intervals(bed_A_as_pybedtool,chrom_len)
    Lr2, Li2, all_chrom2 = read_bed.bed_to_lists_of_intervals(bed_B_as_pybedtool,chrom_len)
    stop = time.time()
    message('BED files read as lists of intervals in '+str(stop-start)+' s', type='DEBUG')


    grand_start = time.time()

    ################################### MINIBATCH  #################################
    # Generate all intersections for a shuffled batch of size n

    minibatches = [minibatch_size for i in range(minibatch_nb)]
    all_intersections = list()
    for k in range(len(minibatches)):
        #message('\n'+'--- Minibatch number : '+str(k)+' ---'+'\n')

        # Display of current progress
        message("--- Minibatch nb. : " + str(k+1)+" / "+str(minibatch_nb))
        #sys.stdout.flush()

        all_intersections = all_intersections + compute_all_intersections_minibatch(Lr1, Li1, Lr2, Li2, all_chrom1, all_chrom2,
                                                                    minibatches[k], use_markov_shuffling, nb_threads)
    message('All intersections have been generated.')



    # ----------------- Compute statistics on the intersections ------------------ #

    start = time.time()
    with Pool(nb_threads) as p:
        stats = p.map(read_bed.compute_stats_for_intersection, all_intersections)
    stop = time.time()
    message('Statistics on overlaps computed in : '+str(stop-start)+' s', type='DEBUG')

    # Unpack the stats
    bp_overlaps = [s[0] for s in stats]
    unlisted_bp_overlaps = [item for sublist in bp_overlaps for item in sublist]
    summed_bp_overlaps = [sum(x) for x in bp_overlaps]
    intersect_nbs = [s[1] for s in stats]



    #### Fitting of a Negative Binomial distribution on the shuffles
    # Only relevant for classical shuffle, not Markov


    # TODO : also do that if mean<100 because it can't be approximated by a normal
    if use_markov_shuffling :
        ps = pn = -1 # TODO explain why -1 is returned !

    else :
        # Renaming esperances and variances
        esperance_fitted_summed_bp_overlaps, variance_fitted_summed_bp_overlaps = np.mean(summed_bp_overlaps), np.var(summed_bp_overlaps)
        esperance_fitted_intersect_nbs, variance_fitted_intersect_nbs = np.mean(intersect_nbs), np.var(intersect_nbs)

        # Check that there is a good adjustment
        # We use a normal law as NB -> norm for large N, but this is irrelevant
        # if the mean is under 50 roughly (should not happen with this kind of
        # data, but just in case)
        if esperance_fitted_summed_bp_overlaps<50 :
            ps = pn = -1
        else :
            ps = nf.check_negbin_adjustment(summed_bp_overlaps,esperance_fitted_summed_bp_overlaps,variance_fitted_summed_bp_overlaps).pvalue
            pn = nf.check_negbin_adjustment(intersect_nbs,esperance_fitted_intersect_nbs,variance_fitted_intersect_nbs).pvalue




    # ---------------------------- True intersections ---------------------------- #
    # Now, calculating the actual p-value for the number of intersections and the
    # total number of overlapping base pairs



    ## True intersection
    # A = pybedtools.BedTool(BED_A) # Useless since bedA and bedB are already pybedtools objects
    # B = pybedtools.BedTool(BED_B)
    true_intersection = bedA.intersect(bedB)

    true_intersect_nb = len(true_intersection)
    true_bp_overlaps = sum([x.length for x in true_intersection])


    # Compute the p-values using the distribution fitted on the shuffles

    # Do not do this for the Markov shuffling, as it is likely a multi-variable fit (see notes)


    # We can only use a Neg Binom p-val if we can fit it, and that is not the case for
    # the Markov shuffle : we must use an empirical p-value
    if use_markov_shuffling :
        pval_intersect_nb = nf.empirical_p_val(true_intersect_nb,intersect_nbs)
        pval_bp_overlaps = nf.empirical_p_val(true_bp_overlaps,summed_bp_overlaps)

    else :
        pval_intersect_nb = 1 - np.exp(nf.log_nb_pval(true_intersect_nb,esperance_fitted_intersect_nbs,variance_fitted_intersect_nbs))
        pval_bp_overlaps = 1 - np.exp(nf.log_nb_pval(true_bp_overlaps,esperance_fitted_summed_bp_overlaps,variance_fitted_summed_bp_overlaps))







    grand_stop = time.time()

    message('--- Total time : '+str(grand_stop-grand_start)+' s ---')
    message('Total time does not include BED reading, as it does not scale with batch size.', type='DEBUG')

























    # result will be a dict or a pandas, see the dataframe d in vanilla to see the format





    # Just adpat the original output_intersect to return it, and see how to merge it
    # with Denis' code


    from collections import OrderedDict

    result = OrderedDict()

    result['nb_intersections_esperance_shuffled'] = np.mean(intersect_nbs)
    result['nb_intersections_variance_shuffled'] = np.var(intersect_nbs)
    result['nb_intersections_fit'] = pn
    result['nb_intersections_true'] = true_intersect_nb
    result['nb_intersections_pvalue'] = pval_intersect_nb
    result['summed_bp_overlaps_esperance_shuffled'] = np.mean(summed_bp_overlaps)
    result['summed_bp_overlaps_variance_shuffled'] = np.var(summed_bp_overlaps)
    result['summed_bp_overlaps_fit'] = ps
    result['summed_bp_overlaps_true'] = true_bp_overlaps
    result['summed_bp_overlaps_pvalue'] = pval_bp_overlaps



    return result








#
#
#
#
#
#
# poub = compute_overlap_stats(bedA, bedB, chrom_len, minibatch_size, minibatch_nb,bed_excl,use_markov_shuffling,nb_threads)
#
#
