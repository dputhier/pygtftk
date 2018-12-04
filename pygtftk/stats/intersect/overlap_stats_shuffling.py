from pytgtftk.stats.intersect import *#
# Of course inmport hem one by one, like 'from pygtftk.stats.intersect import nb_fit as nf'


# PUT ALL OTHER RELEVANT SNIPPERS IN THIS DIRECTORY AND IMPORT THEM !!!


# MUST EXPLAIN IN COMMENTS THAT THIS FILE IS THE ROOT OF ALL THE stats.intersect MODULE !!


from pygtftk.utils import message


def compute_overlap_stats(peak_file=region_mid_point.fn, # bedA
                                     feature_bo=gtf_sub_bed, # bedB

                                     chrom_len=chrom_len): # chromsizes
    """
    Doc. Add link to the paper.

    This is the hub function to compute overlap statsitics based on Monte Carlo shuffling with integration of the inter-region lengths.


    Author : Quentin Ferré <quentin.q.ferre@gmail.com>
    """


    # ALL CODE HERE

    # TODO all my log stuff will be a message instead : message('blabla')


    # result will be a dict or a pandas, see the dataframe d in vanilla to see the format
    return result
