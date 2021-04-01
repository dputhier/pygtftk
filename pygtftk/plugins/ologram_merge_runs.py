#!/usr/bin/env python
"""
Merge a set of OLOGRAM runs into a single run and recalculates statistics based on it.

This treats each run as a "superbatch".
"""

import argparse
import os
import re

import numpy as np
import pandas as pd

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.utils import message

from pygtftk.stats import negbin_fit as nf

__updated__ = ''' 2021-04-01 '''

__notes__ = """
-- Statistics are recalculated by conflating the distributions with a weighting based on the number of runs (see source for the precise formula).
"""


def make_parser():
    """The main argument parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfiles',
                            help="Complete paths to the OLOGRAM output text files",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('txt')),
                            nargs='+')

    parser_grp.add_argument('-o', '--output',
                            help="Merged output file.",
                            default=None,
                            nargs=None,
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('txt')),
                            required=True)

    return parser



def get_conflated_moments(mean1, mean2, var1, var2, n1, n2):
    r"""
    Get the variance and mean of the conflated distribution of 1 and 2. 
    'mean' and 'var' are their individual means and vairances, and 'n' are the sample sizes.

    Example:

    >>> assert get_conflated_moments(0,1,0,1,0,1) == (1,1)
    >>> assert get_conflated_moments(1,1,0,0,1,1) == (1,1)
    >>> assert get_conflated_moments(0,20,0,0,1,1) == (10, 200)

    """

    # If n1 or n2 is zero, disregard them
    if (n1 == 0) & (n2 == 0): return 0,0
    elif n1 == 0: return mean2, var2
    elif n2 == 0: return mean1, var1

    # Weighted average for the mean
    new_mean = (n1*mean1 + n2*mean2)/(n1+n2)

    # Conflated variance
    q1 = (n1-1)*var1 + n1*(mean1**2)
    q2 = (n2-1)*var2 + n2*(mean2**2)
    qc = q1 + q2
    new_var = (qc - (n1+n2)*(new_mean**2)) / (n1+n2-1)

    return new_mean, new_var



def ologram_merge_runs(inputfiles=None,
                        output=None):

    # -------------------------------------------------------------------------
    # Loop over input files
    # -------------------------------------------------------------------------

    runs_to_be_merged = []


    for _, infile in enumerate(inputfiles):
        message("Reading file : " + infile.name)
        # Read the dataset into a temporay dataframe
        run_as_df = pd.read_csv(infile, sep='\t', header=0, index_col='feature_type')
        # Take the feature_type as index for easier later merging
 
        # Add the df to the list to be subsequently merged
        runs_to_be_merged += [run_as_df]


    ## Prepare an empty dataframe for the merged run
    final_index = []
    for run in runs_to_be_merged: final_index += run.index.tolist()
    final_index = sorted(list(set(final_index)))
    final_columns = runs_to_be_merged[0].columns # Should always be the same columns anyways

    
    merged_run = pd.DataFrame(columns = final_columns, index = final_index)

    # -------------------------------------------------------------------------
    # Merging runs and recalculating stats
    # -------------------------------------------------------------------------

    ## For each run, spill its contents inside the merged_run dataframe

    # Initialize to zero
    for combi, _ in merged_run.iterrows():
        merged_run.loc[combi, 'nb_intersections_expectation_shuffled'] = 0
        merged_run.loc[combi, 'nb_intersections_variance_shuffled'] = 0
        merged_run.loc[combi, 'summed_bp_overlaps_expectation_shuffled'] = 0
        merged_run.loc[combi, 'summed_bp_overlaps_variance_shuffled'] = 0

    # Process each run
    runs_already_merged = 0

    for run in runs_to_be_merged:
        message("Treating a run... so far " + str(runs_already_merged) + " are complete.")

        total_combis = run.shape[0]
        i = 0

        for combi, row in run.iterrows():
            
            ## Combine the means and variance with the runs previously merged

            previous_merged_N_mean = merged_run.loc[combi, 'nb_intersections_expectation_shuffled']
            previous_merged_N_var =  merged_run.loc[combi, 'nb_intersections_variance_shuffled']
            previous_merged_S_mean = merged_run.loc[combi, 'summed_bp_overlaps_expectation_shuffled']
            previous_merged_S_var = merged_run.loc[combi, 'summed_bp_overlaps_variance_shuffled'] 

            current_N_mean = row['nb_intersections_expectation_shuffled']
            current_N_var =  row['nb_intersections_variance_shuffled']
            current_S_mean = row['summed_bp_overlaps_expectation_shuffled']
            current_S_var = row['summed_bp_overlaps_variance_shuffled'] 


            # Get the new moments
            new_S_mean, new_S_var = get_conflated_moments(
                mean1 = previous_merged_S_mean, mean2 = current_S_mean,
                var1 = previous_merged_S_var, var2 = current_S_var,
                n1 = runs_already_merged, n2 = 1)
            new_N_mean, new_N_var = get_conflated_moments(
                mean1 = previous_merged_N_mean, mean2 = current_N_mean,
                var1 = previous_merged_N_var, var2 = current_N_var,
                n1 = runs_already_merged, n2 = 1)

            # Overwrite the moments
            merged_run.loc[combi, 'nb_intersections_expectation_shuffled'] = new_N_mean
            merged_run.loc[combi, 'nb_intersections_variance_shuffled'] = new_N_var
            merged_run.loc[combi, 'summed_bp_overlaps_expectation_shuffled'] = new_S_mean
            merged_run.loc[combi, 'summed_bp_overlaps_variance_shuffled'] = new_S_var       


            # True intersections stay the same every time
            merged_run.loc[combi, 'nb_intersections_true'] = row['nb_intersections_true']
            merged_run.loc[combi, 'summed_bp_overlaps_true'] = row['summed_bp_overlaps_true']

            i += 1
            message("Merged combi "+str(i)+" / "+str(total_combis)+" for this run.", type = "DEBUG")

        # Used for the subsequent weighting
        runs_already_merged = runs_already_merged + 1

    message("All runs read. Proceeding to merge statistics.")



    ## At the end, recalculate fold change and p-value    
    i = 0
    for combi, row in merged_run.iterrows():      

        # Do not divide by zero ! Use the true value as fold change if needed
        expectation_fitted_intersect_nbs = merged_run.loc[combi, 'nb_intersections_expectation_shuffled']
        true_intersect_nb = merged_run.loc[combi, 'nb_intersections_true']
        if expectation_fitted_intersect_nbs == 0: ni_fc = true_intersect_nb  
        else: ni_fc = true_intersect_nb / expectation_fitted_intersect_nbs
        if ni_fc != 0: ni_fc = np.log2(ni_fc)
        merged_run.loc[combi, 'nb_intersections_log2_fold_change'] = '{:.5f}'.format(ni_fc) 

        expectation_fitted_summed_bp_overlaps = merged_run.loc[combi, 'summed_bp_overlaps_expectation_shuffled']
        true_bp_overlaps = merged_run.loc[combi, 'summed_bp_overlaps_true']
        if expectation_fitted_summed_bp_overlaps == 0: sbp_fc = true_bp_overlaps
        else: sbp_fc = true_bp_overlaps / expectation_fitted_summed_bp_overlaps
        if sbp_fc != 0: sbp_fc = np.log2(sbp_fc)  # Apply log transformation
        merged_run.loc[combi, 'summed_bp_overlaps_log2_fold_change'] = '{:.5f}'.format(sbp_fc)

        # Fit qualities are not applicable here
        # TODO: Recalculate them somehow ?
        merged_run.loc[combi,'summed_bp_overlaps_negbinom_fit_quality'] = np.nan
        merged_run.loc[combi,'nb_intersections_negbinom_fit_quality'] = np.nan

        
        # Recalculate the p values
        variance_fitted_intersect_nbs = merged_run.loc[combi, 'nb_intersections_variance_shuffled'] 
        variance_fitted_summed_bp_overlaps = merged_run.loc[combi, 'summed_bp_overlaps_variance_shuffled']

        pval_intersect_nb = nf.negbin_pval(true_intersect_nb, expectation_fitted_intersect_nbs,
                                            variance_fitted_intersect_nbs, ft_type=combi)
        pval_bp_overlaps = nf.negbin_pval(true_bp_overlaps, expectation_fitted_summed_bp_overlaps,
                                            variance_fitted_summed_bp_overlaps, ft_type=combi)
                                            
        merged_run.loc[combi, 'nb_intersections_pvalue'] = '{0:.4g}'.format(pval_intersect_nb)
        merged_run.loc[combi, 'summed_bp_overlaps_pvalue'] = '{0:.4g}'.format(pval_bp_overlaps)

        i+=1
        message("Statistics done for combi "+str(i)+" / "+str(total_combis)+" for this run.", type = "DEBUG")


    ## Finally write the merged df to a file
    message("Writing output file.")
    merged_run.insert(0, "feature_type", merged_run.index) # Remember to transform the index back into the feature_type column
    merged_run.to_csv(open(output.name, 'w'), sep="\t", header=True, index=False)







def main():
    """The main function."""

    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    ologram_merge_runs(**args)


if __name__ == '__main__':
    main()
else:

    # 'Bats' tests
    test = '''
    #ologram_merge_runs: get example files
    @test "ologram_merge_runs_0" {
         result=`gtftk get_example -d ologram_1 -f '*'`
      [ "$result" = "" ]
    }

    #ologram_merge_runs: merge them
    # Note that this is an example of what not to do since those files were done for different queries.
    # I treat them as the same query, since this is just to test the principle.
    @test "ologram_merge_runs_1" {
         result=`gtftk ologram_merge_runs -i H3K4me3_ologram_stats.tsv H3K36me3_ologram_stats.tsv H3K79me2_ologram_stats.tsv -o merged_ologram_runs.tsv`
      [ "$result" = "" ]
    }

    #ologram_merge_runs: test value
    @test "ologram_merge_runs_2" {
        result=`cat merged_ologram_runs.tsv | grep "transcript" | cut -f 2`
        [ "$result" = "87.36" ]
    }
    '''

    cmd = CmdObject(name="ologram_merge_runs",
                    message="Merge ologram runs, treating each as a superbatch.",
                    parser=make_parser(),
                    fun=os.path.abspath(__file__),
                    desc=__doc__,
                    group="ologram",
                    notes=__notes__,
                    updated=__updated__,
                    test=test)
    