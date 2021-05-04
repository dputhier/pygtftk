"""
This module contains the MODL algorithm, an itemset mining algorithm that broadly consists of two steps.

Considering an input matrix with transactions/regions.intersections in lines and
elements/items/sets in columns, the algorithm will : (I) learn (MiniBatchDictionaryLearning)
several factorisations with different sparsity constraints to build a library of
candidate words (eg. (1,1,0,0,1)) and (II) select the best words to rebuild the 
original matrix using a greedy algorithm, by trying to rebuild the original 
matrix with a subset of selected words, adding the word that best improves the 
rebuilding to said subset each time.

Dictionary learning is an optimization problem solved by alternatively updating the 
sparse code, as a solution to multiple Lasso problems, considering the dictionary fixed,
and then updating the dictionary to best fit the sparse code.

In OLOGRAM, we use this on the matrix of true overlap flags to find usual common overlaps.

This is mostly useful if there are many files to reduce the number of displayed combinations,
and to counteract the effect of noise on the data. Unlike classical association
rules mining algorithms, this focuses on mining complexes and correlation groups (item sets).

Author : Quentin FERRE <quentin.q.ferre@gmail.com>

This code is available under the GNU GPL license.
"""

import copy, time, os

import numpy as np
import pandas as pd
import random

from pygtftk.stats.intersect.modl import tree   

from sklearn.decomposition import SparseCoder
from sklearn.decomposition import MiniBatchDictionaryLearning
from sklearn.model_selection import KFold

from pygtftk import utils
from pygtftk.utils import message



# ---------------------------------------------------------------------------- #
#                              Utility functions                               #
# ---------------------------------------------------------------------------- #

def test_data_for_modl(nflags = 1000, number_of_sets = 6, noise = 0, cor_groups = [(0,1),(0,1,2,3),(4,5)]):
    """
    Generate testing data for the dictionary learning module.
    """

    # Correlation groups
    if max(max(cor_groups)) > number_of_sets - 1 : raise ValueError('Bad correlation groups')
    # NOTE REMEMBER TO ADD GROUPS OF DIFF SIZES, AND INCLUDING PREVIOUS ONES ? eg : 1,2 and 1,2,3,4

    # Generate some flags
    all_flags = list()

    # For each wanted flag...
    for _ in range(nflags):
        flags = [0] * number_of_sets    # Initialize flags

        # Pick a correlation group
        chosen_group = random.choice(cor_groups)
        for set in chosen_group : flags[set] += 1       # Add the corresponding flags

        # NOISE : Flip each random flag with a probability of `noise`
        for i, f in enumerate(flags):
            cointoss = np.random.uniform()
            if cointoss < noise:
                flags[i] = not f # Flip 0 to 1 and 1 to 0

        # Add to all_flags
        all_flags += [flags]

    # Turn the data into a NumPy array
    result = np.array(all_flags)
    return result



def squish_matrix(x, abundance_threshold = 0, shuffle = True, smother = True):
    r"""
    To reduce redundancy in the matrix lines, take all unique rows of X and 
    build a squished matrix where each line now has the square root of its 
    original abundance divided by sqrt(abundance of most rate), but not lower 
    than abundance_threshold, 1/10000 by default. We use the square root of 
    those abudances instead to dimnish the emphasis on the most frequent 
    combinations by default (smothering).   

    >>> import numpy as np
    >>> from pygtftk.stats.intersect.modl.dict_learning import squish_matrix
    >>> X = np.array([[1,1,0,0]]*1000 + [[0,0,1,1]]*100)
    >>> X_squished = squish_matrix(X, shuffle = False, smother = True)
    >>> np.testing.assert_equal(X_squished, np.array([[0,0,1,1]]*1 + [[1,1,0,0]]*4)) # Note that the rows have been sorted by abundance   
    
    """

    # Get all unique rows and their counts
    all_rows, counts_per_row = np.unique(x, axis=0, return_counts = True)
    min_abundance = abundance_threshold*x.shape[0]

    # Use sqrt to reduce difference between most and least abundant
    # Only if smother is True
    if smother:
        min_abundance = np.sqrt(min_abundance)
        counts_per_row = np.sqrt(counts_per_row)
        
    # Divide counts by lowest count observed (but never divide by under abundance_threshold)
    minimal_count = max(np.min(counts_per_row), min_abundance)
    counts_relative = counts_per_row / minimal_count

    # Round up
    counts_relative = np.ceil(counts_relative).astype(int)

    # Rebuild a new matrix with those counts
    squished_matrix = []     
    for row, count in zip(all_rows, counts_relative):

        squished_matrix += [row] * count


    final_matrix = np.array(squished_matrix)
    if shuffle: np.random.shuffle(final_matrix) # Shuffle it for good measure

    return final_matrix






# ---------------------------------------------------------------------------- #
#                                  Subroutines                                 #
# ---------------------------------------------------------------------------- #

# Tree model for the library of candidate words
from pygtftk.stats.intersect.modl import tree as tree

# Subroutines of the main algorithm
from pygtftk.stats.intersect.modl import subroutines as modl_subroutines




# ---------------------------------------------------------------------------- #
#                                    MAIN                                      #
# ---------------------------------------------------------------------------- #


class Modl:
    r"""
    This class encapsulates the MODL approach :

    Takes as input a matrix of flags, with one flag per intersection,
    and returns the list of interesting combis.

    For more details, please see the module description of pygtftk/stats/intersect/modl/dict_learning.py

    :param flags_matrix: Matrix of overlaps
    :param multiple_overlap_target_combi_size: Candidates longer than this will be discarded
    :param multiple_overlap_max_number_of_combinations: Final desired number of combinations
    :param nb_threads: Number of threads
    :param step_1_factor_allowance: In step 1 of building the candidates, how many words are allowed in the Dictionary Learning as a proportion of multiple_overlap_max_number_of_combinations
    :param error_function: error function used in step 2. Default to manhattan error. 
    :param smother: Should the smothering which reduces each row's abudane to its square root to emphasize rarer combinations be applied ? Default is True
    :param normalize_words: Normalize the words by their summed squares in step 2. Default True.
    :param step_2_alpha: Override the alpha used in step 2.
    :param discretization_threshold: discretization_threshold in step 1. In each atom, elements below D*maximum_for_this_atom will be discarded. Optional.
    :param step_1_alphas: A list to manually override the alphas to be used during step 1. Optional.

    Passing a custom error function, it must have the signature error_function(X_true, X_rebuilt, encoded, dictionary). 
    X_true is the real data, X_rebuilt is the reconstruction to evaluate, code is the encoded version (in our case used 
    to check sparsity), dictionary has one learned atom per row.
    
    >>> from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> flags_matrix = test_data_for_modl(nflags = 1000, number_of_sets = 6, noise = 0.05, cor_groups = [(0,1),(0,1,2,3),(4,5)])
    >>> combi_miner = Modl(flags_matrix, multiple_overlap_max_number_of_combinations = 3)
    >>> interesting_combis = combi_miner.find_interesting_combinations()
    >>> assert set(interesting_combis) == set([(1,1,0,0,0,0),(1,1,1,1,0,0),(0,0,0,0,1,1)])
    
    """

    def __init__(self, flags_matrix,
                 multiple_overlap_target_combi_size = -1,
                 multiple_overlap_max_number_of_combinations = 5,
                 nb_threads = 1,
                 step_1_factor_allowance = 2, 
                 error_function = None,
                 smother = True,
                 normalize_words = True,
                 step_2_alpha = None,
                 discretization_threshold = 0,
                 step_1_alphas = None):

        # Matrix of overlap flags to work with
        self.original_data = flags_matrix

        ## Squishing
        # NOTE Redundancy in lines is not useful and only increases computing time for
        # for a linear change in all errors. So to save time while not losing much information
        # we use a smothered matrix where each line has its original abundance divided 
        # by abundance of most rare (but not lower than abundance_threshold, 1/10000 by default)
        # eg if X = [A * 1000, B * 10], X' = [A * 100, B * 1]
        # Then, to diminish the emphasis on frequent lines, we use in the end
        # the square root of that abundance.
        self.smother = smother  # Use smothering ?
        self.data = squish_matrix(self.original_data,
            abundance_threshold = 1E-4, smother = self.smother)
        # NOTE The original data is kept under self.original_data but is NOT currently used.



        ## Parameters
        # Those are respectively the desired number of signifcant combis, 
        # and the maximum size of the combis (for the filter library step)
        # They are both command line parameters in ologram
        if multiple_overlap_max_number_of_combinations <= 0 : raise ValueError("multiple_overlap_max_number_of_combinations must be greater than 0")
        self.queried_words_nb = multiple_overlap_max_number_of_combinations
        self.max_word_length = multiple_overlap_target_combi_size

        # Step 1
        self.step_1_factor_allowance = step_1_factor_allowance   # In step 1, how many words to ask for each time, as a proportion of multiple_overlap_max_number_of_combinations
        self.step_1_alphas = step_1_alphas                       # Override the list of alphas in step 1?
        self.discretization_threshold = discretization_threshold # When discretizing the words in step 1, what threshold?

        # Step 2
        self.error_function = error_function    # Remember the error function for step 2
        self.normalize_words = normalize_words  # Do we normalize the words by their summed squared in step 2?        
        self.step_2_alpha = step_2_alpha        # Override alpha in step 2?

        self.nb_threads = nb_threads



    # ------------------------ Elementary subroutines ------------------------ #

    def generate_candidate_words(self):
        """
        Generate the library by running dictionary-learning based matrix factorizations 
        on the data matrix
        """

        ## Generate subsamples from self.data
        # Normally 3 (or less if there are less than 3 lines in the matrix)
        N = min(3, self.data.shape[0])
        # Prepare N subsamples each time, to increase learned combi diversity
        # This is done like a cross validation, to ensure no line is forgotten      
        kf = KFold(n_splits=N, shuffle=False)
        subsamples = []
        for indexes, _ in kf.split(self.data):
            subsamples += [self.data[indexes,:]]


        # ---- Candidate words generation
        # Encode with and a variety of alphas

        # In the encoding, to get some diversity, request k* the final queried
        # number of combinations.
        # k defaults to 2 in the initiation.
        n_words_step_one = self.step_1_factor_allowance * self.queried_words_nb

        # Do this for all the subsamples
        self.all_found_words = dict()

        for sub in subsamples:
            words_this_round = modl_subroutines.generate_candidate_words(sub,
                n_words = n_words_step_one, nb_threads = self.nb_threads,
                discretization_threshold = self.discretization_threshold,
                alphas = self.step_1_alphas)

            # Merge the {word:usage} dicts that were found each time by
            # taking the SUM of usages where relevant
            x = self.all_found_words
            y = words_this_round

            self.all_found_words = {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}


    def filter_library(self):

        ## Sort the words by usage (ie. sort the dictionary keys by their values)
        # to produce the final list
        sorted_kv_by_value = sorted(self.all_found_words.items(),
            key=lambda kv: kv[1],   # Sort by usage
            reverse = True) # Get most used words first
        all_candidate_words_ordered = [kv[0] for kv in sorted_kv_by_value]

        ## Remove the words longer than the user wants (with sum higher than multiple_overlap_target_combi_size)
        # Default is -1, meaning no filtering should be applied and all words should be kept
        if self.max_word_length == -1: self.max_word_length = np.inf
        final_words = [tuple(word) for word in all_candidate_words_ordered if sum(word) <= self.max_word_length]


        ## Pre-selection based on usage.
        # To save time on step 2 and not have to perform quite as many sparse encoding,
        # only send the top N candidates from step 1, sorted by total usage across
        # the learned reconstructions in step 1
        # I currently keep 3* queried
        final_words = final_words[:3*self.queried_words_nb]

        message('Keeping only top words from step 1 by usage. Final list of words is '+str(final_words))


        # Finally, record the words : do all the operations to create a Library
        self.library = tree.Library()
        self.library.build_nodes_for_words(final_words)
        self.library.assign_nodes()

        # Remember the number of words in the library for a stop condition later
        self.number_of_words_in_library = len(final_words)

        

    def select_best_words_from_library(self):
        """
        This is step 2. Takes the library of candidates produced at step 1
        and will get the best N words among it that best rebuild the original matrix.
        """

        # You can't request more words than are actually present in the library, 
        # nor than unique elements in the data
        upper_floor_words = min(self.number_of_words_in_library,
            len(np.unique(self.data, axis = 0)))
        if self.queried_words_nb > upper_floor_words: 
            self.queried_words_nb = upper_floor_words
            message("Requesting too many words, reducing to "+str(self.queried_words_nb)) 
        # NOTE It will actually be +1 to make room for the root (0,0,0,...) word, but this is added later
        
        # Read the parameters that were supplied when creating the Modl object
        best_dict = modl_subroutines.build_best_dict_from_library(
            self.data, self.library,    # Data and Library of candidates
            self.queried_words_nb,      # N best words
            self.error_function,        # Potential custom error function
            self.nb_threads,
            self.normalize_words,       # Normalize words by sum of square
            self.step_2_alpha)          # Sparsity control
        
        # Final step : register the best dictionary
        self.best_words = best_dict




    # --- Main function --- #
    def find_interesting_combinations(self):
        """
        MAIN FUNCTION. Will call the others.
        """

        ## Hardcode ignoring of Python (and by extension SKLearn) warnings

        try:
            previous_warning_level = os.environ["PYTHONWARNINGS"]
        except:
            previous_warning_level = 'default'
        
        if utils.VERBOSITY < 2: # Only if not debugging
            os.environ["PYTHONWARNINGS"] = "ignore"
            #warnings.filterwarnings('ignore', module='^{}\.'.format(re.escape("sklearn")))
            message("Filtering out sklearn warnings.")

        ## Now call the functions
        self.generate_candidate_words()
        self.filter_library()
        self.select_best_words_from_library()

        # Re-enable warnings
        os.environ["PYTHONWARNINGS"] = previous_warning_level

        return self.best_words
