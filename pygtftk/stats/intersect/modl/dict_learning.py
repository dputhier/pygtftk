"""

This module contains the MODL algorithm, an itemset mining algorithm that broadly consists of two steps.

Considering an input matrix with transactions/regions.intersections in lines and elements/items/sets in columns,
the algorithm will :
    - learn (MiniBatchDictionaryLearning) several factorisations with different sparsity constraints to build a library of candidate words (eg. (1,1,0,0,1))
    - select the best words to rebuild the original matrix using a greedy algorithm.

Unitary tests present examples of application below.

Describe dict learning Here (or rather in the doc ?)


In OLOGRAM, we use this on the matrix of true overlap flags to find usual common overlaps.


If you instead want to consider exclusive combis, is is equivalent to simly takin the most common rows (will be useful later with intra-sets overlaps.)

EXPLAIN IN DOCS and/or NOTES THAT THIS IS MOSTLY USEFUL IF THERE ARE MANY FILES (dozen+) and/or intraset overlap, and you
specify a custom combi size


Dictionary learning is an optimization problem solved by alternatively updating the sparse code, as a solution to multiple Lasso problems, considering the dictionary fixed, and then updating the dictionary to best fit the sparse code.

MUST SPECIFy that dict learning algo and sparse encoder algo are mostly independant


ADD DETAILS ABOUT MY NEW CUSTOM ALGO :
Goal : find combinations of interest (itemset mininig)
to "compress" the matrix of all the true "overlap flags" to find interesting combinations of custom size (words of the dictionary).
in summary we make various reconstructions at various sparsity constraints to get a set of candidate words of various lengths (!). Then build the final combis from the candidates by getting the best encoding ones.



TODO ADD A BLEEDING UNITARY TEST WITH SIMPLE_02 DATA .1 and .2 THAT I CREATED FOR THIS !


Author : Quentin FERRE <quentin.q.ferre@gmail.com>

This code is available under the GNU GPL license.
"""

import copy, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

# Library tree custom code
from pygtftk.stats.intersect.modl import tree   

# THOSE IMPORTS MIGHT NOT BE NEEDED HERE
from sklearn.decomposition import SparseCoder
from sklearn.decomposition import MiniBatchDictionaryLearning




# TODO USE THIS FOR MESSAGES !!!
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



def squish_matrix(x, abundance_threshold = 0, shuffle = True):
    """
    To reduce redundancy, take all uniques rows and create a matrix where each line has its original
        abundance divided by abundance of most rare (but not lower than abundance_threshold, 1/10000 by default)
            eg if X = [A * 1000, B * 10], X' = [A * 100, B * 1]
            TODO ADD THIS CASE ABOVE AS AN UNITARY TEST

    >>> import numpy as np
    >>> X = np.array([[1,1,0,0]]*1000 + [[0,0,1,1]]*100)
    >>> X_squished = squish_matrix(X, shuffle = False)
    >>> np.testing.assert_equal(X_squished, np.array([[0,0,1,1]]*1 + [[1,1,0,0]]*10)) # Note that the rows have been sorted by abundance   
    """

    # Get all unique rows and their counts
    all_rows, counts_per_row = np.unique(x, axis=0, return_counts = True)


    # Divide counts by lowest count observed (but never divide by under abundance_threshold)
    minimal_count = max(np.min(counts_per_row), abundance_threshold*x.shape[0])
    counts_relative = np.around(counts_per_row / minimal_count).astype(int)

    # Rebuild a new matrix with those counts
    squished_matrix = []     
    for row, count in zip(all_rows, counts_relative):
        squished_matrix += [row] * count

    final_matrix = np.array(squished_matrix)
    if shuffle: np.random.shuffle(final_matrix) # Shuffle it for good measure

    return final_matrix





# TODO May need to not binarize them to be more accurate later (ie. keep ratios ?)
# Not right now, but mark is as a NOTE for improvement
def binarize_words(V):
    return






# ---------------------------------------------------------------------------- #
#                             Tree model for the library                       #
# ---------------------------------------------------------------------------- #
from pygtftk.stats.intersect.modl import tree as tree




# ---------------------------------------------------------------------------- #
#                             SUBROUTINES OF MAIN                              #
# ---------------------------------------------------------------------------- #
from pygtftk.stats.intersect.modl import subroutines as modl_subroutines




# TODO REACTIVATE THOSE AS UNITARY TESTS MAYBE




"""
# TODO Move to subroutines as unitary test ?
if __name__ == "__main__":

    # How big can we go ?
    Xbig = test_data_for_modl(nflags = 1000, number_of_sets = 20,
        cor_groups = [(0,1),(2,3),(4,5),(0,1,2),(3,4,5),(10,11,12),(13,14),(15,16)])

    start_time = time.time()
    U, V, error = modl_subroutines.learn_dictionary_and_encode(Xbig, n_atoms = 7, alpha = 1,
                                            n_jobs = 8, n_iter = 1000)
    stop_time = time.time()
    #print('One DL step took ', str(stop_time-start_time), 'seconds.')

    # NOTE Time scaling with number of lines and columns is only moderate. Even on 1 core, 
    # 1 million (!) lines for 20 sets it was 58 seconds ! 8 cores halve this time.
    






    ## Art of the compromise, only if compromise itself is present in the data !!
    Xt = np.array([[1,0,0],[0,0,1]]*100 + [[1,0,1]])
    #Xt = np.array([[1,0,0],[0,0,1],[0,1,1]]*10)
    r = modl_subroutines.learn_dictionary_and_encode(Xt, n_atoms = 1, alpha = 1E-100)
    #plt.figure(); sns.heatmap(r[1]**2, annot = True)








# KEEP THIS AS UNITARY TEST !!!
if __name__ == "__main__":


    fl = [[1,1,0,0,0,0]]*2 + [[1,1,0,0,1,1]] + [[0,0,0,0,1,1]]*20
    X = np.array(fl*4)
    r = modl_subroutines.generate_candidate_words(X, n_words = 3, nb_threads = 1)



  







    #Not working ! the algo does not want to elarn the (1,1,1,1) word !
 
    fl = [[1,1,1,1,1,1]]*20
    X = np.array(fl)
    r = modl_subroutines.generate_candidate_words(X, n_words = 3, nb_threads = 1)




# ----- STEP 2 OF THE ALGO
# --- HERE KEEP THIS AS UNITARY TEST BELOW
if __name__ == "__main__":
    all_candidate_words = [(1,1,0,0,0,0),(1,1,0,0,1,1),(0,0,0,0,1,1),
                            (1,1,1,1,0,0)] # Useless one

    l = tree.Library()
    l.build_nodes_for_words(all_candidate_words)
    l.assign_nodes()

    fl = [[1,1,0,0,0,0]]*2 + [[1,1,0,0,1,1]] + [[0,0,0,0,1,1]]*2
    X = np.array(fl*4)



    best_dict = modl_subroutines.build_best_dict_from_library(X, l, queried_words_nb = 2)

    assert best_dict == [(0, 0, 0, 0, 1, 1), (1, 1, 0, 0, 0, 0)]





























# Now testing with noisy data.

#keep this code : it will be used in the supplementary/Snakemake github to compare it to apriori

if __name__ == "__main__":
    np.random.seed(42)
    Xmt = test_data_for_modl(nflags = 20000, number_of_sets = 20,
        noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])


    # Manually specify candidates
    all_candidate_words = [
        (1,1,1,1,0,0),(1,1,0,0,0,0),(0,0,0,0,1,1),   # Correct
        (0,0,1,1,0,0),(0,0,1,1,1,1),(1,1,0,0,1,1),   # Incorrect - not observed
        (1,0,0,0,0,0),(0,1,0,0,0,0),(0,0,1,0,0,0)    # Incorrect - Too small
    ]
    l = tree.Library()
    l.build_nodes_for_words(all_candidate_words)
    l.assign_nodes()


    # Try with curated test
    Xmt = np.array([
        [1,1,0,0,0,0],[1,1,0,0,0,0],[1,1,1,1,0,0],[0,0,0,0,1,1]
    ]*10)


    best_dict = modl_subroutines.build_best_dict_from_library(Xmt, l, queried_words_nb = 3)






# Scaling test
if __name__ == "__main__":

    NB_SETS = 6


    Xmt = test_data_for_modl(nflags = 20000, number_of_sets = NB_SETS,
        noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])


    # Generate all possible words as candidates, extreme scaling testing
    import itertools
    all_candidate_words = list(itertools.product([0, 1], repeat=NB_SETS))



    all_candidate_words = [(1,0,0,0,0,0),(1,1,0,0,0,0), (0,1,0,0,0,0)]
    #all_candidate_words = [(0,0,0,0,0,1),(0,0,0,0,1,1), (0,0,0,0,1,0)]
    

    import time 

    start = time.time()
    l = tree.Library()
    l.build_nodes_for_words(all_candidate_words)
    l.assign_nodes()
    stop = time.time()
    print("Took", stop-start)

    #2,5,15,40,64,95,172


    

    # Try with curated test
    start = time.time()
    best_dict = modl_subroutines.build_best_dict_from_library(Xmt, l, queried_words_nb = 1, nb_threads = 8)
    stop = time.time()
    print("Took", stop-start)




"""










# ---------------------------------------------------------------------------- #
#                                MAIN                                          #
# ---------------------------------------------------------------------------- #


import numpy as np

class Modl:
    """
    This class encapsulates the MODL approach :

    Takes as input a matrix of flags, with one flag per intersection,
    and returns the list of interesting combis.


    ADD MUCHO EXPLANATIONS HERE
    """

    def __init__(self, flags_matrix,
                 multiple_overlap_target_combi_size,
                 multiple_overlap_max_number_of_combinations,
                 nb_threads):

        # Matrix of overlap flags to work with
        self.original_data = flags_matrix




        # UPDATE
        # NOTE Redundancy in lines is not useful and only increases computing time for
        # for a linear change in all error. So to save time while not losing much information
        # we use a squished matrix where each line has its original
        # abundance divided by abundance of most rare (but not lower than abundance_threshold, 1/10000 by default)
        # eg if X = [A * 1000, B * 10], X' = [A * 100, B * 1]
        self.data = squish_matrix(self.original_data,
            abundance_threshold = 1E-4)

        # NOTE The original data is kept under self.original_data but is NOT currently used.




        # DEBUG TODO REMOVE THIS
        print("DATA AFTER SQUISHING - UNIQUE AND NUMBER OF ROWS")
        all_combis, counts_per_combi = np.unique(flags_matrix, axis=0, return_counts = True)
        print(all_combis, counts_per_combi)





        # Parameters
        if multiple_overlap_max_number_of_combinations <= 0 : raise ValueError("multiple_overlap_max_number_of_combinations must be greater than 0")
        self.queried_words_nb = multiple_overlap_max_number_of_combinations



        self.max_word_length = multiple_overlap_target_combi_size
        """
        Those are respectively the desired number of signifcant combis, and the maximum size of the combis (for the filter library step); Explain this !
        THOSE WILL BOTH BE COMMAND LINE PARAMETERS
        """




        self.nb_threads = nb_threads
        """
        TODO THIS IS TO BE USED BY SEVERAL SKLEARN FUNCTIONS, MAKE SURE IT IS CORRECTLY DEFINED AND PASSED !
        Default should be equal to the number of threads defined when calling OLOGRAM in command line !!!!!!
        """


    # ------------------------ Elementary subroutines ------------------------ #


    def generate_candidate_words(self):

        """
        Generate the library by running dictionary-learning based matrix factorizations 
        on the data matrix
        """

        """
        If the data matrix is too big, it will be long to factorize. So
        instead we look for interesting combinations in random subsamples of
        the matrix then pool all found combinatinos together.
        """


        ## Generate subsamples from self.data
        N = 3
        # Prepare N subsamples each time, to increase learned combi diversity
        # This is done like a cross validation, to ensure no line is forgotten      
        from sklearn.model_selection import KFold
        kf = KFold(n_splits=N, shuffle=False, random_state=42)

        subsamples = []
        for indexes, _ in kf.split(self.data):
            subsamples += [self.data[indexes,:]]


        # ---- Candidate words generation
        # Encode with Lasso-LARS and a variety of alphas


        # How many words each time ? The query number of words ?
        # I think I wanted to start with half of all unique binarized,
        # NEW IDEA : ask the DL iterations of step 1 to learn twice the number of desired words
        # in the end after step 2? It's an idea. Let us try it.
        n_words_step_one = 2 * self.queried_words_nb

        # TODO MAKE IT A PARAMETER OF THE FUNCTION !!!!!







        # TODO NOTE IN PAPER ALSO Remember that if you allow more words than
        # necessary in the DL step, it simply will not use them.
        # TODO ACTUALLY I have decided to prevent that because it cna cause errosr, if you ask for more
        # words it will be reduced to the number of unique combis

        # TODO SAY IN PAPER random seed can make a word of differenec when not overcomplete ! Different words may be found ! THIS IS one more argument for random sampling !



        # Do this for all the subsamples !!!!!!!!,


        # TODO As has been established we might keep the subsamplings ?

        # Do subsampling
        self.all_found_words = dict()

        for sub in subsamples:
            words_this_round = modl_subroutines.generate_candidate_words(sub,
                n_words = n_words_step_one, nb_threads = self.nb_threads)


            # HOLD ON. I NEED TO MERGE LIBRARIES HERE NO ? MAYBE INSTEAD I COULD RETURN THE WORDS AND THEN ONLY BUILD A TREE AFTER ALL THE SUBSAMPLINGS !
            # Yep. So for now discard the individuals libraries, merge the lists of words, and


            # Merge the {word:usage} dicts that were found each time by
            # taking the SUM of usages where relevant
            """
            # TODO SAY SO IN PAPER !!!!!
            """
            #self.all_found_words += words_this_round

            x = self.all_found_words
            y = words_this_round

            self.all_found_words = {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}


            # DEBUG
            print("ALL CANDIDATE WORDS WITH USAGE NOW:", self.all_found_words)












    def filter_library(self, max_nb_sets_used_in_words):



        ## Sort the words by usage (ie. sort the dictionary keys by their values)
        # to produce the final list
        sorted_kv_by_value = sorted(self.all_found_words.items(),
            key=lambda kv: kv[1],   # Sort by usage
            reverse = True) # Get most used words first
        all_candidate_words_ordered = [kv[0] for kv in sorted_kv_by_value]






        # Remove the words longer than the user wants (with sum higher than multiple_overlap_target_combi_size)
        # DEFAULT is -1, meaning no filtering should be applied and all words should be kept
        if max_nb_sets_used_in_words == -1: max_nb_sets_used_in_words = np.inf

        final_words = [tuple(word) for word in all_candidate_words_ordered if sum(word) <= max_nb_sets_used_in_words]




        # TODO NOTE IN PAPER
        # New step : pre-selection based on usage.
        # To save time on step 2 and not have to perform quite as many sparse encoding,
        # only send the top N candidates from step 1, sorted by total usage across the learned reconstructions in step 1
        
        # TODO I currently keep 3* queried, note it in paper and play with this parameter !!!! MAKE IT A PARAMETER !
        # TODO MAKE IT A PARAMETER OF THE FUNCTION !!!!!
        final_words = final_words[:3*self.queried_words_nb]

        print('Keeping only top words from step 1 by usage. Final list of words is', str(final_words))




        # Finally, record the words : do all the operations to create a Library
        self.library = tree.Library()
        self.library.build_nodes_for_words(final_words)
        self.library.assign_nodes()

        # Remember the number of words in the library for a stop condition later
        self.number_of_words_in_library = len(final_words)

        

    def select_best_words_from_library(self, error_function = None):
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
            message("Requesting too many words, reducing to"+str(self.queried_words_nb)) 
        # NOTE It will actually be+1 to make room for the root (0,0,0,...) word, but this is added later
        

        # Read the parameters that were supplied when creating the Modl object
        best_dict = modl_subroutines.build_best_dict_from_library(
            self.data, self.library,    # Data and Library of candidates
            self.queried_words_nb,      # N best words
            error_function,             # Potential custom error function
            self.nb_threads)
        
        
        # Final step : register the best dictionary
        self.best_words = best_dict









    # --- Main function --- #
    def find_interesting_combinations(self):
        """
        MAIN FUNCTION. Will call the others.
        """





        # TODO Check data and all parameters properly initialized


        """
        Shortly remind here of the meaning of each step
        """
        self.generate_candidate_words()
        self.filter_library(self.max_word_length)

        self.select_best_words_from_library()



        return self.best_words





# Unitary test 
if __name__ == '__main__':


    # Generate some test data
    np.random.seed(42)

    flags_matrix = test_data_for_modl(nflags = 1000, number_of_sets = 6,
        noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])


    start_mining = time.time()
    combi_miner = Modl(flags_matrix,
        multiple_overlap_target_combi_size = -1,
        multiple_overlap_max_number_of_combinations = 3,
        nb_threads = 6)
    
    interesting_combis = combi_miner.find_interesting_combinations()
    stop_mining = time.time()

    print ("----------------")
    print("MODL mined interesting combis in",str(stop_mining-start_mining),"seconds.")    






"""
TODO : MAKE SURE multiple_overlap_target_combi_size AND multiple_overlap_max_number_of_combinations
CAN BE PASSED '-1', DEFAULT TO '-1', AND WILL INDEED HAVE NO MAXIMUM WHEN PASSED '-1'
"""







plt.close('all')