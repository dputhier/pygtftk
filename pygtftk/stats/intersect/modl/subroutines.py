"""
Pure subroutines of the MODL algorithm.

Author : Quentin Ferré <quentin.q.ferre@gmail.com>
"""

import time

from collections import defaultdict

import numpy as np 
import pandas as pd

from sklearn.decomposition import SparseCoder
from sklearn.decomposition import MiniBatchDictionaryLearning


from pygtftk.stats.intersect.modl import tree

from pygtftk.utils import message



def learn_dictionary_and_encode(data, n_atoms = 20, alpha = 0.5,
        n_iter = 200, random_seed = 42, n_jobs = 1,
        fit_algorithm = 'cd',
        transform_algorithm = 'lasso_cd'):
    r"""
    Will learn a dictionary for the data (row-wise) and encode it, with the specified parameters.
    Returns the dictionary, components as a dataframe, and encoded data.

    By default, allows 20 words with alpha of 0.5.

    More info at : https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.MiniBatchDictionaryLearning.html

    >>> from pygtftk.stats.intersect.modl.dict_learning import test_data_for_modl
    >>> from pygtftk.stats.intersect.modl.subroutines import learn_dictionary_and_encode
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> flags_matrix = test_data_for_modl(nflags = 1000, number_of_sets = 6)
    >>> import time
    >>> start_time = time.time()
    >>> U_df, V_df, error = learn_dictionary_and_encode(flags_matrix, n_atoms = 20, alpha = 0.5)
    >>> stop_time = time.time()

    """


    # TODO: Many operations here, such as recreating an object of a pandas 
    # dataframe, might not be necessary ?
    
    data = np.array(data) # Force cast as array

    # Cannot ask for more rows than there are unique lines
    n_atoms = min(n_atoms, len(np.unique(data, axis = 0)))


    dico = MiniBatchDictionaryLearning(n_components=n_atoms, n_iter=n_iter,
        alpha = alpha,
        fit_algorithm = fit_algorithm, 
        transform_algorithm = transform_algorithm, transform_alpha = alpha, 
        positive_dict = True, positive_code = True,
        random_state = random_seed,
        n_jobs = n_jobs 
        ) 

    # NOTE fit_algorithm is used during the learning and transform_algorithm 
    # transforms the data once the estimator has been fitted.
    # We are using coordinate descent (CD) as LARS has troubles with correlated
    # features. Also because we want to be able to enforce positivity in the 
    # dictionary and the code for interpretability.

    dico.fit(data)  # Fit the data

    # Get components (the dictionary).
    # NOTE: Use a try-except for future proofing, as sklearn (v 0.24) seems to 
    # deprecate 'components_' across the board and replaced it with 'dictionary'
    try: V = dico.components_
    except: V = dico.dictionary  
    V_df = pd.DataFrame(V)


    # If the alpha is inadapted, the learned dictionary may have failed to converge
    # and contain NaNs. If that happens, return it now and bypass the next steps
    # NOTE LARS is less vulnerable to it, use it as a fallback before calling it quits.
    if np.isnan(V).any():
        message("Learned dictionary by Coordinate descent contained NaNs. Alpha may have been indadapted. Defaulting to LARS.", type = 'DEBUG')
        
        dico = MiniBatchDictionaryLearning(n_components=n_atoms, n_iter=n_iter,
            alpha = alpha, transform_alpha = alpha, random_state = random_seed, n_jobs = n_jobs,
            fit_algorithm = 'lars',  transform_algorithm = 'lasso_lars') 

        dico.fit(data)
        try: V = dico.components_
        except: V = dico.dictionary  
        V_df = pd.DataFrame(V)
                
        # Only abort if it still does not work
        if np.isnan(V).any():
            message("Fallback still contains NaNs. Aborting.", type = 'DEBUG')
            return None, V_df, None # Return "U", V, "error"


    # Re-encode the data with this dictionary
    encoded_data = dico.transform(data)
    ed_df = pd.DataFrame(encoded_data)

    # Compute associated normalized L2 loss for reference
    reconstructed_features = np.matmul(encoded_data,V)
    error = np.sum((reconstructed_features-data)**2) / np.sum(data**2)

    return ed_df, V_df, error # Return U, V, error






# ---------------------------------------------------------------------------- #
#                        Step 1 of the MODL algorithm                          #
# ---------------------------------------------------------------------------- #

def generate_candidate_words(X, n_words, nb_threads = 1,
    discretization_threshold=0, alphas=None):
    """
    Given a data matrix, will generate a library of candidate words by trying many
    different Dictionary Learnings (different alphas only for now)
    """

    EPSILON = 1E-10

    ## Preliminary steps
    X = np.array(X) # Force cast to array for safety
    nb_features = X.shape[1] # word size

    # NOTE : n_atoms cannot be higher than the number of unique combis or there might be an Exception thrown.
    # It seems to vary and be due to alpha. To be safe, prevent it.
    nb_unique_combis = len(np.unique(X, axis = 0))
    if n_words > nb_unique_combis : n_words = nb_unique_combis

   
    # A too low alpha can result in learning giberrish words, thay may still be used as there is no penalty for it.
    # To counter this, we begin with an alpha of 1/nb_features
    # Rk: if alphas have been manually specified, use them instead
    if alphas is None: alpha = 1/nb_features
    else: alpha = alphas[0]

    # Remember all candidate words found at each step, and use a dictionary so we
    # can remember their total usage for later filtering
    all_candidate_words = defaultdict(lambda: 0)

    stop = False

    misbehaving = 0

    iternb = 0
    while not stop:

        # Encode using DL     
        #start_time = time.time()
        message("> Alpha is "+str(alpha), type = 'DEBUG')

        # NOTE For some reason, the algorithm will never learn the word (1,1,1,...,1,1) even if it the only row in the data.
        # To prevent that, add a useless padding for the last column.
        pad = np.zeros((X.shape[0],1))
        Xpadded = np.append(X, pad, axis=1)

        # Perform the learning
        U, V, error = learn_dictionary_and_encode(Xpadded, n_atoms = n_words, alpha = alpha,
                                               n_jobs = nb_threads)

        #stop_time = time.time() # Careful not to name it 'stop' and overwrite the stop flag :)
        #print('One DL step took ', str(stop_time-start_time), 'seconds.')

        iternb += 1
        # Update alpha so DL will be more sparse in the next step
        # and look for higher-order combinations
        # We add iternb/nb_features to not take too long and not unduly favor high alphas (and longer words)
        
        # Rk: if alphas have been manually specified, use them instead
        if alphas is None: alpha += iternb/nb_features
        else:
            try: alpha = alphas[iternb]
            except: stop = True # Stop if we reached the end of the list
        
        
        # TODO: Also stop regardless after 2*k iterations ?

        # We stop once alpha is too high and preventing any word from being used
        # in the encoding U (which contained NaNs -- handled before -- or has a sum of 0)
        # Rq : can also happen with too low alpha or convergence errors, so stop only if it happens more than 2 times
        if (U is None): misbehaving += 1
        if (misbehaving >=2): stop = True
            
        if U is not None:

            if np.sum(U.values) == 0:
                stop = True
                continue # Do not do the following steps of word extraction !! The stop signal has been sent, and if U is None there is nothing to extract



            # -- Extract interesting words -- #
            # Here extract interesting words from dico and add to all_found_words

            # Binarize them (ie transform 0.9 - 0.9 - 0.001 into 1-1-0)
            # The sum of V_df**2 for each word is always 1.0 so our cutoff to say if a feature is used
            # is V**2 > 1/n_features**2 - epsilon for rounding errors. 
            # We cannot simply binarize, since 0.001 would be binarized to 1 yet it is likely not interesting.
            # WARNING : no longer true if I have intra set overlap (ie 1-1-2 for example) so when looking for the combis I might need to ignore intra set overlap

            # Each learned word is remembered along with its usage in the reconstruction
            # This way when proceeding to step 2 we keep only the most relevant words from step 1 to save time.
            # Gibberish overcomplete words that are found when n_atoms is too high and were not used will have usage of 0
            for i in U.columns.tolist():
                this_word = (V**2).iloc[i,:].values

                # Remove padding
                this_word = this_word[:-1]

                # Above 1/k**2 ?
                first_thres = this_word - (1/nb_features**2) > EPSILON

                # Above discretization_threshold * max ?
                second_thres = this_word - (discretization_threshold * max(this_word)) > EPSILON

                # Above both ?
                this_word_binarized = np.logical_and(first_thres,second_thres).astype(int)


                this_word_binarized = tuple(this_word_binarized.tolist())

                this_word_usage = np.sum(U.values[:,i])

                # Now remember it
                all_candidate_words[this_word_binarized] += this_word_usage

                # Debug print
                message("Word = "+str(this_word_binarized)+" ; Usage = "+str(this_word_usage),
                    type = 'DEBUG')



    # Return the dictionary giving all candidate words and their total usage
    return all_candidate_words




# ---------------------------------------------------------------------------- #
#                        Step 2 of the MODL algorithm                          #
# ---------------------------------------------------------------------------- #

def normalize_and_jitter_matrix_rows(X, normalize = True, jitter = True):
    r"""
    Apply normalization so that sum of suquare of each row is 1,
    Add a slight epsilon so two rows may not have exact same dot product with a given vector
    Sort in lexicographic order

    >>> import numpy as np
    >>> import numpy.testing as npt
    >>> from pygtftk.stats.intersect.modl.subroutines import normalize_and_jitter_matrix_rows
    >>> D = np.array([[1,1,0],[0,1,0],[0,1,1]])
    >>> Dc = normalize_and_jitter_matrix_rows(D)
    >>> Dc_theory = [[0, 1, 0], [7.0736E-5, 7.07107E-1, 7.07107E-1], [7.07107E-1, 7.07107E-1, 9.99859E-05]]
    >>> npt.assert_almost_equal(Dc, Dc_theory, decimal = 5)
    """

    X = np.array(X) # Enforce NumPy array
    X = X[np.lexsort(np.rot90(X))] # Sort the dictionary in lexicographic order BEFORE applying all of this

    ## Applying corrections to the words : normalize and jitter
    new_X = []
    i = 0
    for row in X: 
        
        # To prevent atoms in D from having the exact same dot product with any word
        # that LARS will try to rebuild, add a very small jitter to the words 
        # -> sqrt(i) times 1E-4 to each element of the i-th row of D
        neo = np.array([x + np.sqrt(i)*1E-4 for x in row])
        i += 1        

        # Normalize the words by the sum of their square. Helps counter a tendency to select longer words, which can increase sensitivity to noise.
        if normalize:
            squared_sum = np.sum(neo**2)
            if squared_sum > 0: # Avoid division by zero
                neo = np.sqrt(neo**2/squared_sum)

        new_X += [neo]

    # Now convert to array
    X = np.array(new_X)
    
    return X



def build_best_dict_from_library(data, library, queried_words_nb,
                error_function = None,
                nb_threads = 1, normalize_words = False,
                transform_alpha = None):
    r"""
    Given a data matrix and a library, will select the best n = queried_words_nb words
    with a greedy algorithm from the library to rebuild the data.

    `data` is a matrix with one transaction per row and one element per column, as usual

    The greeedy algorithm works because the problem in this particular instance
    is submodular. In broad strokes, it is assumed that the gain in term of reconstruction
    with the Sparse Coder of adding the word X to a set S is no higher than adding X to a
    set T that contains S. Also, more words cannot make the sparse coder do a worse job so
    the function is monotonous.

    Instead of the reconstruction error, you may pass a different callable of the form
    error_function(data, rebuilt_data, encoded, dictionary) that returns an error value so there
    can be a supervision, but the submodularity might no longer hold.

    >>> import numpy as np
    >>> from pygtftk.stats.intersect.modl import tree
    >>> from pygtftk.stats.intersect.modl import subroutines as modl_subroutines
    >>> X = np.array([[1,1,0,0,0,0],[1,1,1,1,0,0],[0,0,0,0,1,1]]*10)
    >>> candidates = [(1,1,0,0,0,0),(1,1,1,1,0,0),(0,0,0,0,1,1), (1,1,0,0,1,1), (1,0,1,0,1,0),(1,1,0,1,1,0),(1,0,1,1,0,0),(0,1,0,0,1,1)]
    >>> L = tree.Library()
    >>> L.build_nodes_for_words(candidates)
    >>> L.assign_nodes()
    >>> selection = modl_subroutines.build_best_dict_from_library(X, L, queried_words_nb = 3)
    >>> assert set(selection) == set([(1,1,0,0,0,0),(1,1,1,1,0,0),(0,0,0,0,1,1)])

    """


    ## Create objects

    # Initial dictionary
    nb_features = data.shape[1]

    dictionary = [tuple([0] * nb_features)] # Must initialize with a word, might as well be full zeros.

    # Use nonzero alpha for the coder as well, only if not manually specified.
    if transform_alpha is None:
               
        # To encourage use of fewer words, we by default use a rather high alpha (but capped at 0.5)
        transform_alpha = min(1/np.sqrt(nb_features), 0.5)




    ## Coder object, to be updated iteratively
    coder = SparseCoder(dictionary=np.array(dictionary).astype('float64'), # Rq : must make it float64 for some reason
        transform_algorithm='lasso_lars', #transform_algorithm='lasso_cd',
        transform_alpha = transform_alpha,
        positive_code = True, # Very important to ensure we also force positivity here !
        n_jobs = nb_threads)

    # We use lasso_lars, as lasso_cd often had mysterious convergence problems.
    # Lasso_lars can stil enforce positivity, and with precomputed words now the poor handling of degenerate
    # vectors (lars dropping them) is a strength


    # --- Filtering the candidates
    # Iteratively add candidates that improve the reconstruction
    stop = False
    while stop is False:

        ## Generate all_candidate_words for this iteration
        # For now, simply all words of the library that are not already used
        all_candidate_words = tree.get_all_candidates_except(library, dictionary)

        # Test all words for this iteration
        candidates_with_errors = dict()

        message("----------------------", type = 'DEBUG')
        message("Current dictionary : "+str(dictionary), type = 'DEBUG')
        message("> Candidate\tError", type = 'DEBUG')

        for candidate in all_candidate_words:

            # Create a new dictionary for this test
            dict_being_tested = dictionary + [candidate]

            Dt = np.array(dict_being_tested)
            
            # Normalization helps counter a tendency to selcet longer words, and jitter fixes a LARS bug with atoms that have the same dot product
            # Also lexicographic sort.
            Dt_corrected = normalize_and_jitter_matrix_rows(Dt, normalize_words)
            Dt = Dt_corrected
            
            # Update dictionary to the one being currently tested
            coder.dictionary = Dt.astype('float64')

            # Just to be safe
            assert np.allclose(coder.dictionary, Dt.astype('float64'))


            try:

                encoded = coder.transform(data)
            
                rebuilt_data = np.matmul(encoded, Dt)

                # We use Manhattan error so compromises are discouraged : with L2,
                # when rebuilding [(1,0,0)*many,(1,0,1),(0,0,1)*many], (1,0,1) would be picked first
                # Add an alpha that is nonzero so that using longer words is still
                # an improvement even when other words cover it (and as a tiebreaker)
                # but keep it low (1/k) so you don't re-encourage compromise
                def manhattan_dist_with_sparsity(X_true, X_rebuilt, encoded, dictionary):

                    dictionary = None # This particular function ignores the dictionary

                    error_mat = np.abs(X_rebuilt-X_true)
                    error = np.sum(error_mat)

                    alpha = transform_alpha # Fetch the alpha used when calling the main function

                    regul = np.sum(encoded) * alpha

                    final_error = error + regul
                    return final_error


                # Check for error function, default is manhattan dist
                if error_function is None:
                    error_function = manhattan_dist_with_sparsity

                # To compute the error, pass respectively : true data, rebuilt data, encoded data, and current dictionary 
                # TODO: currently, the dictionary passed is BEFORE normalization and all that jazz 
                error = error_function(data, rebuilt_data, encoded, dict_being_tested)

            # On rare occasion, convergence errors can result in a ValueError
            except ValueError:
                error = np.inf

            # Remember error
            candidates_with_errors[tuple(candidate)] = error

            message(str(candidate)+'\t'+str(error), type = 'DEBUG')

      
        # Add best candidate to final dict and remove from candidates.
        try:
            # In case of a tie, the first word is selected.
            best_candidate = min(candidates_with_errors, key=candidates_with_errors.get)
            dictionary += [best_candidate]
        # It should never happen, but if there are no more candidates add an empty word.
        except:
            dictionary += [tuple([0] * nb_features)]


        """ 
        # NOTE for improvement
        # To save time, we may only pass as candidates the "children" (see tree module)
        # of nodes already present under the assumption that if X is better than Y
        # the best parent of X should be better than the best parent of Y ?
        # Then we would need to wait for an optimum, with this draft code

        # If we have more words than the user wants, remove the worst one
        if (len(dictionary) > queried_words_nb):
            # Remove worst word with respect to error (try all, opposite of above)

        # When the word added and the word removed are the same (meaning
        # the word it wants to add does not improve on the dict) stop and
        # return
        if best_candidate == last_word_removed :
            best_words = final_dict
            return best_words
        """

        ## Stop condition : as soon as we have enough words
        # Remark : since the dictionary always contain the (0,0,0,...) word as a supplement, we
        # must stop when it contains the required number of words... plus one !
        if (len(dictionary) >= queried_words_nb + 1):
            stop = True




    message("Final dictionary :"+str(dictionary), type = 'DEBUG')

    # Conclusion of the function : remove the (0,0,0,...) word and make words unique, just in case
    best_words = list(set(dictionary))
    best_words.remove(tuple([0] * nb_features))

    return best_words
