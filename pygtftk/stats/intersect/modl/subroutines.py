"""
Pure subroutines of the MODL algorithm

Author : Quentin Ferré <quentin.q.ferre@gmail.com>
"""


import time
import warnings

from collections import defaultdict

import numpy as np 
import pandas as pd

from sklearn.decomposition import SparseCoder
from sklearn.decomposition import MiniBatchDictionaryLearning
from sklearn.exceptions import ConvergenceWarning

from pygtftk.stats.intersect.modl import tree







# TODO USE THIS FOR MESSAGES !!!
from pygtftk.utils import message



def learn_dictionary_and_encode(data, n_atoms = 20, alpha = 0.5, n_iter = 100, random_seed = 42, n_jobs = 1):
    """
    Will learn a dictionary for the data (row-wise) and encode it, with the specified parameters.
    Returns the dictionary, components as a dataframe, and encoded data.

    By default, allows 20 words with alpha of 0.5.

    More info at : https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.MiniBatchDictionaryLearning.html


    UNITARY TEST : WITH 110000 * 5, 110011, 000011 * 5 --> find 110000 and 000011 as words, use proper alpha
    """


    # TODO : many operations here, suh as recreating an object of a pandas 
    # dataframe, might not be necessary ?

    
    data = np.array(data) # Force cast as array



    # TODO : move the below comments to the notes of this module, and to the paper
    dico = MiniBatchDictionaryLearning(n_components=n_atoms, n_iter=n_iter,
        # transform_n_nonzero_coefs = target_nb_words_per_ex, # Sadly only omp and lars use those, and we cannot use them as we enforce positivity
        alpha = alpha,
        fit_algorithm = 'cd', 
        transform_algorithm = 'lasso_cd', transform_alpha = alpha, # Must use this algo (in fit as well obviously) because we enforce positivity, cannot use omp or lars
        positive_dict = True, positive_code = True,
        random_state = random_seed,
        n_jobs = n_jobs # Multithreading natively
        ) 

    # NOTE fit_algorithm is used during the learning and transform_algorithm transforms the data once the estimator has been fitted
    # We are using coordinate descent (CD) as LARS has troubles with correlated features.
    # Positivity is enforced in the dictionary and the code for interpretability.

    

    dico.fit(data)  # Fit the data

    V = dico.components_ # Get components
    V_df = pd.DataFrame(V)




    # If the alpha is inadapted, the learned dictionary may have failed to converge
    # and contain NaNs. If that happens, return it now and bypass the next steps
    if np.isnan(V).any():
        print("Learned dictionary contained NaNs. Alpha was likeky inadapted. Aborting.")
        return None, V_df, None # Return "U", V, "error"


    # Re-encode the data with this dictionary
    encoded_data = dico.transform(data)
    ed_df = pd.DataFrame(encoded_data)







    # Compute associated normalized L2 loss for reference
    reconstructed_features = np.matmul(encoded_data,V)
    error = np.sum((reconstructed_features-data)**2) / np.sum(data**2)
    #print("Number of 'atoms' : "+str(n_atoms)+" \t ; Alpha = "+str(alpha)+" \t ; Normalized MSE : " + str(error))
    # TODO Don't return the error, calculate it manually when needed ? Or allow passing of an error_function callable ?

    return ed_df, V_df, error # Return U, V, error














# ----------- Step 1 of the MODL algorithm





def generate_candidate_words(X, n_words, nb_threads = 1):
    """
    Given a data matrix, will generate a library of candidate words by trying many
    different Dictionary Learnings (different alphas only for now)
    """

    EPSILON = 1E-10

    # Preliminary steps :
    #target n_words in the rebuilding should be 50% of unique BINARIZED combis by default

    X = np.array(X) # Force cast to array for safety

    nb_features = X.shape[1] # word size


    # NOTE : n_atoms cannot be higher than the number of unique combis or there might be an Exception thrown.
    # It seems to vary and be due to alpha. To be safe, prevent it.
    nb_unique_combis = len(np.unique(X, axis = 0))
    if n_words > nb_unique_combis : n_words = nb_unique_combis


    # print("n_words_in_subroutine",n_words)
    # print("UNIQUE",np.unique(X, axis = 0))
    
    



   
    # A too low alpha can result in learning giberrish word, thay may
    # still be used as there is no penalty for it.
    # To counter this, we begin with an alpha of 1/nb_features
    alpha = 1/nb_features


    # Remember all candidate words found at each step, and use a dictionary so we
    # can remember their total usage for later filtering
   
    all_candidate_words = defaultdict(lambda: 0)

    stop = False


    iternb = 0
    while not stop:

        # Encode using DL     
        start_time = time.time()
        print("> Alpha is", str(alpha))

        # Ignore convergence warnings from faulty alphas TODO explain more, it happens with the NaN that I catch later
        with warnings.catch_warnings():
            
            warnings.filterwarnings("ignore", category=ConvergenceWarning)


            # NOTE For some reason, the algorithm will never learn the word (1,1,1,...,1,1) even if it the only row in the data.
            # To prevent that, add a useless padding for the last column.
            pad = np.zeros((X.shape[0],1))
            Xpadded = np.append(X, pad, axis=1)

            # Perform the learning
            U, V, error = learn_dictionary_and_encode(Xpadded, n_atoms = n_words, alpha = alpha,
                                                n_jobs = nb_threads)



        stop_time = time.time() # Careful not to name it 'stop' and overwrite the stop flag :)
        #print('One DL step took ', str(stop_time-start_time), 'seconds.')
        






        """
        In certain cases, the learned dictionary can be full of NaNs and the transform step would return an exception.
        This usually means the alpha was inadapted.

        # Example where this happens
        Xt = np.array([[1,1,0],[0,0,1],[0,1,1],[0,1,0]]*10) ; r = modl_subroutines.learn_dictionary_and_encode(Xt, n_atoms = 3, alpha = 0.5)

        My current analysis is that it happens when alpha is TOO HIGH, so it is the signal that we should stop.
        """
        # We stop once alpha is too high and preventing any word from being used
        # in the encoding U
        if (U is None) or (np.sum(U.values) == 0) :
            stop = True
            continue # Do not do the following steps of word extraction !! The stop signal has been sent, and if U is None there is nothing to extract
        """
        To decide to stop, actually check if the sum(U) is indeed zero, meaning alpha is too high and no words are used.
        NOTE : It seems in that case the dictionary will often not have been learned and be full of zeroes.
        I remember observing both cases, so the check takes both into account.
        TODO Also automatically stop after, like, 25 alpha steps no matter what.
        """







        # -- Extract interesting words -- #
        # Here extract interesting words from dico and add to all_found_words

        # Do not simply take all words : if the n_atoms is too high you can have leftover gibberish for the words that were not used)
        # Word must have been used at least once in the rebuilding ie. whose coefficient is higher than zero/epsilon in at least one line
        #  (one flag vector) across the true data

        # Now that we know which words are used, we must binarize them (ie transform 0.9 - 0.9 - 0.001 into 1-1-0)
        # The sum of V_df**2 for each word is always 1.0 so the cutoff to say if a feature is used or not
        # is V**2 > 1/n_features - epsilon for rounding errors. We cannot simply binarize, since 0.001 would be binarized to 1 yet it is likely not interesting
        """
        # TODO WRITE IT IN PAPER !
        """
        # WARNING : no longer true if I have intra set overlap (ie 1-1-2 for example) so when looking for the combis I might need to ignore intra set overlap

        # NEW UPDATED VERSION
        # Each learned word is remembered along with its usage in the reconstruction
        # This way when proceeding to step 2 we keep only the most relevant words from step 1 to save time
        # Gibberish overcomplete words that were not used will have usage of 0
            #if the n_atoms is too high you can have leftover gibberish for the words that were not used)
        for i in U.columns.tolist():
            this_word = (V**2).iloc[i,:].values

            # Remove padding
            this_word = this_word[:-1]


            this_word_binarized = ((this_word - 1/nb_features) > EPSILON).astype(int)
            this_word_binarized = tuple(this_word_binarized.tolist())

            this_word_usage = np.sum(U.values[:,i])

            # Now remember it
            all_candidate_words[this_word_binarized] += this_word_usage

            # DEBUG PRINT
            print("word =", this_word_binarized, "; usage =", this_word_usage)

        #print('Step finished')
        

      

        # Update alpha so DL will be more sparse in the next step
        # and look for higher-order combinations
        iternb += 1
        alpha += iternb/nb_features
        # TODO : THIS STEP WILL BE TWEAKED. notably a larger step when nb_of_features > 20 so we will not bog down
        # Which is why I now update it by adding iternb/nb_features
        # This way longer words are not unduly favored
        
        # TODO Also stop regardless after 2*k iterations ?




    # Return the dictionary giving all candidate words and their total usage
    return all_candidate_words




# ----------- Step 2 of the MODL algorithm




def build_best_dict_from_library(data, library, queried_words_nb,
                error_function = None,
                nb_threads = 1,
                transform_alpha = 1E-100):
    """
    Given a data matrix and a library, will select the best n = queried_words_nb words with a greedy algorithm from the library to rebuild the data

    data is a matrix with one transaction per row and one element per column, as usual

    The greeedy algorithm works because the problem in this particular instance is submodular (TODO ADD MORE DETAILS)
    """





    """
    JUNE 2020 : I CHANGED transform_alpha DEFAULT TO 1E-100 ACCORDING TO THE CONCLUSION DRAWN BEFORE
    """





    print("queried_words_nb =",queried_words_nb)



    # --- Subsampling to save time
    # This algo involves making LOTS of reconstructions of the query matrix
    # (self.data) to find the best word at each iteration. To save time
    # while not losing much info (since we are working with overlap flags with
    # high redundancy, do sort-of a minibatch : limit the size of self.data
    # to a random selection of around 1M lines ? (this can also be a parameter depending on the computer available)

    # TODO : for better randomness, maybe redo this subsampling every time ?
    # It's not like sampling takes long anyways...

    # TODO : This varibale is not not used here, hence double check that the subsampling is done properly 
    # when creating the Modl object
    """
    This is done at creating the object, with squishing.
    """





    # --- Create objects
    # INITIAL DICT

    nb_features = data.shape[1]

    dictionary = [tuple([0] * nb_features)] # Must initialize with a word, might as well be full zeros.
    D = np.array(dictionary)

    # Coder object, to be updated iteratively

    # Rq : must make it float64 for some reason
    coder = SparseCoder(dictionary=D.astype('float64'),
        transform_algorithm='lasso_cd',
        transform_alpha = transform_alpha,
        positive_code = True, # Very important to ensure we also force positivity here !
        n_jobs = nb_threads)


    # Like in step 1, I use lasso coordinate descent because LARS had trouble with
    # highly correlated variables, which I saw experimentally
    # For alpha we will likely use an alpha of zero because the manhattan distance ensures
    # longer words are not prioritized DOUBLE CHECK



    # --- Filtering the candidates
    # Iteratively add candidates that improve the reconstruction


    stop = False
    while stop is False:

        # Update coder with the current dictionary, with the word added a last iteration
        D = np.array(dictionary)
        coder.components_ = D




        # Generate all_candidate_words for this iteration
        # Candidate words to be added are simply all words of the library that have not been
        # already used
        all_candidate_words = tree.get_all_candidates_except(library, D)










        # Test all words for this iteration
        candidates_with_errors = dict()
        #previous_error = np.inf
        print("----------------------")
        print("Current dictionary :",dictionary)

        print("> Candidate\tError")

        for candidate in all_candidate_words:


            #print("Testing candidate :",candidate)
            


            dict_being_tested = dictionary + [candidate]
            Dt = np.array(dict_being_tested)
            # Update dictionary to the one being currently tested
            coder.components_ = Dt.astype('float64')


            # Try for this encoding



            # Ignore convergence warnings likely due to zero alpha TODO ascertain that
            with warnings.catch_warnings():
                from sklearn.exceptions import ConvergenceWarning
                warnings.filterwarnings("ignore", category=ConvergenceWarning)

                
                encoded = coder.transform(data)
                rebuilt_data = np.matmul(encoded, Dt)



            # TODO : perhaps there could be a supervision here (not just reconstruction error I mean). Ask Cécile.
            # Yep. This is a potential improvemnt for the user !
            # Can be customized by plugging here !


            """
            On 110000many+110011+000011many It rebuild it very poorly, trying to use 110011 first,
            before 110000 ! Why ? Okay due to NMSE. IMPORTANCE OF USING ANOTHER ERROR IS PROVEN ! CECILE AGREED !
                Furthermore now I have the tree that privileges short words first.
                Just be careful it does not use 110011 as the second word instead. Likely would be due to the error.

            SO WE SWITCHED TO MANHATTAN ERROR. SAY SO IN PAPER.

            # So yeah, manhattan correctly says what we want : it computes a
            # dist(a,c) as half dist(b,c) instead of a quarter like euclidian does !
            """

            # We use manhattan distance by default to discourage using compromises such 
            # as 101 to rebuild 100 or 001 when the word 100 or 001 themselves are available
            def manhattan_dist(a,b):
                error_mat = np.abs(b-a)
                error = np.sum(error_mat)
                return error




            # UPDATE Add an alpha that is nonzero so that using longer words is still
            # an improvement even when other words cover it (and as a tiebreaker)
            # but keep it low (1/k) so you don't re-encourage compromise
            def manhattan_dist_with_sparsity(X_true, X_rebuilt, code):
                error_mat = np.abs(X_rebuilt-X_true)
                error = np.sum(error_mat)

                ALPHA = 1/X_true.shape[1]
                regul = np.sum(code) * ALPHA

                final_error = error + regul
                return final_error



            # Check for error function :
            # Default is manhattan dist
            if error_function is None :
                #error_function = manhattan_dist
                error_function = manhattan_dist_with_sparsity

            #error = error_function(X, rebuilt_X)
            error = error_function(data, rebuilt_data, encoded)


            # TODO : here it should be this

            # Remember error
            candidates_with_errors[tuple(candidate)] = error
            #print("> Error = ",error)
            print(str(candidate)+'\t'+str(error))


      
        # Add best candidate to final dict and remove from candidates
        best_candidate = min(candidates_with_errors, key=candidates_with_errors.get)
        # TODO and note in paper : in case of a tie, I think the first word is selected.

        print("Candidate", str(best_candidate), "was chosen")

        # This line of code was supposed to represent this idea : "once a word
        # is indeed added to the dictionary, remove it from candidates" but since
        # when querying candidates I automatically get only the non-used words,
        # it is redundant. Maybe keep a comment to that effect in case the
        # candidate selection model changes.
        #candidate_words.remove(best_candidate)



        dictionary += [best_candidate]







        # vvvvvvvvvvvvvvvvv PROBABLY DISCARD WHAT IS BELOW vvvvvvvvvvvvvvv #
        # I have proven my algorithm is optimal no ? Heh. In any case bothering
        # with what is below is likely more trouble than it's worth.
        # Since I will no longer use my "use children of present only" stuff
        # it should be good enough


        """
        LEAVE THIS AS A DRAFT ! NECESSARY TO ADDRESS SCALING PROBLEMS !!!!!!!!
        SAME FOR THE "use children of present only" stuff
        """

        """
        # If we have more words than the user wants, remove the worst one
        if (len(dictionary) > queried_words_nb):

            # TODO Remove worst word with respect to error (try all, opposite of above)

            return


        # When the word added and the word removed are the same (meaning
        # the word it wants to add does not improve on the dict) stop and
        # return
        # Cécile and Marina agreed.

        # TODO Also automatically stop after a certain number of iterations !!!!!!

        if best_candidate == last_word_removed :
            # I think final_dict will contain nodes. At the end of the algo,
            # extract the words from the nodes
            best_words = final_dict
            return best_words
        """

        # ^^^^^^^^ END OF STUFF I'LL DISCARD THIS PARTICULAR TIME ^^^^^^^^ #
        # TODO MAKE SURE IT WAS ALSO REMOVED FROM THE PAPER SINCE IT IS SUBMODULAT NOW

        ## Stop condition : as soon as we have enough words
        # Remark : since the dictionary always contain the (0,0,0,...) word as a supplement, we
        # must stop when it contains the required number of words... plus one !
        if (len(dictionary) >= queried_words_nb + 1):
            stop = True


    print("----------------------")
    print("Final dictionary :",dictionary)

    # Conclusion of the function : remove the (0,0,0,...) word and make words unique, just in case
    best_words = list(set(dictionary))
    best_words.remove(tuple([0] * nb_features))
    return best_words
