"""
A module to shuffle BED files and generate new "fake" BED files.
This shuffle keeps the distribution of regions and inter-region legths.
"""

import numpy as np
cimport numpy as np

import pybedtools

from functools import partial
from multiprocessing import Pool


################################################################################
# -------------------------- Shuffling bed files ----------------------------- #
################################################################################


# ------------------------ Custom simple permutation ------------------------- #

def shuffle(arr):
    """
    Given a numpy array, will shuffle its rows independantly.
    """
    x = arr.shape[0]
    y = arr.shape[1]
    rows = np.indices((x, y))[0]
    cols = np.asarray([np.random.permutation(y) for _ in np.arange(x)])
    result = arr[rows, cols]
    return result


# ------------------------ Custom Markov shuffling --------------------------- #

"""
This will shuffle lists using an order 2 Markov model.

WARNING : A Markov shuffling here is not strictly a shuffle, since the resulting
arrays will have different elements. It is more akin to generating a new array
based on the Markovian characteristics of the old one.

This is not recommended in the general case, and should only be used if you
suspect there is an order to the data that you want to keep.

Markov shuffling is very much in BETA-TEST at the moment. You should not rely on
it just yet.

Please note that :
  - This will be very time consuming (hours).
  - A negative binomial cannot be used here, the resulting distribution is often
  multimodal and not a good fit for any of the basic probability distributions.
  - This will result in a higher number of intersections and number of overlapping
  base pairs than the basic shuffle. Even when comparing a BED file against itself,
  you may find the BED file to be "significantly anticorrelated with itself" because
  the Markov shuffle will generate a new BED file based on the characteristics
  of the old one, resulting in a more patterned distribution.
"""

cdef noise(int x, int factor=1000):
    """
    When generating new lengths, add some noise so they are not always multiples
    of 1000 (or the factor). This could have artificially increased (or decreased)
    the number of overlapping nucleotides if the original regions were, for
    example, mostly shorter than 1000.

    Generates random uniform noise between 0 and factor-1, and substract it to x.
    """
    cdef int noise
    noise = int(np.random.uniform() * abs(factor-1))
    return x - noise

cdef roundup(int x, int factor=1000):
    """
    We round up the lengths to the thousands to avoir overfitting the model by
    learning the exact values.
    """
    return x if x % factor == 0 else x + factor - x % factor


def learn_word_dict(corpus):
    """
    Turn a corpus (a list) into a dictionary of all consecutive trios of elements.

    The dictionary gives, for each pair of elements, all the elements that were
    found following this pair in the corpus.
    """

    def make_trios(corpus):
        # Must have at least 3 features
        if len(corpus) < 3 :
          raise ValueError('At least one chromosome in one of the input bed files had fewer than 3 features. You cannot use Markov shuffling (order 2) in this case.')

        for i in range(len(corpus)-2):
            yield (corpus[i], corpus[i+1], corpus[i+2])

    trios = make_trios(corpus)

    word_dict = {}
    for k_1, k_2, k_3 in trios:

        # Values are rounded to the thousands to avoid overfitting
        k_1, k_2, k_3 = roundup(k_1), roundup(k_2), roundup(k_3)

        key = (k_1, k_2)
        if key in word_dict.keys():
            word_dict[key].append(k_3)
        else:
            word_dict[key] = [k_3]

    return word_dict


cdef generate_markov_cython(int dummy, np.ndarray corpus, dict word_dict):
    """
    Given a corpus (list of elements in 1D numpy array format) and the corresponding
    word_dict generated by learn_word_dict from this corpus, create a Markov-shuffled
    corpus of the same length.

    The dummy argument is not used here, it's only so we can multiprocess it easily later.
    """

    cdef int first_word
    cdef int second_word
    cdef list chain
    cdef int n_words
    cdef np.ndarray rangei
    cdef (int, int) last_words
    cdef int result

    # Randomly pick two words to begin
    first_word, second_word = roundup(np.random.choice(corpus)), roundup(np.random.choice(corpus))
    chain = [first_word, second_word]
    n_words = len(corpus) -2 # -2 because we already have 2 words in the chain

    #nb_of_randos = 0
    rangei = np.arange(n_words)
    for i in rangei:
        last_words = roundup(chain[-2]), roundup(chain[-1])

        # If we by chance hit the second-to-last and last element, restart somewhere else
        try:
            result = np.random.choice(word_dict[last_words])
        except KeyError:
            result = np.random.choice(corpus)
            #nb_of_randos = nb_of_randos + 1
        chain.append(noise(result)) # Add noise

    return chain




# Lambda wrapper to be able to use Python multiprocessing
def generate_markov(arguments):
    dummy, corpus, word_dict = arguments
    return generate_markov_cython(dummy, corpus, word_dict)


def markov_shuffle(arr, nb_threads = 8):
    """
    This function takes an array as input, because it is meant to be slotted
    at the place of shuffle(), which shuffles the rows of an array independently.

    As the rows of the array to be fed are all identical at first (it is a tiling of a list),
    this learns a Markov model of order 2 on the first row of the array, then
    replaces each row with an independant realisation of the Markov model.

    WARNING : This is *very* long. Calling this across all chromosomes of a
    modest BED file (20K lines) for a batch size of 50 may take upwards of 1 minute.
    """
    # Take the first row since they are all supposed to be identical
    L = arr[0,] # Remember that L must be a numpy array
    wd = learn_word_dict(arr[0,]) # Learn a dictionary

    # Replace each row with a Markov shuffled version
    ranger = np.arange(arr.shape[0]) # The dummy argument must be a list of integers
    list_arguments = [(i,L,wd) for i in ranger]
    # Using a partial here resulted in heavy Python overhead, so I do this instead

    with Pool(nb_threads) as p:
       rows = p.map(generate_markov,list_arguments)

    result = np.array(rows)
    return result











################################################################################
# ---------------- Generating bed files from shuffled ones ------------------- #
################################################################################

cdef generate_fake_bed(Lr_shuffled, Li_shuffled, chrom):
    """
    Given a list of genomic region lengths and inter-genomic lengths, will
    generate a corresponding BED file as a list of lines.
    """
    cdef list fake_bed
    fake_bed = list()
    cdef int current_position
    current_position = 0
    cdef int k

    cdef np.ndarray rangek
    rangek = np.arange(len(Lr_shuffled))

    cdef int start
    cdef int end

    # We must begin with an inter-region length. Indeed, Li_shuffled has one
    # more element than Lr_shuffled, so we set before the loop current_position = Li_shuffled[0]
    # and use Li_shuffled[k+1] in the loop and carry on from here.
    current_position = Li_shuffled[0]
    for k in rangek :
        start = current_position
        end = start + Lr_shuffled[k]
        current_position = end + Li_shuffled[k+1]

        fake_bed.append((str(chrom),int(start),int(end)))


    # REMARK : as a perspective, we could use the same algorithm but keep the
    # negative inter-region lengths, to generate non-merged beds : this way,
    # we could work with peaks that have some overlp and not juste merge them,
    # and do some statistics on this within-set overlap. Something to consider.

    return fake_bed



def generate_fake_bed_for_i(i,shuffled_Lr1_batches,shuffled_Li1_batches, shuffled_Lr2_batches,shuffled_Li2_batches,all_chroms):
    """
    Partial call to the individual function 'generate_fake_bed'. Used in multiprocessing.
    """

    current_fake_bed_A = list()
    current_fake_bed_B = list()

    for chrom in all_chroms:
        genA = list()
        genB = list()
        genA = generate_fake_bed(
            shuffled_Lr1_batches[chrom][i], shuffled_Li1_batches[chrom][i], chrom)
        genB = generate_fake_bed(
            shuffled_Lr2_batches[chrom][i], shuffled_Li2_batches[chrom][i], chrom)
        current_fake_bed_A = current_fake_bed_A + genA
        current_fake_bed_B = current_fake_bed_B + genB
    return (current_fake_bed_A, current_fake_bed_B)


def batch_to_bedlist(shuffled_Lr1_batches, shuffled_Li1_batches,
                      shuffled_Lr2_batches, shuffled_Li2_batches,
                      all_chroms, minibatch_size, nb_threads = 8):

    """
    Wrapper to multiprocess generate_fake_bed across both bed files.
    """

    generate_for_them = partial(generate_fake_bed_for_i,
            shuffled_Lr1_batches=shuffled_Lr1_batches,shuffled_Li1_batches=shuffled_Li1_batches,
            shuffled_Lr2_batches=shuffled_Lr2_batches,shuffled_Li2_batches=shuffled_Li2_batches,
            all_chroms=all_chroms)

    bedsA = list()
    bedsB = list()

    # Multiprocess
    with Pool(nb_threads) as p:
        AB_tuples = p.map(generate_for_them,np.arange(minibatch_size))

    bedsA = [x[0] for x in AB_tuples]
    bedsB = [x[1] for x in AB_tuples]

    return bedsA,bedsB
