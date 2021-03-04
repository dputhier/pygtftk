"""
Implementation of Apriori algorithm. Use this to compare MODL to apriori.

It mines for association rules and not complexes per se, but can still be useful.

Credits : based on code from https://github.com/coorty/apriori-agorithm-python combined with https://github.com/asaini/Apriori/blob/master/apriori.py

TODO: Implement it as a potential alternative to MODL ? See notes in overlap_stats_compute.py for how to do it easily.
"""

from collections import defaultdict

import numpy as np
import pandas as pd


class Apriori:

    def __init__(self, min_support):
        """
        Apriori-based itemset miner.
        
        min_support is the minimum frequency that an itemset must have to be considered "frequent".

        To run apriori, use the following steps, if X is a matrix with one line per transaction, one column per item and a '1' if the item is presnet in the transaction :

        >>> X = [[1,1,0],[1,0,1]]
        >>> names = ['A','B','C']
        >>> from pygtftk.stats.intersect.modl.apriori import Apriori, matrix_to_list_of_transactions, apriori_results_to_matrix
        >>> transactions = matrix_to_list_of_transactions(X, names)
        >>> myminer = Apriori(min_support = 0)
        >>> myminer.run_apriori(transactions)
        >>> results = myminer.produce_results()
        >>> apriori_results_df = apriori_results_to_matrix(results, names)
        
        """
        self.min_support = min_support

        self.apriori_was_run = False

    # ---- Utility functions ---- #

    ## Initialization
    def get_one_item_set(self, transactions):
        """ Get unique 1-item set in `set` format """
        itemSet = set()
        for line in transactions:
            for item in line: itemSet.add(frozenset([item]))
        return itemSet

    def get_joined_item_set(self, termSet, k):
        """ Generate new k-terms candiate itemset"""
        return set([term1.union(term2) for term1 in termSet for term2 in termSet
                    if len(term1.union(term2)) == k])

    def get_support(self, item):
        """ Get the proportional support of an item"""
        return self.item_count_dict[item] / self.nb_transactions

    def get_items_with_min_support(self, transactions, item_set, frequent_items_set, min_support):
        """ Calculates the support for items in the itemSet and returns a subset
        of the itemSet each of whose elements satisfies the minimum support"""
        result_items_set = set()

        _local_set = defaultdict(int)

        for item in item_set:
            for transac in transactions:
                # Does the transaction contain the item ?
                if item.issubset(transac):
                    frequent_items_set[item] += 1
                    _local_set[item] += 1

        # Only conserve frequent item-set
        n = self.nb_transactions
        for item, count in _local_set.items():
            support = float(count) / len(transactions)

            if support >= self.min_support:
                result_items_set.add(item)

        return result_items_set

        # ---- Main function ---- #

    def run_apriori(self, transactions):
        """
        transactions must have format of lst of transactions ? EXPLAIN
            It is a list that contains sets, eg. [(A,B),(B,C),(A,B,D)]
        Run the apriori algorithm, return the frequent itemsets. 
        """

        ## Initialization of results variables
        self.nb_transactions = len(transactions)  # How many transactions are there ?

        frequent_item_sets = dict()  # a dict store all frequent *-items set
        # Its structure is dict[k] = [all itemsets of length k]

        ## Dictionary to hold itemset counts
        # Key = candidate k_item set ; value = its count
        item_count_dict = defaultdict(int)

        # Begin the algorithm with all 1-item sets
        item_set = self.get_one_item_set(transactions)
        self.unique_items = item_set

        # Get the frequent 1-item sets
        freq_one_item_set = self.get_items_with_min_support(transactions, item_set, item_count_dict, self.min_support)

        # --- Main loop
        # Main idea is to "grow" in length
        k = 1
        current_frequent_term_set = freq_one_item_set

        while current_frequent_term_set != set():
            frequent_item_sets[k] = current_frequent_term_set  # Save result
            k += 1

            # Get new candiate k-terms set
            current_candidate_item_sets = self.get_joined_item_set(current_frequent_term_set, k)

            # Now restrict to only frequent k-terms set
            current_frequent_term_set = self.get_items_with_min_support(transactions, current_candidate_item_sets,
                                                                        item_count_dict, self.min_support)

        # Save results
        self.item_count_dict = item_count_dict
        self.frequent_item_sets = frequent_item_sets
        # Only frequent items(a dict: freqSet[1] indicate frequent 1-term set)

        # Now for what to return
        self.apriori_was_run = True

    def produce_results(self):
        if not self.apriori_was_run: raise Exception("Must run apriori first")

        # Return the items, along with their support
        to_return_items = list()
        for k, itemsets in self.frequent_item_sets.items():
            to_return_items.extend([(tuple(iset), self.get_support(iset))
                                    for iset in itemsets])

        # NOTE To get the association rules, simply go over all itemsets and 
        # compute the confidence of rules.
        # A rule is: "When X is present, Y is also present" where X and Y can be
        # single items or itemsets themselves.
        # For example confidence (X --> Y) = support({X,Y})/support(X)

        return to_return_items


# The apriori algo requires transactions to be in a list of list format, one transaction per list.
# So we must convert our overlap/transaction matrix into such a list of lists.
def matrix_to_list_of_transactions(x, names):
    """
    From a matrix with one line per transaction and one column per element with 1 if present and 0 if absent, 
    returns a list of transaction
    """
    # Enforce type
    names = np.array(names)
    x = np.array(x)

    # Get the list of all nonzero elements in each row
    result = []
    for row in x:
        current_transaction = np.nonzero(row)[0]
        current_transation_items = names[current_transaction]
        result += [current_transation_items.tolist()]
    return result


def apriori_results_to_matrix(results, names):
    """
    Turns back apriori results of the form [(itemset, support)] to a matrix of words
    Does not further filter by support, keeps all itemsets. Filter by support before passing.
    """

    matrix = []

    names = list(names)  # Convert names to a list if an array

    for res in results:
        itemset = res[0]

        itemset_matrix = [0] * len(names)
        for elem in itemset:
            itemset_matrix[names.index(elem)] = 1

        matrix += [itemset_matrix]

    # Dataframe of results
    resdf = pd.DataFrame(matrix, columns=names)
    resdf['support'] = [res[1] for res in results]
    resdf.sort_values(by=['support'], inplace=True, ascending=False)

    return resdf
