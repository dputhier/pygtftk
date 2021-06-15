"""
A tree (or more accurately, graph) based representation of combinations of elements.

Combinations are represented as tuples. For example, if the possible elements are {A,B,C,D}, the combination A+C is represented as (1,0,1,0).

Nodes will be assigned under the following principle : each node X that contains all flags of a node Y.
For example, (1,0,0) is a parent of (1,0,1) but not the other way around.

A Library is a set set of nodes starting with a root of (0,0,0,...)

Author : Quentin Ferré <quentin.q.ferre@gmail.com>
"""

from functools import partial
from itertools import chain

import graphviz as gv
import numpy as np

from pygtftk.stats.intersect.overlap import overlap_regions as oc
from pygtftk.utils import message


class Node:
    """
    Represents a node of the tree. Each node is a word. 
    A word must be a tuple such as (1,1,0,0) meaning A+B if the possibilities are {A,B,C,D}
    """

    def __init__(self, word):
        if not isinstance(word, tuple):
            raise ValueError("`word` must be a tuple, such as (1,0,1)")

        self.word = word

        self.children = []
        self.parents = []

        # Default stats
        self.s = 0
        self.pval = 0
        self.fc = 0

    def add_child(self, to_add):
        self.children += [to_add]
        to_add.parents += [self]

    # Fallback concatenated representation
    def __str__(self):
        return ''.join([str(item) for item in self.word])


def apply_recursively_to_all_nodes(node, function, global_results,
                                   no_duplicates=True,
                                   stop_condition=lambda c, gr: False):
    r"""
    General utility function to recursively apply a function to all nodes of a tree.
    Also pass a global_results dict to be added to.

    Since nodes can have several parents, by default this will remember which nodes
    have been seen and not operate on them twice (`no_duplicates`).

    A stop condition fonction of signature stop_condition(current_node, global_results)
    can be passed : the recursion will stop if it returns True at any point. 

    >>> from pygtftk.stats.intersect.modl.tree import Library, apply_recursively_to_all_nodes
    >>> from pygtftk import utils
    >>> utils.VERBOSITY = 0
    >>> words = [(1,0,0), (1,1,0), (0,1,0), (1,1,1)]
    >>> l = Library()
    >>> l.build_nodes_for_words(words)
    >>> l.assign_nodes()
    >>> manual_print_words = set([str(n) for n in l.assigned_nodes])
    >>> gr = dict()
    >>> apply_recursively_to_all_nodes(l.root_node, str, gr)
    >>> recursive_print_words = set(gr.values())
    >>> assert recursive_print_words == manual_print_words

    """

    # Apply the function to the node
    result = function(node)

    # Record result; need to use a method to work on the reference and send back to the outer scope
    global_results.update({node: result})

    # Then, for all children of the node ...
    for c in node.children:

        # If no_duplicates is True, we proceed only if c is not already in the results
        if not ((c in global_results.keys()) and no_duplicates):

            # We also check for the stop function
            if not stop_condition(c, global_results):
                apply_recursively_to_all_nodes(c, function, global_results)  # Move to the child


class Library:
    r"""
    A tree (or more accurately, graph) based representation of combinations of elements.

    Combinations are represented as tuples. For example, if the possible elements are {A,B,C,D}, the combination A+C is represented as (1,0,1,0).

    Nodes will be assigned under the following principle : each node X that contains all flags of a node Y.
    For example, (1,0,0) is a parent of (1,0,1) but not the other way around.

    A Library is a set set of nodes starting with a root of (0,0,0,...)

    Example of use :

    >>> from pygtftk.stats.intersect.modl.tree import Library
    >>> words = [(1,1,0),(1,0,0),(0,1,0)]
    >>> l = Library()
    >>> l.build_nodes_for_words(words)
    >>> l.assign_nodes()

    """

    def __init__(self):
        self.unassigned_nodes = list()  # Nodes waiting for assignment
        self.nodes_were_assigned = False

        self.assigned_nodes = list()  # Nodes that were assigned

    # ----- Building

    def build_nodes_for_words(self, words):
        """
        Builds unassigned nodes for all words. 

        `words` must be a list of tuples like [(1,0,0),(1,1,0)]
        """

        # Check that the passed list of words is indeed a list of tuples 
        if not (isinstance(words, list) and all(isinstance(word, tuple) for word in words)):
            raise ValueError("`words` must be a list of tuples like [(1,0,0),(1,1,0)]")

        # Ensure all words are unique
        words = sorted(list(set(words)))

        # For each word in words, create a node and add it to the list
        for word in words:
            self.unassigned_nodes += [Node(word)]

            # Remember the length of the word.
            # TODO: Throw exception if there are words of different lengths
            word_size = len(word)  

        # Finally, create root node, all unassigned nodes will branch from it later
        self.root_node = Node(word=tuple([0] * word_size))
        # Add the root to the assigned words
        self.assigned_nodes.append(self.root_node)

    def build_nodes_for_words_from_ologram_result_df(self, result_df, query_name="Query"):
        """
        Used for graphical display later on

        The nodes should also have a pval and fc and s that can be given and that we can query later on
        So the build_nodes_for_words function can also take a pandas dataframe that is an OLOGRAM output !!
        """

        # How many nodes will there be ?
        self.unassigned_nodes = [None] * len(result_df.index)

        combis_in_the_df = [None] * len(result_df.index)

        # Get combis from the df
        for index, row in result_df.iterrows():

            # Remove '[]' and spaces and split combi
            combi_raw = row["feature_type"].translate({ord(i): None for i in ['[', ']', ' ']})

            split_combi = combi_raw.split('+')
            if split_combi[0] == 'Query': split_combi[0] = str(
                query_name)  # Replace query name with potential custom name

            # ALways read ellipsis character as three dots instead
            ELLIPSIS = u'\u2026'
            for i in range(len(split_combi)):
                if split_combi[i] == ELLIPSIS: split_combi[i] = "..."

            combis_in_the_df[index] = tuple(split_combi)

            message("Read this combination : " + str(tuple(split_combi)))

            # NOTE If there are "..." like in "[ TAL1 + MYC + ... ]" for the non-exact combis,
            # it should still work because there is usually no combi with only "..." that would be a shared parent

        # Translate to binary, and remember features_names
        sorted_features = sorted(set(chain.from_iterable(combis_in_the_df)))

        # Ensure the query name is always first and '...' is always last
        if query_name in sorted_features:
            sorted_features.remove(query_name)
            sorted_features.insert(0, str(query_name))

        if '...' in sorted_features:
            sorted_features.remove('...')
            sorted_features.append('...')

        self.features_names = sorted_features

        for index, row in result_df.iterrows():

            # Convert word to tuple, matching the order of features that was computed above
            word_as_strings = combis_in_the_df[index]
            word = [None] * len(self.features_names)
            for i in range(len(self.features_names)):
                if self.features_names[i] in word_as_strings:
                    word[i] = 1
                else:
                    word[i] = 0
            word = tuple(word)

            word_size = len(
                word)  # Remember the length of the word. TODO Throw exception if there are words of different lengths

            # Query relevant info
            self.unassigned_nodes[index] = Node(word)
            self.unassigned_nodes[index].s = row["summed_bp_overlaps_true"]
            self.unassigned_nodes[index].pval = row['summed_bp_overlaps_pvalue']
            self.unassigned_nodes[index].fc = row["summed_bp_overlaps_log2_fold_change"]

        # Finally, create root node, all unassigned nodes will branch from it later
        self.root_node = Node(word=tuple([0] * word_size))
        # Add the root to the assigned words
        self.assigned_nodes.append(self.root_node)

    # ----- Assigning

    def assign_nodes(self, min_inheritance = 0):

        # Sort nodes in this list by total number of nonzero flags, for assignation purposes.
        def nb_nonzero_flags_node(node):
            return sum(tuple([bool(flag) for flag in node.word]))

        self.unassigned_nodes.sort(key=nb_nonzero_flags_node, reverse=True)
        # Each height level of the tree must contain nodes of the same total flag number.
        # So we must add them in increasing order of total flags, otherwise adding a node of
        # total 3 before a total 2 would result in the total 3 attaching to a total 1 or to a root.
        # Order of nodes which have same total length is not important, as nodes do not take as parents nodes of the same length
        # NOTE : this means the tree must be built all at once !

        # The tree must be built all at once. See above.
        if self.nodes_were_assigned: raise TypeError(
            "Cannot assign nodes after nodes were already assigned once. You must rebuild a new Library with all nodes, old and new.")

        # For all unassigned nodes, browse through all other nodes to get it a parent
        while len(self.unassigned_nodes) > 0:

            unode = self.unassigned_nodes.pop()

            # Find closest node with LESS flags, default to root
            new_parent = None

            # Distance is simply number of different flags
            def dist(other_node):

                # Only add it if it is also an exact parent : the child must have
                # all the elements of its parent PLUS potentially others
                # This should also include the fact that is has less flags
                if oc.does_combi_match_query(unode.word, other_node.word, exact=False):
                    dist = [(a - b) ** 2 for a, b in zip(unode.word, other_node.word)]
                    return sum(dist)
                else:
                    return np.inf

            # Get all distances
            all_distances_less_flags = dict()
            # apply_recursively_to_all_nodes(self.root_node, dist, all_distances_less_flags)
            for anode in self.assigned_nodes:
                all_distances_less_flags[anode] = dist(anode)

            # This throws a ValueError if all_distances_less_flags is empty
            # In case of a tie, add all as parents
            try:
                closest_distance = min(all_distances_less_flags.values())
                new_parents_list = [node for node in all_distances_less_flags if
                                    all_distances_less_flags[node] == closest_distance]
            except:
                new_parents_list = []

            # Add current node as child to new_parent
            for new_parent in new_parents_list:

                # Only add the child if it accounts for at a proportion of at least min_inheritance of its parent
                if unode.s/(new_parent.s+1E-100) >= min_inheritance:

                    message('Adding ' + str(unode) + ' to ' + str(new_parent) + ' as distance of ' + str(
                        all_distances_less_flags[new_parent]), type='DEBUG')
                    new_parent.add_child(unode)

            # Now add the node to the list of assigned nodes
            self.assigned_nodes.append(unode)

        self.nodes_were_assigned = True


# def generate_candidates(library, currently_used_words):
#     """
#     From the library, take the children of all nodes that have the flag "used", 
#     unless those children have the flag "used" themselves.
#
#     currently_used_words must be a list of tuples.
#
#     Currently not used, kept for future evolutions.
#     """
#
#     def get_candidates(node):
#         candidates = []
#         if node.word in currently_used_words:
#             for c in node.children:
#                 if c.word in currently_used_words: candidates += get_candidates(c)
#                 else: candidates += [c.word]
#         return candidates
#
#     # Aaaaand... start !
#     candidates = get_candidates(library.root_node)
#     # TODO: because nodes can have several parents, make the candidates list
#     # equal to a set of itself (hence unique words) each time !
#     return candidates


def get_all_candidates_except(library, exclude):
    """
    Returns the words of all nodes, except those words that are in the exclusion list
    """

    # if exclude is an np array, convert it to a list of tuples
    if isinstance(exclude, np.ndarray):
        exclude = [tuple(w) for w in exclude]

    def get_candidates(node):
        candidates = []

        for c in node.children:
            candidates += get_candidates(c)

        if node.word not in exclude:
            candidates += [node.word]

        return candidates

    # And start !
    candidates = get_candidates(library.root_node)

    # Because nodes can have several parents, make the candidates
    # list equal to a set of itself (hence unique words) each time !
    return list(set(candidates))


def colorize(value, maximum=320):
    # Default to 0 if value is NaN
    if np.isnan(value): value = 0

    # RGB values for extremes
    max_neg = np.array([255, 60, 30])
    zero = np.array([255, 255, 255])
    max_pos = np.array([30, 60, 255])

    absolute_value = np.clip(abs(value), -maximum, maximum)

    ratio = absolute_value / maximum

    if value < 0:  color = ratio * max_neg + (1 - ratio) * zero
    if value >= 0: color = ratio * max_pos + (1 - ratio) * zero

    r, g, b = color.astype(int)
    return '#%02x%02x%02x' % (r, g, b)


def output_visualize(tree, output_path, features_names=None):
    """
    Output a visualisation a a given Library (the `tree` argument)

    :param tree: : The library to write
    :param output_path: Path to write to
    :param features_names: Conversion key giving the name of each feature in the vector. For example if features_names = ['A','B','C'], (0,1,1) will be translated as B+C
    """

    # If features names is None, query tree.features_names
    # tree.features_names will be set when building a Library from an OLOGRAM result df
    if features_names is None:
        try:
            features_names = tree.features_names
        except:
            features_names = None

    root_node = tree.root_node  # Get root node

    # Quick workaround : this should be only called by ologram_modl_treeify
    # If the filepath ends in '.pdf', remove it so graphviz can put a '.dot'
    # and a '.pdf' at the end later
    if output_path.endswith(".pdf"): output_path = output_path[:-4]

    s = gv.Digraph('combi_tree', filename=output_path + '.dot',
                   graph_attr={'splines': 'compound'},
                   strict=True,  # Duplicate edges are not allowed (!)
                   node_attr={'shape': 'plaintext', 'fontname': 'Helvetica'})

    # ------ Utility functions

    # Convert a tree node to proper features, giving a string for the display graph
    def node_to_combi_string(node, features_names=None, new_line_every=2):

        if features_names is not None:

            # Root node special handler
            if sum(node.word) == 0: return '-'

            presents = [features_names[i] for i in range(len(node.word)) if node.word[i] != 0]

            # Produce result string
            result = ""
            newline_counter = 0

            for i in range(len(presents)):
                result += str(presents[i])

                if not (i == len(presents) - 1): result += " + "  # add '+' until the last element
                if (newline_counter == new_line_every):
                    result += "<br/>"  # Use an escaped newline
                    newline_counter = 0
                newline_counter += 1

            return result

        # Fallback
        return str(node)

    def format_node_string(combi_string, s_val, p_val, fc_val):

        color_hex = colorize(fc_val, maximum=10)

        res = '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">'
        res += '<TR> <TD PORT="f1" BGCOLOR="' + str(color_hex) + '"><FONT POINT-SIZE="16"><b>' + str(
            combi_string) + '</b></FONT></TD> </TR>'
        res += '<TR> <TD>S = ' + str(s_val) + '<BR/>p-val = ' + '{0:.4g}'.format(
            p_val) + '<BR/>log2(FC) = ' + '{:.5f}'.format(fc_val) + '</TD> </TR> </TABLE>>'
        return res

    # Prepare the function
    def produce_dot_for_node(node, graph):
        for c in node.children:
            node_name = node_to_combi_string(node, features_names)
            child_name = node_to_combi_string(c, features_names)

            message("Drawing " + node_name + ' --> ' + child_name, type='DEBUG')

            ## Add nodes
            # Only add node if not already present of course.
            # If present, the graph's 'body' contains the combi string prefixed with a tab character

            ## Parent
            if not ('\t' + node_name in s.body):
                # print(combi_string, node.s, node.pval, node.fc)
                graph.node(node_name, format_node_string(node_name, node.s, node.pval, node.fc))

            # Child
            if not ('\t' + child_name in s.body):
                graph.node(child_name, format_node_string(child_name, c.s, c.pval, c.fc))

            graph.edge(node_name + ':s', child_name + ':n')

        return 1

    # Iterate over all nodes
    mygraphfunc = partial(produce_dot_for_node, graph=s)
    global_results = {}
    apply_recursively_to_all_nodes(root_node, mygraphfunc, global_results)

    # DEBUG
    # apply_recursively_to_all_nodes(root_node, print, global_results)

    # Now save it
    s.save(output_path + '.dot')
    s.render(output_path, format='pdf', view=False, cleanup=True)
