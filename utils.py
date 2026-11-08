"""
Utility functions used by CASet and DISC.
"""

import hashlib
import math
import re


def jaccard(a, b):
    """
    Compute the Jaccard distance between two sets. Note that Jacc({}, {}) = 0.
    :param a: a set
    :param b: another set
    :return: Jacc(a, b)
    """
    union = len(a | b)
    intersection = len(a & b)
    return (1 - intersection / union) if union > 0 else 0


def ancestor_sets(newick):
    """
    Find the ancestor sets of every mutation in tree represented in Newick format.
    :param newick: a Newick string representation of tree
    :return: a dict of each mutation's set of ancestors
    """
    ancestors = dict()
    _ancestor_sets(newick.replace(' ', ''), ancestors, set())
    return ancestors


def _ancestor_sets(newick, ancestors, current_ancestors):
    """
    Recursive helper function for ancestor_sets().
    :param newick: a Newick string representation of a subtree
    :param ancestors: a pointer to the dict being computed
    :param current_ancestors: the mutations ancestral to every mutation in the current subtree
    """
    newick = newick.strip()
    if not newick:
        raise ValueError("Árvore Newick contém um ramo vazio.")

    last_paren_index = newick.rfind(')')
    label_string = _strip_branch_length(newick[last_paren_index + 1:])

    if label_string == '':
        # ID estável mesmo quando os ramos filhos aparecem em outra ordem.
        leaves = sorted(_leaf_labels(newick))
        digest = hashlib.sha1(','.join(leaves).encode('utf-8')).hexdigest()[:12]
        labels = {f"__internal_{digest}"}
    elif label_string[0] == '{':
        labels = set(label_string[1:-1].split(','))
    else:
        labels = {label_string}

    for label in labels:
        ancestors[label] = labels | current_ancestors

    if last_paren_index != -1:
        for subtree_newick in outer_comma_split(newick[1:last_paren_index]):
            _ancestor_sets(subtree_newick, ancestors, labels | current_ancestors)


def _leaf_labels(newick):
    """Return terminal labels below a subtree, independent of sibling order."""
    newick = newick.strip()
    last_paren_index = newick.rfind(')')

    if last_paren_index == -1:
        label = _strip_branch_length(newick)
        return [label] if label else []

    leaves = []
    for child in outer_comma_split(newick[1:last_paren_index]):
        leaves.extend(_leaf_labels(child))

    return leaves


def _strip_branch_length(label):
    """Remove a standard numeric Newick branch-length suffix from a label."""
    return re.sub(
        r":\s*[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$",
        "",
        label.strip(),
    )


def outer_comma_split(newick):
    """
    Split a Newick substring on commas, ignoring those contained in
    parentheses and braces.
    """
    chunk_start = 0
    parens = 0
    brackets = 0

    for i, char in enumerate(newick):
        if char == '(':
            parens += 1
        elif char == ')':
            parens -= 1
        elif char == '{':
            brackets += 1
        elif char == '}':
            brackets -= 1
        elif char == ',' and parens == 0 and brackets == 0:
            yield newick[chunk_start:i]
            chunk_start = i + 1

    yield newick[chunk_start:]


def get_matrix_string(matrix):
    """
    Construct a human-readable string displaying a distance matrix.
    :param matrix: the matrix to display
    :return: a string displaying the matrix
    """
    string = 'tree\t' + '\t'.join(map(str, range(len(matrix)))) + '\n\n'
    for index, row in enumerate(matrix):
        string += '{}\t'.format(index) + '\t'.join(
            str(round(entry, 4)) for entry in row
        ) + '\n'

    return string


def get_min_max_string(matrix):
    """
    Determine the pairs of trees with min and max distances in a distance matrix and return a string displaying them.
    :param matrix: a matrix of distances
    :return: a string describing the min- and max-distance tree pairs
    """
    n = len(matrix)
    min_value = float('inf')
    min_indices = set()
    max_value = -float('inf')
    max_indices = set()

    for row in range(n):
        for col in range(row + 1, n):
            entry = matrix[row][col]
            if not math.isfinite(entry):
                continue
            if entry == min_value:
                min_indices.add((row, col))
            elif entry < min_value:
                min_indices = {(row, col)}
                min_value = entry

            if entry == max_value:
                max_indices.add((row, col))
            elif entry > max_value:
                max_indices = {(row, col)}
                max_value = entry

    if min_value == float('inf'):
        return 'No defined pairwise distances (fewer than two comparable mutations).\n'

    return (
        'Max distance: {} between tree pairs {}\n'
        'Min distance: {} between tree pairs {}'.format(
            round(max_value, 6),
            max_indices,
            round(min_value, 6),
            min_indices
        )
    )


def get_newicks(input_file):
    """
    Extract Newick strings from an input file with one tree per line.
    Treat anything after a # as a comment.
    :param input_file: the filename from which to grab trees
    :return: a list of Newick strings extracted from the file
    """
    trees = []
    with open(input_file) as f:
        for line in f.readlines():
            line = line.strip().strip(';').strip()
            if '#' in line:
                line = line[:line.index('#')]
            if line != '' and not line.isspace():
                trees.append(line)

    return trees
