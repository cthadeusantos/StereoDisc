import argparse
import pickle
from pathlib import Path

import numpy as np

from utils import *
from plota import save_heatmap, get_filename_with_new_extension

"""
CASet.py calculates pairwise CASet distances between all trees in an input file with one Newick string on each line.
Defaults to intersection distance.
Command line arguments:
[-o --outputFile file] [-u --union] [-p --pickle file] [-t --treePrint] [-m --minmax]
"""


def caset(mutations, t1_ancestors, t2_ancestors):
    """
    Compute the CASet distance between two trees stored as ancestor sets.
    :param mutations: the set of mutations over which to sum
    :param t1_ancestors: a dict storing the t1 ancestor sets of every mutation
    :param t2_ancestors: a dict storing the t2 ancestor sets of every mutation
    :return: the CASet distance between t1 and t2
    """
    m = len(mutations)
    if m < 2:
        return np.nan
    total = 0
    for i in range(m):
        for j in range(i + 1, m):
            t1_ca_set = t1_ancestors[mutations[i]] & t1_ancestors[mutations[j]]
            t2_ca_set = t2_ancestors[mutations[i]] & t2_ancestors[mutations[j]]

            total += jaccard(t1_ca_set, t2_ca_set)

    return total / (m * (m - 1) / 2)


def caset_intersection(t1, t2):
    """
    Compute the CASet distance between two trees, only summing over pairs of mutations in both trees.
    :param t1: a Newick string representation of tree 1
    :param t2: a Newick string representation of tree 2
    :return: the CASet intersection distance between tree 1 and tree 2
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    intersection = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))

    return caset(intersection, t1_ancestors, t2_ancestors)


def caset_union(t1, t2):
    """
    Compute the CASet distance between two trees, summing over all pairs of mutations.
    :param t1: a Newick string representation of tree 1
    :param t2: a Newick string representation of tree 2
    :return: the CASet union distance between tree 1 and tree 2
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    union = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))

    for u in union:
        if u not in t1_ancestors:
            t1_ancestors[u] = set()
        if u not in t2_ancestors:
            t2_ancestors[u] = set()

    return caset(union, t1_ancestors, t2_ancestors)

def matrix_to_numpy(matrix, dtype=float):
    """Convert a matrix-like object (for example, a list of lists) to a NumPy array."""
    return np.asarray(matrix, dtype=dtype)


def save_numpy_matrix(matrix, file_path):
    """Save a matrix as a NumPy ``.npy`` file and return the array that was saved.

    ``file_path`` may omit the ``.npy`` extension; in that case NumPy adds it.
    """
    array = np.asarray(matrix)
    np.save(file_path, array)
    return array

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find pairwise CASet distances of a set of trees.')
    parser.add_argument('inputFile')
    parser.add_argument('--originals', help='Arquivo original usado para nomear o diretório do dataset.')
    parser.add_argument('-c', '--consensusFile', help='Arquivo com a árvore de consenso (Newick).')
    parser.add_argument('-o', '--outputFile')
    parser.add_argument('-u', '--union', action='store_true')
    parser.add_argument('-p', '--pickle')
    parser.add_argument('-t', '--treePrint', action='store_true')
    parser.add_argument('-m', '--minmax', action='store_true')
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()

    package_name = Path(args.originals or args.inputFile).stem
    output_dir = Path("output") / "CASet" / package_name
    output_dir.mkdir(parents=True, exist_ok=True)

    trees = get_newicks(args.inputFile)
    
    if args.treePrint:
        for i, tree in enumerate(trees):
            print('Tree {}:\n{}\n'.format(i, tree))

    distance_measure = caset_union if args.union else caset_intersection
    distance_matrix = [[None for x in range(len(trees))] for y in range(len(trees))]
    for i, t1 in enumerate(trees):
        for j, t2 in enumerate(trees):
            distance_matrix[i][j] = distance_measure(t1, t2)
            distance_matrix[j][i] = distance_matrix[i][j]

    output = get_matrix_string(distance_matrix)

    numpy_matrix = matrix_to_numpy(distance_matrix)

    pairwise_filename = get_filename_with_new_extension(args.inputFile, ".png", "CASet_trees")
    pairwise_numpy_filename = get_filename_with_new_extension(args.inputFile, ".npy", "CASet_trees")
    save_numpy_matrix(numpy_matrix, output_dir / pairwise_numpy_filename)
    save_heatmap(numpy_matrix, output_dir / pairwise_filename, title=pairwise_filename,
                 cmap="RdBu", show_values=args.show)

    if args.consensusFile:
        consensus_trees = get_newicks(args.consensusFile)
        if len(consensus_trees) != 1:
            parser.error('--consensusFile deve conter exatamente uma árvore Newick.')
        consensus = consensus_trees[0]
        comparison = [distance_measure(consensus, tree) for tree in trees]
        consensus_array = matrix_to_numpy(comparison)
        consensus_filename = get_filename_with_new_extension(args.inputFile, ".png", "CASet_consensus")
        consensus_numpy_filename = get_filename_with_new_extension(args.inputFile, ".npy", "CASet_consensus")
        save_numpy_matrix(consensus_array, output_dir / consensus_numpy_filename)
        save_heatmap(consensus_array.reshape(1, -1), output_dir / consensus_filename,
                     title=consensus_filename, cmap="RdBu", show_values=args.show,
                     figsize=(12, 2.6))

        defined_values = consensus_array[np.isfinite(consensus_array)]
        sorted_array = np.sort(defined_values)
        sorted_filename = get_filename_with_new_extension(
            args.inputFile, ".png", "CASet_consensus_sort")
        if sorted_array.size:
            endpoint_labels = [f"{sorted_array[0]:.4f}", f"{sorted_array[-1]:.4f}"]
            sorted_title = f"{sorted_filename} (comparáveis: {sorted_array.size}; indefinidas: {consensus_array.size - sorted_array.size})"
            save_heatmap(sorted_array.reshape(1, -1), output_dir / sorted_filename,
                         title=sorted_title, cmap="RdBu", show_values=args.show,
                         x_tick_labels=endpoint_labels, figsize=(12, 2.6))
        else:
            save_heatmap(np.array([[np.nan]]), output_dir / sorted_filename,
                         title=f"{sorted_filename} (sem pares comparáveis)",
                         cmap="RdBu", show_values=False, figsize=(12, 2.6))

    if args.minmax:
        output += get_min_max_string(distance_matrix)

    if args.pickle:
        with open(output_dir / Path(args.pickle).name, 'wb') as out:
            pickle.dump(distance_matrix, out, pickle.HIGHEST_PROTOCOL)

    if args.outputFile:
        with open(output_dir / Path(args.outputFile).name, 'w') as out:
            out.write(output)
    else:
        print(output)
