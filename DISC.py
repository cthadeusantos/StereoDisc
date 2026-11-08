import argparse
import pickle
from pathlib import Path
import numpy as np

from utils import *
from plota import save_heatmap, get_filename_with_new_extension

"""
DISC.py calculates pairwise DISC distances between all trees in an input file with one Newick string on each line.
Defaults to intersection distance.
Command line arguments:
[-o --outputFile file] [-u --union] [-p --pickle file] [-t --treePrint] [-x --minmax]
"""


def disc(mutations, t1_ancestors, t2_ancestors):
    """
    Compute the DISC distance between two trees stored as ancestor sets.
    :param mutations: the set of mutations over which to sum
    :param t1_ancestors: a dict storing the t1 ancestor sets of every mutation
    :param t2_ancestors: a dict storing the t2 ancestor sets of every mutation
    :return: the DISC distance between t1 and t2
    """
    m = len(mutations)
    if m < 2:
        return np.nan
    total = 0
    for i in range(m):
        for j in range(m):
            if i != j:
                t1_distinct_set = t1_ancestors[mutations[i]] - t1_ancestors[mutations[j]]
                t2_distinct_set = t2_ancestors[mutations[i]] - t2_ancestors[mutations[j]]

                total += jaccard(t1_distinct_set, t2_distinct_set)

    return total / (m * (m - 1))


def disc_intersection(t1, t2):
    """
    Compute the DISC distance between two trees, only summing over pairs of mutations in both trees.
    :param t1: a Newick string representation of tree 1
    :param t2: a Newick string representation of tree 2
    :return: the DISC intersection distance between tree 1 and tree 2
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    intersection = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))

    return disc(intersection, t1_ancestors, t2_ancestors)


def disc_union(t1, t2):
    """
    Compute the DISC distance between two trees, summing over all pairs of mutations.
    :param t1: a Newick string representation of tree 1
    :param t2: a Newick string representation of tree 2
    :return: the DISC union distance between tree 1 and tree 2
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    union = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))

    for u in union:
        if u not in t1_ancestors:
            t1_ancestors[u] = set()
        if u not in t2_ancestors:
            t2_ancestors[u] = set()

    return disc(union, t1_ancestors, t2_ancestors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find pairwise DISC distances of a set of trees.')
    parser.add_argument("inputFile")
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
    output_dir = Path("output") / "DISC" / package_name
    output_dir.mkdir(parents=True, exist_ok=True)

    trees = get_newicks(args.inputFile)

    if args.treePrint:
        for tree in trees:
            print(tree)

    distance_measure = disc_union if args.union else disc_intersection
    distance_matrix = [[None for x in range(len(trees))] for y in range(len(trees))]
    for i, t1 in enumerate(trees):
        for j, t2 in enumerate(trees):
            distance_matrix[i][j] = distance_measure(t1, t2)
            distance_matrix[j][i] = distance_matrix[i][j]

    output = get_matrix_string(distance_matrix)

    pairwise_filename = get_filename_with_new_extension(args.inputFile, ".png", "DISC_trees")
    pairwise_numpy_filename = get_filename_with_new_extension(args.inputFile, ".npy", "DISC_trees")
    pairwise_array = np.asarray(distance_matrix, dtype=float)
    np.save(output_dir / pairwise_numpy_filename, pairwise_array)
    save_heatmap(pairwise_array, output_dir / pairwise_filename, title=pairwise_filename,
                 cmap="RdBu", show_values=args.show)

    if args.consensusFile:
        consensus_trees = get_newicks(args.consensusFile)
        if len(consensus_trees) != 1:
            parser.error('--consensusFile deve conter exatamente uma árvore Newick.')
        comparison = [distance_measure(consensus_trees[0], tree) for tree in trees]
        result_array = np.asarray(comparison, dtype=float)
        consensus_filename = get_filename_with_new_extension(args.inputFile, ".png", "DISC_consensus")
        consensus_numpy_filename = get_filename_with_new_extension(args.inputFile, ".npy", "DISC_consensus")
        np.save(output_dir / consensus_numpy_filename, result_array)
        save_heatmap(result_array.reshape(1, -1), output_dir / consensus_filename,
                     title=consensus_filename, cmap="RdBu", show_values=args.show,
                     figsize=(12, 2.6))

        defined_values = result_array[np.isfinite(result_array)]
        sorted_array = np.sort(defined_values)
        sorted_filename = get_filename_with_new_extension(
            args.inputFile, ".png", "DISC_consensus_sort")
        if sorted_array.size:
            endpoint_labels = [f"{sorted_array[0]:.4f}", f"{sorted_array[-1]:.4f}"]
            sorted_title = f"{sorted_filename} (comparáveis: {sorted_array.size}; indefinidas: {result_array.size - sorted_array.size})"
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
