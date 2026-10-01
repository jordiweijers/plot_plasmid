import numpy as np
from typing import List

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering, leaves_list

from plot_plasmid.models import Contig


def compute_presence_matrix(contigs: List[Contig]) -> np.ndarray:
    """
    Compute a presence/absence matrix of protein families per contig.
    Args:
        contigs (List[Contig]): The contigs to include.
    Returns:
        np.ndarray: A 2D boolean array of shape (n_contigs, n_families) that is True where the contig has at least one feature of the protein family. Features without a protein family are ignored.
    """
    family_sets = [{feature.protein_family for feature in contig.features if feature.protein_family is not None} for contig in contigs]
    column = {family: i for i, family in enumerate(sorted(set().union(*family_sets)))}
    presence_matrix = np.zeros((len(contigs), len(column)), dtype=bool)
    for row, families in enumerate(family_sets):
        presence_matrix[row, [column[family] for family in families]] = True
    return presence_matrix


def compute_jaccard_matrix(presence_matrix: np.ndarray) -> np.ndarray:
    """
    Compute a Jaccard similarity matrix between contigs based on their protein family presence/absence.
    Args:
        presence_matrix (np.ndarray): A 2D boolean array of shape (n_contigs, n_families), as returned by compute_presence_matrix.
    Returns:
        np.ndarray: A 2D array of shape (n_contigs, n_contigs) with Jaccard similarity values between 0 and 1. Two contigs without any protein family have similarity 0.
    """
    intersection = np.logical_and(presence_matrix[:, None, :], presence_matrix[None, :, :]).sum(axis=2)
    union = np.logical_or(presence_matrix[:, None, :], presence_matrix[None, :, :]).sum(axis=2)
    return np.divide(intersection, union, out=np.zeros_like(intersection, dtype=float), where=union > 0)


def order_by_clustering(contigs: List[Contig], method: str = "average") -> List[Contig]:
    """
    Order contigs by hierarchical clustering with optimal leaf ordering on their pairwise Jaccard similarity, so that neighboring contigs in the result are as similar as possible. No contig is anchored to a fixed position.
    Args:
        contigs (List[Contig]): The contigs to order.
        method (str): The scipy linkage method to use.
    Returns:
        List[Contig]: The contigs reordered so neighboring contigs are as similar as possible.
    """
    if len(contigs) < 2:
        return list(contigs)

    presence_matrix = compute_presence_matrix(contigs)
    distance = 1.0 - compute_jaccard_matrix(presence_matrix)
    np.fill_diagonal(distance, 0.0)

    condensed = squareform(distance)
    Z = linkage(condensed, method=method)
    Z_ordered = optimal_leaf_ordering(Z, condensed)
    order = leaves_list(Z_ordered)
    return [contigs[i] for i in order]
