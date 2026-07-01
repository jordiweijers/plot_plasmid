import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering, leaves_list


def canonical_pair(a: int, b: int) -> Tuple[int, int]:
    """
    Sort a pair of contig IDs so the smaller one is always first, fixing a single BLASTP direction per unordered pair.
    Args:
        a (int): The first contig ID.
        b (int): The second contig ID.
    Returns:
        Tuple[int, int]: (a, b) sorted so the smaller ID comes first.
    """
    return (a, b) if a <= b else (b, a)


def swap_query_subject(df: pd.DataFrame) -> pd.DataFrame:
    """
    Swap the query/subject roles in a parsed BLASTP DataFrame, so a query=A vs subject=B result can be reused as query=B vs subject=A.
    Args:
        df (pd.DataFrame): A parsed BLASTP outfmt-6 DataFrame (see parse.parse_blast_output).
    Returns:
        pd.DataFrame: A copy of df with qseqid<->sseqid and qstart/qend<->sstart/send swapped as coordinate pairs. pident, length, mismatch, gapopen, evalue and bitscore are left unchanged.
    """
    swapped = df.rename(columns={
        "qseqid": "sseqid", "sseqid": "qseqid",
        "qstart": "sstart", "sstart": "qstart",
        "qend": "send", "send": "qend",
    })
    return swapped[df.columns]


def score_pair(df: pd.DataFrame, min_identity: float) -> float:
    """
    Compute a similarity score for one contig pair from its parsed BLASTP DataFrame.
    Args:
        df (pd.DataFrame): A parsed BLASTP outfmt-6 DataFrame for the pair (see parse.parse_blast_output).
        min_identity (float): Minimum percent identity for a hit to count towards the score.
    Returns:
        float: The sum of bitscore across rows with pident >= min_identity.
    Raises:
        ValueError: If df is empty or no rows have pident >= min_identity.
    """
    if df.empty:
        raise ValueError("Cannot score an empty BLASTP result.")
    filtered = df[df["pident"] >= min_identity]
    if filtered.empty:
        raise ValueError(f"No BLASTP hits with pident >= {min_identity}.")
    return float(filtered["bitscore"].sum())


def order_by_clustering(
    contig_ids: List[int],
    pair_scores: Dict[Tuple[int, int], float],
    method: str = "average",
) -> List[int]:
    """
    Order contigs by hierarchical clustering with optimal leaf ordering on their pairwise similarity scores, so that neighboring contigs in the result are as similar as possible. No contig is anchored to a fixed position.
    Args:
        contig_ids (List[int]): The contigs to order.
        pair_scores (Dict[Tuple[int, int], float]): A dictionary mapping canonical_pair(a, b) to a similarity score, for every unordered pair among contig_ids. Pairs with no similarity score (e.g. no BLASTP hits above the identity threshold) should be omitted and are treated as similarity 0.
        method (str): The scipy linkage method to use.
    Returns:
        List[int]: contig_ids reordered so neighboring contigs are as similar as possible.
    """
    n = len(contig_ids)
    if n < 2:
        return list(contig_ids)

    index = {contig_id: i for i, contig_id in enumerate(contig_ids)}
    similarity = np.zeros((n, n))
    for (a, b), score in pair_scores.items():
        i, j = index[a], index[b]
        similarity[i, j] = score
        similarity[j, i] = score

    distance = similarity.max() - similarity
    np.fill_diagonal(distance, 0.0)

    condensed = squareform(distance)
    Z = linkage(condensed, method=method)
    Z_ordered = optimal_leaf_ordering(Z, condensed)
    order = leaves_list(Z_ordered)
    return [contig_ids[i] for i in order]
