from __future__ import annotations

from typing import List
import networkx as nx

from .haplotypes import Haplotype
from .distance import hamming_distance


def build_complete_graph(haplotypes: List[Haplotype]) -> nx.Graph:
    G = nx.Graph()
    for i, h in enumerate(haplotypes):
        G.add_node(i, hap_id=h.hap_id)

    for i in range(len(haplotypes)):
        for j in range(i + 1, len(haplotypes)):
            d = hamming_distance(haplotypes[i].sequence, haplotypes[j].sequence, gap_chars=None)
            G.add_edge(i, j, weight=d)

    return G


def build_mst_network(haplotypes: List[Haplotype]) -> nx.Graph:
    """
    MST on the complete haplotype graph using Hamming distances as weights.
    Returns an undirected NetworkX graph with edge attribute 'weight' (mutational steps).
    """
    complete = build_complete_graph(haplotypes)
    mst = nx.minimum_spanning_tree(complete, weight="weight")
    return mst
