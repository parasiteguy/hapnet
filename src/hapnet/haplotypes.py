from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import defaultdict

from .io import SequenceRecord


@dataclass
class Haplotype:
    hap_id: str
    sequence: str
    n_total: int
    counts_by_pop: Dict[str, int]
    members: List[str]  # headers or individuals (we'll store headers for traceability)


def build_haplotypes(records: List[SequenceRecord]) -> Tuple[List[Haplotype], Dict[str, int]]:
    """
    Collapse identical sequences into haplotypes.
    Returns:
      - list of Haplotype objects (H1..Hk)
      - dict mapping header -> haplotype index
    """
    seq_to_members = defaultdict(list)
    seq_to_popcounts = defaultdict(lambda: defaultdict(int))

    for r in records:
        seq_to_members[r.sequence].append(r.header)
        seq_to_popcounts[r.sequence][r.population] += 1

    # Stable ordering: sort by (descending frequency, then sequence) for reproducibility
    seqs_sorted = sorted(seq_to_members.keys(), key=lambda s: (-len(seq_to_members[s]), s))

    haplotypes: List[Haplotype] = []
    header_to_hap_index: Dict[str, int] = {}

    for i, seq in enumerate(seqs_sorted, start=1):
        hap_id = f"H{i}"
        members = seq_to_members[seq]
        counts_by_pop = dict(sorted(seq_to_popcounts[seq].items(), key=lambda x: x[0]))
        hap = Haplotype(
            hap_id=hap_id,
            sequence=seq,
            n_total=len(members),
            counts_by_pop=counts_by_pop,
            members=sorted(members),
        )
        haplotypes.append(hap)

        for h in members:
            header_to_hap_index[h] = i - 1

    return haplotypes, header_to_hap_index
