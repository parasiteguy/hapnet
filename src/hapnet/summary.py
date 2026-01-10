from __future__ import annotations

from typing import List
from pathlib import Path

from .haplotypes import Haplotype


def write_summary(
    haplotypes: List[Haplotype],
    *,
    n_sequences: int,
    out_prefix: str = "hapnet",
):
    """
    Write a summary TSV describing the haplotype network.

    Parameters
    ----------
    haplotypes : list of Haplotype
        Haplotype objects produced by build_haplotypes().
    n_sequences : int
        Total number of input sequences.
    out_prefix : str
        Prefix for output file (e.g., run1 -> run1_summary.tsv).
    """

    n_haplotypes = len(haplotypes)

    # Shared haplotypes = present in >1 population
    n_shared = sum(1 for h in haplotypes if len(h.counts_by_pop) > 1)

    # Private haplotypes = present in exactly one population
    n_private = n_haplotypes - n_shared

    # Largest haplotype size (max number of individuals)
    max_haplotype_size = max((h.n_total for h in haplotypes), default=0)

    out_path = Path(f"{out_prefix}_summary.tsv")

    with out_path.open("w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"n_sequences\t{n_sequences}\n")
        fh.write(f"n_haplotypes\t{n_haplotypes}\n")
        fh.write(f"n_shared_haplotypes\t{n_shared}\n")
        fh.write(f"n_private_haplotypes\t{n_private}\n")
        fh.write(f"max_haplotype_size\t{max_haplotype_size}\n")
