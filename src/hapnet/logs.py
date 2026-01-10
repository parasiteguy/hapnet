from __future__ import annotations

from typing import List
import csv
from .haplotypes import Haplotype


def write_logs(haplotypes: List[Haplotype], out_prefix: str = "hapnet"):
    """
    Writes:
      - {prefix}_haplotypes.tsv: hap_id, n_total, populations, counts_by_pop, sequence
      - {prefix}_membership.tsv: hap_id, member_header
      - {prefix}_shared_haplotypes.tsv: hap_id, populations, counts_by_pop
    """
    hap_path = f"{out_prefix}_haplotypes.tsv"
    mem_path = f"{out_prefix}_membership.tsv"
    shared_path = f"{out_prefix}_shared_haplotypes.tsv"

    # haplotypes.tsv
    with open(hap_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["hap_id", "n_total", "populations", "counts_by_pop", "sequence"])
        for h in haplotypes:
            pops = ",".join(sorted(h.counts_by_pop.keys()))
            counts = ",".join([f"{p}:{h.counts_by_pop[p]}" for p in sorted(h.counts_by_pop.keys())])
            w.writerow([h.hap_id, h.n_total, pops, counts, h.sequence])

    # membership.tsv
    with open(mem_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["hap_id", "member_header"])
        for h in haplotypes:
            for m in h.members:
                w.writerow([h.hap_id, m])

    # shared_haplotypes.tsv
    with open(shared_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["hap_id", "populations", "counts_by_pop"])
        for h in haplotypes:
            if len(h.counts_by_pop) > 1:
                pops = ",".join(sorted(h.counts_by_pop.keys()))
                counts = ",".join([f"{p}:{h.counts_by_pop[p]}" for p in sorted(h.counts_by_pop.keys())])
                w.writerow([h.hap_id, pops, counts])
