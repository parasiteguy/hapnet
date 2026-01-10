from __future__ import annotations

import argparse

from .io import read_fasta
from .haplotypes import build_haplotypes
from .network import build_mst_network
from .plot import plot_network
from .logs import write_logs
from .summary import write_summary


def main():
    parser = argparse.ArgumentParser(
        prog="hapnet",
        description="Build a population-aware haplotype network (MST) from an aligned FASTA file.",
    )

    parser.add_argument(
        "fasta",
        help="Aligned FASTA file. Population must be the last underscore-delimited token in each header.",
    )

    parser.add_argument(
        "--out",
        default="hapnet.png",
        help="Output image file (PNG, PDF, or SVG). Default: hapnet.png",
    )

    parser.add_argument(
        "--log-prefix",
        default="hapnet",
        help="Prefix for TSV log files. Default: hapnet",
    )

    args = parser.parse_args()

    # 1) Read FASTA
    records = read_fasta(args.fasta)

    # 2) Collapse sequences into haplotypes
    haplotypes, _ = build_haplotypes(records)

    # 3) Build MST network
    G = build_mst_network(haplotypes)

    # 4) Plot network
    plot_network(G, haplotypes, out=args.out)

    # 5) Write detailed logs
    write_logs(haplotypes, out_prefix=args.log_prefix)

    # 6) Write summary statistics
    write_summary(
        haplotypes,
        n_sequences=len(records),
        out_prefix=args.log_prefix,
    )

    print(f"Network written to: {args.out}")
    print(f"Logs written with prefix: {args.log_prefix}")
    print(f"Summary written to: {args.log_prefix}_summary.tsv")


if __name__ == "__main__":
    main()
