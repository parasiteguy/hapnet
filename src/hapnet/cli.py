from __future__ import annotations

import argparse
import sys

from . import __version__
from .distance import AMBIGUOUS_CHARS
from .io import read_fasta
from .haplotypes import build_haplotypes
from .network import build_mst_network
from .plot import plot_network
from .logs import write_logs, write_run_metadata
from .summary import write_summary


def main():
    parser = argparse.ArgumentParser(
        prog="hapnet",
        description=(
            "Build a population-aware, MST-based haplotype graph from an aligned FASTA file. "
            "HapNet does not infer median-joining or statistical-parsimony networks."
        ),
    )

    parser.add_argument(
        "fasta",
        help=(
            "Aligned FASTA file. By default, population must be the last "
            "underscore-delimited token in each header."
        ),
    )

    parser.add_argument(
        "--out",
        default="hapnet.png",
        help="Output image file (PNG, PDF, or SVG). Default: hapnet.png",
    )

    parser.add_argument(
        "--log-prefix",
        default="hapnet",
        help="Prefix for TSV output files. Default: hapnet",
    )

    parser.add_argument(
        "--sep",
        default="_",
        help="Header separator used for parsing population/allele fields. Default: '_'",
    )

    parser.add_argument(
        "--phased",
        action="store_true",
        help=(
            "Parse phased diploid headers as individual_allele_population, "
            "e.g. Ind01_a_NK and Ind01_b_NK. HapNet does not infer phase; "
            "it only preserves already phased allele copies."
        ),
    )

    parser.add_argument(
        "--metadata",
        default=None,
        help=(
            "Optional tab-delimited metadata file with columns sequence_id, "
            "individual_id, population, and optional allele. Metadata values "
            "override header parsing."
        ),
    )

    parser.add_argument(
        "--ignore-ambiguous",
        action="store_true",
        help=(
            "Ignore N, ?, -, and . at pairwise-distance sites instead of "
            "treating them as literal character states."
        ),
    )

    parser.add_argument(
        "--show-counts-in-label",
        action="store_true",
        help="Show haplotype counts inside node labels.",
    )

    parser.add_argument(
        "--hide-labels",
        action="store_true",
        help=(
            "Hide haplotype labels inside network nodes. This is useful for "
            "cleaner publication-style figures when haplotype identities are "
            "documented in the output tables."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"hapnet {__version__}",
    )

    args = parser.parse_args()

    print(
        "HapNet constructs a minimum-spanning-tree-based haplotype graph; "
        "it does not infer median-joining, statistical-parsimony, or reticulate networks.",
        file=sys.stderr,
    )

    # 1) Read FASTA
    records = read_fasta(
        args.fasta,
        sep=args.sep,
        phased=args.phased,
        metadata=args.metadata,
    )

    # 2) Collapse sequences into haplotypes
    haplotypes, header_to_hap_index = build_haplotypes(records)

    # 3) Build MST network
    gap_chars = AMBIGUOUS_CHARS if args.ignore_ambiguous else None
    G = build_mst_network(haplotypes, gap_chars=gap_chars)

    # 4) Plot network
    plot_network(
        G,
        haplotypes,
        out=args.out,
        show_labels=not args.hide_labels,
        show_counts_in_label=args.show_counts_in_label,
    )

    # 5) Write detailed logs
    write_logs(
        haplotypes,
        out_prefix=args.log_prefix,
        records=records,
        header_to_hap_index=header_to_hap_index,
        phased=args.phased,
    )

    # 6) Write summary statistics
    write_summary(
        haplotypes,
        n_sequences=len(records),
        out_prefix=args.log_prefix,
    )

    # 7) Write run metadata for reproducibility
    write_run_metadata(
        out_prefix=args.log_prefix,
        hapnet_version=__version__,
        input_file=args.fasta,
        output_file=args.out,
        n_sequences=len(records),
        n_haplotypes=len(haplotypes),
        phased=args.phased,
        metadata_file=args.metadata,
        ignore_ambiguous=args.ignore_ambiguous,
    )

    print(f"Network written to: {args.out}")
    print(f"Logs written with prefix: {args.log_prefix}")
    print(f"Summary written to: {args.log_prefix}_summary.tsv")
    print(f"Run metadata written to: {args.log_prefix}_run_metadata.tsv")

    if args.hide_labels:
        print("Haplotype labels were hidden in the network figure.")
    else:
        print("Haplotype labels were shown in the network figure.")

    if args.phased:
        print(f"Individual genotype table written to: {args.log_prefix}_individual_genotypes.tsv")


if __name__ == "__main__":
    main()
