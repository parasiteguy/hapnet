# hapnet

**hapnet** is a Python package for building **population-aware haplotype networks** from aligned FASTA files.  
It integrates haplotype inference, population membership, and visualization into a reproducible, command-line workflow.

hapnet constructs minimum spanning tree (MST) haplotype networks in which node size reflects haplotype frequency and shared haplotypes are represented as population-colored pie charts.

---

## Features

- Aligned FASTA input
- Population parsed directly from FASTA headers
- Minimum spanning tree (MST) haplotype networks
- Haplotype nodes sized by frequency
- Population-aware pie charts for shared haplotypes
- Mutation tick marks on network edges
- Reproducible TSV logs and summary statistics

---

## Installation

```bash
pip install hapnet

Python 3.9 or newer is required

---

USAGE

hapnet input.fasta --out network.png --log-prefix run1

This command produces:
network.png - haplotype network visualization
run1_haplotypes.tsv - haplotype definitions and frequencies
run1_membership.tsv - sequence-to-haplotype membership
run1_shared_haplotypes.tsv - haplotypes shared among populations
run1_summary.tsv - summary statistics (total haplotypes, shared vs private, etc.)

FASTA HEADER FORMAT

Population identity must be encoded as the final underscore-delimited token in each FASTA header.

Examples:
Ind1_Pop1
Ind7_SiteA_2019_Pop3
MN605578_Pneocaeca_RI

In all cases, the population is interpreted as the final token (Pop1, Pop3, RI).
Sequences must be aligned and of equal length.

SCALABILITY and VISUALIZATION NOTES

Hapnet scales computationally with the number of distinct haplotypes, not the number of input sequences. Static network visualizations are typically most interpretable for datasets containing tens to low hundreds of haplotypes, depending on network structure and frequency distribution. All haplotypes and population memberships are retained in the output logs regardless of visualization complexity. 

LICENSE
Hapnet is released under the MIT License 

