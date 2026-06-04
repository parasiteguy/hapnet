Simulated phased diploid FASTA dataset for HapNet

Files:
- simulated_phased_5pop_50seq.fasta
  Aligned FASTA file with 50 sequences, representing 25 diploid individuals.
  Header format: individual_allele_population.
  Example: >MAInd01_a_MA

- simulated_phased_5pop_50seq_metadata.tsv
  Optional metadata file with sequence_id, individual_id, allele, and population.

- simulated_phased_5pop_50seq_truth.tsv
  Simulation truth table showing the intended haplotype assignment for each sequence.

- simulated_phased_5pop_50seq_expected_haplotype_counts.tsv
  Expected haplotype composition by population before running HapNet.

Dataset design:
- 5 populations: MA, RI, CT, NY, ME
- 5 diploid individuals per population
- 2 phased allele copies per individual
- 50 total sequences
- 320 bp aligned DNA sequences
- 8 simulated haplotypes
- H1 and H2 are broadly shared across populations
- H3-H7 are more population-associated
- H8 is a divergent rare haplotype concentrated in ME

Example HapNet commands:

Header-parsed phased mode:
hapnet simulated_phased_5pop_50seq.fasta --phased --out simulated_phased_network.svg --log-prefix simulated_phased --hide-labels

Metadata-based phased mode:
hapnet simulated_phased_5pop_50seq.fasta --metadata simulated_phased_5pop_50seq_metadata.tsv --phased --out simulated_phased_network.svg --log-prefix simulated_phased --hide-labels
