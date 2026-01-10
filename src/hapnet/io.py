from __future__ import annotations

from dataclasses import dataclass
from typing import List
from Bio import SeqIO


@dataclass(frozen=True)
class SequenceRecord:
    header: str
    individual: str
    population: str
    sequence: str


def parse_header(header: str, sep: str = "_") -> tuple[str, str]:
    """
    Parse FASTA header using the rule:
      population = last underscore-delimited token
      individual = everything before the last token (joined by sep)
    Example:
      Ind7_SiteA_2019_Pop3 -> individual="Ind7_SiteA_2019", population="Pop3"
    """
    parts = header.split(sep)
    if len(parts) < 2:
        raise ValueError(
            f"Header '{header}' does not contain separator '{sep}'. "
            "Expected format like 'Ind1_Pop1' where population is last token."
        )
    pop = parts[-1].strip()
    indiv = sep.join(parts[:-1]).strip()
    if not pop:
        raise ValueError(f"Header '{header}' has an empty population token.")
    if not indiv:
        raise ValueError(f"Header '{header}' has an empty individual token.")
    return indiv, pop


def read_fasta(path: str, sep: str = "_") -> List[SequenceRecord]:
    """
    Read aligned FASTA. Sequences are uppercased and whitespace removed.
    """
    records: List[SequenceRecord] = []
    for rec in SeqIO.parse(path, "fasta"):
        header = rec.id  # rec.description can contain spaces; id is safer default
        indiv, pop = parse_header(header, sep=sep)
        seq = str(rec.seq).upper().replace(" ", "").replace("\n", "")
        records.append(SequenceRecord(header=header, individual=indiv, population=pop, sequence=seq))

    if not records:
        raise ValueError(f"No FASTA records found in '{path}'.")

    # Basic alignment sanity check
    L = len(records[0].sequence)
    for r in records[1:]:
        if len(r.sequence) != L:
            raise ValueError(
                "FASTA must be an aligned set of sequences with equal length. "
                f"Found lengths {L} and {len(r.sequence)} (e.g., record '{r.header}')."
            )

    return records
