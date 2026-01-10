from __future__ import annotations

def hamming_distance(a: str, b: str, gap_chars: set[str] | None = None) -> int:
    """
    Hamming distance across aligned sequences.
    gap_chars: characters to skip if either sequence has that character at a site.
              If None, we count everything (including '-' and 'N') as literal states.
    """
    if len(a) != len(b):
        raise ValueError("Sequences must be equal length (aligned).")

    d = 0
    if gap_chars is None:
        for x, y in zip(a, b):
            if x != y:
                d += 1
        return d

    for x, y in zip(a, b):
        if x in gap_chars or y in gap_chars:
            continue
        if x != y:
            d += 1
    return d
