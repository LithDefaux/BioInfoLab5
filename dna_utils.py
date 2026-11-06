"""Shared helpers for Lab 5: loading a DNA sequence (FASTA from NCBI)."""
import random

MIN_LEN, MAX_LEN = 1000, 3000


def read_fasta(path):
    """Read a FASTA file (as downloaded from NCBI) and return the sequence
    as one uppercase string containing only A, C, G, T."""
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):   # skip header / blank lines
                continue
            seq_parts.append(line.upper())
    seq = "".join(seq_parts)
    return "".join(ch for ch in seq if ch in "ACGT")


def random_sequence(length=2000, seed=None):
    """Fallback if you have no FASTA file: a random DNA sequence."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def load_sequence(path=None, seed=None):
    """Load a sequence from `path` (FASTA) or generate one.
    The lab requires a length between 1000 and 3000 nucleotides."""
    if path:
        seq = read_fasta(path)
        if len(seq) > MAX_LEN:
            seq = seq[:MAX_LEN]                    # trim to the allowed range
        if len(seq) < MIN_LEN:
            raise ValueError(
                f"Sequence has {len(seq)} nt; the lab needs {MIN_LEN}-{MAX_LEN}."
            )
        return seq
    return random_sequence(random.Random(seed).randint(MIN_LEN, MAX_LEN), seed)
