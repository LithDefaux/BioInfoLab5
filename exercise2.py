"""Lab 5 - Exercise 2: detect repetitions in a DNA sequence.

Interpretation used (tandem repeats):
  * the repeated unit has 3-5 bases           (UNIT_MIN..UNIT_MAX)
  * the whole repeated region is 6-10 bases   (REGION_MIN..REGION_MAX)
  and the unit appears at least twice in a row. Examples:
        ATGATG      -> unit ATG, region length 6
        CAGCAGCAG   -> unit CAG, region length 9
        GATTCGATTC  -> unit GATTC, region length 10

Usage:
  python exercise2.py                    # random sequence
  python exercise2.py --fasta seq.fasta  # sequence from NCBI
  python exercise2.py --seed 1
"""
import argparse

from dna_utils import load_sequence

UNIT_MIN, UNIT_MAX = 3, 5
REGION_MIN, REGION_MAX = 6, 10


def is_primitive(unit):
    """True if the unit is not itself a repeat of something shorter
    (so 'AAA' is rejected: it is really a 1-base repeat)."""
    n = len(unit)
    for d in range(1, n):
        if n % d == 0 and unit == unit[:d] * (n // d):
            return False
    return True


def find_tandem_repeats(seq, unit_min=UNIT_MIN, unit_max=UNIT_MAX,
                        region_min=REGION_MIN, region_max=REGION_MAX):
    """Return a list of (start, end, unit, region_length) for every maximal
    tandem repeat whose region length is in [region_min, region_max].

    Idea: a region of period u satisfies seq[i] == seq[i + u] along the
    whole region. So for each u we count runs of consecutive matches
    seq[i] == seq[i + u]; a run of m matches gives a region of m + u bases.
    """
    found = []
    n = len(seq)
    for u in range(unit_min, unit_max + 1):
        i = 0
        while i + u < n:
            if seq[i] != seq[i + u]:
                i += 1
                continue
            j = i
            while j + u < n and seq[j] == seq[j + u]:   # extend the run
                j += 1
            start, end = i, j + u                        # region = seq[start:end]
            length = end - start
            unit = seq[start:start + u]
            # need at least 2 complete copies of the unit, otherwise it is
            # not a repetition (e.g. 'AGCCCA' only has the 1st base repeated)
            if (region_min <= length <= region_max and length >= 2 * u
                    and is_primitive(unit)):
                found.append((start, end, unit, length))
            i = j + 1                                    # continue after the run
    found.sort()
    return found


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", help="FASTA file downloaded from NCBI")
    parser.add_argument("--seed", type=int, help="seed for the random sequence")
    args = parser.parse_args()

    seq = load_sequence(args.fasta, args.seed)
    print(f"Sequence length: {len(seq)} nt")
    print(f"Looking for tandem repeats: unit {UNIT_MIN}-{UNIT_MAX} bases, "
          f"region {REGION_MIN}-{REGION_MAX} bases\n")

    repeats = find_tandem_repeats(seq)
    if not repeats:
        print("No repeats found.")
        return

    print(f"{'start':>6} {'end':>6} {'unit':<6} {'copies':>7}  region")
    for start, end, unit, length in repeats:
        copies = length / len(unit)
        print(f"{start + 1:>6} {end:>6} {unit:<6} {copies:>7.1f}  {seq[start:end]}")
    print(f"\nTotal: {len(repeats)} repeats (positions are 1-based).")


if __name__ == "__main__":
    main()