"""Lab 5 - Exercise 1: simulate gel electrophoresis of random DNA fragments.

Steps (as in the assignment):
  1. take a DNA sequence of 1000-3000 nt (FASTA from NCBI, or a random one)
  2. cut 10 random samples of 100-3000 bases from it
  3. store them in an array (a Python list)
  4. simulate how they migrate through the gel: short = fast, long = slow

Usage:
  python exercise1.py                      # random sequence, shows the gel
  python exercise1.py --fasta seq.fasta    # sequence downloaded from NCBI
  python exercise1.py --animate            # also saves gel_run.gif
  python exercise1.py --seed 42 --no-show  # reproducible, saves PNG only
"""
import argparse
import math
import random

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle

from dna_utils import load_sequence

N_SAMPLES = 10
MIN_FRAGMENT, MAX_FRAGMENT = 100, 3000

# --- migration model -------------------------------------------------------
# In a real gel the distance travelled is roughly linear in log10(size):
# bigger fragment -> more friction -> shorter distance.
# A fragment of SIZE_STUCK bases barely leaves the well; one of SIZE_FASTEST
# bases runs to the bottom of the gel.
SIZE_STUCK = 10_000
SIZE_FASTEST = 50

LADDER = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000]  # reference lane


def final_distance(length):
    """Fraction (0..1) of the gel a fragment covers at the end of the run."""
    num = math.log10(SIZE_STUCK) - math.log10(length)
    den = math.log10(SIZE_STUCK) - math.log10(SIZE_FASTEST)
    return max(0.0, min(1.0, num / den))


def position(length, t):
    """Distance travelled after time t (t from 0 to 1).
    Speed is constant per fragment and bigger for shorter fragments."""
    return final_distance(length) * t


# --- steps 2 and 3 ---------------------------------------------------------
def take_samples(sequence, n=N_SAMPLES, rng=random):
    """Cut n random pieces (100..3000 bases) out of the sequence.
    Returns (samples, positions); `samples` is the array of fragments."""
    samples, positions = [], []
    longest = min(MAX_FRAGMENT, len(sequence))
    for _ in range(n):
        length = rng.randint(MIN_FRAGMENT, longest)
        start = rng.randint(0, len(sequence) - length)
        samples.append(sequence[start:start + length])
        positions.append((start, start + length))
    return samples, positions


# --- step 4: drawing -------------------------------------------------------
def draw_gel(ax, lanes, t):
    """Draw the gel at time t. `lanes` = list of (label, [fragment lengths])."""
    ax.clear()
    ax.set_facecolor("#10151c")
    ax.set_xlim(0, len(lanes))
    ax.set_ylim(1.12, -0.05)          # y grows downwards: wells at the top
    ax.set_xticks([i + 0.5 for i in range(len(lanes))])
    ax.set_xticklabels([name for name, _ in lanes], fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks([])

    for i, (name, lengths) in enumerate(lanes):
        # the well
        ax.add_patch(Rectangle((i + 0.15, -0.04), 0.7, 0.04, color="#3a4452"))
        for L in lengths:
            y = position(L, t)
            brightness = 0.35 + 0.65 * min(L, 3000) / 3000   # more DNA = brighter
            ax.add_patch(Rectangle((i + 0.18, y), 0.64, 0.012,
                                   color=(0.4, 1.0, 0.7), alpha=brightness))
            if name == "Ladder":
                ax.text(i + 0.86, y + 0.006, str(L), color="white",
                        fontsize=7, ha="left", va="center")
    ax.set_title("Gel electrophoresis  (+ electrode at the bottom)", pad=28)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", help="FASTA file downloaded from NCBI")
    parser.add_argument("--seed", type=int, help="seed for reproducible results")
    parser.add_argument("--animate", action="store_true", help="save gel_run.gif")
    parser.add_argument("--no-show", action="store_true", help="do not open a window")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    sequence = load_sequence(args.fasta, args.seed)
    print(f"Sequence length: {len(sequence)} nt")

    samples, positions = take_samples(sequence, rng=rng)   # the array of fragments

    # Report: sorted by migration distance (fastest first)
    print(f"\n{'#':>2} {'start':>6} {'end':>6} {'length':>7} {'distance':>9}")
    order = sorted(range(N_SAMPLES), key=lambda i: -final_distance(len(samples[i])))
    for i in order:
        s, e = positions[i]
        print(f"{i + 1:>2} {s:>6} {e:>6} {len(samples[i]):>7} "
              f"{final_distance(len(samples[i])):>9.3f}")

    # One lane for the ladder + one lane per sample
    lanes = [("Ladder", LADDER)]
    lanes += [(f"S{i + 1}", [len(s)]) for i, s in enumerate(samples)]

    fig, ax = plt.subplots(figsize=(9, 7))

    draw_gel(ax, lanes, t=1.0)
    fig.savefig("gel_result.png", dpi=150, bbox_inches="tight")
    print("\nSaved gel_result.png")

    if args.animate:
        def update(frame):
            draw_gel(ax, lanes, t=frame / 59)
        anim = FuncAnimation(fig, update, frames=60, interval=50)
        anim.save("gel_run.gif", writer=PillowWriter(fps=20))
        print("Saved gel_run.gif")
        draw_gel(ax, lanes, t=1.0)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()