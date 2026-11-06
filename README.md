# Lab 5 – DNA Sequence Analysis

Two small Python programs that work on a DNA sequence (1000–3000 nucleotides):

1. **Gel electrophoresis simulation** – cuts 10 random fragments from the sequence and simulates how they migrate through a gel.
2. **Tandem repeat detection** – finds short repeated units inside the sequence.

## Project structure

```
lab5/
├── dna_utils.py       # loads a FASTA file (or generates a random sequence)
├── exercise1.py       # gel electrophoresis simulation
├── exercise2.py       # tandem repeat finder
├── requirements.txt
└── README.md
```

## Setup

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Getting a sequence from NCBI

1. Open <https://www.ncbi.nlm.nih.gov/nuccore> and pick any nucleotide record.
2. Click **Send to → File → Format: FASTA → Create File**.
3. Save it in this folder, e.g. as `sequence.fasta`.

The sequence must have at least 1000 nucleotides. If it is longer than 3000, only the first 3000 are used.

If no FASTA file is given, both programs generate a random sequence instead.

## Exercise 1 – Gel electrophoresis

```bash
python exercise1.py --fasta sequence.fasta
```

Options:

| Option | Description |
|---|---|
| `--fasta FILE` | FASTA file with the DNA sequence |
| `--seed N` | Seed for reproducible results |
| `--animate` | Also save the run as `gel_run.gif` |
| `--no-show` | Do not open a window, only save files |

**What it does**

- Takes 10 random samples of 100–3000 bases from the sequence and stores them in a list.
- Simulates migration: short fragments meet less friction and travel farther, long fragments travel less.
- The distance depends on `log10(length)`, which is how real gels behave.
- Draws the gel with a ladder lane (reference sizes) and one lane per sample.
- Prints a table of the fragments sorted from fastest to slowest and saves `gel_result.png`.

## Exercise 2 – Tandem repeats

```bash
python exercise2.py --fasta sequence.fasta
```

**What it detects**

A tandem repeat is a short unit repeated back to back:

- the unit has **3–5 bases**
- the whole repeated region is **6–10 bases**
- the unit appears at least twice in a row

| Region | Unit | Copies |
|---|---|---|
| `ATGATG` | `ATG` | 2 |
| `CAGCAGCAG` | `CAG` | 3 |
| `GATTCGATTC` | `GATTC` | 2 |

Units that are themselves repeats of something shorter (for example `AAA`) are ignored. Positions in the output are 1-based.

**How it works**

For each unit length `u` (3, 4, 5) the program scans the sequence and counts runs of consecutive positions where `seq[i] == seq[i + u]`. A run of `m` matches gives a repeated region of `m + u` bases. Regions whose length is between 6 and 10 are reported.

The limits can be changed with the constants at the top of `exercise2.py` (`UNIT_MIN`, `UNIT_MAX`, `REGION_MIN`, `REGION_MAX`).

## Requirements

- Python 3.8+
- matplotlib
- pillow (for the optional GIF animation)
