"""
02_generate_w80a_mutants.py

Reads AKT1.fasta and AKT2.fasta (from step 01), generates the W80A point
mutant of each, and writes ColabFold-ready FASTA files:

  data/AKT1_WT.fasta      data/AKT1_W80A.fasta
  data/AKT2_WT.fasta      data/AKT2_W80A.fasta
  data/all_variants.fasta   (combined, for a single batch ColabFold job)

Position convention: 1-indexed, Met1 = position 1 (matches the numbering
used in the MK-2206/AKT1 structural literature, e.g. Trp-80).
"""

import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

MUTATION_POSITION = 80          # 1-indexed
WT_RESIDUE = "W"
MUT_RESIDUE = "A"

ISOFORMS_TO_MUTATE = ["AKT1", "AKT2"]  # extend to AKT3 if you decide to test it


def read_fasta(path):
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    header = lines[0].lstrip(">")
    seq = "".join(lines[1:])
    return header, seq


def write_fasta(path, header, seq):
    with open(path, "w") as fh:
        fh.write(f">{header}\n{seq}\n")


def make_mutant(seq, position, wt_residue, mut_residue, label):
    idx = position - 1
    actual = seq[idx]
    if actual != wt_residue:
        raise ValueError(
            f"{label}: residue at position {position} is '{actual}', "
            f"expected '{wt_residue}'. Check numbering (canonical vs. "
            f"isoform-2, Met1 offset, etc.) before proceeding."
        )
    return seq[:idx] + mut_residue + seq[idx + 1:]


def main():
    combined_records = []

    for gene in ISOFORMS_TO_MUTATE:
        fasta_path = os.path.join(DATA_DIR, f"{gene}.fasta")
        if not os.path.exists(fasta_path):
            print(f"Skipping {gene}: {fasta_path} not found. Run step 01 first.")
            continue

        header, wt_seq = read_fasta(fasta_path)

        mut_seq = make_mutant(
            wt_seq, MUTATION_POSITION, WT_RESIDUE, MUT_RESIDUE, label=gene
        )

        wt_header = f"{gene}_WT"
        mut_header = f"{gene}_W80A"

        write_fasta(os.path.join(DATA_DIR, f"{gene}_WT.fasta"), wt_header, wt_seq)
        write_fasta(os.path.join(DATA_DIR, f"{gene}_W80A.fasta"), mut_header, mut_seq)

        combined_records.append((wt_header, wt_seq))
        combined_records.append((mut_header, mut_seq))

        print(f"{gene}: WT and W80A FASTA written ({len(wt_seq)} aa each).")

    combined_path = os.path.join(DATA_DIR, "all_variants.fasta")
    with open(combined_path, "w") as fh:
        for header, seq in combined_records:
            fh.write(f">{header}\n{seq}\n")

    print(f"\nCombined batch file for ColabFold: {combined_path}")
    print(
        "\nNext step: submit data/all_variants.fasta (or individual files) to "
        "ColabFold (https://github.com/sokrypton/ColabFold -- the standard "
        "'AlphaFold2.ipynb' notebook on Google Colab). Each header becomes the "
        "job name, so WT and W80A predictions come back clearly labeled."
    )


if __name__ == "__main__":
    main()
