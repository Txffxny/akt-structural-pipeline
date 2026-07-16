"""
03_pocket_conservation_analysis.py

Aligns AKT2 and AKT3 against AKT1 and reports whether the residues known
to form the MK-2206 / miransertib allosteric pocket (PH-kinase domain
interface) are conserved across isoforms.

Pocket residue list (AKT1 numbering, 1-indexed) is taken from the
MK-2206/Inhibitor VIII structural studies:
  Trp-80, Asn-53, Glu-59, Leu-78          (common contact / Trp80 environment)
  Leu-264, Tyr-272                         (shared contacts, both ligands)
  Arg-273, Asp-274, Asp-292, Cys-296       (Inhibitor VIII-specific region)
  Val-201, Val-270                         (MK-2206-specific region)

Output: printed table + data/pocket_conservation.csv
"""

import os
import csv
from Bio.Align import PairwiseAligner, substitution_matrices

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

POCKET_RESIDUES_AKT1 = {
    53: "N", 59: "E", 78: "L", 80: "W",
    201: "V", 264: "L", 270: "V", 272: "Y",
    273: "R", 274: "D", 292: "D", 296: "C",
}


def read_fasta(path):
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    header = lines[0].lstrip(">")
    seq = "".join(lines[1:])
    return header, seq


def build_aligner():
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    return aligner


def map_reference_positions(ref_seq, other_seq, aligner):
    """
    Align ref_seq (AKT1) to other_seq (AKT2/AKT3) and return a dict mapping
    1-indexed AKT1 position -> residue in other_seq at the aligned column
    (or '-' if that column is a gap in other_seq).
    """
    alignment = aligner.align(ref_seq, other_seq)[0]
    aligned_ref, aligned_other = alignment[0], alignment[1]

    mapping = {}
    ref_pos = 0  # 1-indexed position in ungapped ref_seq
    for a_ref, a_other in zip(aligned_ref, aligned_other):
        if a_ref != "-":
            ref_pos += 1
            mapping[ref_pos] = a_other
    return mapping


def main():
    akt1_path = os.path.join(DATA_DIR, "AKT1.fasta")
    if not os.path.exists(akt1_path):
        print("AKT1.fasta not found -- run step 01 first.")
        return

    _, akt1_seq = read_fasta(akt1_path)

    isoform_maps = {}
    for gene in ["AKT2", "AKT3"]:
        path = os.path.join(DATA_DIR, f"{gene}.fasta")
        if not os.path.exists(path):
            print(f"{gene}.fasta not found -- run step 01 first.")
            continue
        _, other_seq = read_fasta(path)
        aligner = build_aligner()
        isoform_maps[gene] = map_reference_positions(akt1_seq, other_seq, aligner)

    rows = []
    print(f"{'Pos':<6}{'AKT1':<6}{'AKT2':<6}{'AKT3':<6}{'Conserved?'}")
    for pos, akt1_res in sorted(POCKET_RESIDUES_AKT1.items()):
        actual_akt1 = akt1_seq[pos - 1] if pos <= len(akt1_seq) else "?"
        akt2_res = isoform_maps.get("AKT2", {}).get(pos, "?")
        akt3_res = isoform_maps.get("AKT3", {}).get(pos, "?")
        conserved = (actual_akt1 == akt2_res == akt3_res)
        print(f"{pos:<6}{actual_akt1:<6}{akt2_res:<6}{akt3_res:<6}{'yes' if conserved else 'NO'}")
        rows.append([pos, actual_akt1, akt2_res, akt3_res, conserved])

    csv_path = os.path.join(DATA_DIR, "pocket_conservation.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["akt1_position", "akt1_residue", "akt2_residue", "akt3_residue", "conserved"])
        writer.writerows(rows)

    print(f"\nSaved: {csv_path}")
    print(
        "\nInterpretation: any 'NO' rows are candidate explanations for "
        "isoform-selective potency (e.g. why AKT3 is ~5-8x less sensitive "
        "to MK-2206/miransertib than AKT1/AKT2) -- these are the residues "
        "worth highlighting on the AlphaFold structure figure."
    )


if __name__ == "__main__":
    main()
