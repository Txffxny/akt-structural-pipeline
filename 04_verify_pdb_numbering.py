"""
04_verify_pdb_numbering.py

Downloads the AKT1 + Inhibitor VIII co-crystal structure (PDB 3O96),
extracts its residue numbering directly from the ATOM records, and checks
it against the UniProt-derived AKT1 sequence (data/AKT1.fasta) at the two
positions we care about: 80 (Trp, expected -- sanity anchor) and 59
(flagged as Gln in UniProt, but expected Glu based on the MK-2206 docking
literature).

Two independent checks are run:
  1. Direct lookup -- is there a residue with author-numbering 59/80 in
     the structure, and what is it?
  2. Alignment-based lookup -- align the resolved crystal sequence to the
     full UniProt sequence and see what maps to UniProt positions 59/80,
     regardless of the structure's internal numbering. This is the one
     that matters if the two disagree.

Requires internet access to files.rcsb.org (run on your own machine, not
inside a locked-down sandbox).
"""

import os
import requests
from Bio.Align import PairwiseAligner, substitution_matrices

PDB_ID = "3O96"
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
POSITIONS_OF_INTEREST = [59, 80]

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def download_pdb(pdb_id: str) -> str:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    out_path = os.path.join(DATA_DIR, f"{pdb_id}.pdb")
    with open(out_path, "w") as fh:
        fh.write(resp.text)
    return out_path


def parse_atom_records(pdb_path: str):
    """
    Parse ATOM records directly (no Bio.PDB structure object needed) and
    return {chain_id: {resseq: one_letter_aa}} using only the CA atom of
    each standard amino acid residue (skips waters, ligands, alt locs
    beyond the first).
    """
    chains = {}
    seen_altloc = set()

    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            resname = line[17:20].strip()
            if resname not in THREE_TO_ONE:
                continue

            chain_id = line[21].strip()
            resseq = int(line[22:26])
            altloc = line[16].strip()

            key = (chain_id, resseq)
            if altloc not in ("", "A") or key in seen_altloc:
                continue
            seen_altloc.add(key)

            chains.setdefault(chain_id, {})[resseq] = THREE_TO_ONE[resname]

    return chains


def read_fasta(path):
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    return "".join(lines[1:])


def direct_lookup(chain_residues, positions):
    print("  Direct lookup (by author residue number):")
    for pos in positions:
        aa = chain_residues.get(pos, None)
        print(f"    position {pos}: {'no residue at this number' if aa is None else aa}")


def alignment_lookup(chain_residues, uniprot_seq, positions):
    resseqs_sorted = sorted(chain_residues.keys())
    pdb_seq = "".join(chain_residues[r] for r in resseqs_sorted)

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"

    alignment = aligner.align(uniprot_seq, pdb_seq)[0]
    aligned_uniprot, aligned_pdb = alignment[0], alignment[1]

    uniprot_pos = 0
    pdb_idx = 0
    mapping = {}
    for a_uni, a_pdb in zip(aligned_uniprot, aligned_pdb):
        if a_uni != "-":
            uniprot_pos += 1
        if a_pdb != "-":
            pdb_idx += 1
            if a_uni != "-":
                mapping[uniprot_pos] = (resseqs_sorted[pdb_idx - 1], a_pdb)

    print("  Alignment-based lookup (UniProt position -> PDB residue, regardless of numbering):")
    for pos in positions:
        if pos in mapping:
            resseq, aa = mapping[pos]
            uni_aa = uniprot_seq[pos - 1]
            match = "match" if aa == uni_aa else "MISMATCH"
            print(f"    UniProt pos {pos} ({uni_aa}) -> PDB residue {resseq} ({aa})  [{match}]")
        else:
            print(f"    UniProt pos {pos}: not resolved in this crystal structure (missing density or not in construct)")


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    akt1_fasta = os.path.join(DATA_DIR, "AKT1.fasta")
    if not os.path.exists(akt1_fasta):
        print("data/AKT1.fasta not found -- run script 01 first.")
        return
    uniprot_seq = read_fasta(akt1_fasta)

    print(f"Downloading {PDB_ID} from RCSB...")
    pdb_path = download_pdb(PDB_ID)
    print(f"  Saved to {pdb_path}")

    chains = parse_atom_records(pdb_path)
    print(f"\nChains found: {list(chains.keys())}")
    for chain_id, residues in chains.items():
        print(f"  Chain {chain_id}: {len(residues)} resolved residues, "
              f"numbering {min(residues)}-{max(residues)}")

    # Use the largest protein chain as the AKT1 chain (co-crystal ligand
    # chains, if any residues sneak through, will be much smaller).
    main_chain_id = max(chains, key=lambda c: len(chains[c]))
    main_chain = chains[main_chain_id]
    print(f"\nUsing chain {main_chain_id} as the AKT1 chain for comparison.\n")

    direct_lookup(main_chain, POSITIONS_OF_INTEREST)
    print()
    alignment_lookup(main_chain, uniprot_seq, POSITIONS_OF_INTEREST)

    print(
        "\nHow to read this: if the direct lookup and alignment-based lookup "
        "agree, the structure's numbering matches UniProt numbering directly "
        "and position 59 (Gln in UniProt) is a genuine feature of AKT1, not "
        "a numbering artifact -- the 'Glu-59' in the docking paper likely "
        "refers to a different construct/numbering convention, or a "
        "different isoform/species, and should be checked against that "
        "paper's own methods before being used in a figure. If the two "
        "checks disagree, trust the alignment-based one."
    )


if __name__ == "__main__":
    main()