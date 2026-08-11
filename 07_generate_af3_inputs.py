"""
07_generate_af3_inputs.py

Generates AlphaFold3 input JSON files for every combination of AKT1
variant (WT, W80A) x candidate ligand (5 total), using the real AKT1
sequences already fetched by the earlier pipeline (01_fetch_akt_sequences.py,
02_generate_w80a_mutants.py) rather than retyping the sequence here.

Ligand set and ground-truth status:
  Inhibitor VIII  - allosteric - PDB 3O96  (verified ground truth)
  Miransertib     - allosteric - PDB 5KCV  (verified ground truth)
  MK-2206         - allosteric - no crystal structure (exploratory only)
  Ipatasertib     - orthosteric - PDB 4EKL (verified ground truth)
  Capivasertib    - orthosteric - PDB 4GV1 (verified ground truth)

All SMILES below are isomeric/stereo-defined, sourced from RCSB ligand
pages (derived from real crystal geometry) except MK-2206, which has no
crystal structure and has no stereocenters to begin with (confirmed by
inspection -- the cyclobutylamine substituent sits at a symmetric ring
position).

Requires: data/AKT1_WT.fasta and data/AKT1_W80A.fasta (from script 02
of the akt-structural-pipeline project) to be present in the same data/
directory as this script.

Output: data/af3_inputs/<variant>_<ligand>.json, one per combination,
plus a manifest.csv documenting which pairs have real ground truth to
validate against and which are exploratory.
"""

import os
import json
import csv

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUT_DIR = os.path.join(DATA_DIR, "af3_inputs")
os.makedirs(OUT_DIR, exist_ok=True)

LIGANDS = {
    "InhibitorVIII": {
        "smiles": "O=C2Nc1ccccc1N2C8CCN(Cc7ccc(c4nc6c(nc4c3ccccc3)cc5ncnc5c6)cc7)CC8",
        "class": "allosteric",
        "ground_truth_pdb": "3O96",
    },
    "Miransertib": {
        "smiles": "Nc1ncccc1c2nc3ccc(nc3n2c4ccc(cc4)C5(N)CCC5)c6ccccc6",
        "class": "allosteric",
        "ground_truth_pdb": "5KCV",
    },
    "MK2206": {
        "smiles": "C1CC(C1)(C2=CC=C(C=C2)C3=C(C=C4C(=N3)C=CN5C4=NNC5=O)C6=CC=CC=C6)N",
        "class": "allosteric",
        "ground_truth_pdb": "",  # no crystal structure -- exploratory
    },
    "Ipatasertib": {
        "smiles": "CC(C)NC[C@@H](C(=O)N1CCN(CC1)c2ncnc3[C@H](O)C[C@@H](C)c23)c4ccc(Cl)cc4",
        "class": "orthosteric",
        "ground_truth_pdb": "4EKL",
    },
    "Capivasertib": {
        "smiles": "NC1(CCN(CC1)c2ncnc3[nH]ccc23)C(=O)N[C@@H](CCO)c4ccc(Cl)cc4",
        "class": "orthosteric",
        "ground_truth_pdb": "4GV1",
    },
}

VARIANTS = {
    "AKT1_WT": "AKT1_WT.fasta",
    "AKT1_W80A": "AKT1_W80A.fasta",
}


def read_fasta_sequence(path):
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    return "".join(lines[1:])


def build_af3_input(job_name, protein_sequence, ligand_smiles):
    return {
        "name": job_name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": protein_sequence,
                }
            },
            {
                "ligand": {
                    "id": ["L"],
                    "smiles": ligand_smiles,
                }
            },
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }


def main():
    manifest_rows = []

    for variant_name, fasta_filename in VARIANTS.items():
        fasta_path = os.path.join(DATA_DIR, fasta_filename)
        if not os.path.exists(fasta_path):
            print(f"MISSING: {fasta_path} -- run the AKT sequence/mutant "
                  f"generation scripts first (01_fetch_akt_sequences.py, "
                  f"02_generate_w80a_mutants.py) before this one.")
            continue

        sequence = read_fasta_sequence(fasta_path)
        print(f"{variant_name}: loaded {len(sequence)} aa sequence from {fasta_path}")

        for ligand_name, ligand_info in LIGANDS.items():
            job_name = f"{variant_name}_{ligand_name}"
            input_json = build_af3_input(job_name, sequence, ligand_info["smiles"])

            out_path = os.path.join(OUT_DIR, f"{job_name}.json")
            with open(out_path, "w") as fh:
                json.dump(input_json, fh, indent=2)

            gt = ligand_info["ground_truth_pdb"] or "NONE (exploratory)"
            print(f"  {ligand_name} ({ligand_info['class']}, ground truth: {gt}) -> {out_path}")

            manifest_rows.append({
                "job_name": job_name,
                "variant": variant_name,
                "ligand": ligand_name,
                "ligand_class": ligand_info["class"],
                "ground_truth_pdb": ligand_info["ground_truth_pdb"],
                "file": out_path,
            })

    manifest_path = os.path.join(OUT_DIR, "manifest.csv")
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nManifest saved: {manifest_path}")
    print(f"Total input files generated: {len(manifest_rows)}")

    verified = [r for r in manifest_rows if r["ground_truth_pdb"]]
    exploratory = [r for r in manifest_rows if not r["ground_truth_pdb"]]
    print(f"\nRECOMMENDED RUN ORDER:")
    print(f"  1. Verified pairs first ({len(verified)} jobs) -- these have real crystal")
    print(f"     structures to check the prediction against:")
    for r in verified:
        print(f"       {r['job_name']}  (compare against PDB {r['ground_truth_pdb']})")
    print(f"  2. Exploratory pairs after ({len(exploratory)} jobs) -- no ground truth,")
    print(f"     treat these results as hypothesis-generating only:")
    for r in exploratory:
        print(f"       {r['job_name']}")


if __name__ == "__main__":
    main()