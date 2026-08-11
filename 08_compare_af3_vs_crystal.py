"""
08_compare_af3_vs_crystal.py

Tests whether AF3's confident predictions (high ipTM, ~0.96-0.97 across
all four verified-ligand WT runs) are actually CORRECT, not just
self-consistent. Confidence scores tell you AF3 committed to a specific
pose; they don't tell you that pose matches reality. This compares the
predicted structure directly against the real crystal structure -- for
BOTH WT and W80A, against the same real crystal reference, so any
deviation in the W80A comparison reflects the mutation's effect relative
to true biology, not just relative to AF3's own WT prediction.

APPROACH: align the predicted structure onto the real crystal using ONLY
the kinase domain CA atoms (residues 130-480) -- this is the one region
reliably well-resolved in both the prediction and every crystal structure
used here, so it's a trustworthy reference frame. After that alignment:

  1. PH-domain RMSD (predicted vs crystal, residues 1-130): if AF3
     correctly predicted the ligand-induced closed conformation, the PH
     domain should land close to where it really is. A large RMSD here
     means AF3's confident prediction is confidently WRONG about domain
     packing, regardless of what ipTM says.
  2. Ligand centroid distance: does the predicted ligand end up in
     roughly the same 3D location as the real, crystallographically-
     resolved ligand, once the kinase domain frame is aligned?

Validated on the WT runs: PH-domain RMSD ~1.9-2.4 A and ligand centroid
distance ~0.64-0.85 A against 3O96 and 5KCV respectively -- AF3's high
confidence for these two ligands is backed by genuinely accurate
structures, not just self-consistent guessing. That validation is what
makes the W80A comparison below meaningful.

CAVEAT on the ligand comparison: this uses ligand CENTROID distance, not
atom-level RMSD -- AF3's predicted ligand atom ordering isn't guaranteed
to match the crystal ligand's PDB atom ordering, and atom-level
correspondence would need substructure matching (RDKit) to do properly.
Centroid distance is a coarser but honest metric: it answers "is the
ligand in the right neighborhood", not "is the exact pose/orientation
correct". Treat it as a first-pass check, not a final verdict.

Requires:
  - Real crystal PDB files: data/3O96.pdb (already fetched by script 04
    of the akt-structural-pipeline project) and data/5KCV.pdb (fetched
    automatically by this script if not already present).
  - AF3 predicted results, downloaded and unzipped into
    data/af3_results/<job_name>/ for each of AKT1_WT_<ligand> and
    AKT1_W80A_<ligand>. Missing folders are skipped cleanly, not errored.

Output: data/af3_vs_crystal/comparison_summary.csv
"""

import os
import glob
import csv
import numpy as np
import requests
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
AF3_RESULTS_DIR = os.path.join(DATA_DIR, "af3_results")
OUT_DIR = os.path.join(DATA_DIR, "af3_vs_crystal")
os.makedirs(OUT_DIR, exist_ok=True)

PH_DOMAIN_REGION = (1, 130)
KINASE_DOMAIN_REGION = (130, 480)

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

CRYSTAL_TARGETS = {
    "InhibitorVIII": {"pdb_id": "3O96", "ligand_resname": "IQO"},
    "Miransertib": {"pdb_id": "5KCV", "ligand_resname": "6S1"},
}

VARIANTS = ["AKT1_WT", "AKT1_W80A"]


def fetch_crystal_structure(pdb_id: str) -> str:
    out_path = os.path.join(DATA_DIR, f"{pdb_id}.pdb")
    if os.path.exists(out_path):
        return out_path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    with open(out_path, "w") as fh:
        fh.write(resp.text)
    print(f"  Downloaded {pdb_id} -> {out_path}")
    return out_path


def get_protein_ca_atoms(structure):
    atoms = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in STANDARD_AA and "CA" in residue:
                    atoms[residue.id[1]] = residue["CA"]
        if atoms:
            break
    return atoms


def get_ligand_atoms(structure, ligand_resname=None):
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in STANDARD_AA or resname == "HOH":
                    continue
                if ligand_resname and resname != ligand_resname:
                    continue
                for atom in residue:
                    if atom.element != "H":
                        coords.append(atom.get_coord())
        if coords:
            break
    return np.array(coords)


def find_predicted_cif(job_name: str) -> str:
    job_dir = os.path.join(AF3_RESULTS_DIR, job_name)
    if not os.path.isdir(job_dir):
        raise FileNotFoundError(
            f"No results folder found at {job_dir} -- download and unzip "
            f"the AF3 result for {job_name} there first."
        )
    cif_files = glob.glob(os.path.join(job_dir, "**", "*.cif"), recursive=True)
    if not cif_files:
        raise FileNotFoundError(f"No .cif files found under {job_dir}")

    top_level = [f for f in cif_files if "seed-" not in os.path.basename(f)]
    if top_level:
        chosen = sorted(top_level, key=len)[0]
    else:
        preferred = [f for f in cif_files if "sample-0" in f]
        chosen = preferred[0] if preferred else sorted(cif_files)[0]
        print(f"  WARNING: no top-level summary .cif found for {job_name}, "
              f"falling back to sample-0 -- verify this is the intended model.")

    if len(cif_files) > 1:
        print(f"  {len(cif_files)} .cif files found for {job_name}, using: {chosen}")
    return chosen


def compute_rmsd(coords_a, coords_b):
    if len(coords_a) == 0 or len(coords_b) == 0:
        return None
    return float(np.sqrt(np.mean(np.sum((coords_a - coords_b) ** 2, axis=1))))


def compare_one(job_name, ligand_name, crystal_info):
    print(f"\n{'='*60}\n{job_name} vs {crystal_info['pdb_id']}\n{'='*60}")

    crystal_path = fetch_crystal_structure(crystal_info["pdb_id"])
    predicted_path = find_predicted_cif(job_name)
    print(f"  Predicted: {predicted_path}")
    print(f"  Crystal:   {crystal_path}")

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)

    crystal_structure = pdb_parser.get_structure("crystal", crystal_path)
    predicted_structure = cif_parser.get_structure("predicted", predicted_path)

    crystal_ca = get_protein_ca_atoms(crystal_structure)
    predicted_ca = get_protein_ca_atoms(predicted_structure)

    fit_resnums = sorted(
        r for r in (set(crystal_ca) & set(predicted_ca))
        if KINASE_DOMAIN_REGION[0] <= r <= KINASE_DOMAIN_REGION[1]
    )
    if len(fit_resnums) < 10:
        raise ValueError(
            f"Only {len(fit_resnums)} common kinase-domain residues between "
            f"prediction and crystal -- numbering mismatch likely, don't "
            f"trust this comparison. Check residue numbering in both files."
        )

    fixed = [crystal_ca[r] for r in fit_resnums]
    moving = [predicted_ca[r] for r in fit_resnums]

    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    kinase_fit_rmsd = superimposer.rms

    all_predicted_atoms = list(predicted_structure.get_atoms())
    superimposer.apply(all_predicted_atoms)

    predicted_ca_transformed = get_protein_ca_atoms(predicted_structure)
    ph_resnums = sorted(
        r for r in (set(crystal_ca) & set(predicted_ca_transformed))
        if PH_DOMAIN_REGION[0] <= r <= PH_DOMAIN_REGION[1]
    )
    if len(ph_resnums) >= 3:
        crystal_ph_coords = np.array([crystal_ca[r].get_coord() for r in ph_resnums])
        predicted_ph_coords = np.array([predicted_ca_transformed[r].get_coord() for r in ph_resnums])
        ph_domain_rmsd = compute_rmsd(crystal_ph_coords, predicted_ph_coords)
        n_ph = len(ph_resnums)
    else:
        ph_domain_rmsd, n_ph = None, len(ph_resnums)

    crystal_ligand_coords = get_ligand_atoms(crystal_structure, crystal_info["ligand_resname"])
    predicted_ligand_coords = get_ligand_atoms(predicted_structure, ligand_resname=None)

    if len(crystal_ligand_coords) > 0 and len(predicted_ligand_coords) > 0:
        crystal_centroid = crystal_ligand_coords.mean(axis=0)
        predicted_centroid = predicted_ligand_coords.mean(axis=0)
        centroid_distance = float(np.linalg.norm(crystal_centroid - predicted_centroid))
        n_crystal_ligand_atoms = len(crystal_ligand_coords)
        n_predicted_ligand_atoms = len(predicted_ligand_coords)
    else:
        centroid_distance = None
        n_crystal_ligand_atoms = len(crystal_ligand_coords)
        n_predicted_ligand_atoms = len(predicted_ligand_coords)
        print(f"  WARNING: could not find ligand atoms in one or both structures "
              f"(crystal: {n_crystal_ligand_atoms}, predicted: {n_predicted_ligand_atoms}) "
              f"-- check ligand_resname / chain identification.")

    print(f"  Kinase-domain fit RMSD (reference frame quality): {kinase_fit_rmsd:.3f} A "
          f"over {len(fit_resnums)} residues")
    if ph_domain_rmsd is not None:
        print(f"  PH-domain RMSD (predicted vs real crystal): {ph_domain_rmsd:.3f} A "
              f"over {n_ph} residues")
    if centroid_distance is not None:
        print(f"  Ligand centroid distance (predicted vs real): {centroid_distance:.3f} A "
              f"({n_predicted_ligand_atoms} predicted heavy atoms vs "
              f"{n_crystal_ligand_atoms} crystal heavy atoms)")

    return {
        "job_name": job_name,
        "crystal_pdb": crystal_info["pdb_id"],
        "kinase_fit_rmsd_A": round(kinase_fit_rmsd, 3),
        "n_kinase_fit_residues": len(fit_resnums),
        "ph_domain_rmsd_A": round(ph_domain_rmsd, 3) if ph_domain_rmsd is not None else "NA",
        "n_ph_domain_residues": n_ph,
        "ligand_centroid_distance_A": round(centroid_distance, 3) if centroid_distance is not None else "NA",
        "n_predicted_ligand_atoms": n_predicted_ligand_atoms,
        "n_crystal_ligand_atoms": n_crystal_ligand_atoms,
    }


def main():
    results = []
    for variant in VARIANTS:
        for ligand_name, crystal_info in CRYSTAL_TARGETS.items():
            job_name = f"{variant}_{ligand_name}"
            try:
                result = compare_one(job_name, ligand_name, crystal_info)
                result["variant"] = variant
                results.append(result)
            except (FileNotFoundError, ValueError) as e:
                print(f"\nSkipping {job_name}: {e}")

    if not results:
        print("\nNo comparisons completed -- check that AF3 results are downloaded "
              "and unzipped into data/af3_results/<job_name>/")
        return

    out_path = os.path.join(OUT_DIR, "comparison_summary.csv")
    fieldnames = ["variant"] + [k for k in results[0].keys() if k != "variant"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\n\nSaved: {out_path}")
    print(
        "\nHow to read this:\n"
        "  - Kinase-fit RMSD is the reference-frame quality -- expect this small "
        "(well under 2 A); it's not the interesting number.\n"
        "  - PH-domain RMSD is the real test: if AF3's confident prediction "
        "(ipTM ~0.96-0.97) actually got the ligand-induced closed conformation "
        "right, this should be small too. A large value here means high "
        "confidence does NOT mean correct -- exactly the gap ipTM/pLDDT can't "
        "reveal on their own.\n"
        "  - Ligand centroid distance is a coarse 'is it in the right "
        "neighborhood' check, not a precise pose comparison -- a small value "
        "is reassuring, a large value is a clear red flag, but don't over-"
        "interpret small differences given the centroid-only approximation.\n"
        "  - Compare WT vs W80A rows for the same ligand/crystal directly: "
        "since both are checked against the SAME real crystal structure, any "
        "increase in RMSD/centroid distance for W80A reflects the mutation's "
        "effect relative to true biology, not just relative to AF3's own WT "
        "prediction."
    )


if __name__ == "__main__":
    main()