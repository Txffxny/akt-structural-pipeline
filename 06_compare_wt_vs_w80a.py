"""
06_compare_wt_vs_w80a.py

Compares each isoform's WT vs W80A ColabFold models by superimposing
them and computing RMSD at the MK-2206/miransertib pocket residues.

FIX -- separate output files, not one combined file: an earlier version
wrote WT and W80A into a single PDB file with a manually-inserted
ENDMDL line between them, but no matching MODEL header lines wrapping
each structure. That's a malformed multi-model PDB -- PyMOL split it
into two objects based on the stray ENDMDL, but the second object's
atoms weren't parsed correctly as a result (confirmed: coloring it
reported no atom count, unlike the first object's clean "3918 atoms").
This version writes WT and W80A as two separate, independently valid
PDB files instead -- load both with separate `load` commands in PyMOL/
ChimeraX rather than expecting one file to contain both.

FIX -- model-matched comparison, not rank-matched: ColabFold runs 5
separately-trained AlphaFold model weights and ranks by pLDDT.
"rank_001" doesn't guarantee the same underlying model was used for
both WT and W80A runs -- this version finds whichever model NUMBER is
shared between the two result sets and uses that pair.

CRITICAL DESIGN POINT -- two separate superpositions, not one: aligns
on the PH domain block (1-130) and kinase domain block (130-480)
independently, since PAE plots showed their relative orientation isn't
trustworthy -- mixing them would contaminate the residue-80 comparison
with that interdomain uncertainty.

AKT2 pocket residue numbering caveat: live pairwise alignment against
AKT1 maps pocket positions correctly (AKT2 is 481 aa vs AKT1's 480 aa).

Output (per isoform):
  - data/structural_comparison/<isoform>_PH_domain_WT.pdb
  - data/structural_comparison/<isoform>_PH_domain_W80A.pdb
  - data/structural_comparison/<isoform>_kinase_domain_WT.pdb
  - data/structural_comparison/<isoform>_kinase_domain_W80A.pdb
  - data/structural_comparison/rmsd_summary.csv
"""

import os
import re
import glob
import csv
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.Align import PairwiseAligner, substitution_matrices

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
COLABFOLD_DIR = os.path.join(DATA_DIR, "colabfold_results")
OUT_DIR = os.path.join(DATA_DIR, "structural_comparison")
os.makedirs(OUT_DIR, exist_ok=True)

PH_DOMAIN_REGION = (1, 130)
KINASE_DOMAIN_REGION = (130, 480)

POCKET_RESIDUES_PH = [53, 59, 78, 80]
POCKET_RESIDUES_KINASE = [201, 264, 270, 272, 273, 274, 292, 296]

ISOFORMS = ["AKT1", "AKT2"]

MODEL_NUM_RE = re.compile(r"_model_(\d+)_seed")


def find_all_models(isoform: str, variant: str) -> dict:
    pattern = os.path.join(COLABFOLD_DIR, f"{isoform}_{variant}_*", "**",
                            f"{isoform}_{variant}_*_unrelaxed_rank_*.pdb")
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        raise FileNotFoundError(
            f"No PDBs found for {isoform}_{variant} under {COLABFOLD_DIR}. "
            f"Check the folder was actually extracted."
        )
    models = {}
    for path in matches:
        m = MODEL_NUM_RE.search(os.path.basename(path))
        if m:
            models[int(m.group(1))] = path
    return models


def find_matched_pair(isoform: str):
    wt_models = find_all_models(isoform, "WT")
    w80a_models = find_all_models(isoform, "W80A")

    shared = sorted(set(wt_models) & set(w80a_models))
    if not shared:
        raise ValueError(
            f"{isoform}: no model number is shared between WT {sorted(wt_models)} "
            f"and W80A {sorted(w80a_models)} -- cannot do an apples-to-apples comparison."
        )
    chosen = shared[0]
    if len(shared) > 1:
        print(f"  Note: models {shared} shared between WT and W80A for {isoform}; "
              f"using model_{chosen} (lowest shared number).")
    return wt_models[chosen], w80a_models[chosen], chosen


def get_all_ca_atoms(structure):
    atoms = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    atoms[residue.id[1]] = residue["CA"]
        break
    return atoms


def read_fasta(path):
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    return "".join(lines[1:])


def map_akt1_pocket_to_akt2():
    akt1_seq = read_fasta(os.path.join(DATA_DIR, "AKT1.fasta"))
    akt2_seq = read_fasta(os.path.join(DATA_DIR, "AKT2.fasta"))

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"

    alignment = aligner.align(akt1_seq, akt2_seq)[0]
    aligned_akt1, aligned_akt2 = alignment[0], alignment[1]

    mapping = {}
    akt1_pos, akt2_pos = 0, 0
    for a1, a2 in zip(aligned_akt1, aligned_akt2):
        if a1 != "-":
            akt1_pos += 1
        if a2 != "-":
            akt2_pos += 1
        if a1 != "-" and a2 != "-":
            mapping[akt1_pos] = akt2_pos

    return mapping


def compute_rmsd(coords_a: dict, coords_b: dict, residue_numbers):
    ca, cb = [], []
    for resnum in residue_numbers:
        if resnum in coords_a and resnum in coords_b:
            ca.append(coords_a[resnum])
            cb.append(coords_b[resnum])
    if len(ca) < 1:
        return None, 0
    ca, cb = np.array(ca), np.array(cb)
    rmsd = np.sqrt(np.mean(np.sum((ca - cb) ** 2, axis=1)))
    return rmsd, len(ca)


def save_structure(structure, out_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path)


def superimpose_on_region(wt_path, w80a_path, region_start, region_end,
                           wt_out_path, w80a_out_path):
    parser = PDBParser(QUIET=True)
    wt_structure = parser.get_structure("WT", wt_path)
    w80a_structure = parser.get_structure("W80A", w80a_path)

    wt_all = get_all_ca_atoms(wt_structure)
    w80a_all = get_all_ca_atoms(w80a_structure)

    fit_resnums = sorted(
        r for r in (set(wt_all) & set(w80a_all)) if region_start <= r <= region_end
    )
    if len(fit_resnums) < 10:
        raise ValueError(
            f"Only {len(fit_resnums)} common residues in region {region_start}-"
            f"{region_end} -- too few to fit a superposition, don't trust this."
        )

    fixed = [wt_all[r] for r in fit_resnums]
    moving = [w80a_all[r] for r in fit_resnums]

    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    region_baseline_rmsd = superimposer.rms

    all_w80a_atoms = list(w80a_structure.get_atoms())
    superimposer.apply(all_w80a_atoms)

    w80a_all_transformed = get_all_ca_atoms(w80a_structure)

    wt_coords = {r: a.get_coord() for r, a in wt_all.items()}
    w80a_coords = {r: a.get_coord() for r, a in w80a_all_transformed.items()}

    save_structure(wt_structure, wt_out_path)
    save_structure(w80a_structure, w80a_out_path)

    return wt_coords, w80a_coords, region_baseline_rmsd, len(fit_resnums)


def main():
    summary_rows = []

    print("Mapping AKT1 pocket residue numbers onto AKT2 (live pairwise alignment)...")
    akt1_to_akt2_map = map_akt1_pocket_to_akt2()

    for isoform in ISOFORMS:
        print(f"\n{'='*60}\n{isoform}\n{'='*60}")

        wt_path, w80a_path, model_used = find_matched_pair(isoform)
        print(f"  Using model_{model_used} for BOTH WT and W80A (model-matched):")
        print(f"  WT:   {wt_path}")
        print(f"  W80A: {w80a_path}")

        if isoform == "AKT1":
            ph_pocket = POCKET_RESIDUES_PH
            kinase_pocket = POCKET_RESIDUES_KINASE
        else:
            ph_pocket = [akt1_to_akt2_map.get(p) for p in POCKET_RESIDUES_PH]
            ph_pocket = [p for p in ph_pocket if p is not None]
            kinase_pocket = [akt1_to_akt2_map.get(p) for p in POCKET_RESIDUES_KINASE]
            kinase_pocket = [p for p in kinase_pocket if p is not None]
            print(f"  AKT2-mapped PH pocket residues: {ph_pocket}")
            print(f"  AKT2-mapped kinase pocket residues: {kinase_pocket}")

        ph_wt_out = os.path.join(OUT_DIR, f"{isoform}_PH_domain_WT.pdb")
        ph_w80a_out = os.path.join(OUT_DIR, f"{isoform}_PH_domain_W80A.pdb")
        wt_c, w80a_c, ph_baseline_rmsd, n_ph_fit = superimpose_on_region(
            wt_path, w80a_path, *PH_DOMAIN_REGION, ph_wt_out, ph_w80a_out
        )
        print(f"\n  [PH domain alignment] baseline RMSD over {n_ph_fit} fitted residues: "
              f"{ph_baseline_rmsd:.3f} A")
        ph_pocket_rmsd, n_ph_pocket = compute_rmsd(wt_c, w80a_c, ph_pocket)
        if ph_pocket_rmsd is not None:
            print(f"  [PH domain alignment] pocket-residue RMSD "
                  f"(residues {ph_pocket}, includes mutation site): "
                  f"{ph_pocket_rmsd:.3f} A over {n_ph_pocket} residues")
        print(f"  Saved: {ph_wt_out}")
        print(f"  Saved: {ph_w80a_out}")

        kinase_wt_out = os.path.join(OUT_DIR, f"{isoform}_kinase_domain_WT.pdb")
        kinase_w80a_out = os.path.join(OUT_DIR, f"{isoform}_kinase_domain_W80A.pdb")
        wt_c2, w80a_c2, kinase_baseline_rmsd, n_kinase_fit = superimpose_on_region(
            wt_path, w80a_path, *KINASE_DOMAIN_REGION, kinase_wt_out, kinase_w80a_out
        )
        print(f"\n  [Kinase domain alignment] baseline RMSD over {n_kinase_fit} fitted residues: "
              f"{kinase_baseline_rmsd:.3f} A")
        kinase_pocket_rmsd, n_kinase_pocket = compute_rmsd(wt_c2, w80a_c2, kinase_pocket)
        if kinase_pocket_rmsd is not None:
            print(f"  [Kinase domain alignment] pocket-residue RMSD "
                  f"(residues {kinase_pocket}): "
                  f"{kinase_pocket_rmsd:.3f} A over {n_kinase_pocket} residues")
        print(f"  Saved: {kinase_wt_out}")
        print(f"  Saved: {kinase_w80a_out}")

        summary_rows.append({
            "isoform": isoform,
            "model_used": model_used,
            "ph_domain_baseline_rmsd_A": round(ph_baseline_rmsd, 3),
            "ph_pocket_rmsd_A": round(ph_pocket_rmsd, 3) if ph_pocket_rmsd is not None else "NA",
            "n_ph_pocket_residues": n_ph_pocket,
            "kinase_domain_baseline_rmsd_A": round(kinase_baseline_rmsd, 3),
            "kinase_pocket_rmsd_A": round(kinase_pocket_rmsd, 3) if kinase_pocket_rmsd is not None else "NA",
            "n_kinase_pocket_residues": n_kinase_pocket,
        })

    summary_path = os.path.join(OUT_DIR, "rmsd_summary.csv")
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\n\nSaved summary: {summary_path}")
    print(
        "\nTo view in PyMOL, load the WT and W80A files for a given domain "
        "as two SEPARATE load commands, e.g.:\n"
        "  load .../AKT1_PH_domain_WT.pdb\n"
        "  load .../AKT1_PH_domain_W80A.pdb\n"
        "then color/select as before -- no split_states needed, each file "
        "is already a single clean object."
    )


if __name__ == "__main__":
    main()