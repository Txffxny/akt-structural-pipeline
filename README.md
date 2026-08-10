# AKT Structural Side Project — Pipeline

Three scripts, run in order, on a machine with normal internet access
(needs to reach `rest.uniprot.org`).

```
pip install biopython requests
python 01_fetch_akt_sequences.py
python 02_generate_w80a_mutants.py
python 03_pocket_conservation_analysis.py
```

## What each step does

**01_fetch_akt_sequences.py**
Pulls canonical AKT1 (P31749), AKT2 (P31751), AKT3 (Q9Y243) sequences
from UniProt, saves as `data/AKT1.fasta` etc. Prints a sanity check that
residue 80 is Trp (W) in each — if that fails, stop and check numbering
before going further (canonical vs. isoform-2 numbering differs by a few
residues in some databases).

**02_generate_w80a_mutants.py**
Reads AKT1/AKT2 FASTA, confirms residue 80 == W, writes W80A mutants.
Outputs `data/AKT1_WT.fasta`, `data/AKT1_W80A.fasta`, same for AKT2, plus
`data/all_variants.fasta` as a single batch file. Submit these directly to
ColabFold (https://github.com/sokrypton/ColabFold — the AlphaFold2.ipynb
notebook on Google Colab is the easiest route, free GPU, no local install).
Each FASTA header becomes the job name so WT vs W80A predictions come back
clearly labeled.

**03_pocket_conservation_analysis.py**
Aligns AKT2 and AKT3 against AKT1 and checks whether the residues known to
form the MK-2206/miransertib allosteric pocket (Trp-80 and its structural
neighbors: Asn-53, Glu-59, Leu-78, Leu-264, Val-270, Tyr-272, Arg-273,
Asp-274, Asp-292, Cys-296, Val-201 — from the MK-2206 and Inhibitor VIII
co-crystal structures) are conserved or diverge across isoforms. Outputs
`data/pocket_conservation.csv` and a printed table. Any residue flagged
"NO" is a candidate explanation for AKT3's lower potency to these
allosteric drugs, and worth highlighting on the eventual structure figure.

## After this: the actual fold

1. Take `data/all_variants.fasta` (or the individual WT/W80A files) into
   ColabFold on Google Colab.
2. Download the resulting PDB files.
3. Load them alongside the existing experimental structures for
   comparison — PDB 3O96 (AKT1 + Inhibitor VIII) and the ΔHM-AKT1 +
   Inhibitor VIII co-crystal described in Wu et al., PLOS ONE 2010 — in
   PyMOL or ChimeraX, overlay WT vs W80A, and colour the pocket residues
   from step 03's output.

## Notes

- These scripts were tested end-to-end on synthetic placeholder sequences
  to confirm the logic (mutation, alignment, conservation flagging) works
  correctly. They have **not** been run against the real UniProt data —
  that fetch needs to happen from your own machine.
- If you later add AKT3 to the mutagenesis plan, just add `"AKT3"` to
  `ISOFORMS_TO_MUTATE` in script 02.

  ## WT vs W80A structural comparison (script 06)

Compared ColabFold-predicted WT and W80A models for both AKT1 and AKT2,
using two separate domain-restricted superpositions (PH domain, kinase
domain) rather than one whole-chain alignment -- necessary because PAE
plots showed AlphaFold is confident about each domain's internal fold
but not their relative orientation, so a whole-chain alignment would
let that interdomain uncertainty contaminate the pocket-residue RMSD.

Comparisons are model-matched (same underlying AlphaFold model number
used for WT and W80A), not rank-matched -- an initial rank_001-vs-
rank_001 comparison for AKT2 silently compared two different underlying
models (model_4 vs model_3), inflating the apparent baseline RMSD
roughly 6-fold (6.933 A vs the corrected 1.211 A). Always verify the
model number, not just the rank, before trusting a ColabFold structural
comparison.

**Result:** in all four comparisons (AKT1 PH domain, AKT1 kinase domain,
AKT2 PH domain, AKT2 kinase domain), the pocket-residue RMSD was equal
to or lower than that domain's own baseline RMSD -- i.e. AlphaFold's
static single-state prediction does not show a detectable localized
structural perturbation at the pocket residues (including residue 80
itself) as a result of the W80A mutation, in either isoform.

This is a legitimate negative result, not a failed experiment: it's
consistent with Trp80's proposed role being about stabilizing a specific
*conformation* (the closed, inhibitor-competent PH-in state) via a
stacking interaction with the inhibitor, rather than being a rigid
structural scaffold residue whose loss should visibly deform the fold.
A static AlphaFold prediction cannot resolve a conformational-
equilibrium-shifting mutation -- only real drug-response data (the
4-drug panel) can actually test that mechanism.
