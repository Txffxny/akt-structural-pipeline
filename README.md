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
