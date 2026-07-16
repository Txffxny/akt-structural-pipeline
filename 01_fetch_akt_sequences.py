"""
01_fetch_akt_sequences.py

Pulls the canonical human AKT1, AKT2, AKT3 sequences from UniProt and
saves them as individual FASTA files. Run this first -- everything
downstream reads from data/*.fasta.

Requires internet access to www.uniprot.org (works fine on a normal
machine/HPC login node; will not work inside network-locked sandboxes).
"""

import os
import requests

# Canonical UniProt accessions, human
ACCESSIONS = {
    "AKT1": "P31749",
    "AKT2": "P31751",
    "AKT3": "Q9Y243",
}

OUT_DIR = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(OUT_DIR, exist_ok=True)


def fetch_fasta(accession: str) -> str:
    """Fetch canonical isoform FASTA text from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.text


def parse_fasta(fasta_text: str):
    """Return (header, sequence) from a single-record FASTA string."""
    lines = fasta_text.strip().splitlines()
    header = lines[0].lstrip(">")
    seq = "".join(lines[1:])
    return header, seq


def main():
    summary = []
    for gene, acc in ACCESSIONS.items():
        print(f"Fetching {gene} ({acc}) from UniProt...")
        fasta_text = fetch_fasta(acc)
        header, seq = parse_fasta(fasta_text)

        out_path = os.path.join(OUT_DIR, f"{gene}.fasta")
        with open(out_path, "w") as fh:
            fh.write(f">{gene}|{acc}\n{seq}\n")

        summary.append((gene, acc, len(seq), seq[79] if len(seq) >= 80 else "N/A"))
        print(f"  -> {len(seq)} aa, saved to {out_path}")

    print("\nSanity check -- residue at position 80 (1-indexed) for each isoform:")
    print(f"{'Gene':<6}{'Accession':<12}{'Length':<8}{'Res80'}")
    for gene, acc, length, res80 in summary:
        flag = "" if res80 == "W" else "  <-- expected Trp (W), check numbering/isoform!"
        print(f"{gene:<6}{acc:<12}{length:<8}{res80}{flag}")


if __name__ == "__main__":
    main()
