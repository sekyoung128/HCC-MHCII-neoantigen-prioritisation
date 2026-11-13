#!/usr/bin/env python3
"""Phase 1: derive and validate mutation-centred WT/MUT peptide pairs.

Run from inside the mhcii-neoantigen-discovery project after downloading
canonical UniProt FASTA files to raw/.

This creates an initial *one-window-per-length* design, NOT an exhaustive
scan of all possible mutation-spanning MHC-II peptide windows.
"""

import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
RESULTS = ROOT / "results"
RAW.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

MUTATIONS = [
    {
        "gene": "CTNNB1", "variant": "S45P", "position": 45,
        "wt_aa": "S", "mut_aa": "P", "cdna": "c.133T>C",
        "uniprot": "P35222", "refseq_transcript": "NM_001904.4",
        "refseq_protein": "NP_001895.1", "protein_length": 781,
        "expected_21mer_wt": "IHSGATTTAPSLSGKGNPEEE",
        "source_url": "https://ckb.genomenon.com/geneVariant/show?geneVariantId=786",
    },
    {
        "gene": "TP53", "variant": "R249S", "position": 249,
        "wt_aa": "R", "mut_aa": "S", "cdna": "c.747G>T",
        "uniprot": "P04637", "refseq_transcript": "NM_000546.6",
        "refseq_protein": "NP_000537.3", "protein_length": 393,
        "expected_21mer_wt": "NSSCMGGMNRRPILTIITLED",
        "source_url": "https://www.ncbi.nlm.nih.gov/clinvar/RCV000013145/",
    },
]

LENGTHS = (15, 18, 21)
VALID_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path, accession):
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}. Download canonical UniProt FASTA for {accession} first."
        )
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
             if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Not a FASTA file: {path}")
    if accession not in lines[0]:
        raise ValueError(f"Unexpected FASTA accession in {path}: {lines[0]}")
    if any(line.startswith(">") for line in lines[1:]):
        raise ValueError(f"Expected exactly one reference sequence in {path}")
    sequence = "".join(lines[1:]).upper()
    if not sequence or set(sequence) - VALID_AAS:
        raise ValueError(f"Noncanonical or missing amino acids in {path}")
    return sequence


def save_csv(path, rows, columns):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


peptide_rows = []
validation_rows = []

for m in MUTATIONS:
    fasta_path = RAW / f"{m['gene']}_{m['uniprot']}.fasta"
    seq = read_fasta(fasta_path, m["uniprot"])
    pos = m["position"]

    if len(seq) != m["protein_length"]:
        raise ValueError(f"{m['gene']}: wrong protein length {len(seq)}")
    if seq[pos - 1] != m["wt_aa"]:
        raise ValueError(
            f"{m['gene']}: residue {pos} is {seq[pos - 1]}, "
            f"not expected {m['wt_aa']}"
        )
    check_start = pos - 10  # 21-mer: mutation is at position 11
    if seq[check_start - 1: check_start - 1 + 21] != m["expected_21mer_wt"]:
        raise ValueError(f"{m['gene']}: 21-mer context differs from reference check")

    validation_rows.append({
        "gene": m["gene"],
        "uniprot": m["uniprot"],
        "reference_length": len(seq),
        "confirmed_wt_residue": seq[pos - 1],
        "fasta_sha256": hashlib.sha256(fasta_path.read_bytes()).hexdigest(),
        "status": "PASS",
    })

    for length in LENGTHS:
        # 1-based protein coordinates. For 18mers mutation is residue 9 of 18.
        start = pos - (length - 1) // 2
        end = start + length - 1
        if start < 1 or end > len(seq):
            raise ValueError(f"{m['gene']}: peptide outside reference sequence")
        wt = seq[start - 1:end]
        mutation_index = pos - start  # zero-based peptide index
        mut = wt[:mutation_index] + m["mut_aa"] + wt[mutation_index + 1:]
        if len(wt) != length or len(mut) != length:
            raise ValueError("Incorrect peptide length")
        if sum(a != b for a, b in zip(wt, mut)) != 1:
            raise ValueError("WT and MUT must differ at exactly one residue")
        if wt[mutation_index] != m["wt_aa"] or mut[mutation_index] != m["mut_aa"]:
            raise ValueError("Mutation not at expected peptide position")

        for label, peptide in (("WT", wt), ("MUT", mut)):
            peptide_rows.append({
                "peptide_id": f"{m['gene']}_{m['variant']}_{length}_{label}",
                "gene": m["gene"],
                "variant": m["variant"],
                "uniprot": m["uniprot"],
                "reference_protein": m["refseq_protein"],
                "type": label,
                "length": length,
                "protein_start_1based": start,
                "protein_end_1based": end,
                "mutation_position_1based": pos,
                "mutation_position_in_peptide_1based": mutation_index + 1,
                "peptide": peptide,
            })

assert len(peptide_rows) == 12, "Expected 12 WT/MUT peptide rows"
assert len({row["peptide"] for row in peptide_rows}) == 12, \
    "Expected 12 distinct peptide sequences"

save_csv(
    RAW / "phase1_mutations.csv", MUTATIONS,
    list(MUTATIONS[0].keys()),
)
save_csv(
    RESULTS / "phase1_peptides.csv", peptide_rows,
    list(peptide_rows[0].keys()),
)
save_csv(
    RESULTS / "phase1_reference_validation.csv", validation_rows,
    list(validation_rows[0].keys()),
)
(RAW / "phase1_netmhciipan_input.txt").write_text(
    "\n".join(row["peptide"] for row in peptide_rows) + "\n",
    encoding="utf-8",
)

print("PHASE 1 REFERENCE VALIDATION PASSED")
for row in validation_rows:
    print(f"  {row['gene']}: {row['reference_length']} aa; "
          f"WT residue {row['confirmed_wt_residue']}; PASS")
print("Created raw/phase1_mutations.csv (2 variants)")
print("Created results/phase1_peptides.csv (12 sequences; 6 WT/MUT pairs)")
print("Created results/phase1_reference_validation.csv")
print("Created raw/phase1_netmhciipan_input.txt (12 lines)")
print("Next: 12 peptide sequences x 2 HLA-DR alleles = 24 prediction rows")
