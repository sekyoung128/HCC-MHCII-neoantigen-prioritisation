#!/usr/bin/env python3
"""Phase 1 exploratory extension: every mutation-spanning 15/18/21-aa window.

Run from an existing mhcii-neoantigen-discovery project containing:
  raw/phase1_mutations.csv
  raw/CTNNB1_P35222.fasta
  raw/TP53_P04637.fasta
  results/phase1_peptides.csv

Does NOT call NetMHCIIpan, modify prior Phase 0/1 files, or claim immunogenicity.
Creates independently auditable peptide manifests and web-server input files.
"""
from __future__ import annotations

import csv
import hashlib
from collections import Counter
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
RESULTS = ROOT / "results"
RAW_EXT = RAW / "phase1_extension"
RESULTS_EXT = RESULTS / "phase1_extension"
LENGTHS = (15, 18, 21)
GENES = {"CTNNB1": "P35222", "TP53": "P04637"}
HLA = ("DRB1_0101", "DRB1_0401")
ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing required file: {path}")
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def save_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot save empty data: {path}")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path, accession: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"Missing reference: {path}")
    lines = [s.strip() for s in path.read_text(encoding="utf-8").splitlines()
             if s.strip()]
    if len(lines) < 2 or not lines[0].startswith(">") or accession not in lines[0]:
        raise ValueError(f"Wrong FASTA header/accession: {path}")
    if any(s.startswith(">") for s in lines[1:]):
        raise ValueError(f"Multiple sequences in FASTA: {path}")
    seq = "".join(lines[1:]).upper()
    if not seq or (set(seq) - ALPHABET):
        raise ValueError(f"Invalid protein sequence: {path}")
    return seq


def main() -> None:
    mutations = read_csv(RAW / "phase1_mutations.csv")
    pilot = read_csv(RESULTS / "phase1_peptides.csv")
    if len(mutations) != 2 or {m["gene"] for m in mutations} != set(GENES):
        raise ValueError("Expected exactly the original CTNNB1 and TP53 variants")
    if len(pilot) != 12:
        raise ValueError("Original Phase 1 design must contain 12 peptide rows")

    rows: list[dict[str, object]] = []
    pairs: list[dict[str, object]] = []
    fasta_hashes: list[str] = []
    for m in sorted(mutations, key=lambda r: (r["gene"] != "CTNNB1", r["gene"])):
        gene = m["gene"]
        acc = GENES[gene]
        if m["uniprot"] != acc:
            raise ValueError(f"{gene}: UniProt mismatch")
        fasta_path = RAW / f"{gene}_{acc}.fasta"
        seq = read_fasta(fasta_path, acc)
        protein_length = int(m["protein_length"])
        pos = int(m["position"])  # one-based protein coordinate
        old = m["wt_aa"]
        new = m["mut_aa"]
        if len(seq) != protein_length or seq[pos - 1] != old or old == new:
            raise ValueError(f"{gene}: reference length/residue mismatch")
        fasta_sha256 = hashlib.sha256(fasta_path.read_bytes()).hexdigest()
        fasta_hashes.append(f"{gene} {acc}: {fasta_sha256}")

        for length in LENGTHS:
            # All possible windows of exactly `length` that include residue `pos`.
            first_start = max(1, pos - length + 1)
            last_start = min(pos, protein_length - length + 1)
            if last_start - first_start + 1 != length:
                raise ValueError(f"{gene}: mutation too close to an end for {length} windows")
            for start in range(first_start, last_start + 1):
                end = start + length - 1
                wt = seq[start - 1:end]
                index = pos - start  # zero-based within peptide
                mut = wt[:index] + new + wt[index + 1:]
                if not (len(wt) == len(mut) == length):
                    raise ValueError("Incorrect peptide length")
                diffs = [i for i, (a, b) in enumerate(zip(wt, mut)) if a != b]
                if diffs != [index] or wt[index] != old or mut[index] != new:
                    raise ValueError(f"{gene} {length} {start}: invalid mutant pairing")
                pair_id = f"{gene}_{m['variant']}_{length}aa_start{start:04d}"
                common: dict[str, object] = {
                    "pair_id": pair_id,
                    "gene": gene,
                    "variant": m["variant"],
                    "uniprot": acc,
                    "reference_protein": m["refseq_protein"],
                    "length": length,
                    "protein_start_1based": start,
                    "protein_end_1based": end,
                    "mutation_position_1based": pos,
                    "mutation_position_in_peptide_1based": index + 1,
                }
                for kind, peptide in (("WT", wt), ("MUT", mut)):
                    rows.append({
                        "peptide_id": f"{pair_id}_{kind}",
                        **common,
                        "type": kind,
                        "peptide": peptide,
                    })
                pairs.append({**common, "wt_peptide": wt, "mut_peptide": mut})

    assert len(pairs) == 108, f"Expected 108 WT/MUT windows, got {len(pairs)}"
    assert len(rows) == 216, f"Expected 216 peptide entries, got {len(rows)}"
    peptide_counts = Counter(str(r["peptide"]) for r in rows)
    duplicates = [s for s, n in peptide_counts.items() if n > 1]
    if duplicates:
        raise ValueError(
            f"Detected {len(duplicates)} repeated peptide sequence(s) across window IDs. "
            "Resolve sequence-to-window mapping before server submission. "
            f"Examples: {duplicates[:3]}"
        )

    # Every original Phase 1 pilot peptide must appear exactly in this full scan.
    full_set = {(r["gene"], r["variant"], str(r["length"]),
                 str(r["protein_start_1based"]), r["type"], r["peptide"])
                for r in rows}
    for p in pilot:
        key = (p["gene"], p["variant"], str(p["length"]),
               str(p["protein_start_1based"]), p["type"], p["peptide"])
        if key not in full_set:
            raise ValueError(f"Original pilot peptide absent from extension: {key}")

    # Write only after all reference/sequence/pilot validations have passed.
    RAW_EXT.mkdir(parents=True, exist_ok=True)
    RESULTS_EXT.mkdir(parents=True, exist_ok=True)
    save_csv(RESULTS_EXT / "peptide_manifest.csv", rows)
    save_csv(RESULTS_EXT / "pair_manifest.csv", pairs)
    for gene in GENES:
        seqs = [str(row["peptide"]) for row in rows if row["gene"] == gene]
        if len(seqs) != 108 or len(set(seqs)) != 108:
            raise ValueError(f"{gene}: expected 108 unique inputs")
        (RAW_EXT / f"{gene}_input.txt").write_text(
            "\n".join(seqs) + "\n", encoding="utf-8"
        )

    summary = [
        "PHASE 1 EXTENSION DESIGN VALIDATION PASSED",
        "Study status: exploratory extension (planned after inspecting pilot)",
        "Cancer context: HCC; mutations: CTNNB1 S45P, TP53 R249S",
        "Source: validated canonical UniProt FASTAs from original Phase 1",
        "Peptide window lengths: 15, 18, 21 aa",
        "Window rule: all peptide windows of given length spanning mutated residue",
        "Candidate interpretation: individual windows overlap; not independent mutations",
        "Expected per gene: 15+18+21=54 WT/MUT pairs, 108 unique sequences",
        "Total: 108 pairs; 216 peptide sequences; 2 HLA = 432 predictions",
        "Original 12 pilot peptide sequences: all preserved within extension",
        "HLA scenarios: DRB1_0101, DRB1_0401 (not patient HLA typing)",
        "Core filtering: none; later examine whether predicted 9-aa core contains mutation",
        "Outputs are prediction inputs only; experimental validation not performed",
        "Reference FASTA SHA256:",
        *fasta_hashes,
        "Files:",
        " raw/phase1_extension/CTNNB1_input.txt (108 sequences)",
        " raw/phase1_extension/TP53_input.txt (108 sequences)",
        " results/phase1_extension/peptide_manifest.csv (216 rows)",
        " results/phase1_extension/pair_manifest.csv (108 rows)",
    ]
    (RESULTS_EXT / "design_qc.txt").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )
    (RESULTS_EXT / "run_plan.md").write_text(
        "# Phase 1 exploratory extension — presubmission record\n\n"
        f"Prepared on: {date.today().isoformat()}\n\n"
        "The pilot was inspected before deciding on this exhaustive window extension. "
        "Do not describe this selection as prospectively prespecified before Phase 1 pilot.\n\n"
        "- Model: NetMHCIIpan-4.3k; input PEPTIDE; context encoding OFF\n"
        "- Predicted endpoints: EL+BA; report raw ranks without filtering\n"
        "- Alleles: DRB1_0101, DRB1_0401 (hypothetical)\n"
        "- Strong/weak threshold: 1% / 5%; report continuous ranks too\n"
        "- Inversion for non-DP: OFF; filtering: OFF; sort by score: OFF; XLS: ON\n"
        "- Two separate submissions, one per gene (108 peptides each)\n"
        "- Analyse ALL windows; do not choose only favourable windows post hoc\n"
        "- Group/interpret by mutation and allele; overlapping windows are dependent\n"
        "- Check mutation's presence in predicted 9-aa core in follow-up analysis\n"
        "- Candidate selection: consider absolute MUT ranks and WT/MUT deltas separately\n"
        "- Save full results as raw/phase1_extension/CTNNB1_output.txt, TP53_output.txt\n"
        "- Save XLS files in same directory; keep originals unchanged\n"
        "- Do not infer experimental binding, natural presentation, or CD4 response\n",
        encoding="utf-8",
    )
    print("\n".join(summary))


if __name__ == "__main__":
    main()
