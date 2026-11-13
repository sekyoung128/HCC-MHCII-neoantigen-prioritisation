#!/usr/bin/env python3
"""Phase 1 HCC exploratory extension: audit and analyse 432 NetMHCIIpan predictions.

Usage (from anywhere): python phase1_extension_analyse.py
Requires only existing project files from phase1_expand_windows.py plus the
CTNNB1_output.txt and TP53_output.txt copied from the NetMHCIIpan web results.
Outputs are saved under results/phase1_extension/analysis/.
No network access, experiments, or new inference is performed by this script.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
EXT_RAW = RAW / "phase1_extension"
EXT_RESULTS = ROOT / "results" / "phase1_extension"
OUT = EXT_RESULTS / "analysis"
FIG = OUT / "figures"
GENES = ("CTNNB1", "TP53")
ALLELES = ("DRB1_0101", "DRB1_0401")
LENGTHS = (15, 18, 21)
AA = set("ACDEFGHIKLMNPQRSTVWY")
EXPECTED = {"CTNNB1": ("S45P", 45, "S", "P", "P35222", 781),
            "TP53": ("R249S", 249, "R", "S", "P04637", 393)}


def require(condition, reason):
    if not condition:
        raise ValueError("QC FAILED: " + str(reason))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_fasta(path, accession):
    require(path.is_file(), f"Reference FASTA missing: {path}")
    lines = [s.strip() for s in path.read_text(encoding="utf-8").splitlines() if s.strip()]
    require(len(lines) > 1 and lines[0].startswith(">") and accession in lines[0],
            f"Wrong FASTA header/accession: {path}")
    require(not any(x.startswith(">") for x in lines[1:]), "Multiple FASTA entries")
    seq = "".join(lines[1:])
    require(set(seq).issubset(AA), f"Invalid amino-acid symbols in {path}")
    return seq


def get_design():
    p = EXT_RESULTS / "peptide_manifest.csv"
    pp = EXT_RESULTS / "pair_manifest.csv"
    require(p.is_file() and pp.is_file(), "Missing manifest(s); run phase1_expand_windows.py")
    d = pd.read_csv(p, keep_default_na=False)
    pairs = pd.read_csv(pp, keep_default_na=False)
    needed = {"peptide_id", "pair_id", "gene", "variant", "uniprot", "length",
              "protein_start_1based", "protein_end_1based", "mutation_position_1based",
              "mutation_position_in_peptide_1based", "type", "peptide"}
    require(needed.issubset(d.columns), f"Manifest missing: {needed - set(d.columns)}")
    require(len(d) == 216 and len(pairs) == 108, "Expected 216 peptide / 108 pair rows")
    require(not d.peptide_id.duplicated().any(), "Duplicated peptide ID")
    require(not d.peptide.duplicated().any(),
            "Repeated sequence: cannot uniquely attribute NetMHCIIpan output to windows")
    require(set(d.gene) == set(GENES) and set(d.type) == {"WT", "MUT"},
            "Unexpected gene or WT/MUT label")
    require(set(d.length) == set(LENGTHS), "Unexpected peptide length")
    require(set(pairs.pair_id) == set(d.pair_id) and not pairs.pair_id.duplicated().any(),
            "pair_manifest IDs inconsistent with peptide_manifest")

    for gene in GENES:
        variant, pos, old, new, accession, size = EXPECTED[gene]
        seq = load_fasta(RAW / f"{gene}_{accession}.fasta", accession)
        require(len(seq) == size and seq[pos-1] == old,
                f"{gene} canonical sequence length / residue mismatch")
        dg = d[d.gene == gene]
        require(len(dg) == 108, f"{gene}: wrong manifest sequence count")
        require(set(dg.variant) == {variant} and set(dg.uniprot) == {accession},
                f"{gene}: mutation / UniProt metadata mismatch")
        for length in LENGTHS:
            dd = dg[dg.length == length]
            require(len(dd) == length * 2 and dd.pair_id.nunique() == length,
                    f"{gene} {length}-mer: expected {length} overlapping windows")
        for row in dg.itertuples(index=False):
            start = int(row.protein_start_1based)
            end = int(row.protein_end_1based)
            within = int(row.mutation_position_in_peptide_1based)
            require(end-start+1 == int(row.length) and 1 <= start <= pos <= end <= size,
                    f"Invalid reference window: {row.peptide_id}")
            require(int(row.mutation_position_1based) == pos and within == pos-start+1,
                    f"Incorrect mutation offset: {row.peptide_id}")
            wt = seq[start-1:end]
            expected_pep = wt if row.type == "WT" else wt[:within-1]+new+wt[within:]
            require(row.peptide == expected_pep and len(row.peptide) == row.length,
                    f"Reference/peptide mismatch: {row.peptide_id}")

    for pair_id, g in d.groupby("pair_id"):
        require(len(g) == 2 and set(g.type) == {"WT", "MUT"},
                f"Incomplete design WT/MUT pair: {pair_id}")
        wt = g.set_index("type").loc["WT"]
        mut = g.set_index("type").loc["MUT"]
        mutation_i = int(wt.mutation_position_in_peptide_1based)-1
        diffs = [i for i, (x, y) in enumerate(zip(wt.peptide, mut.peptide)) if x != y]
        require(diffs == [mutation_i], f"Pair not a single mutation: {pair_id}")
        p = pairs.loc[pairs.pair_id == pair_id].iloc[0]
        require(p.wt_peptide == wt.peptide and p.mut_peptide == mut.peptide,
                f"Pair manifest mismatch: {pair_id}")

    for gene in GENES:
        input_file = EXT_RAW / f"{gene}_input.txt"
        require(input_file.is_file(), f"Input list missing: {input_file}")
        input_seqs = [line.strip() for line in input_file.read_text(encoding="utf-8").splitlines()
                      if line.strip()]
        expected_seqs = d.loc[d.gene == gene, "peptide"].tolist()
        require(input_seqs == expected_seqs,
                f"{gene} input file differs from manifest (including order)")
    return d, pairs


def parse_output(gene):
    path = EXT_RAW / f"{gene}_output.txt"
    require(path.is_file(), f"Server result missing: {path}")
    raw = path.read_text(encoding="utf-8")
    require("NetMHCIIpan version 4.3k" in raw,
            f"{gene}: wrong or missing predictor version")
    require("Input is in PEPTIDE format" in raw,
            f"{gene}: wrong input format / context encoding setting")
    require("Prediction Mode: EL+BA" in raw,
            f"{gene}: EL+BA not confirmed")
    seen_alleles = re.findall(r"^# Allele:\s*(\S+)", raw, re.MULTILINE)
    require(set(seen_alleles) == set(ALLELES) and len(seen_alleles) == 2,
            f"{gene}: expected exactly two allele sections, got {seen_alleles}")
    records = []
    for ln, line in enumerate(raw.splitlines(), start=1):
        if not re.match(r"^\s*\d+\s+DRB1_", line):
            continue
        f = line.split()
        require(len(f) >= 14, f"{gene} line {ln}: incomplete prediction row: {line[:110]}")
        level = " ".join(f[14:])
        require(level in ("", "SB", "WB", "<= SB", "<= WB"),
                f"{gene} line {ln}: unknown trailing BindLevel: {level}")
        try:
            row = {"Pos": int(f[0]), "MHC": f[1], "Peptide": f[2],
                   "Of": int(f[3]), "Core": f[4], "Core_Rel": float(f[5]),
                   "Inverted": int(f[6]), "Identity": f[7], "Score_EL": float(f[8]),
                   "Rank_EL": float(f[9]), "Exp_Bind": f[10], "Score_BA": float(f[11]),
                   "Rank_BA": float(f[12]), "Affinity_nM": float(f[13]),
                   "BindLevel": level, "source_gene": gene, "source_line": ln}
        except (ValueError, IndexError) as e:
            raise ValueError(f"QC FAILED: {gene} line {ln}: cannot parse {line}") from e
        records.append(row)
    p = pd.DataFrame.from_records(records)
    require(len(p) == 216, f"{gene}: expected 216 predictions; found {len(p)}")
    require(p.groupby("MHC").size().to_dict() == {a: 108 for a in ALLELES},
            f"{gene}: expected 108 rows for each allele")
    require(not p[["MHC", "Peptide"]].duplicated().any(),
            f"{gene}: duplicated peptide–HLA predictions")
    for fld in ("Score_BA", "Score_EL", "Rank_BA", "Rank_EL", "Affinity_nM", "Core_Rel"):
        require(np.isfinite(p[fld].to_numpy()).all(),
                f"{gene}: invalid/missing {fld}")
    for fld in ("Rank_BA", "Rank_EL"):
        require(p[fld].between(0, 100).all(), f"{gene}: percentile rank outside 0–100")
    require((p.Affinity_nM > 0).all(), f"{gene}: nonpositive predicted affinity")
    require(p.Core_Rel.between(0, 1).all(), f"{gene}: Core_Rel outside [0,1]")
    require(p.Core.str.len().eq(9).all(), f"{gene}: core must be 9 aa")
    for r in p.itertuples(index=False):
        require(0 <= r.Of <= len(r.Peptide)-9 and r.Peptide[r.Of:r.Of+9] == r.Core,
                f"{gene} line {r.source_line}: predicted core does not match Of / sequence")
    return p, path


def make_figures(comp):
    FIG.mkdir(parents=True, exist_ok=True)
    files = []
    for gene in GENES:
        g = comp[comp.gene == gene]
        for endpoint in ("BA", "EL"):
            # Each heatmap covers every 15-, 18- and 21-mer window, allele by allele.
            grid = np.full((6, 21), np.nan)
            row_labels = []
            for idx, (allele, length) in enumerate((a, L) for a in ALLELES for L in LENGTHS):
                row_labels.append(f"{allele} | {length} aa")
                group = g[(g.MHC == allele) & (g.length == length)]
                require(len(group) == length, f"Heatmap grid incomplete: {gene} {allele} {length}")
                for r in group.itertuples(index=False):
                    col = int(r.mutation_position_in_peptide_1based) - 1
                    grid[idx, col] = getattr(r, f"Rank_{endpoint}_MUT")
            fig, ax = plt.subplots(figsize=(10, 4.5))
            mat = ax.imshow(np.ma.masked_invalid(grid), aspect="auto", vmin=0, vmax=100)
            ax.set_xticks(range(21), labels=[str(x) for x in range(1, 22)])
            ax.set_yticks(range(6), labels=row_labels)
            ax.set_xlabel("Position of mutated residue within peptide (1-based)")
            ax.set_title(f"{gene} | mutant %Rank_{endpoint} | all mutation-spanning windows")
            fig.colorbar(mat, ax=ax, label=f"Mutant %Rank_{endpoint} (lower = more favourable)")
            fig.tight_layout()
            file = FIG / f"{gene.lower()}_mutant_{endpoint.lower()}_heatmap.png"
            fig.savefig(file, dpi=200)
            plt.close(fig)
            files.append(file)

        fig, ax = plt.subplots(figsize=(7, 5.5))
        for allele, group in g.groupby("MHC"):
            ax.scatter(group.Delta_BA_WT_minus_MUT, group.Delta_EL_WT_minus_MUT,
                       label=allele, s=25, alpha=0.65)
        ax.axhline(0, linestyle="--", linewidth=1)
        ax.axvline(0, linestyle="--", linewidth=1)
        ax.set_xlabel("WT %Rank_BA − MUT %Rank_BA")
        ax.set_ylabel("WT %Rank_EL − MUT %Rank_EL")
        ax.set_title(f"{gene} | mutation-associated BA vs EL shifts")
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        file = FIG / f"{gene.lower()}_delta_ba_vs_el.png"
        fig.savefig(file, dpi=200)
        plt.close(fig)
        files.append(file)

        fig, ax = plt.subplots(figsize=(7, 5.5))
        for allele, group in g.groupby("MHC"):
            ax.scatter(group.Rank_BA_MUT, group.Delta_BA_WT_minus_MUT,
                       label=allele, s=25, alpha=0.65)
        ax.axhline(0, linestyle="--", linewidth=1)
        ax.axvline(5, linestyle=":", linewidth=1)
        ax.set_xlabel("Absolute mutant %Rank_BA (lower = more favourable)")
        ax.set_ylabel("WT %Rank_BA − MUT %Rank_BA (positive = improved mutant rank)")
        ax.set_title(f"{gene} | absolute BA rank vs mutation effect")
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        file = FIG / f"{gene.lower()}_absolute_vs_delta_ba.png"
        fig.savefig(file, dpi=200)
        plt.close(fig)
        files.append(file)
    return files


def main():
    design, _ = get_design()
    frames = []
    hashes = {}
    for gene in GENES:
        frame, path = parse_output(gene)
        frames.append(frame)
        hashes[path.name] = sha(path)
        observed = set(frame.Peptide)
        expected = set(design.loc[design.gene == gene, "peptide"])
        require(observed == expected,
                f"{gene}: output peptide set differs from submitted manifest; "
                f"missing={len(expected-observed)}, extra={len(observed-expected)}")

    predictions = pd.concat(frames, ignore_index=True)
    require(len(predictions) == 432, "Wrong combined output row count")
    annotated = predictions.merge(design, left_on="Peptide", right_on="peptide",
                                  how="left", indicator=True, validate="many_to_one")
    require(annotated._merge.eq("both").all(), "Unmapped model result")
    require(annotated.source_gene.eq(annotated.gene).all(),
            "Gene label inconsistent between file and manifest")
    annotated = annotated.drop(columns=["_merge", "peptide"])
    require(annotated.peptide_id.nunique() == 216 and annotated.MHC.nunique() == 2,
            "Incomplete peptide/allele coverage")
    require(not annotated[["peptide_id", "MHC"]].duplicated().any(),
            "Duplicated peptide–allele mapping")

    keys = ["pair_id", "gene", "variant", "length", "protein_start_1based",
            "protein_end_1based", "mutation_position_1based",
            "mutation_position_in_peptide_1based", "MHC"]
    values = ["peptide_id", "Peptide", "Rank_BA", "Rank_EL", "Score_BA",
              "Score_EL", "Affinity_nM", "Core", "Of", "Core_Rel"]
    wt = annotated.loc[annotated.type == "WT", keys+values]
    mut = annotated.loc[annotated.type == "MUT", keys+values]
    comp = wt.merge(mut, on=keys, how="outer", suffixes=("_WT", "_MUT"),
                    validate="one_to_one", indicator=True)
    require(len(comp) == 216 and comp._merge.eq("both").all(),
            "Expected 216 complete WT/MUT–HLA pairs")
    comp = comp.drop(columns="_merge")
    for endpoint in ("BA", "EL"):
        comp[f"Delta_{endpoint}_WT_minus_MUT"] = (
            comp[f"Rank_{endpoint}_WT"]-comp[f"Rank_{endpoint}_MUT"])
        comp[f"MUT_{endpoint}_rank_le_1"] = comp[f"Rank_{endpoint}_MUT"].le(1)
        comp[f"MUT_{endpoint}_rank_le_5"] = comp[f"Rank_{endpoint}_MUT"].le(5)
    comp["Mutation_in_core_WT"] = (
        (comp.mutation_position_in_peptide_1based-1 >= comp.Of_WT) &
        (comp.mutation_position_in_peptide_1based-1 < comp.Of_WT+9))
    comp["Mutation_in_core_MUT"] = (
        (comp.mutation_position_in_peptide_1based-1 >= comp.Of_MUT) &
        (comp.mutation_position_in_peptide_1based-1 < comp.Of_MUT+9))
    comp["Predicted_core_offset_changed"] = comp.Of_WT.ne(comp.Of_MUT)
    ba = comp.Delta_BA_WT_minus_MUT
    el = comp.Delta_EL_WT_minus_MUT
    comp["BA_EL_shift_direction"] = np.select(
        [(ba>0)&(el>0), (ba<0)&(el<0), (ba>0)&(el<0), (ba<0)&(el>0)],
        ["both favour MUT", "both favour WT", "BA favours MUT; EL favours WT",
         "BA favours WT; EL favours MUT"], default="one or both unchanged")
    comp = comp.sort_values(["gene", "MHC", "length", "protein_start_1based"])

    # Optional reproducibility check: the 12 original pilot peptides should receive
    # consistent scores under the same model build, alleles and settings.
    pilot_note = "Pilot prediction cross-check: unavailable (no annotated pilot CSV)"
    pilot_file = ROOT / "results" / "phase1_predictions_annotated.csv"
    pilot_diff = None
    if pilot_file.is_file():
        pilot = pd.read_csv(pilot_file)
        pilot_needed = {"Peptide", "MHC", "Rank_BA", "Rank_EL"}
        require(pilot_needed.issubset(pilot.columns) and len(pilot) == 24,
                "Original pilot prediction CSV missing expected columns / rows")
        pilot_diff = pilot[["Peptide", "MHC", "Rank_BA", "Rank_EL"]].merge(
            annotated[["Peptide", "MHC", "Rank_BA", "Rank_EL"]],
            on=["Peptide", "MHC"], how="left", validate="one_to_one",
            indicator=True, suffixes=("_pilot", "_extension"))
        require(pilot_diff._merge.eq("both").all(),
                "Some Phase 1 pilot peptide–HLA rows missing from expanded model results")
        for endpoint in ("BA", "EL"):
            pilot_diff[f"abs_difference_{endpoint}"] = (
                pilot_diff[f"Rank_{endpoint}_pilot"] -
                pilot_diff[f"Rank_{endpoint}_extension"]).abs()
        max_diff = pilot_diff[["abs_difference_BA", "abs_difference_EL"]].max().max()
        require(max_diff <= 0.05,
                f"Pilot/extension rank predictions disagree (max difference {max_diff:.3f}). "
                "Check that both server runs used the same build and settings.")
        pilot_note = ("Pilot prediction cross-check: 24/24 overlapping peptide–HLA rows "
                      f"consistent (largest rank difference {max_diff:.3f})")
        pilot_diff = pilot_diff.drop(columns="_merge")

    summary_rows = []
    for (gene, hla, length), group in comp.groupby(["gene", "MHC", "length"], sort=True):
        summary_rows.append({
            "gene": gene, "MHC": hla, "length": length, "window_pairs": len(group),
            "best_mutant_BA_rank": group.Rank_BA_MUT.min(),
            "best_mutant_EL_rank": group.Rank_EL_MUT.min(),
            "median_mutant_BA_rank": group.Rank_BA_MUT.median(),
            "median_mutant_EL_rank": group.Rank_EL_MUT.median(),
            "largest_improvement_delta_BA": group.Delta_BA_WT_minus_MUT.max(),
            "largest_improvement_delta_EL": group.Delta_EL_WT_minus_MUT.max(),
            "n_mutant_BA_rank_le_1": int(group.MUT_BA_rank_le_1.sum()),
            "n_mutant_BA_rank_le_5": int(group.MUT_BA_rank_le_5.sum()),
            "n_mutant_EL_rank_le_1": int(group.MUT_EL_rank_le_1.sum()),
            "n_mutant_EL_rank_le_5": int(group.MUT_EL_rank_le_5.sum()),
            "n_both_BA_EL_rank_le_5": int((group.MUT_BA_rank_le_5 & group.MUT_EL_rank_le_5).sum()),
            "n_mutation_in_mutant_core": int(group.Mutation_in_core_MUT.sum()),
            "n_core_offset_changed": int(group.Predicted_core_offset_changed.sum()),
        })
    summary = pd.DataFrame(summary_rows)
    require(len(summary) == 12 and summary.window_pairs.sum() == 216,
            "Summary groups inconsistent")

    # Only write outputs after input and prediction QC have passed.
    OUT.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(OUT / "predictions_annotated_432.csv", index=False)
    comp.to_csv(OUT / "wt_mut_comparison_216.csv", index=False)
    summary.to_csv(OUT / "window_summary_by_gene_hla_length.csv", index=False)
    if pilot_diff is not None:
        pilot_diff.to_csv(OUT / "pilot_vs_extension_consistency.csv", index=False)
    ranked_columns = ["gene", "variant", "MHC", "length", "protein_start_1based",
                      "mutation_position_in_peptide_1based", "Peptide_MUT", "Peptide_WT",
                      "Rank_BA_MUT", "Rank_EL_MUT", "Delta_BA_WT_minus_MUT",
                      "Delta_EL_WT_minus_MUT", "Mutation_in_core_MUT",
                      "Predicted_core_offset_changed"]
    for name, sorter, ascending in [
        ("top10_lowest_mutant_ba_rank.csv", "Rank_BA_MUT", True),
        ("top10_lowest_mutant_el_rank.csv", "Rank_EL_MUT", True),
        ("top10_largest_ba_rank_improvement.csv", "Delta_BA_WT_minus_MUT", False),
    ]:
        comp.sort_values(sorter, ascending=ascending).head(10)[ranked_columns].to_csv(
            OUT / name, index=False)
    figures = make_figures(comp)

    overall = {
        "MUT BA <= 1%": int(comp.MUT_BA_rank_le_1.sum()),
        "MUT BA <= 5%": int(comp.MUT_BA_rank_le_5.sum()),
        "MUT EL <= 1%": int(comp.MUT_EL_rank_le_1.sum()),
        "MUT EL <= 5%": int(comp.MUT_EL_rank_le_5.sum()),
        "Both MUT BA and EL <= 5%": int((comp.MUT_BA_rank_le_5 & comp.MUT_EL_rank_le_5).sum()),
        "BA delta positive": int((ba > 0).sum()),
        "EL delta positive": int((el > 0).sum()),
        "Predicted core offset changed": int(comp.Predicted_core_offset_changed.sum()),
        "Mutation inside predicted MUT core": int(comp.Mutation_in_core_MUT.sum()),
    }
    report = [
        "PHASE 1 EXTENSION PREDICTION VALIDATION PASSED",
        "Status: post-pilot exploratory extension; NOT independent predictor validation",
        "Predictor: NetMHCIIpan-4.3k | PEPTIDE | EL+BA (per output header)",
        "Source variants: HCC-associated CTNNB1 S45P, TP53 R249S",
        "HLA: DRB1_0101; DRB1_0401 (hypothetical; no patient HLA typing)",
        f"Model rows: {len(predictions)}; unique peptides: {annotated.peptide_id.nunique()}",
        f"WT/MUT-HLA matched comparisons: {len(comp)}",
        f"Summary strata (gene x HLA x length): {len(summary)}",
        pilot_note,
        *[f"{label}: {n} of 216" for label, n in overall.items()],
        f"Minimum absolute MUT BA %Rank: {comp.Rank_BA_MUT.min():.2f}",
        f"Minimum absolute MUT EL %Rank: {comp.Rank_EL_MUT.min():.2f}",
        "Input SHA-256:",
        *[f"  {gene}_input.txt: {sha(EXT_RAW / f'{gene}_input.txt')}" for gene in GENES],
        "Output SHA-256:",
        *[f"  {k}: {v}" for k, v in hashes.items()],
        "Interpretation: lower percentile rank means more favourable predictor score;",
        "  delta = WT rank minus MUT rank (positive = lower MUT rank).",
        "  BA and EL endpoints differ; do not equate shifts with immunogenicity.",
        "  Many windows overlap; 216 pairs are NOT 216 independent mutations.",
        "  These windows/HLA contexts were explored after inspecting the pilot.",
        "  No measured binding, natural presentation or CD4 functional confirmation.",
        f"Figures created: {len(figures)}",
    ]
    (OUT / "analysis_qc_and_results.txt").write_text("\n".join(report)+"\n", encoding="utf-8")
    print("\n".join(report[:16]))
    print("\nSaved outputs to:", OUT)


if __name__ == "__main__":
    main()
