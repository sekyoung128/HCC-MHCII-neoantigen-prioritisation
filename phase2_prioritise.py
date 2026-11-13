#!/usr/bin/env python3
"""Locked Phase 2 descriptive candidate-decision comparison.

Run in the mhcii-neoantigen-discovery project root:
    python phase2_prioritise.py
Optionally for testing: --input PATH --rules PATH --out DIRECTORY

Only uses pre-existing Phase 1 predictions. No experimental labels or internet.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "results/phase1_extension/analysis/wt_mut_comparison_216.csv"
DEFAULT_RULES = ROOT / "phase2_rules.json"
DEFAULT_OUT = ROOT / "results/phase2"


def check(ok, message):
    if not bool(ok):
        raise ValueError("PHASE 2 QC FAILED: " + message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return p.parse_args()


def validate(df, rules):
    mandatory = [
        "pair_id", "gene", "variant", "MHC", "length", "protein_start_1based",
        "protein_end_1based", "mutation_position_1based", "mutation_position_in_peptide_1based",
        "Peptide_WT", "Peptide_MUT", "Core_WT", "Core_MUT", "Of_WT", "Of_MUT",
        "Rank_BA_WT", "Rank_BA_MUT", "Rank_EL_WT", "Rank_EL_MUT",
        "Delta_BA_WT_minus_MUT", "Delta_EL_WT_minus_MUT",
        "Mutation_in_core_MUT", "Predicted_core_offset_changed",
    ]
    check(all(x in df for x in mandatory), "Required columns missing")
    check(len(df) == rules["expected_rows"], "Unexpected input row count")
    check(df[mandatory].notna().all().all(), "Missing value in required columns")
    check(set(df.gene) == set(rules["expected_gene_variants"]), "Unexpected gene set")
    for gene, variant in rules["expected_gene_variants"].items():
        check(set(df.loc[df.gene == gene, "variant"]) == {variant}, "Incorrect mutation " + gene)
    check(set(df.MHC) == set(rules["expected_HLA"]), "Unexpected HLA alleles")
    check(set(df.length) == set(rules["expected_lengths"]), "Unexpected peptide lengths")
    for (gene, hla, leng), sub in df.groupby(["gene", "MHC", "length"]):
        check(len(sub) == int(leng), f"Wrong number of windows: {gene}/{hla}/{leng}")
    check(not df.duplicated(["pair_id", "MHC"]).any(), "Duplicate candidate keys")
    check((df.groupby("pair_id").MHC.nunique() == len(rules["expected_HLA"])).all(), "Missing HLA partner for a window")
    for col in ["Rank_BA_WT", "Rank_BA_MUT", "Rank_EL_WT", "Rank_EL_MUT"]:
        check(pd.to_numeric(df[col], errors="coerce").between(0, 100).all(), "Invalid rank: " + col)
    check(np.allclose(df.Delta_BA_WT_minus_MUT, df.Rank_BA_WT-df.Rank_BA_MUT, atol=1e-8), "Wrong BA delta")
    check(np.allclose(df.Delta_EL_WT_minus_MUT, df.Rank_EL_WT-df.Rank_EL_MUT, atol=1e-8), "Wrong EL delta")
    for r in df.itertuples():
        check(len(r.Peptide_WT) == len(r.Peptide_MUT) == r.length, "Invalid peptide length")
        check(r.protein_end_1based-r.protein_start_1based+1 == r.length, "Invalid source-protein window")
        at = r.mutation_position_in_peptide_1based - 1
        check(0 <= at < r.length, "Invalid mutation location")
        check(r.protein_start_1based + at == r.mutation_position_1based, "Mutation position mismatch")
        diffs = [j for j,(a,b) in enumerate(zip(r.Peptide_WT, r.Peptide_MUT)) if a != b]
        check(diffs == [at], "WT/MUT mismatch not limited to mutation position")
        for version in ("WT", "MUT"):
            seq = getattr(r, "Peptide_" + version)
            core = getattr(r, "Core_" + version)
            off = int(getattr(r, "Of_" + version))
            check(0 <= off <= len(seq) - 9 and seq[off:off+9] == core,
                  "Core and core-offset mismatch")
        expected_in_core = int(r.Of_MUT) <= at < int(r.Of_MUT) + 9
        expected_offset_change = int(r.Of_MUT) != int(r.Of_WT)
        check(str(r.Mutation_in_core_MUT).lower() == str(expected_in_core).lower(), "Mutation-in-core flag mismatch")
        check(str(r.Predicted_core_offset_changed).lower() == str(expected_offset_change).lower(), "Core offset change flag mismatch")


def add_features(df):
    d = df.copy()
    d["candidate_id"] = d.pair_id + "|" + d.MHC
    d["deltaBA"] = d.Delta_BA_WT_minus_MUT
    d["deltaEL"] = d.Delta_EL_WT_minus_MUT
    d["improved_BA"] = d.deltaBA > 0
    d["improved_EL"] = d.deltaEL > 0
    d["opposite_BA_EL_delta_sign"] = (d.deltaBA * d.deltaEL) < 0
    d["abs_delta_endpoint_gap"] = (d.deltaBA-d.deltaEL).abs()
    d["min_mutant_rank"] = d[["Rank_BA_MUT", "Rank_EL_MUT"]].min(axis=1)
    d["core_group"] = d.gene + "|" + d.MHC + "|" + d.Core_MUT
    # HLA contrasts are within exactly the same mutation/peptide window.
    pivot_rank = d.pivot(index="pair_id", columns="MHC", values="Rank_BA_MUT")
    pivot_delta = d.pivot(index="pair_id", columns="MHC", values="deltaBA")
    a,b = sorted(pivot_rank.columns)
    hla_sensitivity = pd.DataFrame({
        "pair_id":pivot_rank.index,
        "HLA_mut_BA_rank_range":(pivot_rank[a]-pivot_rank[b]).abs().to_numpy(),
        "HLA_BA_delta_sign_opposite":((pivot_delta[a] * pivot_delta[b]) < 0).to_numpy(),
    })
    d = d.merge(hla_sensitivity,on="pair_id",validate="many_to_one",how="left")
    check(d.HLA_mut_BA_rank_range.notna().all(), "HLA sensitivity join failed")
    return d


def get_orders(d):
    # Locked lexicographic algorithms; no 'vaccine suitability' score.
    a = d.sort_values(["Rank_BA_MUT", "Rank_EL_MUT", "candidate_id"],
                      ascending=[True, True, True], kind="mergesort")
    b = d.sort_values(["improved_BA", "Rank_BA_MUT", "deltaBA", "candidate_id"],
                      ascending=[False, True, False, True], kind="mergesort")
    q1 = d[d.opposite_BA_EL_delta_sign].sort_values(
        ["abs_delta_endpoint_gap", "min_mutant_rank", "candidate_id"],
        ascending=[False, True, True], kind="mergesort")
    q2 = d[d.Predicted_core_offset_changed].sort_values(
        ["Mutation_in_core_MUT", "Rank_BA_MUT", "abs_delta_endpoint_gap", "candidate_id"],
        ascending=[False, True, False, True], kind="mergesort")
    q3 = d.sort_values(
        ["HLA_BA_delta_sign_opposite", "HLA_mut_BA_rank_range", "Rank_BA_MUT", "candidate_id"],
        ascending=[False, False, True, True], kind="mergesort")
    fallback = d.sort_values(["Rank_EL_MUT", "candidate_id"],ascending=[True, True],kind="mergesort")
    queues = [("BA_EL_mutation_shift_conflict", list(q1.candidate_id)),
              ("WT_MUT_core_offset_shift", list(q2.candidate_id)),
              ("cross_HLA_BA_sensitivity", list(q3.candidate_id))]
    seen, c_order, reason = set(), [], {}
    positions = [0,0,0]
    while len(c_order) < len(d):
        progressed = False
        for i,(label,items) in enumerate(queues):
            while positions[i] < len(items) and items[positions[i]] in seen:
                positions[i]+=1
            if positions[i] < len(items):
                cid = items[positions[i]]
                positions[i]+=1
                if cid not in seen:
                    seen.add(cid);c_order.append(cid);reason[cid]=label
                    progressed=True
        if not progressed: break
    for cid in fallback.candidate_id:
        if cid not in seen:
            seen.add(cid);c_order.append(cid);reason[cid]="EL_rank_fallback"
    check(len(c_order) == len(d) and len(seen) == len(d), "C ranking does not cover candidate pool")
    return {"A": list(a.candidate_id), "B": list(b.candidate_id), "C":c_order},reason


def save_plots(d, ranks, out, ks):
    # Each figure is standalone; use matplotlib default palette and markers.
    jrows=[]
    for k in ks:
        selections={s:set(v[:k]) for s,v in ranks.items()}
        for s1,s2 in combinations("ABC",2):
            union=selections[s1]|selections[s2]
            jrows.append({"K":k,"strategies":s1+" vs "+s2,
                          "Jaccard":len(selections[s1]&selections[s2])/len(union)})
    j=pd.DataFrame(jrows)
    fig,ax=plt.subplots(figsize=(7,4.5))
    for lab,g in j.groupby("strategies"):
        ax.plot(g.K,g.Jaccard,marker="o",label=lab)
    ax.set(ylim=(-0.02,1.02),xlabel="Top-K",ylabel="Jaccard overlap",title="Candidate decision overlap by strategy")
    ax.set_xticks(ks);ax.grid(alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(out/"fig1_topk_jaccard.png",dpi=220);plt.close(fig)

    fig,ax=plt.subplots(figsize=(6,6))
    ax.scatter(d.rank_A,d.rank_B,s=22,alpha=.55)
    ax.plot([1,len(d)],[1,len(d)],linestyle="--",linewidth=1,alpha=.6)
    ax.set(xlabel="A: BA-first ranking",ylabel="B: mutation-aware ranking",
           title="Global candidate rank changes (216 windows x HLA)")
    ax.grid(alpha=.18);fig.tight_layout()
    fig.savefig(out/"fig2_rank_A_vs_B.png",dpi=220);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,5.5))
    ax.scatter(d.Rank_BA_MUT,d.deltaBA,marker=".",label="all candidates",alpha=.45)
    for strat,marker in [("A","o"),("B","^"),("C","s")]:
        sub=d[d[f"rank_{strat}"] <= 10]
        ax.scatter(sub.Rank_BA_MUT,sub.deltaBA,marker=marker,
                   facecolors="none",s=80,label=f"Top-10 {strat}")
    ax.axhline(0,linestyle="--",linewidth=.9)
    ax.set(xlabel="Mutant BA percentile rank (lower = more favourable)",
           ylabel="WT rank minus MUT rank (BA; positive = mutant-favourable)",
           title="Absolute BA rank versus mutation-associated shift")
    ax.legend(fontsize=9);ax.grid(alpha=.18);fig.tight_layout()
    fig.savefig(out/"fig3_absolute_ba_vs_delta.png",dpi=220);plt.close(fig)

    union=[]
    for strat in "ABC":
        union.extend(ranks[strat][:10])
    union=list(dict.fromkeys(union))
    display = d.set_index("candidate_id").loc[union]
    fig,ax=plt.subplots(figsize=(8,max(5,len(union)*.32+1)))
    selection=np.array([[int(display.loc[c,f"rank_{s}"]<=10) for s in "ABC"] for c in union])
    ax.imshow(selection,aspect="auto",interpolation="nearest",vmin=0,vmax=1)
    ax.set_xticks(range(3),["A · BA", "B · mutation", "C · uncertainty"])
    ax.set_yticks(range(len(union)),[f"{display.loc[c,'gene']} {display.loc[c,'MHC']} "
               f"{int(display.loc[c,'length'])}aa start{int(display.loc[c,'protein_start_1based'])}" for c in union])
    ax.set_title("Top-10 selection membership (bright = selected)")
    for i in range(len(union)):
        for j in range(3):
            ax.text(j,i,"●" if selection[i,j] else "",ha="center",va="center",fontsize=10)
    fig.tight_layout();fig.savefig(out/"fig4_top10_membership.png",dpi=220);plt.close(fig)


def run(args):
    check(args.input.is_file(),f"Missing input CSV: {args.input}")
    check(args.rules.is_file(),f"Missing rule file: {args.rules}")
    rules=json.loads(args.rules.read_text(encoding="utf-8"))
    current_hash=digest(args.input)
    check(current_hash==rules["source_sha256"],
          f"Dataset SHA-256 has changed. Expected {rules['source_sha256']}; got {current_hash}. Investigate before running.")
    ks=[int(x) for x in rules["primary_topk"]]
    check(ks==[5,10,20],"Top-K parameters differ from locked method")
    df=pd.read_csv(args.input)
    validate(df,rules)
    d=add_features(df)
    ranks, reasons=get_orders(d)
    for s,ids in ranks.items():
        rankmap={cid:i+1 for i,cid in enumerate(ids)}
        d[f"rank_{s}"]=d.candidate_id.map(rankmap).astype(int)
    d["C_selection_question"]=d.candidate_id.map(reasons)
    check(len(d)==216,"Candidate count changed during joins")
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    d.sort_values("rank_A").to_csv(out/"phase2_all_216_ranked.csv",index=False)

    toplists=[]
    for s in "ABC":
        top=d[d[f"rank_{s}"]<=10].sort_values(f"rank_{s}").copy()
        top.insert(0,"selected_rank",top[f"rank_{s}"])
        top.insert(0,"strategy",s)
        top.to_csv(out/f"phase2_top10_{s}.csv",index=False)
        toplists.append(top)
    pd.concat(toplists,ignore_index=True).to_csv(out/"phase2_top10_all_strategies.csv",index=False)

    overlaps=[];all_membership=[]
    for k in ks:
        sets={s:set(ranks[s][:k]) for s in "ABC"}
        for s1,s2 in combinations("ABC",2):
            i=len(sets[s1]&sets[s2]); u=len(sets[s1]|sets[s2])
            overlaps.append({"K":k,"comparison":s1+"_vs_"+s2,
                             "intersection":i,"union":u,"Jaccard":i/u})
        overlaps.append({"K":k,"comparison":"A_B_C_all_three",
                         "intersection":len(sets["A"]&sets["B"]&sets["C"]),
                         "union":len(sets["A"]|sets["B"]|sets["C"]),
                         "Jaccard":len(sets["A"]&sets["B"]&sets["C"])/len(sets["A"]|sets["B"]|sets["C"])})
        for row in d.itertuples(index=False):
            all_membership.append({"K":k,"candidate_id":row.candidate_id,
                        "A":row.rank_A<=k,"B":row.rank_B<=k,"C":row.rank_C<=k})
    pd.DataFrame(overlaps).to_csv(out/"phase2_topK_overlap.csv",index=False)
    pd.DataFrame(all_membership).to_csv(out/"phase2_topK_membership.csv",index=False)

    corrs=[]
    for s1,s2 in combinations("ABC",2):
        corrs.append({"comparison":s1+"_vs_"+s2,
                      "Spearman_rho_descriptive":d[f"rank_{s1}"].corr(d[f"rank_{s2}"],method="spearman")})
    pd.DataFrame(corrs).to_csv(out/"phase2_rank_correlation.csv",index=False)

    div=[]
    for s in "ABC":
        for k in ks:
            g=d[d[f"rank_{s}"]<=k]
            div.append({"strategy":s,"K":k,"CTNNB1":int((g.gene=="CTNNB1").sum()),
                "TP53":int((g.gene=="TP53").sum()),
                "DRB1_0101":int((g.MHC=="DRB1_0101").sum()),
                "DRB1_0401":int((g.MHC=="DRB1_0401").sum()),
                "unique_core_groups":g.core_group.nunique(),
                "core_redundancy_windows":len(g)-g.core_group.nunique(),
                "positive_delta_BA":int(g.improved_BA.sum()),
                "BA_rank_le_5":int((g.Rank_BA_MUT<=5).sum()),
                "EL_rank_le_5":int((g.Rank_EL_MUT<=5).sum()),
                "mutation_in_core":int(g.Mutation_in_core_MUT.sum()),
                "core_offset_changed":int(g.Predicted_core_offset_changed.sum())})
    diversity=pd.DataFrame(div)
    diversity.to_csv(out/"phase2_topK_diversity_and_caveats.csv",index=False)

    # A core collapse is a sensitivity audit, not another independent validation.
    collapsed=[]
    for s in "ABC":
        g=d.sort_values(f"rank_{s}").drop_duplicates("core_group")
        collapsed.append({"strategy":s,"unique_core_groups_in_pool":len(g),
            "top10_distinct_mutant_cores":d[d[f"rank_{s}"]<=10].core_group.nunique(),
            "top10_core_collapsed_representatives":min(10,len(g)),
            "first10_distinct_core_ids":";".join(g.head(10).candidate_id)})
    pd.DataFrame(collapsed).to_csv(out/"phase2_core_collapse_sensitivity.csv",index=False)
    save_plots(d,ranks,out,ks)

    top10={s:set(ranks[s][:10]) for s in "ABC"}
    lines=[
        "PHASE 2 DECISION ANALYSIS PASSED",
        "Status: post-pilot method lock; descriptive, NOT independently validated or prospectively blinded",
        "Input: "+str(args.input),
        "Input SHA-256: "+current_hash,
        "Rules SHA-256: "+digest(args.rules),
        f"Candidates: {len(d)} (2 mutations; 2 HLA; overlapping windows)",
        "Primary Top-K: "+str(ks),
        "No experimental CD4/binding outcomes were used.",
        "NOTE: Strategies optimise different goals; overlap and correlations do NOT establish biological superiority.",
    ]
    for s in "ABC":
        g=d[d[f"rank_{s}"]<=10]
        lines.append(f"Top10 {s}: {len(g)} candidates, {g.core_group.nunique()} unique (gene,HLA,MUT core) groups; "
                     f"{int(g.improved_BA.sum())} BA improvements; {int((g.Rank_BA_MUT<=5).sum())} BA ranks <=5")
        for r in g.sort_values(f"rank_{s}").head(3).itertuples():
            lines.append(f"  rank {getattr(r,'rank_'+s)} {r.candidate_id} "
                         f"MUT_BA={r.Rank_BA_MUT:.2f}, MUT_EL={r.Rank_EL_MUT:.2f}, dBA={r.deltaBA:.2f}, "
                         f"question={r.C_selection_question if s=='C' else s}")
    for s1,s2 in combinations("ABC",2):
        inter=len(top10[s1]&top10[s2]);union=len(top10[s1]|top10[s2]);
        lines.append(f"Top10 overlap {s1}-{s2}: {inter} shared; Jaccard={inter/union:.3f}")
    lines.append("All output generated solely from the supplied Phase 1 prediction comparisons.")
    (out/"phase2_QC_and_results.txt").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines))


if __name__=="__main__":
    run(parse_args())
