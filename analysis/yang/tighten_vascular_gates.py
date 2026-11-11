"""
Tighten Yang endo/peri gates on top of soft vascular QC, then re-run Aim3 + H4.

Gate v2 (marker-aware):
  - Require canonical lineage markers (expression > 0)
  - Require module-score margins vs other vascular lineages + parenchyma
  - Exclude neuronal/astro/oligo/micro contaminant markers (incl. LINGO1)
  - Drop ambiguous cells (not clean endo or peri for Aim3/H4)

Writes:
  data/processed/yang/yang_vascular_qc_tight.h5ad
  results/yang/gate_v2_*.{csv,md}
  results/yang/aim3_lineage_stratified_de_tight.csv
  results/h4/yang_territory_classifier_h4_tight.{csv,md}
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis"))
from continue_existing_data import run_aim3_lineage, run_h4_classifier  # noqa: E402

IN_H5AD = ROOT / "data/processed/yang/yang_vascular_qc.h5ad"
OUT_H5AD = ROOT / "data/processed/yang/yang_vascular_qc_tight.h5ad"
OUT_GATE = ROOT / "results/yang/gate_v2_summary.csv"
OUT_MD = ROOT / "results/yang/gate_v2_report.md"
OUT_AIM3 = ROOT / "results/yang/aim3_lineage_stratified_de_tight.csv"
OUT_H4 = ROOT / "results/h4/yang_territory_classifier_h4_tight.csv"
OUT_H4_MD = ROOT / "results/h4/yang_territory_classifier_h4_tight_report.md"
OUT_COMPARE = ROOT / "results/yang/aim3_h4_soft_vs_tight.md"

# Canonical positive markers (any one sufficient)
ENDO_MARKERS = ["CLDN5", "PECAM1", "FLT1", "VWF", "CDH5"]
PERI_MARKERS = ["PDGFRB", "RGS5", "KCNJ8", "ABCC9"]

# Contaminant genes excluded when above threshold (log1p-normalized)
CONTAMINANT_MAX = {
    "LINGO1": 0.5,
    "SNAP25": 0.5,
    "SLC17A7": 0.5,
    "RBFOX3": 0.5,
    "GAD1": 0.5,
    "MBP": 1.0,
    "MOBP": 1.0,
    "PLP1": 1.0,
    "AQP4": 1.0,
    "GFAP": 1.5,
    "PTPRC": 0.5,
}

# Module-score thresholds
ENDO_MIN_SCORE = 0.5
PERI_MIN_SCORE = 0.25
LINEAGE_MARGIN = 0.40  # vs next-best vascular lineage
PARENCHYMA_MAX = 0.0
NEURON_MAX = 0.05
ASTRO_MAX = 0.25
OLIGO_MAX = 0.0
MICRO_MAX = 0.05


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def expr_matrix(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return pd.DataFrame(index=adata.obs_names)
    X = adata[:, present].X
    if sparse.issparse(X):
        X = X.toarray()
    return pd.DataFrame(X, index=adata.obs_names, columns=present)


def apply_tight_gates(adata: ad.AnnData) -> ad.AnnData:
    obs = adata.obs
    all_genes = list(CONTAMINANT_MAX) + ENDO_MARKERS + PERI_MARKERS
    expr = expr_matrix(adata, all_genes)

    # Contaminant mask (True = contaminated)
    contam = pd.Series(False, index=obs.index)
    for g, thr in CONTAMINANT_MAX.items():
        if g in expr.columns:
            contam |= expr[g] > thr

    endo_pos = pd.Series(False, index=obs.index)
    for g in ENDO_MARKERS:
        if g in expr.columns:
            endo_pos |= expr[g] > 0

    peri_pos = pd.Series(False, index=obs.index)
    for g in PERI_MARKERS:
        if g in expr.columns:
            peri_pos |= expr[g] > 0

    vasc_cols = ["score_endo", "score_peri", "score_smc", "score_fib"]
    endo_margin = obs["score_endo"] - obs[["score_peri", "score_smc", "score_fib"]].max(axis=1)
    peri_margin = obs["score_peri"] - obs[["score_endo", "score_smc", "score_fib"]].max(axis=1)

    clean_shared = (
        (~contam)
        & (obs["score_parenchyma"] <= PARENCHYMA_MAX)
        & (obs["score_neuron"] <= NEURON_MAX)
        & (obs["score_astro"] <= ASTRO_MAX)
        & (obs["score_oligo"] <= OLIGO_MAX)
        & (obs["score_micro"] <= MICRO_MAX)
    )

    is_endo = (
        clean_shared
        & endo_pos
        & (obs["score_endo"] >= ENDO_MIN_SCORE)
        & (endo_margin >= LINEAGE_MARGIN)
    )
    is_peri = (
        clean_shared
        & peri_pos
        & (obs["score_peri"] >= PERI_MIN_SCORE)
        & (peri_margin >= LINEAGE_MARGIN)
        & (~is_endo)  # prefer endo if both pass (rare)
    )

    lineage = pd.Series("ambiguous", index=obs.index, dtype=object)
    lineage[is_endo] = "endo"
    lineage[is_peri] = "peri"
    # Keep smc/fib from soft gate only if they also pass clean_shared + their score dominates
    soft = obs["vascular_lineage"].astype(str)
    for lin, col, amin in [("smc", "score_smc", 0.25), ("fib", "score_fib", 0.25)]:
        marg = obs[col] - obs[[c for c in vasc_cols if c != col]].max(axis=1)
        keep = clean_shared & (soft == lin) & (obs[col] >= amin) & (marg >= LINEAGE_MARGIN)
        lineage[keep & (lineage == "ambiguous")] = lin

    out = adata.copy()
    out.obs["vascular_lineage_soft"] = soft.astype("category")
    out.obs["vascular_lineage"] = lineage.astype("category")
    out.obs["pass_gate_v2"] = (lineage != "ambiguous").astype(bool)
    out.obs["contam_gate_v2"] = contam.astype(bool)
    out.obs["gate_version"] = "v2_marker_aware"

    # Restrict to cells that pass as any vascular lineage (drop ambiguous)
    out = out[out.obs["pass_gate_v2"]].copy()
    return out


def lingo1_rank(de: pd.DataFrame, lineage: str, disease: str) -> str:
    sub = de[(de["lineage"] == lineage) & (de["disease"] == disease)]
    if sub.empty or "names" not in sub.columns:
        return "n/a"
    hit = sub.reset_index(drop=True)
    idxs = hit.index[hit["names"] == "LINGO1"].tolist()
    if not idxs:
        return "not in top list"
    i = idxs[0]
    row = hit.iloc[i]
    return f"rank {i+1}, logFC={row['logfoldchanges']:.2f}, padj={row['pvals_adj']:.2e}"


def main():
    log(f"Loading {IN_H5AD}")
    soft = sc.read_h5ad(IN_H5AD)
    n_soft = soft.n_obs
    soft_counts = soft.obs["vascular_lineage"].astype(str).value_counts().to_dict()
    log(f"Soft-gated cells: {n_soft:,} | {soft_counts}")

    tight = apply_tight_gates(soft)
    log(f"After tight gate: {tight.n_obs:,}")
    tight_counts = tight.obs["vascular_lineage"].astype(str).value_counts().to_dict()
    log(f"Tight lineages: {tight_counts}")

    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    tight.write_h5ad(OUT_H5AD, compression="gzip")
    log(f"Wrote {OUT_H5AD}")

    # Gate funnel summary
    rows = []
    for stage, obj, col in [
        ("soft", soft, "vascular_lineage"),
        ("tight", tight, "vascular_lineage"),
    ]:
        g = (
            obj.obs.groupby([col, "territory_contrast", "disease"], observed=True)
            .size()
            .reset_index(name="n_cells")
        )
        g = g.rename(columns={col: "vascular_lineage"})
        g["gate"] = stage
        rows.append(g)
    gate_df = pd.concat(rows, ignore_index=True)
    OUT_GATE.parent.mkdir(parents=True, exist_ok=True)
    gate_df.to_csv(OUT_GATE, index=False)

    # Aim3 + H4 on tight object
    log("Aim3 lineage DE (tight)")
    aim3 = run_aim3_lineage(tight)
    aim3.to_csv(OUT_AIM3, index=False)
    log(f"Wrote {OUT_AIM3} ({len(aim3)} rows)")

    log("H4 classifier (tight)")
    h4_df, h4_sum = run_h4_classifier(tight)
    h4_df.to_csv(OUT_H4, index=False)
    h4_md = [
        "# H4 first-pass — Yang endothelial territory classifier (gate v2 tight)",
        "",
        "Trained on **control** endo (tight): PCA hippocampus vs ACA/MCA borderzone SFG.",
        "Transferred to **AD** endo without retraining.",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for k, v in h4_sum.items():
        if isinstance(v, float):
            h4_md.append(f"| {k} | {v:.4f} |")
        else:
            h4_md.append(f"| {k} | {v} |")
    h4_md += [
        "",
        f"**Delta AUROC (AD - control CV) = {h4_sum['delta_auroc_ad_minus_control']:.4f}**",
        "",
        "Negative delta is consistent with H4 (degradation); still exploratory until donor-aware CV.",
    ]
    OUT_H4_MD.write_text("\n".join(h4_md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_H4}")

    # Soft Aim3 for LINGO1 comparison if available
    soft_aim3_path = ROOT / "results/yang/aim3_lineage_stratified_de.csv"
    soft_h4_path = ROOT / "results/h4/yang_territory_classifier_h4.csv"
    soft_aim3 = pd.read_csv(soft_aim3_path) if soft_aim3_path.exists() else pd.DataFrame()
    soft_h4 = pd.read_csv(soft_h4_path) if soft_h4_path.exists() else pd.DataFrame()

    lines = [
        "# Yang gate v2 (tight) report",
        "",
        "## Gate rules",
        "",
        f"- Endo: any of {ENDO_MARKERS} > 0; score_endo ≥ {ENDO_MIN_SCORE}; "
        f"margin vs other vascular ≥ {LINEAGE_MARGIN}",
        f"- Peri: any of {PERI_MARKERS} > 0; score_peri ≥ {PERI_MIN_SCORE}; "
        f"margin ≥ {LINEAGE_MARGIN}",
        f"- Shared: parenchyma ≤ {PARENCHYMA_MAX}; neuron ≤ {NEURON_MAX}; "
        f"astro ≤ {ASTRO_MAX}; oligo ≤ {OLIGO_MAX}; micro ≤ {MICRO_MAX}",
        f"- Contaminant expression caps: {CONTAMINANT_MAX}",
        "",
        "## Cell counts",
        "",
        f"- Soft vascular gate: **{n_soft:,}** -> `{soft_counts}`",
        f"- Tight gate v2: **{tight.n_obs:,}** -> `{tight_counts}`",
        f"- Retention: **{100.0 * tight.n_obs / n_soft:.1f}%** of soft-gated cells",
        f"- Object: `{OUT_H5AD.relative_to(ROOT).as_posix()}`",
        "",
        "## Lineage x contrast x disease (tight)",
        "",
        gate_df[gate_df["gate"] == "tight"]
        .pivot_table(
            index=["vascular_lineage", "territory_contrast"],
            columns="disease",
            values="n_cells",
            aggfunc="sum",
            fill_value=0,
        )
        .to_string(),
        "",
        "## Aim3 top hits (tight endo/peri)",
        "",
    ]
    for lineage in ["endo", "peri"]:
        for disease in ["control", "AD"]:
            sub = aim3[(aim3["lineage"] == lineage) & (aim3["disease"] == disease)].head(10)
            lines += [
                f"### {lineage} / {disease}",
                "",
                sub[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False),
                "",
            ]

    lines += [
        "## H4 (tight)",
        "",
        f"- Control CV AUROC: **{h4_sum['control_cv_auroc']:.3f}**",
        f"- AD transfer AUROC: **{h4_sum['ad_transfer_auroc']:.3f}**",
        f"- Delta: **{h4_sum['delta_auroc_ad_minus_control']:.3f}**",
        "",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")

    # Soft vs tight comparison
    cmp_lines = [
        "# Soft vs tight gate — Aim3 / H4 comparison",
        "",
        "## Retention",
        f"- Soft: {n_soft:,} -> Tight: {tight.n_obs:,} ({100.0 * tight.n_obs / n_soft:.1f}%)",
        f"- Soft endo/peri: {soft_counts.get('endo', 0):,} / {soft_counts.get('peri', 0):,}",
        f"- Tight endo/peri: {tight_counts.get('endo', 0):,} / {tight_counts.get('peri', 0):,}",
        "",
        "## LINGO1 rank in Aim3 DE (lower/absent = cleaner)",
        "",
        "| Stratum | Soft gate | Tight gate |",
        "|---|---|---|",
    ]
    for lineage in ["endo", "peri"]:
        for disease in ["control", "AD"]:
            s = lingo1_rank(soft_aim3, lineage, disease) if len(soft_aim3) else "n/a"
            t = lingo1_rank(aim3, lineage, disease)
            cmp_lines.append(f"| {lineage}/{disease} | {s} | {t} |")

    cmp_lines += ["", "## H4 AUROC", ""]
    if len(soft_h4):
        srow = soft_h4.iloc[0]
        cmp_lines += [
            "| Gate | Control CV | AD transfer | Delta |",
            "|---|---:|---:|---:|",
            f"| Soft | {srow['control_cv_auroc']:.3f} | {srow['ad_transfer_auroc']:.3f} | "
            f"{srow['delta_auroc_ad_minus_control']:.3f} |",
            f"| Tight | {h4_sum['control_cv_auroc']:.3f} | {h4_sum['ad_transfer_auroc']:.3f} | "
            f"{h4_sum['delta_auroc_ad_minus_control']:.3f} |",
        ]
    else:
        cmp_lines.append(
            f"Tight only: control={h4_sum['control_cv_auroc']:.3f}, "
            f"AD={h4_sum['ad_transfer_auroc']:.3f}, "
            f"delta={h4_sum['delta_auroc_ad_minus_control']:.3f}"
        )

    OUT_COMPARE.write_text("\n".join(cmp_lines) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_COMPARE}")
    try:
        print("\n".join(cmp_lines))
    except UnicodeEncodeError:
        print(OUT_COMPARE.read_text(encoding="utf-8").encode("ascii", "replace").decode("ascii"))


if __name__ == "__main__":
    main()
