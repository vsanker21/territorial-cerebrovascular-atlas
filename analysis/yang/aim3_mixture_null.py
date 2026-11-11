"""
H2 / Aim 3 — convex mixture null (hardened).

Primary (within-Allen):
  Borderzone = MCA_leaning Endo; flanks = ACA + MCA.
  Features = genes with strongest |ACA−MCA| difference (flank-separating).
  Fit border ≈ w·ACA+(1−w)·MCA; residual pathway modules vs gene-permutation null.

Aim3 exploratory:
  Yang control SFG vascular cells ≈ Allen ACC Endo + Garcia MCA Endo (cross-dataset).

LOCKED interpretation rules in report.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
ALLEN = ROOT / "data/processed/allen/allen_endo_multiregion.h5ad"
YANG = ROOT / "data/processed/yang/yang_vascular_qc_tight.h5ad"
GARCIA = ROOT / "data/processed/garcia/garcia_temporal_mca_qc_tight.h5ad"
OUT = ROOT / "results/yang"

MODULES = {
    "BBB_tight_junction": ["CLDN5", "OCLN", "TJP1", "TJP2", "F11R", "ESAM"],
    "angiogenesis": ["VEGFA", "FLT1", "KDR", "PECAM1", "ENG", "VWF", "CDH5"],
    "ECM_basement": ["COL4A1", "COL4A2", "LAMA4", "LAMB1", "HSPG2", "NID1"],
    "hypoxia_stress": ["HIF1A", "VEGFA", "BNIP3", "SLC2A1", "DDIT3", "ATF4"],
    "inflammation": ["ICAM1", "VCAM1", "SELE", "CCL2", "IL6", "NFKBIA"],
    "autoregulation_proxy": ["KCNJ8", "ABCC9", "NOTCH3", "PDGFRB", "ACTA2", "TAGLN"],
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X, dtype=float)


def donor_means(a: ad.AnnData, donor_col: str, genes: list[str], min_cells=3) -> pd.DataFrame:
    present = [g for g in genes if g in a.var_names]
    X = dense(a[:, present].X)
    donors = a.obs[donor_col].astype(str).to_numpy()
    rows = []
    for d in pd.unique(donors):
        m = donors == d
        if int(m.sum()) < min_cells:
            continue
        rows.append({"donor": d, "n": int(m.sum()), **{g: float(X[m, i].mean()) for i, g in enumerate(present)}})
    return pd.DataFrame(rows)


def select_flank_genes(aca: np.ndarray, mca: np.ndarray, genes: list[str], k=500) -> list[str]:
    delta = np.abs(aca - mca)
    order = np.argsort(delta)[::-1]
    keep = [genes[i] for i in order[:k] if np.isfinite(delta[i]) and delta[i] > 0]
    return keep


def fit_mixture(border: np.ndarray, aca: np.ndarray, mca: np.ndarray) -> dict:
    ws = np.linspace(0, 1, 101)
    best = None
    b2 = float(np.sum(border**2)) + 1e-12
    for w in ws:
        pred = w * aca + (1.0 - w) * mca
        sse = float(np.sum((border - pred) ** 2))
        if best is None or sse < best["sse"]:
            cos = float(np.dot(border, pred) / (np.linalg.norm(border) * np.linalg.norm(pred) + 1e-12))
            best = {
                "w_ACA": float(w),
                "w_MCA": float(1.0 - w),
                "sse": sse,
                "cosine": cos,
                "pred": pred,
                "frac_var_explained": float(1.0 - sse / b2),
            }
    best["residual"] = border - best["pred"]
    return best


def module_mean(vec: np.ndarray, genes: list[str], module_genes: list[str]) -> float:
    idx = [i for i, g in enumerate(genes) if g in set(module_genes)]
    if len(idx) < 2:
        return np.nan
    return float(np.mean(vec[idx]))


def analyze(border_pb, aca_pb, mca_pb, genes_all, label: str, k=500, n_perm=300) -> dict:
    genes = [g for g in genes_all if g in border_pb.columns and g in aca_pb.columns and g in mca_pb.columns]
    if len(border_pb) < 2 or len(aca_pb) < 2 or len(mca_pb) < 2 or len(genes) < 200:
        return {
            "contrast": label,
            "status": "insufficient",
            "n_genes": len(genes),
            "n_border_donors": int(len(border_pb)),
            "n_aca_donors": int(len(aca_pb)),
            "n_mca_donors": int(len(mca_pb)),
        }
    aca_v = aca_pb[genes].mean(axis=0).to_numpy(float)
    mca_v = mca_pb[genes].mean(axis=0).to_numpy(float)
    border_v = border_pb[genes].mean(axis=0).to_numpy(float)
    sel = select_flank_genes(aca_v, mca_v, genes, k=k)
    ix = [genes.index(g) for g in sel]
    aca_s = aca_v[ix]
    mca_s = mca_v[ix]
    border_s = border_v[ix]
    # standardize features by flank pooled std
    sd = np.std(np.vstack([aca_s, mca_s]), axis=0) + 1e-6
    aca_z = (aca_s - aca_s.mean()) / sd
    mca_z = (mca_s - mca_s.mean()) / sd
    border_z = (border_s - border_s.mean()) / sd

    fit = fit_mixture(border_z, aca_z, mca_z)
    # Module residuals on full shared-gene space with locked mixture weights
    # (selected genes may omit sparse pathway genes).
    w = fit["w_ACA"]
    full_pred = w * aca_v + (1.0 - w) * mca_v
    full_resid = border_v - full_pred
    obs_mod = {name: module_mean(full_resid, genes, mg) for name, mg in MODULES.items()}

    rng = np.random.default_rng(0)
    null_fve = []
    null_mod = {name: [] for name in MODULES}
    for _ in range(n_perm):
        b = border_z.copy()
        rng.shuffle(b)
        f = fit_mixture(b, aca_z, mca_z)
        null_fve.append(f["frac_var_explained"])
        # Permute residual gene labels for module null (keeps mixture structure)
        resid_perm = full_resid.copy()
        rng.shuffle(resid_perm)
        for name, mg in MODULES.items():
            null_mod[name].append(module_mean(resid_perm, genes, mg))

    mod_stats = {}
    for name in MODULES:
        o = obs_mod[name]
        null = np.asarray(null_mod[name], float)
        null = null[np.isfinite(null)]
        if not np.isfinite(o) or len(null) < 50:
            mod_stats[name] = {"obs_residual_module": o, "p_two_sided": np.nan}
            continue
        p = float((np.abs(null - null.mean()) >= abs(o - null.mean())).mean())
        mod_stats[name] = {
            "obs_residual_module": float(o),
            "null_mean": float(null.mean()),
            "p_two_sided": p,
        }

    return {
        "contrast": label,
        "status": "ok",
        "n_genes_selected": len(sel),
        "n_border_donors": int(len(border_pb)),
        "n_aca_donors": int(len(aca_pb)),
        "n_mca_donors": int(len(mca_pb)),
        "w_ACA": fit["w_ACA"],
        "w_MCA": fit["w_MCA"],
        "frac_var_explained": fit["frac_var_explained"],
        "cosine": fit["cosine"],
        "p_fve_ge_obs": float((np.asarray(null_fve) >= fit["frac_var_explained"]).mean()),
        "residual_modules": mod_stats,
    }


def verdict_from(primary: dict) -> str:
    if primary.get("status") != "ok":
        return "INCONCLUSIVE_insufficient_data"
    sig = [
        k
        for k, v in primary.get("residual_modules", {}).items()
        if isinstance(v, dict) and np.isfinite(v.get("p_two_sided", np.nan)) and v["p_two_sided"] < 0.05
    ]
    fve = primary["frac_var_explained"]
    # High FVE: flanks explain most variance. Significant residual modules: structured leftover.
    if fve >= 0.7 and len(sig) == 0:
        return "MIXTURE_NULL_NOT_REJECTED_within_Allen"
    if fve >= 0.7 and len(sig) >= 1:
        return "MOSTLY_MIXTURE_WITH_STRUCTURED_RESIDUAL_exploratory"
    if fve < 0.7 and len(sig) >= 1:
        return "STRUCTURED_RESIDUAL_H2_COMPATIBLE_exploratory"
    return "INCONCLUSIVE_mixed_signals"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log("Allen primary mixture…")
    allen = ad.read_h5ad(ALLEN)
    genes = list(allen.var_names)
    aca = allen[allen.obs["tv_arm"].astype(str) == "ACA"]
    mca = allen[allen.obs["tv_arm"].astype(str) == "MCA"]
    bz = allen[allen.obs["tv_arm"].astype(str) == "MCA_leaning"]
    primary = analyze(
        donor_means(bz, "donor_id", genes),
        donor_means(aca, "donor_id", genes),
        donor_means(mca, "donor_id", genes),
        genes,
        "within_Allen_MCA_leaning_vs_ACA_MCA",
    )

    log("Yang Aim3 exploratory…")
    yang = ad.read_h5ad(YANG)
    y = yang[
        (yang.obs["disease"].astype(str).str.lower() == "control")
        & (yang.obs["region"].astype(str) == "superior_frontal_cortex")
    ].copy()
    # Prefer endo; else all vascular (Aim3 stream often uses vascular)
    if "vascular_lineage" in y.obs.columns:
        ye = y[y.obs["vascular_lineage"].astype(str).str.lower().str.contains("endo")]
        if ye.n_obs >= 50:
            y = ye
    garcia = ad.read_h5ad(GARCIA)
    if "SampleID" in garcia.obs.columns:
        garcia.obs["donor_id"] = garcia.obs["SampleID"].astype(str)
    else:
        garcia.obs["donor_id"] = "garcia_" + garcia.obs.index.astype(str)
    # Use symbol-aligned Allen Endo ACA (not ENSG-indexed allen_acc CELLxGENE export)
    allen_aca = allen[allen.obs["tv_arm"].astype(str) == "ACA"].copy()
    genes2 = sorted(set(y.var_names) & set(allen_aca.var_names) & set(garcia.var_names))
    exploratory = analyze(
        donor_means(y, "donor_id", genes2, min_cells=5),
        donor_means(allen_aca, "donor_id", genes2, min_cells=3),
        donor_means(garcia, "donor_id", genes2, min_cells=3),
        genes2,
        "Yang_SFG_ctrl_vs_AllenACA_GarciaMCA",
        k=400,
    )

    verdict = verdict_from(primary)
    summary = {
        "hypothesis": "H2_convex_mixture_null",
        "verdict_primary": verdict,
        "primary_within_allen": primary,
        "aim3_yang_exploratory": exploratory,
        "interpretation_locks": [
            "Primary H2 is within-Allen MCA_leaning as mixture of ACA+MCA (same technology).",
            "Yang SFG vs Allen ACA + Garcia MCA is exploratory (batch + enrichment asymmetry).",
            "Yang SFG is Aim3 borderzone, not confirmatory Class B ACA.",
            "MIXTURE_NULL_NOT_REJECTED => borderzone largely explained by flanks.",
            "STRUCTURED_RESIDUAL_H2_COMPATIBLE => leftover module signal after mixture (exploratory).",
        ],
    }
    (OUT / "aim3_mixture_null_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    def mods(d):
        rm = d.get("residual_modules") or {}
        return "\n".join(
            f"  - {k}: resid={v.get('obs_residual_module')}, p={v.get('p_two_sided')}" for k, v in rm.items()
        )

    md = [
        "# H2 / Aim 3 — convex mixture null",
        "",
        f"**Primary verdict:** `{verdict}`",
        "",
        "## Primary — within-Allen (flank-separating genes)",
        "Borderzone: MCA_leaning Endo. Flanks: ACA (ACC) + MCA (MTG+ANG).",
        f"- status={primary.get('status')}",
        f"- w_ACA={primary.get('w_ACA')}, w_MCA={primary.get('w_MCA')}",
        f"- frac_var_explained={primary.get('frac_var_explained')}, cosine={primary.get('cosine')}",
        f"- p(FVE ≥ obs | gene permutation)={primary.get('p_fve_ge_obs')}",
        f"- donors border/ACA/MCA={primary.get('n_border_donors')}/{primary.get('n_aca_donors')}/{primary.get('n_mca_donors')}",
        f"- n_genes_selected={primary.get('n_genes_selected')}",
        "Residual modules:",
        mods(primary),
        "",
        "## Exploratory — Yang SFG (Aim3) vs Allen ACA + Garcia MCA",
        "**Caveats:** cross-dataset; enrichment asymmetry; not confirmatory H1/H2.",
        f"- status={exploratory.get('status')}",
        f"- w_ACA={exploratory.get('w_ACA')}, w_MCA={exploratory.get('w_MCA')}",
        f"- frac_var_explained={exploratory.get('frac_var_explained')}, cosine={exploratory.get('cosine')}",
        f"- p(FVE ≥ obs)={exploratory.get('p_fve_ge_obs')}",
        f"- donors border/ACA/MCA={exploratory.get('n_border_donors')}/{exploratory.get('n_aca_donors')}/{exploratory.get('n_mca_donors')}",
        "Residual modules:",
        mods(exploratory),
        "",
        "## Locked interpretation",
        "- Primary within-Allen: flanks explain most borderzone variance (high FVE), with residual pathway modules that can remain non-null after mixture.",
        "- Verdict `MOSTLY_MIXTURE_WITH_STRUCTURED_RESIDUAL_exploratory` = mixture-dominant + structured leftover (exploratory H2-compatible residual).",
        "- Verdict `STRUCTURED_RESIDUAL_H2_COMPATIBLE_exploratory` = low FVE + coherent residual (stronger H2-compatible).",
        "- Do not call Yang SFG clean ACA; Yang arm is exploratory (cross-dataset / enrichment asymmetry).",
        "",
    ]
    (OUT / "aim3_mixture_null_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps({"verdict_primary": verdict, "primary": {k: primary.get(k) for k in ["w_ACA", "w_MCA", "frac_var_explained", "p_fve_ge_obs", "status"]}, "yang": {k: exploratory.get(k) for k in ["status", "w_ACA", "frac_var_explained"]}}, indent=2))


if __name__ == "__main__":
    main()
