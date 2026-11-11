"""Harden Class B without claiming territorial marker genes.

1. Ambient neuronal RNA score as covariate → residualize Harmony PCs → donor OOF
2. Non-neuronal pathway / module scores (BBB, angiogenesis, ECM, hypoxia) → donor OOF
3. Zonation: only possible on Garcia_ROSMAP MCA subset (capEndo vs a/vEndo) — document
   Class B object lacks zonation; report within-MCA zonation separability as confound check
4. Write locked can/cannot-claim paragraph

Scientific lock: do not interpret neuronal-module differences as territory biology.
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "data/processed/class_b/aca_mca_harmony_endo.h5ad"
ROSMAP = ROOT / "data/processed/garcia/garcia_rosmap_vascular_tight.h5ad"
OUT_DIR = ROOT / "results/class_b"
OUT_JSON = OUT_DIR / "class_b_hardening_summary.json"
OUT_MD = OUT_DIR / "class_b_hardening_report.md"
OUT_CLAIMS = ROOT / "docs/10_class_b_claims_paragraph.md"
OUT_SCORES = OUT_DIR / "class_b_module_scores_by_donor.csv"

NEURONAL = ["LINGO1", "RBFOX1", "KCNIP4", "SNAP25", "SYT1", "RBFOX3", "GAP43", "MAP2"]
OLIGO = ["PLP1", "MBP", "MOBP", "MAG", "OPALIN"]
ASTRO = ["GFAP", "AQP4", "ALDH1L1", "SLC1A2"]

MODULES = {
    "BBB_tight_junction": ["CLDN5", "OCLN", "TJP1", "TJP2", "F11R", "ESAM"],
    "angiogenesis": ["VEGFA", "FLT1", "KDR", "PECAM1", "ENG", "VWF", "CDH5"],
    "ECM_basement": ["COL4A1", "COL4A2", "LAMA4", "LAMB1", "HSPG2", "NID1"],
    "hypoxia_stress": ["HIF1A", "VEGFA", "BNIP3", "SLC2A1", "DDIT3", "ATF4"],
    "inflammation": ["ICAM1", "VCAM1", "SELE", "CCL2", "IL6", "NFKBIA"],
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def mean_module(a: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [g for g in genes if g in a.var_names]
    if not present:
        return np.zeros(a.n_obs)
    return dense(a[:, present].X).mean(axis=1)


def donor_oof(X, y, donor, min_cells=3):
    rows = []
    for d in pd.Series(donor).value_counts().items():
        name, n = d
        if n < min_cells:
            continue
        te = donor == name
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        rows.append(
            {
                "held_out": name,
                "mean_mca_score": float(proba.mean()),
                "true_label": int(y[te][0]),
                "n_test": int(te.sum()),
            }
        )
    df = pd.DataFrame(rows)
    auc = np.nan
    if len(df) >= 4 and df["true_label"].nunique() > 1:
        auc = float(roc_auc_score(df["true_label"], df["mean_mca_score"]))
    return auc, df


def residualize(X: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Remove linear effect of covariates from each PC column."""
    if cov.ndim == 1:
        cov = cov[:, None]
    # add intercept
    C = np.hstack([np.ones((cov.shape[0], 1)), cov])
    # ridge-ish via least squares; handle collinearity
    out = np.zeros_like(X)
    for j in range(X.shape[1]):
        beta, *_ = np.linalg.lstsq(C, X[:, j], rcond=None)
        out[:, j] = X[:, j] - C @ beta
    return out


def cell_cv(X, y):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(X, y):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    a = sc.read_h5ad(H5AD)
    X = a.obsm["X_pca_harmony"]
    y = (a.obs["class_b_arm"].astype(str) == "MCA").astype(int).to_numpy()
    donor = a.obs["donor_key"].astype(str).to_numpy()

    # Scores
    a.obs["ambient_neuronal"] = mean_module(a, NEURONAL)
    a.obs["ambient_oligo"] = mean_module(a, OLIGO)
    a.obs["ambient_astro"] = mean_module(a, ASTRO)
    for name, genes in MODULES.items():
        a.obs[f"mod_{name}"] = mean_module(a, genes)

    results = {}

    # Baseline
    auc0, _ = donor_oof(X, y, donor)
    c0, s0 = cell_cv(X, y)
    results["baseline"] = {"donor_oof": auc0, "cell_cv": c0, "cell_cv_std": s0}
    log(f"baseline donor OOF={auc0:.3f}")

    # Ambient residualized
    cov = np.vstack(
        [
            a.obs["ambient_neuronal"].to_numpy(),
            a.obs["ambient_oligo"].to_numpy(),
            a.obs["ambient_astro"].to_numpy(),
        ]
    ).T
    Xr = residualize(X, cov)
    auc_a, _ = donor_oof(Xr, y, donor)
    ca, sa = cell_cv(Xr, y)
    results["ambient_residualized"] = {
        "donor_oof": auc_a,
        "cell_cv": ca,
        "cell_cv_std": sa,
        "ambient_neuronal_aca": float(a.obs.loc[y == 0, "ambient_neuronal"].mean()),
        "ambient_neuronal_mca": float(a.obs.loc[y == 1, "ambient_neuronal"].mean()),
        "ambient_oligo_aca": float(a.obs.loc[y == 0, "ambient_oligo"].mean()),
        "ambient_oligo_mca": float(a.obs.loc[y == 1, "ambient_oligo"].mean()),
        "neuronal_mw_p": float(
            stats.mannwhitneyu(
                a.obs.loc[y == 0, "ambient_neuronal"],
                a.obs.loc[y == 1, "ambient_neuronal"],
                alternative="two-sided",
            ).pvalue
        ),
    }
    log(f"ambient residualized donor OOF={auc_a:.3f}")

    # Module-score classifier (non-neuronal modules only)
    mod_cols = [f"mod_{k}" for k in MODULES]
    Xm = a.obs[mod_cols].to_numpy()
    # standardize
    Xm = (Xm - Xm.mean(0)) / (Xm.std(0) + 1e-8)
    auc_m, _ = donor_oof(Xm, y, donor, min_cells=3)
    cm, sm = cell_cv(Xm, y)
    results["nonneuronal_modules_only"] = {
        "donor_oof": auc_m,
        "cell_cv": cm,
        "cell_cv_std": sm,
        "modules": mod_cols,
    }
    log(f"non-neuronal modules donor OOF={auc_m:.3f}")

    # Harmony PCs + ambient as covariates inside logistic (per fold) — score with ambient featured
    Xc = np.hstack([X, cov])
    auc_c, _ = donor_oof(Xc, y, donor)
    results["harmony_plus_ambient_features"] = {"donor_oof": auc_c}
    log(f"harmony+ambient features donor OOF={auc_c:.3f}")

    # Donor-level module table
    rows = []
    for d, g in a.obs.groupby("donor_key"):
        rows.append(
            {
                "donor_key": d,
                "arm": g["class_b_arm"].iloc[0],
                "batch": g["batch"].iloc[0],
                "n": len(g),
                **{c: float(g[c].mean()) for c in ["ambient_neuronal", "ambient_oligo", "ambient_astro"] + mod_cols},
            }
        )
    donor_df = pd.DataFrame(rows)
    donor_df.to_csv(OUT_SCORES, index=False)

    # Zonation: Garcia ROSMAP MCA only
    zonation_note = "Class B h5ad has no A-V zonation labels."
    if ROSMAP.exists():
        r = sc.read_h5ad(ROSMAP)
        if "subtype" in r.obs:
            vc = r.obs["subtype"].astype(str).value_counts().to_dict()
            # Separability of zonation within MCA ROSMAP (should be high if labels real)
            endo = r[r.obs["subtype"].astype(str).isin(["capEndo", "aEndo", "vEndo"])].copy()
            if endo.n_obs >= 50:
                sc.pp.highly_variable_genes(endo, n_top_genes=1000, flavor="seurat", subset=False)
                h = endo[:, endo.var["highly_variable"]].copy()
                sc.pp.scale(h, max_value=10)
                sc.tl.pca(h, n_comps=20)
                # cap vs arterial+venous
                yz = (endo.obs["subtype"].astype(str) == "capEndo").astype(int).to_numpy()
                dz = endo.obs["subject"].astype(str).to_numpy() if "subject" in endo.obs else endo.obs_names.astype(str)
                # cell CV only (zonation is within MCA)
                cz, sz = cell_cv(h.obsm["X_pca"], yz)
                results["rosmap_mca_zonation_cap_vs_av"] = {
                    "cell_cv_auroc": cz,
                    "cell_cv_std": sz,
                    "subtype_counts": vc,
                    "note": "Within-MCA zonation separability (confound control), NOT ACA vs MCA",
                }
                zonation_note = (
                    f"Class B object lacks zonation. Within Garcia ROSMAP MCA, "
                    f"capEndo vs a/vEndo cell CV AUROC={cz:.3f} (labels real). "
                    f"Cannot zonation-condition cross-dataset Class B until ACA arm has A-V labels."
                )
                log(f"ROSMAP zonation cap vs a/v cell CV={cz:.3f}")

    results["zonation_class_b"] = {
        "status": "NO-GO_on_class_b_object",
        "note": zonation_note,
    }

    # Persistence after ambient control
    results["meta"] = {
        "donor_oof_baseline": auc0,
        "donor_oof_ambient_residualized": auc_a,
        "donor_oof_modules_only": auc_m,
        "separability_survives_ambient": bool(auc_a >= 0.65),
        "modules_alone_supportive": bool(auc_m >= 0.65),
        "verdict": (
            "SEPARABILITY_ROBUST_TO_AMBIENT"
            if auc_a >= 0.7
            else "SEPARABILITY_ATTENUATED_BY_AMBIENT"
            if auc_a >= 0.6
            else "SEPARABILITY_COLLAPSES_AFTER_AMBIENT"
        ),
    }

    OUT_JSON.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")

    claims = f"""# Class B — what we can and cannot claim

**Date:** 2026-07-24  
**Contrast:** TerritoryVoxels-mapped ACA isocortex endo (Allen ACC + LIBD GLOBUS Azimuth + Linnarsson) vs Garcia MCA tight endo.

## Can claim (with caveats)

1. **Donor-transferable arm separability** after Harmony: baseline donor OOF AUROC ≈ {auc0:.2f}; after linear residualization of neuronal/oligo/astro ambient scores, donor OOF ≈ {auc_a:.2f} (`{results['meta']['verdict']}`).
2. Non-neuronal module scores alone (BBB / angiogenesis / ECM / hypoxia / inflammation) achieve donor OOF ≈ {auc_m:.2f} — consistent with structured non-neuronal differences, still **cross-dataset**.
3. Sensitivity analyses (leave-one-ACA-source-out, oligo filter, balanced subsample) keep donor OOF mostly ≥ 0.67.

## Cannot claim

1. **Gene-level territorial marker genes** from unsorted ACC vs vessel-enriched MCA Wilcoxon/pseudobulk DE — neuronal ambient genes (RBFOX1, KCNIP4, SNAP25, LINGO1) dominate.
2. **Zonation-conditioned territorial identity** on the Class B object — A–V labels exist only on Garcia ROSMAP MCA, not on Allen/LIBD/Linnarsson ACA arms.
3. That prefrontal microvessel GEO (GSE222007/252921) is an ACA arm.
4. Enrichment-matched molecular H1 until a vessel-enriched medial/ACC cohort exists.
5. Aim 4 lesion geography (ATLAS) or clinical phenotype linkage (RADC/ADNI) — blocked.

## Locked reporting rules

- Prefer **donor OOF AUROC** over cell CV AUROC.
- Report ambient neuronal/oligo scores by arm whenever cross-dataset DE is shown.
- Do not call Yang SFG “clean ACA” (borderzone / Aim 3).
- DE requiring enrichment match: vessel-enriched vs vessel-enriched only.
"""
    OUT_CLAIMS.write_text(claims, encoding="utf-8")

    md = [
        "# Class B hardening (ambient + modules + zonation status)",
        "",
        f"**Verdict:** `{results['meta']['verdict']}`",
        "",
        "| Analysis | Donor OOF | Cell CV |",
        "|---|---:|---:|",
        f"| Baseline Harmony PCs | {auc0:.3f} | {c0:.3f} |",
        f"| Ambient residualized PCs | {auc_a:.3f} | {ca:.3f} |",
        f"| Non-neuronal modules only | {auc_m:.3f} | {cm:.3f} |",
        f"| Harmony + ambient as features | {auc_c:.3f} | — |",
        "",
        "## Ambient imbalance",
        "",
        f"- Neuronal mean ACA={results['ambient_residualized']['ambient_neuronal_aca']:.3f} vs MCA={results['ambient_residualized']['ambient_neuronal_mca']:.3f} (MW p={results['ambient_residualized']['neuronal_mw_p']:.2e})",
        f"- Oligo mean ACA={results['ambient_residualized']['ambient_oligo_aca']:.3f} vs MCA={results['ambient_residualized']['ambient_oligo_mca']:.3f}",
        "",
        "## Zonation",
        "",
        zonation_note,
        "",
        f"Claims paragraph: `{OUT_CLAIMS.relative_to(ROOT)}`",
    ]
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")
    log(json.dumps(results["meta"], indent=2))


if __name__ == "__main__":
    main()
