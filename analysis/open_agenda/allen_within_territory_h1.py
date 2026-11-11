"""Within-atlas Allen Endo territory contrasts (confirmatory-style H1 path).

Same tech + published Endo taxonomy; map tissues → TerritoryVoxels priors.
Primary: ACC (ACA) vs MTG/ANG (MCA-leaning temporal) vs V1 (PCA-leaning) vs S1/A1.
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
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
NN = ROOT / "data/raw/aca_cohort_hunt/allen_cortical_acc/NonNeuronal.h5ad"
TV_ROSMAP = ROOT / "results/class_b/garcia_rosmap_territory_by_region.json"
OUT_DIR = ROOT / "results/h1"
OUT_H5AD = ROOT / "data/processed/allen/allen_endo_multiregion.h5ad"
OUT_JSON = OUT_DIR / "allen_within_territory_summary.json"
OUT_MD = OUT_DIR / "allen_within_territory_report.md"
OUT_DONOR = OUT_DIR / "allen_within_territory_donor_oof.csv"

# Tissue → operational territory label (locked rationale in report)
TISSUE_TV = {
    "anterior cingulate cortex": {"tv_arm": "ACA", "tv_p_ACA": 1.0, "tv_p_MCA": 0.0, "tv_p_PCA": 0.0, "assign": 0.72},
    "middle temporal gyrus": {"tv_arm": "MCA", "tv_p_ACA": 0.0, "tv_p_MCA": 1.0, "tv_p_PCA": 0.0, "assign": 0.55},
    "angular gyrus": {"tv_arm": "MCA", "tv_p_ACA": 0.0, "tv_p_MCA": 1.0, "tv_p_PCA": 0.0, "assign": 0.50},
    "primary visual cortex": {"tv_arm": "PCA", "tv_p_ACA": 0.0, "tv_p_MCA": 0.0, "tv_p_PCA": 1.0, "assign": 0.60},
    "primary somatosensory cortex": {"tv_arm": "MCA", "tv_p_ACA": 0.15, "tv_p_MCA": 0.85, "tv_p_PCA": 0.0, "assign": 0.40},
    "primary auditory cortex": {"tv_arm": "MCA", "tv_p_ACA": 0.0, "tv_p_MCA": 1.0, "tv_p_PCA": 0.0, "assign": 0.45},
    "dorsolateral prefrontal cortex": {"tv_arm": "MCA_leaning", "tv_p_ACA": 0.35, "tv_p_MCA": 0.65, "tv_p_PCA": 0.0, "assign": 0.30},
    "primary motor cortex": {"tv_arm": "MCA_leaning", "tv_p_ACA": 0.25, "tv_p_MCA": 0.75, "tv_p_PCA": 0.0, "assign": 0.35},
}

NEURONAL = ["LINGO1", "RBFOX1", "KCNIP4", "SNAP25"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def to_symbols(a: ad.AnnData) -> ad.AnnData:
    a = a.copy()
    if "feature_name" not in a.var:
        return a
    sym = a.var["feature_name"].astype(str).to_numpy()
    ok = np.array([s not in {"nan", "NA", ""} and not s.startswith("ENSG") for s in sym])
    a = a[:, ok].copy()
    a.var_names = pd.Index(a.var["feature_name"].astype(str))
    a.var_names_make_unique()
    return a


def donor_oof(X, y, donor, min_cells=3):
    rows = []
    for name, n in pd.Series(donor).value_counts().items():
        if n < min_cells:
            continue
        te = donor == name
        tr = ~te
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2 and True:
            if len(np.unique(y[tr])) < 2:
                continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        rows.append(
            {
                "held_out": name,
                "mean_pos_score": float(proba.mean()),
                "true_label": int(y[te][0]),
                "n_test": int(te.sum()),
            }
        )
    df = pd.DataFrame(rows)
    auc = np.nan
    if len(df) >= 4 and df["true_label"].nunique() > 1:
        auc = float(roc_auc_score(df["true_label"], df["mean_pos_score"]))
    return auc, df


def run_contrast(a: ad.AnnData, pos_arm: str, neg_arm: str, label: str):
    sub = a[a.obs["tv_arm"].isin([pos_arm, neg_arm])].copy()
    if sub.n_obs < 40:
        return None
    # mild ambient note only — keep published Endo
    sc.pp.highly_variable_genes(sub, n_top_genes=1500, flavor="seurat", subset=False)
    h = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=20)
    y = (sub.obs["tv_arm"].astype(str) == pos_arm).astype(int).to_numpy()
    donor = sub.obs["donor_id"].astype(str).to_numpy()
    auc, ddf = donor_oof(h.obsm["X_pca"], y, donor)
    ddf["contrast"] = label
    # ambient imbalance
    present = [g for g in NEURONAL if g in sub.var_names]
    amb = dense(sub[:, present].X).mean(1) if present else np.zeros(sub.n_obs)
    return {
        "contrast": label,
        "pos": pos_arm,
        "neg": neg_arm,
        "n_pos": int((y == 1).sum()),
        "n_neg": int((y == 0).sum()),
        "n_donors_pos": int(sub.obs.loc[y == 1, "donor_id"].nunique()),
        "n_donors_neg": int(sub.obs.loc[y == 0, "donor_id"].nunique()),
        "donor_oof_auroc": auc,
        "ambient_neuronal_pos": float(amb[y == 1].mean()) if (y == 1).any() else None,
        "ambient_neuronal_neg": float(amb[y == 0].mean()) if (y == 0).any() else None,
        "folds": ddf,
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    log(f"Reading {NN.name}")
    a = sc.read_h5ad(NN)
    a = a[a.obs["cell_type"].astype(str).str.lower().str.contains("endo")].copy()
    a = to_symbols(a)
    if a.var.index.name == "feature_name" or "feature_name" in a.var.columns:
        a.var = a.var.copy()
        if "feature_name" in a.var.columns:
            # avoid AnnData write conflict when index name == column
            a.var = a.var.drop(columns=["feature_name"], errors="ignore")
        a.var.index.name = None
    # attach TV
    tissue = a.obs["tissue"].astype(str).str.lower()
    rows = []
    for t in tissue:
        meta = TISSUE_TV.get(t)
        if meta is None:
            rows.append({"tv_arm": "other", "tv_p_ACA": np.nan, "tv_p_MCA": np.nan, "tv_p_PCA": np.nan, "tv_assign": np.nan})
        else:
            rows.append(
                {
                    "tv_arm": meta["tv_arm"],
                    "tv_p_ACA": meta["tv_p_ACA"],
                    "tv_p_MCA": meta["tv_p_MCA"],
                    "tv_p_PCA": meta["tv_p_PCA"],
                    "tv_assign": meta["assign"],
                }
            )
    tv = pd.DataFrame(rows, index=a.obs_names)
    for c in tv.columns:
        a.obs[c] = tv[c]
    a.obs["source_dataset"] = "Allen_Cortical_NonNeuronal"
    a.obs["vessel_enriched"] = False
    a.write_h5ad(OUT_H5AD, compression="gzip")
    log(f"Endo by tissue:\n{a.obs.groupby(['tissue','tv_arm']).size()}")

    contrasts = [
        ("ACA", "MCA", "ACC_ACA_vs_temporal_MCA"),  # ACC vs MTG+ANG (both MCA)
        ("ACA", "PCA", "ACC_ACA_vs_V1_PCA"),
        ("MCA", "PCA", "temporal_MCA_vs_V1_PCA"),
    ]
    # Build MCA pool = MTG+ANG only for cleaner Class B analog (exclude S1/A1/DLPFC lean)
    a2 = a.copy()
    a2.obs["tv_arm"] = a2.obs["tv_arm"].astype(str)
    mask_mca_clean = a2.obs["tissue"].astype(str).str.lower().isin(
        ["middle temporal gyrus", "angular gyrus"]
    )
    a2.obs.loc[~mask_mca_clean & (a2.obs["tv_arm"] == "MCA"), "tv_arm"] = "MCA_other"

    results = []
    fold_frames = []
    for pos, neg, lab in contrasts:
        # for ACA vs MCA use clean temporal MCA
        src = a2 if "MCA" in (pos, neg) and lab.startswith("ACC") else a
        if lab == "ACC_ACA_vs_temporal_MCA":
            sub_src = a2[a2.obs["tv_arm"].isin(["ACA", "MCA"])].copy()
        elif lab == "temporal_MCA_vs_V1_PCA":
            sub_src = a2[a2.obs["tv_arm"].isin(["MCA", "PCA"])].copy()
        else:
            sub_src = a[a.obs["tv_arm"].isin([pos, neg])].copy()
        # inline contrast on subset
        if sub_src.n_obs < 40:
            log(f"skip {lab}")
            continue
        sc.pp.highly_variable_genes(sub_src, n_top_genes=1500, flavor="seurat", subset=False)
        h = sub_src[:, sub_src.var["highly_variable"]].copy()
        sc.pp.scale(h, max_value=10)
        sc.tl.pca(h, n_comps=20)
        y = (sub_src.obs["tv_arm"].astype(str) == pos).astype(int).to_numpy()
        donor = sub_src.obs["donor_id"].astype(str).to_numpy()
        auc, ddf = donor_oof(h.obsm["X_pca"], y, donor)
        ddf["contrast"] = lab
        fold_frames.append(ddf)
        present = [g for g in NEURONAL if g in sub_src.var_names]
        amb = dense(sub_src[:, present].X).mean(1) if present else np.zeros(sub_src.n_obs)
        rec = {
            "contrast": lab,
            "pos": pos,
            "neg": neg if lab != "ACC_ACA_vs_temporal_MCA" else "MCA_temporal_MTG_ANG",
            "n_pos": int((y == 1).sum()),
            "n_neg": int((y == 0).sum()),
            "n_donors_pos": int(sub_src.obs.loc[y == 1, "donor_id"].nunique()),
            "n_donors_neg": int(sub_src.obs.loc[y == 0, "donor_id"].nunique()),
            "donor_oof_auroc": auc,
            "ambient_neuronal_pos": float(amb[y == 1].mean()),
            "ambient_neuronal_neg": float(amb[y == 0].mean()),
        }
        results.append(rec)
        log(f"{lab}: donor OOF={auc} n={rec['n_pos']}/{rec['n_neg']}")

    if fold_frames:
        pd.concat(fold_frames, ignore_index=True).to_csv(OUT_DONOR, index=False)

    summary = {
        "object": str(OUT_H5AD.relative_to(ROOT)),
        "n_endo_total": int(a.n_obs),
        "by_tissue": a.obs["tissue"].astype(str).value_counts().to_dict(),
        "by_tv_arm": a.obs["tv_arm"].astype(str).value_counts().to_dict(),
        "contrasts": results,
        "scientific_note": (
            "Within-dataset Allen Endo; unsorted cortex; TV arms are operational priors "
            "(ACC→ACA, MTG/ANG→MCA, V1→PCA). Not vessel-enriched. Donor n is small — "
            "interpret OOF cautiously."
        ),
        "verdict": (
            "WITHIN_ATLAS_SIGNAL"
            if any((r.get("donor_oof_auroc") or 0) >= 0.7 for r in results)
            else "WITHIN_ATLAS_WEAK_OR_UNDERPOWERED"
        ),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")

    lines = [
        "# Allen within-atlas Endo — territory-mapped H1 path",
        "",
        f"**Verdict:** `{summary['verdict']}`",
        "",
        "Same atlas, published Endo labels, no cross-dataset enrichment asymmetry.",
        "",
        "| Contrast | n pos/neg | donors pos/neg | Donor OOF | Ambient neur pos/neg |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in results:
        lines.append(
            f"| {r['contrast']} | {r['n_pos']}/{r['n_neg']} | {r['n_donors_pos']}/{r['n_donors_neg']} | "
            f"{r['donor_oof_auroc'] if r['donor_oof_auroc']==r['donor_oof_auroc'] else float('nan'):.3f} | "
            f"{r['ambient_neuronal_pos']:.2f}/{r['ambient_neuronal_neg']:.2f} |"
        )
    lines += [
        "",
        "## Mapping rationale",
        "- ACC → ACA (TerritoryVoxels cingulate anterior)",
        "- MTG + angular gyrus → MCA (Garcia ROSMAP Class B precedent)",
        "- V1 → PCA-leaning occipital",
        "- DLPFC/M1/S1 treated as MCA-leaning / other (not primary Class B)",
        "",
        "## Caveats",
        "- Unsorted snRNA; ambient neuronal RNA still present",
        "- Few donors (≤6); donor OOF unstable if <4 usable donors per arm",
        "- Not enrichment-matched to Garcia BVE",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
