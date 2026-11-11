"""
Next molecular pass (no ATLAS/RADC):

1. Within-ROSMAP endothelial contrasts (published Endo subtypes + neuronal caps):
   - MCA isocortex (Midtemporal+Angular) vs PFC (same-territory region null)
   - MCA isocortex vs PCA medial temporal (territory+tissue)
   - Capillary-endo (capEndo) only versions to condition on A-V zone (H1 confound control)

2. Donor/subject-aware Class B on tight Yang×Garcia Harmony object:
   - Leave-one-batch-out and subject-level score summaries

Writes under results/h1/, results/class_b/, docs/
"""
from __future__ import annotations

import re
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
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
ROSMAP = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
HARM = ROOT / "data/processed/class_b/yang_garcia_harmony_endo_tight.h5ad"

OUT_H1_DE = ROOT / "results/h1/rosmap_endo_region_territory_de.csv"
OUT_H1_MD = ROOT / "results/h1/rosmap_within_dataset_h1_report.md"
OUT_CB = ROOT / "results/class_b/yang_garcia_tight_donor_aware.csv"
OUT_CB_MD = ROOT / "results/class_b/yang_garcia_tight_donor_aware_report.md"
OUT_MEMO = ROOT / "docs/09_scientific_status_update_2026-07-24.md"
OUT_SUMMARY = ROOT / "results/analysis_progress_summary.md"

NEURONAL_CAPS = {"LINGO1": 0.5, "RBFOX1": 0.5, "KCNIP4": 0.5, "SNAP25": 0.5}
ENDO_MARKERS = ["CLDN5", "PECAM1", "FLT1", "VWF", "CDH5"]
ENDO_SUBTYPES = {"aEndo", "capEndo", "vEndo"}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def expr_ok(adata: ad.AnnData) -> pd.Series:
    genes = list(NEURONAL_CAPS) + ENDO_MARKERS
    present = [g for g in genes if g in adata.var_names]
    X = dense(adata[:, present].X)
    df = pd.DataFrame(X, index=adata.obs_names, columns=present)
    contam = pd.Series(False, index=adata.obs_names)
    for g, thr in NEURONAL_CAPS.items():
        if g in df.columns:
            contam |= df[g] > thr
    endo_pos = pd.Series(False, index=adata.obs_names)
    for g in ENDO_MARKERS:
        if g in df.columns:
            endo_pos |= df[g] > 0
    return (~contam) & endo_pos


def subsample(adata: ad.AnnData, groupby: str, groups: list[str], n_max: int, seed: int):
    rng = np.random.default_rng(seed)
    parts = []
    for g in groups:
        idx = adata.obs_names[adata.obs[groupby].astype(str) == g]
        if len(idx) == 0:
            continue
        if len(idx) > n_max:
            idx = rng.choice(idx.to_numpy(), n_max, replace=False)
        parts.append(adata[idx])
    return ad.concat(parts, join="inner", merge="same") if parts else adata[[]].copy()


def run_de(adata, groupby, group, reference, n_genes=150):
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group],
        reference=reference,
        method="wilcoxon",
        n_genes=n_genes,
    )
    return sc.get.rank_genes_groups_df(adata, group=group)


def within_rosmap_h1(rosmap: ad.AnnData) -> tuple[pd.DataFrame, dict]:
    """Published Endo subtypes + neuronal caps; region/territory contrasts."""
    endo = rosmap[rosmap.obs["subtype"].astype(str).isin(ENDO_SUBTYPES)].copy()
    keep = expr_ok(endo)
    endo = endo[keep].copy()
    endo.obs["arm"] = endo.obs["class_b_arm"].astype(str)
    endo.obs["zonation"] = endo.obs["subtype"].astype(str)

    frames = []
    counts = {
        "n_endo_raw_subtypes": int(rosmap.obs["subtype"].astype(str).isin(ENDO_SUBTYPES).sum()),
        "n_endo_after_neuronal_caps": int(endo.n_obs),
        "by_arm": endo.obs["arm"].value_counts().to_dict(),
        "by_zonation": endo.obs["zonation"].value_counts().to_dict(),
    }

    contrasts = [
        (
            "MCA_isocortex",
            "PFC_MCA_leaning",
            "MCA_isocortex_vs_PFC_same_territory_region_null",
            "same_territory_region_contrast",
        ),
        (
            "MCA_isocortex",
            "PCA_medial_temporal",
            "MCA_isocortex_vs_PCA_medial_temporal",
            "territory_and_tissue_differ",
        ),
    ]

    # All endo subtypes pooled
    for group, ref, name, status in contrasts:
        sub = endo[endo.obs["arm"].isin([group, ref])].copy()
        if sub.obs["arm"].nunique() < 2:
            continue
        n_max = 800
        sub = subsample(sub, "arm", [group, ref], n_max, seed=abs(hash(name)) % 10_000)
        sub.obs["arm"] = sub.obs["arm"].astype("category")
        df = run_de(sub, "arm", group, ref)
        df["contrast"] = name
        df["stratum"] = "all_endo_subtypes"
        df["h1_status"] = status
        df["n_cells"] = sub.n_obs
        frames.append(df)
        log(f"ROSMAP DE {name} all-endo: {sub.n_obs}")

    # Capillary-only (zonation-conditioned)
    cap = endo[endo.obs["zonation"] == "capEndo"].copy()
    for group, ref, name, status in contrasts:
        sub = cap[cap.obs["arm"].isin([group, ref])].copy()
        if sub.obs["arm"].nunique() < 2:
            continue
        if min((sub.obs["arm"] == group).sum(), (sub.obs["arm"] == ref).sum()) < 30:
            log(f"Skip capEndo {name}: too few cells")
            continue
        sub = subsample(sub, "arm", [group, ref], 500, seed=abs(hash("cap" + name)) % 10_000)
        sub.obs["arm"] = sub.obs["arm"].astype("category")
        df = run_de(sub, "arm", group, ref)
        df["contrast"] = name
        df["stratum"] = "capEndo_only"
        df["h1_status"] = status
        df["n_cells"] = sub.n_obs
        frames.append(df)
        log(f"ROSMAP DE {name} capEndo: {sub.n_obs}")

    # Simple variance: PCA on HVGs, R2 of arm after subject
    h1_var = {}
    for stratum_name, subset in [("all_endo", endo), ("capEndo", cap)]:
        use = subset[subset.obs["arm"].isin(["MCA_isocortex", "PFC_MCA_leaning", "PCA_medial_temporal"])]
        if use.n_obs < 100:
            continue
        use = subsample(
            use,
            "arm",
            ["MCA_isocortex", "PFC_MCA_leaning", "PCA_medial_temporal"],
            600,
            seed=9,
        )
        sc.pp.highly_variable_genes(use, n_top_genes=1000, flavor="seurat", subset=False)
        X = dense(use[:, use.var["highly_variable"]].X)
        X = X - X.mean(0, keepdims=True)
        u, s, _ = np.linalg.svd(X, full_matrices=False)
        pcs = u[:, :20] * s[:20]
        subj = pd.get_dummies(use.obs["subject"].astype(str), drop_first=True).to_numpy(float)
        arm = pd.get_dummies(use.obs["arm"].astype(str), drop_first=True).to_numpy(float)

        def mean_r2(Y, Z):
            from sklearn.linear_model import LinearRegression

            r2s = []
            for j in range(Y.shape[1]):
                lr = LinearRegression().fit(Z, Y[:, j])
                pred = lr.predict(Z)
                ss_res = ((Y[:, j] - pred) ** 2).sum()
                ss_tot = ((Y[:, j] - Y[:, j].mean()) ** 2).sum()
                r2s.append(1 - ss_res / ss_tot if ss_tot > 0 else 0)
            return float(np.mean(r2s))

        r2_s = mean_r2(pcs, subj) if subj.size else 0.0
        r2_sa = mean_r2(pcs, np.hstack([subj, arm]))
        h1_var[stratum_name] = {
            "n_cells": int(use.n_obs),
            "r2_subject": r2_s,
            "r2_subject_arm": r2_sa,
            "delta_arm_after_subject": r2_sa - r2_s,
        }
        log(f"Variance {stratum_name}: delta_arm={r2_sa - r2_s:.4f}")

    counts["variance_partition"] = h1_var
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return out, counts


def recover_donor_keys(harm: ad.AnnData) -> pd.Series:
    """Recover donor/library IDs from cell barcodes when donor_key missing."""
    keys = []
    for name, batch in zip(harm.obs_names.astype(str), harm.obs["batch"].astype(str), strict=True):
        if batch == "Yang":
            m = re.match(r"(GSM\d+)", name)
            keys.append(m.group(1) if m else name.split("_")[0])
        elif batch == "Garcia_exvivo":
            keys.append(name.split(".")[0])  # D20-125061
        else:
            # ROSMAP: Region.BARCODE-lib-...
            parts = name.split(".")
            keys.append(parts[0] + "_" + name.rsplit("-", 2)[-2] if "-" in name else name)
    return pd.Series(keys, index=harm.obs_names, name="donor_key")


def donor_aware_class_b(harm: ad.AnnData) -> tuple[pd.DataFrame, dict]:
    """
    On Harmony embedding with out-of-fold scores:
    - 5-fold cell CV
    - Leave-one-donor-out AUROC (donors recovered from barcodes)
    - Donor-mean OOF scores -> donor-level AUROC
    """
    if "X_pca_harmony" not in harm.obsm:
        raise SystemExit("Harmony embedding missing on tight object")
    X = harm.obsm["X_pca_harmony"]
    y = (harm.obs["class_b_group"].astype(str) == "Garcia_MCA").astype(int).to_numpy()
    batch = harm.obs["batch"].astype(str).to_numpy()
    if "donor_key" in harm.obs:
        donor = harm.obs["donor_key"].astype(str).to_numpy()
    else:
        donor = recover_donor_keys(harm).to_numpy()
        harm.obs["donor_key"] = donor

    rows = []
    # Cell-level CV
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    cell_aucs = []
    for i, (tr, te) in enumerate(skf.split(X, y)):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        auc = roc_auc_score(y[te], proba)
        cell_aucs.append(auc)
        rows.append({"fold_type": "cell_cv", "held_out": f"fold{i}", "auroc": auc, "n_test": len(te)})

    # Leave-one-donor-out (only donors with >=5 cells)
    donors = pd.Series(donor)
    vc = donors.value_counts()
    usable = [d for d in vc.index if vc[d] >= 5]
    oof = np.full(len(y), np.nan)
    lodo_aucs = []
    for d in usable:
        te = donor == d
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        oof[te] = proba
        if len(np.unique(y[te])) > 1:
            auc = roc_auc_score(y[te], proba)
            lodo_aucs.append(auc)
            rows.append(
                {
                    "fold_type": "leave_one_donor",
                    "held_out": d,
                    "auroc": auc,
                    "n_test": int(te.sum()),
                    "batch": batch[te][0],
                }
            )
        else:
            rows.append(
                {
                    "fold_type": "leave_one_donor",
                    "held_out": d,
                    "auroc": np.nan,
                    "mean_garcia_score": float(proba.mean()),
                    "n_test": int(te.sum()),
                    "batch": batch[te][0],
                    "true_label": int(y[te][0]),
                    "note": "single-class donor (OOF score only)",
                }
            )

    # Donor-level AUROC from OOF mean scores
    drows = []
    for d in usable:
        m = donor == d
        if np.isnan(oof[m]).all():
            continue
        drows.append(
            {
                "fold_type": "donor_oof_score",
                "held_out": d,
                "mean_garcia_score": float(np.nanmean(oof[m])),
                "true_label": int(y[m][0]),
                "n_test": int(m.sum()),
                "batch": batch[m][0],
            }
        )
    rows.extend(drows)
    ddf = pd.DataFrame(drows)
    donor_auc = np.nan
    if len(ddf) >= 4 and ddf["true_label"].nunique() > 1:
        donor_auc = float(roc_auc_score(ddf["true_label"], ddf["mean_garcia_score"]))

    summary = {
        "cell_cv_auroc_mean": float(np.mean(cell_aucs)) if cell_aucs else np.nan,
        "cell_cv_auroc_std": float(np.std(cell_aucs)) if cell_aucs else np.nan,
        "lodo_within_donor_auroc_mean": float(np.mean(lodo_aucs)) if lodo_aucs else np.nan,
        "n_lodo_multiclass_donors": len(lodo_aucs),
        "donor_level_oof_auroc": donor_auc,
        "n_donors_scored": int(len(ddf)),
        "n_yang_donors": int((ddf["true_label"] == 0).sum()) if len(ddf) else 0,
        "n_garcia_donors": int((ddf["true_label"] == 1).sum()) if len(ddf) else 0,
    }
    return pd.DataFrame(rows), summary


def top_genes(de: pd.DataFrame, contrast: str, stratum: str, n=10) -> str:
    sub = de[(de["contrast"] == contrast) & (de["stratum"] == stratum)].head(n)
    if sub.empty:
        return "_none_"
    return sub[["names", "logfoldchanges", "pvals_adj"]].to_string(index=False)


def main():
    OUT_H1_DE.parent.mkdir(parents=True, exist_ok=True)
    OUT_CB.parent.mkdir(parents=True, exist_ok=True)

    log("Loading ROSMAP")
    rosmap = sc.read_h5ad(ROSMAP)
    de, counts = within_rosmap_h1(rosmap)
    de.to_csv(OUT_H1_DE, index=False)

    var = counts.get("variance_partition", {})
    h1_md = [
        "# Within-ROSMAP H1 path — published Endo subtypes + neuronal caps",
        "",
        "Garcia ROSMAP is the only **within-dataset** multi-region vascular object on disk.",
        "No ACA cortex arm; contrasts are:",
        "1. MCA isocortex (Midtemporal+Angular) vs PFC (both MCA-leaning) — **region null**",
        "2. MCA isocortex vs PCA medial temporal — territory differs but tissue class also differs",
        "",
        "Neuronal caps: LINGO1/RBFOX1/KCNIP4/SNAP25 ≤ 0.5; require canonical endo marker.",
        "capEndo stratum conditions on A–V zonation (H1 confound control).",
        "",
        "## Counts",
        "",
        f"- Endo subtypes raw: **{counts['n_endo_raw_subtypes']:,}**",
        f"- After neuronal caps: **{counts['n_endo_after_neuronal_caps']:,}**",
        f"- By arm: `{counts['by_arm']}`",
        f"- By zonation: `{counts['by_zonation']}`",
        "",
        "## Variance partition (arm after subject)",
        "",
    ]
    for k, v in var.items():
        h1_md.append(
            f"- **{k}**: delta R2 arm after subject = **{v['delta_arm_after_subject']:.4f}** "
            f"(subject R2={v['r2_subject']:.4f}; n={v['n_cells']})"
        )
    h1_md += ["", "## Top DE genes", ""]
    if len(de):
        for contrast in de["contrast"].unique():
            for stratum in de.loc[de["contrast"] == contrast, "stratum"].unique():
                h1_md += [
                    f"### {contrast} | {stratum}",
                    "",
                    top_genes(de, contrast, stratum),
                    "",
                ]
    OUT_H1_MD.write_text("\n".join(h1_md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_H1_MD}")

    log("Donor-aware Class B on tight Harmony object")
    harm = sc.read_h5ad(HARM)
    cb_df, cb_sum = donor_aware_class_b(harm)
    cb_df.to_csv(OUT_CB, index=False)
    cb_md = f"""# Donor-aware Class B — tight Harmony object (out-of-fold)

Donors recovered from barcodes (Yang GSM; Garcia D20-*; ROSMAP region_lib).
Scores are **leave-one-donor-out** (not in-sample).

| Metric | Value |
|---|---:|
| Cell-level CV AUROC | {cb_sum['cell_cv_auroc_mean']:.3f} +/- {cb_sum['cell_cv_auroc_std']:.3f} |
| Donor-level OOF AUROC | {cb_sum['donor_level_oof_auroc']:.3f} |
| N Yang / Garcia donors scored | {cb_sum['n_yang_donors']} / {cb_sum['n_garcia_donors']} |
| Multiclass within-donor LODO folds | {cb_sum['n_lodo_multiclass_donors']} |

Donor-level OOF AUROC tests whether a held-out donor's mean score ranks Yang vs Garcia correctly.
Near-chance = little transferable territorial signal beyond batch/donor structure.

## Per-fold / donor table

{cb_df.to_string(index=False)}
"""
    OUT_CB_MD.write_text(cb_md, encoding="utf-8")
    log(f"Wrote {OUT_CB_MD}")

    # Scientific memo
    delta_all = var.get("all_endo", {}).get("delta_arm_after_subject", float("nan"))
    delta_cap = var.get("capEndo", {}).get("delta_arm_after_subject", float("nan"))
    memo = f"""# Scientific status update — 2026-07-24

## Centre of gravity (unchanged)
Vascular cells may carry territorial molecular identity partly independent of anatomy and
A-V zonation; it may degrade in disease and predict lesion geography (H1-H4).

## What open data now support / refute

### Yang (PCA hipp vs SFG borderzone)
- Tight endo/peri gates removed LINGO1 contamination from Aim3 DE.
- **H4 not supported on Yang alone:** within-donor AUROC high in both control and AD;
  cross-donor LODO ~chance. Region=territory confound remains.

### Garcia Class B
- **No ACA isocortex arm** within Garcia (PFC/MT/AG all MCA by TerritoryVoxels).
- Soft YangxGarcia Harmony DE was **neuron-contamination artifact** (LINGO1/KCNIP4/RBFOX1).
- After Garcia endo gate v2b: neuronal hits gone; Harmony AUROC 0.85 -> **0.68**;
  n Garcia MCA endo = 319. Residual signal weaker; still exploratory.

### Within-ROSMAP (this pass) — best within-dataset H1 path on disk
- Published Endo subtypes + neuronal caps.
- Arm-after-subject delta R2: all-endo **{delta_all:.4f}**; capEndo **{delta_cap:.4f}**.
- MCA vs PFC is a **same-territory region null**, not true ACA<->MCA Class B.
- Report: `results/h1/rosmap_within_dataset_h1_report.md`

### Donor-aware tight Class B (OOF)
- Cell CV AUROC **{cb_sum['cell_cv_auroc_mean']:.3f}**; donor-level OOF AUROC **{cb_sum['donor_level_oof_auroc']:.3f}**.
- Report: `results/class_b/yang_garcia_tight_donor_aware_report.md`

## Still blocked
- ATLAS passphrase (lesion x TerritoryVoxels Aim 4)
- RADC ROSMAP clinical-pathologic (link Garcia subject / ROSMAP IDs)
- ADNI

## Design implications locked
1. Yang = Aim3 borderzone contrast, not clean Class B.
2. True ACA<->MCA Class B needs an ACA cortex cohort not present in Garcia open data.
3. Soft cross-dataset DE without neuronal gates is not interpretable.
4. Prefer within-dataset + zonation-conditioned tests for confirmatory H1.
"""
    OUT_MEMO.write_text(memo, encoding="utf-8")
    log(f"Wrote {OUT_MEMO}")

    summary = f"""# Analysis progress — existing data

## This pass
- Within-ROSMAP H1: `results/h1/rosmap_within_dataset_h1_report.md`
  - Arm-after-subject delta R2: all-endo {delta_all:.4f}; capEndo {delta_cap:.4f}
- Donor-aware tight Class B (OOF): cell CV {cb_sum['cell_cv_auroc_mean']:.3f}; donor AUROC {cb_sum['donor_level_oof_auroc']:.3f}
- Status memo: `docs/09_scientific_status_update_2026-07-24.md`

## Prior (tight Class B)
- Garcia temporal 207 + ROSMAP MCA 112 tight endo; Harmony AUROC 0.682
- Soft Class B neuronal DE was artifact

## Blocked
- ATLAS passphrase / RADC / ADNI
"""
    OUT_SUMMARY.write_text(summary, encoding="utf-8")
    try:
        print(summary)
        print(memo)
    except UnicodeEncodeError:
        print(summary.encode("ascii", "replace").decode("ascii"))


if __name__ == "__main__":
    main()
