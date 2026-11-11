"""
Class B path on disk:

1. Build Garcia ROSMAP vascular h5ad; attach TerritoryVoxels priors by region.
2. Within-Garcia contrasts (honest Class B status — no clean ACA cortex arm).
3. Harmonized Yang SFG (tight endo) x Garcia MCA endo via Harmony + DE.
4. ATLAS passphrase gate status note.

Outputs under results/class_b/ and data/processed/garcia/
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import anndata as ad
import nibabel as nib
import numpy as np
import pandas as pd
import scanpy as sc
from nilearn import datasets, image
from scipy import sparse
from scipy.io import mmread
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
MTX_DIR = ROOT / "data/processed/garcia/rosmap_mtx"
OUT_ROSMAP = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
YANG = ROOT / "data/processed/yang/yang_vascular_qc_tight.h5ad"
GARCIA_TEMP = ROOT / "data/processed/garcia/garcia_temporal_mca_qc.h5ad"
TV_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/TerritoryVoxels_BMM.nii"
)

OUT_DIR = ROOT / "results/class_b"
OUT_WITHIN = OUT_DIR / "within_garcia_class_b_de.csv"
OUT_WITHIN_MD = OUT_DIR / "within_garcia_class_b_report.md"
OUT_HARM = OUT_DIR / "yang_garcia_harmony_class_b_de.csv"
OUT_HARM_MD = OUT_DIR / "yang_garcia_harmony_class_b_report.md"
OUT_HARM_H5AD = ROOT / "data/processed/class_b/yang_garcia_harmony_endo.h5ad"
OUT_ATLAS = ROOT / "results/atlas/atlas_passphrase_gate.md"
OUT_SUMMARY = ROOT / "results/analysis_progress_summary.md"

# TerritoryVoxels labels 1-4 among assigned voxels
TV_MAP = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}

# Harvard-Oxford cortical labels pooled per Garcia brain_region
REGION_HO = {
    "Midtemporal_cortex": [
        "Middle Temporal Gyrus, anterior division",
        "Middle Temporal Gyrus, posterior division",
        "Middle Temporal Gyrus, temporooccipital part",
    ],
    "Angular_gyrus": ["Angular Gyrus"],
    "Prefrontal_cortex": ["Middle Frontal Gyrus"],  # prior proxy; MCA-leaning
    "Hippocampus": None,  # subcortical — use prior from Yang probe / skip HO cort
    "Entorhinal_cortex": [
        "Parahippocampal Gyrus, anterior division",
        "Parahippocampal Gyrus, posterior division",
    ],
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def safe_print(text: str):
    try:
        print(text)
    except UnicodeEncodeError:
        print(text.encode("ascii", "replace").decode("ascii"))


def tv_roi_vector(ho_labels: list[str] | None) -> dict:
    """TerritoryVoxels fraction among labels 1-4 inside HO ROI (pooled)."""
    if not ho_labels:
        return {
            "tv_dominant": "unknown",
            "tv_p_ACA": np.nan,
            "tv_p_MCA": np.nan,
            "tv_p_PCA": np.nan,
            "tv_p_VB": np.nan,
            "tv_p_dominant": np.nan,
            "tv_assignment_fraction": np.nan,
        }
    ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-1mm")
    labels = list(ho["labels"])
    ho_data = np.asanyarray(image.load_img(ho["maps"]).dataobj)
    mask = np.zeros(ho_data.shape, dtype=bool)
    for lab in ho_labels:
        if lab not in labels:
            raise KeyError(lab)
        mask |= ho_data == labels.index(lab)

    tv_img = nib.load(str(TV_PATH))
    # resample TV to HO grid if needed
    if tv_img.shape[:3] != ho_data.shape:
        tv_img = image.resample_to_img(tv_img, ho["maps"], interpolation="nearest")
    tv = np.asanyarray(tv_img.dataobj)
    # round float-scaled labels
    tv = np.rint(tv).astype(np.int32)
    vals = tv[mask]
    assigned = vals[np.isin(vals, [1, 2, 3, 4])]
    frac_asg = float(len(assigned) / max(1, len(vals)))
    if len(assigned) == 0:
        return {
            "tv_dominant": "unassigned",
            "tv_p_ACA": 0.0,
            "tv_p_MCA": 0.0,
            "tv_p_PCA": 0.0,
            "tv_p_VB": 0.0,
            "tv_p_dominant": 0.0,
            "tv_assignment_fraction": frac_asg,
        }
    counts = {TV_MAP[i]: float((assigned == i).sum()) for i in [1, 2, 3, 4]}
    total = sum(counts.values())
    probs = {f"tv_p_{k}": counts[k] / total for k in ["ACA", "MCA", "PCA", "VB"]}
    dom = max(counts, key=counts.get)
    return {
        "tv_dominant": dom,
        **probs,
        "tv_p_dominant": probs[f"tv_p_{dom}"],
        "tv_assignment_fraction": frac_asg,
    }


def load_rosmap() -> ad.AnnData:
    log("Loading ROSMAP MTX")
    mat = mmread(MTX_DIR / "matrix.mtx").T.tocsr()  # cells x genes
    genes = pd.read_csv(MTX_DIR / "genes.tsv", header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(MTX_DIR / "barcodes.tsv", header=None)[0].astype(str).tolist()
    meta = pd.read_csv(MTX_DIR / "meta.csv", index_col=0)
    meta.index = meta.index.astype(str)
    # align
    meta = meta.reindex(barcodes)
    adata = ad.AnnData(X=mat, obs=meta, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs["dataset"] = "Garcia_ROSMAP"
    adata.obs["celltype"] = adata.obs["celltype"].astype(str)
    adata.obs["vascular_lineage"] = (
        adata.obs["celltype"]
        .map({"Endo": "endo", "Per": "peri", "SMC": "smc", "Fib": "fib"})
        .fillna("other")
    )
    # Territory attach
    region_tv = {}
    for region, holabs in REGION_HO.items():
        if region == "Hippocampus":
            # from prior Yang TerritoryVoxels hippocampus
            region_tv[region] = {
                "tv_dominant": "PCA",
                "tv_p_ACA": 0.0,
                "tv_p_MCA": 0.25,
                "tv_p_PCA": 0.75,
                "tv_p_VB": 0.0,
                "tv_p_dominant": 0.75,
                "tv_assignment_fraction": 0.31,
                "tv_source": "prior_Yang_hippocampus_TerritoryVoxels",
            }
        else:
            log(f"TerritoryVoxels for {region}")
            region_tv[region] = {**tv_roi_vector(holabs), "tv_source": "TerritoryVoxels_HO"}
    # defaults for other regions
    for r in adata.obs["brain_region"].astype(str).unique():
        if r not in region_tv:
            region_tv[r] = {
                "tv_dominant": "unknown",
                "tv_p_ACA": np.nan,
                "tv_p_MCA": np.nan,
                "tv_p_PCA": np.nan,
                "tv_p_VB": np.nan,
                "tv_p_dominant": np.nan,
                "tv_assignment_fraction": np.nan,
                "tv_source": "unmapped",
            }
    tv_df = pd.DataFrame(region_tv).T
    tv_df.index.name = "brain_region"
    br = adata.obs["brain_region"].astype(str)
    for col in tv_df.columns:
        mapped = br.map(tv_df[col])
        if col.startswith("tv_p_") or col == "tv_assignment_fraction":
            adata.obs[col] = pd.to_numeric(mapped, errors="coerce")
        else:
            adata.obs[col] = mapped.astype(str)

    # Class B arm labels
    def arm(br: str, dom: str) -> str:
        if br in ("Midtemporal_cortex", "Angular_gyrus"):
            return "MCA_isocortex"
        if br == "Prefrontal_cortex":
            return "PFC_MCA_leaning"  # not ACA
        if br in ("Hippocampus", "Entorhinal_cortex") and dom == "PCA":
            return "PCA_medial_temporal"
        return "other"

    adata.obs["class_b_arm"] = [
        arm(str(b), str(d))
        for b, d in zip(adata.obs["brain_region"], adata.obs["tv_dominant"], strict=True)
    ]
    # h5ad-safe dtypes
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str)
    return adata, region_tv


def normalize(adata: ad.AnnData) -> ad.AnnData:
    a = adata.copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.filter_cells(a, min_genes=200)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def subsample(adata: ad.AnnData, groupby: str, groups: list[str], n_max: int, seed: int):
    rng = np.random.default_rng(seed)
    parts = []
    for g in groups:
        idx = adata.obs_names[adata.obs[groupby].astype(str) == g]
        if len(idx) == 0:
            continue
        if len(idx) > n_max:
            idx = rng.choice(idx.to_numpy(), size=n_max, replace=False)
        parts.append(adata[idx])
    return ad.concat(parts, join="inner", merge="same") if parts else adata[[]].copy()


def run_de(adata: ad.AnnData, groupby: str, group: str, reference: str, n_genes: int = 150):
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group],
        reference=reference,
        method="wilcoxon",
        n_genes=n_genes,
    )
    return sc.get.rank_genes_groups_df(adata, group=group)


def within_garcia(rosmap: ad.AnnData) -> tuple[pd.DataFrame, str]:
    """
    Best within-Garcia contrasts on ROSMAP endo.
    True ACA↔MCA Class B is unavailable (no ACA cortex arm).
    """
    endo = rosmap[rosmap.obs["vascular_lineage"].astype(str) == "endo"].copy()
    frames = []
    notes = []

    # Contrast 1: MCA isocortex (MT+AG) vs PFC (both MCA-leaning) — territory-matched region contrast
    c1 = endo[endo.obs["class_b_arm"].isin(["MCA_isocortex", "PFC_MCA_leaning"])].copy()
    if c1.obs["class_b_arm"].nunique() == 2:
        c1 = subsample(c1, "class_b_arm", ["MCA_isocortex", "PFC_MCA_leaning"], 2000, 1)
        c1.obs["class_b_arm"] = c1.obs["class_b_arm"].astype("category")
        df = run_de(c1, "class_b_arm", "MCA_isocortex", "PFC_MCA_leaning")
        df["contrast"] = "MCA_isocortex_vs_PFC_MCA_leaning"
        df["class_b_status"] = "NOT_true_ClassB_same_territory_region_contrast"
        frames.append(df)
        notes.append(
            f"MCA isocortex vs PFC: {c1.n_obs} endo cells — both MCA-leaning "
            "(region contrast, not ACA↔MCA)."
        )

    # Contrast 2: MCA isocortex vs PCA medial temporal (tissue-class differing)
    c2 = endo[endo.obs["class_b_arm"].isin(["MCA_isocortex", "PCA_medial_temporal"])].copy()
    if c2.obs["class_b_arm"].nunique() == 2:
        c2 = subsample(c2, "class_b_arm", ["MCA_isocortex", "PCA_medial_temporal"], 2000, 2)
        c2.obs["class_b_arm"] = c2.obs["class_b_arm"].astype("category")
        df = run_de(c2, "class_b_arm", "MCA_isocortex", "PCA_medial_temporal")
        df["contrast"] = "MCA_isocortex_vs_PCA_medial_temporal"
        df["class_b_status"] = "exploratory_ClassA_like_tissue_differing"
        frames.append(df)
        notes.append(
            f"MCA isocortex vs PCA medial temporal: {c2.n_obs} endo cells — "
            "territory differs but tissue class also differs."
        )

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return out, "\n".join(notes)


def harmony_yang_garcia(yang: ad.AnnData, garcia_temp: ad.AnnData, rosmap: ad.AnnData):
    """
    True Class B proxy: Yang SFG borderzone endo (control) vs Garcia MCA endo
    (ex vivo temporal + ROSMAP MT/AG), Harmony-integrated, then DE on corrected PCA.
    Also territory classifier AUROC on Harmony embedding.
    """
    y = yang[
        (yang.obs["vascular_lineage"].astype(str) == "endo")
        & (yang.obs["disease"].astype(str) == "control")
        & (yang.obs["territory_contrast"].astype(str) == "ACA_MCA_borderzone_SFG")
    ].copy()
    y.obs["class_b_group"] = "Yang_SFG_borderzone"
    y.obs["batch"] = "Yang"
    y.obs["donor_key"] = y.obs["donor_id"].astype(str)

    gt = garcia_temp[garcia_temp.obs["vascular_lineage"].astype(str) == "endo"].copy()
    if "score_endo" in gt.obs.columns:
        others = gt.obs[["score_peri", "score_smc", "score_fib"]].max(axis=1)
        gt = gt[(gt.obs["score_endo"] >= 0.25) & (gt.obs["score_endo"] - others >= 0.2)].copy()
    gt.obs["class_b_group"] = "Garcia_MCA"
    gt.obs["batch"] = "Garcia_exvivo"
    donor_col = "PatientID" if "PatientID" in gt.obs.columns else "SampleID"
    gt.obs["donor_key"] = gt.obs[donor_col].astype(str)

    gr = rosmap[
        (rosmap.obs["vascular_lineage"].astype(str) == "endo")
        & (rosmap.obs["class_b_arm"].astype(str) == "MCA_isocortex")
    ].copy()
    gr.obs["class_b_group"] = "Garcia_MCA"
    gr.obs["batch"] = "Garcia_ROSMAP"
    gr.obs["donor_key"] = gr.obs["subject"].astype(str)

    # Shared genes
    shared = y.var_names.intersection(gt.var_names).intersection(gr.var_names)
    log(f"Shared genes: {len(shared)}")
    parts = []
    for a, nmax, seed in [(y, 4000, 10), (gt, 4000, 11), (gr, 2000, 12)]:
        a = a[:, shared].copy()
        if a.n_obs > nmax:
            rng = np.random.default_rng(seed)
            a = a[rng.choice(a.obs_names.to_numpy(), nmax, replace=False)].copy()
        parts.append(a)
    comb = ad.concat(parts, join="inner", merge="same", index_unique="-")
    # Ensure log-normalized (Yang/Garcia QC already are; ROSMAP normalized in load)
    sc.pp.highly_variable_genes(comb, n_top_genes=2000, flavor="seurat", subset=False)
    comb_hvg = comb[:, comb.var["highly_variable"]].copy()
    sc.pp.scale(comb_hvg, max_value=10)
    sc.tl.pca(comb_hvg, n_comps=30, svd_solver="arpack")
    log("Running Harmony")
    sc.external.pp.harmony_integrate(comb_hvg, key="batch", basis="X_pca", adjusted_basis="X_pca_harmony")
    sc.pp.neighbors(comb_hvg, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(comb_hvg)

    # DE on raw log expression but report Harmony-aware classifier
    comb.obs["class_b_group"] = comb.obs["class_b_group"].astype("category")
    # transfer harmony coords
    comb.obsm["X_pca_harmony"] = comb_hvg.obsm["X_pca_harmony"]
    comb.obsm["X_umap"] = comb_hvg.obsm["X_umap"]
    sc.tl.rank_genes_groups(
        comb,
        groupby="class_b_group",
        groups=["Garcia_MCA"],
        reference="Yang_SFG_borderzone",
        method="wilcoxon",
        n_genes=200,
    )
    de = sc.get.rank_genes_groups_df(comb, group="Garcia_MCA")
    de["contrast"] = "Garcia_MCA_vs_Yang_SFG_borderzone_harmony"
    de["caveat"] = "Harmony batch correction; DE still on expression (not Harmony space)"

    # Classifier on Harmony PCs
    X = comb.obsm["X_pca_harmony"]
    ybin = (comb.obs["class_b_group"].astype(str) == "Garcia_MCA").astype(int).to_numpy()
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs = []
    for tr, te in skf.split(X, ybin):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], ybin[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(ybin[te], proba))
    clf_summary = {
        "harmony_cv_auroc_mean": float(np.mean(aucs)),
        "harmony_cv_auroc_std": float(np.std(aucs)),
        "n_yang": int((comb.obs["class_b_group"] == "Yang_SFG_borderzone").sum()),
        "n_garcia_mca": int((comb.obs["class_b_group"] == "Garcia_MCA").sum()),
        "n_batches": int(comb.obs["batch"].nunique()),
        "n_cells": int(comb.n_obs),
    }
    return comb, de, clf_summary


def atlas_gate_note():
    enc = ROOT / "data/raw/atlas3_training_preprocessed.tar.gz"
    head = enc.read_bytes()[:8] if enc.exists() else b""
    salted = head.startswith(b"Salted__") or head.startswith(b"U2FsdGVk")
    md = f"""# ATLAS passphrase gate

| Item | Status |
|---|---|
| Encrypted package | `{"present" if enc.exists() else "MISSING"}` ({enc}) |
| OpenSSL Salted__ | **{salted}** |
| Passphrase | **NOT SET** — awaiting ATLAS Google-form email (per-person, non-shareable) |
| Decrypt script | `analysis/atlas/decrypt_atlas3.py` |

When passphrase arrives:

```bash
set ATLAS_PASSPHRASE=your_emailed_code
python analysis/atlas/decrypt_atlas3.py
```

Then: extract -> resample lesions to MNI152 / TerritoryVoxels -> per-subject
territory lesion fractions -> nest geometry (MRA) into lesion-probability models.
"""
    OUT_ATLAS.parent.mkdir(parents=True, exist_ok=True)
    OUT_ATLAS.write_text(md, encoding="utf-8")
    return md


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_HARM_H5AD.parent.mkdir(parents=True, exist_ok=True)

    rosmap, region_tv = load_rosmap()
    rosmap = normalize(rosmap)
    rosmap.write_h5ad(OUT_ROSMAP, compression="gzip")
    log(f"Wrote {OUT_ROSMAP} ({rosmap.n_obs} cells)")
    (OUT_DIR / "garcia_rosmap_territory_by_region.json").write_text(
        json.dumps(region_tv, indent=2, default=float), encoding="utf-8"
    )

    # 1) Within-Garcia
    log("Within-Garcia contrasts")
    within_df, within_notes = within_garcia(rosmap)
    within_df.to_csv(OUT_WITHIN, index=False)
    arm_counts = (
        rosmap.obs.groupby(["class_b_arm", "vascular_lineage"], observed=True)
        .size()
        .reset_index(name="n")
    )
    within_md = [
        "# Within-Garcia Class B path — ROSMAP vascular",
        "",
        "## Critical finding",
        "",
        "Garcia has **no clean ACA isocortex arm** on disk.",
        "- Midtemporal + Angular -> MCA (TerritoryVoxels)",
        "- Prefrontal (MFG proxy) -> MCA-leaning",
        "- Hippocampus / entorhinal -> PCA-leaning (tissue class differs)",
        "",
        "True ACA↔MCA Class B therefore requires **Yang SFG × Garcia MCA** ",
        "(harmonized), not within-Garcia alone.",
        "",
        "## TerritoryVoxels by region",
        "",
        "```json",
        json.dumps(region_tv, indent=2, default=float),
        "```",
        "",
        "## Arm × lineage counts",
        "",
        arm_counts.to_string(index=False),
        "",
        "## Contrasts run",
        "",
        within_notes,
        "",
    ]
    if len(within_df):
        for contrast in within_df["contrast"].unique():
            sub = within_df[within_df["contrast"] == contrast].head(12)
            within_md += [
                f"### {contrast}",
                "",
                sub[["names", "scores", "logfoldchanges", "pvals_adj", "class_b_status"]].to_string(
                    index=False
                ),
                "",
            ]
    OUT_WITHIN_MD.write_text("\n".join(within_md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_WITHIN_MD}")

    # 2) Harmonized Yang × Garcia
    log("Loading Yang + Garcia temporal for Harmony Class B")
    yang = sc.read_h5ad(YANG)
    garcia_temp = sc.read_h5ad(GARCIA_TEMP)
    comb, de, clf_sum = harmony_yang_garcia(yang, garcia_temp, rosmap)
    comb.write_h5ad(OUT_HARM_H5AD, compression="gzip")
    de.to_csv(OUT_HARM, index=False)
    top = de.head(20)
    bot = de.tail(15)
    harm_md = [
        "# Harmonized Class B — Yang SFG borderzone vs Garcia MCA (Harmony)",
        "",
        "Batches: Yang | Garcia ex vivo temporal | Garcia ROSMAP Midtemporal+Angular.",
        "Endothelial only; Yang controls; Harmony on PCA.",
        "",
        "## Separability on Harmony embedding",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| 5-fold CV AUROC | {clf_sum['harmony_cv_auroc_mean']:.3f} +/- {clf_sum['harmony_cv_auroc_std']:.3f} |",
        f"| N Yang SFG | {clf_sum['n_yang']} |",
        f"| N Garcia MCA | {clf_sum['n_garcia_mca']} |",
        f"| N cells / batches | {clf_sum['n_cells']} / {clf_sum['n_batches']} |",
        "",
        "High AUROC after Harmony = residual territorial/regional signal beyond batch.",
        "Yang arm is ACA/MCA **borderzone** (not pure ACA) — Aim3/Class B hybrid.",
        "",
        "## Top genes higher in Garcia MCA",
        "",
        top[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False),
        "",
        "## Higher in Yang SFG (bottom ranks)",
        "",
        bot[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string(index=False),
        "",
        f"Object: `{OUT_HARM_H5AD.relative_to(ROOT).as_posix()}`",
    ]
    OUT_HARM_MD.write_text("\n".join(harm_md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_HARM_MD}")

    # 3) ATLAS gate
    atlas_md = atlas_gate_note()
    log(f"Wrote {OUT_ATLAS}")

    summary = f"""# Analysis progress — existing data

## Class B (this pass)

### Within-Garcia
- **No clean ACA cortex arm** — true ACA↔MCA Class B blocked within Garcia alone
- ROSMAP object: `data/processed/garcia/garcia_rosmap_vascular.h5ad`
- MCA isocortex (Midtemporal+Angular) available; PFC is MCA-leaning
- Report: `results/class_b/within_garcia_class_b_report.md`

### Harmonized Yang × Garcia (primary Class B path)
- Harmony CV AUROC: **{clf_sum['harmony_cv_auroc_mean']:.3f}** +/- {clf_sum['harmony_cv_auroc_std']:.3f}
- Cells: Yang SFG {clf_sum['n_yang']} vs Garcia MCA {clf_sum['n_garcia_mca']}
- DE: `results/class_b/yang_garcia_harmony_class_b_de.csv`
- Report: `results/class_b/yang_garcia_harmony_class_b_report.md`
- Object: `data/processed/class_b/yang_garcia_harmony_endo.h5ad`

## ATLAS
- Still blocked on passphrase — `results/atlas/atlas_passphrase_gate.md`
- Decrypt ready: `analysis/atlas/decrypt_atlas3.py`

## Still waiting
- ATLAS encryption code email
- RADC clinical-pathologic approval
- ADNI approval
"""
    OUT_SUMMARY.write_text(summary, encoding="utf-8")
    safe_print(summary)
    safe_print("\n" + atlas_md)


if __name__ == "__main__":
    main()
