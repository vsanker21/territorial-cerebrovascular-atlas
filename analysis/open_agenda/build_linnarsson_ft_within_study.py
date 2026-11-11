"""
Build Linnarsson/Siletti within-study frontal/temporal vascular object +
TerritoryVoxels-mapped Endo contrasts (ACA cingulate vs lateral F vs temporal).

Scientific locks:
- Same study / chemistry / annotation ontology (Siletti HBCA CELLxGENE).
- ACA pool: A24, A32, A25, A23 (A14 excluded: MCA orbital).
- Lateral frontal: MFG A46 + IFG A44–45 (MCA-dominant HO ROIs).
- Temporal: MTG + STG (MCA-dominant).
- Prefer donor OOF AUROC; with n_donors≈3 this is underpowered — report honestly.
- Not vessel-enriched; ambient neuronal scores reported by arm.
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/aca_cohort_hunt/linnarsson_acc"
OUT_H5AD = ROOT / "data/processed/aca_cohort_hunt/linnarsson_ft_vascular.h5ad"
OUT_DIR = ROOT / "results/h1"
OUT_JSON = OUT_DIR / "linnarsson_within_study_summary.json"
OUT_MD = OUT_DIR / "linnarsson_within_study_report.md"
OUT_DONOR = OUT_DIR / "linnarsson_within_study_donor_oof.csv"
OUT_TV = ROOT / "results/aca_cohort_hunt/linnarsson_ft_tv_priors.json"
TV_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/TerritoryVoxels_BMM.nii"
)

DISSECTIONS = {
    "A24": {
        "file": "Cx_CgGr_A24.h5ad",
        "arm": "ACA",
        "ho_labels": ["Cingulate Gyrus, anterior division"],
        "role": "primary_ACA_cingulate",
    },
    "A32": {
        "file": "Cx_A32_RoG.h5ad",
        "arm": "ACA",
        "ho_labels": ["Paracingulate Gyrus", "Cingulate Gyrus, anterior division"],
        "role": "ACA_dorsal_MFC",
    },
    "A25": {
        "file": "Cx_A25_SCG.h5ad",
        "arm": "ACA",
        "ho_labels": ["Frontal Medial Cortex", "Cingulate Gyrus, anterior division"],
        "role": "ACA_subcallosal_borderline",
    },
    "A23": {
        "file": "Cx_A23_CgGC.h5ad",
        "arm": "ACA",
        "ho_labels": ["Cingulate Gyrus, posterior division"],
        "role": "ACA_caudal_cingulate",
    },
    "MFG_A46": {
        "file": "Cx_MFG_A46.h5ad",
        "arm": "lateral_frontal",
        "ho_labels": ["Middle Frontal Gyrus"],
        "role": "MCA_lateral_frontal_A46",
    },
    "IFG": {
        "file": "Cx_IFG_A44A45.h5ad",
        "arm": "lateral_frontal",
        "ho_labels": ["Inferior Frontal Gyrus, pars opercularis", "Inferior Frontal Gyrus, pars triangularis"],
        "role": "MCA_lateral_frontal_IFG",
    },
    "MTG": {
        "file": "Cx_MTG.h5ad",
        "arm": "temporal",
        "ho_labels": [
            "Middle Temporal Gyrus, anterior division",
            "Middle Temporal Gyrus, posterior division",
            "Middle Temporal Gyrus, temporooccipital part",
        ],
        "role": "MCA_temporal_MTG",
    },
    "STG": {
        "file": "Cx_STG.h5ad",
        "arm": "temporal",
        "ho_labels": ["Superior Temporal Gyrus, anterior division", "Superior Temporal Gyrus, posterior division"],
        "role": "MCA_temporal_STG",
    },
}

VASC_TYPES = {
    "endothelial cell",
    "pericyte",
    "vascular associated smooth muscle cell",
    "fibroblast",
}
TV_MAP = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
NEURONAL = ["LINGO1", "RBFOX1", "KCNIP4", "SNAP25"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def tv_roi_vector(ho_labels: list[str]) -> dict:
    ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-1mm")
    labels = list(ho["labels"])
    ho_data = np.asanyarray(image.load_img(ho["maps"]).dataobj)
    mask = np.zeros(ho_data.shape, dtype=bool)
    missing = []
    for lab in ho_labels:
        if lab not in labels:
            missing.append(lab)
            continue
        mask |= ho_data == labels.index(lab)
    if not mask.any():
        return {
            "tv_dominant": "unassigned",
            "tv_p_ACA": 0.0,
            "tv_p_MCA": 0.0,
            "tv_p_PCA": 0.0,
            "tv_p_VB": 0.0,
            "tv_p_dominant": 0.0,
            "tv_assignment_fraction": 0.0,
            "ho_missing": missing,
        }
    tv_img = nib.load(str(TV_PATH))
    if tv_img.shape[:3] != ho_data.shape:
        tv_img = image.resample_to_img(tv_img, ho["maps"], interpolation="nearest")
    tv = np.rint(np.asanyarray(tv_img.dataobj)).astype(np.int32)
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
            "ho_missing": missing,
        }
    counts = {TV_MAP[k]: float((assigned == k).sum()) for k in TV_MAP}
    total = float(len(assigned))
    probs = {f"tv_p_{k}": counts[k] / total for k in counts}
    dom = max(counts, key=counts.get)
    return {
        "tv_dominant": dom,
        **probs,
        "tv_p_dominant": probs[f"tv_p_{dom}"],
        "tv_assignment_fraction": frac_asg,
        "ho_missing": missing,
        "tv_source": "TerritoryVoxels_HO",
    }


def to_symbols(a: ad.AnnData) -> ad.AnnData:
    """Map Ensembl → symbol using feature_name / Gene if present, else leave IDs."""
    a = a.copy()
    sym_col = None
    for c in ("feature_name", "Gene", "gene_symbols", "symbol"):
        if c in a.var.columns:
            sym_col = c
            break
    if sym_col is None:
        # recover from a raw sibling slice
        raw = RAW / "Cx_MFG_A46.h5ad"
        if raw.exists():
            ref = sc.read_h5ad(raw, backed="r")
            # align by ensembl index
            idx = a.var_names.astype(str)
            if "feature_name" in ref.var.columns:
                m = ref.var["feature_name"].astype(str).reindex(idx)
                a.var["feature_name"] = m.values
                sym_col = "feature_name"
            del ref
    if sym_col is None:
        return a
    sym = a.var[sym_col].astype(str).to_numpy()
    ok = np.array([s not in {"nan", "NA", ""} and not str(s).startswith("ENSG") for s in sym])
    a = a[:, ok].copy()
    a.var_names = pd.Index(sym[ok])
    a.var_names_make_unique()
    return a


def donor_lodo_cell_auroc(X, y, donor, min_cells=5):
    """Leave-one-donor-out with cell-level AUROC.

    Siletti donors contribute cells to *all* arms, so donor-mean single-label
    OOF (as used for cross-cohort Class B) is invalid here.
    """
    rows = []
    y_all, p_all = [], []
    for name, n in pd.Series(donor).value_counts().items():
        if n < min_cells:
            continue
        te = donor == name
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        fold_auc = float("nan")
        if len(np.unique(y[te])) > 1 and te.sum() >= 4:
            fold_auc = float(roc_auc_score(y[te], proba))
        rows.append(
            {
                "held_out": name,
                "n_test": int(te.sum()),
                "n_pos_test": int(y[te].sum()),
                "n_neg_test": int((1 - y[te]).sum()),
                "cell_auroc": fold_auc,
                "mean_pos_score": float(proba.mean()),
            }
        )
        y_all.extend(y[te].tolist())
        p_all.extend(proba.tolist())
    df = pd.DataFrame(rows)
    pooled = float("nan")
    if len(y_all) >= 8 and len(set(y_all)) > 1:
        pooled = float(roc_auc_score(y_all, p_all))
    mean_fold = float(np.nanmean(df["cell_auroc"])) if len(df) and df["cell_auroc"].notna().any() else float("nan")
    return pooled, mean_fold, df


def run_contrast(a: ad.AnnData, pos: str, neg: str, label: str):
    sub = a[a.obs["study_arm"].isin([pos, neg])].copy()
    if sub.n_obs < 20:
        return None, None
    sc.pp.highly_variable_genes(sub, n_top_genes=min(1500, sub.n_vars - 1), flavor="seurat", subset=False)
    h = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(h, max_value=10)
    n_comp = min(20, h.n_obs - 2, h.n_vars - 1)
    if n_comp < 5:
        return None, None
    sc.tl.pca(h, n_comps=n_comp)
    y = (sub.obs["study_arm"].astype(str) == pos).astype(int).to_numpy()
    donor = sub.obs["donor_id"].astype(str).to_numpy()
    pooled, mean_fold, ddf = donor_lodo_cell_auroc(h.obsm["X_pca"], y, donor)
    ddf["contrast"] = label
    present = [g for g in NEURONAL if g in sub.var_names]
    amb = dense(sub[:, present].X).mean(1) if present else np.zeros(sub.n_obs)
    return {
        "contrast": label,
        "pos": pos,
        "neg": neg,
        "n_pos": int((y == 1).sum()),
        "n_neg": int((y == 0).sum()),
        "n_donors_pos": int(sub.obs.loc[y == 1, "donor_id"].nunique()),
        "n_donors_neg": int(sub.obs.loc[y == 0, "donor_id"].nunique()),
        "donor_lodo_pooled_cell_auroc": pooled,
        "donor_lodo_mean_fold_auroc": mean_fold,
        "n_donor_folds": int(len(ddf)),
        "ambient_neuronal_mean_pos": float(amb[y == 1].mean()) if (y == 1).any() else None,
        "ambient_neuronal_mean_neg": float(amb[y == 0].mean()) if (y == 0).any() else None,
        "n_ambient_genes": len(present),
        "scientific_note": (
            "Same 3 Siletti donors contribute to all arms; metric = LODO pooled "
            "cell-level AUROC (not cross-cohort donor-mean OOF)."
        ),
    }, ddf


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)

    missing = [c for c, m in DISSECTIONS.items() if not (RAW / m["file"]).exists()]
    if missing:
        raise FileNotFoundError(f"Missing CELLxGENE slices: {missing}")

    tv_priors = {}
    for code, meta in DISSECTIONS.items():
        log(f"TerritoryVoxels for {code}")
        tv_priors[code] = {
            **tv_roi_vector(meta["ho_labels"]),
            "ho_labels": meta["ho_labels"],
            "arm": meta["arm"],
            "role": meta["role"],
        }
    OUT_TV.write_text(json.dumps(tv_priors, indent=2), encoding="utf-8")

    pieces = []
    per = {}
    for code, meta in DISSECTIONS.items():
        path = RAW / meta["file"]
        log(f"Extract vascular {path.name}")
        a = sc.read_h5ad(path)
        keep = a.obs["cell_type"].astype(str).isin(VASC_TYPES)
        b = a[keep].copy()
        # keep symbol map through concat
        if "feature_name" in a.var.columns and "feature_name" not in b.var.columns:
            b.var["feature_name"] = a.var["feature_name"]
        elif "feature_name" in a.var.columns:
            b.var["feature_name"] = a.var.loc[b.var_names, "feature_name"].values
        b.obs["dissection_code"] = code
        b.obs["study_arm"] = meta["arm"]
        b.obs["dissection_role"] = meta["role"]
        b.obs["source_dataset"] = "linnarsson_siletti_hbca_cellxgene"
        tv = tv_priors[code]
        for k, v in tv.items():
            if k.startswith("tv_") or k in {"role", "arm"}:
                b.obs[k] = json.dumps(v) if isinstance(v, list) else v
        per[code] = {
            "n_vascular": int(b.n_obs),
            "n_endo": int((b.obs["cell_type"].astype(str) == "endothelial cell").sum()),
            "n_donors": int(b.obs["donor_id"].nunique()) if b.n_obs else 0,
            "arm": meta["arm"],
            "tv_dominant": tv["tv_dominant"],
            "tv_p_ACA": tv.get("tv_p_ACA"),
            "tv_p_MCA": tv.get("tv_p_MCA"),
            "tv_assignment_fraction": tv.get("tv_assignment_fraction"),
            "ho_missing": tv.get("ho_missing"),
        }
        pieces.append(b)
        del a

    all_v = ad.concat(pieces, join="outer", merge="same", index_unique="-")
    all_v.write_h5ad(OUT_H5AD)
    log(f"Wrote {OUT_H5AD} n={all_v.n_obs}")

    endo = all_v[all_v.obs["cell_type"].astype(str) == "endothelial cell"].copy()
    endo = to_symbols(endo)
    sc.pp.normalize_total(endo, target_sum=1e4)
    sc.pp.log1p(endo)

    contrasts = [
        ("ACA", "lateral_frontal", "ACA_vs_lateral_frontal"),
        ("ACA", "temporal", "ACA_vs_temporal"),
        ("lateral_frontal", "temporal", "lateral_frontal_vs_temporal_MCA_null"),
    ]
    results = []
    donor_rows = []
    for pos, neg, lab in contrasts:
        log(f"Contrast {lab}")
        r, ddf = run_contrast(endo, pos, neg, lab)
        if r:
            results.append(r)
            if ddf is not None and len(ddf):
                donor_rows.append(ddf)

    if donor_rows:
        pd.concat(donor_rows, ignore_index=True).to_csv(OUT_DONOR, index=False)

    # Verdict logic (LODO pooled cell AUROC; n_donors fixed at 3 in Siletti)
    def _ok(r):
        a = r.get("donor_lodo_pooled_cell_auroc")
        return a == a and a is not None and a >= 0.7

    territory = [r for r in results if r["contrast"].startswith("ACA_")]
    null = [r for r in results if "null" in r["contrast"]]
    terr_ok = any(_ok(r) for r in territory)
    null_low = all(
        (r.get("donor_lodo_pooled_cell_auroc") != r.get("donor_lodo_pooled_cell_auroc"))
        or r.get("donor_lodo_pooled_cell_auroc", 1) < 0.65
        for r in null
    ) if null else False
    n_donors = int(endo.obs["donor_id"].nunique())
    if n_donors <= 3 and terr_ok:
        verdict = "SUPPORTIVE_but_UNDERPOWERED_n_donors_3__slices_acquired"
    elif n_donors <= 3 and not terr_ok:
        verdict = "UNDERPOWERED_n_donors_3__slices_acquired"
    elif terr_ok and null_low:
        verdict = "SUPPORTIVE_within_Siletti_territory_separability"
    else:
        verdict = "INCONCLUSIVE_within_Siletti"

    summary = {
        "object": str(OUT_H5AD.relative_to(ROOT)),
        "n_vascular": int(all_v.n_obs),
        "n_endo": int(endo.n_obs),
        "n_donors": int(endo.obs["donor_id"].nunique()),
        "endo_by_arm": endo.obs["study_arm"].astype(str).value_counts().to_dict(),
        "endo_by_dissection": endo.obs["dissection_code"].astype(str).value_counts().to_dict(),
        "per_dissection": per,
        "tv_priors": tv_priors,
        "contrasts": results,
        "lateral_frontal_temporal_on_disk": True,
        "verdict": verdict,
        "locked_interpretation": [
            "Lateral frontal + temporal are MCA-dominant on TerritoryVoxels; ACA arm is cingulate pool.",
            "Frontal-vs-temporal is a same-territory (MCA) null within Siletti, not ACA vs MCA.",
            "Siletti adult census uses ~3 donors — donor OOF folds are few; do not overclaim.",
            "Not vessel-enriched; ambient neuronal means are reported per contrast.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    lines = [
        "# Linnarsson / Siletti within-study vascular (F/T acquired)",
        "",
        f"**Verdict:** `{verdict}`",
        "",
        "## Design",
        "- Same study: Siletti et al. HBCA CELLxGENE collection `283d65eb-…`",
        "- **ACA:** A24, A32, A25, A23 (A14 excluded)",
        "- **Lateral frontal (MCA):** MFG A46 + IFG A44–45",
        "- **Temporal (MCA):** MTG + STG",
        "- Metric: leave-one-**donor**-out **pooled cell-level** AUROC on Endo PCA",
        "- Note: the same 3 adult donors contribute cells to ACA, lateral frontal, and temporal — cross-cohort donor-mean OOF is invalid here.",
        "",
        "## Counts",
        f"- Vascular object: `{OUT_H5AD.relative_to(ROOT)}` — n={all_v.n_obs}",
        f"- Endo n={endo.n_obs}; donors={summary['n_donors']}",
        f"- Endo by arm: `{summary['endo_by_arm']}`",
        f"- Endo by dissection: `{summary['endo_by_dissection']}`",
        "",
        "## TerritoryVoxels priors (dominant)",
    ]
    for code, tv in tv_priors.items():
        lines.append(
            f"- {code}: dominant={tv.get('tv_dominant')} "
            f"pACA={tv.get('tv_p_ACA'):.2f} pMCA={tv.get('tv_p_MCA'):.2f} "
            f"assign={tv.get('tv_assignment_fraction'):.2f}"
        )
    lines += [
        "",
        "## Contrasts (LODO pooled cell AUROC)",
        "| contrast | n_pos | n_neg | donors | pooled AUROC | mean fold AUROC | ambient pos/neg |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in results:
        auc = r.get("donor_lodo_pooled_cell_auroc")
        mf = r.get("donor_lodo_mean_fold_auroc")
        auc_s = f"{auc:.3f}" if auc == auc else "NA"
        mf_s = f"{mf:.3f}" if mf == mf else "NA"
        lines.append(
            f"| {r['contrast']} | {r['n_pos']} | {r['n_neg']} | "
            f"{r['n_donors_pos']}+{r['n_donors_neg']} | {auc_s} | {mf_s} | "
            f"{r.get('ambient_neuronal_mean_pos'):.3f}/{r.get('ambient_neuronal_mean_neg'):.3f} |"
        )
    lines += [
        "",
        "## Interpretation locks",
        "- Do **not** treat lateral frontal as ACA (TerritoryVoxels MCA-dominant).",
        "- F vs T = MCA vs MCA null (regional, not territorial).",
        "- Only 3 adult donors in the Siletti census → underpowered as sole H1 proof even if AUROC is high.",
        "- Class B still requires external ACA arms; enrichment-matched DE needs vessel-enriched ACC (no public download).",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    log(json.dumps({"verdict": verdict, "contrasts": results}, indent=2, default=str))


if __name__ == "__main__":
    main()
