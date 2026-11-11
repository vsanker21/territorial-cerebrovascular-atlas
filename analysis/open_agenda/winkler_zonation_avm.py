"""Winkler UCSC endothelial: A-V zonation programs; AVM vs CTRL exploratory (no territory)."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "data/raw/winkler/ucsc/endothelial/meta.tsv"
EXPR = ROOT / "data/raw/winkler/ucsc/endothelial/exprMatrix.tsv.gz"
AM_META = ROOT / "data/raw/winkler/ucsc/am-endothelial/meta.tsv"
AM_EXPR = ROOT / "data/raw/winkler/ucsc/am-endothelial/exprMatrix.tsv.gz"
OUT_DIR = ROOT / "results/h1"
OUT_JSON = OUT_DIR / "winkler_zonation_avm_summary.json"
OUT_MD = OUT_DIR / "winkler_zonation_avm_report.md"


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def load_ucsc(meta_path: Path, expr_path: Path, max_genes: int = 2000):
    meta = pd.read_csv(meta_path, sep="\t")
    if "Cell" in meta.columns:
        meta = meta.set_index("Cell")
    log(f"Loading matrix {expr_path.name} (may take a minute)")
    # UCSC exprMatrix: genes x cells, first col gene
    # Read in chunks for memory — sample HV-ish by variance on subset
    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    # columns may be cells
    common = [c for c in df.columns if c in meta.index]
    if len(common) < 100:
        # try without suffix mismatches
        meta.index = meta.index.astype(str)
        df.columns = df.columns.astype(str)
        common = [c for c in df.columns if c in meta.index]
    df = df[common]
    meta = meta.loc[common]
    # pick top variance genes
    v = df.var(axis=1).to_numpy()
    top = np.argsort(v)[::-1][:max_genes]
    df = df.iloc[top]
    X = df.T.to_numpy(dtype=np.float32)
    return meta, X, df.index.to_list()


def cell_cv(X, y):
    y = np.asarray(y)
    if len(np.unique(y)) < 2:
        return np.nan, np.nan
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs = []
    Xs = StandardScaler().fit_transform(X)
    for tr, te in skf.split(Xs, y):
        if len(np.unique(y[tr])) < 2:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
        clf.fit(Xs[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(Xs[te])[:, 1]))
    return (float(np.mean(aucs)), float(np.std(aucs))) if aucs else (np.nan, np.nan)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {
        "region_metadata": False,
        "territory_test": "NO-GO",
        "note": "Winkler UCSC has A-V clusters but no brain-region/TerritoryVoxels labels.",
    }

    # Control endothelial zonation
    meta, X, genes = load_ucsc(META, EXPR)
    cl = meta["Cluster"].astype(str) if "Cluster" in meta.columns else meta["clusters"].astype(str)
    # Cap vs Arterial
    y_cap = cl.str.contains("Cap").astype(int).to_numpy()
    y_art = cl.str.startswith("Art").astype(int).to_numpy()
    # Cap vs Art
    keep = (y_cap == 1) | (y_art == 1)
    auc_z, sd_z = cell_cv(X[keep], y_cap[keep])
    results["ctrl_endothelial"] = {
        "n": int(len(meta)),
        "cluster_counts": cl.value_counts().to_dict(),
        "cap_vs_arterial_cell_cv_auroc": auc_z,
        "cap_vs_arterial_cell_cv_std": sd_z,
        "samples": meta["Sample"].astype(str).value_counts().to_dict() if "Sample" in meta else {},
    }
    log(f"CTRL Cap vs Art AUROC={auc_z}")

    # AVM endothelial if present
    if AM_META.exists() and AM_EXPR.exists():
        am_meta, am_X, _ = load_ucsc(AM_META, AM_EXPR)
        # infer AVM vs CTRL from Sample
        samp = am_meta["Sample"].astype(str) if "Sample" in am_meta.columns else pd.Series(["unk"] * len(am_meta))
        is_avm = samp.str.contains("AVM|avm|nidus", case=False, regex=True) | (
            am_meta.get("AVM", pd.Series([0] * len(am_meta))).astype(str).isin(["1", "True", "true", "AVM"])
        )
        # also Cluster Nidus
        cl2 = am_meta["Cluster"].astype(str) if "Cluster" in am_meta.columns else am_meta.get("clusters", pd.Series(["x"] * len(am_meta))).astype(str)
        if is_avm.sum() < 50:
            # try Sample naming
            is_avm = ~samp.str.upper().str.startswith("CTRL")
        if is_avm.sum() >= 50 and (~is_avm).sum() >= 50:
            auc_d, sd_d = cell_cv(am_X, is_avm.astype(int).to_numpy())
            results["am_endothelial_avm_vs_ctrl"] = {
                "n": int(len(am_meta)),
                "n_avm": int(is_avm.sum()),
                "n_ctrl": int((~is_avm).sum()),
                "cell_cv_auroc": auc_d,
                "cell_cv_std": sd_d,
                "cluster_counts": cl2.value_counts().to_dict(),
                "claim": "EXPLORATORY_disease_vs_control_NOT_territory",
            }
            log(f"AVM vs CTRL AUROC={auc_d}")
        else:
            results["am_endothelial_avm_vs_ctrl"] = {
                "status": "could_not_parse_AVM_labels",
                "sample_top": samp.value_counts().head(15).to_dict(),
            }

    results["verdict"] = "ZONATION_PROGRAMS_PRESENT_TERRITORY_NO_GO"
    OUT_JSON.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    md = [
        "# Winkler UCSC — zonation vs territory",
        "",
        f"**Verdict:** `{results['verdict']}`",
        "",
        "- A–V-like EC clusters present (Art/Cap/Venous/Venule).",
        "- **No region / TerritoryVoxels metadata** → cannot test residual territory signal within Winkler.",
        "- Use as zonation covariate library / AVM disease exploratory only.",
        "",
        f"CTRL Cap vs Arterial cell CV AUROC: **{auc_z:.3f}**",
    ]
    if "am_endothelial_avm_vs_ctrl" in results and "cell_cv_auroc" in results["am_endothelial_avm_vs_ctrl"]:
        md.append(
            f"AVM vs CTRL (am-endothelial) cell CV AUROC: **{results['am_endothelial_avm_vs_ctrl']['cell_cv_auroc']:.3f}** (exploratory H4-adjacent, not territory)."
        )
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    log(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
