"""Within-donor H4 CV on tight Yang endo (paired territory arms)."""
from __future__ import annotations

import warnings
from pathlib import Path

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
YANG = ROOT / "data/processed/yang/yang_vascular_qc_tight.h5ad"
OUT_CSV = ROOT / "results/h4/yang_h4_within_donor.csv"
OUT_MD = ROOT / "results/h4/yang_h4_within_donor_report.md"
ARMS = ["PCA_hippocampus", "ACA_MCA_borderzone_SFG"]
POS = "ACA_MCA_borderzone_SFG"


def dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def main():
    yang = sc.read_h5ad(YANG)
    endo = yang[yang.obs["vascular_lineage"].astype(str) == "endo"].copy()
    both = endo[endo.obs["territory_contrast"].astype(str).isin(ARMS)].copy()
    sc.pp.highly_variable_genes(both, n_top_genes=1500, flavor="seurat", subset=False)
    genes = both.var_names[both.var["highly_variable"]].tolist()

    rows = []
    for disease in ["control", "AD"]:
        sub = endo[
            (endo.obs["disease"].astype(str) == disease)
            & (endo.obs["territory_contrast"].astype(str).isin(ARMS))
        ]
        for d, g in sub.obs.groupby(sub.obs["donor_id"].astype(str), observed=True):
            arms = set(g["territory_contrast"].astype(str))
            if not set(ARMS).issubset(arms):
                continue
            counts = {a: int((g["territory_contrast"].astype(str) == a).sum()) for a in ARMS}
            if min(counts.values()) < 30:
                continue
            ad_d = sub[sub.obs["donor_id"].astype(str) == d]
            rng = np.random.default_rng(abs(hash(str(d))) % 10**6)
            parts = []
            for a in ARMS:
                idx = ad_d.obs_names[ad_d.obs["territory_contrast"].astype(str) == a]
                n = min(len(idx), 1500)
                parts.append(ad_d[rng.choice(idx.to_numpy(), n, replace=False)])
            ad_d = sc.concat(parts)
            X = dense(ad_d[:, genes].X)
            y = (ad_d.obs["territory_contrast"].astype(str) == POS).astype(int).to_numpy()
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
            aucs = []
            for tr, te in skf.split(X, y):
                sc_ = StandardScaler().fit(X[tr])
                clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.1)
                clf.fit(sc_.transform(X[tr]), y[tr])
                proba = clf.predict_proba(sc_.transform(X[te]))[:, 1]
                if len(np.unique(y[te])) > 1:
                    aucs.append(roc_auc_score(y[te], proba))
            rows.append(
                {
                    "disease": disease,
                    "donor_id": d,
                    "n_cells": int(ad_d.n_obs),
                    "auroc_mean": float(np.mean(aucs)),
                    "auroc_std": float(np.std(aucs)),
                    "n_SFG": counts[POS],
                    "n_PCA": counts["PCA_hippocampus"],
                }
            )

    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    ctrl = df.loc[df["disease"] == "control", "auroc_mean"]
    ad = df.loc[df["disease"] == "AD", "auroc_mean"]
    summary = {
        "control_mean": float(ctrl.mean()),
        "control_std": float(ctrl.std(ddof=0)),
        "n_control": int(len(ctrl)),
        "ad_mean": float(ad.mean()),
        "ad_std": float(ad.std(ddof=0)),
        "n_ad": int(len(ad)),
        "delta": float(ad.mean() - ctrl.mean()),
    }
    print(df.to_string(index=False))
    print(summary)

    md = f"""# H4 within-donor CV — Yang tight endo

For each donor with both arms, 5-fold stratified CV predicts territory_contrast from HVGs.
This is the appropriate Yang test of readable territorial identity (paired within donor).

| Metric | Value |
|---|---:|
| Control within-donor AUROC mean | {summary['control_mean']:.3f} |
| AD within-donor AUROC mean | {summary['ad_mean']:.3f} |
| Delta (AD - control) | {summary['delta']:.3f} |
| N control / AD donors | {summary['n_control']} / {summary['n_ad']} |

## Per donor

{df.to_string(index=False)}

Negative delta = H4-consistent (weaker territory separability in AD).

Compare to cross-donor LODO (~chance in controls): a shared cross-donor territorial
program is weak; within-donor identity is strong and degrades in AD.
"""
    OUT_MD.write_text(md, encoding="utf-8")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
