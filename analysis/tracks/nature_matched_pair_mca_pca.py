"""
Nature-standard rigor: matched Midtemporal (MCA) vs Entorhinal (PCA)
donor-aware endothelial and mural contrasts with bootstrap CIs.

Locks:
  - Donor-level inference only for primary claims
  - No soft DE marker gene lists
  - Ambient residualization for endo modules
  - Report regardless of significance; power note if under MDE
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"F:\Cerebrovascular SCS")
H5AD = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"
OUT = ROOT / "results/tracks/nature_rigor"
OUT.mkdir(parents=True, exist_ok=True)

MCA_REG = "Midtemporal_cortex"
PCA_REG = "Entorhinal_cortex"
PAIR = (MCA_REG, PCA_REG)

MODULES = {
    "BBB_tight_junction": ["CLDN5", "OCLN", "TJP1", "TJP2", "CDH5", "ESAM", "PECAM1"],
    "angiogenesis": ["VEGFA", "FLT1", "KDR", "ENG", "VWF", "ANGPT2", "TEK"],
    "ECM_basement": ["COL4A1", "COL4A2", "LAMA4", "LAMB1", "NID1", "HSPG2"],
    "hypoxia_stress": ["HIF1A", "VEGFA", "SLC2A1", "BNIP3", "DDIT4"],
    "inflammation": ["ICAM1", "VCAM1", "SELE", "CCL2", "IL6", "NFKB1"],
    "flow_shear": ["KLF2", "KLF4", "NOS3", "THBD", "PLAT", "MEF2C"],
    "mural_contractile": ["ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "PDGFRB", "RGS5"],
}
AMBIENT = ["RBFOX1", "SNAP25", "SYT1", "SLC17A7", "GAD1", "MBP", "PLP1", "MOBP", "AQP4", "GFAP"]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def log1p_norm(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.max() > 50:
        lib = X.sum(axis=1, keepdims=True)
        lib[lib == 0] = 1.0
        X = np.log1p(X / lib * 1e4)
    return X


def module_mat(X: np.ndarray, var_names: list[str], genes: list[str]) -> np.ndarray:
    idx = [var_names.index(g) for g in genes if g in var_names]
    if not idx:
        return np.full(X.shape[0], np.nan)
    return X[:, idx].mean(axis=1)


def bootstrap_mean_diff(a: np.ndarray, b: np.ndarray, B: int = 2000, seed: int = 0):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5:
        return {
            "n_a": len(a),
            "n_b": len(b),
            "mean_a": float(np.nanmean(a)) if len(a) else np.nan,
            "mean_b": float(np.nanmean(b)) if len(b) else np.nan,
            "delta_b_minus_a": np.nan,
            "ci95": [np.nan, np.nan],
            "cohens_d": np.nan,
            "mannwhitney_p": np.nan,
        }
    delta = float(b.mean() - a.mean())
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(B):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        boots.append(bb.mean() - aa.mean())
    boots = np.asarray(boots)
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / max(len(a) + len(b) - 2, 1))
    d = float(delta / pooled) if pooled > 0 else np.nan
    st = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
        "delta_b_minus_a": delta,
        "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
        "cohens_d": d,
        "mannwhitney_p": float(st.pvalue),
    }


def donor_means(scores: pd.DataFrame, score_col: str, celltype: str) -> pd.DataFrame:
    sub = scores[scores["celltype"] == celltype].copy()
    return (
        sub.groupby(["donor", "region"], observed=True)[score_col]
        .mean()
        .reset_index()
        .pivot(index="donor", columns="region", values=score_col)
    )


def main():
    log("Loading Garcia ROSMAP vascular")
    a = ad.read_h5ad(H5AD)
    obs = a.obs.copy()
    obs["donor"] = obs["subject"].astype(str)
    obs["region"] = obs["brain_region"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)

    # Restrict to pair regions
    mask = obs["region"].isin(PAIR).values
    idx = np.where(mask)[0]
    X = a.X[idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = log1p_norm(np.asarray(X, dtype=float))
    var_names = list(a.var_names)
    meta = obs.iloc[idx].copy()

    # Ambient score + residualize modules later
    amb = module_mat(X, var_names, AMBIENT)
    meta["ambient"] = amb

    rows = []
    module_cols = {}
    for name, genes in MODULES.items():
        s = module_mat(X, var_names, genes)
        # residualize on ambient (linear)
        m = np.isfinite(s) & np.isfinite(amb)
        resid = s.copy()
        if m.sum() > 10 and np.nanstd(amb[m]) > 0:
            b1 = np.polyfit(amb[m], s[m], 1)
            resid = s - (b1[0] * amb + b1[1])
        meta[name] = s
        meta[f"{name}_resid"] = resid
        module_cols[name] = f"{name}_resid"

    # Cell counts by donor × region × type
    counts = (
        meta.groupby(["donor", "region", "celltype"], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    counts.to_csv(OUT / "pair_cell_counts.csv", index=False)

    results = {
        "pair": {"MCA_region": MCA_REG, "PCA_region": PCA_REG},
        "n_cells_pair": int(len(meta)),
        "n_donors_any": int(meta["donor"].nunique()),
        "design": "donor-mean module scores; Midtemporal(MCA) vs Entorhinal(PCA); ambient-residualized",
        "MDE": {"cohens_d": 0.5, "note": "from Nature amendment v2"},
        "tests": {},
    }

    for ctype in ["Endo", "Per", "SMC"]:
        results["tests"][ctype] = {}
        for mod, col in module_cols.items():
            if ctype == "Endo" and mod == "mural_contractile":
                continue
            if ctype in ("Per", "SMC") and mod in ("BBB_tight_junction", "flow_shear"):
                # still compute but lower priority
                pass
            piv = donor_means(meta, col, ctype)
            if MCA_REG not in piv.columns or PCA_REG not in piv.columns:
                results["tests"][ctype][mod] = {"error": "missing region columns"}
                continue
            # paired donors with both regions
            both = piv[[MCA_REG, PCA_REG]].dropna()
            if len(both) >= 5:
                d = both[PCA_REG].values - both[MCA_REG].values
                w = stats.wilcoxon(both[PCA_REG], both[MCA_REG])
                rng = np.random.default_rng(1)
                boots = [float(rng.choice(d, size=len(d), replace=True).mean()) for _ in range(2000)]
                pooled = float(np.std(d, ddof=1)) if len(d) > 1 else np.nan
                cohend = float(np.mean(d) / pooled) if pooled and pooled > 0 else np.nan
                rec = {
                    "design": "paired_donors_both_regions",
                    "n_paired_donors": int(len(both)),
                    "mean_MCA": float(both[MCA_REG].mean()),
                    "mean_PCA": float(both[PCA_REG].mean()),
                    "delta_PCA_minus_MCA": float(np.mean(d)),
                    "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                    "cohens_d_paired": cohend,
                    "wilcoxon_p": float(w.pvalue),
                    "meets_MDE_d05": bool(abs(cohend) >= 0.5) if np.isfinite(cohend) else False,
                    "ci_excludes_0": bool(np.percentile(boots, 2.5) > 0 or np.percentile(boots, 97.5) < 0),
                }
            else:
                # unpaired among region-present donors
                a_vals = piv[MCA_REG].dropna().values
                b_vals = piv[PCA_REG].dropna().values
                rec = bootstrap_mean_diff(a_vals, b_vals)
                rec["design"] = "unpaired_region_present_donors"
                rec["meets_MDE_d05"] = bool(abs(rec.get("cohens_d") or 0) >= 0.5)
                ci = rec.get("ci95", [np.nan, np.nan])
                rec["ci_excludes_0"] = bool(
                    np.isfinite(ci[0]) and np.isfinite(ci[1]) and (ci[0] > 0 or ci[1] < 0)
                )
            # Nature label
            if rec.get("n_paired_donors", rec.get("n_a", 0)) < 10:
                rec["reporting_label"] = "EXPLORATORY_UNDERPOWERED"
            elif rec.get("meets_MDE_d05") and rec.get("ci_excludes_0") and rec.get("wilcoxon_p", rec.get("mannwhitney_p", 1)) < 0.05:
                rec["reporting_label"] = "PASSES_PRESPEC_MDE_EXPLORATORY_UNTIL_REPLICATION"
            else:
                rec["reporting_label"] = "NULL_OR_BELOW_MDE"
            results["tests"][ctype][mod] = rec

    # Separability LODO AUROC endo MCA vs PCA pair (donor-aware)
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    endo = meta[meta["celltype"] == "Endo"].copy()
    feat_cols = [module_cols[m] for m in MODULES if m != "mural_contractile"]
    # donor means
    dm = (
        endo.groupby(["donor", "region"], observed=True)[feat_cols]
        .mean()
        .reset_index()
    )
    dm = dm[dm["region"].isin(PAIR)]
    # need donors with both for paired; for AUROC use all donor-region rows
    y = (dm["region"] == PCA_REG).astype(int).values
    Xd = dm[feat_cols].values
    donors = dm["donor"].values
    aucs = []
    uniq = np.unique(donors)
    for d in uniq:
        te = donors == d
        tr = ~te
        if len(np.unique(y[tr])) < 2 or te.sum() == 0:
            continue
        sc = StandardScaler()
        Xtr = sc.fit_transform(Xd[tr])
        Xte = sc.transform(Xd[te])
        clf = LogisticRegression(max_iter=2000)
        clf.fit(Xtr, y[tr])
        if len(np.unique(y[te])) < 2:
            # single region for held-out donor — skip AUROC
            continue
        aucs.append(float(roc_auc_score(y[te], clf.predict_proba(Xte)[:, 1])))
    # Cell-pooled LODO as secondary only
    # Donor-mean LODO with leave-one-donor: if donor has one region only, score that point
    aucs2 = []
    for d in uniq:
        te = donors == d
        tr = ~te
        if len(np.unique(y[tr])) < 2:
            continue
        sc = StandardScaler()
        Xtr = sc.fit_transform(Xd[tr])
        Xte = sc.transform(Xd[te])
        clf = LogisticRegression(max_iter=2000).fit(Xtr, y[tr])
        proba = clf.predict_proba(Xte)[:, 1]
        # accuracy-style: mean predicted prob aligned with label
        aucs2.append(float(np.mean(np.abs(proba - y[te]) < 0.5)))  # balanced accuracy proxy per donor
    # Proper: train on donors, evaluate all held-out donor-region points' ranking needs ≥2 classes in test — use pooled OOF scores
    oof = np.full(len(y), np.nan)
    for d in uniq:
        te = donors == d
        tr = ~te
        if len(np.unique(y[tr])) < 2 or te.sum() == 0:
            continue
        sc = StandardScaler()
        clf = LogisticRegression(max_iter=2000)
        clf.fit(sc.fit_transform(Xd[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(Xd[te]))[:, 1]
    m = np.isfinite(oof)
    if m.sum() >= 10 and len(np.unique(y[m])) == 2:
        auc = float(roc_auc_score(y[m], oof[m]))
        # bootstrap donors
        rng = np.random.default_rng(2)
        ba = []
        for _ in range(1000):
            take = rng.choice(uniq, size=len(uniq), replace=True)
            idx = np.concatenate([np.where(donors == dd)[0] for dd in take])
            if len(np.unique(y[idx])) < 2:
                continue
            oo = np.full(len(idx), np.nan)
            # recompute quickly with same global oof resampled
            ba.append(float(roc_auc_score(y[idx], oof[idx])))
        results["endo_donor_region_OOF_AUROC"] = {
            "auroc": auc,
            "ci95": [float(np.percentile(ba, 2.5)), float(np.percentile(ba, 97.5))] if ba else [np.nan, np.nan],
            "n_donor_region_rows": int(m.sum()),
            "n_donors": int(len(uniq)),
            "label": "EXPLORATORY_MATCHED_PAIR_SEPARABILITY",
        }
    else:
        results["endo_donor_region_OOF_AUROC"] = {"auroc": np.nan, "note": "insufficient class balance"}

    with open(OUT / "nature_matched_pair_summary.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    # Flat table
    flat = []
    for ctype, mods in results["tests"].items():
        for mod, rec in mods.items():
            if not isinstance(rec, dict) or "error" in rec:
                continue
            flat.append({"celltype": ctype, "module": mod, **{k: v for k, v in rec.items() if k != "ci95"}, "ci95_lo": rec.get("ci95", [None, None])[0], "ci95_hi": rec.get("ci95", [None, None])[1]})
    pd.DataFrame(flat).to_csv(OUT / "nature_matched_pair_tests.csv", index=False)

    _report(results)
    log(f"Wrote {OUT}")


def _report(r: dict):
    lines = [
        "# Nature rigor — Midtemporal (MCA) vs Entorhinal (PCA)",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d')}",
        f"**Design:** {r['design']}",
        f"**n cells (pair):** {r['n_cells_pair']}; donors: {r['n_donors_any']}",
        "",
        "## Endo donor-region OOF AUROC",
        f"`{json.dumps(r.get('endo_donor_region_OOF_AUROC', {}))}`",
        "",
        "## Module tests (ambient-residualized; report all)",
        "",
    ]
    for ctype, mods in r["tests"].items():
        lines.append(f"### {ctype}")
        for mod, rec in mods.items():
            if not isinstance(rec, dict):
                continue
            lab = rec.get("reporting_label", "")
            lines.append(f"- **{mod}**: `{lab}` — {json.dumps({k: rec[k] for k in rec if k != 'reporting_label'})}")
        lines.append("")
    lines += [
        "## Nature locks",
        "- No marker-gene DE lists.",
        "- Labels EXPLORATORY until enrichment-matched replication.",
        "- MDE: |d|≥0.5 and CI excludes 0 for ‘passes pre-spec’ exploratory flag.",
        "",
    ]
    (OUT / "nature_matched_pair_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
