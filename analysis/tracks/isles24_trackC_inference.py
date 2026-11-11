#!/usr/bin/env python3
"""ISLES'24 Track C inference: bootstrap ΔR², nested CV, multi-outcome.

Uses existing subject table (TopCoW + LVO vol + lesion shares).
Exploratory until larger n; not confirmatory Nature claim alone.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ROOT = Path(r"F:\Cerebrovascular SCS")
OUT = ROOT / "results/nature_gates/isles24_trackC"
TABLE = OUT / "isles24_trackC_subject_table.csv"
OUT.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def design_matrix(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    parts = []
    for c in cols:
        if c == "lvo_volume_log":
            parts.append(np.log1p(pd.to_numeric(df["lvo_volume_mm3"], errors="coerce").fillna(0).values))
        elif c == "cow_vol_MCA_log":
            parts.append(np.log1p(pd.to_numeric(df["cow_vol_MCA"], errors="coerce").fillna(0).values))
        else:
            parts.append(pd.to_numeric(df[c], errors="coerce").values)
    return np.column_stack(parts)


SPECS = {
    "age": ["age"],
    "age_LVOvol": ["age", "lvo_volume_log"],
    "age_CoWcomplete": ["age", "cow_completeness"],
    "age_CoWMCA": ["age", "cow_vol_MCA_log"],
    "age_angio": ["age", "lvo_volume_log", "cow_completeness", "cow_vol_MCA_log"],
    "age_full": [
        "age",
        "lvo_volume_log",
        "cow_completeness",
        "cow_vol_MCA_log",
        "cow_laterality_MCA",
    ],
}


def ols_r2(y: np.ndarray, X: np.ndarray) -> float:
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    if len(y) < max(10, X.shape[1] + 3):
        return np.nan
    Xd = np.column_stack([np.ones(len(y)), X])
    coef, *_ = lstsq(Xd, y, rcond=None)
    pred = Xd @ coef
    sse = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2)) + 1e-12
    return float(1 - sse / sst)


def cv_r2(y: np.ndarray, X: np.ndarray, k: int = 5, seed: int = 0) -> float:
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    n = len(y)
    if n < max(15, X.shape[1] + 5):
        return np.nan
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    folds = np.array_split(idx, k)
    preds = np.full(n, np.nan)
    for i in range(k):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(k) if j != i])
        if len(np.unique(y[tr])) < 1 or len(tr) <= X.shape[1] + 1:
            continue
        Xtr = np.column_stack([np.ones(len(tr)), X[tr]])
        Xte = np.column_stack([np.ones(len(te)), X[te]])
        coef, *_ = lstsq(Xtr, y[tr], rcond=None)
        preds[te] = Xte @ coef
    m2 = np.isfinite(preds)
    if m2.sum() < 10:
        return np.nan
    sse = float(np.sum((y[m2] - preds[m2]) ** 2))
    sst = float(np.sum((y[m2] - y[m2].mean()) ** 2)) + 1e-12
    return float(1 - sse / sst)


def bootstrap_delta(
    y: np.ndarray,
    X_base: np.ndarray,
    X_full: np.ndarray,
    B: int = 2000,
    seed: int = 1,
) -> dict:
    m = np.isfinite(y) & np.all(np.isfinite(X_base), axis=1) & np.all(np.isfinite(X_full), axis=1)
    y, X_base, X_full = y[m], X_base[m], X_full[m]
    n = len(y)
    if n < 12:
        return {"n": n, "delta": np.nan, "ci95": [np.nan, np.nan], "p_boot": np.nan}
    r0 = ols_r2(y, X_base)
    r1 = ols_r2(y, X_full)
    delta = float(r1 - r0) if np.isfinite(r0) and np.isfinite(r1) else np.nan
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(B):
        take = rng.choice(n, size=n, replace=True)
        d = ols_r2(y[take], X_full[take]) - ols_r2(y[take], X_base[take])
        if np.isfinite(d):
            boots.append(d)
    boots = np.asarray(boots)
    if len(boots) < 100:
        return {"n": n, "delta": delta, "ci95": [np.nan, np.nan], "p_boot": np.nan}
    # one-sided: P(delta <= 0)
    p = float(np.mean(boots <= 0))
    return {
        "n": int(n),
        "R2_base": float(r0),
        "R2_full": float(r1),
        "delta": delta,
        "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
        "p_boot_delta_le_0": p,
        "ci_excludes_0": bool(np.percentile(boots, 2.5) > 0),
    }


def spearman_block(df: pd.DataFrame, ycol: str, xs: list[str]) -> list[dict]:
    from scipy import stats

    rows = []
    y = pd.to_numeric(df[ycol], errors="coerce")
    for x in xs:
        if x.endswith("_log"):
            # handled via design; skip raw
            continue
        if x not in df.columns and x not in ("lvo_volume_log", "cow_vol_MCA_log"):
            continue
        if x == "lvo_volume_log":
            xv = np.log1p(pd.to_numeric(df["lvo_volume_mm3"], errors="coerce"))
            name = "lvo_volume_log"
        elif x == "cow_vol_MCA_log":
            xv = np.log1p(pd.to_numeric(df["cow_vol_MCA"], errors="coerce"))
            name = "cow_vol_MCA_log"
        else:
            xv = pd.to_numeric(df[x], errors="coerce")
            name = x
        m = y.notna() & np.isfinite(xv)
        if m.sum() < 10:
            continue
        rho, p = stats.spearmanr(y[m], xv[m])
        rows.append({"outcome": ycol, "feature": name, "spearman_rho": float(rho), "p": float(p), "n": int(m.sum())})
    return rows


def main() -> None:
    df = pd.read_csv(TABLE)
    if "qc_pass" in df.columns:
        sub = df[(df["status"] == "ok") & (df["qc_pass"] == 1)].copy()
    else:
        sub = df[df["status"] == "ok"].copy()
    log(f"Inference n={len(sub)}")

    outcomes = [
        c
        for c in [
            "lesion_share_MCA",
            "lesion_share_ACA",
            "lesion_share_PCA",
            "lesion_share_VB",
            "lesion_vox_total",
            "nihss_admission",
            "mrs_3month",
        ]
        if c in sub.columns
    ]

    model_rows = []
    delta_rows = []
    spearman_rows = []

    for ycol in outcomes:
        y = pd.to_numeric(sub[ycol], errors="coerce").values
        # log for voxel counts
        if ycol == "lesion_vox_total":
            y = np.log1p(y)
            yname = "log1p_lesion_vox_total"
        else:
            yname = ycol

        for name, cols in SPECS.items():
            X = design_matrix(sub, cols)
            r2 = ols_r2(y, X)
            r2cv = cv_r2(y, X, k=5, seed=hash(yname) % 10000)
            model_rows.append(
                {
                    "outcome": yname,
                    "model": name,
                    "R2_in": r2,
                    "R2_cv5": r2cv,
                    "n": int(np.isfinite(y).sum()),
                    "p_features": len(cols),
                }
            )

        X0 = design_matrix(sub, SPECS["age"])
        for full_name in ("age_LVOvol", "age_angio", "age_full"):
            X1 = design_matrix(sub, SPECS[full_name])
            boot = bootstrap_delta(y, X0, X1, B=2000, seed=7)
            delta_rows.append({"outcome": yname, "contrast": f"{full_name}_minus_age", **boot})

        spearman_rows.extend(
            spearman_block(
                sub.assign(**{yname: y}) if yname != ycol else sub,
                ycol if yname == ycol else ycol,  # spearman on original scale for shares
                ["age", "lvo_volume_mm3", "cow_completeness", "cow_vol_MCA", "cow_laterality_MCA", "nihss_admission"],
            )
        )

    # Fix spearman for lesion_vox on log — recompute cleanly
    spearman_rows = []
    for ycol in outcomes:
        spearman_rows.extend(
            spearman_block(
                sub,
                ycol,
                ["age", "lvo_volume_mm3", "cow_completeness", "cow_vol_MCA", "cow_laterality_MCA"],
            )
        )

    pd.DataFrame(model_rows).to_csv(OUT / "isles24_trackC_inference_models.csv", index=False)
    pd.DataFrame(delta_rows).to_csv(OUT / "isles24_trackC_inference_deltaR2.csv", index=False)
    pd.DataFrame(spearman_rows).to_csv(OUT / "isles24_trackC_inference_spearman.csv", index=False)

    # Primary Nature-relevant: MCA share angio delta
    primary = [r for r in delta_rows if r["outcome"] == "lesion_share_MCA" and "age_full" in r["contrast"]]
    primary = primary[0] if primary else {}
    mca_models = [r for r in model_rows if r["outcome"] == "lesion_share_MCA"]

    verdict = "EXPLORATORY_INFERENCE"
    if primary.get("ci_excludes_0") and primary.get("delta", 0) > 0:
        verdict = "EXPLORATORY_ANGIO_DELTA_CI_EXCLUDES_0__NOT_CONFIRMATORY_ALONE"

    status = {
        "date": time.strftime("%Y-%m-%d"),
        "n": int(len(sub)),
        "primary_outcome": "lesion_share_MCA",
        "primary_delta_age_full_minus_age": primary,
        "mca_models": mca_models,
        "verdict": verdict,
        "locks": [
            "ArterialAtlas_level2 not TV_BMM",
            "Subject-level CoW/LVO not atlas-constant MRA",
            "Not confirmatory Nature Track C alone (n, registration QC audit, replication)",
            "TRCI/BCI excluded",
        ],
    }
    (OUT / "isles24_trackC_inference_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")

    lines = [
        "# ISLES'24 Track C — formal inference (exploratory)",
        "",
        f"**Date:** {status['date']}",
        f"**n (QC-pass):** {status['n']}",
        f"**Verdict:** `{verdict}`",
        "",
        "## Primary: MCA lesion share — ΔR² (full angio − age)",
        f"`{json.dumps(primary)}`",
        "",
        "## MCA share model R² (in-sample / 5-fold CV)",
    ]
    for r in mca_models:
        lines.append(f"- **{r['model']}**: R²_in={r['R2_in']:.3f}; R²_cv5={r['R2_cv5']:.3f}")
    lines += [
        "",
        "## Files",
        "- `isles24_trackC_inference_models.csv`",
        "- `isles24_trackC_inference_deltaR2.csv`",
        "- `isles24_trackC_inference_spearman.csv`",
        "",
        "## Interpretation",
        "Positive angio ΔR² with CI excluding 0 supports Track C architecture mediation as an "
        "**exploratory** mechanism signal. Confirmatory package still needs larger n, registration "
        "audit, and (for Nature) linkage to molecular gates.",
        "",
    ]
    (OUT / "isles24_trackC_inference_report.md").write_text("\n".join(lines), encoding="utf-8")
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
