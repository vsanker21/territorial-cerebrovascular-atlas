"""Moonshot #4 — CoW topology as instrument for territorial lesion share.

Nested CV + split conformal. Only OOS metrics count.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from _paths import HONESTY, ISLES24_TABLE, RES

OUT = RES / "04_cow_iv_conformal"
OUT.mkdir(parents=True, exist_ok=True)


def _conformal_intervals(y_cal, pred_cal, pred_test, alpha=0.1):
    """Split conformal absolute residual intervals."""
    resid = np.abs(y_cal - pred_cal)
    q = np.quantile(resid, 1 - alpha, method="higher")
    return pred_test - q, pred_test + q, float(q)


def run() -> dict:
    df = pd.read_csv(ISLES24_TABLE)
    df = df[df["qc_pass"] == 1].copy()
    results = {"honesty_lock": HONESTY["4"], "n_qc": int(len(df))}

    feature_cols = [
        c
        for c in df.columns
        if c.startswith("cow_present_")
        or c in [
            "cow_completeness",
            "cow_n_labels",
            "cow_laterality_MCA",
            "cow_vol_ACA",
            "cow_vol_MCA",
            "cow_vol_PCA",
            "cow_vol_Acom",
            "cow_vol_Pcom",
            "lvo_present",
            "lvo_volume_mm3",
            "age",
        ]
    ]
    # Derived instruments
    df["missing_Acom"] = 1 - df.get("cow_present_Acom", 0)
    df["missing_Pcom_any"] = 1 - df.get("cow_present_Pcom", ((df.get("cow_present_R_Pcom", 0) + df.get("cow_present_L_Pcom", 0)) > 0).astype(int))
    df["aca_asym"] = (df["cow_vol_R_ACA"] - df["cow_vol_L_ACA"]).abs() / (df["cow_vol_ACA"] + 1e-6)
    feature_cols += ["missing_Acom", "missing_Pcom_any", "aca_asym"]
    feature_cols = [c for c in feature_cols if c in df.columns]

    outcomes = {
        "lesion_share_MCA": "lesion_share_MCA",
        "lesion_share_ACA": "lesion_share_ACA",
        "mca_minus_aca_share": None,
    }
    df["mca_minus_aca_share"] = df["lesion_share_MCA"] - df["lesion_share_ACA"]
    outcomes["mca_minus_aca_share"] = "mca_minus_aca_share"

    X_full = df[feature_cols].fillna(0).values.astype(float)
    # Demographics-only null
    demo_cols = [c for c in ["age", "sex"] if c in df.columns]
    X_demo = pd.get_dummies(df[demo_cols], drop_first=True).fillna(0).values.astype(float) if demo_cols else np.ones((len(df), 1))

    pipe = Pipeline(
        [
            ("sc", StandardScaler()),
            ("ridge", RidgeCV(alphas=np.logspace(-2, 3, 20))),
        ]
    )

    rows = []
    for oname, ocol in outcomes.items():
        y = df[ocol].values.astype(float)
        # Nested-ish: outer 5-fold OOS predictions
        kf = KFold(n_splits=5, shuffle=True, random_state=0)
        oos_full = cross_val_predict(clone(pipe), X_full, y, cv=kf)
        oos_demo = cross_val_predict(clone(pipe), X_demo, y, cv=kf)

        def r2(y_true, y_pred):
            ss_res = ((y_true - y_pred) ** 2).sum()
            ss_tot = ((y_true - y_true.mean()) ** 2).sum()
            return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

        r2_full = r2(y, oos_full)
        r2_demo = r2(y, oos_demo)
        # Spearman
        from scipy.stats import spearmanr

        sp_full = spearmanr(y, oos_full)
        sp_demo = spearmanr(y, oos_demo)

        # Split conformal on last fold style: use 60% train, 20% cal, 20% test once for coverage check
        n = len(y)
        idx = np.arange(n)
        rng = np.random.default_rng(1)
        rng.shuffle(idx)
        n_tr = int(0.6 * n)
        n_cal = int(0.2 * n)
        tr, cal, te = idx[:n_tr], idx[n_tr : n_tr + n_cal], idx[n_tr + n_cal :]
        model = clone(pipe).fit(X_full[tr], y[tr])
        pred_cal = model.predict(X_full[cal])
        pred_te = model.predict(X_full[te])
        lo, hi, q = _conformal_intervals(y[cal], pred_cal, pred_te, alpha=0.1)
        cover = float(np.mean((y[te] >= lo) & (y[te] <= hi)))

        rows.append(
            {
                "outcome": oname,
                "oos_r2_cow_features": r2_full,
                "oos_r2_demo_only": r2_demo,
                "delta_oos_r2": r2_full - r2_demo,
                "oos_spearman_cow": float(sp_full.correlation) if sp_full.correlation is not None else np.nan,
                "oos_spearman_p": float(sp_full.pvalue),
                "oos_spearman_demo": float(sp_demo.correlation) if sp_demo.correlation is not None else np.nan,
                "conformal_q90": q,
                "conformal_coverage_90": cover,
                "n": n,
            }
        )
        pd.DataFrame(
            {
                "y": y,
                "oos_pred_cow": oos_full,
                "oos_pred_demo": oos_demo,
            }
        ).to_csv(OUT / f"oos_predictions_{oname}.csv", index=False)

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "cow_iv_nested_cv_conformal.csv", index=False)

    # Verdict: any outcome with delta_oos_r2 > 0.05 and positive oos_r2?
    success = bool(((tab["oos_r2_cow_features"] > 0) & (tab["delta_oos_r2"] > 0.05)).any())
    results["models"] = tab.to_dict(orient="records")
    results["features"] = feature_cols
    results["track_C_topology_support"] = "PASS_EXPLORATORY" if success else "FAIL_OOS"
    results["note"] = (
        "PASS only if OOS R²>0 and improves on demographics by >0.05 on any territorial share outcome."
    )

    (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(
        "\n".join(
            [
                "# Moonshot 04 — CoW topology IV + conformal",
                "",
                f"**Honesty lock:** {HONESTY['4']}",
                "",
                f"- n QC-pass: {results['n_qc']}",
                f"- Track C topology support: **{results['track_C_topology_support']}**",
                "",
                tab.to_markdown(index=False),
            ]
        ),
        encoding="utf-8",
    )
    print("DONE 04", results["track_C_topology_support"])
    return results


if __name__ == "__main__":
    run()
