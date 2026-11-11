"""Moonshot #6 — Hierarchical Bayesian ecological territorial vulnerability.

ECOLOGICAL — not confirmatory molecular H3.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from _paths import ADNI_WMH, ATLAS_FRAC, HONESTY, ISLES_FRAC, MODULE_BY_DONOR, RES

OUT = RES / "06_bayesian_ecological"
OUT.mkdir(parents=True, exist_ok=True)

TERR = ["ACA", "MCA", "PCA", "VB"]


def _mean_shares(df: pd.DataFrame, prefix: str) -> dict:
    """prefix e.g. lesion_frac_of_lesion_ or wmh_frac_of_wmh_"""
    out = {}
    for t in TERR:
        col = f"{prefix}{t}"
        if col in df.columns:
            out[t] = float(df[col].mean())
        else:
            # atlas file uses lesion_frac_of_lesion_
            alt = f"lesion_frac_of_lesion_{t}"
            out[t] = float(df[alt].mean()) if alt in df.columns else np.nan
    return out


def _gibbs_hierarchical(y, x, terr_idx, n_terr, n_iter=4000, burn=1000, seed=0):
    """
    y_i = alpha[terr] + beta * x_i + eps
    alpha[t] ~ N(mu, tau^2); semi-conjugate Gibbs (known noise via IG).
    Ecological: few points — still returns posterior samples.
    """
    rng = np.random.default_rng(seed)
    n = len(y)
    # init
    alpha = np.zeros(n_terr)
    beta = 0.0
    mu = 0.0
    tau2 = 1.0
    sig2 = 1.0
    samples = []
    for it in range(n_iter):
        # alpha_t | rest
        for t in range(n_terr):
            idx = np.where(terr_idx == t)[0]
            if len(idx) == 0:
                alpha[t] = rng.normal(mu, np.sqrt(tau2))
                continue
            prec = len(idx) / sig2 + 1 / tau2
            mean = (np.sum(y[idx] - beta * x[idx]) / sig2 + mu / tau2) / prec
            alpha[t] = rng.normal(mean, 1 / np.sqrt(prec))
        # beta
        r = y - alpha[terr_idx]
        prec_b = np.sum(x ** 2) / sig2 + 1e-6
        mean_b = np.sum(x * r) / sig2 / prec_b
        beta = rng.normal(mean_b, 1 / np.sqrt(prec_b))
        # mu
        prec_m = n_terr / tau2 + 1e-6
        mu = rng.normal(alpha.mean() * (n_terr / tau2) / prec_m, 1 / np.sqrt(prec_m))  # approx
        mu = rng.normal(alpha.mean(), np.sqrt(tau2 / max(n_terr, 1)))
        # tau2 IG
        ss = np.sum((alpha - mu) ** 2) + 1e-3
        tau2 = 1 / rng.gamma(2 + n_terr / 2, 1 / (1 + 0.5 * ss))
        # sig2 IG
        resid = y - alpha[terr_idx] - beta * x
        sig2 = 1 / rng.gamma(2 + n / 2, 1 / (1 + 0.5 * np.sum(resid ** 2)))
        if it >= burn and it % 5 == 0:
            samples.append({"beta": beta, "mu": mu, "tau2": tau2, "sig2": sig2, **{f"a_{t}": alpha[t] for t in range(n_terr)}})
    return pd.DataFrame(samples)


def run() -> dict:
    results = {"honesty_lock": HONESTY["6"], "label": "ECOLOGICAL_NOT_CONFIRMATORY_MOLECULAR_H3"}

    # Imaging geography means
    atlas = pd.read_csv(ATLAS_FRAC)
    isles = pd.read_csv(ISLES_FRAC)
    adni = pd.read_csv(ADNI_WMH)
    # QC filters if present
    if "status" in atlas.columns:
        atlas = atlas[atlas["status"].astype(str).str.contains("ok|pass", case=False, na=True)]
    shares = {
        "ATLAS": _mean_shares(atlas, "lesion_frac_of_lesion_"),
        "ISLES22": _mean_shares(isles, "lesion_frac_of_lesion_"),
        "ADNI_WMH": {
            t: float(adni[f"wmh_frac_of_wmh_{t}"].mean()) if f"wmh_frac_of_wmh_{t}" in adni.columns else np.nan
            for t in TERR
        },
    }
    pd.DataFrame(shares).T.to_csv(OUT / "mean_territory_shares_by_cohort.csv")

    # Molecular modules by territory from Class B donor table
    mod = pd.read_csv(MODULE_BY_DONOR)
    # expect columns with territory and module scores
    terr_col = None
    for c in ["tv_dominant", "territory", "arm", "class_b_arm"]:
        if c in mod.columns:
            terr_col = c
            break
    score_cols = [
        c
        for c in mod.columns
        if any(k in c.lower() for k in ["bbb", "angio", "ecm", "hypox", "inflam", "module", "score"])
        or c.lower().startswith("mod")
    ]
    if not score_cols:
        # numeric columns excluding ids
        score_cols = [
            c
            for c in mod.columns
            if mod[c].dtype != object and c.lower() not in {"donor", "n_cells", "age"}
        ]

    # Build ecological table: for each territory, mean module + mean imaging share (ATLAS primary)
    eco_rows = []
    if terr_col:
        for t in TERR:
            sub = mod[mod[terr_col].astype(str).str.contains(t, case=False, na=False)]
            if len(sub) == 0:
                # try exact
                sub = mod[mod[terr_col].astype(str) == t]
            row = {"territory": t, "n_donors_mol": int(len(sub))}
            for c in score_cols:
                row[f"mol_{c}"] = float(sub[c].mean()) if len(sub) else np.nan
            for cohort, sh in shares.items():
                row[f"share_{cohort}"] = sh.get(t, np.nan)
            eco_rows.append(row)
    else:
        # territory-level from imaging only + placeholder molecular from within-garcia if needed
        for t in TERR:
            eco_rows.append(
                {
                    "territory": t,
                    "n_donors_mol": 0,
                    **{f"share_{k}": v.get(t, np.nan) for k, v in shares.items()},
                }
            )

    eco = pd.DataFrame(eco_rows)
    eco.to_csv(OUT / "ecological_territory_table.csv", index=False)

    # Hierarchical model: predict ATLAS share from first molecular PC / mean module
    mol_cols = [c for c in eco.columns if c.startswith("mol_")]
    y_cohorts = ["share_ATLAS", "share_ISLES22", "share_ADNI_WMH"]
    model_summaries = []

    for ycol in y_cohorts:
        if ycol not in eco.columns:
            continue
        y = eco[ycol].values.astype(float)
        if np.isnan(y).all():
            continue
        if mol_cols:
            Xmat = eco[mol_cols].apply(pd.to_numeric, errors="coerce").values.astype(float)
            # fill nan with col mean
            for j in range(Xmat.shape[1]):
                col = Xmat[:, j]
                if np.isnan(col).any():
                    Xmat[np.isnan(col), j] = np.nanmean(col)
            # use mean molecular score as x
            x = np.nanmean(Xmat, axis=1)
        else:
            # null: use ordinal territory code only (no molecular) — still ecological geography
            x = np.arange(len(eco), dtype=float)
        ok = np.isfinite(y) & np.isfinite(x)
        y, x = y[ok], x[ok]
        terr_idx = np.arange(len(y))  # each territory one obs — hierarchical prior still regularizes
        if len(y) < 3:
            continue
        post = _gibbs_hierarchical(y, x, terr_idx, n_terr=len(y), n_iter=6000, burn=2000)
        beta = post["beta"]
        model_summaries.append(
            {
                "outcome": ycol,
                "n_territories": int(len(y)),
                "beta_mean": float(beta.mean()),
                "beta_sd": float(beta.std()),
                "beta_ci95": [float(np.quantile(beta, 0.025)), float(np.quantile(beta, 0.975))],
                "P_beta_gt_0": float((beta > 0).mean()),
                "spearman_raw": float(stats.spearmanr(x, y).correlation) if len(y) >= 3 else np.nan,
            }
        )
        post.to_csv(OUT / f"posterior_{ycol}.csv", index=False)

        # Posterior predictive check: replicate y
        # simple: y_rep = a_t + beta * x
        y_rep_mean = []
        for t in range(len(y)):
            a_t = post[f"a_{t}"].values
            y_rep_mean.append(float(np.mean(a_t + post["beta"].values * x[t])))
        ppc = pd.DataFrame({"territory": eco.loc[ok, "territory"].values, "y_obs": y, "y_ppc_mean": y_rep_mean, "x": x})
        ppc.to_csv(OUT / f"ppc_{ycol}.csv", index=False)

    ms = pd.DataFrame(model_summaries)
    ms.to_csv(OUT / "hierarchical_model_summary.csv", index=False)

    results["shares"] = shares
    results["models"] = model_summaries
    results["molecular_cols_used"] = mol_cols
    results["verdict"] = (
        "Ecological hierarchical posteriors estimated; do NOT interpret as confirmatory molecular H3. "
        "Use PPC and beta CI for honest uncertainty at n_territories≈3–4."
    )

    (OUT / "summary.json").write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    (OUT / "report.md").write_text(
        "\n".join(
            [
                "# Moonshot 06 — Hierarchical Bayesian ecological vulnerability",
                "",
                f"**Honesty lock:** {HONESTY['6']}",
                "",
                "## Mean territory shares",
                pd.DataFrame(shares).T.to_markdown(),
                "",
                "## Hierarchical models",
                ms.to_markdown(index=False) if len(ms) else "_no model_",
                "",
                results["verdict"],
            ]
        ),
        encoding="utf-8",
    )
    print("DONE 06")
    return results


if __name__ == "__main__":
    run()
