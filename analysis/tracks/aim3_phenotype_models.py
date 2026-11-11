#!/usr/bin/env python3
"""Aim 3 phenotype models — runs only after RADC Standard extract lands.

Pre-specified models match Nature Tier-2 phenotype linkage (not Track B identity).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/tracks/aim3_radc"
OUT.mkdir(parents=True, exist_ok=True)

RADC_CANDIDATES = [
    ROOT / "data/raw/radc/Project7508_Standard_Clinical_Pathologic.csv",
    ROOT / "data/raw/radc/standard_clinical_pathologic_variables.csv",
    ROOT / "data/raw/radc/rosmap_standard.csv",
]

# Expected variable *classes* (exact RADC names harmonized on receipt)
EXPECTED_CLASSES = {
    "id": ["projid", "individualID", "proj_id"],
    "demographics": ["age_death", "age_bl", "msex", "sex", "educ"],
    "diagnosis": ["cogdx", "dcfdx", "diagnosis", "clin_dx"],
    "ad_path": ["braaksc", "ceradsc", "niareagansc", "np_reagan"],
    "cvd_path": [
        "arteriolosclerosis",
        "cvda",
        "caa_4gp",
        "amyloid",
        "ci_num2_gct",
        "ci_num2_mct",
        "atherosclerosis",
    ],
    "genetics": ["apoe_genotype", "apoe4"],
}


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # fuzzy contains
    for cand in candidates:
        for c in df.columns:
            if cand.lower() in c.lower():
                return c
    return None


def donor_molecular_scores() -> pd.DataFrame:
    """Load precomputed donor-level territorial molecular scores if present; else stub keys."""
    # Prefer existing Class B / module donor tables when available
    candidates = [
        ROOT / "results/tracks/nature_rigor/matched_mca_pca_donor_scores.csv",
        ROOT / "results/metadata_audit/garcia_rosmap_samples.csv",
    ]
    for p in candidates:
        if p.exists():
            df = pd.read_csv(p)
            if "projid" in df.columns or "subject" in df.columns:
                return df
    return pd.DataFrame()


def fit_assoc(y: pd.Series, X: pd.DataFrame) -> dict:
    """Simple OLS association with age/sex covariates when available — exploratory scaffold."""
    from numpy.linalg import lstsq

    mask = y.notna()
    for c in X.columns:
        mask &= X[c].notna()
    yy = y[mask].astype(float).values
    XX = X.loc[mask].astype(float).values
    if len(yy) < 20 or XX.shape[1] == 0:
        return {"n": int(mask.sum()), "status": "underpowered_or_empty"}
    XX = np.column_stack([np.ones(len(yy)), XX])
    coef, *_ = lstsq(XX, yy, rcond=None)
    resid = yy - XX @ coef
    sse = float(np.sum(resid**2))
    sst = float(np.sum((yy - yy.mean()) ** 2)) + 1e-12
    return {
        "n": int(len(yy)),
        "r2": float(1 - sse / sst),
        "coef": coef.tolist(),
        "status": "fit",
    }


def main() -> None:
    radc_path = next((p for p in RADC_CANDIDATES if p.exists()), None)
    status: dict = {
        "date": time.strftime("%Y-%m-%d"),
        "rra": "RRA-539487",
        "project": "RADC #7508",
        "expected_classes": EXPECTED_CLASSES,
        "not_identity_proof": True,
        "trci_bci_excluded": True,
    }

    if radc_path is None:
        status["verdict"] = "GATED_AWAITING_RADC_DUA_DELIVERY"
        (OUT / "aim3_phenotype_models_status.json").write_text(
            json.dumps(status, indent=2), encoding="utf-8"
        )
        (OUT / "aim3_phenotype_models_report.md").write_text(
            "\n".join(
                [
                    "# Aim 3 phenotype models",
                    "",
                    f"**Verdict:** `{status['verdict']}`",
                    "",
                    "Place Standard Clinical Pathologic Variables extract under `data/raw/radc/`, then re-run:",
                    "`python analysis/tracks/aim3_phenotype_models.py`",
                    "",
                    "Pre-specified outcomes: diagnosis, AD path scores, CVD neuropathology.",
                    "Covariates: age, sex, educ, APOE when available.",
                    "This is phenotype linkage — not developmental territorial identity.",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        log(json.dumps(status, indent=2))
        return

    radc = pd.read_csv(radc_path, low_memory=False)
    idcol = pick_col(radc, EXPECTED_CLASSES["id"])
    mol = donor_molecular_scores()
    status["radc_path"] = str(radc_path)
    status["radc_n"] = int(len(radc))
    status["id_col"] = idcol
    status["mapped_columns"] = {
        k: pick_col(radc, v) for k, v in EXPECTED_CLASSES.items() if k != "id"
    }

    if idcol is None:
        status["verdict"] = "RADC_RECEIVED__NO_ID_COLUMN"
        (OUT / "aim3_phenotype_models_status.json").write_text(
            json.dumps(status, indent=2), encoding="utf-8"
        )
        log(json.dumps(status, indent=2))
        return

    # Link
    if "projid" in mol.columns:
        linked = mol.merge(radc, left_on="projid", right_on=idcol, how="inner")
    elif "subject" in mol.columns:
        # try parse ROSMAP numeric from subject
        linked = mol.copy()
        linked["_join"] = linked["subject"].astype(str).str.extract(r"(\d+)", expand=False)
        radc["_join"] = radc[idcol].astype(str)
        linked = linked.merge(radc, on="_join", how="inner")
    else:
        linked = pd.DataFrame()

    linked.to_csv(OUT / "aim3_phenotype_linked.csv", index=False)
    status["n_linked"] = int(len(linked))

    # Prefer a molecular score column if present
    score_col = next(
        (
            c
            for c in linked.columns
            if c.lower() in {"donor_oof", "territorial_score", "endo_score", "module_hypoxia"}
            or "score" in c.lower()
        ),
        None,
    )
    results = {}
    if score_col and len(linked) >= 20:
        y = pd.to_numeric(linked[score_col], errors="coerce")
        cov_cols = [
            c
            for c in [
                status["mapped_columns"].get("demographics"),
            ]
            if c
        ]
        # build covariate matrix from mapped demography columns present
        X_parts = []
        names = []
        for cls in ("demographics", "genetics"):
            col = status["mapped_columns"].get(cls)
            # demographics mapping returns one col; expand known names present
        for c in linked.columns:
            cl = c.lower()
            if cl in {"age_death", "age_bl", "msex", "sex", "educ", "apoe4"} or "apoe" in cl:
                X_parts.append(pd.to_numeric(linked[c], errors="coerce"))
                names.append(c)
        X = pd.concat(X_parts, axis=1) if X_parts else pd.DataFrame(index=linked.index)
        if not X.empty:
            X.columns = names[: X.shape[1]]
        for outcome_cls in ("diagnosis", "ad_path", "cvd_path"):
            ocol = status["mapped_columns"].get(outcome_cls)
            if ocol and ocol in linked.columns:
                yy = pd.to_numeric(linked[ocol], errors="coerce")
                # association: molecular score ~ outcome + covariates (swap if outcome categorical ordinal)
                XX = X.copy()
                XX["outcome"] = yy
                results[outcome_cls] = fit_assoc(y, XX)

    status["models"] = results
    status["score_col"] = score_col
    status["verdict"] = (
        "LINKED__MODELS_FIT" if results else "LINKED__AWAITING_SCORE_OR_COLUMNS"
    )
    (OUT / "aim3_phenotype_models_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    (OUT / "aim3_phenotype_models_report.md").write_text(
        "\n".join(
            [
                "# Aim 3 phenotype models",
                "",
                f"**Verdict:** `{status['verdict']}`",
                f"- Linked n={status.get('n_linked')}",
                f"- Score column: `{score_col}`",
                f"- Mapped columns: `{json.dumps(status['mapped_columns'])}`",
                "",
                "Phenotype linkage only — not Track B identity; TRCI/BCI excluded.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
