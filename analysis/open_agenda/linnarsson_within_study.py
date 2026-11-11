"""Linnarsson within-study vascular note + limited within-ACA dissection contrast.

Lateral frontal / temporal Siletti slices are NOT on disk → powered F/T comparison NO-GO.
Document A24/A32/A25/A23 endo counts and A14 MCA-leaning exclusion.
"""
from __future__ import annotations

import json
from pathlib import Path

import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
LIN = ROOT / "data/processed/aca_cohort_hunt/linnarsson_aca_vascular.h5ad"
PRIORS = ROOT / "results/aca_cohort_hunt/linnarsson_aca_tv_priors.json"
OUT_MD = ROOT / "results/h1/linnarsson_within_study_report.md"
OUT_JSON = ROOT / "results/h1/linnarsson_within_study_summary.json"


def main():
    priors = json.loads(PRIORS.read_text(encoding="utf-8"))
    a = sc.read_h5ad(LIN)
    endo = a[a.obs["cell_type"].astype(str) == "endothelial cell"].copy()
    # dissection column?
    dcol = None
    for c in ["dissection", "tissue", "region", "subregion"]:
        if c in endo.obs:
            dcol = c
            break
    by = endo.obs[dcol].astype(str).value_counts().to_dict() if dcol else {}
    summary = {
        "n_vascular": int(a.n_obs),
        "n_endo": int(endo.n_obs),
        "endo_by_dissection": by,
        "donors": int(endo.obs["donor_id"].nunique()) if "donor_id" in endo.obs else None,
        "priors": priors,
        "lateral_frontal_temporal_on_disk": False,
        "verdict": "NO-GO_powered_frontal_vs_temporal",
        "limited_go": "ACA cingulate pool seed only; within-dissection endo n too sparse for donor OOF",
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = f"""# Linnarsson / Siletti within-study vascular

**Verdict:** `{summary['verdict']}`

## On disk
- Processed ACA vascular pool: `{LIN.relative_to(ROOT)}` — n vascular={a.n_obs}, endo={endo.n_obs}
- Dissections present: A24, A32, A25, A23 (A14 excluded as MCA orbital)
- Endo by dissection: `{by}`
- Donors: {summary['donors']}

## NO-GO
Lateral frontal and temporal cortex Siletti slices were **not downloaded**.  
Cannot run A24/A32 vs lateral frontal vs temporal within-study Class B analog without acquiring those CELLxGENE dissections.

## Limited use
TerritoryVoxels-mapped ACA seed for external Class B arm only (already used). Within-cingulate endo n is too small for powered donor-OOF H1.
"""
    OUT_MD.write_text(md, encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
