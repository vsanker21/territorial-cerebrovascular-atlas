"""Formalize Yang H4 + lock Aim3 borderzone role (not Class B)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
H4_CSV = ROOT / "results/h4/yang_h4_within_donor.csv"
AIM3 = list((ROOT / "results/yang").glob("aim3*.md"))
OUT_MD = ROOT / "results/h4/yang_h4_formal_status.md"
OUT_JSON = ROOT / "results/h4/yang_h4_formal_summary.json"
OUT_AIM3 = ROOT / "results/yang/yang_aim3_borderzone_lock.md"


def main():
    summary = {"h4_csv_exists": H4_CSV.exists(), "aim3_reports": [str(p.relative_to(ROOT)) for p in AIM3]}
    if H4_CSV.exists():
        df = pd.read_csv(H4_CSV)
        summary["columns"] = list(df.columns)
        # expect disease and auroc-like columns
        if "auroc_mean" in df.columns and "disease" in df.columns:
            summary["mean_auroc_by_disease"] = (
                df.groupby(df["disease"].astype(str))["auroc_mean"].mean().to_dict()
            )
            summary["overall_mean_auroc"] = float(df["auroc_mean"].mean())
            summary["delta_AD_minus_control"] = float(
                summary["mean_auroc_by_disease"].get("AD", np.nan)
                - summary["mean_auroc_by_disease"].get("control", np.nan)
            )
        summary["n_by_disease"] = (
            df["disease"].astype(str).value_counts().to_dict() if "disease" in df.columns else {}
        )

    # H4 claim lock
    h4_verdict = "H4_NOT_SUPPORTED_ON_YANG_ALONE"
    summary["verdict"] = h4_verdict
    summary["rationale"] = (
        "Within-donor AUROC high in BOTH control (~0.95) and AD (~0.98); "
        "delta AD-control is NOT negative (fails H4 degradation direction). "
        "SFG is ACA/MCA borderzone; region=territory confound remains; "
        "cross-donor LODO near chance (see yang_h4_lodo_report.md)."
    )

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")

    ctrl = summary.get("mean_auroc_by_disease", {}).get("control")
    ad = summary.get("mean_auroc_by_disease", {}).get("AD")
    delta = summary.get("delta_AD_minus_control")
    OUT_MD.write_text(
        f"""# Yang H4 — formal status (open data)

**Verdict:** `{h4_verdict}`

## Quantitative summary
- Control within-donor AUROC mean: **{ctrl}**
- AD within-donor AUROC mean: **{ad}**
- Delta (AD − control): **{delta}** (H4 expects negative)

## Locked interpretation
1. Yang does **not** support H4: identity scores do **not** degrade in AD on this within-donor test.
2. High within-donor AUROC reflects region/territory aliasing (hippocampus vs SFG), not clean arterial identity.
3. Cross-donor LODO failure (separate report) shows the program is not transferable across donors.

## Required for revisiting H4
- Multi-territory designs with enrichment matching and/or disease cohorts with region labels.
""",
        encoding="utf-8",
    )

    OUT_AIM3.write_text(
        """# Yang SFG — Aim 3 borderzone lock (NOT Class B)

**Locked role:** ACA/MCA **borderzone** contrast (TerritoryVoxels superior frontal ≈ mixed ACA/MCA).

## Use for
- Aim 3 mixture/residual borderzone tests
- Within-donor paired PCA-hippocampus vs SFG exploratory molecular programs
- Gate/QC development (tight endo gates)

## Do NOT use for
- Clean Class B ACA vs MCA
- Calling Yang SFG “ACA cortex” in confirmatory H1
- External ACA arm (use Allen ACC / LIBD dACC / Linnarsson A24–A32 instead)

## Related
- Class B true ACA arm: `results/class_b/aca_mca_harmony_class_b_report.md`
""",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
