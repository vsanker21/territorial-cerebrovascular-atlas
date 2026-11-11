"""ADNI Aim4 WMH geography stratified by DIAGNOSIS (CN/MCI/Dementia)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
NEST = ROOT / "results/aim4/adni/aim4_nest_wmh_mra_subject_table.csv"
BASE = ROOT / "results/aim4/adni/adni_aim4_wmh_volume_baseline.csv"
OUT = ROOT / "results/aim4/adni"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nest = pd.read_csv(NEST)
    base = pd.read_csv(BASE)
    # join keys
    if "RID" in nest.columns and "RID" in base.columns:
        d = nest.merge(base[["RID", "DIAGNOSIS"]].drop_duplicates("RID"), on="RID", how="left")
    else:
        raise SystemExit("RID missing")
    d["dx"] = d["DIAGNOSIS"].astype(str).str.strip()
    # normalize
    map_dx = {
        "CN": "CN",
        "MCI": "MCI",
        "Dementia": "Dementia",
        "AD": "Dementia",
        "nan": "unknown",
        "None": "unknown",
    }
    d["dx_group"] = d["dx"].map(lambda x: map_dx.get(x, x if x in ("CN", "MCI", "Dementia") else "other"))

    def summarize(df, label):
        if len(df) == 0:
            return {"stratum": label, "n": 0}
        mca = df["wmh_frac_of_wmh_MCA"]
        aca = df["wmh_frac_of_wmh_ACA"]
        return {
            "stratum": label,
            "n": int(len(df)),
            "mean_share_MCA": float(mca.mean()),
            "mean_share_ACA": float(aca.mean()),
            "frac_subjects_MCA_gt_ACA": float((mca > aca).mean()),
            "median_TOTAL_WMH_ucd": float(df["TOTAL_WMH_ucd"].median()),
            "median_n_wmh": float(df["n_wmh_voxels"].median()),
        }

    rows = [summarize(d, "all_primary")]
    for g in ["CN", "MCI", "Dementia", "other", "unknown"]:
        rows.append(summarize(d[d["dx_group"] == g], g))

    tests = {}
    cn = d.loc[d["dx_group"] == "CN", "wmh_frac_of_wmh_MCA"]
    mci = d.loc[d["dx_group"] == "MCI", "wmh_frac_of_wmh_MCA"]
    if len(cn) and len(mci):
        u = stats.mannwhitneyu(cn, mci, alternative="two-sided")
        tests["MCA_share_CN_vs_MCI"] = {"U": float(u.statistic), "p": float(u.pvalue)}

    out = {
        "strata": rows,
        "tests": tests,
        "dx_match_rate": float(d["DIAGNOSIS"].notna().mean()),
        "interpretation_locks": [
            "Descriptive sensitivity only.",
            "Primary geography lock remains mean MCA share > mean ACA share on primary QC set.",
            "Subject-level MCA>ACA fraction may be near 50% even when mean share is MCA-dominant.",
            "Dementia n is exploratory if small.",
        ],
    }
    pd.DataFrame(rows).to_csv(OUT / "aim4_wmh_dx_strata.csv", index=False)
    (OUT / "aim4_wmh_dx_summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    md = [
        "# ADNI Aim4 WMH geography — DIAGNOSIS strata",
        "",
        f"Primary QC n={len(d)}; DIAGNOSIS match rate={out['dx_match_rate']:.2f}",
        "",
        pd.DataFrame(rows).to_string(index=False),
        "",
        f"Tests: {tests}",
        "",
        "Locked: mean MCA share > mean ACA share remains the primary geography claim;",
        "report subject-level MCA>ACA fractions as secondary descriptive stats.",
        "",
    ]
    (OUT / "aim4_wmh_dx_report.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
