"""ADNI Aim4 geography stratified by sex and APOE (primary QC set).

APOE genotypes come from ADNIMERGE2/data/APOERES.rda (has RID+GENOTYPE).
The flat RMT_APOERES_*.csv only has ADNIOnlineID and cannot join to RID.
"""
from __future__ import annotations

import json
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
NEST = ROOT / "results/aim4/adni/aim4_nest_wmh_mra_subject_table.csv"
ADNI = ROOT / "data/raw/ADNI"
OUT = ROOT / "results/aim4/adni"
APOE_CACHE = OUT / "adni_apoeres_rid_genotype.csv"


def load_apoe_by_rid() -> pd.DataFrame:
    if APOE_CACHE.exists():
        return pd.read_csv(APOE_CACHE)
    tar_path = ADNI / "ADNIMERGE2.tar.gz"
    try:
        import pyreadr
    except ImportError as e:
        raise RuntimeError("pyreadr required to extract APOERES.rda") from e
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extract("ADNIMERGE2/data/APOERES.rda", path=td)
        obj = pyreadr.read_r(str(td / "ADNIMERGE2/data/APOERES.rda"))
        ap = list(obj.values())[0]
    keep = ap[["RID", "GENOTYPE"]].dropna().copy()
    keep["RID"] = pd.to_numeric(keep["RID"], errors="coerce")
    keep = keep.dropna(subset=["RID"]).drop_duplicates("RID")
    keep.to_csv(APOE_CACHE, index=False)
    return keep


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nest = pd.read_csv(NEST)
    nest["RID"] = pd.to_numeric(nest["RID"], errors="coerce")
    ap = load_apoe_by_rid()
    ap = ap.rename(columns={"GENOTYPE": "APOE"})
    d = nest.merge(ap, on="RID", how="left")
    d["APOE_e4"] = d["APOE"].astype(str).str.contains("4").astype(float)
    d.loc[d["APOE"].isna(), "APOE_e4"] = np.nan

    def summarize(df, label):
        return {
            "stratum": label,
            "n": int(len(df)),
            "mean_share_MCA": float(df["wmh_frac_of_wmh_MCA"].mean()) if len(df) else None,
            "mean_share_ACA": float(df["wmh_frac_of_wmh_ACA"].mean()) if len(df) else None,
            "median_n_wmh": float(df["n_wmh_voxels"].median()) if len(df) else None,
            "median_TOTAL_WMH_ucd": float(df["TOTAL_WMH_ucd"].median()) if len(df) else None,
        }

    rows = [
        summarize(d, "all_primary"),
        summarize(d[d["sex_f"] == 1], "female"),
        summarize(d[d["sex_f"] == 0], "male"),
        summarize(d[d["APOE_e4"] == 1], "APOE_e4_carrier"),
        summarize(d[d["APOE_e4"] == 0], "APOE_e4_noncarrier"),
    ]
    tests = {}
    f = d.loc[d["sex_f"] == 1, "wmh_frac_of_wmh_MCA"]
    m = d.loc[d["sex_f"] == 0, "wmh_frac_of_wmh_MCA"]
    if len(f) and len(m):
        u = stats.mannwhitneyu(f, m, alternative="two-sided")
        tests["MCA_share_female_vs_male"] = {"U": float(u.statistic), "p": float(u.pvalue)}
    e4 = d.loc[d["APOE_e4"] == 1, "wmh_frac_of_wmh_MCA"]
    non = d.loc[d["APOE_e4"] == 0, "wmh_frac_of_wmh_MCA"]
    if len(e4) and len(non):
        u = stats.mannwhitneyu(e4, non, alternative="two-sided")
        tests["MCA_share_e4_vs_non"] = {"U": float(u.statistic), "p": float(u.pvalue)}

    out = {
        "strata": rows,
        "tests": tests,
        "apoe_match_rate": float(d["APOE"].notna().mean()),
        "apoe_source": "ADNIMERGE2/data/APOERES.rda (RID+GENOTYPE)",
        "note": "Descriptive sensitivity only; does not change primary MCA-dominant geography lock.",
    }
    pd.DataFrame(rows).to_csv(OUT / "aim4_wmh_apoe_sex_strata.csv", index=False)
    (OUT / "aim4_wmh_apoe_sex_summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    md = [
        "# ADNI Aim4 WMH geography — sex / APOE strata",
        "",
        f"Primary QC n={len(d)}; APOE match rate={out['apoe_match_rate']:.2f}",
        f"APOE source: `{out['apoe_source']}`",
        "",
        pd.DataFrame(rows).to_string(index=False),
        "",
        f"Tests: {tests}",
        "",
        "Locked: descriptive sensitivity only; does not change primary MCA-dominant geography lock.",
        "",
    ]
    (OUT / "aim4_wmh_apoe_sex_report.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
