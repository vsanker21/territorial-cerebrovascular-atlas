"""
Build ADNI Aim-4 analysis inventory from on-disk tables (no full unzip).

Joins:
  - All_Subjects_Key_MRI (FLAIR / T1w series)
  - UCD_WMH + ADSP_PHC_FLAIR volumes
  - PTDEMOG + RMT_APOERES
  - ADNIMERGE2 ADSL/DXSUM if extractable via pyreadr / rpy2; else demog only

Outputs under results/aim4/adni/
"""
from __future__ import annotations

import json
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ADNI = ROOT / "data/raw/ADNI"
OUT = ROOT / "results/aim4/adni"
OUT.mkdir(parents=True, exist_ok=True)


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def rid_from_ptid(s: pd.Series) -> pd.Series:
    return s.astype(str).str.extract(r"_S_0*(\d+)$", expand=False).astype(float)


def load_adsl_dx():
    """Try to pull ADSL + DXSUM from ADNIMERGE2.tar.gz via pyreadr."""
    tar_path = ADNI / "ADNIMERGE2.tar.gz"
    out = {}
    try:
        import pyreadr
    except ImportError:
        log("pyreadr not installed — clinical limited to PTDEMOG/APOE CSVs")
        return out
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        with tarfile.open(tar_path, "r:gz") as tar:
            for member in ("ADNIMERGE2/data/ADSL.rda", "ADNIMERGE2/data/DXSUM.rda"):
                try:
                    tar.extract(member, path=td)
                except KeyError:
                    continue
        for key, rel in [("ADSL", "ADNIMERGE2/data/ADSL.rda"), ("DXSUM", "ADNIMERGE2/data/DXSUM.rda")]:
            p = td / rel
            if not p.exists():
                continue
            try:
                obj = pyreadr.read_r(str(p))
                name = list(obj.keys())[0]
                out[key] = obj[name]
                log(f"Loaded {key} n={len(out[key])}")
            except Exception as e:
                log(f"Skip {key} ({type(e).__name__}: {e})")
    return out


def main():
    key = pd.read_csv(ADNI / "All_Subjects_Key_MRI_24Jul2026.csv")
    key["RID"] = rid_from_ptid(key["subject_id"])
    flair = key[key["series_type"].astype(str).str.upper() == "FLAIR"].copy()
    t1 = key[key["series_type"].astype(str).str.upper() == "T1W"].copy()
    flair = flair[flair["magnetic_field_strength"].astype(float) >= 2.5]
    t1 = t1[pd.to_numeric(t1["magnetic_field_strength"], errors="coerce").fillna(3) >= 2.5]

    wmh = pd.read_csv(ADNI / "UCD_WMH_24Jul2026.csv")
    phc = pd.read_csv(ADNI / "ADSP_PHC_FLAIR_24Jul2026.csv")
    dem = pd.read_csv(ADNI / "PTDEMOG_24Jul2026.csv")
    apoe = pd.read_csv(ADNI / "RMT_APOERES_24Jul2026.csv")

    # earliest FLAIR visit per subject (prefer sc/bl)
    vis_rank = {v: i for i, v in enumerate(
        ["sc", "scmri", "bl", "init", "4_sc", "v02", "m03", "m06", "m12", "m24"]
    )}
    flair["vis_rank"] = flair["image_visit"].astype(str).str.lower().map(vis_rank).fillna(99)
    flair_sorted = flair.sort_values(["subject_id", "vis_rank", "image_date"])
    flair_bl = flair_sorted.groupby("subject_id", as_index=False).first()

    # match UCD WMH nearest visit by RID (use first UCD row by examdate)
    wmh = wmh.sort_values(["RID", "EXAMDATE"])
    wmh_bl = wmh.groupby("RID", as_index=False).first()

    cohort = flair_bl.merge(wmh_bl, on="RID", how="left", suffixes=("_flair", "_wmh"))
    # demog: one row per RID (first)
    dem_bl = dem.sort_values(["RID", "VISCODE"]).groupby("RID", as_index=False).first()
    keep_dem = [c for c in ["RID", "PTGENDER", "PTDOBYY", "PTEDUCAT", "PTID"] if c in dem_bl.columns]
    cohort = cohort.merge(dem_bl[keep_dem], on="RID", how="left")

    # APOE via ADNIOnlineID if possible — RMT table may not have RID; leave as side table
    r = {
        "n_flair_series_3t": int(len(flair)),
        "n_flair_subjects_3t": int(flair["subject_id"].nunique()),
        "n_t1w_series_3t": int(len(t1)),
        "n_cohort_primary_flair": int(len(cohort)),
        "n_with_ucd_wmh": int(cohort["TOTAL_WMH"].notna().sum()) if "TOTAL_WMH" in cohort.columns else 0,
        "n_with_sex": int(cohort["PTGENDER"].notna().sum()) if "PTGENDER" in cohort.columns else 0,
        "flair_acq_type": flair["acquisition_type"].astype(str).value_counts().to_dict(),
        "spatial_wmh_masks_on_disk": False,
        "path": "local_FLAIR_segmentation_for_territory_fractions",
    }

    # PHC harmonized WMH join on RID+nearest
    phc_bl = phc.sort_values(["RID", "PHC_SCANDATE"]).groupby("RID", as_index=False).first()
    cohort = cohort.merge(
        phc_bl[["RID", "PHC_WMHVOL_Harmonized", "PHC_LONIUID", "PHC_SCANDATE"]],
        on="RID",
        how="left",
    )
    r["n_with_phc_wmh"] = int(cohort["PHC_WMHVOL_Harmonized"].notna().sum())

    clinical = load_adsl_dx()
    r["adsl_joined"] = False
    r["dxsum_joined"] = False
    if "DXSUM" in clinical:
        dx = clinical["DXSUM"].copy()
        if "RID" in dx.columns:
            dx["RID"] = pd.to_numeric(dx["RID"], errors="coerce")
            # baseline-ish diagnosis
            if "VISCODE" in dx.columns:
                dx["_rank"] = dx["VISCODE"].astype(str).str.lower().map(
                    {"sc": 0, "bl": 1, "scmri": 0}
                ).fillna(50)
                dx = dx.sort_values(["RID", "_rank"])
            dx_bl = dx.groupby("RID", as_index=False).first()
            dx_cols = [c for c in dx_bl.columns if c.upper() in {"RID", "DIAGNOSIS", "DXCURREN", "DXCHANGE", "DX"}
                       or "DIAG" in c.upper()]
            if "RID" not in dx_cols:
                dx_cols = ["RID"] + [c for c in dx_bl.columns if c != "RID"][:6]
            cohort = cohort.merge(dx_bl[dx_cols].drop_duplicates("RID"), on="RID", how="left")
            r["dxsum_joined"] = True
    if "ADSL" in clinical:
        r["adsl_joined"] = True

    out_csv = OUT / "adni_aim4_cohort_inventory.csv"
    cohort.to_csv(out_csv, index=False)
    (OUT / "adni_aim4_inventory_summary.json").write_text(json.dumps(r, indent=2), encoding="utf-8")

    # volume-only baseline table for immediate Aim 4
    vol_cols = [
        c
        for c in [
            "subject_id",
            "RID",
            "image_visit",
            "image_date",
            "image_id",
            "acquisition_type",
            "TOTAL_WMH",
            "TOTAL_BRAIN",
            "CEREBRUM_WHITE",
            "PHC_WMHVOL_Harmonized",
            "PTGENDER",
            "PTDOBYY",
            "PTEDUCAT",
        ]
        if c in cohort.columns
    ]
    # add any dx cols
    for c in cohort.columns:
        if c.upper().startswith("DX") or c.upper() in {"DIAGNOSIS", "DX_BL"}:
            if c not in vol_cols:
                vol_cols.append(c)
    vol = cohort[vol_cols].copy()
    if "TOTAL_WMH" in vol.columns and "TOTAL_BRAIN" in vol.columns:
        vol["wmh_frac_brain"] = vol["TOTAL_WMH"] / vol["TOTAL_BRAIN"].replace(0, np.nan)
    vol.to_csv(OUT / "adni_aim4_wmh_volume_baseline.csv", index=False)

    md = [
        "# ADNI Aim 4 inventory — existing data (no IDA WMH masks)",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d')}",
        "",
        "## Locked path",
        "- Spatial WMH masks **not** available on LONI for this account.",
        "- Proceed with: UCD/ADSP **volume** outcomes + **local FLAIR→WMH segmentation** for TerritoryVoxels fractions.",
        "",
        "## Counts",
        f"- 3T FLAIR series: {r['n_flair_series_3t']} ({r['n_flair_subjects_3t']} subjects)",
        f"- Primary cohort (1 FLAIR/subject): {r['n_cohort_primary_flair']}",
        f"- With UCD TOTAL_WMH: {r['n_with_ucd_wmh']}",
        f"- With ADSP PHC harmonized WMH: {r['n_with_phc_wmh']}",
        f"- ADSL join: {r['adsl_joined']}",
        "",
        "## On-disk imaging packs",
        "- T1: `Clinical_T1w_Imaging_Cohort_MRI*.zip`",
        "- FLAIR: `2D T2 FLAIR.zip`, `2D T2 FLAIR_dataset.zip`",
        "",
        "## Outputs",
        f"- `{out_csv.relative_to(ROOT)}`",
        f"- `results/aim4/adni/adni_aim4_wmh_volume_baseline.csv`",
        "",
    ]
    (OUT / "adni_aim4_inventory_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()
