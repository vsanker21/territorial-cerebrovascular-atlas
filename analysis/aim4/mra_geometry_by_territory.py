"""
MRA vessel-geometry covariates summarized within Liu TerritoryVoxels.

Provides the Aim-4 nested-null 'geometry' layer without requiring ATLAS lesions yet.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import image

ROOT = Path(__file__).resolve().parents[2]
NITRC = ROOT / "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/Atlas_182_MNI152"
TV_PATH = NITRC / "TerritoryVoxels_BMM.nii"
MRA_DIR = ROOT / "data/raw/mra_atlas/figshare"
OUT_CSV = ROOT / "results/aim4/mra_geometry_by_territory.csv"
OUT_MD = ROOT / "results/aim4/mra_geometry_by_territory.md"

TV_LABELS = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
MRA_FILES = {
    "vessel_probability": "vesselProbabilities.nii.gz",
    "vessel_radius": "vesselRadius.nii.gz",
    "vessel_radius_std": "vesselRadiusStd.nii.gz",
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    if not TV_PATH.exists():
        raise SystemExit(f"Missing {TV_PATH}")
    tv_img = nib.load(str(TV_PATH))
    tv = np.rint(np.asanyarray(tv_img.dataobj, dtype=float)).astype(np.int32)
    if tv.ndim == 4:
        tv = tv[..., 0]

    rows = []
    for key, fname in MRA_FILES.items():
        path = MRA_DIR / fname
        if not path.exists():
            log(f"Missing {path}")
            continue
        log(f"Loading {fname}")
        img = nib.load(str(path))
        if img.shape[:3] != tv.shape[:3]:
            log(f"Resampling {fname} to TerritoryVoxels grid")
            img = image.resample_to_img(img, tv_img, interpolation="linear")
        data = np.asanyarray(img.dataobj, dtype=float)
        if data.ndim == 4:
            data = data[..., 0]
        for lab, name in TV_LABELS.items():
            mask = tv == lab
            vals = data[mask]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            rows.append(
                {
                    "metric": key,
                    "territory": name,
                    "tv_label": lab,
                    "n_voxels": int(vals.size),
                    "mean": float(vals.mean()),
                    "std": float(vals.std()),
                    "median": float(np.median(vals)),
                    "p90": float(np.percentile(vals, 90)),
                }
            )

    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    # pivot for readability
    if len(df):
        piv = df.pivot_table(index="territory", columns="metric", values="mean")
        md = [
            "# MRA geometry by Liu TerritoryVoxels",
            "",
            "Mean vessel metrics inside TerritoryVoxels labels 1–4 (Aim-4 geometry null).",
            "",
            piv.to_string(),
            "",
            f"Wrote `{OUT_CSV.relative_to(ROOT).as_posix()}`",
        ]
        OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
        print(piv.to_string())
    log("Done")


if __name__ == "__main__":
    main()
