"""ATLAS laterality sensitivity: L/R parcels before major-territory collapse."""
from __future__ import annotations

import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import image

ROOT = Path(__file__).resolve().parents[2]
ATLAS_ROOT = ROOT / "data/raw/atlas_stroke/extracted/ATLAS3_Training_Preprocessed"
LEVEL2 = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/ArterialAtlas_level2.nii"
)
QC = ROOT / "results/aim4/atlas/atlas_lesion_tv_fractions.csv"
OUT = ROOT / "results/aim4/atlas"

# level2: 1 ACAL 2 ACAR 3 MCAL 4 MCAR 5 PCAL 6 PCAR 7 VBL 8 VBR
PARCELS = {
    1: ("ACA", "L"),
    2: ("ACA", "R"),
    3: ("MCA", "L"),
    4: ("MCA", "R"),
    5: ("PCA", "L"),
    6: ("PCA", "R"),
    7: ("VB", "L"),
    8: ("VB", "R"),
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(QC)
    primary = set(qc.loc[qc["qc_primary"] == True, "subject_id"].astype(str))  # noqa: E712
    masks = sorted(ATLAS_ROOT.rglob("*label-lesion*mask*.nii.gz"))
    ref = nib.load(str(masks[0]))
    l2 = image.resample_to_img(LEVEL2, ref, interpolation="nearest", force_resample=True, copy_header=True)
    l2a = np.rint(np.asanyarray(l2.dataobj)).astype(np.int32)
    if l2a.ndim == 4:
        l2a = l2a[..., 0]

    rows = []
    for i, mpath in enumerate(masks, 1):
        sid = None
        for p in mpath.parts:
            if p.startswith("sub-"):
                sid = p
                break
        if sid not in primary:
            continue
        if mpath.stat().st_size < 100:
            continue
        if i % 200 == 0:
            log(f"{i} masks scanned; rows={len(rows)}")
        les = np.asanyarray(nib.load(str(mpath)).dataobj) > 0
        n = int(les.sum())
        if n <= 0:
            continue
        row = {"subject_id": sid, "n_lesion_voxels": n}
        for lab, (maj, hemi) in PARCELS.items():
            v = int(((l2a == lab) & les).sum())
            row[f"vox_{maj}_{hemi}"] = v
            row[f"share_{maj}_{hemi}"] = v / n
        row["share_L"] = sum(row[f"share_{m}_L"] for m in ["ACA", "MCA", "PCA", "VB"])
        row["share_R"] = sum(row[f"share_{m}_R"] for m in ["ACA", "MCA", "PCA", "VB"])
        # dominant major from L+R
        maj_share = {
            m: row[f"share_{m}_L"] + row[f"share_{m}_R"] for m in ["ACA", "MCA", "PCA", "VB"]
        }
        row["dominant_major"] = max(maj_share, key=maj_share.get)
        row["dominant_hemi"] = "L" if row["share_L"] >= row["share_R"] else "R"
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "atlas_laterality_subject_table.csv", index=False)

    # Does MCA dominance hold within each hemisphere?
    summary = {
        "n": int(len(df)),
        "mean_share_L": float(df["share_L"].mean()),
        "mean_share_R": float(df["share_R"].mean()),
        "dominant_hemi_counts": df["dominant_hemi"].value_counts().to_dict(),
        "dominant_major_counts": df["dominant_major"].value_counts().to_dict(),
        "mean_share_MCA_L": float(df["share_MCA_L"].mean()),
        "mean_share_MCA_R": float(df["share_MCA_R"].mean()),
        "mean_share_ACA_L": float(df["share_ACA_L"].mean()),
        "mean_share_ACA_R": float(df["share_ACA_R"].mean()),
        "frac_dom_MCA_when_dom_L": float(
            (df.loc[df["dominant_hemi"] == "L", "dominant_major"] == "MCA").mean()
        ),
        "frac_dom_MCA_when_dom_R": float(
            (df.loc[df["dominant_hemi"] == "R", "dominant_major"] == "MCA").mean()
        ),
        "interpretation": (
            "If MCA remains the modal dominant major territory within L-dominant and R-dominant "
            "subjects, collapse L/R→major is not an artifact of unilateral atlas error."
        ),
    }
    (OUT / "atlas_laterality_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# ATLAS laterality sensitivity",
        "",
        f"n={summary['n']}",
        f"- Mean lesion share L/R: {summary['mean_share_L']:.3f} / {summary['mean_share_R']:.3f}",
        f"- Dominant hemisphere counts: {summary['dominant_hemi_counts']}",
        f"- Dominant major counts: {summary['dominant_major_counts']}",
        f"- Mean MCA share L/R: {summary['mean_share_MCA_L']:.3f} / {summary['mean_share_MCA_R']:.3f}",
        f"- Frac MCA-dominant | L-dom / R-dom: {summary['frac_dom_MCA_when_dom_L']:.3f} / {summary['frac_dom_MCA_when_dom_R']:.3f}",
        "",
        summary["interpretation"],
        "",
    ]
    (OUT / "atlas_laterality_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
