"""
ISLES'22 Aim 4 stroke arm: native lesion masks -> MNI -> ArterialAtlas_level2.

Steps per subject:
  1. Affine-register DWI (mask space) -> MNI152 2 mm
  2. Apply transform to lesion mask (nearest)
  3. Overlay ArterialAtlas_level2 collapsed to ACA/MCA/PCA/VB

Same atlas rule as ATLAS: do NOT use TerritoryVoxels_BMM for stroke geography.

Usage:
  python analysis/aim4/isles_lesion_tv_pipeline.py --n 20
  python analysis/aim4/isles_lesion_tv_pipeline.py --resume
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import SimpleITK as sitk
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
ISLES = ROOT / "data/raw/isles/ISLES-2022/ISLES-2022"
LEVEL2_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/ArterialAtlas_level2.nii"
)
OUT = ROOT / "results/aim4/isles"
WORK = OUT / "_work"
MAJOR = ["ACA", "MCA", "PCA", "VB"]
LEVEL2_TO_MAJOR = {1: "ACA", 2: "ACA", 3: "MCA", 4: "MCA", 5: "PCA", 6: "PCA", 7: "VB", 8: "VB"}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def find_cases() -> list[dict]:
    rows = []
    for msk in sorted(ISLES.glob("derivatives/sub-*/ses-*/sub-*_msk.nii.gz")):
        sub = msk.parts[-3]  # sub-strokecaseXXXX under derivatives
        # path: derivatives/sub-X/ses-Y/file
        sub = [p for p in msk.parts if p.startswith("sub-")][0]
        ses = [p for p in msk.parts if p.startswith("ses-")][0]
        dwi = ISLES / sub / ses / "dwi" / f"{sub}_{ses}_dwi.nii.gz"
        flair = ISLES / sub / ses / "anat" / f"{sub}_{ses}_FLAIR.nii.gz"
        if not dwi.exists():
            continue
        rows.append({"subject_id": sub, "session": ses, "mask_path": str(msk), "dwi_path": str(dwi), "flair_path": str(flair) if flair.exists() else None})
    return rows


def load_participants() -> pd.DataFrame:
    p = ISLES / "participants.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    df["age_years"] = (
        df["age"].astype(str).str.extract(r"(\d+)", expand=False).astype(float)
    )
    df["sex_f"] = df["sex"].astype(str).str.upper().map({"F": 1.0, "M": 0.0})
    return df


def nib_to_sitk(path: Path) -> sitk.Image:
    img = sitk.ReadImage(str(path))
    img = sitk.Cast(img, sitk.sitkFloat32)
    return sitk.DICOMOrient(img, "RAS")


def affine_to_mni(moving: sitk.Image, mni: sitk.Image) -> sitk.Transform:
    fixed_n = sitk.Normalize(mni)
    moving_n = sitk.Normalize(moving)
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.2)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    init = sitk.CenteredTransformInitializer(
        fixed_n, moving_n, sitk.AffineTransform(3), sitk.CenteredTransformInitializerFilter.GEOMETRY
    )
    reg.SetInitialTransform(init, inPlace=False)
    return reg.Execute(fixed_n, moving_n)


def apply_tx(moving: sitk.Image, reference: sitk.Image, tx: sitk.Transform, interp) -> sitk.Image:
    return sitk.Resample(moving, reference, tx, interp, 0.0, moving.GetPixelID())


def mni_template_sitk() -> tuple[sitk.Image, nib.Nifti1Image]:
    WORK.mkdir(parents=True, exist_ok=True)
    # Match ADNI 2 mm MNI152
    mni_nib = datasets.load_mni152_template(resolution=2)
    tmp = WORK / "_mni152_2mm.nii.gz"
    nib.save(mni_nib, str(tmp))
    mni_sitk = sitk.Cast(sitk.DICOMOrient(sitk.ReadImage(str(tmp)), "RAS"), sitk.sitkFloat32)
    return mni_sitk, mni_nib


def collapse_level2(level2: np.ndarray) -> np.ndarray:
    out = np.zeros(level2.shape, dtype=np.int16)
    code = {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}
    for lab, maj in LEVEL2_TO_MAJOR.items():
        out[level2 == lab] = code[maj]
    return out


def tv_fractions(lesion: np.ndarray, tv: np.ndarray) -> dict:
    lesion = lesion.astype(bool)
    row = {"n_lesion_voxels": int(lesion.sum())}
    assigned = 0
    code = {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}
    for name in MAJOR:
        m = tv == code[name]
        v = int((lesion & m).sum())
        assigned += v
        row[f"lesion_vox_{name}"] = v
        row[f"tv_vox_{name}"] = int(m.sum())
        row[f"lesion_frac_of_lesion_{name}"] = v / row["n_lesion_voxels"] if row["n_lesion_voxels"] else 0.0
        row[f"lesion_burden_in_{name}"] = v / row[f"tv_vox_{name}"] if row[f"tv_vox_{name}"] else 0.0
    row["lesion_vox_unassigned"] = row["n_lesion_voxels"] - assigned
    row["lesion_frac_assigned_tv"] = assigned / row["n_lesion_voxels"] if row["n_lesion_voxels"] else 0.0
    shares = {n: row[f"lesion_frac_of_lesion_{n}"] for n in MAJOR}
    if row["n_lesion_voxels"] > 0 and any(shares.values()):
        dom = max(shares, key=shares.get)
        row["dominant_territory"] = dom
        row["dominant_share"] = shares[dom]
    else:
        row["dominant_territory"] = None
        row["dominant_share"] = 0.0
    return row


def apply_qc(row: dict) -> dict:
    reasons, flags = [], []
    n = int(row.get("n_lesion_voxels") or 0)
    frac = float(row.get("lesion_frac_assigned_tv") or 0.0)
    if n <= 0:
        reasons.append("FAIL_ZERO_LESION")
    if n > 0 and frac < 0.10:
        flags.append("FLAG_LOW_TV")
    dice = float(row.get("reg_brain_overlap_proxy") or 1.0)
    if dice < 0.05:
        reasons.append("FAIL_REGISTRATION")
    row["qc_reasons"] = ";".join(reasons)
    row["qc_flags"] = ";".join(flags)
    row["qc_primary"] = len(reasons) == 0
    row["qc_sensitivity_strict"] = len(reasons) == 0 and "FLAG_LOW_TV" not in flags
    return row


def process_one(case: dict, mni_sitk: sitk.Image, mni_nib: nib.Nifti1Image, tv: np.ndarray) -> dict:
    dwi = nib_to_sitk(Path(case["dwi_path"]))
    # DWI may be 4D — take first volume
    if dwi.GetDimension() == 4:
        dwi = dwi[:, :, :, 0]
    msk = sitk.ReadImage(case["mask_path"])
    msk = sitk.Cast(sitk.DICOMOrient(msk, "RAS"), sitk.sitkUInt8)
    if msk.GetDimension() == 4:
        msk = msk[:, :, :, 0]

    tx = affine_to_mni(dwi, mni_sitk)
    dwi_mni = apply_tx(dwi, mni_sitk, tx, sitk.sitkLinear)
    msk_mni = apply_tx(msk, mni_sitk, tx, sitk.sitkNearestNeighbor)

    # save mask for audit (optional lightweight)
    sub_dir = WORK / case["subject_id"]
    sub_dir.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(msk_mni, str(sub_dir / "lesion_mni.nii.gz"))

    # registration proxy: fraction of brain template voxels with DWI signal
    dwi_arr = sitk.GetArrayFromImage(dwi_mni)
    mni_arr = sitk.GetArrayFromImage(mni_sitk)
    brain = mni_arr > np.percentile(mni_arr[mni_arr > 0], 10) if np.any(mni_arr > 0) else mni_arr > 0
    signal = dwi_arr > np.percentile(dwi_arr, 50)
    overlap = float((brain & signal).sum() / max(int(brain.sum()), 1))

    lesion = sitk.GetArrayFromImage(msk_mni) > 0
    # SimpleITK array is z,y,x; nibabel tv is x,y,z — convert via nibabel path for consistency
    tmp_msk = sub_dir / "lesion_mni.nii.gz"
    lesion_nib = nib.load(str(tmp_msk))
    # resample to exact mni_nib grid if needed
    if lesion_nib.shape != mni_nib.shape:
        lesion_nib = image.resample_to_img(
            lesion_nib, mni_nib, interpolation="nearest", force_resample=True, copy_header=True
        )
    lesion = np.asanyarray(lesion_nib.dataobj) > 0
    if lesion.shape != tv.shape:
        raise RuntimeError(f"shape mismatch lesion{lesion.shape} tv{tv.shape}")

    fr = tv_fractions(lesion, tv)
    row = {
        **case,
        **fr,
        "reg_brain_overlap_proxy": overlap,
        "status": "ok",
        "space": "MNI152_2mm",
        "atlas_map": "ArterialAtlas_level2",
        "registration": "DWI_affine_to_MNI152_2mm",
    }
    return apply_qc(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=0, help="max subjects (0=all)")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    out_csv = OUT / "isles_lesion_tv_fractions.csv"
    cases = find_cases()
    if not cases:
        raise SystemExit(f"No ISLES masks under {ISLES}")
    if args.n and args.n > 0:
        cases = cases[: args.n]
    log(f"ISLES cases: {len(cases)}")

    done = {}
    if args.resume and out_csv.exists():
        prev = pd.read_csv(out_csv)
        for _, r in prev.iterrows():
            if str(r.get("status", "")) == "ok":
                done[str(r["subject_id"])] = r.to_dict()
        log(f"Resume: keeping {len(done)} completed")

    dem = load_participants()
    dem_map = dem.set_index("participant_id").to_dict("index") if len(dem) else {}

    log("Preparing MNI + ArterialAtlas_level2 ...")
    mni_sitk, mni_nib = mni_template_sitk()
    l2_img = nib.load(str(LEVEL2_PATH))
    l2_r = image.resample_to_img(
        l2_img, mni_nib, interpolation="nearest", force_resample=True, copy_header=True
    )
    l2 = np.rint(np.asanyarray(l2_r.dataobj, dtype=float)).astype(np.int32)
    if l2.ndim == 4:
        l2 = l2[..., 0]
    tv = collapse_level2(l2)
    tv_out = OUT / "ArterialAtlas_level2_major_on_MNI152_2mm.nii.gz"
    nib.save(nib.Nifti1Image(tv.astype(np.int16), mni_nib.affine), str(tv_out))
    log("TV major counts: " + ", ".join(f"{n}={(tv == i).sum()}" for i, n in enumerate(MAJOR, 1)))

    rows = list(done.values())
    todo = [c for c in cases if c["subject_id"] not in done]
    log(f"To process: {len(todo)}")

    for i, case in enumerate(todo, 1):
        sid = case["subject_id"]
        if i == 1 or i % 10 == 0 or i == len(todo):
            log(f"[{i}/{len(todo)}] {sid}")
        try:
            row = process_one(case, mni_sitk, mni_nib, tv)
            meta = dem_map.get(sid, {})
            row["age_years"] = meta.get("age_years")
            row["sex_f"] = meta.get("sex_f")
            row["weight"] = meta.get("weight")
        except Exception as e:
            row = {**case, "status": "error", "error": str(e), "qc_primary": False, "atlas_map": "ArterialAtlas_level2"}
            log(f"  ERROR {sid}: {e}")
        rows.append(row)
        if i % 10 == 0 or i == len(todo):
            pd.DataFrame(rows).to_csv(out_csv, index=False)

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    ok = df[df["status"] == "ok"].copy()
    primary = ok[ok["qc_primary"] == True]  # noqa: E712
    strict = ok[ok["qc_sensitivity_strict"] == True]  # noqa: E712

    def share_means(d):
        return {t: float(d[f"lesion_frac_of_lesion_{t}"].mean()) if len(d) else None for t in MAJOR}

    summary = {
        "dataset": "ISLES-2022",
        "doi": "10.5281/zenodo.7960856",
        "atlas_map": "ArterialAtlas_level2 (L/R collapsed to ACA/MCA/PCA/VB)",
        "registration": "DWI affine -> MNI152 2mm; nearest for mask",
        "n_cases": int(len(cases)),
        "n_ok": int(len(ok)),
        "n_error": int((df["status"] != "ok").sum()),
        "n_primary_qc": int(len(primary)),
        "n_strict_qc": int(len(strict)),
        "median_lesion_voxels_primary": float(primary["n_lesion_voxels"].median()) if len(primary) else None,
        "median_frac_assigned_tv_primary": float(primary["lesion_frac_assigned_tv"].median()) if len(primary) else None,
        "mean_lesion_share_primary": share_means(primary),
        "dominant_territory_counts_primary": {
            str(k): int(v) for k, v in primary["dominant_territory"].value_counts(dropna=False).items()
        }
        if len(primary)
        else {},
        "interpretation_locks": [
            "ISLES'22 is an independent stroke replication arm (not Liu construction).",
            "Use ArterialAtlas_level2 bilateral; do not use TerritoryVoxels_BMM.",
            "Native-space masks require DWI->MNI registration before territory overlay.",
            "Descriptive geography + MRA null; not molecular H3.",
        ],
    }
    (OUT / "isles_lesion_tv_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# ISLES'22 lesion x ArterialAtlas_level2",
        "",
        f"**n_ok / primary / strict:** {summary['n_ok']} / {summary['n_primary_qc']} / {summary['n_strict_qc']}",
        f"- Mean lesion share ACA/MCA/PCA/VB: {summary['mean_lesion_share_primary']}",
        f"- Dominant counts: {summary['dominant_territory_counts_primary']}",
        f"- Median assigned frac: {summary['median_frac_assigned_tv_primary']}",
        "",
        "## Locked interpretation",
        "- Independent stroke outcome arm for Aim 4 replication vs ATLAS.",
        "- Bilateral level2 majors only; TV_BMM NO-GO for stroke geography.",
        "- Do not claim molecular H3 from lesion shares alone.",
        "",
        f"Table: `{out_csv.relative_to(ROOT)}`",
        "",
    ]
    (OUT / "isles_lesion_tv_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
