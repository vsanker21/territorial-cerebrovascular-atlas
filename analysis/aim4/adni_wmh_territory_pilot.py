"""
ADNI Aim-4 pilot: local FLAIR WMH segmentation → MNI → TerritoryVoxels fractions.

Scientific locks:
- No IDA spatial WMH masks; segmentation is local and exploratory.
- Method: N4 bias-correct FLAIR; affine to MNI152 via nilearn; WM prior from
  MNI tissue probability; hyperintensity = FLAIR z>z_thr within WM ∩ brain.
- Pilot N subjects only (default 5) extracted selectively from zip packs.

Usage:
  python analysis/aim4/adni_wmh_territory_pilot.py --n 5
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import tempfile
import time
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import SimpleITK as sitk
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
ADNI = ROOT / "data/raw/ADNI"
TV_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/TerritoryVoxels_BMM.nii"
)
OUT = ROOT / "results/aim4/adni"
WORK = ROOT / "data/processed/adni/pilot"
TV_MAP = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def find_flair_zips():
    return [
        p
        for p in [
            ADNI / "2D T2 FLAIR_dataset.zip",
            ADNI / "2D T2 FLAIR.zip",
        ]
        if p.exists()
    ]


def cache_flair_members(zips: list[Path], cache: Path) -> pd.DataFrame:
    if cache.exists():
        return pd.read_parquet(cache) if cache.suffix == ".parquet" else pd.read_csv(cache)
    rows = []
    for zp in zips:
        log(f"Indexing {zp.name} …")
        with zipfile.ZipFile(zp) as zf:
            for n in zf.namelist():
                if n.endswith("/") or "flair" not in n.lower():
                    continue
                # ADNI/PTID/Series/date/Ixxxx/...
                parts = n.replace("\\", "/").split("/")
                if len(parts) < 5 or parts[0] != "ADNI":
                    continue
                ptid = parts[1]
                series = parts[2]
                m = re.search(r"/I(\d+)/", "/" + n.replace("\\", "/"))
                rows.append(
                    {
                        "zip": str(zp),
                        "member": n,
                        "ptid": ptid,
                        "series": series,
                        "image_id": int(m.group(1)) if m else None,
                    }
                )
    df = pd.DataFrame(rows)
    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(cache, index=False)
    except Exception:
        df.to_csv(cache.with_suffix(".csv"), index=False)
    log(f"Cached {len(df)} FLAIR members / {df['ptid'].nunique()} subjects")
    return df


def extract_subject_dicoms(zip_path: Path, members: list[str], dest: Path):
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        for m in members:
            target = dest / Path(m).name
            with zf.open(m) as src, open(target, "wb") as out:
                shutil.copyfileobj(src, out)


def dicom_dir_to_nifti(dicom_dir: Path, out_nii: Path) -> Path:
    reader = sitk.ImageSeriesReader()
    series_ids = reader.GetGDCMSeriesIDs(str(dicom_dir))
    if not series_ids:
        raise RuntimeError(f"No DICOM series in {dicom_dir}")
    # pick largest series
    best, best_n = None, -1
    for sid in series_ids:
        files = reader.GetGDCMSeriesFileNames(str(dicom_dir), sid)
        if len(files) > best_n:
            best, best_n = files, len(files)
    reader.SetFileNames(best)
    img = reader.Execute()
    # N4
    mask = sitk.OtsuThreshold(img, 0, 1, 200)
    img = sitk.Cast(img, sitk.sitkFloat32)
    corrector = sitk.N4BiasFieldCorrectionImageFilter()
    try:
        img = corrector.Execute(img, mask)
    except Exception:
        pass
    out_nii.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(img, str(out_nii))
    return out_nii


def flair_wmh_to_mni(flair_nii: Path, z_thr: float = 2.0, wm_pct: float = 99.0):
    """Exploratory WMH: FLAIR hyperintensity within MNI WM prior after affine to MNI.

    Uses max(mean + z_thr*sd, wm_pct percentile) inside WM∩brain. Affine-only
    FLAIR→MNI is coarse; treat as pilot QC, not confirmatory lesion maps.
    """
    mni = datasets.load_mni152_template(resolution=2)
    wm = datasets.load_mni152_wm_mask(threshold=0.5, resolution=2)
    brain_m = datasets.load_mni152_brain_mask(threshold=0.5, resolution=2)
    # Reorient to RAS before nilearn resample
    img_sitk = sitk.ReadImage(str(flair_nii))
    img_sitk = sitk.DICOMOrient(img_sitk, "RAS")
    tmp = flair_nii.with_name(flair_nii.stem + "_ras.nii.gz")
    sitk.WriteImage(img_sitk, str(tmp))
    flair_mni = image.resample_to_img(
        image.load_img(str(tmp)),
        mni,
        interpolation="continuous",
        force_resample=True,
        copy_header=True,
    )
    f = np.asanyarray(flair_mni.dataobj, dtype=float)
    w = np.asanyarray(wm.dataobj, dtype=bool)
    b = np.asanyarray(brain_m.dataobj, dtype=bool)
    if w.shape != f.shape:
        w = np.asanyarray(
            image.resample_to_img(wm, flair_mni, interpolation="nearest", force_resample=True, copy_header=True).dataobj
        ) > 0
        b = np.asanyarray(
            image.resample_to_img(brain_m, flair_mni, interpolation="nearest", force_resample=True, copy_header=True).dataobj
        ) > 0
    wm_mask = w & b & np.isfinite(f) & (f > 0)
    vals = f[wm_mask]
    if vals.size < 100:
        raise RuntimeError(f"Too few WM voxels after registration (n={vals.size})")
    mu, sd = float(vals.mean()), float(vals.std() + 1e-6)
    thr_z = mu + z_thr * sd
    thr_p = float(np.percentile(vals, wm_pct))
    thr = min(thr_z, thr_p)  # more sensitive of the two
    wmh = (f >= thr) & wm_mask
    # remove tiny speckles
    wmh_img = sitk.GetImageFromArray(wmh.astype(np.uint8).transpose(2, 1, 0))
    wmh_img = sitk.BinaryMorphologicalOpening(wmh_img, [1, 1, 1])
    wmh = sitk.GetArrayFromImage(wmh_img).transpose(2, 1, 0).astype(bool)
    return flair_mni, wmh, {"thr": thr, "thr_z": thr_z, "thr_p": thr_p, "n_wm": int(vals.size)}


def territory_fractions(wmh: np.ndarray, tv_img: nib.Nifti1Image):
    tv = np.rint(np.asanyarray(tv_img.dataobj, dtype=float)).astype(np.int32)
    if tv.ndim == 4:
        tv = tv[..., 0]
    # resample TV to wmh grid if needed
    wmh_img = nib.Nifti1Image(wmh.astype(np.uint8), tv_img.affine)
    # assume both MNI 2mm-ish — resample TV to wmh shape via nilearn using flair_mni affine
    return tv, wmh


def compute_tv_fracs(wmh: np.ndarray, flair_mni: nib.Nifti1Image):
    tv_img = nib.load(str(TV_PATH))
    tv_r = image.resample_to_img(tv_img, flair_mni, interpolation="nearest")
    tv = np.rint(np.asanyarray(tv_r.dataobj, dtype=float)).astype(np.int32)
    if tv.ndim == 4:
        tv = tv[..., 0]
    if wmh.shape != tv.shape:
        raise RuntimeError(f"shape mismatch wmh{wmh.shape} tv{tv.shape}")
    row = {"n_wmh_voxels": int(wmh.sum())}
    for lab, name in TV_MAP.items():
        m = tv == lab
        row[f"wmh_vox_{name}"] = int((wmh & m).sum())
        row[f"tv_vox_{name}"] = int(m.sum())
        row[f"wmh_frac_of_wmh_{name}"] = (
            row[f"wmh_vox_{name}"] / row["n_wmh_voxels"] if row["n_wmh_voxels"] else 0.0
        )
        row[f"wmh_burden_in_{name}"] = (
            row[f"wmh_vox_{name}"] / row[f"tv_vox_{name}"] if row[f"tv_vox_{name}"] else 0.0
        )
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5, help="pilot subjects")
    ap.add_argument("--z-thr", type=float, default=2.0)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    inv = pd.read_csv(OUT / "adni_aim4_cohort_inventory.csv")
    inv = inv[inv["TOTAL_WMH"].notna()].copy()
    # prefer subjects present in FLAIR zip index
    zips = find_flair_zips()
    if not zips:
        raise SystemExit("No FLAIR zips found")
    cache = WORK / "flair_zip_members.parquet"
    members = cache_flair_members(zips, cache)
    # pick pilot ptids that appear in both inventory and zip
    inv["ptid"] = inv["subject_id"].astype(str)
    common = sorted(set(inv["ptid"]) & set(members["ptid"].astype(str)))
    if not common:
        raise SystemExit("No overlapping PTIDs between inventory and FLAIR zip index")
    pilot_ids = common[: args.n]
    log(f"Pilot PTIDs: {pilot_ids}")

    rows = []
    for ptid in pilot_ids:
        log(f"=== {ptid} ===")
        sub_m = members[members["ptid"] == ptid]
        # one series: most members
        series = sub_m.groupby("series").size().sort_values(ascending=False).index[0]
        sel = sub_m[sub_m["series"] == series]
        zip_path = Path(sel.iloc[0]["zip"])
        mems = sel["member"].tolist()
        # limit dicoms
        dcm_dir = WORK / ptid / "flair_dcm"
        if dcm_dir.exists():
            shutil.rmtree(dcm_dir)
        log(f"Extracting {len(mems)} DICOM files from {zip_path.name} / {series}")
        extract_subject_dicoms(zip_path, mems, dcm_dir)
        nii = WORK / ptid / "flair_n4.nii.gz"
        try:
            dicom_dir_to_nifti(dcm_dir, nii)
            flair_mni, wmh, thr_info = flair_wmh_to_mni(nii, z_thr=args.z_thr)
            wmh_path = WORK / ptid / "wmh_mni.nii.gz"
            nib.save(nib.Nifti1Image(wmh.astype(np.uint8), flair_mni.affine), str(wmh_path))
            fr = compute_tv_fracs(wmh, flair_mni)
            meta = inv[inv["ptid"] == ptid].iloc[0]
            fr.update(
                {
                    "ptid": ptid,
                    "RID": int(meta["RID"]) if pd.notna(meta["RID"]) else None,
                    "TOTAL_WMH_ucd": float(meta["TOTAL_WMH"]) if pd.notna(meta["TOTAL_WMH"]) else None,
                    "series": series,
                    "z_thr": args.z_thr,
                    **{f"seg_{k}": v for k, v in thr_info.items()},
                    "status": "ok",
                }
            )
            rows.append(fr)
            log(f"OK {ptid} wmh_vox={fr['n_wmh_voxels']}")
        except Exception as e:
            log(f"FAIL {ptid}: {e}")
            rows.append({"ptid": ptid, "status": f"error: {e}"})
        # free dicoms to save space
        shutil.rmtree(dcm_dir, ignore_errors=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "adni_wmh_territory_pilot.csv", index=False)
    summary = {
        "n_requested": args.n,
        "n_ok": int((df["status"] == "ok").sum()) if "status" in df.columns else 0,
        "method": "N4 FLAIR + affine MNI + WM-prior z-score hyperintensity (exploratory)",
        "ida_spatial_masks": False,
        "tv_atlas": str(TV_PATH.relative_to(ROOT)),
        "note": "Pilot only; not confirmatory Aim 4 until QC + fuller n + optional T1-based WM segmentation",
    }
    (OUT / "adni_wmh_territory_pilot_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# ADNI WMH → TerritoryVoxels pilot (local segmentation)",
        "",
        f"**Verdict:** `PILOT_LOCAL_WMH__NO_IDA_MASKS`",
        "",
        f"- Method: {summary['method']}",
        f"- OK: {summary['n_ok']} / {summary['n_requested']}",
        "",
        "```",
        df.to_string(index=False),
        "```",
        "",
        "## Locked",
        "- Not using IDA spatial WMH (unavailable).",
        "- Territory fractions are exploratory until scaled + QC'd.",
        "- UCD TOTAL_WMH remains the validated volume endpoint for full-cohort analyses.",
        "",
    ]
    (OUT / "adni_wmh_territory_pilot_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
