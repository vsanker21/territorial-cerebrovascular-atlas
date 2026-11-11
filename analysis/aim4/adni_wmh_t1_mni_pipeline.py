"""
Hardened ADNI Aim-4 WMH pipeline: FLAIR → T1 (rigid) → MNI (affine) → TerritoryVoxels.

Uses SimpleITK registration (ANTs not required). No IDA spatial masks.

Usage:
  python analysis/aim4/adni_wmh_t1_mni_pipeline.py --n 15 --min-wmh 3.0
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
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
WORK = ROOT / "data/processed/adni/t1_mni"
TV_MAP = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def load_member_cache(cache: Path, builder):
    if cache.exists():
        if cache.suffix == ".parquet":
            return pd.read_parquet(cache)
        return pd.read_csv(cache)
    df = builder()
    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_parquet(cache, index=False)
    except Exception:
        df.to_csv(cache.with_suffix(".csv"), index=False)
    return df


def index_zip(zips: list[Path], name_filter, label: str) -> pd.DataFrame:
    rows = []
    for zp in zips:
        if not zp.exists():
            continue
        log(f"Indexing {label}: {zp.name}")
        with zipfile.ZipFile(zp) as zf:
            for n in zf.namelist():
                if n.endswith("/"):
                    continue
                nl = n.lower()
                if not name_filter(nl):
                    continue
                parts = n.replace("\\", "/").split("/")
                if len(parts) < 5 or parts[0] != "ADNI":
                    continue
                ptid, series, date = parts[1], parts[2], parts[3]
                m = re.search(r"/I(\d+)/", "/" + n.replace("\\", "/"))
                rows.append(
                    {
                        "zip": str(zp),
                        "member": n,
                        "ptid": ptid,
                        "series": series,
                        "acq_date": date[:10] if date else None,
                        "image_id": int(m.group(1)) if m else None,
                    }
                )
    df = pd.DataFrame(rows)
    log(f"{label}: {len(df)} files / {df['ptid'].nunique() if len(df) else 0} subjects")
    return df


def extract_members(zip_path: Path, members: list[str], dest: Path, alt_zips: list[Path] | None = None):
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    zips = [zip_path] + [z for z in (alt_zips or []) if z.resolve() != zip_path.resolve() and z.exists()]
    # verify members exist; fall back across zips
    extracted = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as zf:
            names = set(zf.namelist())
            present = [m for m in members if m in names]
            if len(present) < max(10, int(0.5 * len(members))):
                continue
            for m in present:
                out = dest / Path(m).name
                with zf.open(m) as src, open(out, "wb") as f:
                    shutil.copyfileobj(src, f)
                extracted += 1
            if extracted:
                return zp
    raise FileNotFoundError(f"Could not extract series from {[z.name for z in zips]} (n_members={len(members)})")


def dicom_to_sitk(dicom_dir: Path) -> sitk.Image:
    reader = sitk.ImageSeriesReader()
    sids = reader.GetGDCMSeriesIDs(str(dicom_dir))
    if not sids:
        raise RuntimeError(f"No DICOM series in {dicom_dir}")
    best, nbest = None, -1
    for sid in sids:
        files = reader.GetGDCMSeriesFileNames(str(dicom_dir), sid)
        if len(files) > nbest:
            best, nbest = files, len(files)
    reader.SetFileNames(best)
    img = reader.Execute()
    img = sitk.Cast(img, sitk.sitkFloat32)
    mask = sitk.OtsuThreshold(img, 0, 1, 200)
    try:
        img = sitk.N4BiasFieldCorrectionImageFilter().Execute(img, mask)
    except Exception:
        pass
    return sitk.DICOMOrient(img, "RAS")


def write_sitk(img: sitk.Image, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(img, str(path))


def rigid_register(fixed: sitk.Image, moving: sitk.Image) -> tuple[sitk.Image, sitk.Transform]:
    fixed_n = sitk.Normalize(fixed)
    moving_n = sitk.Normalize(moving)
    # multi-res rigid
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(50)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.2)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsGradientDescent(
        learningRate=1.0, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=10
    )
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    init = sitk.CenteredTransformInitializer(
        fixed_n, moving_n, sitk.Euler3DTransform(), sitk.CenteredTransformInitializerFilter.GEOMETRY
    )
    reg.SetInitialTransform(init, inPlace=False)
    tx = reg.Execute(fixed_n, moving_n)
    resampled = sitk.Resample(moving, fixed, tx, sitk.sitkLinear, 0.0, moving.GetPixelID())
    return resampled, tx


def affine_to_mni(moving: sitk.Image, mni_sitk: sitk.Image) -> tuple[sitk.Image, sitk.Transform]:
    fixed_n = sitk.Normalize(mni_sitk)
    moving_n = sitk.Normalize(moving)
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.15)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsGradientDescent(
        learningRate=1.0, numberOfIterations=120, convergenceMinimumValue=1e-6, convergenceWindowSize=12
    )
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([6, 4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([3, 2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    init = sitk.CenteredTransformInitializer(
        fixed_n, moving_n, sitk.AffineTransform(3), sitk.CenteredTransformInitializerFilter.GEOMETRY
    )
    reg.SetInitialTransform(init, inPlace=False)
    tx = reg.Execute(fixed_n, moving_n)
    resampled = sitk.Resample(moving, mni_sitk, tx, sitk.sitkLinear, 0.0, moving.GetPixelID())
    return resampled, tx


def apply_tx(moving: sitk.Image, reference: sitk.Image, tx: sitk.Transform, interp=sitk.sitkLinear) -> sitk.Image:
    return sitk.Resample(moving, reference, tx, interp, 0.0, moving.GetPixelID())


def mni_template_sitk() -> sitk.Image:
    mni = datasets.load_mni152_template(resolution=2)
    tmp = WORK / "_mni152_2mm.nii.gz"
    if not tmp.exists():
        nib.save(mni, str(tmp))
    img = sitk.ReadImage(str(tmp))
    return sitk.Cast(sitk.DICOMOrient(img, "RAS"), sitk.sitkFloat32)


def segment_wmh_flair_mni(flair_mni: sitk.Image, z_thr: float = 1.5, wm_pct: float = 99.0):
    """WMH on FLAIR already in MNI space using MNI WM∩brain priors."""
    tmp = WORK / "_tmp_flair_mni.nii.gz"
    write_sitk(flair_mni, tmp)
    mni = datasets.load_mni152_template(resolution=2)
    wm = datasets.load_mni152_wm_mask(threshold=0.5, resolution=2)
    brain = datasets.load_mni152_brain_mask(threshold=0.5, resolution=2)
    f_img = image.load_img(str(tmp))
    f_img = image.resample_to_img(f_img, mni, interpolation="continuous", force_resample=True, copy_header=True)
    f = np.asanyarray(f_img.dataobj, dtype=float)
    w = np.asanyarray(wm.dataobj, dtype=bool)
    b = np.asanyarray(brain.dataobj, dtype=bool)
    if w.shape != f.shape:
        w = (
            np.asanyarray(
                image.resample_to_img(wm, f_img, interpolation="nearest", force_resample=True, copy_header=True).dataobj
            )
            > 0
        )
        b = (
            np.asanyarray(
                image.resample_to_img(brain, f_img, interpolation="nearest", force_resample=True, copy_header=True).dataobj
            )
            > 0
        )
    mask = w & b & np.isfinite(f) & (f > 0)
    vals = f[mask]
    if vals.size < 500:
        raise RuntimeError(f"WM voxels too few after T1-MNI path: {vals.size}")
    thr = min(float(vals.mean() + z_thr * vals.std()), float(np.percentile(vals, wm_pct)))
    wmh = (f >= thr) & mask
    # opening
    arr = wmh.astype(np.uint8).transpose(2, 1, 0)
    simg = sitk.GetImageFromArray(arr)
    simg = sitk.BinaryMorphologicalOpening(simg, [1, 1, 1])
    wmh = sitk.GetArrayFromImage(simg).transpose(2, 1, 0).astype(bool)
    return f_img, wmh, thr, int(vals.size)


def tv_fractions(wmh: np.ndarray, flair_mni_nib: nib.Nifti1Image) -> dict:
    tv_img = nib.load(str(TV_PATH))
    tv_r = image.resample_to_img(
        tv_img, flair_mni_nib, interpolation="nearest", force_resample=True, copy_header=True
    )
    tv = np.rint(np.asanyarray(tv_r.dataobj, dtype=float)).astype(np.int32)
    if tv.ndim == 4:
        tv = tv[..., 0]
    row = {"n_wmh_voxels": int(wmh.sum())}
    assigned = 0
    for lab, name in TV_MAP.items():
        m = tv == lab
        v = int((wmh & m).sum())
        assigned += v
        row[f"wmh_vox_{name}"] = v
        row[f"tv_vox_{name}"] = int(m.sum())
        row[f"wmh_frac_of_wmh_{name}"] = v / row["n_wmh_voxels"] if row["n_wmh_voxels"] else 0.0
        row[f"wmh_burden_in_{name}"] = v / row[f"tv_vox_{name}"] if row[f"tv_vox_{name}"] else 0.0
    row["wmh_vox_unassigned"] = row["n_wmh_voxels"] - assigned
    row["wmh_frac_assigned_tv"] = assigned / row["n_wmh_voxels"] if row["n_wmh_voxels"] else 0.0
    return row


def pick_series(df: pd.DataFrame, prefer_substr: list[str] | None = None) -> pd.DataFrame:
    """Choose one series with most files; prefer name substrings."""
    if prefer_substr:
        for s in prefer_substr:
            hit = df[df["series"].astype(str).str.contains(s, case=False, na=False)]
            if len(hit):
                df = hit
                break
    counts = df.groupby("series").size().sort_values(ascending=False)
    series = counts.index[0]
    return df[df["series"] == series]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=15, help="max subjects (use large n for full paired set)")
    ap.add_argument("--min-wmh", type=float, default=3.0, help="min UCD TOTAL_WMH cc for pilot")
    ap.add_argument("--z-thr", type=float, default=1.5)
    ap.add_argument("--resume", action="store_true", help="skip ptids already status=ok with wmh_mni.nii.gz")
    ap.add_argument("--run-qc", action="store_true", help="run qc_adni_wmh_t1_mni.py after finish")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    out_csv = OUT / "adni_wmh_t1_mni_scaled.csv"

    inv = pd.read_csv(OUT / "adni_aim4_cohort_inventory.csv")
    inv["ptid"] = inv["subject_id"].astype(str)
    inv = inv[inv["TOTAL_WMH"].notna() & (inv["TOTAL_WMH"] >= args.min_wmh)].copy()

    flair_zips = [ADNI / "2D T2 FLAIR_dataset.zip", ADNI / "2D T2 FLAIR.zip"]
    # Index smaller T1 pack first; add large pack only if pairing is sparse
    t1_small = ADNI / "Clinical_T1w_Imaging_Cohort_MRI.zip"
    t1_large = ADNI / "Clinical_T1w_Imaging_Cohort_MRI_1.zip"

    flair_m = load_member_cache(
        WORK / "flair_zip_members.parquet",
        lambda: index_zip(flair_zips, lambda n: "flair" in n, "FLAIR"),
    )
    # also reuse pilot cache if present
    alt = ROOT / "data/processed/adni/pilot/flair_zip_members.parquet"
    if not (WORK / "flair_zip_members.parquet").exists() and alt.exists():
        shutil.copy2(alt, WORK / "flair_zip_members.parquet")
        flair_m = pd.read_parquet(WORK / "flair_zip_members.parquet")

    def _t1_filt(n: str) -> bool:
        return ("mprage" in n or "mp-rage" in n or "spgr" in n or "ir-spgr" in n) and "flair" not in n

    t1_m = load_member_cache(
        WORK / "t1_zip_members_small.parquet",
        lambda: index_zip([t1_small], _t1_filt, "T1_small"),
    )
    common = sorted(set(inv["ptid"]) & set(flair_m["ptid"]) & set(t1_m["ptid"]))
    if len(common) < min(args.n, 267) and t1_large.exists():
        log(f"Only {len(common)} pairs in small T1 zip — indexing large T1 pack")
        t1_big = load_member_cache(
            WORK / "t1_zip_members_large.parquet",
            lambda: index_zip([t1_large], _t1_filt, "T1_large"),
        )
        t1_m = pd.concat([t1_m, t1_big], ignore_index=True)
        common = sorted(set(inv["ptid"]) & set(flair_m["ptid"]) & set(t1_m["ptid"]))

    log(f"Paired FLAIR+T1+UCD(>= {args.min_wmh} cc): {len(common)} subjects")
    inv_c = inv[inv["ptid"].isin(common)].sort_values("TOTAL_WMH", ascending=False)
    pilot = inv_c["ptid"].head(args.n).tolist()
    if len(pilot) < args.n:
        log(f"WARNING: only {len(pilot)} paired subjects available")

    prior_rows: list[dict] = []
    done = set()
    if args.resume and out_csv.exists():
        prior = pd.read_csv(out_csv)
        for _, r in prior.iterrows():
            ptid = str(r["ptid"])
            if str(r.get("status", "")) == "ok" and (WORK / ptid / "wmh_mni.nii.gz").exists():
                prior_rows.append(r.to_dict())
                done.add(ptid)
        log(f"Resume: keeping {len(done)} completed subjects")
    pilot = [p for p in pilot if p not in done]
    log(f"To process n={len(pilot)} (plus {len(done)} resumed); first: {pilot[:5]}")

    mni_sitk = mni_template_sitk()
    rows = list(prior_rows)
    for i, ptid in enumerate(pilot, 1):
        log(f"[{i}/{len(pilot)}] {ptid}")
        try:
            fsub = pick_series(
                flair_m[flair_m["ptid"] == ptid],
                prefer_substr=["3D_FLAIR", "FLAIR", "2D_FLAIR"],
            )
            tsub = pick_series(
                t1_m[t1_m["ptid"] == ptid],
                prefer_substr=["MPRAGE", "MP-RAGE", "IR-SPGR", "SPGR"],
            )
            fzip, tzip = Path(fsub.iloc[0]["zip"]), Path(tsub.iloc[0]["zip"])
            fdcm, tdcm = WORK / ptid / "flair_dcm", WORK / ptid / "t1_dcm"
            extract_members(fzip, fsub["member"].tolist(), fdcm, alt_zips=flair_zips)
            extract_members(tzip, tsub["member"].tolist(), tdcm, alt_zips=[t1_small, t1_large])

            flair = dicom_to_sitk(fdcm)
            t1 = dicom_to_sitk(tdcm)
            write_sitk(flair, WORK / ptid / "flair_n4.nii.gz")
            write_sitk(t1, WORK / ptid / "t1_n4.nii.gz")

            flair_t1, tx_f2t = rigid_register(t1, flair)
            write_sitk(flair_t1, WORK / ptid / "flair_in_t1.nii.gz")

            t1_mni, tx_t2m = affine_to_mni(t1, mni_sitk)
            flair_mni = apply_tx(flair_t1, mni_sitk, tx_t2m)
            write_sitk(t1_mni, WORK / ptid / "t1_mni.nii.gz")
            write_sitk(flair_mni, WORK / ptid / "flair_mni.nii.gz")

            f_nib, wmh, thr, n_wm = segment_wmh_flair_mni(flair_mni, z_thr=args.z_thr)
            nib.save(nib.Nifti1Image(wmh.astype(np.uint8), f_nib.affine), str(WORK / ptid / "wmh_mni.nii.gz"))
            fr = tv_fractions(wmh, f_nib)
            meta = inv[inv["ptid"] == ptid].iloc[0]
            fr.update(
                {
                    "ptid": ptid,
                    "RID": int(meta["RID"]) if pd.notna(meta["RID"]) else None,
                    "TOTAL_WMH_ucd": float(meta["TOTAL_WMH"]),
                    "flair_series": fsub.iloc[0]["series"],
                    "t1_series": tsub.iloc[0]["series"],
                    "seg_thr": thr,
                    "seg_n_wm": n_wm,
                    "method": "FLAIR->T1_rigid + T1->MNI_affine + WM_zscore",
                    "status": "ok",
                }
            )
            # crude QC: Dice-like overlap of brain after reg via nonzero
            fr["qc_flair_mni_nonzero"] = int(np.count_nonzero(sitk.GetArrayFromImage(flair_mni) > 0))
            fr["qc_t1_mni_nonzero"] = int(np.count_nonzero(sitk.GetArrayFromImage(t1_mni) > 0))
            rows.append(fr)
            log(
                f"  OK wmh={fr['n_wmh_voxels']} assigned_frac={fr['wmh_frac_assigned_tv']:.2f} "
                f"ACA/MCA/PCA/VB={fr['wmh_vox_ACA']}/{fr['wmh_vox_MCA']}/{fr['wmh_vox_PCA']}/{fr['wmh_vox_VB']}"
            )
        except Exception as e:
            log(f"  FAIL {ptid}: {e}")
            rows.append({"ptid": ptid, "status": f"error: {e}"})
        finally:
            shutil.rmtree(WORK / ptid / "flair_dcm", ignore_errors=True)
            shutil.rmtree(WORK / ptid / "t1_dcm", ignore_errors=True)
            # checkpoint after each subject
            pd.DataFrame(rows).drop_duplicates("ptid", keep="last").to_csv(out_csv, index=False)

    df = pd.DataFrame(rows)
    # dedupe by ptid keeping last
    if len(df) and "ptid" in df.columns:
        df = df.drop_duplicates("ptid", keep="last")
    df.to_csv(out_csv, index=False)
    ok = df[df["status"] == "ok"] if "status" in df.columns else df
    summary = {
        "n_requested": args.n,
        "n_rows": int(len(df)),
        "n_ok": int(len(ok)),
        "n_newly_processed": int(len(pilot)),
        "n_resumed": int(len(done)),
        "min_ucd_wmh": args.min_wmh,
        "method": "FLAIR->T1 rigid -> T1->MNI affine (SimpleITK) -> WM-prior hyperintensity",
        "median_wmh_voxels": float(ok["n_wmh_voxels"].median()) if len(ok) and "n_wmh_voxels" in ok else None,
        "median_frac_assigned_tv": float(ok["wmh_frac_assigned_tv"].median())
        if len(ok) and "wmh_frac_assigned_tv" in ok
        else None,
        "mean_frac_ACA": float(ok["wmh_frac_of_wmh_ACA"].mean()) if len(ok) and "wmh_frac_of_wmh_ACA" in ok else None,
        "mean_frac_MCA": float(ok["wmh_frac_of_wmh_MCA"].mean()) if len(ok) and "wmh_frac_of_wmh_MCA" in ok else None,
        "mean_frac_PCA": float(ok["wmh_frac_of_wmh_PCA"].mean()) if len(ok) and "wmh_frac_of_wmh_PCA" in ok else None,
        "mean_frac_VB": float(ok["wmh_frac_of_wmh_VB"].mean()) if len(ok) and "wmh_frac_of_wmh_VB" in ok else None,
    }
    corr = None
    if len(ok) >= 4 and "TOTAL_WMH_ucd" in ok.columns:
        from scipy.stats import spearmanr

        rho, p = spearmanr(ok["TOTAL_WMH_ucd"], ok["n_wmh_voxels"])
        corr = {"spearman_rho": float(rho), "p": float(p)}
        summary["ucd_vs_local_wmh_voxels"] = corr

    (OUT / "adni_wmh_t1_mni_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = [
        "# ADNI WMH — T1-based registration (scaled)",
        "",
        f"**n_ok / n_rows:** {summary['n_ok']} / {summary['n_rows']} "
        f"(new={summary['n_newly_processed']}, resumed={summary['n_resumed']})",
        f"- Median local WMH voxels: {summary.get('median_wmh_voxels')}",
        f"- Median fraction of WMH in labeled TV: {summary.get('median_frac_assigned_tv')}",
        f"- Mean WMH share ACA/MCA/PCA/VB: "
        f"{summary.get('mean_frac_ACA')}/{summary.get('mean_frac_MCA')}/"
        f"{summary.get('mean_frac_PCA')}/{summary.get('mean_frac_VB')}",
        f"- UCD vs local Spearman: `{corr}`",
        "",
        "Run `python analysis/aim4/qc_adni_wmh_t1_mni.py` then "
        "`python analysis/aim4/aim4_nest_wmh_mra.py` for nested Aim 4 tables.",
        "",
    ]
    (OUT / "adni_wmh_t1_mni_report.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(summary, indent=2))

    if args.run_qc:
        import subprocess
        import sys

        log("Running QC…")
        subprocess.check_call([sys.executable, str(ROOT / "analysis/aim4/qc_adni_wmh_t1_mni.py")])
        log("Running Aim 4 nest (TV WMH + MRA geometry)…")
        subprocess.check_call([sys.executable, str(ROOT / "analysis/aim4/aim4_nest_wmh_mra.py")])


if __name__ == "__main__":
    main()
