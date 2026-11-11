#!/usr/bin/env python3
"""Selective ISLES'24 acquisition (n=15) + Track C angio/lesion pilot.

Uses SimpleITK affine DWI→MNI (same pattern as ISLES'22 Aim 4) and
ArterialAtlas_level2 majors — NOT TerritoryVoxels_BMM.
"""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import SimpleITK as sitk
from nilearn import datasets, image

ROOT = Path(r"F:\Cerebrovascular SCS")
DATA = ROOT / "data/raw/isles24"
CASES = DATA / "cases"
OUT = ROOT / "results/nature_gates/isles24_trackC"
OUT.mkdir(parents=True, exist_ok=True)
WORK = OUT / "_work"
WORK.mkdir(parents=True, exist_ok=True)

LEVEL2_PATH = ROOT / (
    "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/"
    "Atlas_182_MNI152/ArterialAtlas_level2.nii"
)
MAJOR_ON_MNI = ROOT / "results/aim4/isles/ArterialAtlas_level2_major_on_MNI152_2mm.nii.gz"
LEVEL2_TO_MAJOR = {1: "ACA", 2: "ACA", 3: "MCA", 4: "MCA", 5: "PCA", 6: "PCA", 7: "VB", 8: "VB"}
N_CASES = 15


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def write_nii_bytes(blob: dict, dest: Path) -> Path | None:
    if not isinstance(blob, dict) or blob.get("bytes") is None:
        return None
    raw = blob["bytes"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if raw[:2] == b"\x1f\x8b":
        data = gzip.decompress(raw)
        out = dest if str(dest).endswith(".nii") else Path(str(dest) + ".nii")
        out.write_bytes(data)
        return out
    dest.write_bytes(raw)
    return dest


def download_cases(n: int = N_CASES) -> list[Path]:
    from huggingface_hub import hf_hub_download
    import pyarrow.parquet as pq

    CASES.mkdir(parents=True, exist_ok=True)
    extracted = []
    for i in range(n):
        fname = f"data/train-{i:05d}-of-00149.parquet"
        sid_guess = f"case_{i:05d}"
        # Skip re-extract if phenotype exists with files
        existing = list(CASES.glob("sub-stroke*/phenotype.json"))
        # Prefer deterministic extract by parquet index
        log(f"HF download {fname}")
        local = hf_hub_download(
            "hugging-science/isles24-stroke",
            fname,
            repo_type="dataset",
            local_dir=str(DATA),
        )
        row = pq.read_table(local).to_pydict()
        sid = row["subject_id"][0]
        cdir = CASES / sid
        if (cdir / "phenotype.json").exists() and (cdir / f"{sid}_lesion_mask.nii").exists():
            log(f"  skip extract (exists): {sid}")
            extracted.append(cdir)
            continue
        cdir.mkdir(parents=True, exist_ok=True)
        meta = {
            "subject_id": sid,
            "age": row.get("age", [None])[0],
            "sex": row.get("sex", [None])[0],
            "nihss_admission": row.get("nihss_admission", [None])[0],
            "mrs_admission": row.get("mrs_admission", [None])[0],
            "mrs_3month": row.get("mrs_3month", [None])[0],
            "files": {},
        }
        for field in ("cta", "dwi", "lesion_mask", "lvo_mask", "cow_segmentation"):
            dest = cdir / f"{sid}_{field}.nii"
            p = write_nii_bytes(row[field][0], dest)
            meta["files"][field] = str(p) if p else None
            if p:
                log(f"  {sid} {field}: {p.stat().st_size/1e6:.1f} MB")
        (cdir / "phenotype.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
        extracted.append(cdir)
    return extracted


def nib_to_sitk(path: Path) -> sitk.Image:
    img = sitk.Cast(sitk.ReadImage(str(path)), sitk.sitkFloat32)
    return sitk.DICOMOrient(img, "RAS")


def affine_to_mni(moving: sitk.Image, mni: sitk.Image) -> sitk.Transform:
    fixed_n = sitk.Normalize(mni)
    moving_n = sitk.Normalize(moving)
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.2)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=100,
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10,
    )
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    init = sitk.CenteredTransformInitializer(
        fixed_n,
        moving_n,
        sitk.AffineTransform(3),
        sitk.CenteredTransformInitializerFilter.GEOMETRY,
    )
    reg.SetInitialTransform(init, inPlace=False)
    return reg.Execute(fixed_n, moving_n)


def apply_tx(moving: sitk.Image, reference: sitk.Image, tx: sitk.Transform, interp) -> sitk.Image:
    return sitk.Resample(moving, reference, tx, interp, 0.0, moving.GetPixelID())


def mni_template_sitk() -> tuple[sitk.Image, nib.Nifti1Image]:
    mni_nib = datasets.load_mni152_template(resolution=2)
    tmp = WORK / "_mni152_2mm.nii.gz"
    if not tmp.exists():
        nib.save(mni_nib, str(tmp))
    mni_sitk = sitk.Cast(sitk.DICOMOrient(sitk.ReadImage(str(tmp)), "RAS"), sitk.sitkFloat32)
    return mni_sitk, nib.load(str(tmp))


def collapse_level2(level2: np.ndarray) -> np.ndarray:
    out = np.zeros(level2.shape, dtype=np.int16)
    code = {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}
    for lab, maj in LEVEL2_TO_MAJOR.items():
        out[level2 == lab] = code[maj]
    return out


def ensure_major_atlas(mni_nib: nib.Nifti1Image) -> tuple[np.ndarray, np.ndarray]:
    if MAJOR_ON_MNI.exists():
        img = nib.load(str(MAJOR_ON_MNI))
        return np.asanyarray(img.dataobj).astype(np.int16), img.affine
    l2 = nib.load(str(LEVEL2_PATH))
    l2_r = image.resample_to_img(
        l2, mni_nib, interpolation="nearest", force_resample=True, copy_header=True
    )
    major = collapse_level2(np.asanyarray(l2_r.dataobj))
    MAJOR_ON_MNI.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(major.astype(np.int16), l2_r.affine), str(MAJOR_ON_MNI))
    return major, l2_r.affine


def register_lesion(dwi_path: Path, lesion_path: Path, out_dir: Path, mni_sitk, mni_nib, major, major_aff) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        dwi = nib_to_sitk(dwi_path)
        lesion = nib_to_sitk(lesion_path)
        tx = affine_to_mni(dwi, mni_sitk)
        les_mni_sitk = apply_tx(lesion, mni_sitk, tx, sitk.sitkNearestNeighbor)
        out_les = out_dir / "lesion_mni.nii.gz"
        sitk.WriteImage(les_mni_sitk, str(out_les))
        les_nib = nib.load(str(out_les))
        # Resample to atlas grid if needed
        ref = nib.Nifti1Image(major.astype(np.float32), major_aff)
        les_on = image.resample_to_img(
            les_nib, ref, interpolation="nearest", force_resample=True, copy_header=True
        )
        les = np.asanyarray(les_on.dataobj) > 0
        names = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
        total = int(les.sum())
        shares = {
            "status": "ok",
            "method": "SimpleITK_affine_DWI_to_MNI_nearest_lesion",
            "atlas": "ArterialAtlas_level2_major",
            "lesion_vox_total": total,
            "lesion_mni": str(out_les),
        }
        for code, name in names.items():
            vox = int(((major == code) & les).sum())
            shares[f"lesion_vox_{name}"] = vox
            shares[f"lesion_share_{name}"] = float(vox / total) if total else np.nan
        return shares
    except Exception as e:
        return {"status": "reg_failed", "error": str(e)}


def mask_volume_mm3(path: Path | None) -> float | None:
    if path is None or not Path(path).exists():
        return None
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    vox = float(np.prod(img.header.get_zooms()[:3]))
    return float((data > 0).sum() * vox)


def lvo_present(path: Path | None) -> int | None:
    if path is None or not Path(path).exists():
        return None
    data = np.asanyarray(nib.load(str(path)).dataobj)
    return int((data > 0).sum() > 0)


def analyze_cases(case_dirs: list[Path]) -> pd.DataFrame:
    log("Preparing MNI + ArterialAtlas_level2 major")
    mni_sitk, mni_nib = mni_template_sitk()
    major, major_aff = ensure_major_atlas(mni_nib)
    rows = []
    for cdir in case_dirs:
        meta = json.loads((cdir / "phenotype.json").read_text(encoding="utf-8"))
        files = meta["files"]
        dwi = Path(files["dwi"]) if files.get("dwi") else None
        lesion = Path(files["lesion_mask"]) if files.get("lesion_mask") else None
        lvo = Path(files["lvo_mask"]) if files.get("lvo_mask") else None
        cow = Path(files["cow_segmentation"]) if files.get("cow_segmentation") else None
        row = {
            "subject_id": meta["subject_id"],
            "age": meta.get("age"),
            "sex": meta.get("sex"),
            "nihss_admission": meta.get("nihss_admission"),
            "mrs_3month": meta.get("mrs_3month"),
            "lvo_present": lvo_present(lvo),
            "lvo_volume_mm3": mask_volume_mm3(lvo),
            "cow_volume_mm3": mask_volume_mm3(cow),
            "lesion_volume_native_mm3": mask_volume_mm3(lesion),
        }
        if dwi and lesion and dwi.exists() and lesion.exists():
            log(f"Register {meta['subject_id']}")
            row.update(register_lesion(dwi, lesion, cdir / "mni", mni_sitk, mni_nib, major, major_aff))
        else:
            row["status"] = "missing_dwi_or_lesion"
        rows.append(row)
    return pd.DataFrame(rows)


def nested_pilot(df: pd.DataFrame) -> dict:
    from numpy.linalg import lstsq

    if "status" not in df.columns:
        return {"verdict": "NO_STATUS_COLUMN", "n": int(len(df))}
    sub = df[df["status"] == "ok"].copy()
    ycol = "lesion_share_MCA"
    if ycol not in sub.columns or len(sub) < 8:
        return {"verdict": "UNDERPOWERED_OR_INCOMPLETE", "n": int(len(sub))}

    def fit_r2(y, X):
        y = np.asarray(y, float)
        X = np.asarray(X, float)
        m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        y, X = y[m], X[m]
        if len(y) < 8:
            return np.nan, int(len(y))
        Xd = np.column_stack([np.ones(len(y)), X])
        coef, *_ = lstsq(Xd, y, rcond=None)
        pred = Xd @ coef
        sse = float(np.sum((y - pred) ** 2))
        sst = float(np.sum((y - y.mean()) ** 2)) + 1e-12
        return float(1 - sse / sst), int(len(y))

    y = sub[ycol].values
    age = pd.to_numeric(sub["age"], errors="coerce").values
    lvo = pd.to_numeric(sub["lvo_present"], errors="coerce").values
    cow = np.log1p(pd.to_numeric(sub["cow_volume_mm3"], errors="coerce").fillna(0).values)
    r0, n0 = fit_r2(y, age.reshape(-1, 1))
    r1, n1 = fit_r2(y, np.column_stack([age, lvo]))
    r2v, n2 = fit_r2(y, np.column_stack([age, lvo, cow]))
    return {
        "outcome": ycol,
        "n": int(min(n0, n1, n2)),
        "R2_age": r0,
        "R2_age_LVO": r1,
        "R2_age_LVO_CoW": r2v,
        "delta_R2_LVO": float(r1 - r0) if np.isfinite(r1) and np.isfinite(r0) else None,
        "delta_R2_CoW": float(r2v - r1) if np.isfinite(r2v) and np.isfinite(r1) else None,
        "verdict": "EXPLORATORY_PILOT__NOT_CONFIRMATORY",
        "note": "Subject-level angio (LVO/CoW) + ArterialAtlas_level2 lesion shares; enlarge n + QC for confirmatory.",
    }


def main() -> None:
    log(f"Selective ISLES'24 n={N_CASES}")
    case_dirs = download_cases(N_CASES)
    df = analyze_cases(case_dirs)
    df.to_csv(OUT / "isles24_trackC_subject_table.csv", index=False)
    nest = nested_pilot(df)
    status = {
        "date": time.strftime("%Y-%m-%d"),
        "n_requested": N_CASES,
        "n_cases_on_disk": len(case_dirs),
        "n_reg_ok": int((df["status"] == "ok").sum()) if "status" in df.columns else 0,
        "atlas_lock": "ArterialAtlas_level2 major (not TV_BMM)",
        "nested_pilot": nest,
        "verdict": "ISLES24_PILOT_ACQUIRED",
        "trci_bci_excluded": True,
        "table": str(OUT / "isles24_trackC_subject_table.csv"),
    }
    (OUT / "isles24_trackC_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    (ROOT / "results/nature_gates/isles24_trackC_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    (OUT / "isles24_trackC_report.md").write_text(
        "\n".join(
            [
                "# ISLES'24 Track C selective pilot",
                "",
                f"**Verdict:** `{status['verdict']}`",
                f"- Cases: {status['n_cases_on_disk']}; reg OK: {status['n_reg_ok']}",
                f"- Nested: `{json.dumps(nest)}`",
                "",
                "## Locks",
                "- ArterialAtlas_level2 (not TV_BMM); TRCI/BCI excluded",
                "- Exploratory until larger n + registration QC",
                "",
            ]
        ),
        encoding="utf-8",
    )
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
