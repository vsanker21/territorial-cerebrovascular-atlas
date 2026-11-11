#!/usr/bin/env python3
"""Expand ISLES'24 Track C: more cases, TopCoW multi-label features, registration QC.

TopCoW labels (Yang et al.):
  1 BA, 2 R-PCA, 3 L-PCA, 4 R-ICA, 5 R-MCA, 6 L-ICA, 7 L-MCA,
  8 R-Pcom, 9 L-Pcom, 10 Acom, 11 R-ACA, 12 L-ACA, 15 3rd-A2

Locks: ArterialAtlas_level2 (not TV_BMM); exploratory until larger n + QC; no TRCI/BCI.
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
from numpy.linalg import lstsq

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

# Expand to 80 total (indices 0..79)
N_TOTAL = 80

TOPCOW = {
    1: "BA",
    2: "R_PCA",
    3: "L_PCA",
    4: "R_ICA",
    5: "R_MCA",
    6: "L_ICA",
    7: "L_MCA",
    8: "R_Pcom",
    9: "L_Pcom",
    10: "Acom",
    11: "R_ACA",
    12: "L_ACA",
    15: "A2_3rd",
}
TERR_AGG = {
    "ACA": [11, 12],
    "MCA": [5, 7],
    "PCA": [2, 3],
    "VB": [1],
    "ICA": [4, 6],
    "Pcom": [8, 9],
    "Acom": [10],
}


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


def download_cases(n_total: int = N_TOTAL) -> list[Path]:
    from huggingface_hub import hf_hub_download
    import pyarrow.parquet as pq

    CASES.mkdir(parents=True, exist_ok=True)
    extracted = []
    for i in range(n_total):
        fname = f"data/train-{i:05d}-of-00149.parquet"
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
            log(f"exists {sid}")
            extracted.append(cdir)
            continue
        log(f"extract {fname} -> {sid}")
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
        # Skip full CTA volume (hundreds of MB); CoW/LVO masks encode subject angio.
        for field in ("dwi", "lesion_mask", "lvo_mask", "cow_segmentation"):
            dest = cdir / f"{sid}_{field}.nii"
            p = write_nii_bytes(row[field][0], dest)
            meta["files"][field] = str(p) if p else None
            if p:
                log(f"  {field}: {p.stat().st_size/1e6:.1f} MB")
        meta["files"]["cta"] = None
        (cdir / "phenotype.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
        extracted.append(cdir)
    return extracted


def nib_to_sitk(path: Path) -> sitk.Image:
    return sitk.DICOMOrient(sitk.Cast(sitk.ReadImage(str(path)), sitk.sitkFloat32), "RAS")


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
        numberOfIterations=120,
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


def apply_tx(moving, reference, tx, interp):
    return sitk.Resample(moving, reference, tx, interp, 0.0, moving.GetPixelID())


def mni_template_sitk():
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


def ensure_major_atlas(mni_nib):
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


def cow_features(path: Path | None) -> dict:
    out = {f"cow_present_{v}": 0 for v in TOPCOW.values()}
    out.update({f"cow_vol_{k}": 0.0 for k in TERR_AGG})
    out["cow_n_labels"] = 0
    out["cow_completeness"] = np.nan
    out["cow_volume_mm3"] = np.nan
    out["cow_laterality_MCA"] = np.nan  # (R-L)/(R+L)
    if path is None or not Path(path).exists():
        return out
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    vox = float(np.prod(img.header.get_zooms()[:3]))
    out["cow_volume_mm3"] = float((data > 0).sum() * vox)
    present = []
    for lab, name in TOPCOW.items():
        n = int((data == lab).sum())
        out[f"cow_present_{name}"] = int(n > 0)
        if n > 0:
            present.append(name)
        out[f"cow_vol_{name}"] = float(n * vox)
    for terr, labs in TERR_AGG.items():
        out[f"cow_vol_{terr}"] = float(sum(int((data == lab).sum()) for lab in labs) * vox)
        out[f"cow_present_{terr}"] = int(any((data == lab).sum() > 0 for lab in labs))
    # Completeness: fraction of core CoW labels present (BA, bilat ICA/MCA/ACA/PCA, Acom, bilat Pcom)
    core = ["BA", "R_ICA", "L_ICA", "R_MCA", "L_MCA", "R_ACA", "L_ACA", "R_PCA", "L_PCA", "Acom", "R_Pcom", "L_Pcom"]
    out["cow_n_labels"] = int(sum(out[f"cow_present_{n}"] for n in TOPCOW.values()))
    out["cow_completeness"] = float(sum(out[f"cow_present_{n}"] for n in core) / len(core))
    r = out.get("cow_vol_R_MCA", 0.0)
    l = out.get("cow_vol_L_MCA", 0.0)
    out["cow_laterality_MCA"] = float((r - l) / (r + l + 1e-6))
    return out


def mask_volume_mm3(path: Path | None) -> float | None:
    if path is None or not Path(path).exists():
        return None
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    return float((data > 0).sum() * float(np.prod(img.header.get_zooms()[:3])))


def lvo_stats(path: Path | None) -> dict:
    if path is None or not Path(path).exists():
        return {"lvo_present": None, "lvo_volume_mm3": None, "lvo_centroid_x": None}
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj) > 0
    vox = float(np.prod(img.header.get_zooms()[:3]))
    if not data.any():
        return {"lvo_present": 0, "lvo_volume_mm3": 0.0, "lvo_centroid_x": np.nan}
    coords = np.column_stack(np.where(data))
    # axis 0 ~ i; use world x from affine
    ijk = coords.mean(axis=0)
    xyz = img.affine.dot([ijk[0], ijk[1], ijk[2], 1.0])
    return {
        "lvo_present": 1,
        "lvo_volume_mm3": float(data.sum() * vox),
        "lvo_centroid_x": float(xyz[0]),
    }


def register_lesion(dwi_path, lesion_path, out_dir, mni_sitk, major, major_aff) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        dwi = nib_to_sitk(dwi_path)
        lesion = nib_to_sitk(lesion_path)
        tx = affine_to_mni(dwi, mni_sitk)
        dwi_mni = apply_tx(dwi, mni_sitk, tx, sitk.sitkLinear)
        les_mni = apply_tx(lesion, mni_sitk, tx, sitk.sitkNearestNeighbor)
        out_les = out_dir / "lesion_mni.nii.gz"
        out_dwi = out_dir / "dwi_mni.nii.gz"
        sitk.WriteImage(les_mni, str(out_les))
        sitk.WriteImage(dwi_mni, str(out_dwi))

        # Correlation QC: brain overlap proxy — fraction of lesion inside major atlas FOV
        les_nib = nib.load(str(out_les))
        ref = nib.Nifti1Image(major.astype(np.float32), major_aff)
        les_on = image.resample_to_img(
            les_nib, ref, interpolation="nearest", force_resample=True, copy_header=True
        )
        les = np.asanyarray(les_on.dataobj) > 0
        atlas_fg = major > 0
        total = int(les.sum())
        in_atlas = int((les & atlas_fg).sum())
        names = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
        shares = {
            "status": "ok",
            "method": "SimpleITK_affine_DWI_to_MNI_nearest_lesion",
            "atlas": "ArterialAtlas_level2_major",
            "lesion_vox_total": total,
            "lesion_frac_in_atlas": float(in_atlas / total) if total else np.nan,
            "qc_lesion_empty": int(total == 0),
            "qc_low_atlas_overlap": int((in_atlas / total) < 0.5) if total else 1,
            "lesion_mni": str(out_les),
        }
        for code, name in names.items():
            vox = int(((major == code) & les).sum())
            shares[f"lesion_vox_{name}"] = vox
            shares[f"lesion_share_{name}"] = float(vox / total) if total else np.nan
        # QC pass: non-empty lesion and ≥50% overlap with atlas
        shares["qc_pass"] = int(
            (not shares["qc_lesion_empty"]) and (not shares["qc_low_atlas_overlap"])
        )
        return shares
    except Exception as e:
        return {"status": "reg_failed", "error": str(e), "qc_pass": 0}


def analyze_cases(case_dirs: list[Path]) -> pd.DataFrame:
    log("MNI + atlas")
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
            "lesion_volume_native_mm3": mask_volume_mm3(lesion),
        }
        row.update(lvo_stats(lvo))
        row.update(cow_features(cow))
        if dwi and lesion and dwi.exists() and lesion.exists():
            # Skip re-reg if qc already present unless forced
            mni_dir = cdir / "mni"
            prev = mni_dir / "lesion_mni.nii.gz"
            if prev.exists() and (cdir / "reg_features.json").exists():
                reg = json.loads((cdir / "reg_features.json").read_text(encoding="utf-8"))
                log(f"reuse reg {meta['subject_id']}")
            else:
                log(f"register {meta['subject_id']}")
                reg = register_lesion(dwi, lesion, mni_dir, mni_sitk, major, major_aff)
                (cdir / "reg_features.json").write_text(json.dumps(reg, indent=2), encoding="utf-8")
            row.update(reg)
        else:
            row["status"] = "missing_dwi_or_lesion"
            row["qc_pass"] = 0
        rows.append(row)
    return pd.DataFrame(rows)


def fit_r2(y, X):
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    if len(y) < 10 or X.shape[1] == 0:
        return np.nan, int(len(y))
    Xd = np.column_stack([np.ones(len(y)), X])
    coef, *_ = lstsq(Xd, y, rcond=None)
    pred = Xd @ coef
    sse = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2)) + 1e-12
    return float(1 - sse / sst), int(len(y))


def nested_models(df: pd.DataFrame) -> dict:
    sub = df[(df.get("status") == "ok") & (df.get("qc_pass") == 1)].copy()
    if len(sub) < 10:
        sub = df[df.get("status") == "ok"].copy()
        qc_note = "QC filter too strict — using all reg-ok"
    else:
        qc_note = "QC-pass only"
    ycol = "lesion_share_MCA"
    if ycol not in sub.columns or len(sub) < 10:
        return {"verdict": "UNDERPOWERED", "n": int(len(sub))}

    age = pd.to_numeric(sub["age"], errors="coerce").values
    lvo_vol = np.log1p(pd.to_numeric(sub["lvo_volume_mm3"], errors="coerce").fillna(0).values)
    cow_c = pd.to_numeric(sub["cow_completeness"], errors="coerce").values
    cow_mca = np.log1p(pd.to_numeric(sub["cow_vol_MCA"], errors="coerce").fillna(0).values)
    cow_lat = pd.to_numeric(sub["cow_laterality_MCA"], errors="coerce").values
    y = sub[ycol].values

    models = {}
    specs = {
        "age": age.reshape(-1, 1),
        "age_LVOvol": np.column_stack([age, lvo_vol]),
        "age_CoWcomplete": np.column_stack([age, cow_c]),
        "age_CoWMCA": np.column_stack([age, cow_mca]),
        "age_LVOvol_CoWcomplete_CoWMCA": np.column_stack([age, lvo_vol, cow_c, cow_mca]),
        "age_full_angio": np.column_stack([age, lvo_vol, cow_c, cow_mca, cow_lat]),
    }
    for name, X in specs.items():
        r2, n = fit_r2(y, X)
        models[name] = {"R2": r2, "n": n}

    # deltas vs age
    base = models["age"]["R2"]
    deltas = {
        k: (float(v["R2"] - base) if np.isfinite(v["R2"]) and np.isfinite(base) else None)
        for k, v in models.items()
        if k != "age"
    }
    return {
        "outcome": ycol,
        "qc_note": qc_note,
        "n": int(models["age"]["n"]),
        "models": models,
        "delta_R2_vs_age": deltas,
        "mean_lesion_shares": {
            t: float(pd.to_numeric(sub[f"lesion_share_{t}"], errors="coerce").mean())
            for t in ("ACA", "MCA", "PCA", "VB")
            if f"lesion_share_{t}" in sub.columns
        },
        "verdict": "EXPLORATORY_EXPANDED__NOT_CONFIRMATORY",
        "note": "Subject-level CoW/LVO features; confirmatory needs larger n + registration QC audit.",
    }


def main() -> None:
    log(f"Expand ISLES'24 to n={N_TOTAL}")
    case_dirs = download_cases(N_TOTAL)
    # Force recompute CoW features even if reg reused — always re-analyze
    # Clear stale reg_features only for new cases; for old cases without cow in table, re-run full analyze
    # Delete reg_features to refresh QC fields on all? Better: always register if reg_features missing qc_pass
    for cdir in case_dirs:
        rf = cdir / "reg_features.json"
        if rf.exists():
            reg = json.loads(rf.read_text(encoding="utf-8"))
            if "qc_pass" not in reg:
                rf.unlink()  # force re-reg with QC
    df = analyze_cases(case_dirs)
    df.to_csv(OUT / "isles24_trackC_subject_table.csv", index=False)
    nest = nested_models(df)
    n_ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    n_qc = int((df.get("qc_pass") == 1).sum()) if "qc_pass" in df.columns else 0
    status = {
        "date": time.strftime("%Y-%m-%d"),
        "n_requested": N_TOTAL,
        "n_cases_on_disk": len(case_dirs),
        "n_reg_ok": n_ok,
        "n_qc_pass": n_qc,
        "atlas_lock": "ArterialAtlas_level2 major (not TV_BMM)",
        "cow_label_scheme": "TopCoW",
        "nested": nest,
        "verdict": "ISLES24_EXPANDED_WITH_COW_FEATURES",
        "trci_bci_excluded": True,
        "table": str(OUT / "isles24_trackC_subject_table.csv"),
    }
    (OUT / "isles24_trackC_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    (ROOT / "results/nature_gates/isles24_trackC_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    lines = [
        "# ISLES'24 Track C expanded pilot",
        "",
        f"**Verdict:** `{status['verdict']}`",
        f"- Cases: {status['n_cases_on_disk']}; reg OK: {n_ok}; QC pass: {n_qc}",
        f"- Nested: `{json.dumps(nest, indent=2)}`",
        "",
        "## Features",
        "- TopCoW multi-label presence/volumes; MCA laterality; completeness",
        "- LVO volume + centroid x (binary LVO often invariant)",
        "- Lesion shares on ArterialAtlas_level2 majors",
        "",
        "## Locks",
        "- Not TV_BMM; exploratory; TRCI/BCI excluded",
        "",
    ]
    (OUT / "isles24_trackC_report.md").write_text("\n".join(lines), encoding="utf-8")
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
