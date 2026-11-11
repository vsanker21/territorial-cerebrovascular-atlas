#!/usr/bin/env python3
"""Recompute CTA features after SimpleITK affine CTA → NCCT (CoW/LVO space).

Fixes empty sampling caused by raw CTA vs space-ncct CoW grid mismatch.
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
from numpy.linalg import lstsq

ROOT = Path(r"F:\Cerebrovascular SCS")
DATA = ROOT / "data/raw/isles24"
CASES = DATA / "cases"
OUT = ROOT / "results/nature_gates/isles24_trackC"
TABLE = OUT / "isles24_trackC_subject_table.csv"
OUT.mkdir(parents=True, exist_ok=True)

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
}


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def write_nii_bytes(blob, dest: Path) -> Path | None:
    if not isinstance(blob, dict) or blob.get("bytes") is None:
        return None
    raw = blob["bytes"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if raw[:2] == b"\x1f\x8b":
        out = dest if str(dest).endswith(".nii") else Path(str(dest) + ".nii")
        out.write_bytes(gzip.decompress(raw))
        return out
    dest.write_bytes(raw)
    return dest


def ensure_ncct(cdir: Path, sid: str) -> Path | None:
    dest = cdir / f"{sid}_ncct.nii"
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    from huggingface_hub import hf_hub_download
    import pyarrow.parquet as pq

    # Prefer already-cached parquets (train-00000..00079 cover CTA subset)
    for i in list(range(80)) + list(range(80, 149)):
        fname = f"data/train-{i:05d}-of-00149.parquet"
        cached = DATA / fname
        if cached.exists():
            local = str(cached)
        else:
            if i >= 80:
                continue  # don't pull new archives just for NCCT beyond cached set
            local = hf_hub_download(
                "hugging-science/isles24-stroke",
                fname,
                repo_type="dataset",
                local_dir=str(DATA),
            )
        row = pq.read_table(local).to_pydict()
        if row["subject_id"][0] == sid:
            return write_nii_bytes(row["ncct"][0], dest)
    return None


def nib_to_sitk(path: Path) -> sitk.Image:
    return sitk.DICOMOrient(sitk.Cast(sitk.ReadImage(str(path)), sitk.sitkFloat32), "RAS")


def affine_register(moving: sitk.Image, fixed: sitk.Image) -> sitk.Transform:
    fixed_n = sitk.Normalize(fixed)
    moving_n = sitk.Normalize(moving)
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.15)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=80,
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


def resample(moving, fixed, tx, interp):
    return sitk.Resample(moving, fixed, tx, interp, 0.0, moving.GetPixelID())


def cta_features_registered(cdir: Path) -> dict:
    meta = json.loads((cdir / "phenotype.json").read_text(encoding="utf-8"))
    sid = meta["subject_id"]
    files = meta["files"]
    out = {"subject_id": sid, "cta_status": "ok"}
    cta_p = Path(files["cta"]) if files.get("cta") else cdir / f"{sid}_cta.nii"
    cow_p = Path(files["cow_segmentation"]) if files.get("cow_segmentation") else None
    lvo_p = Path(files["lvo_mask"]) if files.get("lvo_mask") else None
    if not cta_p.exists() or cow_p is None or not cow_p.exists():
        out["cta_status"] = "missing"
        return out

    ncct_p = ensure_ncct(cdir, sid)
    # Fixed image: prefer NCCT; else DWI (same grid as CoW)
    if ncct_p and ncct_p.exists():
        fixed_path = ncct_p
        out["cta_reg_fixed"] = "ncct"
    else:
        dwi = Path(files.get("dwi") or "")
        if not dwi.exists():
            out["cta_status"] = "no_fixed"
            return out
        fixed_path = dwi
        out["cta_reg_fixed"] = "dwi"

    mni_dir = cdir / "cta_ncct"
    mni_dir.mkdir(parents=True, exist_ok=True)
    out_cta = mni_dir / "cta_on_ncct.nii.gz"
    try:
        if out_cta.exists():
            cta_on = nib_to_sitk(out_cta)
        else:
            log(f"  register CTA→{out['cta_reg_fixed']} {sid}")
            fixed = nib_to_sitk(fixed_path)
            moving = nib_to_sitk(cta_p)
            tx = affine_register(moving, fixed)
            cta_on = resample(moving, fixed, tx, sitk.sitkLinear)
            sitk.WriteImage(cta_on, str(out_cta))
        # Load arrays on CoW grid
        cta_nib = nib.load(str(out_cta))
        cow_img = nib.load(str(cow_p))
        cta = np.asanyarray(cta_nib.dataobj, dtype=np.float32)
        cow = np.asanyarray(cow_img.dataobj)
        if cta.shape != cow.shape:
            # resample cow to cta_on using sitk
            cow_s = nib_to_sitk(cow_p)
            cow_r = sitk.Resample(
                cow_s,
                cta_on,
                sitk.Transform(),
                sitk.sitkNearestNeighbor,
                0.0,
                cow_s.GetPixelID(),
            )
            tmp = mni_dir / "cow_on_cta.nii.gz"
            sitk.WriteImage(cow_r, str(tmp))
            cow = np.asanyarray(nib.load(str(tmp)).dataobj)

        fg = cow > 0
        if fg.sum() < 10:
            out["cta_status"] = "empty_cow_overlap"
            return out
        vals = cta[fg]
        out["cta_cow_mean"] = float(vals.mean())
        out["cta_cow_median"] = float(np.median(vals))
        out["cta_cow_p10"] = float(np.percentile(vals, 10))
        out["cta_cow_p90"] = float(np.percentile(vals, 90))
        out["cta_cow_range_p90_p10"] = float(out["cta_cow_p90"] - out["cta_cow_p10"])

        for lab, name in TOPCOW.items():
            m = cow == lab
            out[f"cta_mean_{name}"] = float(cta[m].mean()) if m.sum() >= 5 else np.nan
        for terr, labs in {
            "MCA": [5, 7],
            "ACA": [11, 12],
            "PCA": [2, 3],
            "ICA": [4, 6],
            "BA": [1],
        }.items():
            m = np.isin(cow, labs)
            out[f"cta_mean_{terr}"] = float(cta[m].mean()) if m.sum() >= 5 else np.nan

        for a, b, key in [
            ("cta_mean_R_MCA", "cta_mean_L_MCA", "cta_asym_MCA"),
            ("cta_mean_R_ICA", "cta_mean_L_ICA", "cta_asym_ICA"),
        ]:
            va, vb = out.get(a, np.nan), out.get(b, np.nan)
            out[key] = float((va - vb) / (va + vb + 1e-6)) if np.isfinite(va) and np.isfinite(vb) else np.nan

        if lvo_p and lvo_p.exists():
            lvo = np.asanyarray(nib.load(str(lvo_p)).dataobj) > 0
            if lvo.shape != cta.shape:
                lvo_s = nib_to_sitk(lvo_p)
                lvo_r = sitk.Resample(
                    lvo_s, cta_on, sitk.Transform(), sitk.sitkNearestNeighbor, 0.0, lvo_s.GetPixelID()
                )
                tmp = mni_dir / "lvo_on_cta.nii.gz"
                sitk.WriteImage(lvo_r, str(tmp))
                lvo = np.asanyarray(nib.load(str(tmp)).dataobj) > 0
            if lvo.any():
                out["cta_mean_LVO"] = float(cta[lvo].mean())
                out["cta_LVO_over_CoW"] = float(out["cta_mean_LVO"] / (out["cta_cow_median"] + 1e-6))
            else:
                out["cta_mean_LVO"] = np.nan
                out["cta_LVO_over_CoW"] = np.nan
        return out
    except Exception as e:
        out["cta_status"] = f"reg_failed:{type(e).__name__}:{e}"
        return out


def ols_r2(y, X):
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    if len(y) < max(12, X.shape[1] + 4):
        return np.nan, int(len(y))
    Xd = np.column_stack([np.ones(len(y)), X])
    coef, *_ = lstsq(Xd, y, rcond=None)
    pred = Xd @ coef
    sse = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2)) + 1e-12
    return float(1 - sse / sst), int(len(y))


def cv_r2(y, X, k=5, seed=0):
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    n = len(y)
    if n < max(15, X.shape[1] + 5):
        return np.nan
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    folds = np.array_split(idx, k)
    preds = np.full(n, np.nan)
    for i in range(k):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(k) if j != i])
        if len(tr) <= X.shape[1] + 1:
            continue
        coef, *_ = lstsq(np.column_stack([np.ones(len(tr)), X[tr]]), y[tr], rcond=None)
        preds[te] = np.column_stack([np.ones(len(te)), X[te]]) @ coef
    m2 = np.isfinite(preds)
    if m2.sum() < 10:
        return np.nan
    sse = float(np.sum((y[m2] - preds[m2]) ** 2))
    sst = float(np.sum((y[m2] - y[m2].mean()) ** 2)) + 1e-12
    return float(1 - sse / sst)


def bootstrap_delta(y, X0, X1, B=2000, seed=5):
    m = np.isfinite(y) & np.all(np.isfinite(X0), axis=1) & np.all(np.isfinite(X1), axis=1)
    y, X0, X1 = y[m], X0[m], X1[m]
    n = len(y)
    r0, _ = ols_r2(y, X0)
    r1, _ = ols_r2(y, X1)
    delta = float(r1 - r0) if np.isfinite(r0) and np.isfinite(r1) else np.nan
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(B):
        take = rng.choice(n, size=n, replace=True)
        a, _ = ols_r2(y[take], X0[take])
        b, _ = ols_r2(y[take], X1[take])
        if np.isfinite(a) and np.isfinite(b):
            boots.append(b - a)
    boots = np.asarray(boots)
    return {
        "n": int(n),
        "R2_base": float(r0) if np.isfinite(r0) else None,
        "R2_full": float(r1) if np.isfinite(r1) else None,
        "delta": delta,
        "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
        if len(boots)
        else [None, None],
        "ci_excludes_0": bool(len(boots) and np.percentile(boots, 2.5) > 0),
        "p_boot_delta_le_0": float(np.mean(boots <= 0)) if len(boots) else None,
    }


def main():
    # CTA cases = those with cta.nii
    cdirs = sorted([p.parent for p in CASES.glob("*/sub-stroke*_cta.nii")])
    log(f"CTA cases {len(cdirs)}")
    rows = []
    for cdir in cdirs:
        log(f"CTA-reg features {cdir.name}")
        rows.append(cta_features_registered(cdir))
    feat = pd.DataFrame(rows)
    feat.to_csv(OUT / "isles24_cta_features.csv", index=False)
    log(f"finite cta_mean_MCA {feat['cta_mean_MCA'].notna().sum()}/{len(feat)}")

    base = pd.read_csv(TABLE)
    merged = base.merge(feat, on="subject_id", how="inner")
    if "qc_pass" in merged.columns:
        merged = merged[merged["qc_pass"] == 1].copy()
    merged = merged[merged["cta_status"] == "ok"].copy()
    merged.to_csv(OUT / "isles24_trackC_with_cta_features.csv", index=False)

    y = pd.to_numeric(merged["lesion_share_MCA"], errors="coerce").values
    age = pd.to_numeric(merged["age"], errors="coerce").values
    lvo = np.log1p(pd.to_numeric(merged["lvo_volume_mm3"], errors="coerce").fillna(0).values)
    cow_m = np.log1p(pd.to_numeric(merged["cow_vol_MCA"], errors="coerce").fillna(0).values)
    cta_mca = pd.to_numeric(merged["cta_mean_MCA"], errors="coerce").values
    cta_asym = pd.to_numeric(merged["cta_asym_MCA"], errors="coerce").values
    cta_lvo = pd.to_numeric(merged["cta_LVO_over_CoW"], errors="coerce").values
    cta_range = pd.to_numeric(merged["cta_cow_range_p90_p10"], errors="coerce").values

    specs = {
        "age": np.column_stack([age]),
        "age_CoW_LVO": np.column_stack([age, lvo, cow_m]),
        "age_CTA": np.column_stack([age, cta_mca, cta_asym, cta_lvo, cta_range]),
        "age_CoW_LVO_CTA": np.column_stack([age, lvo, cow_m, cta_mca, cta_asym, cta_lvo, cta_range]),
    }
    models = []
    for name, X in specs.items():
        r2, n = ols_r2(y, X)
        models.append({"model": name, "R2_in": r2, "R2_cv5": cv_r2(y, X, seed=11), "n": n, "p": X.shape[1]})
    deltas = {
        "CTA_minus_age": bootstrap_delta(y, specs["age"], specs["age_CTA"]),
        "CTA_minus_CoWLVO": bootstrap_delta(y, specs["age_CoW_LVO"], specs["age_CoW_LVO_CTA"]),
        "full_minus_age": bootstrap_delta(y, specs["age"], specs["age_CoW_LVO_CTA"]),
    }
    from scipy import stats

    spear = []
    for col in ["cta_mean_MCA", "cta_asym_MCA", "cta_LVO_over_CoW", "cta_cow_range_p90_p10", "cta_mean_LVO"]:
        xv = pd.to_numeric(merged[col], errors="coerce")
        m = np.isfinite(y) & xv.notna()
        if m.sum() < 10:
            continue
        rho, p = stats.spearmanr(y[m], xv[m])
        spear.append({"feature": col, "rho": float(rho), "p": float(p), "n": int(m.sum())})

    full = next(r for r in models if r["model"] == "age_CoW_LVO_CTA")
    d = deltas["full_minus_age"]
    if np.isfinite(full["R2_cv5"] or np.nan) and full["R2_cv5"] > 0 and d.get("ci_excludes_0"):
        verdict = "CTA_CV_POSITIVE_EXPLORATORY__REPLICATE"
    elif d.get("ci_excludes_0") and (full.get("R2_cv5") or -1) <= 0:
        verdict = "CTA_INSAMPLE_DELTA_BUT_CV_NONPOSITIVE__NOT_CONFIRMATORY"
    else:
        verdict = "CTA_FEATURES_NULL_OR_WEAK__NOT_CONFIRMATORY"

    status = {
        "date": time.strftime("%Y-%m-%d"),
        "n": int(len(merged)),
        "models": models,
        "delta_R2": deltas,
        "spearman": spear,
        "verdict": verdict,
        "method": "SimpleITK_affine_CTA_to_NCCT_then_TopCoW_intensity",
        "trci_bci_excluded": True,
    }
    (OUT / "isles24_cta_features_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    lines = [
        "# ISLES'24 CTA features (registered to NCCT)",
        "",
        f"**Verdict:** `{verdict}`",
        f"- n={status['n']}",
        f"- method: {status['method']}",
        "",
        "## Models",
    ]
    for r in models:
        lines.append(f"- **{r['model']}**: R²_in={r['R2_in']}; R²_cv5={r['R2_cv5']}; n={r['n']}")
    lines += ["", "## ΔR²", f"```\n{json.dumps(deltas, indent=2)}\n```", "", "## Spearman", f"```\n{json.dumps(spear, indent=2)}\n```", ""]
    (OUT / "isles24_cta_features_report.md").write_text("\n".join(lines), encoding="utf-8")
    # Update Yang email with confirmed address
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
