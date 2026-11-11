#!/usr/bin/env python3
"""Stronger subject-level CTA features for ISLES'24 Track C.

Uses CTA intensity sampled inside TopCoW labels and LVO mask as patency /
hypoattenuation proxies (no full vesselness CNN).

Steps:
  1. Ensure CTA on disk for TARGET_CTA cases (download from HF parquet if needed)
  2. Extract intensity / asymmetry features
  3. Merge with existing subject table; nested + CV inference on CTA subset
"""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from numpy.linalg import lstsq

ROOT = Path(r"F:\Cerebrovascular SCS")
DATA = ROOT / "data/raw/isles24"
CASES = DATA / "cases"
OUT = ROOT / "results/nature_gates/isles24_trackC"
OUT.mkdir(parents=True, exist_ok=True)
TABLE = OUT / "isles24_trackC_subject_table.csv"

TARGET_CTA = 40  # QC subset with full CTA
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


def ensure_cta_for_subset(n: int = TARGET_CTA) -> list[Path]:
    """Download CTA only for first n parquet indices that have cases on disk / extract."""
    from huggingface_hub import hf_hub_download
    import pyarrow.parquet as pq

    have = []
    for i in range(n):
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
        cdir.mkdir(parents=True, exist_ok=True)
        cta_path = cdir / f"{sid}_cta.nii"
        if cta_path.exists() and cta_path.stat().st_size > 1_000_000:
            log(f"CTA exists {sid}")
        else:
            log(f"CTA download/extract {sid}")
            p = write_nii_bytes(row["cta"][0], cta_path)
            if not p:
                continue
            log(f"  wrote {p.stat().st_size/1e6:.1f} MB")
        # Ensure CoW/LVO exist for feature sampling
        meta_path = cdir / "phenotype.json"
        if meta_path.exists():
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        else:
            meta = {"subject_id": sid, "files": {}}
        meta.setdefault("files", {})
        meta["files"]["cta"] = str(cta_path)
        for field in ("lvo_mask", "cow_segmentation", "dwi", "lesion_mask"):
            dest = cdir / f"{sid}_{field}.nii"
            if not dest.exists():
                p = write_nii_bytes(row[field][0], dest)
                meta["files"][field] = str(p) if p else None
            else:
                meta["files"][field] = str(dest)
        for k in ("age", "sex", "nihss_admission", "mrs_admission", "mrs_3month"):
            if k in row and (k not in meta or meta.get(k) is None):
                meta[k] = row[k][0]
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        have.append(cdir)
    return have


def robust_stats(x: np.ndarray) -> dict:
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return {"n": int(len(x)), "mean": np.nan, "median": np.nan, "p10": np.nan, "p90": np.nan}
    return {
        "n": int(len(x)),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "p10": float(np.percentile(x, 10)),
        "p90": float(np.percentile(x, 90)),
    }


def cta_features(cdir: Path) -> dict:
    meta = json.loads((cdir / "phenotype.json").read_text(encoding="utf-8"))
    sid = meta["subject_id"]
    files = meta["files"]
    cta_p = Path(files["cta"]) if files.get("cta") else None
    cow_p = Path(files["cow_segmentation"]) if files.get("cow_segmentation") else None
    lvo_p = Path(files["lvo_mask"]) if files.get("lvo_mask") else None
    out = {"subject_id": sid, "cta_status": "ok"}
    if cta_p is None or not cta_p.exists() or cow_p is None or not cow_p.exists():
        out["cta_status"] = "missing"
        return out

    cta_img = nib.load(str(cta_p))
    cow_img = nib.load(str(cow_p))
    cta = np.asanyarray(cta_img.dataobj, dtype=np.float32)
    cow = np.asanyarray(cow_img.dataobj)

    # Resample CoW to CTA grid if needed
    if cta.shape != cow.shape:
        from nilearn.image import resample_to_img

        cow_r = resample_to_img(
            cow_img, cta_img, interpolation="nearest", force_resample=True, copy_header=True
        )
        cow = np.asanyarray(cow_r.dataobj)

    # Global CoW CTA intensity
    fg = cow > 0
    st = robust_stats(cta[fg])
    out["cta_cow_mean"] = st["mean"]
    out["cta_cow_median"] = st["median"]
    out["cta_cow_p10"] = st["p10"]
    out["cta_cow_p90"] = st["p90"]

    # Per-label means (patency proxy)
    for lab, name in TOPCOW.items():
        m = cow == lab
        if m.sum() >= 5:
            out[f"cta_mean_{name}"] = float(cta[m].mean())
        else:
            out[f"cta_mean_{name}"] = np.nan

    # Territory aggregates
    for terr, labs in {
        "MCA": [5, 7],
        "ACA": [11, 12],
        "PCA": [2, 3],
        "ICA": [4, 6],
        "BA": [1],
    }.items():
        m = np.isin(cow, labs)
        out[f"cta_mean_{terr}"] = float(cta[m].mean()) if m.sum() >= 5 else np.nan

    # Asymmetry (R-L)/(R+L) on MCA and ICA means
    for side_pair, key in [
        (("cta_mean_R_MCA", "cta_mean_L_MCA"), "cta_asym_MCA"),
        (("cta_mean_R_ICA", "cta_mean_L_ICA"), "cta_asym_ICA"),
    ]:
        a, b = out.get(side_pair[0], np.nan), out.get(side_pair[1], np.nan)
        if np.isfinite(a) and np.isfinite(b):
            out[key] = float((a - b) / (a + b + 1e-6))
        else:
            out[key] = np.nan

    # LVO region CTA intensity vs contralateral MCA mean
    if lvo_p is not None and lvo_p.exists():
        lvo_img = nib.load(str(lvo_p))
        lvo = np.asanyarray(lvo_img.dataobj) > 0
        if lvo.shape != cta.shape:
            from nilearn.image import resample_to_img

            lvo_r = resample_to_img(
                lvo_img, cta_img, interpolation="nearest", force_resample=True, copy_header=True
            )
            lvo = np.asanyarray(lvo_r.dataobj) > 0
        if lvo.any():
            out["cta_mean_LVO"] = float(cta[lvo].mean())
            # hypoattenuation index: LVO / CoW median
            out["cta_LVO_over_CoW"] = float(out["cta_mean_LVO"] / (out["cta_cow_median"] + 1e-6))
        else:
            out["cta_mean_LVO"] = np.nan
            out["cta_LVO_over_CoW"] = np.nan
    else:
        out["cta_mean_LVO"] = np.nan
        out["cta_LVO_over_CoW"] = np.nan

    # Contrast dynamic range in CoW
    out["cta_cow_range_p90_p10"] = (
        float(out["cta_cow_p90"] - out["cta_cow_p10"])
        if np.isfinite(out["cta_cow_p90"]) and np.isfinite(out["cta_cow_p10"])
        else np.nan
    )
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
        Xtr = np.column_stack([np.ones(len(tr)), X[tr]])
        Xte = np.column_stack([np.ones(len(te)), X[te]])
        coef, *_ = lstsq(Xtr, y[tr], rcond=None)
        preds[te] = Xte @ coef
    m2 = np.isfinite(preds)
    if m2.sum() < 10:
        return np.nan
    sse = float(np.sum((y[m2] - preds[m2]) ** 2))
    sst = float(np.sum((y[m2] - y[m2].mean()) ** 2)) + 1e-12
    return float(1 - sse / sst)


def bootstrap_delta(y, X0, X1, B=2000, seed=3):
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
        "R2_base": float(r0),
        "R2_full": float(r1),
        "delta": delta,
        "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
        if len(boots)
        else [np.nan, np.nan],
        "ci_excludes_0": bool(len(boots) and np.percentile(boots, 2.5) > 0),
        "p_boot_delta_le_0": float(np.mean(boots <= 0)) if len(boots) else np.nan,
    }


def main() -> None:
    log(f"Ensure CTA for n={TARGET_CTA}")
    cdirs = ensure_cta_for_subset(TARGET_CTA)
    rows = []
    for cdir in cdirs:
        log(f"CTA features {cdir.name}")
        rows.append(cta_features(cdir))
    feat = pd.DataFrame(rows)
    feat.to_csv(OUT / "isles24_cta_features.csv", index=False)

    base = pd.read_csv(TABLE) if TABLE.exists() else pd.DataFrame()
    if len(base):
        merged = base.merge(feat, on="subject_id", how="inner")
    else:
        merged = feat
    # Keep QC-pass if present
    if "qc_pass" in merged.columns:
        merged = merged[merged["qc_pass"] == 1].copy()
    merged.to_csv(OUT / "isles24_trackC_with_cta_features.csv", index=False)
    log(f"CTA-feature table n={len(merged)}")

    y = pd.to_numeric(merged["lesion_share_MCA"], errors="coerce").values
    age = pd.to_numeric(merged["age"], errors="coerce").values
    # Prior CoW/LVO
    lvo = np.log1p(pd.to_numeric(merged.get("lvo_volume_mm3"), errors="coerce").fillna(0).values)
    cow_m = np.log1p(pd.to_numeric(merged.get("cow_vol_MCA"), errors="coerce").fillna(0).values)
    # New CTA
    cta_mca = pd.to_numeric(merged.get("cta_mean_MCA"), errors="coerce").values
    cta_asym = pd.to_numeric(merged.get("cta_asym_MCA"), errors="coerce").values
    cta_lvo = pd.to_numeric(merged.get("cta_LVO_over_CoW"), errors="coerce").values
    cta_range = pd.to_numeric(merged.get("cta_cow_range_p90_p10"), errors="coerce").values

    specs = {
        "age": np.column_stack([age]),
        "age_CoW_LVO": np.column_stack([age, lvo, cow_m]),
        "age_CTA": np.column_stack([age, cta_mca, cta_asym, cta_lvo, cta_range]),
        "age_CoW_LVO_CTA": np.column_stack([age, lvo, cow_m, cta_mca, cta_asym, cta_lvo, cta_range]),
    }
    model_rows = []
    for name, X in specs.items():
        r2, n = ols_r2(y, X)
        model_rows.append(
            {
                "model": name,
                "R2_in": r2,
                "R2_cv5": cv_r2(y, X, seed=11),
                "n": n,
                "p": X.shape[1],
            }
        )
    deltas = {
        "CTA_minus_age": bootstrap_delta(y, specs["age"], specs["age_CTA"]),
        "CTA_minus_CoWLVO": bootstrap_delta(y, specs["age_CoW_LVO"], specs["age_CoW_LVO_CTA"]),
        "full_minus_age": bootstrap_delta(y, specs["age"], specs["age_CoW_LVO_CTA"]),
    }

    # Spearman of CTA features vs MCA share
    from scipy import stats

    spear = []
    for col in [
        "cta_mean_MCA",
        "cta_asym_MCA",
        "cta_LVO_over_CoW",
        "cta_cow_range_p90_p10",
        "cta_mean_LVO",
        "cta_mean_ICA",
    ]:
        if col not in merged.columns:
            continue
        xv = pd.to_numeric(merged[col], errors="coerce")
        m = np.isfinite(y) & xv.notna()
        if m.sum() < 10:
            continue
        rho, p = stats.spearmanr(y[m], xv[m])
        spear.append({"feature": col, "rho": float(rho), "p": float(p), "n": int(m.sum())})

    status = {
        "date": time.strftime("%Y-%m-%d"),
        "n_cta_subset": int(len(merged)),
        "models": model_rows,
        "delta_R2": deltas,
        "spearman": spear,
        "verdict": "CTA_FEATURES_EXPLORATORY",
        "note": (
            "CTA intensity-in-CoW / LVO proxies. Confirmatory only if CV R² improves "
            "and holds at larger n; still not molecular Nature package alone."
        ),
        "trci_bci_excluded": True,
    }
    # Upgrade verdict if CV positive and delta CI excludes 0
    full = next(r for r in model_rows if r["model"] == "age_CoW_LVO_CTA")
    d = deltas["full_minus_age"]
    if (
        np.isfinite(full.get("R2_cv5") or np.nan)
        and full["R2_cv5"] > 0
        and d.get("ci_excludes_0")
    ):
        status["verdict"] = "CTA_CV_POSITIVE_EXPLORATORY__REPLICATE"
    elif d.get("ci_excludes_0") and (full.get("R2_cv5") or -1) <= 0:
        status["verdict"] = "CTA_INSAMPLE_DELTA_BUT_CV_NONPOSITIVE__NOT_CONFIRMATORY"

    (OUT / "isles24_cta_features_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    lines = [
        "# ISLES'24 stronger CTA features (Track C)",
        "",
        f"**Verdict:** `{status['verdict']}`",
        f"- CTA subset n={status['n_cta_subset']}",
        "",
        "## Models (MCA lesion share)",
    ]
    for r in model_rows:
        lines.append(f"- **{r['model']}**: R²_in={r['R2_in']:.3f}; R²_cv5={r['R2_cv5']:.3f}; n={r['n']}")
    lines += ["", "## ΔR² bootstrap", f"`{json.dumps(deltas, indent=2)}`", "", "## Spearman", f"`{json.dumps(spear, indent=2)}`", ""]
    (OUT / "isles24_cta_features_report.md").write_text("\n".join(lines), encoding="utf-8")
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
