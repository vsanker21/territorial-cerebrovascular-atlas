"""
ADNI WMH z-threshold sensitivity on existing FLAIR→T1→MNI warps.

Re-segments from flair_mni.nii.gz (no re-registration). Applies the same QC rules
as qc_adni_wmh_t1_mni.py. Primary lock remains z=1.5; this report is sensitivity only.

Does not claim molecular H1 / H3.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import SimpleITK as sitk
from nilearn import datasets, image
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis/aim4"))
from adni_wmh_t1_mni_pipeline import TV_MAP, TV_PATH  # noqa: E402

OUT = ROOT / "results/aim4/adni"
WORK = ROOT / "data/processed/adni/t1_mni"
SCALED = OUT / "adni_wmh_t1_mni_scaled.csv"
QC = OUT / "adni_wmh_t1_mni_qc.csv"
INV = OUT / "adni_aim4_cohort_inventory.csv"
MRA = ROOT / "results/aim4/mra_geometry_by_territory.csv"

Z_GRID = (1.0, 1.25, 1.5, 1.75, 2.0)
PRIMARY_Z = 1.5
WM_PCT = 99.0


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def r2_ols(X, y):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(X).all(axis=1) & np.isfinite(y)
    if mask.sum() < 8:
        return np.nan, int(mask.sum())
    Xs = StandardScaler().fit_transform(X[mask])
    return float(LinearRegression().fit(Xs, y[mask]).score(Xs, y[mask])), int(mask.sum())


def apply_qc_row(n_wmh, frac_tv, ucd, flair_nz, t1_nz, brain_dice, brain_nvox):
    reasons = []
    flags = []
    if n_wmh <= 0:
        reasons.append("FAIL_ZERO_WMH")
    if flair_nz < 0.15 * brain_nvox or t1_nz < 0.15 * brain_nvox:
        reasons.append("FAIL_EMPTY_REG")
    if brain_dice == brain_dice and brain_dice < 0.25:
        reasons.append("FAIL_BRAIN_DICE")
    if ucd == ucd and ucd >= 10.0 and n_wmh < 20:
        reasons.append("FAIL_UCD_DISCORD")
    if n_wmh > 0 and frac_tv < 0.10:
        flags.append("FLAG_LOW_TV")
    return reasons, flags


def morph_open(wmh: np.ndarray) -> np.ndarray:
    arr = wmh.astype(np.uint8).transpose(2, 1, 0)
    simg = sitk.GetImageFromArray(arr)
    simg = sitk.BinaryMorphologicalOpening(simg, [1, 1, 1])
    return sitk.GetArrayFromImage(simg).transpose(2, 1, 0).astype(bool)


def tv_fractions_cached(wmh: np.ndarray, tv: np.ndarray) -> dict:
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


def nest_metrics(df: pd.DataFrame, dem: pd.DataFrame):
    d = df.merge(dem, on="ptid", how="left")
    scan_year = pd.to_datetime(d["image_date"], errors="coerce").dt.year
    dob = pd.to_datetime(d["PTDOBYY"], errors="coerce")
    birth_year = dob.dt.year.fillna(pd.to_numeric(d["PTDOBYY"], errors="coerce"))
    d["age_approx"] = scan_year - birth_year
    d["sex_f"] = (pd.to_numeric(d["PTGENDER"], errors="coerce") == 2).astype(float)
    y = d["TOTAL_WMH_ucd"].astype(float).to_numpy()
    burden = [f"wmh_burden_in_{t}" for t in ["ACA", "MCA", "PCA", "VB"]]
    r0, _ = r2_ols(d[["age_approx", "sex_f"]].to_numpy(), y)
    r3, n3 = r2_ols(d[["age_approx", "sex_f"] + burden].to_numpy(), y)
    rho, p = stats.spearmanr(d["TOTAL_WMH_ucd"], d["n_wmh_voxels"])
    shares = {t: float(d[f"wmh_frac_of_wmh_{t}"].mean()) for t in ["ACA", "MCA", "PCA", "VB"]}
    return {
        "n": int(len(d)),
        "spearman_rho": float(rho) if rho == rho else None,
        "spearman_p": float(p) if p == p else None,
        "r2_M0": r0,
        "r2_M3": r3,
        "delta_M3_M0": (r3 - r0) if (r3 == r3 and r0 == r0) else None,
        "mean_share": shares,
        "n_ols": n3,
    }


def main():
    scaled = pd.read_csv(SCALED)
    qc0 = pd.read_csv(QC) if QC.exists() else None
    dice_map = {}
    if qc0 is not None and "brain_dice_flair" in qc0.columns:
        dice_map = dict(zip(qc0["ptid"].astype(str), qc0["brain_dice_flair"].astype(float)))

    inv = pd.read_csv(INV)
    inv["ptid"] = inv["subject_id"].astype(str)
    dem = inv[["ptid", "PTGENDER", "PTDOBYY", "image_date", "TOTAL_WMH"]].drop_duplicates("ptid")

    mra = pd.read_csv(MRA)
    geo = {}
    for _, r in mra.iterrows():
        geo[f"mra_{r['metric']}_{r['territory']}"] = float(r["mean"])
    geo_feat = [c for c in geo if c.startswith("mra_vessel_")]

    mni = datasets.load_mni152_template(resolution=2)
    wm = datasets.load_mni152_wm_mask(threshold=0.5, resolution=2)
    brain = datasets.load_mni152_brain_mask(threshold=0.5, resolution=2)
    brain_nvox = int(np.asanyarray(brain.dataobj).astype(bool).sum())
    tv_img = nib.load(str(TV_PATH))

    detail_rows = []
    n_subj = len(scaled)
    for i, r in scaled.iterrows():
        ptid = str(r["ptid"])
        flair_p = WORK / ptid / "flair_mni.nii.gz"
        if not flair_p.exists():
            log(f"SKIP missing flair_mni {ptid}")
            continue
        flair_nz = float(r.get("qc_flair_mni_nonzero") or 0)
        t1_nz = float(r.get("qc_t1_mni_nonzero") or 0)
        ucd = float(r["TOTAL_WMH_ucd"]) if pd.notna(r.get("TOTAL_WMH_ucd")) else np.nan
        dice = float(dice_map.get(ptid, np.nan))
        if (i + 1) % 20 == 0 or i == 0:
            log(f"[{i+1}/{n_subj}] {ptid}")
        try:
            f_img = image.load_img(str(flair_p))
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
                        image.resample_to_img(
                            brain, f_img, interpolation="nearest", force_resample=True, copy_header=True
                        ).dataobj
                    )
                    > 0
                )
            mask = w & b & np.isfinite(f) & (f > 0)
            vals = f[mask]
            if vals.size < 500:
                raise RuntimeError(f"WM voxels too few: {vals.size}")
            tv_r = image.resample_to_img(tv_img, f_img, interpolation="nearest", force_resample=True, copy_header=True)
            tv = np.rint(np.asanyarray(tv_r.dataobj, dtype=float)).astype(np.int32)
            if tv.ndim == 4:
                tv = tv[..., 0]
            mu, sd = float(vals.mean()), float(vals.std())
            p99 = float(np.percentile(vals, WM_PCT))
        except Exception as e:
            for z in Z_GRID:
                detail_rows.append(
                    {"ptid": ptid, "z_thr": z, "status": "error", "error": str(e), "TOTAL_WMH_ucd": ucd}
                )
            continue

        for z in Z_GRID:
            thr = min(mu + float(z) * sd, p99)
            wmh = morph_open((f >= thr) & mask)
            fr = tv_fractions_cached(wmh, tv)
            reasons, flags = apply_qc_row(
                fr["n_wmh_voxels"], fr["wmh_frac_assigned_tv"], ucd, flair_nz, t1_nz, dice, brain_nvox
            )
            detail_rows.append(
                {
                    "ptid": ptid,
                    "RID": r.get("RID"),
                    "z_thr": z,
                    "status": "ok",
                    "seg_thr": thr,
                    "TOTAL_WMH_ucd": ucd,
                    "brain_dice_flair": dice,
                    "qc_flair_mni_nonzero": flair_nz,
                    "qc_t1_mni_nonzero": t1_nz,
                    "qc_reasons": ";".join(reasons),
                    "qc_flags": ";".join(flags),
                    "qc_primary": len(reasons) == 0,
                    "qc_sensitivity_strict": len(reasons) == 0 and "FLAG_LOW_TV" not in flags,
                    **fr,
                    **{k: geo[k] for k in geo_feat},
                }
            )

    detail = pd.DataFrame(detail_rows)
    detail_path = OUT / "adni_wmh_zthr_robustness_detail.csv"
    detail.to_csv(detail_path, index=False)

    summary_rows = []
    nest_by_z = {}
    for z in Z_GRID:
        sub = detail[(detail["z_thr"] == z) & (detail["status"] == "ok")].copy()
        n_all = len(sub)
        n_zero = int((sub["n_wmh_voxels"] <= 0).sum())
        primary = sub[sub["qc_primary"] == True].copy()  # noqa: E712
        strict = sub[sub["qc_sensitivity_strict"] == True].copy()  # noqa: E712
        rho_all, p_all = stats.spearmanr(sub["TOTAL_WMH_ucd"], sub["n_wmh_voxels"])
        nest_p = nest_metrics(primary, dem) if len(primary) >= 8 else None
        nest_s = nest_metrics(strict, dem) if len(strict) >= 8 else None
        summary_rows.append(
            {
                "z_thr": z,
                "is_primary_lock": z == PRIMARY_Z,
                "n_subjects": n_all,
                "n_zero_wmh": n_zero,
                "frac_zero_wmh": n_zero / n_all if n_all else None,
                "n_pass_primary": int(len(primary)),
                "n_pass_strict": int(len(strict)),
                "median_wmh_voxels_all": float(sub["n_wmh_voxels"].median()) if n_all else None,
                "median_wmh_voxels_primary": float(primary["n_wmh_voxels"].median()) if len(primary) else None,
                "spearman_all_rho": float(rho_all) if rho_all == rho_all else None,
                "spearman_all_p": float(p_all) if p_all == p_all else None,
                "nest_primary": nest_p,
                "nest_strict": nest_s,
            }
        )
        nest_by_z[str(z)] = {"primary": nest_p, "strict": nest_s}

    flat = []
    for row in summary_rows:
        flat.append(
            {
                "z_thr": row["z_thr"],
                "is_primary_lock": row["is_primary_lock"],
                "n_subjects": row["n_subjects"],
                "n_zero_wmh": row["n_zero_wmh"],
                "frac_zero_wmh": row["frac_zero_wmh"],
                "n_pass_primary": row["n_pass_primary"],
                "n_pass_strict": row["n_pass_strict"],
                "median_wmh_voxels_all": row["median_wmh_voxels_all"],
                "median_wmh_voxels_primary": row["median_wmh_voxels_primary"],
                "spearman_all_rho": row["spearman_all_rho"],
                "spearman_primary_rho": (row["nest_primary"] or {}).get("spearman_rho"),
                "delta_M3_M0_primary": (row["nest_primary"] or {}).get("delta_M3_M0"),
                "delta_M3_M0_strict": (row["nest_strict"] or {}).get("delta_M3_M0"),
                "mean_share_MCA_primary": ((row["nest_primary"] or {}).get("mean_share") or {}).get("MCA"),
                "mean_share_ACA_primary": ((row["nest_primary"] or {}).get("mean_share") or {}).get("ACA"),
            }
        )
    flat_df = pd.DataFrame(flat)
    flat_df.to_csv(OUT / "adni_wmh_zthr_robustness_summary.csv", index=False)

    lock = next(x for x in summary_rows if x["z_thr"] == PRIMARY_Z)
    payload = {
        "primary_z_lock": PRIMARY_Z,
        "z_grid": list(Z_GRID),
        "n_subjects": int(len(scaled)),
        "summary": flat,
        "nest_by_z": nest_by_z,
        "interpretation_locks": [
            "Primary analysis remains z_thr=1.5 + QC primary PASS (locked before scaling).",
            "Lower z increases sensitivity (fewer zero-WMH) but risks false-positive hyperintensities.",
            "Higher z decreases sensitivity; more FAIL_ZERO_WMH / UCD discordance.",
            "M3−M0 remains a partition/concordance check of local TV WMH vs UCD TOTAL_WMH, not H3 molecular proof.",
            "Does not unlock enrichment-matched Class B DE or ATLAS stroke Aim 4.",
            "Class B claim unchanged: donor OOF separability yes; gene-level territorial markers not from unsorted ACC vs BVE MCA.",
        ],
        "crosscheck_z15_vs_locked_qc": {
            "robustness_n_pass_primary": lock["n_pass_primary"],
            "robustness_n_pass_strict": lock["n_pass_strict"],
            "robustness_n_zero": lock["n_zero_wmh"],
        },
    }
    (OUT / "adni_wmh_zthr_robustness_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "# ADNI WMH z-threshold robustness",
        "",
        f"**Primary lock:** z_thr = {PRIMARY_Z} (unchanged).",
        f"**Grid:** {list(Z_GRID)} on existing `flair_mni.nii.gz` (no re-registration).",
        "",
        "## Summary",
        flat_df.to_string(index=False),
        "",
        "## Locked interpretation",
        "- Primary Aim 4 non-stroke table remains QC-pass at z=1.5.",
        "- Zero-WMH fraction vs UCD≥3 cc is a known sensitivity limit of the exploratory WM-prior z-score segmenter.",
        "- Concordance (Spearman local vs UCD) and MCA-dominant share should be checked for qualitative stability across z.",
        "- M3−M0 ΔR² is partition/concordance, **not** molecular H3.",
        "- Gene-level Class B territorial markers remain blocked (unsorted ACC vs BVE MCA); do not re-run soft unmatched DE.",
        "",
        f"Detail: `{detail_path.relative_to(ROOT)}`",
        "",
    ]
    (OUT / "adni_wmh_zthr_robustness_report.md").write_text("\n".join(lines), encoding="utf-8")
    log(json.dumps({"wrote": str(OUT / "adni_wmh_zthr_robustness_summary.csv"), **payload["crosscheck_z15_vs_locked_qc"]}, indent=2))


if __name__ == "__main__":
    main()
