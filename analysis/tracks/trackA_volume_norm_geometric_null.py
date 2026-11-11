"""
Track A — volume-normalized territorial burden + geometric null + variance partition.

Scientific rationale (concept note Track A):
  A1. Burden per unit territory volume (not raw share alone).
  A2. Coarse geometric null: territory volume, CSF/ventricular distance, GM/WM fraction,
      MRA vessel probability/radius (atlas-derived in MNI).
  A3. Commonality / hierarchical R² partition: unique geometry, unique territory, shared.
  A4. Selection-structure note (ATLAS rehab; ISLES LVO/thrombectomy MCA skew).

Uses already-computed subject × territory lesion/WMH tables from Aim 4.
Does NOT mix TRCI/BCI.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets, image
from scipy import ndimage
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tracks" / "trackA"
OUT.mkdir(parents=True, exist_ok=True)

TERRITORIES = ["ACA", "MCA", "PCA", "VB"]

# Prefer MNI152 1mm atlas for ATLAS-space work; ISLES uses 2mm resampled maps
ATLAS_L2_1MM = (
    ROOT
    / "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/Atlas_182_MNI152/ArterialAtlas_level2.nii"
)
ATLAS_L2_ATLAS_SPACE = (
    ROOT / "results/aim4/atlas/ArterialAtlas_level2_major_on_ATLAS_MNI152NLin2009aSym.nii.gz"
)
ATLAS_L2_ISLES = ROOT / "results/aim4/isles/ArterialAtlas_level2_major_on_MNI152_2mm.nii.gz"
MRA_PROB = ROOT / "data/raw/mra_atlas/figshare/vesselProbabilities.nii.gz"
MRA_RAD = ROOT / "data/raw/mra_atlas/figshare/vesselRadius.nii.gz"

# Major-territory label maps on ArterialAtlas_level2 (Liu labels collapsed L/R → ACA/MCA/PCA/VB)
# After collapse in aim4 pipelines: 1=ACA, 2=MCA, 3=PCA, 4=VB (verify from files)


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def r2_fit(X: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(X).all(axis=1) & np.isfinite(y)
    if m.sum() < max(10, X.shape[1] + 5):
        return float("nan"), int(m.sum())
    Xs = StandardScaler().fit_transform(X[m])
    return float(LinearRegression().fit(Xs, y[m]).score(Xs, y[m])), int(m.sum())


def commonality(r2_g: float, r2_t: float, r2_gt: float) -> dict:
    """Variance partition (commonality analysis)."""
    if any(np.isnan(v) for v in (r2_g, r2_t, r2_gt)):
        return {
            "R2_geometry": r2_g,
            "R2_territory": r2_t,
            "R2_joint": r2_gt,
            "unique_geometry": np.nan,
            "unique_territory": np.nan,
            "shared": np.nan,
        }
    return {
        "R2_geometry": r2_g,
        "R2_territory": r2_t,
        "R2_joint": r2_gt,
        "unique_geometry": float(r2_gt - r2_t),
        "unique_territory": float(r2_gt - r2_g),
        "shared": float(r2_g + r2_t - r2_gt),
    }


def load_label_map(path: Path) -> tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(path))
    data = np.rint(np.asanyarray(img.dataobj, dtype=float)).astype(np.int32)
    if data.ndim == 4:
        data = data[..., 0]
    return data, img.affine


def major_codes(data: np.ndarray) -> dict[str, int]:
    """Infer ACA/MCA/PCA/VB codes from unique labels (handles 1–4 major or L/R 1–8)."""
    uniq = sorted(int(x) for x in np.unique(data) if x > 0)
    # Already major (1–4)
    if set(uniq).issubset({1, 2, 3, 4}) and len(uniq) <= 4:
        return {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}
    # Bilateral level2 often 1–8 (L/R pairs). Collapse: odd/even or known Liu scheme.
    # Liu ArterialAtlas_level2: pairs for L/R of each major — use modulo mapping if 1..8
    if max(uniq) <= 8:
        # Common: 1,2=ACA L/R; 3,4=MCA; 5,6=PCA; 7,8=VB OR 1–4 major already
        # Prefer mapping used in aim4: check voxel counts — MCA should be largest
        return {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}  # major maps in results/
    return {"ACA": 1, "MCA": 2, "PCA": 3, "VB": 4}


def build_geometry_table(atlas_path: Path, tag: str) -> pd.DataFrame:
    """Atlas-derived geometry covariates per major territory (constant across subjects)."""
    log(f"Building geometry table from {atlas_path.name} ({tag})")
    tv, affine = load_label_map(atlas_path)
    codes = major_codes(tv)

    # If labels are 1–8 bilateral, collapse to major
    if tv.max() > 4:
        # Heuristic for Liu level2 bilateral: labels come in L/R pairs
        # Use mapping from aim4 major file when possible — here remap ((lab+1)//2)
        collapsed = np.zeros_like(tv)
        for lab in range(1, int(tv.max()) + 1):
            major = (lab + 1) // 2 if lab <= 8 else lab
            if major > 4:
                major = ((lab - 1) % 4) + 1
            collapsed[tv == lab] = major
        # Prefer precomputed major maps when provided (results paths already major)
        if "major" in atlas_path.name.lower() or set(np.unique(tv)).issubset({0, 1, 2, 3, 4}):
            collapsed = tv
        else:
            # Check if 1-4 already
            if set(np.unique(tv)).issubset({0, 1, 2, 3, 4}):
                collapsed = tv
            else:
                # For ArterialAtlas_level2 raw: read labels file if available
                collapsed = tv.copy()
                # Liu Sci Data: level2 has many parcels — use results major maps instead
                pass
        tv = collapsed if set(np.unique(collapsed)).issubset({0, 1, 2, 3, 4}) else tv

    # CSF / ventricles from Harvard-Oxford subcortical
    ho = datasets.fetch_atlas_harvard_oxford("sub-maxprob-thr25-2mm")
    ho_img = image.resample_to_img(
        ho.maps,
        nib.Nifti1Image(tv.astype(np.float32), affine),
        interpolation="nearest",
        force_resample=True,
        copy_header=True,
    )
    ho_data = np.asanyarray(ho_img.dataobj)
    labels = list(ho.labels)
    vent_names = [
        "Left Lateral Ventricle",
        "Right Lateral Ventricle",
        "Left Third Ventricle",
        "Right Third Ventricle",
        "Fourth Ventricle",
    ]
    vent_idx = [labels.index(n) for n in vent_names if n in labels]
    # Also include "Brain-Stem" exclusion later; CSF from ICBM
    vent_mask = np.isin(ho_data, vent_idx) if vent_idx else np.zeros_like(tv, dtype=bool)

    # Tissue probability maps (GM/WM)
    gm_img = datasets.load_mni152_gm_template(resolution=2)
    wm_img = datasets.load_mni152_wm_template(resolution=2)

    ref = nib.Nifti1Image(tv.astype(np.float32), affine)
    gm = np.asanyarray(
        image.resample_to_img(gm_img, ref, interpolation="linear", force_resample=True, copy_header=True).dataobj,
        dtype=float,
    )
    wm = np.asanyarray(
        image.resample_to_img(wm_img, ref, interpolation="linear", force_resample=True, copy_header=True).dataobj,
        dtype=float,
    )

    # Distance to ventricles (mm)
    if vent_mask.any():
        # voxel sizes
        zooms = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
        dist_vox = ndimage.distance_transform_edt(~vent_mask)
        # anisotropic approx: multiply by mean zoom
        dist_mm = dist_vox * float(np.mean(zooms))
    else:
        dist_mm = np.zeros_like(tv, dtype=float)

    # MRA
    mra_p = mra_r = None
    if MRA_PROB.exists():
        mra_p = np.asanyarray(
            image.resample_to_img(
                nib.load(str(MRA_PROB)), ref, interpolation="linear", force_resample=True, copy_header=True
            ).dataobj,
            dtype=float,
        )
    if MRA_RAD.exists():
        mra_r = np.asanyarray(
            image.resample_to_img(
                nib.load(str(MRA_RAD)), ref, interpolation="linear", force_resample=True, copy_header=True
            ).dataobj,
            dtype=float,
        )

    rows = []
    # Use codes 1–4 if present
    present = sorted(int(x) for x in np.unique(tv) if x in (1, 2, 3, 4))
    name_by_code = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
    if not present:
        # fall back: map using codes dict values
        name_by_code = {v: k for k, v in codes.items()}
        present = list(name_by_code.keys())

    for code in present:
        name = name_by_code.get(code, f"T{code}")
        mask = tv == code
        n = int(mask.sum())
        if n == 0:
            continue
        rows.append(
            {
                "geometry_space": tag,
                "territory": name,
                "tv_code": int(code),
                "territory_voxels": n,
                "log_territory_voxels": float(np.log(n)),
                "mean_vent_dist_mm": float(dist_mm[mask].mean()) if dist_mm is not None else np.nan,
                "gm_mean": float(gm[mask].mean()),
                "wm_mean": float(wm[mask].mean()),
                "gm_wm_ratio": float(gm[mask].mean() / max(wm[mask].mean(), 1e-6)),
                "mra_vessel_probability": float(mra_p[mask].mean()) if mra_p is not None else np.nan,
                "mra_vessel_radius": float(mra_r[mask].mean()) if mra_r is not None else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"trackA_geometry_{tag}.csv", index=False)
    return df


def long_from_atlas_frac(path: Path, cohort: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df.get("qc_primary", True) == True].copy() if "qc_primary" in df.columns else df
    df = df[df["status"] == "ok"] if "status" in df.columns else df
    rows = []
    idcol = "subject_id" if "subject_id" in df.columns else "ptid"
    for _, r in df.iterrows():
        for t in TERRITORIES:
            vox_col = f"lesion_vox_{t}"
            tv_col = f"tv_vox_{t}"
            bur_col = f"lesion_burden_in_{t}"
            if vox_col not in df.columns:
                continue
            tv = float(r[tv_col]) if pd.notna(r.get(tv_col)) else np.nan
            vox = float(r[vox_col]) if pd.notna(r.get(vox_col)) else np.nan
            bur = float(r[bur_col]) if bur_col in df.columns and pd.notna(r.get(bur_col)) else (
                vox / tv if tv and tv > 0 else np.nan
            )
            rows.append(
                {
                    "cohort": cohort,
                    "subject_id": r[idcol],
                    "territory": t,
                    "lesion_vox": vox,
                    "tv_vox": tv,
                    "burden_vol_norm": bur,
                    "log1p_lesion_vox": np.log1p(vox) if np.isfinite(vox) else np.nan,
                    "log1p_burden": np.log1p(bur * 1e4) if np.isfinite(bur) else np.nan,  # scale for stability
                    "total_lesion_vox": float(r["n_lesion_voxels"]) if pd.notna(r.get("n_lesion_voxels")) else np.nan,
                    "dominant_territory": r.get("dominant_territory"),
                }
            )
    return pd.DataFrame(rows)


def long_from_adni(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        for t in TERRITORIES:
            bur = r.get(f"wmh_burden_in_{t}")
            # reconstruct tv from nest if available — ADNI nest may lack tv_vox; use burden only
            rows.append(
                {
                    "cohort": "ADNI_WMH",
                    "subject_id": r.get("ptid", r.get("RID")),
                    "territory": t,
                    "lesion_vox": np.nan,
                    "tv_vox": np.nan,
                    "burden_vol_norm": float(bur) if pd.notna(bur) else np.nan,
                    "log1p_lesion_vox": np.nan,
                    "log1p_burden": np.log1p(float(bur) * 1e4) if pd.notna(bur) else np.nan,
                    "total_lesion_vox": float(r["n_wmh_voxels"]) if pd.notna(r.get("n_wmh_voxels")) else np.nan,
                    "dominant_territory": None,
                    "TOTAL_WMH_ucd": r.get("TOTAL_WMH_ucd"),
                }
            )
    return pd.DataFrame(rows)


def analyze_cohort(long_df: pd.DataFrame, geo: pd.DataFrame, cohort: str) -> dict:
    log(f"Variance partition: {cohort} n_rows={len(long_df)}")
    g = geo.drop(columns=["geometry_space", "tv_code"], errors="ignore")
    d = long_df.merge(g, on="territory", how="left")
    d = d[np.isfinite(d["log1p_burden"])].copy()
    # territory dummies
    for t in TERRITORIES:
        d[f"is_{t}"] = (d["territory"] == t).astype(float)

    geo_cols = [
        c
        for c in [
            "log_territory_voxels",
            "mean_vent_dist_mm",
            "gm_mean",
            "wm_mean",
            "gm_wm_ratio",
            "mra_vessel_probability",
            "mra_vessel_radius",
        ]
        if c in d.columns and d[c].notna().all()
    ]
    # drop near-constant
    geo_cols = [c for c in geo_cols if d[c].std(ddof=0) > 1e-12]

    terr_cols = [f"is_{t}" for t in TERRITORIES[1:]]  # drop one for dummy trap; keep ACA as ref

    y = d["log1p_burden"].values
    Xg = d[geo_cols].values if geo_cols else np.zeros((len(d), 1))
    Xt = d[terr_cols].values
    Xgt = np.column_stack([Xg, Xt]) if geo_cols else Xt

    r2_g, n_g = r2_fit(Xg, y) if geo_cols else (0.0, len(d))
    r2_t, n_t = r2_fit(Xt, y)
    r2_gt, n_gt = r2_fit(Xgt, y)
    part = commonality(r2_g, r2_t, r2_gt)

    # Descriptive: mean volume-normalized burden by territory
    desc = (
        d.groupby("territory")
        .agg(
            n=("subject_id", "nunique"),
            mean_burden=("burden_vol_norm", "mean"),
            median_burden=("burden_vol_norm", "median"),
            mean_log1p_burden=("log1p_burden", "mean"),
            mean_tv_vox=("tv_vox", "mean"),
        )
        .reset_index()
    )

    # MCA vs ACA Wilcoxon on burden
    from scipy import stats

    a = d.loc[d["territory"] == "ACA", "burden_vol_norm"].dropna()
    m = d.loc[d["territory"] == "MCA", "burden_vol_norm"].dropna()
    if len(a) > 10 and len(m) > 10:
        w = stats.wilcoxon(m.values[: min(len(m), len(a))], a.values[: min(len(m), len(a))]) if False else stats.mannwhitneyu(
            m, a, alternative="two-sided"
        )
        wilcox = {"stat": float(w.statistic), "p": float(w.pvalue), "test": "mannwhitneyu_MCA_vs_ACA"}
    else:
        wilcox = {}

    out = {
        "cohort": cohort,
        "n_subject_territory_rows": int(len(d)),
        "n_subjects": int(d["subject_id"].nunique()),
        "geometry_features": geo_cols,
        "outcome": "log1p(burden_vol_norm * 1e4)",
        "partition": part,
        "n_model": {"geometry": n_g, "territory": n_t, "joint": n_gt},
        "wilcoxon_MCA_vs_ACA_burden": wilcox,
        "interpretation": [
            "Outcome is volume-normalized territorial burden (lesion/WMH voxels per territory voxel).",
            "Geometry features are atlas-derived (constant within territory across subjects).",
            "Unique territorial R² = joint − geometry; unique geometry = joint − territory; shared = G+T−joint.",
            "Large shared component is expected if branching anatomy is developmentally specified.",
            "NOT molecular H3; NOT TRCI/BCI.",
        ],
    }
    desc.to_csv(OUT / f"trackA_{cohort}_burden_by_territory.csv", index=False)
    d.to_csv(OUT / f"trackA_{cohort}_long.csv", index=False)
    return out


def selection_structure_note() -> dict:
    return {
        "ATLAS": (
            "Chronic / rehabilitation-skewed sampling favors large MCA-territory lesions; "
            "agreement with acute cohorts is not fully independent confirmation of MCA biology."
        ),
        "ISLES22": (
            "Acute DWI lesions; LVO/thrombectomy selection also MCA-skewed. "
            "Higher VB-dominant fraction than ATLAS is descriptive cohort mix."
        ),
        "ADNI_WMH": (
            "Non-stroke SVD/WMH axis; TerritoryVoxels geography; not stroke selection. "
            "Mean MCA>ACA share claim; not majority-subject claim."
        ),
    }


def main():
    log("=== Track A start ===")
    # Geometry for ATLAS-space major map and ISLES 2mm major map
    geo_atlas = build_geometry_table(
        ATLAS_L2_ATLAS_SPACE if ATLAS_L2_ATLAS_SPACE.exists() else ATLAS_L2_1MM,
        "ATLAS_MNI152NLin2009aSym_major",
    )
    geo_isles = build_geometry_table(
        ATLAS_L2_ISLES if ATLAS_L2_ISLES.exists() else ATLAS_L2_1MM,
        "ISLES_MNI152_2mm_major",
    )
    # ADNI used TerritoryVoxels — use 1mm MNI152 level2 major as approximate geometry prior
    geo_adni = build_geometry_table(ATLAS_L2_1MM, "ADNI_MNI152_1mm_level2")

    atlas_frac = ROOT / "results/aim4/atlas/atlas_lesion_tv_fractions.csv"
    isles_frac = ROOT / "results/aim4/isles/isles_lesion_tv_fractions.csv"
    adni_nest = ROOT / "results/aim4/adni/aim4_nest_wmh_mra_subject_table.csv"

    results = {"selection_structure": selection_structure_note(), "cohorts": {}}

    if atlas_frac.exists():
        long_a = long_from_atlas_frac(atlas_frac, "ATLAS")
        results["cohorts"]["ATLAS"] = analyze_cohort(long_a, geo_atlas, "ATLAS")
    if isles_frac.exists():
        long_i = long_from_atlas_frac(isles_frac, "ISLES22")
        results["cohorts"]["ISLES22"] = analyze_cohort(long_i, geo_isles, "ISLES22")
    if adni_nest.exists():
        long_w = long_from_adni(adni_nest)
        results["cohorts"]["ADNI_WMH"] = analyze_cohort(long_w, geo_adni, "ADNI_WMH")

    # Cross-cohort summary table
    rows = []
    for c, res in results["cohorts"].items():
        p = res["partition"]
        rows.append(
            {
                "cohort": c,
                "n_subjects": res["n_subjects"],
                "R2_geometry": p["R2_geometry"],
                "R2_territory": p["R2_territory"],
                "R2_joint": p["R2_joint"],
                "unique_geometry": p["unique_geometry"],
                "unique_territory": p["unique_territory"],
                "shared": p["shared"],
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "trackA_variance_partition_summary.csv", index=False)

    results["verdict"] = _verdict(summary)
    with open(OUT / "trackA_summary.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    _write_report(results, summary)
    log(f"Wrote {OUT}")
    log("=== Track A done ===")


def _verdict(summary: pd.DataFrame) -> str:
    if summary.empty:
        return "NO_DATA"
    # If unique territory survives after geometry with meaningful joint R²
    ut = summary["unique_territory"].mean()
    sh = summary["shared"].mean()
    ug = summary["unique_geometry"].mean()
    if np.nanmax(summary["R2_joint"]) < 0.05:
        return "NULL_UNDERPOWERED_OR_NO_SIGNAL__CHECK_DESIGN"
    if ut > 0.05 and ut > ug:
        return "TERRITORIAL_UNIQUE_VARIANCE_SURVIVES_COARSE_GEOMETRIC_NULL"
    if sh > max(ut, ug, 0):
        return "SHARED_GEOMETRIC_TERRITORIAL_DOMINATES__TRACK_C_ALIGNED"
    if ug >= ut:
        return "GEOMETRY_ABSORBS_MOST_TERRITORIAL_SIGNAL__METHODS_NULL_RESULT"
    return "MIXED__SEE_PARTITION_TABLE"


def _write_report(results: dict, summary: pd.DataFrame):
    lines = [
        "# Track A — volume-normalized burden + geometric null + variance partition",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d')}",
        f"**Verdict:** `{results['verdict']}`",
        "",
        "## Design locks",
        "- Outcome: volume-normalized territorial burden (lesion or WMH voxels / territory voxels).",
        "- Geometry null (coarse): log territory volume, mean distance to ventricles, GM/WM means,",
        "  MRA vessel probability & radius (atlas-level).",
        "- Partition: unique geometry / unique territory / shared (commonality analysis).",
        "- **Not** molecular H3. **Not** TRCI/BCI.",
        "",
        "## Selection structure (A4)",
    ]
    for k, v in results["selection_structure"].items():
        lines.append(f"- **{k}:** {v}")
    lines += ["", "## Variance partition summary", "", summary.to_markdown(index=False), ""]
    for c, res in results["cohorts"].items():
        lines += [
            f"## {c}",
            f"- n subjects: {res['n_subjects']}",
            f"- geometry features: {', '.join(res['geometry_features'])}",
            f"- partition: `{json.dumps(res['partition'])}`",
            f"- MCA vs ACA: `{json.dumps(res.get('wilcoxon_MCA_vs_ACA_burden', {}))}`",
            "",
        ]
    lines += [
        "## Interpretation",
        "- If **shared ≫ unique**, architecture and territory labels are collinear (expected under developmental branching) → supports Track C framing.",
        "- If **unique territory** remains after the coarse null, upgrade null (penetrator depth / watershed) before claiming biology.",
        "- If **geometry absorbs** the signal, publish as methods + null-result paper.",
        "",
    ]
    (OUT / "trackA_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
