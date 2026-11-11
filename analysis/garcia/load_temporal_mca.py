"""
Load Garcia Nature 2022 ex vivo temporal (MCA Class B) cerebrovascular nuclei.

Scientific framing
------------------
- Ex vivo epilepsy temporal tip/lobe samples are TerritoryVoxels MCA (P≈1.0).
- Join GEO/MIT barcode meta (SampleID) to MOESM2 Ex Vivo Sample ID table.
- Attach TerritoryVoxels MCA stream; mark analysis_stream = ClassB_MCA.

Outputs
-------
data/processed/garcia/garcia_temporal_mca.h5ad
data/processed/garcia/garcia_temporal_mca_sample_table.csv
results/garcia/garcia_temporal_mca_report.md
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyreadr
import rpy2.robjects as ro
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "feasibility"))
sys.path.insert(0, str(ROOT / "analysis" / "yang"))

MOESM2 = ROOT / "data/raw/garcia/esm/41586_2022_4521_MOESM2_ESM.xlsx"
COUNTS = ROOT / "data/raw/garcia/mit/brain.BBB.human.counts.GEOsubmission.rds"
META = ROOT / "data/raw/garcia/mit/brain.BBB.human.meta.GEOsubmission.rds"
# Fallback gzipped GEO
COUNTS_GZ = ROOT / "data/raw/garcia/GSE173731_brain.BBB.human.counts.rds.gz"
META_GZ = ROOT / "data/raw/garcia/GSE173731_brain.BBB.human.meta.rds.gz"

OUT_H5AD = ROOT / "data/processed/garcia/garcia_temporal_mca.h5ad"
OUT_SAMPLES = ROOT / "data/processed/garcia/garcia_temporal_mca_sample_table.csv"
OUT_MD = ROOT / "results/garcia/garcia_temporal_mca_report.md"
TV_LOOKUP = ROOT / "data/processed/garcia/garcia_class_b_territory_lookup.json"

# Temporal pole TerritoryVoxels (from confirmatory interpretation): MCA 100%
TEMPORAL_TV = {
    "roi": "Temporal Pole",
    "method": "TerritoryVoxels_BMM_labels_1to4",
    "analysis_stream": "ClassB_MCA",
    "territory_contrast": "MCA_temporal",
    "tv_dominant": "MCA",
    "tv_p_ACA": 0.0,
    "tv_p_MCA": 1.0,
    "tv_p_PCA": 0.0,
    "tv_p_VB": 0.0,
    "tv_p_dominant": 1.0,
    "tv_assignment_fraction": 0.43,
    "tv_aca_mca_borderzone_flag": False,
    "source_note": "NITRC confirmatory ROI probe; temporal pole among assigned = MCA 100%",
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def load_exvivo_table() -> pd.DataFrame:
    ex = pd.read_excel(MOESM2, sheet_name="Ex Vivo Sample ID")
    # Normalize column names
    ex.columns = [str(c).strip() for c in ex.columns]
    # Expected: BVE Sequencing ID, Tissue location, ...
    id_col = [c for c in ex.columns if "BVE" in c.upper() or c == "SampleID"][0]
    loc_col = [c for c in ex.columns if "tissue" in c.lower() or "location" in c.lower()][0]
    ex = ex.rename(columns={id_col: "SampleID", loc_col: "tissue_location"})
    ex["SampleID"] = ex["SampleID"].astype(str).str.strip()
    ex["tissue_location"] = ex["tissue_location"].astype(str)
    # Keep temporal samples only
    mask = ex["tissue_location"].str.contains("temporal", case=False, na=False)
    ex = ex.loc[mask].copy()
    ex["region"] = "temporal_cortex"
    ex["class_b_role"] = "MCA_primary"
    return ex


def read_meta() -> pd.DataFrame:
    path = META if META.exists() else None
    if path is None and META_GZ.exists():
        import gzip
        import shutil

        tmp = ROOT / "data/raw/garcia/_tmp_meta.rds"
        with gzip.open(META_GZ, "rb") as f, open(tmp, "wb") as o:
            shutil.copyfileobj(f, o)
        path = tmp
    meta = pyreadr.read_r(str(path))[None]
    if meta.index.name or meta.index.dtype != "int64":
        meta = meta.copy()
        meta["barcode"] = meta.index.astype(str)
    else:
        meta["barcode"] = meta.index.astype(str)
    meta["SampleID"] = meta["SampleID"].astype(str)
    return meta


def read_counts_subset(keep_barcodes: list[str]) -> ad.AnnData:
    """Load dgCMatrix RDS via R and subset to keep_barcodes; return AnnData."""
    path = COUNTS if COUNTS.exists() else None
    if path is None:
        raise SystemExit(f"Missing counts RDS at {COUNTS}")
    rp = path.resolve().as_posix()
    log(f"Reading counts RDS via R: {path.name}")
    ro.r("library(Matrix)")
    ro.r(f'x <- readRDS("{rp}")')
    dims = list(ro.r("dim(x)"))
    log(f"Full matrix dim: {dims}")

    # Pass keep barcodes to R
    keep = ro.StrVector(keep_barcodes)
    ro.globalenv["keep_bc"] = keep
    ro.r(
        """
        cn <- colnames(x)
        # barcodes in meta may be rownames like Sample.cell; try exact then SampleID prefix
        hit <- cn %in% keep_bc
        if (sum(hit) == 0) {
          # try without transformation
          hit <- cn %in% keep_bc
        }
        xs <- x[, hit, drop=FALSE]
        """
    )
    nkeep = int(ro.r("ncol(xs)")[0])
    log(f"Subset columns: {nkeep}")
    if nkeep == 0:
        # Diagnostic: show colnames style
        print(ro.r("head(colnames(x), 10)"))
        print(ro.r("head(keep_bc, 10)"))
        raise SystemExit("No overlapping barcodes between meta subset and counts matrix")

    # Extract sparse components
    ro.r(
        """
        xs <- as(xs, 'dgCMatrix')
        i <- xs@i
        p <- xs@p
        y <- xs@x
        rn <- rownames(xs)
        cn2 <- colnames(xs)
        nr <- nrow(xs)
        nc <- ncol(xs)
        """
    )
    i = np.array(ro.r("as.integer(i)"), dtype=np.int32)
    p = np.array(ro.r("as.integer(p)"), dtype=np.int32)
    y = np.array(ro.r("as.numeric(y)"), dtype=np.float64)
    rn = list(ro.r("rn"))
    cn = list(ro.r("cn2"))
    nr = int(ro.r("nr")[0])
    nc = int(ro.r("nc")[0])
    # R dgCMatrix is gene x cell; AnnData wants cell x gene
    mat = sparse.csc_matrix((y, i, p), shape=(nr, nc))
    X = mat.T.tocsr()
    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index(cn, name="barcode")
    adata.var_names = pd.Index(rn, name="gene")
    adata.var_names_make_unique()
    return adata


def attach_territory(adata: ad.AnnData, sample_table: pd.DataFrame) -> ad.AnnData:
    for k, v in TEMPORAL_TV.items():
        if k in ("roi", "source_note"):
            continue
        adata.obs[k] = v
    adata.obs["dataset"] = "garcia_nature_2022"
    adata.obs["region"] = "temporal_cortex"
    # map SampleID-level covariates if present on obs
    return adata


def compute_tv_live() -> dict:
    """Recompute Temporal Pole TerritoryVoxels for provenance (optional live check)."""
    try:
        import nibabel as nib
        from nilearn import datasets, image
        from yang_territory_probe import resample_mask_to_atlas, roi_mask_from_labels

        nitrc = (
            ROOT
            / "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/Atlas_182_MNI152"
        )
        tv_path = nitrc / "TerritoryVoxels_BMM.nii"
        if not tv_path.exists():
            return TEMPORAL_TV
        cort = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
        art = nib.load(str(nitrc / "ArterialAtlas_level2.nii"))
        mask = roi_mask_from_labels(cort.maps, list(cort.labels), ["Temporal Pole"])
        mask_rs = resample_mask_to_atlas(mask, cort.maps.affine, art)
        tv = np.rint(np.asanyarray(nib.load(str(tv_path)).dataobj, dtype=float)).astype(np.int32)
        if tv.ndim == 4:
            tv = tv[..., 0]
        vals = tv[mask_rs]
        assigned = vals[(vals >= 1) & (vals <= 4)]
        labels = {1: "ACA", 2: "MCA", 3: "PCA", 4: "VB"}
        probs = {n: float((assigned == lab).sum() / len(assigned)) if len(assigned) else 0.0 for lab, n in labels.items()}
        out = dict(TEMPORAL_TV)
        out.update(
            {
                "tv_p_ACA": probs["ACA"],
                "tv_p_MCA": probs["MCA"],
                "tv_p_PCA": probs["PCA"],
                "tv_p_VB": probs["VB"],
                "tv_dominant": max(probs, key=probs.get) if assigned.size else None,
                "tv_p_dominant": max(probs.values()) if assigned.size else float("nan"),
                "tv_assignment_fraction": float(assigned.size / mask_rs.sum()) if mask_rs.sum() else float("nan"),
                "live_recomputed": True,
            }
        )
        return out
    except Exception as e:  # noqa: BLE001
        log(f"Live TV recompute skipped: {e}")
        return TEMPORAL_TV


def main():
    global TEMPORAL_TV
    TEMPORAL_TV = compute_tv_live()
    TV_LOOKUP.parent.mkdir(parents=True, exist_ok=True)
    TV_LOOKUP.write_text(json.dumps(TEMPORAL_TV, indent=2), encoding="utf-8")

    exvivo = load_exvivo_table()
    log(f"Ex vivo temporal samples in MOESM2: {len(exvivo)}")
    meta = read_meta()
    # Restrict to ex vivo temporal SampleIDs
    keep_ids = set(exvivo["SampleID"])
    meta_sub = meta[meta["SampleID"].isin(keep_ids)].copy()
    # Prefer BVP (brain vascular prep) over whole-tissue libraries when labeled
    if "snRNAPreparation" in meta_sub.columns:
        bvp = meta_sub["snRNAPreparation"].astype(str).str.upper().isin(["BVP", "BVE"])
        if bvp.any():
            log(f"Filtering to BVP/BVE prep: {int(bvp.sum())}/{len(meta_sub)} barcodes")
            meta_sub = meta_sub.loc[bvp].copy()
    log(f"Barcodes in meta for temporal SampleIDs: {len(meta_sub)}")
    if meta_sub.empty:
        print("MOESM2 IDs:", sorted(keep_ids)[:20])
        print("Meta SampleIDs:", sorted(meta["SampleID"].unique())[:30])
        raise SystemExit("No meta barcodes matched ex vivo temporal SampleIDs")

    meta_sub = meta_sub.merge(exvivo, on="SampleID", how="left")
    barcodes = meta_sub["barcode"].astype(str).tolist()

    adata = read_counts_subset(barcodes)
    # Align obs
    meta_sub = meta_sub.set_index("barcode")
    common = adata.obs_names.intersection(meta_sub.index)
    adata = adata[common].copy()
    for col in meta_sub.columns:
        adata.obs[col] = meta_sub.loc[adata.obs_names, col].values

    adata = attach_territory(adata, exvivo)
    # QC
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    import scanpy as sc

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata.uns["territory_annotation"] = TEMPORAL_TV
    adata.uns["load_timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    log(f"Writing {OUT_H5AD} ({adata.n_obs} cells × {adata.n_vars} genes)")
    adata.write_h5ad(OUT_H5AD, compression="gzip")

    # Sample summary
    samp = (
        adata.obs.groupby("SampleID", observed=True)
        .agg(
            n_cells=("SampleID", "size"),
            PatientID=("PatientID", "first"),
            tissue_location=("tissue_location", "first"),
            snRNAPreparation=("snRNAPreparation", "first"),
            median_genes=("n_genes_by_counts", "median"),
            median_counts=("total_counts", "median"),
        )
        .reset_index()
    )
    for k, v in TEMPORAL_TV.items():
        if k.startswith("tv_") or k in ("analysis_stream", "territory_contrast", "method"):
            samp[k] = v
    samp.to_csv(OUT_SAMPLES, index=False)

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Garcia temporal MCA (Class B) load",
        "",
        f"- Cells: **{adata.n_obs:,}**",
        f"- Genes: **{adata.n_vars:,}**",
        f"- Samples: **{samp['SampleID'].nunique()}**",
        f"- Object: `{OUT_H5AD.relative_to(ROOT).as_posix()}`",
        f"- Territory: **{TEMPORAL_TV['tv_dominant']}** "
        f"(P={TEMPORAL_TV['tv_p_dominant']:.3f}, assign frac={TEMPORAL_TV['tv_assignment_fraction']:.3f})",
        "",
        "## Samples",
        "",
        samp.to_string(index=False),
        "",
        "Class B lock: temporal tip/lobe → MCA via TerritoryVoxels; pairs with Yang PCA hippocampus "
        "and/or SFG borderzone depending on contrast.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log("Done")
    print(samp.to_string(index=False))


if __name__ == "__main__":
    main()
