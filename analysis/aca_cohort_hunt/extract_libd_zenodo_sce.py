"""Extract sce_azimuth / vascular labels from completed Zenodo spatialdACC zip."""
from __future__ import annotations

import json
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ZIP_CANDIDATES = [
    ROOT / "data/raw/aca_cohort_hunt/zenodo_spatialdACC/spatialdACC.zip",
    ROOT / "data/raw/aca_cohort_hunt/zenodo_spatialdACC/LieberInstitute_spatialdACC-v1.0.0-preprint.zip",
]
PART = ROOT / "data/raw/aca_cohort_hunt/zenodo_spatialdACC/spatialdACC.zip.part"
OUTDIR = ROOT / "data/raw/aca_cohort_hunt/libd_annotations/zenodo_extract"
TARGETS = (
    "sce_azimuth.Rdata",
    "sce_Final.Rdata",
    "sce_harmony.Rdata",
)


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def find_zip() -> Path:
    for z in ZIP_CANDIDATES:
        if z.exists() and z.stat().st_size > 1_000_000_000:
            return z
    if PART.exists() and PART.stat().st_size >= 13_510_296_335 - 1024:
        # rename completed part
        dest = ZIP_CANDIDATES[0]
        PART.rename(dest)
        return dest
    raise FileNotFoundError("Zenodo zip not complete yet")


def main():
    zpath = find_zip()
    log(f"Opening {zpath} ({zpath.stat().st_size/1e9:.2f} GB)")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    extracted = []
    with zipfile.ZipFile(zpath) as zf:
        names = zf.namelist()
        log(f"n entries: {len(names)}")
        hits = [n for n in names if any(t in n.replace("\\", "/") for t in TARGETS)]
        log(f"target hits: {len(hits)}")
        for n in hits:
            log(f"extract {n}")
            zf.extract(n, OUTDIR)
            extracted.append(n)
        # also list any snRNA-seq/05_azimuth paths
        az = [n for n in names if "05_azimuth" in n and n.lower().endswith((".rdata", ".rds", ".csv"))]
        (OUTDIR / "azimuth_member_list.json").write_text(json.dumps(az, indent=2), encoding="utf-8")
    summary = {"zip": str(zpath), "extracted": extracted}
    (OUTDIR / "extract_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
