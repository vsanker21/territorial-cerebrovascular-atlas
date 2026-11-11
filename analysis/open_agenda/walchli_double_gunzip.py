"""Stream double-gunzip Wälchli GEO Seurat .rds.gz (still gzip after one pass)."""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

ROOT = Path(r"F:\Cerebrovascular SCS")
SUPPL = ROOT / "data/raw/walchli/geo_GSE256493/suppl"
OUTDIR = SUPPL / "_decomp"
OUTDIR.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "adult_temporal_sorted": "GSE256493_Adult_control_brain_temporal_lobe_sorted_endothelial_cells_seurat_object.rds.gz",
    "avm_sorted": "GSE256493_Arteriovenous_malformations_sorted_endothelial_cells_seurat_object.rds.gz",
    "lgg_sorted": "GSE256493_Lower_grade_glioma_sorted_endothelial_cells_seurat_object.rds.gz",
}


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def double_gunzip(src: Path, dest: Path):
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log(f"SKIP {dest.name} exists ({dest.stat().st_size/1e9:.2f} GB)")
        return dest
    mid = dest.with_suffix(".pass1.gz")
    log(f"pass1 {src.name}")
    t0 = time.time()
    with gzip.open(src, "rb") as fi, mid.open("wb") as fo:
        while True:
            chunk = fi.read(16 * 1024 * 1024)
            if not chunk:
                break
            fo.write(chunk)
    magic = mid.read_bytes()[:2]
    log(f"  pass1 done {mid.stat().st_size/1e9:.2f} GB magic={magic.hex()} in {time.time()-t0:.0f}s")
    if magic != b"\x1f\x8b":
        # already a plain rds after one pass
        mid.replace(dest)
        log(f"  single-gzip only → {dest.name}")
        return dest
    log(f"pass2 {mid.name}")
    t1 = time.time()
    with gzip.open(mid, "rb") as fi, dest.open("wb") as fo:
        n = 0
        while True:
            chunk = fi.read(16 * 1024 * 1024)
            if not chunk:
                break
            fo.write(chunk)
            n += len(chunk)
            if n % (512 * 1024 * 1024) < 16 * 1024 * 1024:
                log(f"  pass2 {n/1e9:.2f} GB")
    mid.unlink(missing_ok=True)
    log(f"OK {dest.name} {dest.stat().st_size/1e9:.2f} GB in {time.time()-t1:.0f}s")
    return dest


def main():
    report = {}
    for key, bn in TARGETS.items():
        src = SUPPL / bn
        if not src.exists():
            report[key] = {"status": "missing"}
            continue
        dest = OUTDIR / f"{key}.rds"
        try:
            double_gunzip(src, dest)
            report[key] = {"status": "ok", "rds": str(dest), "bytes": dest.stat().st_size}
        except Exception as e:
            report[key] = {"status": "error", "error": str(e)}
            log(f"ERROR {key}: {e}")
    out = ROOT / "results/h4/walchli_double_gunzip_report.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    log(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
