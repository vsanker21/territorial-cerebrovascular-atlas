"""Probe Garcia RDS metadata / matrix shapes."""
from __future__ import annotations

from pathlib import Path

import pyreadr
import rpy2.robjects as ro

ROOT = Path(__file__).resolve().parents[2]


def main():
    for p in [
        ROOT / "data/raw/garcia/mit/brain.BBB.humanHD.meta.GEOsubmission.rds",
        ROOT / "data/raw/garcia/mit/ROSMAP.VascularCells.meta.rds",
    ]:
        print("===", p.name)
        try:
            v = pyreadr.read_r(str(p))[None]
            print(v.shape, list(v.columns))
            print(v.head(2).to_string()[:900])
            for c in v.columns:
                if any(
                    x in c.lower()
                    for x in ["region", "tissue", "loc", "cell", "cluster", "class", "zone", "subtype"]
                ):
                    print(c, ":", v[c].astype(str).value_counts().head(12).to_dict())
        except Exception as e:  # noqa: BLE001
            print("fail", e)

    p = ROOT / "data/raw/garcia/mit/brain.BBB.human.counts.GEOsubmission.rds"
    print("reading dims via R...", p)
    ro.r("library(Matrix)")
    # forward slashes for R on Windows
    rp = str(p).replace("\\", "/")
    ro.r(f'x <- readRDS("{rp}")')
    print(ro.r("class(x)"))
    print(ro.r("dim(x)"))
    print(ro.r("head(colnames(x), 5)"))
    print(ro.r("head(rownames(x), 5)"))


if __name__ == "__main__":
    main()
