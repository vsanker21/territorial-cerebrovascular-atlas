"""Quick structural spot-check of a few wmh_mni.nii.gz maps."""
from __future__ import annotations

import csv
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "results/aim4/adni/adni_wmh_t1_mni_scaled.csv"
PROC = ROOT / "data/processed/adni/t1_mni"


def main() -> None:
    rows = [r for r in csv.DictReader(CSV.open()) if r.get("status") == "ok"]
    by_n = sorted(rows, key=lambda r: int(float(r["n_wmh_voxels"])))
    pick = [rows[0], rows[len(rows) // 2], rows[-1], by_n[0], by_n[-1]]
    seen: set[str] = set()
    print(f"spot-check n_ok={len(rows)}")
    for r in pick:
        ptid = r["ptid"]
        if ptid in seen:
            continue
        seen.add(ptid)
        path = PROC / ptid / "wmh_mni.nii.gz"
        if not path.exists():
            print(f"{ptid}: MISSING {path}")
            continue
        img = nib.load(str(path))
        d = np.asanyarray(img.dataobj) > 0
        nz = int(d.sum())
        com = np.round(np.column_stack(np.where(d)).mean(axis=0), 1) if nz else None
        frac = float(r["wmh_frac_assigned_tv"])
        print(
            f"{ptid}: n_wmh={nz} csv={r['n_wmh_voxels']} "
            f"frac_tv={frac:.3f} com={com} shape={d.shape} "
            f"det={np.linalg.det(img.affine):.3f}"
        )


if __name__ == "__main__":
    main()
