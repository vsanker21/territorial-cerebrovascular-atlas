#!/usr/bin/env python3
"""OpenNeuro ds005096 TOF-MRA methods pilot — subject-level vessel metrics.

Scientific role: calibrate subject-level angio feature extraction.
NOT confirmatory Track C (no territorial stroke/WMH lesions in this aneurysm cohort).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/raw/openneuro_tof_mra/ds005096"
OUT = ROOT / "results/nature_gates"
OUT.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def mask_stats(path: Path) -> dict:
    try:
        import nibabel as nib
    except ImportError:
        return {"path": str(path), "status": "nibabel_missing"}
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    vox = img.header.get_zooms()[:3]
    vox_vol = float(np.prod(vox))
    pos = data > 0
    n = int(pos.sum())
    return {
        "path": str(path.relative_to(DATA)),
        "status": "ok",
        "shape": list(data.shape),
        "voxel_mm": [float(x) for x in vox],
        "n_foreground_voxels": n,
        "volume_mm3": float(n * vox_vol),
        "intensity_mean_fg": float(data[pos].mean()) if n else None,
    }


def main() -> None:
    parts = DATA / "participants.tsv"
    participants = pd.read_csv(parts, sep="\t") if parts.exists() else pd.DataFrame()

    rows = []
    for sub in ["sub-000", "sub-001", "sub-002"]:
        # Prefer parent-artery nifti segmentations
        hits = list((DATA / "derivatives" / sub).rglob("*Parent*Artery*.nii*"))
        hits = [h for h in hits if h.suffix.lower() in {".nii", ".gz"} or h.name.endswith(".nii")]
        # also incomplete downloads ending oddly
        hits += [
            h
            for h in (DATA / "derivatives" / sub).rglob("*.nii*")
            if "Parent" in h.name or "parent" in h.name
        ]
        hits = sorted(set(hits))
        if not hits:
            rows.append({"subject": sub, "status": "no_parent_artery_nifti"})
            continue
        # pick first complete-looking file (skip partial .C7Be... temps if tiny sibling exists)
        chosen = None
        for h in hits:
            if h.stat().st_size > 1_000_000 and not h.name.endswith((".C7Be527D",)):
                chosen = h
                break
        if chosen is None:
            chosen = max(hits, key=lambda p: p.stat().st_size)
        st = mask_stats(chosen)
        st["subject"] = sub
        rows.append(st)

    status = {
        "date": time.strftime("%Y-%m-%d"),
        "dataset": "OpenNeuro ds005096 (RBWH TOF-MRA intracranial aneurysm)",
        "n_participants_tsv": int(len(participants)),
        "pilot_subjects": ["sub-000", "sub-001", "sub-002"],
        "role": "METHODS_CALIBRATION_ONLY",
        "not_confirmatory_track_c": (
            "No territorial stroke/WMH lesion maps; cannot test burden ~ angio + territory."
        ),
        "subject_metrics": rows,
        "verdict": "OPENNEURO_TOF_PILOT_PARTIAL__METHODS_ONLY",
        "next_confirmatory": "ISLES'24 CTA + CoW/LVO masks + follow-up lesion masks (subject-level angio + lesions)",
    }
    (OUT / "openneuro_tof_mra_pilot.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(OUT / "openneuro_tof_mra_pilot_metrics.csv", index=False)

    md = [
        "# OpenNeuro TOF-MRA methods pilot (ds005096)",
        "",
        f"**Verdict:** `{status['verdict']}`",
        "",
        "## Role",
        status["not_confirmatory_track_c"],
        "",
        f"Participants.tsv n={status['n_participants_tsv']}. Pilot subjects: {status['pilot_subjects']}.",
        "",
        "## Metrics",
        "See `openneuro_tof_mra_pilot_metrics.csv`.",
        "",
        "## Confirmatory Track C path",
        status["next_confirmatory"],
        "",
    ]
    (OUT / "openneuro_tof_mra_pilot.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
