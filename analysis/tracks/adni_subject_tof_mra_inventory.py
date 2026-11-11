#!/usr/bin/env python3
"""ADNI subject-level TOF/MRA inventory for Track C mechanism gate.

Scientific question: do on-disk ADNI tables contain subject-level angiography
usable for nested architecture models (burden ~ angio + territory)?

Atlas-constant MRA covariates in Aim 4 are NOT subject-level and remain locked
as collinear with territory labels.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ADNI = ROOT / "data/raw/ADNI"
OUT = ROOT / "results/nature_gates"
OUT.mkdir(parents=True, exist_ok=True)

PAT = re.compile(r"\bMRA\b|TOF|ANGIO|TIME[\s_-]*OF[\s_-]*FLIGHT|ARTERIOGRAM", re.I)


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def main() -> None:
    key = pd.read_csv(ADNI / "All_Subjects_Key_MRI_24Jul2026.csv", low_memory=False)
    series_types = key["series_type"].astype(str).str.upper().value_counts().to_dict()

    text_cols = [c for c in key.columns if key[c].dtype == object]
    hit_mask = pd.Series(False, index=key.index)
    hit_cols: dict[str, int] = {}
    for c in text_cols:
        m = key[c].astype(str).str.contains(PAT)
        if m.any():
            hit_cols[c] = int(m.sum())
            hit_mask |= m

    nest = ROOT / "results/aim4/adni/aim4_nest_wmh_mra_subject_table.csv"
    wmh_n = int(len(pd.read_csv(nest))) if nest.exists() else None

    # Literature / protocol lock: ADNI core MRI does not include TOF-MRA
    status = {
        "date": time.strftime("%Y-%m-%d"),
        "source_key": "All_Subjects_Key_MRI_24Jul2026.csv",
        "n_series_rows": int(len(key)),
        "n_unique_subjects": int(key["subject_id"].nunique()),
        "series_type_counts": series_types,
        "mra_tof_angio_text_hits": hit_cols,
        "n_mra_tof_angio_rows": int(hit_mask.sum()),
        "wmh_nest_n_subjects": wmh_n,
        "wmh_nest_mra_covariates": "atlas-constant only (vessel_probability/radius by territory)",
        "adni_protocol_tof_mra": False,
        "verdict": "ADNI_SUBJECT_TOF_MRA__NO_GO_ON_DISK_AND_PROTOCOL",
        "implication": (
            "Track C subject-level architecture mediation cannot be completed on ADNI WMH. "
            "Keep ADNI for territory-normalized WMH geography only. "
            "Pursue stroke/aneurysm cohorts with subject TOF+lesions, or OpenNeuro TOF methods pilot."
        ),
        "next_sources": [
            "OpenNeuro ds005096 TOF-MRA (methods/calibration; no stroke lesions)",
            "Stroke center cohorts with TOF/CTA + MNI lesions (request/collaborate)",
            "Do not re-interpret atlas MRA null as subject architecture",
        ],
    }

    (OUT / "adni_subject_tof_mra_inventory.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )

    md = [
        "# ADNI subject-level TOF/MRA inventory — Track C gate",
        "",
        f"**Date:** {status['date']}",
        f"**Verdict:** `{status['verdict']}`",
        "",
        "## Finding",
        f"- MRI key rows scanned: **{status['n_series_rows']:,}** "
        f"({status['n_unique_subjects']:,} subjects).",
        f"- Text hits for MRA/TOF/angio/TOF synonyms: **{status['n_mra_tof_angio_rows']}**.",
        f"- Present `series_type` values include T1W, FLAIR, ASL, dMRI, BOLD, etc. — **no MRA/TOF class**.",
        f"- Aim 4 WMH nest n={wmh_n}: MRA columns are **atlas-constant**, not per-subject angio.",
        "",
        "## Protocol context",
        "ADNI core MRI (ADNI1–4 documentation) emphasizes T1, FLAIR, T2*, dMRI, ASL, tf-fMRI — "
        "**not** TOF-MRA. Local CSV search across all ADNI tables found no TOF/MRA modality fields "
        "(only `MRACQUISITIONTYPEFLAIR`, unrelated).",
        "",
        "## Scientific implication (Nature Track C)",
        status["implication"],
        "",
        "## Do not claim",
        "- That atlas-constant MRA null proves or refutes subject-level architecture mediation.",
        "- That ADNI WMH geography alone is a complete Nature mechanism package.",
        "",
        "## Next",
        *[f"- {x}" for x in status["next_sources"]],
        "",
    ]
    (OUT / "adni_subject_tof_mra_inventory.md").write_text("\n".join(md), encoding="utf-8")
    log(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
