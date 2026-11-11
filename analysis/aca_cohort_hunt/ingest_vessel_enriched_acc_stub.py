"""Stub pipeline for vessel-enriched ACC / medial ACA vascular object.

Run only after a lab/RADC share lands under data/raw/aca_vessel_enriched/.
Does not invent public substitutes (GSE222007, Garcia temporal, Wälchli TL).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/aca_vessel_enriched"
OUT = ROOT / "data/processed/aca_cohort_hunt"
RESULTS = ROOT / "results/aca_cohort_hunt"
MANIFEST = RAW / "provenance.json"

REQUIRED_KEYS = [
    "source_lab",
    "dua_or_accession",
    "anatomy_label",  # must be ACC / BA24/32/33 / medial ACA-dominant
    "enrichment_method",  # BVE / microvessel / FACS CD31 etc.
    "n_donors_expected",
    "file_path",
]


def main():
    if not MANIFEST.exists():
        print(
            "NO-GO: missing data/raw/aca_vessel_enriched/provenance.json\n"
            "Create it when a vessel-enriched ACC share arrives. "
            "Do not point this at MCA temporal or lateral prefrontal GEO.",
            file=sys.stderr,
        )
        return 1
    meta = json.loads(MANIFEST.read_text(encoding="utf-8"))
    missing = [k for k in REQUIRED_KEYS if k not in meta]
    if missing:
        print(f"NO-GO: provenance missing keys {missing}", file=sys.stderr)
        return 1
    anatomy = str(meta["anatomy_label"]).lower()
    if any(x in anatomy for x in ("temporal", "mtg", "stg", "mca")) and "acc" not in anatomy:
        print("NO-GO: anatomy looks MCA/temporal — not an ACA arm", file=sys.stderr)
        return 1
    if not any(x in anatomy for x in ("acc", "cingulate", "ba24", "ba32", "medial")):
        print("NO-GO: anatomy_label must clearly indicate ACC/medial ACA", file=sys.stderr)
        return 1

    RESULTS.mkdir(parents=True, exist_ok=True)
    plan = {
        "status": "READY_TO_INGEST",
        "next_steps": [
            "Load counts + published vascular labels (prefer published over marker seed)",
            "TerritoryVoxels HO prior; require tv_dominant==ACA or high tv_p_ACA",
            "Endo gate; donor OOF vs Garcia MCA tight endo",
            "Enrichment-matched pseudobulk DE; report ambient scores",
            "Update docs/10_class_b_claims_paragraph.md only if DE survives ambient",
        ],
        "provenance": meta,
        "outputs_planned": {
            "h5ad": str((OUT / "vessel_enriched_acc_endo.h5ad").relative_to(ROOT)),
            "class_b_refresh": "data/processed/class_b/aca_mca_harmony_endo_enrichment_matched.h5ad",
            "report": "results/class_b/enrichment_matched_aca_mca_report.md",
        },
    }
    outp = RESULTS / "vessel_enriched_acc_ingest_plan.json"
    outp.write_text(json.dumps(plan, indent=2), encoding="utf-8")
    print(json.dumps(plan, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
