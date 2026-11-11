"""
Metadata audit: can Winkler / Garcia support Class B (ACA vs MCA isocortex)?

Garcia Nature 2022 Supplementary Table 1 is ingested when present.
Winkler Science 2022 regional detail is literature-curated (normal tissue =
epilepsy lobectomy cortex; lobe often temporal — verify against table S1 / dbGaP).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GARCIA_XLSX = ROOT / "data" / "raw" / "garcia" / "Garcia_SuppTable1.xlsx"
OUT_DIR = ROOT / "results" / "metadata_audit"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# A priori territorial priors for Garcia ROSMAP regions (hard-atlas / clinical anatomy)
GARCIA_REGION_PRIOR = {
    "prefrontal cortex": {
        "ho_roi": "middle_frontal_gyrus / frontal_medial_cortex (split)",
        "a_priori_territory": "mixed ACA/MCA — requires subregional split for Class B",
        "class_b_candidate": "conditional",
        "notes": "PFC spans ACA medial and MCA lateral; do not assign single territory without finer label.",
    },
    "mid-temporal cortex": {
        "ho_roi": "Middle Temporal Gyrus",
        "a_priori_territory": "MCA",
        "class_b_candidate": "yes_mca",
        "notes": "Strong MCA prior; excellent Class B partner vs medial frontal ACA.",
    },
    "angular gyrus": {
        "ho_roi": "Angular Gyrus",
        "a_priori_territory": "MCA (inferior parietal)",
        "class_b_candidate": "yes_mca",
        "notes": "Classic MCA parietal territory.",
    },
    "entorhinal cortex": {
        "ho_roi": "not in HO cortical maxprob as named; near temporal/para-hippocampal",
        "a_priori_territory": "PCA / posterior circulation influence",
        "class_b_candidate": "no_use_for_class_a_pca",
        "notes": "Better as PCA-leaning cortical ROI than MCA Class B.",
    },
    "thalamus": {
        "ho_roi": "Thalamus",
        "a_priori_territory": "mixed perforators (PCA/ACA/VB-associated)",
        "class_b_candidate": "no",
        "notes": "Perforator territories; high entropy expected — Aim 3 / special case.",
    },
    "hippocampus": {
        "ho_roi": "Hippocampus",
        "a_priori_territory": "PCA",
        "class_b_candidate": "no_use_for_class_a",
        "notes": "Pairs with Yang hippocampus for cross-dataset PCA replication.",
    },
    "mammillary body": {
        "ho_roi": "none standard",
        "a_priori_territory": "posterior perforator / uncertain",
        "class_b_candidate": "no",
        "notes": "Too small / atypical for primary territorial tests.",
    },
}

WINKLER_CURATED = {
    "dataset": "Winkler EA et al. Science 2022;375:eabi7377",
    "normal_tissue": {
        "description": "Normal cerebral cortex from tailored lobectomies for epilepsy (n donors with scRNA ≈ 5 for normal arm contributing 74,535 cells before AVM joint atlas).",
        "region_granularity": "Reported as cerebral cortex; epilepsy lobectomies are frequently temporal — exact lobe per donor requires table S1 / dbGaP phs002624.",
        "territorial_implication": (
            "If predominantly temporal neocortex → MCA-dominant prior. "
            "If frontal/parietal mix without labels → cannot support clean Class B."
        ),
        "class_b_status": "PENDING_TABLE_S1",
        "access": "dbGaP phs002624.v2.p1 (controlled); UCSC cell browser adult-brain-vasc.cells.ucsc.edu; Zenodo cited in Science for objects",
    },
    "avm_tissue": {
        "description": "Intraoperative angiographically confirmed AVMs (5 donors; 106,853 cells reported for AVM arm).",
        "region_granularity": "Lesion-centric (nidus), not cytoarchitectonic region — territorial indexing may use AVM arterial feeder territory if recorded clinically.",
        "class_b_status": "not_applicable_use_for_H4",
        "role": "H4 disease de-differentiation / zonation loss; not naive territorial identity discovery.",
    },
    "recommendation": (
        "Do not assign Winkler normal cells to MCA until lobe-level metadata are verified. "
        "Garcia mid-temporal + angular gyrus are the current best open Class B MCA candidates."
    ),
}


def normalize_region(x: str) -> str:
    if not isinstance(x, str):
        return ""
    x = x.strip().lower()
    x = x.replace("_", " ").replace("-", " ")
    x = re.sub(r"\s+", " ", x)
    aliases = {
        "pfc": "prefrontal cortex",
        "dlpfc": "prefrontal cortex",
        "prefrontal": "prefrontal cortex",
        "mtc": "mid-temporal cortex",
        "middle temporal cortex": "mid-temporal cortex",
        "mid temporal cortex": "mid-temporal cortex",
        "ag": "angular gyrus",
        "erc": "entorhinal cortex",
        "entorhinal": "entorhinal cortex",
        "hippocampal": "hippocampus",
        "mb": "mammillary body",
        "mammillary": "mammillary body",
    }
    return aliases.get(x, x)


def audit_garcia() -> dict:
    result = {
        "file": str(GARCIA_XLSX),
        "present": GARCIA_XLSX.exists(),
        "sheets": [],
        "region_counts": [],
        "class_b_assessment": {},
    }
    if not GARCIA_XLSX.exists():
        result["error"] = "Garcia supplementary XLSX missing"
        return result

    xl = pd.ExcelFile(GARCIA_XLSX)
    result["sheets"] = xl.sheet_names
    frames = []
    for sheet in xl.sheet_names:
        df = pd.read_excel(GARCIA_XLSX, sheet_name=sheet)
        df.columns = [str(c).strip() for c in df.columns]
        df["_sheet"] = sheet
        frames.append(df)
        # save each sheet
        df.to_csv(OUT_DIR / f"garcia_{re.sub(r'[^A-Za-z0-9]+', '_', sheet)}.csv", index=False)

    big = pd.concat(frames, ignore_index=True)
    # Heuristic: find a region-like column
    region_col = None
    for c in big.columns:
        cl = c.lower()
        if any(k in cl for k in ("region", "brain", "area", "structure", "tissue")):
            region_col = c
            break
    result["detected_region_column"] = region_col
    if region_col:
        vc = (
            big[region_col]
            .dropna()
            .astype(str)
            .map(normalize_region)
            .value_counts()
            .rename_axis("region")
            .reset_index(name="n_rows")
        )
        result["region_counts"] = vc.to_dict(orient="records")
        vc.to_csv(OUT_DIR / "garcia_region_counts.csv", index=False)

        assessments = []
        for region in vc["region"]:
            prior = None
            for key, val in GARCIA_REGION_PRIOR.items():
                if key in region or region in key:
                    prior = {"region_matched": key, **val}
                    break
            assessments.append({"region_in_table": region, "prior": prior})
        result["class_b_assessment"] = assessments

    # Explicit priors table for paper methods
    prior_df = pd.DataFrame(
        [{"region": k, **v} for k, v in GARCIA_REGION_PRIOR.items()]
    )
    prior_df.to_csv(OUT_DIR / "garcia_region_territory_priors.csv", index=False)

    mca_yes = [k for k, v in GARCIA_REGION_PRIOR.items() if v["class_b_candidate"] == "yes_mca"]
    result["class_b_mca_regions"] = mca_yes
    result["verdict"] = (
        f"Garcia multi-region design can support Class B if mid-temporal and/or angular gyrus "
        f"samples are present with usable cell counts. Prioritized MCA ROIs: {mca_yes}. "
        f"Prefrontal cortex must be split or treated as mixed."
    )
    return result


def main():
    garcia = audit_garcia()
    payload = {
        "garcia": garcia,
        "winkler": WINKLER_CURATED,
        "overall_recommendation": {
            "class_b_primary_source": "Garcia ROSMAP mid-temporal cortex and angular gyrus",
            "class_b_aca_source": "Yang superior frontal (empirically ACA) and/or Garcia medial-leaning PFC if split; HO frontal medial cortex as anatomical target",
            "class_a_pca": "Yang hippocampus + Garcia hippocampus + Garcia entorhinal (secondary)",
            "winkler_role": "H4 AVM + zonation reference; Class B only after lobe metadata verification",
            "blocking_items": [
                "NITRC login for BMM/borderzone maps",
                "Winkler table S1 / dbGaP lobe labels",
                "Confirm Garcia sheet region field parsing against actual column names",
            ],
        },
    }
    out = OUT_DIR / "garcia_winkler_metadata_audit.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# Metadata audit — Garcia & Winkler (Class B readiness)",
        "",
        "## Verdict",
        "",
        payload["overall_recommendation"]["class_b_primary_source"] + ".",
        "",
        f"- ACA arm: {payload['overall_recommendation']['class_b_aca_source']}",
        f"- PCA Class A: {payload['overall_recommendation']['class_a_pca']}",
        f"- Winkler: {payload['overall_recommendation']['winkler_role']}",
        "",
        "## Garcia",
        "",
        garcia.get("verdict", ""),
        "",
        f"File present: {garcia.get('present')}; sheets: {garcia.get('sheets')}",
        f"Detected region column: {garcia.get('detected_region_column')}",
        f"MCA Class B candidate regions: {garcia.get('class_b_mca_regions')}",
        "",
        "## Winkler (curated from publication text)",
        "",
        WINKLER_CURATED["recommendation"],
        "",
        "Normal tissue class_b_status:",
        WINKLER_CURATED["normal_tissue"]["class_b_status"],
        "",
        "## Blocking items",
        "",
    ]
    for b in payload["overall_recommendation"]["blocking_items"]:
        md.append(f"- {b}")
    (OUT_DIR / "garcia_winkler_metadata_audit.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(payload["overall_recommendation"], indent=2))
    print("Garcia region column:", garcia.get("detected_region_column"))
    print("Wrote", out)


if __name__ == "__main__":
    main()
