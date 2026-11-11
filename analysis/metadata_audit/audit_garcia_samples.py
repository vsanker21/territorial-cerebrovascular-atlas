"""Audit Garcia Nature 2022 MOESM2 sample tables for Class B territorial readiness."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "feasibility"))
from yang_territory_probe import (  # noqa: E402
    LEVEL2_MAJOR,
    load_arterial_atlas,
    load_harvard_oxford,
    resample_mask_to_atlas,
    roi_mask_from_labels,
    score_roi,
)

MOESM2 = ROOT / "data" / "raw" / "garcia" / "esm" / "41586_2022_4521_MOESM2_ESM.xlsx"
ATLAS = ROOT / "data" / "raw" / "arterial_atlas" / "Atlas_182"
OUT = ROOT / "results" / "metadata_audit"
OUT.mkdir(parents=True, exist_ok=True)

# Map Garcia labels → HO ROIs for territorial scoring
REGION_MAP = {
    "Prefrontal_cortex": {
        "atlas": "cortical",
        "labels": ["Middle Frontal Gyrus"],
        "role": "mixed_or_MCA_leaning_proxy — PFC not medial-split; MFG used as lateral proxy",
        "class_b": "conditional",
    },
    "Mid_temporal_cortex": {
        "atlas": "cortical",
        "labels": ["Middle Temporal Gyrus"],
        "role": "MCA Class B primary",
        "class_b": "yes_mca",
    },
    "Angular_gyrus": {
        "atlas": "cortical",
        "labels": ["Angular Gyrus"],
        "role": "MCA Class B primary",
        "class_b": "yes_mca",
    },
    "Entorhinal_cortex": {
        "atlas": "cortical",
        "labels": ["Parahippocampal Gyrus"],
        "role": "PCA-leaning Class A partner",
        "class_b": "no",
    },
    "Hippocampus": {
        "atlas": "subcortical",
        "labels": ["Left Hippocampus", "Right Hippocampus"],
        "role": "PCA Class A / Yang replication",
        "class_b": "no",
    },
    "Thalamus": {
        "atlas": "subcortical",
        "labels": ["Left Thalamus", "Right Thalamus"],
        "role": "perforator / high entropy",
        "class_b": "no",
    },
}

# Ex vivo temporal tip ≈ temporal pole / middle temporal
EXVIVO_TEMPORAL = {
    "atlas": "cortical",
    "labels": ["Temporal Pole"],
    "role": "Garcia ex vivo epilepsy cortex — MCA prior",
    "class_b": "yes_mca",
}


def main():
    exvivo = pd.read_excel(MOESM2, sheet_name="Ex Vivo Sample ID")
    hd = pd.read_excel(MOESM2, sheet_name="Control and HD Samples")
    rosmap = pd.read_excel(MOESM2, sheet_name="ROSMAP")

    exvivo.to_csv(OUT / "garcia_exvivo_samples.csv", index=False)
    hd.to_csv(OUT / "garcia_hd_control_samples.csv", index=False)
    rosmap.to_csv(OUT / "garcia_rosmap_samples.csv", index=False)

    rosmap_counts = (
        rosmap["brain_region"].astype(str).value_counts().rename_axis("brain_region").reset_index(name="n_libraries")
    )
    rosmap_counts.to_csv(OUT / "garcia_rosmap_region_counts.csv", index=False)

    hd_counts = (
        hd["Brain Region"].astype(str).value_counts().rename_axis("brain_region").reset_index(name="n_samples")
    )
    hd_counts.to_csv(OUT / "garcia_hd_region_counts.csv", index=False)

    art_l1 = load_arterial_atlas(ATLAS / "ArterialAtlas.nii")
    art_l2 = load_arterial_atlas(ATLAS / "ArterialAtlas_level2.nii")
    atlases = load_harvard_oxford()

    rows = []
    # Score ROSMAP regions present
    for region, n in zip(rosmap_counts["brain_region"], rosmap_counts["n_libraries"]):
        key = region
        spec = REGION_MAP.get(key)
        if not spec:
            # try fuzzy
            for k, v in REGION_MAP.items():
                if k.lower().replace("_", "") in region.lower().replace("_", "").replace(" ", ""):
                    spec = v
                    key = k
                    break
        if not spec:
            rows.append({"source": "ROSMAP", "label": region, "n": int(n), "status": "unmapped_roi"})
            continue
        scored = score_roi(key, spec, atlases, art_l1, art_l2)
        scored.update({"source": "ROSMAP", "label": region, "n_libraries": int(n), "class_b": spec["class_b"], "role": spec["role"]})
        rows.append(scored)

    # Ex vivo temporal
    scored = score_roi("exvivo_temporal", EXVIVO_TEMPORAL, atlases, art_l1, art_l2)
    scored.update(
        {
            "source": "ex_vivo_BVE",
            "label": "temporal tip/lobe (from Tissue location)",
            "n_libraries": int(exvivo.shape[0]),
            "class_b": "yes_mca",
            "role": EXVIVO_TEMPORAL["role"],
            "tissue_locations": exvivo["Tissue location"].astype(str).value_counts().to_dict(),
        }
    )
    rows.append(scored)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "garcia_region_territory_scores.csv", index=False)

    mca_ok = df[(df.get("class_b") == "yes_mca") & (df.get("dominant_l2") == "MCA") & (df.get("high_specificity_pdom_ge_0.50") == True)]
    verdict = {
        "sample_table": str(MOESM2.name),
        "exvivo_n": int(exvivo.shape[0]),
        "exvivo_locations": exvivo["Tissue location"].astype(str).value_counts().to_dict(),
        "rosmap_regions": rosmap_counts.to_dict(orient="records"),
        "hd_regions": hd_counts.to_dict(orient="records"),
        "class_b_mca_ready_rows": mca_ok[["source", "label", "dominant_l2", "p_dominant_l2", "entropy_l2_bits"]].to_dict(orient="records")
        if len(mca_ok)
        else [],
        "verdict": (
            "Garcia supports Class B MCA arm via (i) ex vivo temporal epilepsy cortex and "
            "(ii) ROSMAP mid-temporal / angular gyrus if present in region counts. "
            "Pair with Yang superior frontal (ACA) for territory-differing isocortex contrast. "
            "Winkler normal cortex remains pending lobe labels but is likely temporal-enriched by clinical context."
        ),
    }
    (OUT / "garcia_class_b_verdict.json").write_text(json.dumps(verdict, indent=2), encoding="utf-8")

    md = [
        "# Garcia sample audit (MOESM2) — Class B readiness",
        "",
        "## Verdict",
        "",
        verdict["verdict"],
        "",
        "## Ex vivo tissue locations",
        "",
        "```json",
        json.dumps(verdict["exvivo_locations"], indent=2),
        "```",
        "",
        "## ROSMAP brain_region counts",
        "",
        rosmap_counts.to_markdown(index=False),
        "",
        "## Territory scores",
        "",
        df.to_markdown(index=False),
        "",
    ]
    (OUT / "garcia_class_b_verdict.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
