"""Verify Yang + Garcia Class B objects and write combined contrast summary."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
YANG = ROOT / "data/processed/yang/yang_vascular.h5ad"
GARCIA = ROOT / "data/processed/garcia/garcia_temporal_mca.h5ad"
OUT = ROOT / "results/class_b_contrast_ready.md"


def main():
    y = ad.read_h5ad(YANG, backed="r")
    g = ad.read_h5ad(GARCIA, backed="r")

    y_stream = y.obs.groupby("analysis_stream", observed=True).size()
    y_contrast = y.obs.groupby("territory_contrast", observed=True).size()
    g_n = g.n_obs
    g_samples = g.obs["SampleID"].nunique()

    lines = [
        "# Class B / territorial contrast — objects ready",
        "",
        "## Yang (`yang_vascular.h5ad`)",
        f"- Cells: **{y.n_obs:,}** (raw filtered 10x barcodes; published vascular n is post-QC)",
        f"- Genes: **{y.n_vars:,}**",
        "- TerritoryVoxels streams attached per library",
        "",
        "### analysis_stream",
        y_stream.to_string(),
        "",
        "### territory_contrast",
        y_contrast.to_string(),
        "",
        "## Garcia temporal MCA (`garcia_temporal_mca.h5ad`)",
        f"- Cells: **{g_n:,}** (BVP prep, temporal tissue only)",
        f"- Samples: **{g_samples}**",
        f"- Genes: **{g.n_vars:,}**",
        f"- tv_dominant: **{g.obs['tv_dominant'].iloc[0]}** (P={g.obs['tv_p_MCA'].iloc[0]})",
        f"- analysis_stream: **{g.obs['analysis_stream'].iloc[0]}**",
        "",
        "## Locked contrasts for next DE / variance analyses",
        "",
        "| Contrast | Arm A | Arm B | Hypothesis |",
        "|---|---|---|---|",
        "| Within-Yang Aim3 | PCA_hippocampus | ACA_MCA_borderzone_SFG | H2 / Aim3 borderzone |",
        "| Class B MCA vs Yang PCA | Garcia MCA_temporal | Yang PCA_hippocampus | H1 territory (region-confounded; treat carefully) |",
        "| Class B MCA vs Yang SFG | Garcia MCA_temporal | Yang ACA_MCA_borderzone_SFG | MCA vs watershed |",
        "",
        "Note: Cross-dataset Class B tests require donor/batch harmonization and vascular-cell QC before DE.",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT.read_text(encoding="utf-8"))
    y.file.close()
    g.file.close()


if __name__ == "__main__":
    main()
