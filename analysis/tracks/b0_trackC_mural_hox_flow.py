"""
B0 gates on available multiregion vascular snRNA-seq (Garcia ROSMAP object =
ROSMAP multi-region vascular nuclei already local; Sun 2023 full release
gated on Dropbox/Synapse if MIT RDS incomplete).

Gates:
  B0i  mural (Per+SMC) nuclei per region per donor — fail if median < 5
  B0ii HOX anterior vs posterior enrichment in mural cells

Also runs Track C shear/flow module scores on endothelium by territory.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tracks" / "b0_trackC"
OUT.mkdir(parents=True, exist_ok=True)
H5AD = ROOT / "data/processed/garcia/garcia_rosmap_vascular.h5ad"

# Territory priors for Garcia ROSMAP regions (Liu TerritoryVoxels locked elsewhere)
REGION_TERRITORY = {
    "Prefrontal_cortex": "ACA_MCA_border_or_MCA",  # PFC often MCA-leaning lateral; medial=ACA — flag ambiguous
    "Angular_gyrus": "MCA",
    "Midtemporal_cortex": "MCA",
    "Entorhinal_cortex": "PCA",
    "Hippocampus": "PCA",
    "Anterior_thalamus": "PCA_VB_perforator",
    "Mammillary_body": "PCA_VB_perforator",
}

ANTERIOR_REGIONS = {"Prefrontal_cortex", "Angular_gyrus", "Midtemporal_cortex"}
POSTERIOR_REGIONS = {"Entorhinal_cortex", "Hippocampus", "Anterior_thalamus", "Mammillary_body"}

# Matched-cytoarchitecture pair from concept note
MATCHED_PAIR = ("Midtemporal_cortex", "Entorhinal_cortex")  # MCA vs PCA

HOX_GENES = [
    "HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA7", "HOXA9", "HOXA10",
    "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9",
    "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXC10",
    "HOXD1", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10",
]
# Shear / flow endothelial programs (Track C)
FLOW_GENES = [
    "KLF2", "KLF4", "NOS3", "THBD", "PLAT", "CYBB",  # flow / atheroprotective
    "ICAM1", "VCAM1", "SELE", "CCL2",  # inflammatory / disturbed flow
    "ASS1", "ARG2", "MEF2C", "FOXP1",
]


def log(msg: str):
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def module_score(X: np.ndarray, gene_idx: list[int]) -> np.ndarray:
    if not gene_idx:
        return np.full(X.shape[0], np.nan)
    return np.asarray(X[:, gene_idx].mean(axis=1)).ravel()


def main():
    log("Loading Garcia ROSMAP vascular h5ad (multiregion; local)")
    a = ad.read_h5ad(H5AD)
    # densify gene subset only
    genes = list(a.var_names)
    gene_to_i = {g: i for i, g in enumerate(genes)}

    obs = a.obs.copy()
    obs["donor"] = obs["subject"].astype(str)
    obs["region"] = obs["brain_region"].astype(str)
    obs["celltype"] = obs["celltype"].astype(str)
    obs["is_mural"] = obs["celltype"].isin(["Per", "SMC"])
    obs["is_endo"] = obs["celltype"] == "Endo"
    obs["pole"] = obs["region"].map(
        lambda r: "anterior" if r in ANTERIOR_REGIONS else ("posterior" if r in POSTERIOR_REGIONS else "other")
    )
    obs["territory_prior"] = obs["region"].map(REGION_TERRITORY)

    # ----- B0i mural counts -----
    # Denominator = donors with ANY vascular cells in that region (not all 214 donors).
    mural = obs[obs["is_mural"]].copy()
    region_donors = (
        obs.groupby(["donor", "region"], observed=True)
        .size()
        .rename("n_vascular_any")
        .reset_index()
    )
    counts = (
        mural.groupby(["donor", "region", "celltype"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    counts_tot = (
        mural.groupby(["donor", "region"], observed=True)
        .size()
        .rename("n_mural")
        .reset_index()
    )
    # Left-join onto region-present donors; missing mural → 0
    base = region_donors.merge(counts_tot, on=["donor", "region"], how="left")
    base["n_mural"] = base["n_mural"].fillna(0).astype(int)
    pivot = base.pivot_table(index="donor", columns="region", values="n_mural", fill_value=np.nan)
    pivot.to_csv(OUT / "b0_mural_counts_donor_by_region.csv")
    counts.to_csv(OUT / "b0_mural_counts_long.csv", index=False)
    base.to_csv(OUT / "b0_mural_counts_region_present_donors.csv", index=False)

    medians = {c: float(np.nanmedian(pivot[c].values)) for c in pivot.columns}
    pair_stats = {}
    for reg in list(MATCHED_PAIR) + ["Prefrontal_cortex", "Hippocampus", "Anterior_thalamus"]:
        if reg not in pivot.columns:
            continue
        s = pivot[reg].dropna()
        pair_stats[reg] = {
            "n_donors_with_region": int(len(s)),
            "n_donors_with_mural": int((s > 0).sum()),
            "median_mural": float(s.median()),
            "mean_mural": float(s.mean()),
            "pct_donors_ge_5": float((s >= 5).mean()),
            "pct_donors_ge_10": float((s >= 10).mean()),
        }

    b0i_pass = all(
        pair_stats.get(r, {}).get("median_mural", 0) >= 5
        and pair_stats.get(r, {}).get("pct_donors_ge_5", 0) >= 0.5
        for r in MATCHED_PAIR
        if r in pair_stats
    )
    b0i_verdict = (
        "PASS_MURAL_COUNTS_ADEQUATE_AMONG_REGION_PRESENT_DONORS"
        if b0i_pass
        else (
            "FAIL_SINGLE_DIGIT_MURAL_COUNTS__TRACK_B_BLOCKED"
            if any(pair_stats.get(r, {}).get("median_mural", 0) < 5 for r in MATCHED_PAIR if r in pair_stats)
            else "MARGINAL_MURAL_COUNTS"
        )
    )

    # ----- B0ii HOX test -----
    # Need expression matrix — load into memory for mural cells only
    log("HOX module: subset mural cells")
    mural_idx = np.where(obs["is_mural"].values)[0]
    # Get X
    if mural_idx.size == 0:
        hox_res = {"verdict": "NO_MURAL_CELLS"}
    else:
        Xm = a.X[mural_idx]
        if hasattr(Xm, "toarray"):
            Xm = Xm.toarray()
        Xm = np.asarray(Xm, dtype=float)
        # log1p normalize per cell if looks like counts
        if Xm.max() > 50:
            lib = Xm.sum(axis=1, keepdims=True)
            lib[lib == 0] = 1
            Xm = np.log1p(Xm / lib * 1e4)

        hox_present = [g for g in HOX_GENES if g in gene_to_i]
        hox_i = [gene_to_i[g] for g in hox_present]
        hox_score = module_score(Xm, hox_i)
        mural_meta = obs.iloc[mural_idx].copy()
        mural_meta["hox_score"] = hox_score

        ant = mural_meta.loc[mural_meta["pole"] == "anterior", "hox_score"]
        post = mural_meta.loc[mural_meta["pole"] == "posterior", "hox_score"]
        # Donor-level means
        donor_hox = (
            mural_meta.groupby(["donor", "pole"], observed=True)["hox_score"]
            .mean()
            .reset_index()
            .pivot(index="donor", columns="pole", values="hox_score")
        )
        # Matched pair Midtemporal vs Entorhinal
        pair_hox = (
            mural_meta[mural_meta["region"].isin(MATCHED_PAIR)]
            .groupby(["donor", "region"], observed=True)["hox_score"]
            .mean()
            .reset_index()
            .pivot(index="donor", columns="region", values="hox_score")
        )

        def _mw(a, b):
            a = a.dropna()
            b = b.dropna()
            if len(a) < 5 or len(b) < 5:
                return {"n_a": len(a), "n_b": len(b), "p": np.nan, "stat": np.nan}
            st = stats.mannwhitneyu(post if False else b, a, alternative="two-sided")
            # Expect posterior > anterior HOX under lineage model
            return {
                "n_a": int(len(a)),
                "n_b": int(len(b)),
                "mean_a": float(a.mean()),
                "mean_b": float(b.mean()),
                "stat": float(st.statistic),
                "p": float(st.pvalue),
            }

        cell_mw = _mw(ant, post)
        # donor-level: posterior vs anterior
        if {"anterior", "posterior"}.issubset(donor_hox.columns):
            # paired donors present in both
            both = donor_hox.dropna(subset=["anterior", "posterior"])
            if len(both) >= 5:
                w = stats.wilcoxon(both["posterior"], both["anterior"])
                donor_test = {
                    "n_paired_donors": int(len(both)),
                    "mean_anterior": float(both["anterior"].mean()),
                    "mean_posterior": float(both["posterior"].mean()),
                    "wilcoxon_stat": float(w.statistic),
                    "p": float(w.pvalue),
                    "direction_expected": "posterior_HOX_gt_anterior",
                    "direction_observed": (
                        "posterior_gt_anterior"
                        if both["posterior"].mean() > both["anterior"].mean()
                        else "anterior_ge_posterior"
                    ),
                }
            else:
                donor_test = {"n_paired_donors": int(len(both)), "p": np.nan}
        else:
            donor_test = {"error": "missing poles"}

        pair_test = {}
        if set(MATCHED_PAIR).issubset(pair_hox.columns):
            bothp = pair_hox.dropna()
            if len(bothp) >= 5:
                w = stats.wilcoxon(bothp[MATCHED_PAIR[1]], bothp[MATCHED_PAIR[0]])
                pair_test = {
                    "pair": f"{MATCHED_PAIR[0]}(MCA) vs {MATCHED_PAIR[1]}(PCA)",
                    "n_paired_donors": int(len(bothp)),
                    "mean_midtemporal": float(bothp[MATCHED_PAIR[0]].mean()),
                    "mean_entorhinal": float(bothp[MATCHED_PAIR[1]].mean()),
                    "wilcoxon_p": float(w.pvalue),
                    "direction_observed": (
                        "entorhinal_gt_midtemporal"
                        if bothp[MATCHED_PAIR[1]].mean() > bothp[MATCHED_PAIR[0]].mean()
                        else "midtemporal_ge_entorhinal"
                    ),
                }

        # Detection: adult mural HOX RNA is often near-zero (lineage may be chromatin-only)
        hox_mat = Xm[:, hox_i] if hox_i else np.zeros((Xm.shape[0], 1))
        det_rate = float((hox_mat > 0).any(axis=1).mean()) if hox_i else 0.0
        mean_hox = float(hox_mat.mean()) if hox_i else 0.0

        hox_support = (
            donor_test.get("direction_observed") == "posterior_gt_anterior"
            and donor_test.get("p", 1) < 0.05
            and mean_hox > 1e-4
        )
        if det_rate < 0.05 or mean_hox < 1e-4:
            hox_verdict = (
                "HOX_RNA_UNDETECTABLE_IN_ADULT_MURAL__LINEAGE_NEEDS_CHROMATIN__B0ii_RNA_NOGO"
            )
        elif hox_support:
            hox_verdict = "PASS_HOX_POSTERIOR_ENRICHED__LINEAGE_SIGNAL"
        else:
            hox_verdict = "FAIL_OR_NULL_HOX_LINEAGE_TEST__TRACK_B_WEAK"

        hox_res = {
            "hox_genes_present": hox_present,
            "n_hox_genes_present": len(hox_present),
            "mural_cell_detection_rate_any_HOX": det_rate,
            "mural_mean_HOX_module": mean_hox,
            "cell_level_MW_anterior_vs_posterior": cell_mw,
            "donor_level_wilcoxon": donor_test,
            "matched_pair_midtemporal_vs_entorhinal": pair_test,
            "verdict": hox_verdict,
            "scientific_note": (
                "Concept note: lineage may persist as chromatin state, not adult RNA label. "
                "Near-zero HOX RNA is compatible with that and blocks RNA-only B0ii."
            ),
        }
        mural_meta.to_csv(OUT / "b0_mural_hox_cell_scores.csv", index=False)
        donor_hox.to_csv(OUT / "b0_mural_hox_donor_by_pole.csv")

    # ----- Track C: endothelial flow programs by territory prior -----
    log("Track C: endothelial flow module by region/territory")
    endo_idx = np.where(obs["is_endo"].values)[0]
    Xe = a.X[endo_idx]
    if hasattr(Xe, "toarray"):
        Xe = Xe.toarray()
    Xe = np.asarray(Xe, dtype=float)
    if Xe.max() > 50:
        lib = Xe.sum(axis=1, keepdims=True)
        lib[lib == 0] = 1
        Xe = np.log1p(Xe / lib * 1e4)

    flow_present = [g for g in FLOW_GENES if g in gene_to_i]
    flow_i = [gene_to_i[g] for g in flow_present]
    flow_score = module_score(Xe, flow_i)
    endo_meta = obs.iloc[endo_idx].copy()
    endo_meta["flow_score"] = flow_score
    # KLF2/KLF4 alone
    for g in ["KLF2", "KLF4", "NOS3", "ICAM1", "VCAM1"]:
        if g in gene_to_i:
            endo_meta[g] = Xe[:, gene_to_i[g]]

    # Donor-region means
    endo_dr = (
        endo_meta.groupby(["donor", "region", "territory_prior", "pole"], observed=True)
        .agg(n=("flow_score", "size"), flow_score=("flow_score", "mean"), **{
            g: (g, "mean") for g in ["KLF2", "KLF4", "NOS3"] if g in endo_meta.columns
        })
        .reset_index()
    )
    endo_dr.to_csv(OUT / "trackC_endo_flow_donor_region.csv", index=False)

    # Anterior vs posterior flow (expect difference if haemodynamic)
    ant_f = endo_dr.loc[endo_dr["pole"] == "anterior", "flow_score"]
    post_f = endo_dr.loc[endo_dr["pole"] == "posterior", "flow_score"]
    if len(ant_f) >= 5 and len(post_f) >= 5:
        st = stats.mannwhitneyu(ant_f, post_f, alternative="two-sided")
        flow_pole = {
            "mean_anterior": float(ant_f.mean()),
            "mean_posterior": float(post_f.mean()),
            "p": float(st.pvalue),
            "n_ant": int(len(ant_f)),
            "n_post": int(len(post_f)),
        }
    else:
        flow_pole = {}

    # MCA (midtemporal+angular) vs PCA (hippocampus+entorhinal) endo flow
    mca_regs = {"Midtemporal_cortex", "Angular_gyrus"}
    pca_regs = {"Hippocampus", "Entorhinal_cortex"}
    mca_f = endo_dr.loc[endo_dr["region"].isin(mca_regs), "flow_score"]
    pca_f = endo_dr.loc[endo_dr["region"].isin(pca_regs), "flow_score"]
    if len(mca_f) >= 5 and len(pca_f) >= 5:
        st = stats.mannwhitneyu(mca_f, pca_f, alternative="two-sided")
        flow_terr = {
            "mean_MCA_regions": float(mca_f.mean()),
            "mean_PCA_regions": float(pca_f.mean()),
            "p": float(st.pvalue),
        }
    else:
        flow_terr = {}

    trackC_verdict = (
        "FLOW_PROGRAM_DIFFERS_BY_POLE_OR_TERRITORY__ARCHITECTURE_READOUT_PLAUSIBLE"
        if (flow_pole.get("p", 1) < 0.05 or flow_terr.get("p", 1) < 0.05)
        else "FLOW_PROGRAM_NULL_BETWEEN_POLES__TRACK_C_INCONCLUSIVE_ON_THIS_OBJECT"
    )

    summary = {
        "data_source": {
            "path": str(H5AD),
            "note": (
                "Local Garcia ROSMAP multiregion vascular object used for B0/Track C. "
                "Sun 2023 MIT counts_full.rds download was incomplete/unreadable (exact 32MiB); "
                "full Sun release requires Dropbox/Synapse (syn51015750). "
                "Regions overlap Sun design; cell n is smaller than published 22,514."
            ),
            "n_cells": int(a.n_obs),
            "n_donors": int(obs["donor"].nunique()),
            "celltype_counts": obs["celltype"].value_counts().to_dict(),
            "region_counts": obs["region"].value_counts().to_dict(),
        },
        "B0i_mural_counts": {
            "medians_by_region": medians,
            "matched_pair_stats": pair_stats,
            "verdict": b0i_verdict,
        },
        "B0ii_HOX": hox_res,
        "TrackC_flow": {
            "flow_genes_present": flow_present,
            "pole_contrast": flow_pole,
            "MCA_vs_PCA_region_contrast": flow_terr,
            "verdict": trackC_verdict,
        },
        "TrackB_gate": {
            "vessel_enriched_ACC": "CLOSED_Heiman_Matthews__PENDING_Walchli",
            "sun_full_release": "BLOCKED_incomplete_MIT_RDS__need_Dropbox_or_Synapse",
            "chromatin": "NOT_AVAILABLE",
            "verdict": "TRACK_B_CLOSED__DO_NOT_SPEND",
        },
        "overall": None,
    }

    # Overall routing
    b0i_v = b0i_verdict
    hox_v = hox_res.get("verdict", "")
    if "FAIL_SINGLE_DIGIT" in b0i_v:
        summary["overall"] = "B0i_FAIL__HOLD_TRACK_B__TRACK_A_AND_C_PRIMARY"
    elif "HOX_RNA_UNDETECTABLE" in hox_v and "PASS" in b0i_v:
        summary["overall"] = (
            "B0i_PASS_COUNTS__B0ii_HOX_RNA_NOGO__TRACK_B_NEEDS_CHROMATIN__TRACK_C_PRIMARY"
        )
    elif "PASS" in b0i_v and "PASS_HOX" in hox_v:
        summary["overall"] = "B0_PASS__TRACK_B_ELIGIBLE_IF_ACC_AND_CHROMATIN_OPEN"
    else:
        summary["overall"] = "B0_MIXED__PREFER_TRACK_A_AND_C__HOLD_TRACK_B"

    with open(OUT / "b0_trackC_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=str)

    _report(summary)
    log(f"Wrote {OUT}")
    return summary


def _report(s: dict):
    lines = [
        "# B0 gates + Track C (flow programs)",
        "",
        f"**Date:** {time.strftime('%Y-%m-%d')}",
        f"**Overall routing:** `{s['overall']}`",
        "",
        "## Data",
        f"- Source: `{s['data_source']['path']}`",
        f"- Note: {s['data_source']['note']}",
        f"- n_cells={s['data_source']['n_cells']}, n_donors={s['data_source']['n_donors']}",
        "",
        "## B0i — mural counts (Per+SMC) per region per donor",
        f"- Verdict: `{s['B0i_mural_counts']['verdict']}`",
        f"- Medians: `{json.dumps(s['B0i_mural_counts']['medians_by_region'])}`",
        f"- Matched pair: `{json.dumps(s['B0i_mural_counts']['matched_pair_stats'])}`",
        "",
        "## B0ii — HOX anterior vs posterior (lineage proxy)",
        f"- Verdict: `{s['B0ii_HOX'].get('verdict')}`",
        f"- Donor Wilcoxon: `{json.dumps(s['B0ii_HOX'].get('donor_level_wilcoxon'))}`",
        f"- Matched pair: `{json.dumps(s['B0ii_HOX'].get('matched_pair_midtemporal_vs_entorhinal'))}`",
        "",
        "## Track C — endothelial flow/shear programs",
        f"- Verdict: `{s['TrackC_flow']['verdict']}`",
        f"- Pole: `{json.dumps(s['TrackC_flow']['pole_contrast'])}`",
        f"- MCA vs PCA regions: `{json.dumps(s['TrackC_flow']['MCA_vs_PCA_region_contrast'])}`",
        "",
        "## Track B gate",
        f"- `{json.dumps(s['TrackB_gate'])}`",
        "",
        "## Locks",
        "- TRCI/BCI not mixed into these tracks.",
        "- Track A imaging partition reported separately.",
        "- RADC Aim 3 still DUA-gated.",
        "",
    ]
    (OUT / "b0_trackC_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
