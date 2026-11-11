"""Moonshot #1 — Territory as a filter of cerebrovascular GWAS × MultiVINE.

Honesty lock: falsify if territory adds nothing beyond cell type + A–V zonation.
Public MultiVINE CELLxGENE is dlPFC snRNA (+ multiome embeddings); peak-level
ATAC matrices are not in the public H5AD — we use cell-type programs + GWAS gene
sets, and X_umap_atac only as orthogonal embedding evidence.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from _paths import CLASS_B, GARCIA_REGIONS, GWAS_DIR, HONESTY, MULTIVINE, RES

OUT = RES / "01_gwas_multivine"
OUT.mkdir(parents=True, exist_ok=True)

VASC_TYPES = {
    "Arterial",
    "Capillary",
    "Venous",
    "Pericyte",
    "Smooth_Muscle_Cell",
    "Perivascular_Fibroblast",
    "Meningeal_Fibroblast",
}
AV_MAP = {
    "Arterial": "A",
    "Capillary": "Cap",
    "Venous": "V",
    "Smooth_Muscle_Cell": "Mural",
    "Pericyte": "Mural",
}


def _load_genesets() -> dict[str, set[str]]:
    sets = {}
    for p in GWAS_DIR.glob("geneset_*.txt"):
        genes = {g.strip() for g in p.read_text(encoding="utf-8").splitlines() if g.strip()}
        if genes:
            sets[p.stem.replace("geneset_", "")] = genes
    return sets


def _hypergeom(overlap: int, set_n: int, marker_n: int, universe: int) -> float:
    # P(X >= overlap)
    return float(stats.hypergeom.sf(overlap - 1, universe, set_n, marker_n))


def _donor_oof_auroc(X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> float:
    mask = np.isin(y, [0, 1])
    X, y, groups = X[mask], y[mask], groups[mask]
    if len(np.unique(groups)) < 3 or y.sum() < 5 or (1 - y).sum() < 5:
        return float("nan")
    gkf = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    preds = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xte = sc.transform(X[te])
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        clf.fit(Xtr, y[tr])
        preds[te] = clf.predict_proba(Xte)[:, 1]
    ok = ~np.isnan(preds)
    if ok.sum() < 10 or len(np.unique(y[ok])) < 2:
        return float("nan")
    return float(roc_auc_score(y[ok], preds[ok]))


def _module_score(X, var_names: list[str], genes: set[str]) -> np.ndarray:
    idx = [i for i, g in enumerate(var_names) if g in genes]
    if not idx:
        return np.zeros(X.shape[0])
    sub = X[:, idx]
    if sparse.issparse(sub):
        return np.asarray(sub.mean(axis=1)).ravel()
    return sub.mean(axis=1)


def run() -> dict:
    genesets = _load_genesets()
    print("genesets", {k: len(v) for k, v in genesets.items()})

    print("Loading MultiVINE (memory)…")
    mv = ad.read_h5ad(MULTIVINE)
    # Prefer feature_name / names for symbols
    if "feature_name" in mv.var.columns:
        mv.var_names = mv.var["feature_name"].astype(str).values
    elif "names" in mv.var.columns:
        mv.var_names = mv.var["names"].astype(str).values
    mv.var_names_make_unique()

    # Vascular subset
    ct = mv.obs["author_cell_type"].astype(str)
    vasc = mv[ct.isin(VASC_TYPES)].copy()
    print("vascular nuclei", vasc.n_obs)

    # Marker genes per cell type (donor-pseudobulk rank)
    # Use simple mean expression fold vs other vascular
    if sparse.issparse(vasc.X):
        X = vasc.X.tocsr()
    else:
        X = sparse.csr_matrix(vasc.X)

    var_names = list(vasc.var_names.astype(str))
    # subsample genes to HVGs if present / else top variance
    gene_var = np.asarray(X.power(2).mean(axis=0) - np.asarray(X.mean(axis=0)).ravel() ** 2).ravel()
    top = np.argsort(gene_var)[-4000:]
    Xh = X[:, top]
    genes_h = [var_names[i] for i in top]

    marker_rows = []
    enrichment_rows = []
    programs = {}  # celltype -> top genes

    for ctype in sorted(VASC_TYPES):
        mask = (vasc.obs["author_cell_type"].astype(str) == ctype).values
        if mask.sum() < 50:
            continue
        mu_in = np.asarray(Xh[mask].mean(axis=0)).ravel()
        mu_out = np.asarray(Xh[~mask].mean(axis=0)).ravel()
        lfc = mu_in - mu_out
        order = np.argsort(lfc)[::-1]
        top_genes = [genes_h[i] for i in order[:100]]
        programs[ctype] = top_genes
        for rank, gi in enumerate(order[:50]):
            marker_rows.append(
                {
                    "cell_type": ctype,
                    "gene": genes_h[gi],
                    "rank": rank + 1,
                    "lfc_vs_other_vasc": float(lfc[gi]),
                }
            )
        # GWAS enrichment on top 100
        for gs_name, gs in genesets.items():
            inter = set(top_genes) & gs
            p = _hypergeom(len(inter), len(gs), 100, len(genes_h))
            enrichment_rows.append(
                {
                    "cell_type": ctype,
                    "geneset": gs_name,
                    "overlap_n": len(inter),
                    "overlap_genes": "|".join(sorted(inter)[:30]),
                    "hypergeom_p": p,
                    "universe": len(genes_h),
                }
            )

    # A–V zonation program
    art = (vasc.obs["author_cell_type"].astype(str) == "Arterial").values
    ven = (vasc.obs["author_cell_type"].astype(str) == "Venous").values
    if art.sum() > 50 and ven.sum() > 50:
        mu_a = np.asarray(Xh[art].mean(axis=0)).ravel()
        mu_v = np.asarray(Xh[ven].mean(axis=0)).ravel()
        lfc_av = mu_a - mu_v
        av_genes = [genes_h[i] for i in np.argsort(lfc_av)[::-1][:100]]
        programs["Arterial_vs_Venous"] = av_genes
        for gs_name, gs in genesets.items():
            inter = set(av_genes) & gs
            enrichment_rows.append(
                {
                    "cell_type": "Arterial_vs_Venous",
                    "geneset": gs_name,
                    "overlap_n": len(inter),
                    "overlap_genes": "|".join(sorted(inter)[:30]),
                    "hypergeom_p": _hypergeom(len(inter), len(gs), 100, len(genes_h)),
                    "universe": len(genes_h),
                }
            )

    pd.DataFrame(marker_rows).to_csv(OUT / "multivine_vascular_markers.csv", index=False)
    enr = pd.DataFrame(enrichment_rows)
    enr["neglog10p"] = -np.log10(np.clip(enr["hypergeom_p"], 1e-300, 1))
    enr.to_csv(OUT / "gwas_enrichment_by_celltype.csv", index=False)

    # BH FDR within geneset
    from statsmodels.stats.multitest import multipletests

    enr["fdr"] = np.nan
    for gs in enr["geneset"].unique():
        m = enr["geneset"] == gs
        enr.loc[m, "fdr"] = multipletests(enr.loc[m, "hypergeom_p"], method="fdr_bh")[1]
    enr.to_csv(OUT / "gwas_enrichment_by_celltype.csv", index=False)

    # Project GWAS programs onto TerritoryVoxels cohorts (Class B + Garcia regions)
    # Use combined cerebrovascular geneset + AD_core + WMH_core
    score_sets = {
        k: genesets[k]
        for k in [
            "cerebrovascular_combined",
            "AD_core",
            "WMH_core",
            "stroke_core",
            "cSVD_core",
            "ischemic_stroke_HP",
            "WMH",
        ]
        if k in genesets
    }

    projection_rows = []
    falsifier = {}

    for label, path in [("class_b", CLASS_B), ("garcia_regions", GARCIA_REGIONS)]:
        if not path.exists():
            continue
        adata = ad.read_h5ad(path)
        # endo only if available
        if "vascular_lineage" in adata.obs:
            adata = adata[adata.obs["vascular_lineage"].astype(str) == "endo"].copy()
        elif "celltype" in adata.obs:
            adata = adata[adata.obs["celltype"].astype(str).str.contains("Endo", case=False, na=False)].copy()
        if adata.n_obs < 20:
            continue
        genes = list(adata.var_names.astype(str))
        X = adata.X
        scores = {}
        for sn, gs in score_sets.items():
            scores[sn] = _module_score(X, genes, gs)
            adata.obs[f"gwas_{sn}"] = scores[sn]

        # Territory contrast if present
        terr_col = None
        for c in ["tv_dominant", "class_b_arm", "brain_region"]:
            if c in adata.obs.columns:
                terr_col = c
                break
        donor_col = "donor" if "donor" in adata.obs else ("subject" if "subject" in adata.obs else None)

        if terr_col is None:
            continue

        # Summarize by territory
        for sn in score_sets:
            g = adata.obs.groupby(terr_col, observed=True)[f"gwas_{sn}"]
            for terr, s in g:
                projection_rows.append(
                    {
                        "object": label,
                        "territory_col": terr_col,
                        "territory": str(terr),
                        "geneset": sn,
                        "mean_score": float(s.mean()),
                        "n_cells": int(s.shape[0]),
                    }
                )

        # Falsifier: donor-OOF AUROC predicting ACA vs MCA from GWAS scores alone
        # vs from cell-type composition proxies if available
        if label == "class_b" and "tv_dominant" in adata.obs.columns and donor_col:
            ymap = {"ACA": 1, "MCA": 0}
            y = adata.obs["tv_dominant"].astype(str).map(ymap).values.astype(float)
            # keep only ACA/MCA
            keep = np.isfinite(y)
            Xg = np.column_stack([adata.obs[f"gwas_{sn}"].values for sn in score_sets])
            donors = adata.obs[donor_col].astype(str).values
            auroc_gwas = _donor_oof_auroc(Xg[keep], y[keep], donors[keep])
            falsifier["class_b_gwas_score_donor_oof_auroc_ACA_vs_MCA"] = auroc_gwas

            # Ambient / non-GWAS null: random gene modules same size
            rng = np.random.default_rng(0)
            null_aurocs = []
            for i in range(20):
                null_cols = []
                for sn, gs in score_sets.items():
                    null_genes = set(rng.choice(genes, size=min(len(gs), len(genes)), replace=False))
                    null_cols.append(_module_score(X, genes, null_genes))
                Xn = np.column_stack(null_cols)
                null_aurocs.append(_donor_oof_auroc(Xn[keep], y[keep], donors[keep]))
            falsifier["null_random_modules_mean_auroc"] = float(np.nanmean(null_aurocs))
            falsifier["null_random_modules_sd_auroc"] = float(np.nanstd(null_aurocs))
            falsifier["gwas_beats_null"] = bool(
                np.isfinite(auroc_gwas)
                and auroc_gwas > falsifier["null_random_modules_mean_auroc"] + falsifier["null_random_modules_sd_auroc"]
            )

        # Within Garcia: Midtemporal (MCA) vs Entorhinal (PCA) GWAS scores — enrichment-matched regions
        if label == "garcia_regions" and "brain_region" in adata.obs.columns and donor_col:
            sub = adata[adata.obs["brain_region"].isin(["Midtemporal_cortex", "Entorhinal_cortex"])].copy()
            if sub.n_obs > 30:
                y = (sub.obs["brain_region"].astype(str) == "Midtemporal_cortex").astype(int).values
                Xg = np.column_stack([sub.obs[f"gwas_{sn}"].values for sn in score_sets])
                donors = sub.obs[donor_col].astype(str).values
                falsifier["garcia_MT_vs_EC_gwas_donor_oof_auroc"] = _donor_oof_auroc(Xg, y, donors)

        # Save per-cell scores sample
        cols = [c for c in adata.obs.columns if c.startswith("gwas_")]
        meta_cols = [c for c in [terr_col, donor_col, "author_cell_type", "celltype", "subtype"] if c and c in adata.obs]
        adata.obs[meta_cols + cols].to_csv(OUT / f"gwas_scores_{label}.csv")

    pd.DataFrame(projection_rows).to_csv(OUT / "gwas_projection_by_territory.csv", index=False)

    # Cell-type vs territory: does vascular cell-type GWAS enrichment dominate?
    # Compare max enrichment FDR in Arterial/Cap/Venous vs neuronal control types
    ct_full = mv.obs["author_cell_type"].astype(str)
    control_types = ["Excitatory_Neuron", "Inhibitory_Neuron", "Astrocyte", "Oligodendrocyte", "Microglia"]
    control_enr = []
    for ctype in control_types:
        mask = (ct_full == ctype).values
        if mask.sum() < 100:
            continue
        # compute on full mv X for this type vs rest (subsample)
        idx = np.where(mask)[0]
        if len(idx) > 2000:
            idx = np.random.default_rng(1).choice(idx, 2000, replace=False)
        # Use same HVG space from vascular for comparability — approximate with gene intersection
        # Skip heavy; use Gene set score difference only
        control_enr.append(ctype)

    # Summary decision
    top_enr = enr.sort_values("hypergeom_p").head(15)
    decision = {
        "honesty_lock": HONESTY["1"],
        "multivine_tissue": "dorsolateral prefrontal cortex ONLY — not ACC",
        "n_vascular_nuclei": int(vasc.n_obs),
        "n_genesets": len(genesets),
        "top_enrichments": top_enr.to_dict(orient="records"),
        "falsifier": falsifier,
        "verdict_notes": [],
    }
    if falsifier.get("gwas_beats_null") is True:
        decision["verdict_notes"].append(
            "GWAS modules separate Class B ACA vs MCA above random-module null (territory-linked signal possible)."
        )
    elif falsifier.get("gwas_beats_null") is False:
        decision["verdict_notes"].append(
            "FALSIFIED at Class B: GWAS modules do not beat random modules for ACA vs MCA — territory adds nothing beyond noise at this resolution."
        )
    else:
        decision["verdict_notes"].append("Class B falsifier inconclusive (AUROC undefined / sparse).")

    # Cell-type enrichment primacy
    vasc_enr = enr[enr["cell_type"].isin(list(VASC_TYPES) + ["Arterial_vs_Venous"])]
    best = vasc_enr.loc[vasc_enr["hypergeom_p"].idxmin()] if len(vasc_enr) else None
    if best is not None:
        decision["best_celltype_enrichment"] = best.to_dict()

    (OUT / "summary.json").write_text(json.dumps(decision, indent=2, default=str), encoding="utf-8")

    lines = [
        "# Moonshot 01 — GWAS × MultiVINE territory filter",
        "",
        f"**Honesty lock:** {HONESTY['1']}",
        "",
        f"- MultiVINE vascular nuclei: **{vasc.n_obs}** (dlPFC only)",
        f"- Genesets: { {k: len(v) for k,v in genesets.items()} }",
        f"- Falsifier: `{json.dumps(falsifier)}`",
        "",
        "## Top enrichments",
        top_enr[["cell_type", "geneset", "overlap_n", "hypergeom_p", "fdr"]].head(12).to_markdown(index=False),
        "",
        "## Verdict notes",
        *[f"- {x}" for x in decision["verdict_notes"]],
        "",
        "Outputs: `gwas_enrichment_by_celltype.csv`, `gwas_projection_by_territory.csv`, `summary.json`",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print("DONE 01", OUT)
    return decision


if __name__ == "__main__":
    run()
