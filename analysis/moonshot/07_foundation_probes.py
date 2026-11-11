"""Moonshot #7 — Foundation-model-style territorial probes aligned to GWAS programs.

Primary: donor-OOF linear probes on Class B Harmony embeddings + MultiVINE Harmony
gene-program projection. Secondary: gene-embedding probe (GenePT-style mean pooling
from a downloadable gene vector table if available; else skip with clear status).

Honesty: methods lever; with #1 as supplement — not standalone Nat Neurosci claim.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from _paths import CLASS_B, GWAS_DIR, HONESTY, MULTIVINE, RES

OUT = RES / "07_foundation_probes"
OUT.mkdir(parents=True, exist_ok=True)
EMB_DIR = Path(r"F:\Cerebrovascular SCS\data\raw\gene_embeddings")
EMB_DIR.mkdir(parents=True, exist_ok=True)


def _donor_oof(X, y, groups) -> dict:
    mask = np.isfinite(y)
    X, y, groups = X[mask], y[mask], groups[mask]
    if len(np.unique(groups)) < 3 or len(np.unique(y)) < 2:
        return {"auroc": float("nan"), "n": int(len(y))}
    gkf = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    preds = np.full(len(y), np.nan)
    for tr, te in gkf.split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        sc = StandardScaler()
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        clf.fit(sc.fit_transform(X[tr]), y[tr])
        preds[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    ok = ~np.isnan(preds)
    return {
        "auroc": float(roc_auc_score(y[ok], preds[ok])) if ok.sum() > 10 and len(np.unique(y[ok])) > 1 else float("nan"),
        "n": int(ok.sum()),
        "n_donors": int(len(np.unique(groups[ok]))) if ok.any() else 0,
    }


def _try_download_genept() -> Path | None:
    """Attempt small public gene embedding download; return path or None."""
    # GenePT ChatGPT-3.5 gene embeddings (parquet/csv mirrors vary). Try a known raw github release.
    candidates = [
        (
            "https://raw.githubusercontent.com/yiqunchen/GenePT/main/data/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle",
            EMB_DIR / "GenePT_gene_embedding_ada_text.pickle",
        ),
    ]
    for url, dest in candidates:
        if dest.exists() and dest.stat().st_size > 1000:
            return dest
        try:
            print("Downloading", url)
            urllib.request.urlretrieve(url, dest)
            if dest.stat().st_size > 1000:
                return dest
        except Exception as e:
            print("gene emb download fail", e)
    return None


def _map_terr(obs: pd.DataFrame) -> np.ndarray:
    if "tv_dominant" in obs.columns:
        t = obs["tv_dominant"].astype(str)
        y = np.full(len(t), np.nan)
        y[t == "ACA"] = 1
        y[t == "MCA"] = 0
        # expand
        for i, v in enumerate(t):
            if np.isnan(y[i]):
                if any(k in v for k in ["ACA", "Allen", "LIBD", "ACC"]):
                    y[i] = 1
                elif any(k in v for k in ["MCA", "Garcia"]):
                    y[i] = 0
        return y
    return np.full(len(obs), np.nan)


def run() -> dict:
    results = {"honesty_lock": HONESTY["7"]}
    cb = ad.read_h5ad(CLASS_B)
    y = _map_terr(cb.obs)
    donor = cb.obs["donor"].astype(str).values if "donor" in cb.obs else np.array([f"d{i//10}" for i in range(cb.n_obs)])

    # Probe A: Harmony PCA embeddings
    if "X_pca_harmony" in cb.obsm:
        Zh = np.asarray(cb.obsm["X_pca_harmony"])
    else:
        X = cb.X.toarray() if sparse.issparse(cb.X) else np.asarray(cb.X)
        from sklearn.decomposition import PCA

        Zh = PCA(n_components=min(30, X.shape[1]), random_state=0).fit_transform(X)

    probe_h = _donor_oof(Zh, y, donor)
    results["probe_harmony_pca"] = probe_h

    # Probe B: GWAS gene-module scores as "foundation program" features
    genesets = {}
    for p in GWAS_DIR.glob("geneset_*.txt"):
        genesets[p.stem.replace("geneset_", "")] = {
            g.strip() for g in p.read_text(encoding="utf-8").splitlines() if g.strip()
        }
    genes = list(cb.var_names.astype(str))
    X = cb.X
    cols = []
    used = []
    for name, gs in genesets.items():
        idx = [i for i, g in enumerate(genes) if g in gs]
        if len(idx) < 3:
            continue
        sub = X[:, idx]
        cols.append(np.asarray(sub.mean(axis=1)).ravel() if sparse.issparse(sub) else sub.mean(axis=1))
        used.append(name)
    if cols:
        Xg = np.column_stack(cols)
        probe_g = _donor_oof(Xg, y, donor)
        results["probe_gwas_modules"] = probe_g
        results["gwas_modules_used"] = used
        # Align: correlation of probe coefficients with MultiVINE marker overlap — done after MV
    else:
        probe_g = {"auroc": float("nan")}
        results["probe_gwas_modules"] = probe_g

    # Probe C: MultiVINE vascular Harmony space transfer — train on MV arterial vs venous,
    # apply gene weights... (different cells). Instead: project Class B genes onto MultiVINE
    # cell-type centroids for Arterial/Cap/Venous and use distances as features.
    print("Loading MultiVINE backed for centroids…")
    mv = ad.read_h5ad(MULTIVINE, backed="r")
    if "feature_name" in mv.var.columns:
        mv_genes = mv.var["feature_name"].astype(str).values
    else:
        mv_genes = mv.var_names.astype(str).values
    shared = [g for g in genes if g in set(mv_genes)]
    shared = shared[:3000]
    results["n_shared_genes_mv"] = len(shared)
    # Use precomputed X_harmony on MultiVINE for vascular cells — centroid by author_cell_type
    # Can't easily slice backed for gene means without loading; load obs + harmony only
    ct = mv.obs["author_cell_type"].astype(str)
    H = np.asarray(mv.obsm["X_harmony"])
    centroids = {}
    for t in ["Arterial", "Capillary", "Venous", "Pericyte", "Smooth_Muscle_Cell"]:
        m = (ct == t).values
        if m.sum() > 50:
            centroids[t] = H[m].mean(axis=0)
    # For Class B, we don't have MultiVINE harmony. Use GWAS+Harmony concat probe instead.
    if cols:
        Xcat = np.hstack([Zh, Xg])
        results["probe_harmony_plus_gwas"] = _donor_oof(Xcat, y, donor)

    # Alignment metric: gene importance from logistic on GWAS modules vs enrichment from #1
    enr_path = RES / "01_gwas_multivine" / "gwas_enrichment_by_celltype.csv"
    alignment = {}
    if enr_path.exists() and cols:
        enr = pd.read_csv(enr_path)
        # Fit full-data logistic for coefficients on modules (descriptive)
        mask = np.isfinite(y)
        sc = StandardScaler()
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        clf.fit(sc.fit_transform(Xg[mask]), y[mask])
        coef = pd.Series(clf.coef_.ravel(), index=used)
        # Compare absolute coef to -log10p best vascular enrichment per geneset
        best_p = enr.groupby("geneset")["hypergeom_p"].min()
        common = [g for g in used if g in best_p.index]
        if len(common) >= 3:
            a = coef.loc[common].abs().values
            b = (-np.log10(np.clip(best_p.loc[common].values, 1e-300, 1)))
            r, p = stats.spearmanr(a, b)
            alignment = {"spearman_coef_vs_multivine_enrichment": float(r), "p": float(p), "n": len(common)}
        results["alignment_to_multivine_gwas"] = alignment

    # Optional GenePT
    emb_path = _try_download_genept()
    if emb_path is not None:
        try:
            import pickle

            with open(emb_path, "rb") as f:
                gene_emb = pickle.load(f)
            # gene_emb: dict gene -> vector
            dim = len(next(iter(gene_emb.values())))
            # cell embedding = mean of expressed gene embeddings (binarize top genes)
            Xbin = X
            if sparse.issparse(Xbin):
                # use top 500 expressed genes per cell is expensive; use gene-weight by mean expr
                mean_expr = np.asarray(Xbin.mean(axis=0)).ravel()
            else:
                mean_expr = Xbin.mean(axis=0)
            # Build gene embedding matrix for shared genes present in gene_emb
            g_idx = []
            E = []
            for i, g in enumerate(genes):
                if g in gene_emb:
                    g_idx.append(i)
                    E.append(np.asarray(gene_emb[g], dtype=float))
            if len(E) > 100:
                E = np.vstack(E)
                # Weighted cell embeddings: X[:,g_idx] @ E
                Xs = X[:, g_idx]
                if sparse.issparse(Xs):
                    cell_emb = Xs @ E
                    cell_emb = np.asarray(cell_emb)
                else:
                    cell_emb = Xs @ E
                results["probe_genept_meanpool"] = _donor_oof(cell_emb, y, donor)
                results["genept_n_genes"] = len(E)
            else:
                results["genept_status"] = "too_few_genes_in_embedding"
        except Exception as e:
            results["genept_status"] = f"failed:{e}"
    else:
        results["genept_status"] = "DOWNLOAD_UNAVAILABLE — Harmony+GWAS probes primary"

    # Null: shuffle territory labels across donors (within-donor shuffle is a no-op
    # when each donor is territory-pure)
    rng = np.random.default_rng(0)
    donors_u = np.unique(donor)
    # donor-level majority label
    donor_lab = {}
    for d in donors_u:
        m = donor == d
        vals = y[m]
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            donor_lab[d] = np.nan
        else:
            donor_lab[d] = float(np.round(np.mean(vals)))
    labs = np.array([donor_lab[d] for d in donors_u])
    perm = rng.permutation(labs)
    donor_lab_shuf = {d: perm[i] for i, d in enumerate(donors_u)}
    y_shuf = np.array([donor_lab_shuf[d] for d in donor], dtype=float)
    results["probe_harmony_donor_label_shuffle_null"] = _donor_oof(Zh, y_shuf, donor)
    if cols:
        results["probe_gwas_donor_label_shuffle_null"] = _donor_oof(Xg, y_shuf, donor)

    (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(
        "\n".join(
            [
                "# Moonshot 07 — Foundation-style territorial probes",
                "",
                f"**Honesty lock:** {HONESTY['7']}",
                "",
                f"```json\n{json.dumps(results, indent=2)}\n```",
            ]
        ),
        encoding="utf-8",
    )
    print("DONE 07", results)
    return results


if __name__ == "__main__":
    run()
