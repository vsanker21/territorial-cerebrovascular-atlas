"""Moonshot #5 — Wasserstein territorial manifold + borderzone geodesics.

Upgrade H2 mixture claim into OT geometry; conformal territory labels.
"""
from __future__ import annotations

import json

import anndata as ad
import numpy as np
import ot
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from _paths import CLASS_B, HONESTY, RES, YANG

OUT = RES / "05_wasserstein_manifold"
OUT.mkdir(parents=True, exist_ok=True)


def _empirical_measure(Z: np.ndarray, max_n: int = 200, rng=None):
    rng = rng or np.random.default_rng(0)
    if len(Z) > max_n:
        Z = Z[rng.choice(len(Z), max_n, replace=False)]
    a = np.ones(len(Z)) / len(Z)
    return Z, a


def _w2(Z1, Z2, rng) -> float:
    X, a = _empirical_measure(Z1, rng=rng)
    Y, b = _empirical_measure(Z2, rng=rng)
    M = ot.dist(X, Y, metric="euclidean")
    # squared euclidean for W2^2 then sqrt
    M2 = M ** 2
    return float(np.sqrt(ot.emd2(a, b, M2)))


def _barycenter_interp(Z_a, Z_b, t: float, rng) -> np.ndarray:
    """Approximate geodesic midpoint via McCann interpolation on EMD coupling (paired)."""
    X, a = _empirical_measure(Z_a, max_n=150, rng=rng)
    Y, b = _empirical_measure(Z_b, max_n=150, rng=rng)
    # equalize sizes
    n = min(len(X), len(Y))
    X, Y = X[:n], Y[:n]
    a = np.ones(n) / n
    b = np.ones(n) / n
    M = ot.dist(X, Y, metric="euclidean")
    G = ot.emd(a, b, M)
    # sample pairs from coupling
    flat = G.ravel()
    flat = flat / flat.sum()
    idx = rng.choice(n * n, size=n, p=flat)
    i = idx // n
    j = idx % n
    return (1 - t) * X[i] + t * Y[j]


def run() -> dict:
    rng = np.random.default_rng(0)
    cb = ad.read_h5ad(CLASS_B)
    if "X_pca_harmony" in cb.obsm:
        Z = np.asarray(cb.obsm["X_pca_harmony"])[:, :20]
    else:
        X = cb.X.toarray() if sparse.issparse(cb.X) else np.asarray(cb.X)
        Z = PCA(n_components=20, random_state=0).fit_transform(StandardScaler().fit_transform(X))

    # Territory labels
    if "tv_dominant" in cb.obs:
        terr = cb.obs["tv_dominant"].astype(str).values
    else:
        terr = cb.obs["class_b_arm"].astype(str).values

    def map_t(t):
        t = str(t)
        if t in ("ACA", "MCA", "PCA", "VB"):
            return t
        if any(k in t for k in ["ACA", "Allen", "LIBD", "Linn", "ACC"]):
            return "ACA"
        if any(k in t for k in ["MCA", "Garcia", "MTG", "temporal"]):
            return "MCA"
        return "other"

    labels = np.array([map_t(t) for t in terr])
    territories = [t for t in ["ACA", "MCA", "PCA"] if (labels == t).sum() >= 20]
    if "ACA" not in territories or "MCA" not in territories:
        # force from class_b_arm binary
        territories = ["ACA", "MCA"]
        # already mapped

    # Pairwise W2
    dist_rows = []
    centers = {}
    for t in territories:
        centers[t] = Z[labels == t]
    for i, t1 in enumerate(territories):
        for t2 in territories[i:]:
            w = _w2(centers[t1], centers[t2], rng)
            dist_rows.append({"t1": t1, "t2": t2, "W2": w})
    dist_df = pd.DataFrame(dist_rows)
    dist_df.to_csv(OUT / "wasserstein_pairwise.csv", index=False)

    # H2 borderzone geodesic test with Yang SFG if available
    h2 = {}
    if YANG.exists():
        yang = ad.read_h5ad(YANG)
        # subsample endo
        if "cell_type" in yang.obs.columns:
            mask = yang.obs["cell_type"].astype(str).str.contains("Endo|endothelial", case=False, na=False)
            yang = yang[mask].copy() if mask.any() else yang
        elif "vascular_lineage" in yang.obs:
            yang = yang[yang.obs["vascular_lineage"].astype(str) == "endo"].copy()
        # Project Yang into Class B PCA space is hard without shared model;
        # Instead: use Class B if SFG-like arm exists, else skip gene alignment
        # Build shared PCA on concatenated expression for endo genes intersection
        shared = sorted(set(cb.var_names.astype(str)) & set(yang.var_names.astype(str)))
        if len(shared) >= 500:
            rng2 = np.random.default_rng(1)
            def sub_expr(a, n=400):
                if a.n_obs > n:
                    a = a[rng2.choice(a.n_obs, n, replace=False)].copy()
                X = a[:, shared].X
                X = X.toarray() if sparse.issparse(X) else np.asarray(X)
                return X
            Xa = sub_expr(cb[labels == "ACA"])
            Xm = sub_expr(cb[labels == "MCA"])
            Xy = sub_expr(yang)
            Xall = np.vstack([Xa, Xm, Xy])
            P = PCA(n_components=15, random_state=0).fit(StandardScaler().fit_transform(Xall))
            Za = P.transform(StandardScaler().fit_transform(Xa))
            Zm = P.transform(StandardScaler().fit_transform(Xm))
            Zy = P.transform(StandardScaler().fit_transform(Xy))
            # Geodesic: W2(SFG, midpoint ACA-MCA) vs W2(SFG, ACA) and W2(SFG, MCA)
            mid = _barycenter_interp(Za, Zm, 0.5, rng)
            w_mid = _w2(Zy, mid, rng)
            w_aca = _w2(Zy, Za, rng)
            w_mca = _w2(Zy, Zm, rng)
            w_am = _w2(Za, Zm, rng)
            # Mixture FVE-like: 1 - W2(SFG, mid)^2 / mean W2 to flanks^2
            fve = 1 - (w_mid ** 2) / (((w_aca ** 2 + w_mca ** 2) / 2) + 1e-12)
            h2 = {
                "W2_SFG_ACA": w_aca,
                "W2_SFG_MCA": w_mca,
                "W2_SFG_geodesic_mid": w_mid,
                "W2_ACA_MCA": w_am,
                "OT_mixture_FVE": float(fve),
                "prior_H2_FVE": 0.94,
                "interpretation": (
                    "If OT_mixture_FVE high, SFG sits near ACA–MCA geodesic (mixture/borderzone), "
                    "consistent with H2; SFG is NOT clean ACA."
                ),
            }

    # Conformal territory labels: distance to territory barycenters in Z
    bary = {t: centers[t].mean(axis=0) for t in territories}
    # calibration split
    idx = np.arange(len(Z))
    rng.shuffle(idx)
    n_cal = int(0.3 * len(idx))
    cal, te = idx[:n_cal], idx[n_cal:]
    # nonconformity = distance to claimed class barycenter
    scores = {}
    for t in territories:
        d = np.linalg.norm(Z[cal][labels[cal] == t] - bary[t], axis=1) if (labels[cal] == t).any() else np.array([np.inf])
        scores[t] = float(np.quantile(d, 0.9))
    # predict
    conf_rows = []
    correct = 0
    covered = 0
    n_te = 0
    for i in te:
        if labels[i] not in territories:
            continue
        n_te += 1
        dists = {t: float(np.linalg.norm(Z[i] - bary[t])) for t in territories}
        pred = min(dists, key=dists.get)
        # prediction set: all t with dist <= score_t
        pset = [t for t in territories if dists[t] <= scores[t]]
        if not pset:
            pset = [pred]
        conf_rows.append({"true": labels[i], "pred": pred, "set": "|".join(pset), "set_size": len(pset)})
        correct += int(pred == labels[i])
        covered += int(labels[i] in pset)
    conformal = {
        "accuracy": correct / max(n_te, 1),
        "set_coverage": covered / max(n_te, 1),
        "mean_set_size": float(np.mean([r["set_size"] for r in conf_rows])) if conf_rows else None,
        "n_test": n_te,
    }
    pd.DataFrame(conf_rows).to_csv(OUT / "conformal_territory_labels.csv", index=False)

    results = {
        "honesty_lock": HONESTY["5"],
        "territories": territories,
        "pairwise_W2": dist_df.to_dict(orient="records"),
        "H2_borderzone_OT": h2,
        "conformal": conformal,
    }
    (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUT / "report.md").write_text(
        "\n".join(
            [
                "# Moonshot 05 — Wasserstein territorial manifold",
                "",
                f"**Honesty lock:** {HONESTY['5']}",
                "",
                "## Pairwise W2",
                dist_df.to_markdown(index=False),
                "",
                f"## H2 OT borderzone\n```json\n{json.dumps(h2, indent=2)}\n```",
                "",
                f"## Conformal labels\n```json\n{json.dumps(conformal, indent=2)}\n```",
            ]
        ),
        encoding="utf-8",
    )
    print("DONE 05", results)
    return results


if __name__ == "__main__":
    run()
