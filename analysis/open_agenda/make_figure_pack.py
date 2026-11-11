"""Figure pack for open-data Class B / H1 status (matplotlib only; paper polish)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures" / "open_agenda"
FIG.mkdir(parents=True, exist_ok=True)

# Restrained scientific palette (avoid purple/glow AI defaults)
C_ACA = "#1B4F72"
C_MCA = "#B9770E"
C_NULL = "#5D6D7E"
C_BAR = "#2E86AB"
C_WARN = "#922B21"


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fig_tv_priors():
    paths = [
        ROOT / "results/aca_cohort_hunt/linnarsson_ft_tv_priors.json",
        ROOT / "results/aca_cohort_hunt/linnarsson_aca_tv_priors.json",
    ]
    p = next((x for x in paths if x.exists()), None)
    if p is None:
        return
    d = json.loads(p.read_text(encoding="utf-8"))
    labels, p_aca, p_mca = [], [], []
    for k, v in d.items():
        labels.append(k)
        p_aca.append(float(v.get("tv_p_ACA", 0) or 0))
        p_mca.append(float(v.get("tv_p_MCA", 0) or 0))
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.bar(x - 0.18, p_aca, 0.36, label="p(ACA)", color=C_ACA)
    ax.bar(x + 0.18, p_mca, 0.36, label="p(MCA)", color=C_MCA)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("TerritoryVoxels fraction")
    ax.set_title("TerritoryVoxels priors — Siletti dissections (HO ROIs)")
    ax.axhline(0.5, color=C_NULL, ls="--", lw=0.8)
    ax.legend(frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "01_territoryvoxels_linnarsson_priors.png", dpi=200)
    fig.savefig(FIG / "01_territoryvoxels_linnarsson_priors.pdf")
    plt.close(fig)


def fig_class_b_umap():
    h5 = ROOT / "data/processed/class_b/aca_mca_harmony_endo.h5ad"
    if not h5.exists():
        return
    a = sc.read_h5ad(h5)
    if "X_umap" not in a.obsm:
        use = "X_pca_harmony" if "X_pca_harmony" in a.obsm else "X_pca"
        sc.pp.neighbors(a, use_rep=use)
        sc.tl.umap(a)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    color_maps = {
        "class_b_arm": { "ACA": C_ACA, "MCA": C_MCA },
    }
    for ax, key, title in [
        (axes[0], "class_b_arm", "Class B arm"),
        (axes[1], "batch", "Batch / source"),
    ]:
        cats = a.obs[key].astype(str)
        uniq = sorted(cats.unique())
        cmap = color_maps.get(key, {})
        for i, c in enumerate(uniq):
            m = cats == c
            ax.scatter(
                a.obsm["X_umap"][m, 0],
                a.obsm["X_umap"][m, 1],
                s=10,
                alpha=0.75,
                color=cmap.get(c, plt.cm.tab10(i % 10)),
                label=str(c),
                linewidths=0,
            )
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(markerscale=1.8, fontsize=8, frameon=False, loc="best")
        _style(ax)
    fig.suptitle("ACA vs MCA Harmony UMAP (Endo)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "02_class_b_umap.png", dpi=200)
    fig.savefig(FIG / "02_class_b_umap.pdf")
    plt.close(fig)


def fig_sensitivity():
    p = ROOT / "results/class_b/class_b_hardening_summary.json"
    alt = ROOT / "results/class_b/aca_mca_sensitivity_summary.json"
    rows = []
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        # flatten known keys
        mapping = [
            ("baseline", "donor_oof_baseline", "Baseline"),
            ("ambient_residualized", "donor_oof_ambient_residualized", "Ambient residualized"),
            ("modules_only", "donor_oof_modules_only", "Modules only"),
        ]
        for _, k, lab in mapping:
            if k in d and d[k] is not None:
                rows.append((lab, float(d[k])))
        # nested
        for k, v in d.items():
            if isinstance(v, dict) and "donor_oof" in v:
                rows.append((k.replace("_", " "), float(v["donor_oof"])))
    elif alt.exists():
        d = json.loads(alt.read_text(encoding="utf-8"))
        for k, v in d.items():
            if isinstance(v, dict) and "donor_oof_auroc" in v:
                rows.append((k, float(v["donor_oof_auroc"])))
    # dedupe keeping first
    seen = set()
    uniq = []
    for lab, val in rows:
        if lab in seen or val != val:
            continue
        seen.add(lab)
        uniq.append((lab, val))
    if not uniq:
        return
    uniq = sorted(uniq, key=lambda x: x[1])
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.barh([r[0] for r in uniq], [r[1] for r in uniq], color=C_BAR)
    ax.axvline(0.7, color=C_WARN, ls="--", lw=1.0, label="0.70 supportive")
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Donor OOF AUROC")
    ax.set_title("Class B sensitivity — donor leave-one-out AUROC")
    ax.legend(frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "03_class_b_sensitivity_donor_oof.png", dpi=200)
    fig.savefig(FIG / "03_class_b_sensitivity_donor_oof.pdf")
    plt.close(fig)


def fig_ambient():
    scores = ROOT / "results/class_b/class_b_module_scores_by_donor.csv"
    if not scores.exists():
        return
    df = pd.read_csv(scores)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    colors = {"ACA": C_ACA, "MCA": C_MCA}
    for ax, col, title in [
        (axes[0], "ambient_neuronal", "Ambient neuronal score"),
        (axes[1], "ambient_oligo", "Ambient oligo score"),
    ]:
        if col not in df.columns:
            continue
        for arm, g in df.groupby("arm"):
            ax.hist(g[col], bins=12, alpha=0.55, label=str(arm), color=colors.get(str(arm), C_NULL))
        ax.set_title(title)
        ax.set_xlabel("Donor mean score")
        ax.legend(frameon=False)
        _style(ax)
    fig.suptitle("Class B ambient contamination (donor means)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "04_class_b_ambient_panel.png", dpi=200)
    fig.savefig(FIG / "04_class_b_ambient_panel.pdf")
    plt.close(fig)


def fig_linnarsson_ft():
    p = ROOT / "results/h1/linnarsson_within_study_summary.json"
    if not p.exists():
        return
    d = json.loads(p.read_text(encoding="utf-8"))
    contrasts = d.get("contrasts") or []
    if not contrasts:
        return
    labels = [c["contrast"].replace("_", "\n") for c in contrasts]
    aucs = []
    for c in contrasts:
        a = c.get("donor_lodo_pooled_cell_auroc", c.get("donor_oof_auroc"))
        aucs.append(float(a) if a is not None and a == a else np.nan)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cols = [C_ACA if "ACA" in c["contrast"] else C_NULL for c in contrasts]
    ax.bar(range(len(labels)), aucs, color=cols)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.7, color=C_WARN, ls="--", lw=1.0)
    ax.set_ylabel("Donor OOF AUROC")
    ax.set_title("Within-Siletti Endo contrasts (TerritoryVoxels arms)")
    _style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "05_linnarsson_within_study_oof.png", dpi=200)
    fig.savefig(FIG / "05_linnarsson_within_study_oof.pdf")
    plt.close(fig)


def main():
    fig_tv_priors()
    fig_class_b_umap()
    fig_sensitivity()
    fig_ambient()
    fig_linnarsson_ft()
    (FIG / "README.md").write_text(
        "# Open-agenda figure pack (paper)\n\n"
        "| File | Content |\n|---|---|\n"
        "| 01_territoryvoxels_linnarsson_priors | TerritoryVoxels p(ACA)/p(MCA) for Siletti dissections |\n"
        "| 02_class_b_umap | Class B ACA vs MCA Harmony UMAP |\n"
        "| 03_class_b_sensitivity_donor_oof | Donor OOF AUROC sensitivity |\n"
        "| 04_class_b_ambient_panel | Ambient neuronal/oligo by arm |\n"
        "| 05_linnarsson_within_study_oof | Within-Siletti F/T/ACA donor OOF |\n\n"
        "PNG (200 dpi) + PDF for each panel. Claims text: `docs/10_class_b_claims_paragraph.md`.\n",
        encoding="utf-8",
    )
    print("figures in", FIG)


if __name__ == "__main__":
    main()
