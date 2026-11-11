"""Build ADNI + ATLAS Aim4 figure panels from frozen tables."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures/open_agenda"
ADNI = ROOT / "results/aim4/adni"
ATLAS = ROOT / "results/aim4/atlas"
ISLES = ROOT / "results/aim4/isles"


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 200,
            "savefig.dpi": 200,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )


def save(fig, stem: str):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


def panel_adni():
    nest = pd.read_csv(ADNI / "aim4_nest_wmh_mra_subject_table.csv")
    summary = json.loads((ADNI / "aim4_nest_wmh_mra_summary.json").read_text(encoding="utf-8"))
    terr = ["ACA", "MCA", "PCA", "VB"]
    means = [float(nest[f"wmh_frac_of_wmh_{t}"].mean()) for t in terr]
    sems = [float(nest[f"wmh_frac_of_wmh_{t}"].sem()) for t in terr]
    ols = summary["ols_r2_TOTAL_WMH"]
    r2s = [
        float(ols["M0_demographics"]["r2"]),
        float(ols["M1_mra_geometry_null"]["r2"]),
        float(ols["M3_demographics_plus_tv_burden"]["r2"]),
    ]
    labels = ["M0 demog", "M1 MRA null", "M3 + TV WMH"]

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    ax = axes[0]
    x = np.arange(len(terr))
    ax.bar(x, means, yerr=sems, color=["#2a6f97", "#c44536", "#3d5a40", "#6c757d"], width=0.7, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(terr)
    ax.set_ylabel("Mean WMH share (of local WMH)")
    ax.set_title(f"ADNI primary QC n={len(nest)}")
    ax.set_ylim(0, max(means) * 1.35)

    ax = axes[1]
    ax.bar(np.arange(len(r2s)), r2s, color=["#495057", "#adb5bd", "#2a6f97"], width=0.65)
    ax.set_xticks(np.arange(len(r2s)))
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("Nested OLS R² (TOTAL_WMH)")
    ax.set_title("Geography nest (not molecular H3)")
    ax.set_ylim(0, 1.05)
    fig.suptitle("Aim 4 non-stroke — ADNI WMH × TerritoryVoxels", fontsize=11, y=1.02)
    fig.tight_layout()
    save(fig, "06_adni_aim4_wmh_nest")


def panel_atlas():
    nest = pd.read_csv(ATLAS / "aim4_nest_atlas_mra_subject_table.csv")
    summary = json.loads((ATLAS / "aim4_nest_atlas_mra_summary.json").read_text(encoding="utf-8"))
    terr = ["ACA", "MCA", "PCA", "VB"]
    means = [float(nest[f"lesion_frac_of_lesion_{t}"].mean()) for t in terr]
    sems = [float(nest[f"lesion_frac_of_lesion_{t}"].sem()) for t in terr]
    vc = nest["dominant_territory"].astype(str).value_counts().reindex(terr, fill_value=0)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    ax = axes[0]
    x = np.arange(len(terr))
    ax.bar(x, means, yerr=sems, color=["#2a6f97", "#c44536", "#3d5a40", "#6c757d"], width=0.7, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(terr)
    ax.set_ylabel("Mean lesion share")
    ax.set_title(f"ATLAS primary QC n={len(nest)}")
    ax.set_ylim(0, max(means) * 1.25)

    ax = axes[1]
    ax.bar(np.arange(len(terr)), vc.values, color=["#2a6f97", "#c44536", "#3d5a40", "#6c757d"], width=0.7)
    ax.set_xticks(np.arange(len(terr)))
    ax.set_xticklabels(terr)
    ax.set_ylabel("n subjects (dominant territory)")
    ax.set_title(f"Dominant major (MCA={summary['dominant_territory_counts'].get('MCA', int(vc['MCA']))})")
    fig.suptitle("Aim 4 stroke — ATLAS × ArterialAtlas_level2", fontsize=11, y=1.02)
    fig.tight_layout()
    save(fig, "07_atlas_aim4_lesion_nest")


def panel_isles():
    nest = pd.read_csv(ISLES / "aim4_nest_isles_mra_subject_table.csv")
    summary = json.loads((ISLES / "aim4_nest_isles_mra_summary.json").read_text(encoding="utf-8"))
    terr = ["ACA", "MCA", "PCA", "VB"]
    means = [float(nest[f"lesion_frac_of_lesion_{t}"].mean()) for t in terr]
    sems = [float(nest[f"lesion_frac_of_lesion_{t}"].sem()) for t in terr]
    vc = nest["dominant_territory"].astype(str).value_counts().reindex(terr, fill_value=0)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    ax = axes[0]
    x = np.arange(len(terr))
    ax.bar(x, means, yerr=sems, color=["#2a6f97", "#c44536", "#3d5a40", "#6c757d"], width=0.7, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(terr)
    ax.set_ylabel("Mean lesion share")
    ax.set_title(f"ISLES'22 primary QC n={len(nest)}")
    ax.set_ylim(0, max(means) * 1.25)

    ax = axes[1]
    ax.bar(np.arange(len(terr)), vc.values, color=["#2a6f97", "#c44536", "#3d5a40", "#6c757d"], width=0.7)
    ax.set_xticks(np.arange(len(terr)))
    ax.set_xticklabels(terr)
    ax.set_ylabel("n subjects (dominant territory)")
    ax.set_title(f"Dominant major (MCA={summary['dominant_territory_counts'].get('MCA', int(vc['MCA']))})")
    fig.suptitle("Aim 4 stroke replication — ISLES'22 × ArterialAtlas_level2", fontsize=11, y=1.02)
    fig.tight_layout()
    save(fig, "08_isles_aim4_lesion_nest")


def main():
    style()
    panel_adni()
    panel_atlas()
    panel_isles()


if __name__ == "__main__":
    main()
