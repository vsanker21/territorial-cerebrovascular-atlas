"""
Build a *labeled methodological proxy* for Liu probabilistic / borderzone maps.

Scientific honesty
------------------
NITRC hosts ProbArterialAtlas_BMM / BorderZone_ProbAve / TerritoryVoxels_BMM and
requires a logged-in account. Until those files are ingested, confirmatory Aim 1/3
analyses MUST NOT treat this proxy as BMM.

What this builds
----------------
From Liu hard level-2 labels (ACA/MCA/PCA/VB), construct a 4-channel soft
assignment by Gaussian-smoothed one-hot maps, renormalized to a simplex per voxel.
Borderzone fraction proxy = 1 - max_c p_c (high where territories compete).
MCA/ACA and MCA/PCA ratio channels mimic BorderZone_ProbAve structure.

These are geometric softenings of the hard atlas, NOT lesion-derived Bernoulli
mixture probabilities. File names and JSON sidecars mark them PROXY.
"""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[2]
ATLAS_DIR = ROOT / "data" / "raw" / "arterial_atlas" / "Atlas_182"
OUT_DIR = ROOT / "data" / "processed" / "arterial_proxy"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FSL_MNI152_1MM_AFFINE = np.array(
    [[-1.0, 0.0, 0.0, 90.0], [0.0, 1.0, 0.0, -126.0], [0.0, 0.0, 1.0, -72.0], [0.0, 0.0, 0.0, 1.0]]
)
MAJOR = {"ACA": (1, 2), "MCA": (3, 4), "PCA": (5, 6), "VB": (7, 8)}
ORDER = ["ACA", "MCA", "PCA", "VB"]
SIGMA_MM = 2.0  # isotropic smoothing in voxels at 1 mm


def load_level2() -> tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(ATLAS_DIR / "ArterialAtlas_level2.nii"))
    data = np.asanyarray(img.dataobj)
    if data.ndim == 4:
        data = data[..., 0]
    return data.astype(np.int16), FSL_MNI152_1MM_AFFINE


def main():
    hard, affine = load_level2()
    brain = (hard > 0) & (hard < 9)
    onehots = []
    for name in ORDER:
        m = np.isin(hard, MAJOR[name]).astype(np.float32)
        soft = gaussian_filter(m, sigma=SIGMA_MM)
        onehots.append(soft)
    stack = np.stack(onehots, axis=-1)  # X,Y,Z,4
    # Zero outside brain support of hard atlas before renormalize
    stack[~brain] = 0
    denom = stack.sum(axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        probs = np.divide(stack, denom, out=np.zeros_like(stack), where=denom > 0)

    pmax = probs.max(axis=-1)
    border_frac = np.where(brain, 1.0 - pmax, 0.0)

    eps = 1e-6
    mca = probs[..., 1]
    aca = probs[..., 0]
    pca = probs[..., 2]
    ratio = np.stack(
        [
            np.where(brain, mca / (aca + eps), 0.0),
            np.where(brain, mca / (pca + eps), 0.0),
        ],
        axis=-1,
    )

    meta = {
        "type": "PROXY_NOT_BMM",
        "warning": (
            "Geometric soft assignment from hard Liu labels via Gaussian smoothing. "
            "Replace with NITRC ProbArterialAtlas_BMM.nii and BorderZone_ProbAve.nii "
            "before confirmatory analyses or pre-registration lock."
        ),
        "sigma_mm": SIGMA_MM,
        "channels_prob": ORDER,
        "channels_ratio": ["MCA_over_ACA", "MCA_over_PCA"],
        "nitrc_required": [
            "https://www.nitrc.org/frs/download.php/12271/ArterialAtlas.zip",
            "https://www.nitrc.org/frs/download.php/12382/Atlas_MNI152.zip",
            "https://www.nitrc.org/frs/download.php/12244/Vascular_Probabilistic_Maps.nii",
        ],
    }

    nib.save(nib.Nifti1Image(probs.astype(np.float32), affine), str(OUT_DIR / "PROXY_soft_territory_probs.nii.gz"))
    nib.save(nib.Nifti1Image(border_frac.astype(np.float32), affine), str(OUT_DIR / "PROXY_borderzone_fraction.nii.gz"))
    nib.save(nib.Nifti1Image(ratio.astype(np.float32), affine), str(OUT_DIR / "PROXY_borderzone_ratios.nii.gz"))
    (OUT_DIR / "PROXY_README.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("Wrote proxy maps to", OUT_DIR)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
