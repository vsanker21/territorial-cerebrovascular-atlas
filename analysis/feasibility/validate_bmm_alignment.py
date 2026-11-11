"""Sanity-check BMM channel order and MNI alignment against hard atlas + anatomy."""
from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from nilearn import datasets, image

ROOT = Path(__file__).resolve().parents[2]
BMM_DIR = ROOT / "data/raw/arterial_atlas/nitrc/extracted/Atlas_MNI152/Atlas_182_MNI152"
HARD_DIR = ROOT / "data/raw/arterial_atlas/Atlas_182"

ORDER = ["ACA", "MCA", "PCA", "VB"]
# Hard level2 intensities → major (from labels file)
HARD_TO_MAJOR = {1: "ACA", 2: "ACA", 3: "MCA", 4: "MCA", 5: "PCA", 6: "PCA", 7: "VB", 8: "VB"}


def load3(path):
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj)
    if data.ndim == 4 and data.shape[-1] == 1:
        data = data[..., 0]
    return img, data


def main():
    bmm_img, bmm = load3(BMM_DIR / "ProbArterialAtlas_BMM.nii")
    hard_img = nib.load(str(HARD_DIR / "ArterialAtlas_level2.nii"))
    hard = np.asanyarray(hard_img.dataobj)
    if hard.ndim == 4:
        hard = hard[..., 0]
    # Restore FSL affine on GitHub hard map for comparison space
    hard_img = nib.Nifti1Image(hard.astype(np.int16), bmm_img.affine)

    terr_img, terr = load3(BMM_DIR / "TerritoryVoxels_BMM.nii")
    border_img, border = load3(BMM_DIR / "BorderZone_ProbAve.nii")

    print("BMM affine:\n", bmm_img.affine)
    print("BMM shape", bmm.shape, "dtype", bmm.dtype, "min/max", float(np.nanmin(bmm)), float(np.nanmax(bmm)))
    print("Border nan fraction", float(np.isnan(border).mean()), "min/max", float(np.nanmin(border)), float(np.nanmax(border)))
    print("TerritoryVoxels unique", np.unique(terr))

    # Channel dominance vs hard major label agreement
    # Restrict to hard-labeled brain
    hard_rs = image.resample_to_img(
        hard_img, bmm_img, interpolation="nearest", force_resample=True, copy_header=True
    )
    hard_a = np.asanyarray(hard_rs.dataobj).astype(int)
    if hard_a.ndim == 4:
        hard_a = hard_a[..., 0]
    mask = (hard_a > 0) & (hard_a < 9) & (bmm.sum(-1) > 0)
    dom = bmm[mask].argmax(axis=-1)
    hard_maj = np.array([HARD_TO_MAJOR.get(int(v), "X") for v in hard_a[mask]])
    bmm_maj = np.array([ORDER[i] for i in dom])
    print("\nAgreement hard-major vs BMM-argmax (same space):")
    for maj in ORDER:
        sel = hard_maj == maj
        if sel.sum() == 0:
            continue
        agree = (bmm_maj[sel] == maj).mean()
        # distribution of BMM argmax within hard maj
        dist = {ORDER[i]: float((dom[sel] == i).mean()) for i in range(4)}
        print(f"  hard={maj}: n={sel.sum()}, agree={agree:.3f}, bmm_dist={dist}")

    # Anatomical landmark probes (MNI mm)
    landmarks = {
        "medial_PFC_approx": (0, 50, 0),       # expect ACA
        "lateral_IFG_approx": (50, 20, 10),    # expect MCA
        "occipital_pole_approx": (0, -90, 0),  # expect PCA
        "hippocampus_L_approx": (-30, -20, -15),  # expect PCA
    }
    inv = np.linalg.inv(bmm_img.affine)
    print("\nLandmark BMM vectors (channel order assumed ACA,MCA,PCA,VB):")
    for name, xyz in landmarks.items():
        ijk = inv @ np.array([xyz[0], xyz[1], xyz[2], 1.0])
        i, j, k = [int(round(x)) for x in ijk[:3]]
        if min(i, j, k) < 0 or i >= bmm.shape[0] or j >= bmm.shape[1] or k >= bmm.shape[2]:
            print(f"  {name}: out of bounds ijk={(i,j,k)}")
            continue
        p = bmm[i, j, k]
        print(f"  {name} xyz={xyz} ijk={(i,j,k)} p={dict(zip(ORDER, map(float, p)))} argmax={ORDER[int(p.argmax())]} hard={HARD_TO_MAJOR.get(int(hard_a[i,j,k]), hard_a[i,j,k])}")

    # HO SFG mean under BMM vs hard
    cort = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    labels = list(cort.labels)
    data = np.asanyarray(cort.maps.dataobj)
    sfg = data == labels.index("Superior Frontal Gyrus")
    fmed = data == labels.index("Frontal Medial Cortex")
    for roi_name, mask2 in [("SFG", sfg), ("FrontalMedial", fmed)]:
        mimg = nib.Nifti1Image(mask2.astype(np.uint8), cort.maps.affine)
        mrs = image.resample_to_img(mimg, bmm_img, interpolation="nearest", force_resample=True, copy_header=True)
        m = np.asanyarray(mrs.dataobj).astype(bool) & (bmm.sum(-1) > 0)
        mean_p = bmm[m].mean(0)
        mean_p = mean_p / mean_p.sum()
        hard_counts = {}
        for v in hard_a[m]:
            maj = HARD_TO_MAJOR.get(int(v), "X")
            hard_counts[maj] = hard_counts.get(maj, 0) + 1
        tot = sum(hard_counts.values()) or 1
        hard_frac = {k: v / tot for k, v in hard_counts.items()}
        print(f"\n{roi_name}: BMM_mean={dict(zip(ORDER, map(float, mean_p)))} hard_frac={hard_frac} n={m.sum()}")


if __name__ == "__main__":
    main()
