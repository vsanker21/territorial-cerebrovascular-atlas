"""Inspect Liu arterial atlas geometry and label inventory."""
from pathlib import Path
import numpy as np
import nibabel as nib

ROOT = Path(__file__).resolve().parents[2]
ATLAS_DIR = ROOT / "data" / "raw" / "arterial_atlas" / "Atlas"

for name in ("ArterialAtlas.nii", "ArterialAtlas_level2.nii"):
    img = nib.load(str(ATLAS_DIR / name))
    data = np.asanyarray(img.dataobj)
    labels = np.unique(data.astype(int))
    print(f"=== {name} ===")
    print("shape:", img.shape)
    print("zooms:", img.header.get_zooms())
    print("affine:\n", img.affine)
    print("n_labels:", len(labels))
    print("labels:", labels.tolist())
    print()
