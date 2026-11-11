"""Unzip NITRC Arterial Atlas packages and inventory probabilistic / borderzone files."""
from __future__ import annotations

import json
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
NITRC = ROOT / "data" / "raw" / "arterial_atlas" / "nitrc"
EXTRACTED = NITRC / "extracted"
OUT_JSON = ROOT / "results" / "feasibility" / "nitrc_ingest_inventory.json"

REQUIRED = [
    "ProbArterialAtlas_BMM.nii",
    "BorderZone_ProbAve.nii",
    "TerritoryVoxels_BMM.nii",
    "ArterialAtlas.nii",
    "ArterialAtlas_level2.nii",
]
OPTIONAL = [
    "ProbArterialAtlas_average.nii",
    "BorderZone_ProbAve2.nii",
]
WANTED = REQUIRED + OPTIONAL


def is_real_zip(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 1000:
        return False
    with path.open("rb") as f:
        magic = f.read(4)
    return magic.startswith(b"PK")


def main():
    EXTRACTED.mkdir(parents=True, exist_ok=True)
    inventory = {"zips": {}, "found": {}, "missing_wanted": [], "notes": []}

    for zname in ("ArterialAtlas.zip", "Atlas_MNI152.zip", "ArterialAtlas_zip.zip"):
        zpath = NITRC / zname
        inventory["zips"][zname] = {
            "exists": zpath.exists(),
            "size": zpath.stat().st_size if zpath.exists() else 0,
            "is_zip": is_real_zip(zpath) if zpath.exists() else False,
        }
        if not is_real_zip(zpath):
            inventory["notes"].append(
                f"{zname} missing or is not a ZIP (likely NITRC login HTML). See docs/05_NITRC_ingest_instructions.md"
            )
            continue
        with zipfile.ZipFile(zpath, "r") as zf:
            zf.extractall(EXTRACTED / zpath.stem)
            inventory["zips"][zname]["members"] = zf.namelist()

    # Inventory wanted files anywhere under nitrc/
    for want in WANTED:
        hits = list(NITRC.rglob(want)) + list(NITRC.rglob(want + ".gz"))
        if not hits:
            inventory["missing_wanted"].append(want)
            continue
        path = hits[0]
        info = {"path": str(path.relative_to(ROOT)), "size": path.stat().st_size}
        try:
            img = nib.load(str(path))
            data = np.asanyarray(img.dataobj)
            info.update({"shape": list(data.shape), "zooms": [float(x) for x in img.header.get_zooms()[:4]]})
        except Exception as e:  # noqa: BLE001
            info["load_error"] = str(e)
        inventory["found"][want] = info

    missing_required = [w for w in REQUIRED if w not in inventory["found"]]
    inventory["missing_required"] = missing_required
    inventory["missing_optional"] = [w for w in OPTIONAL if w not in inventory["found"]]
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(inventory, indent=2), encoding="utf-8")
    print(json.dumps(inventory, indent=2))
    if missing_required:
        raise SystemExit(
            "Required NITRC maps missing: "
            + ", ".join(missing_required)
            + ". See docs/05_NITRC_ingest_instructions.md"
        )
    print("NITRC ingest OK — required BMM/borderzone maps present.")


if __name__ == "__main__":
    main()
