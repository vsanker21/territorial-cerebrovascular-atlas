from pathlib import Path
from collections import Counter
import zipfile
import re

root = Path(r"F:/Cerebrovascular SCS/data/raw/ADNI")
files = [p for p in root.rglob("*") if p.is_file()]
print(f"files={len(files)} size_GB={sum(f.stat().st_size for f in files)/1e9:.2f}")
print()
for p in sorted(root.iterdir(), key=lambda x: x.name.lower()):
    if p.is_file():
        print(f"FILE {p.stat().st_size/1e9:8.3f} GB  {p.name}")
    else:
        sub = [f for f in p.rglob("*") if f.is_file()]
        sz = sum(f.stat().st_size for f in sub)
        print(f"DIR  {sz/1e9:8.3f} GB  {p.name}/  ({len(sub)} files)")

# classify by name keywords
print("\n=== keyword hits among filenames ===")
keys = ["flair", "wmh", "mask", "mprage", "adnimerge", "ucd", "nifti", "nii"]
for k in keys:
    hits = [p for p in files if k in p.name.lower()]
    if hits:
        print(k, len(hits), "e.g.", hits[0].name)

# sample zip contents for flair
print("\n=== zip/tar peeks for FLAIR ===")
for p in sorted(root.glob("*.zip")) + sorted(root.glob("*.tar*")):
    name = p.name.lower()
    if p.stat().st_size < 1e6 and "metadata" in name:
        continue
    try:
        if p.suffix.lower() == ".zip" or p.name.endswith(".zip"):
            with zipfile.ZipFile(p) as zf:
                names = zf.namelist()
            fl = sum(1 for n in names if "flair" in n.lower())
            t1 = sum(1 for n in names if re.search(r"mprage|mp-rage|spgr", n, re.I))
            print(f"{p.name}: members={len(names)} flair_name_hits={fl} t1_name_hits={t1} size_GB={p.stat().st_size/1e9:.2f}")
    except Exception as e:
        print(p.name, "ERR", e)
