"""Quick anatomy audit of GEO series matrices + on-disk Walchli-like paths."""
from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HUNT = ROOT / "data/raw/aca_cohort_hunt"

PATS = [
    "cingulate",
    "anterior",
    "frontal",
    "Brodmann",
    "prefrontal",
    "temporal",
    "hippocamp",
    "region",
    "tissue",
    "source_name",
    "characteristics",
]


def main():
    for acc in ["GSE222007", "GSE252921", "GSE226606"]:
        p = HUNT / acc / f"{acc}_series_matrix.txt.gz"
        if not p.exists():
            print(acc, "no matrix")
            continue
        txt = gzip.open(p, "rt", errors="replace").read()
        print("====", acc, "====")
        for pat in PATS:
            hits = [line for line in txt.splitlines() if pat.lower() in line.lower()]
            if not hits:
                continue
            print(pat, "n=", len(hits))
            for h in hits[:8]:
                print(" ", h[:220])

    print("==== WALCHLI / VASC PATHS ====")
    raw = ROOT / "data/raw"
    for p in raw.rglob("*"):
        name = p.name.lower()
        if any(k in name for k in ["walch", "waelch", "gse256493", "winkler", "crouch"]):
            print(p)


if __name__ == "__main__":
    main()
