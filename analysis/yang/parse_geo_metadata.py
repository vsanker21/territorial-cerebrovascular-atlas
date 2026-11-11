"""Parse GSE163577 series matrix into clean sample metadata + territory stubs."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
IN_PATH = ROOT / "data" / "raw" / "yang" / "GSE163577_series_matrix.txt.gz"
OUT_DIR = ROOT / "data" / "processed" / "yang"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    rows: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []
    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = line.split("\t")
        key = parts[0][len("!Sample_") :]
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key == "characteristics_ch1":
            char_rows.append(vals)
        else:
            rows[key] = vals
    n = len(next(iter(rows.values())))
    df = pd.DataFrame({k: v for k, v in rows.items() if len(v) == n})
    # Expand characteristics into columns when patterned as "key: value"
    for i, vals in enumerate(char_rows):
        if len(vals) != n:
            continue
        keys = []
        parsed = []
        for v in vals:
            m = re.match(r"\s*([^:]+):\s*(.*)\s*$", v)
            if m:
                keys.append(m.group(1).strip().lower().replace(" ", "_"))
                parsed.append(m.group(2).strip())
            else:
                keys.append(f"characteristic_{i}")
                parsed.append(v)
        # If all keys identical, use that name
        if len(set(keys)) == 1:
            df[keys[0]] = parsed
        else:
            df[f"characteristic_{i}"] = vals
    return df


def infer_region(title: str, source: str = "", chars: str = "") -> str:
    blob = f"{title} {source} {chars}".lower()
    if "hpc" in blob or "hippocamp" in blob:
        return "hippocampus"
    if "ctx" in blob or "cortex" in blob or "frontal" in blob or "sfc" in blob:
        return "superior_frontal_cortex"
    return "unknown"


def infer_disease(title: str) -> str:
    t = title.lower()
    if t.startswith("ad") or " alzheimer" in t or re.search(r"\bad\b", t):
        return "AD"
    if "control" in t or "nci" in t or "healthy" in t:
        return "control"
    return "unknown"


def main():
    df = parse_series_matrix(IN_PATH)
    title_col = "title" if "title" in df.columns else df.columns[0]
    source_col = "source_name_ch1" if "source_name_ch1" in df.columns else ""
    df["region"] = [
        infer_region(t, df[source_col].iloc[i] if source_col else "", "")
        for i, t in enumerate(df[title_col])
    ]
    df["disease"] = [infer_disease(t) for t in df[title_col]]
    # donor id from title patterns like "Control 1 Hpc" / "AD 1 Ctx"
    donor = []
    for t in df[title_col]:
        m = re.search(r"(Control|AD)\s*(\d+)", t, re.I)
        if not m:
            donor.append("unknown")
        else:
            group = "AD" if m.group(1).upper() == "AD" else "Control"
            donor.append(f"{group}_{m.group(2)}")
    df["donor_id"] = donor

    out = OUT_DIR / "yang_sample_metadata.csv"
    df.to_csv(out, index=False)
    summary = (
        df.groupby(["region", "disease"], dropna=False)
        .size()
        .reset_index(name="n_samples")
        .sort_values(["region", "disease"])
    )
    summary.to_csv(OUT_DIR / "yang_sample_summary.csv", index=False)
    print(df[[title_col, "geo_accession", "region", "disease", "donor_id"]].to_string(index=False))
    print("\nSummary:")
    print(summary.to_string(index=False))
    print("Wrote", out)


if __name__ == "__main__":
    main()
