"""Resolve and download Allen Cortical ACC + probe assets."""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

ROOT = Path(r"F:\Cerebrovascular SCS")
OUT = ROOT / "data/raw/aca_cohort_hunt/allen_cortical_acc"
OUT.mkdir(parents=True, exist_ok=True)

CID = "d17249d2-0e6e-4500-abb8-e6c93fa1ac6f"
WANT = {
    "be401db3-d732-408a-b0c4-71af0458b8ab": "Dissection_ACC",
    "5346f9c6-755e-4336-94cc-38706ec00c2f": "Supercluster_NonNeuronal",
}


def get(url: str):
    req = urllib.request.Request(url, headers={"User-Agent": "aca-hunt"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def main():
    c = get(f"https://api.cellxgene.cziscience.com/curation/v1/collections/{CID}")
    (OUT / "collection.json").write_text(json.dumps(c, indent=2), encoding="utf-8")
    urls = {}
    for d in c.get("datasets", []):
        did = d.get("dataset_id")
        if did not in WANT:
            continue
        print("DATASET", did, d.get("title"), "keys", list(d.keys())[:30])
        assets = d.get("assets") or []
        print(" assets", assets)
        for a in assets:
            if a.get("filetype") == "H5AD" or str(a.get("filename", "")).endswith(".h5ad"):
                urls[did] = a.get("url") or a.get("s3_uri")
                print(" H5AD", urls[did])
    # fallback patterns
    for did, name in WANT.items():
        if did in urls and urls[did]:
            continue
        for u in [
            f"https://datasets.cellxgene.cziscience.com/{did}.h5ad",
            f"https://cellxgene-data-public.s3.amazonaws.com/cell-census/2024-07-01/h5ads/{did}.h5ad",
            f"https://cellxgene-data-public.s3.amazonaws.com/cell-census/2023-12-15/h5ads/{did}.h5ad",
            f"https://cellxgene-data-public.s3.amazonaws.com/cell-census/2023-07-25/h5ads/{did}.h5ad",
        ]:
            try:
                req = urllib.request.Request(u, method="HEAD", headers={"User-Agent": "aca"})
                with urllib.request.urlopen(req, timeout=30) as r:
                    print("OK HEAD", u, r.status, r.headers.get("Content-Length"))
                    urls[did] = u
                    break
            except Exception as e:
                print("FAIL", u, type(e).__name__)
    (OUT / "download_urls.json").write_text(json.dumps(urls, indent=2), encoding="utf-8")
    print("URLS", urls)


if __name__ == "__main__":
    main()
