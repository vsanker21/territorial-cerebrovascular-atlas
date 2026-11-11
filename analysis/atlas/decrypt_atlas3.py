"""
Decrypt ATLAS R3.0 preprocessed package (OpenSSL Salted__, base64).

Official NITRC command:
  openssl aes-256-cbc -md sha256 -d -a -in atlas3_training_preprocessed.tar.gz \\
    -out ATLAS_R3.0_preprocessed.tar.gz

Passphrase: from ATLAS Google-form email / local key file (do NOT commit or share).

Usage:
  set ATLAS_PASSPHRASE=...
  python analysis/atlas/decrypt_atlas3.py

  # or read from key file without echoing:
  python analysis/atlas/decrypt_atlas3.py --key-file "data/raw/ATLAS Encryption Key.txt"
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENC = ROOT / "data/raw/atlas3_training_preprocessed.tar.gz"
OUT_TAR_GZ = ROOT / "data/raw/atlas_stroke/ATLAS_R3.0_preprocessed.tar.gz"
DEFAULT_KEY_FILE = ROOT / "data/raw/ATLAS Encryption Key.txt"


def is_openssl_salted(path: Path) -> bool:
    with path.open("rb") as f:
        head = f.read(16)
    return head.startswith(b"Salted__") or head.startswith(b"U2FsdGVk")


def parse_key_file(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"encryption key is:\s*(\S+)", text, flags=re.IGNORECASE)
    if not m:
        raise SystemExit(f"Could not parse encryption key from {path}")
    return m.group(1).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--passphrase", default=os.environ.get("ATLAS_PASSPHRASE", ""))
    ap.add_argument("--key-file", default="", help="ATLAS Encryption Key.txt (preferred)")
    ap.add_argument("--in", dest="inp", default=str(ENC))
    ap.add_argument("--out", dest="out", default=str(OUT_TAR_GZ))
    args = ap.parse_args()
    inp = Path(args.inp)
    out = Path(args.out)
    if not inp.exists():
        raise SystemExit(f"Missing {inp}")
    if not is_openssl_salted(inp):
        print("File does not look OpenSSL-Salted/base64; try opening as tar.gz directly.")
        return

    passphrase = args.passphrase
    key_file = Path(args.key_file) if args.key_file else DEFAULT_KEY_FILE
    if not passphrase and key_file.exists():
        passphrase = parse_key_file(key_file)
        print(f"Loaded passphrase from {key_file.name} (not printed).")
    if not passphrase:
        raise SystemExit(
            "Provide passphrase via --passphrase, ATLAS_PASSPHRASE, or --key-file "
            f"(expected {DEFAULT_KEY_FILE})."
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    pass_file = out.parent / "_pass.tmp"
    try:
        pass_file.write_text(passphrase, encoding="utf-8")
        # Official ATLAS/NITRC flags: aes-256-cbc -md sha256 -d -a
        cmd = [
            "openssl",
            "aes-256-cbc",
            "-md",
            "sha256",
            "-d",
            "-a",
            "-in",
            str(inp),
            "-out",
            str(out),
            "-pass",
            f"file:{pass_file}",
        ]
        print("Running openssl aes-256-cbc -md sha256 -d -a ...")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr)
            raise SystemExit("Decryption failed — check passphrase / cipher.")
    finally:
        if pass_file.exists():
            pass_file.unlink()

    print(f"Wrote {out} ({out.stat().st_size} bytes)")
    with out.open("rb") as f:
        magic = f.read(3)
    print("output magic", magic)
    if magic != b"\x1f\x8b\x08":
        print("WARNING: expected gzip magic 1f8b08; verify before extract.")


if __name__ == "__main__":
    main()
