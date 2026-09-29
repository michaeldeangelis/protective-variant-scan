#!/usr/bin/env python
"""Download gnomAD v4.1 constraint metrics (anonymous HTTPS, ~91 MiB) and write <data>/constraint.csv."""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from protscan.adapters import common, gnomad  # noqa: E402

README_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/README.txt"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(ROOT / "data"))
    ap.add_argument("--force", action="store_true", help="re-download even if the file exists")
    a = ap.parse_args()
    data = Path(a.data)
    raw = data / "raw" / "gnomad"
    manifest = data / "raw" / "MANIFEST.tsv"

    for url in (gnomad.URL, README_URL):
        dest = raw / url.rsplit("/", 1)[1]
        if a.force or not dest.exists():
            rec = common.download(url, dest)
            rec["path"] = str(dest.relative_to(data))
            common.record_manifest(manifest, rec)
            print(f"downloaded {dest.name}: {rec['bytes']} bytes sha256={rec['sha256'][:16]}...")

    df = pd.read_csv(raw / gnomad.URL.rsplit("/", 1)[1], sep="\t", low_memory=False, na_values=["NA"])
    out = gnomad.normalize_constraint(df)
    out.to_csv(data / "constraint.csv", index=False)
    print(f"wrote {data / 'constraint.csv'}: {len(out)} genes")


if __name__ == "__main__":
    main()
