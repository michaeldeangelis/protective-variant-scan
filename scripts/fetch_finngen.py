#!/usr/bin/env python
"""Download FinnGen R13 LoF burden results (public bucket, anonymous HTTPS, ~550 MB) and write <data>/burden_finngen.csv.gz.

The gz is streamed and filtered to the declared proxy and trade-off endpoints; the full table is kept only as the raw download.
FinnGen asks users to submit its online form (docs/data-sources.md); this script does not submit it.
"""
import argparse
import sys
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402
import requests  # noqa: E402

from protscan.adapters import common, finngen  # noqa: E402

BUCKET = "finngen-public-data-r13"
OBJECTS = [
    "lof/finngen_R13_lof_readme",
    "lof/data/finngen_R13_lof_variants.txt",
    "summary_stats/finngen_R13_manifest.tsv",
    "lof/data/finngen_R13_lof.txt.gz",
]
LOF = "lof/data/finngen_R13_lof.txt.gz"


def gcs_md5(name: str) -> str:
    r = requests.get(f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o/{quote(name, safe='')}", timeout=60)
    r.raise_for_status()
    return r.json()["md5Hash"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(ROOT / "data"))
    ap.add_argument("--force", action="store_true", help="re-download even if the file exists")
    a = ap.parse_args()
    data = Path(a.data)
    raw = data / "raw" / "finngen"
    manifest = data / "raw" / "MANIFEST.tsv"

    for name in OBJECTS:
        dest = raw / name.rsplit("/", 1)[1]
        if a.force or not dest.exists():
            url = f"https://storage.googleapis.com/{BUCKET}/{name}"
            rec = common.download(url, dest)
            if rec["md5_b64"] != gcs_md5(name):
                raise SystemExit(f"md5 mismatch for {name}")
            rec["path"] = str(dest.relative_to(data))
            common.record_manifest(manifest, rec)
            print(f"downloaded {dest.name}: {rec['bytes']} bytes, md5 ok, sha256={rec['sha256'][:16]}...")

    raw_lof = finngen.read_lof(raw / LOF.rsplit("/", 1)[1])
    present = sorted(set(raw_lof["PHENO"]))
    frames, missing = [], []
    for trait, eps in finngen.ENDPOINTS.items():
        gone = [e for e in eps if e not in present]
        if gone:
            missing.append((trait, gone))
        df = finngen.build_trait(trait, raw_lof)
        print(f"{trait}: endpoints {eps}, {len(df)} genes")
        if len(df):
            frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(data / "burden_finngen.csv.gz", index=False)
    print(f"wrote {data / 'burden_finngen.csv.gz'}: {len(out)} rows")
    if missing:
        print("ENDPOINTS NOT IN THE LOF TABLE (recorded as absent, not substituted):", missing)


if __name__ == "__main__":
    main()
