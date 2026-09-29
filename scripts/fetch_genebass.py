#!/usr/bin/env python
"""Fetch Genebass gene-based burden results via the open API behind app.genebass.org and write <data>/burden_genebass.csv.gz.

No login, no billing project (the Hail tables in gs://ukbb-exome-public are requester-pays and are not used).
One request at a time with a delay: ~19 analyses x 3 burden sets, cached under <data>/raw/genebass/.
"""
import argparse
import json
import sys
import time
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402
import requests  # noqa: E402

from protscan.adapters import common, genebass  # noqa: E402

HEADERS = {"User-Agent": "protscan-fetch/0.1 (research; polite, single-threaded)"}
DELAY_S = 1.0


def fetch(session, url: str, dest: Path, data: Path, manifest: Path, force: bool) -> None:
    if dest.exists() and not force:
        return
    rec = common.download(url, dest, session=session)
    rec["path"] = str(dest.relative_to(data))
    common.record_manifest(manifest, rec)
    print(f"downloaded {dest.name}: {rec['bytes']} bytes sha256={rec['sha256'][:16]}...", flush=True)
    time.sleep(DELAY_S)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(ROOT / "data"))
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    a = ap.parse_args()
    data = Path(a.data)
    raw = data / "raw" / "genebass"
    manifest = data / "raw" / "MANIFEST.tsv"
    s = requests.Session()
    s.headers.update(HEADERS)

    fetch(s, f"{genebass.API}/phenotypes", raw / "phenotypes.json", data, manifest, a.force)
    meta = {m["analysis_id"]: m for m in json.loads((raw / "phenotypes.json").read_text())}
    wanted = [aid for ids in genebass.TRAITS.values() for aid in ids]
    missing = [aid for aid in wanted if aid not in meta]
    if missing:
        raise SystemExit(f"analysis ids not in Genebass /phenotypes: {missing}")

    qc = {}
    for bset in genebass.BURDEN_SET_TO_MASK:
        dest = raw / f"gene_qc_{genebass.BURDEN_SET_TO_MASK[bset]}.json"
        fetch(s, f"{genebass.API}/gene-qc/burden-set/{quote(bset, safe='')}", dest, data, manifest, a.force)
        qc[bset] = json.loads(dest.read_text())

    frames = []
    for trait, ids in genebass.TRAITS.items():
        per = []
        for aid in ids:
            parts = []
            for bset, mask in genebass.BURDEN_SET_TO_MASK.items():
                dest = raw / f"gene_manhattan_{aid}_{mask}.json"
                url = f"{genebass.API}/analysis/{quote(aid, safe='')}/gene-manhattan?burdenSet={quote(bset, safe='')}"
                fetch(s, url, dest, data, manifest, a.force)
                recs = json.loads(dest.read_text())
                if not isinstance(recs, list) or not recs or "BETA_Burden" not in recs[0]:
                    raise SystemExit(f"unexpected response for {aid} / {bset}")
                parts.append(genebass.normalize_analysis(recs, qc[bset], mask, genebass.analysis_n_total(meta[aid])))
            per.append(pd.concat(parts, ignore_index=True))
        df = genebass.build_trait(trait, per)
        print(f"{trait}: {ids} -> {len(df)} rows", flush=True)
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(data / "burden_genebass.csv.gz", index=False)
    print(f"wrote {data / 'burden_genebass.csv.gz'}: {len(out)} rows")


if __name__ == "__main__":
    main()
