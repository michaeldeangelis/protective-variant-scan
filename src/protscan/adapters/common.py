"""Shared helpers for the source adapters: streamed download with checksum log, SE reconstruction, meta-analysis."""
from __future__ import annotations

import base64
import csv
import datetime as dt
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

COLUMNS = ["gene", "trait", "mask", "cohort", "beta", "se", "p", "n_carriers", "n_total", "source"]
MANIFEST_COLUMNS = ["url", "path", "bytes", "sha256", "downloaded_utc"]
# float64 underflows near 1e-308; reported p-values below this are floored, not dropped.
P_FLOOR = 1e-300


def sha256_file(path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest, session=None, timeout: int = 300) -> dict:
    """Stream url to dest (atomic rename) and return a manifest record. Anonymous HTTPS only.

    md5_b64 is returned for comparison with the GCS objects' md5Hash; it is not written to the manifest.
    """
    import requests

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    h = hashlib.sha256()
    m = hashlib.md5()
    n = 0
    with (session or requests).get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(part, "wb") as fh:
            for block in r.iter_content(1 << 20):
                fh.write(block)
                h.update(block)
                m.update(block)
                n += len(block)
    part.replace(dest)
    return {
        "url": url,
        "path": str(dest),
        "bytes": n,
        "sha256": h.hexdigest(),
        "md5_b64": base64.b64encode(m.digest()).decode(),
        "downloaded_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def record_manifest(manifest_path, rec: dict) -> None:
    """Upsert one download record (keyed by path) into a TSV manifest."""
    manifest_path = Path(manifest_path)
    rows = {}
    if manifest_path.exists():
        with open(manifest_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows[row["path"]] = row
    rows[rec["path"]] = {k: rec[k] for k in MANIFEST_COLUMNS}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows.values())


def clip_p(p):
    return np.clip(np.asarray(p, dtype=float), P_FLOOR, 1.0)


def se_from_beta_p(beta, p):
    """Wald reconstruction se = |beta| / z, z = Phi^-1(1 - p/2). NaN where beta == 0 (se is unidentifiable)."""
    beta = np.asarray(beta, dtype=float)
    z = norm.isf(clip_p(p) / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.abs(beta) / z
    return np.where((beta == 0) | ~np.isfinite(se) | (se <= 0), np.nan, se)


def approx_carriers(af, n_total):
    """Expected carriers from a cumulative allele frequency: 2 * af * N, capped at N. Approximate."""
    n_total = np.asarray(n_total, dtype=float)
    return np.minimum(np.round(2.0 * np.asarray(af, dtype=float) * n_total), n_total)


def ivw_combine(frames, n_rule: str) -> pd.DataFrame:
    """Fixed-effect inverse-variance meta-analysis of several results for one trait, per (gene, mask).

    frames: DataFrames with gene, mask, beta, se, n_carriers, n_total (one per component analysis).
    n_rule: 'sum' when the components use disjoint people (sexes); 'max' when they reuse the same people.
    A gene present in only one component keeps that component's estimate.
    """
    if n_rule not in ("sum", "max"):
        raise ValueError("n_rule must be 'sum' or 'max'")
    df = pd.concat(frames, ignore_index=True)
    df = df[df["se"] > 0].copy()
    df["w"] = 1.0 / df["se"] ** 2
    df["wb"] = df["w"] * df["beta"]
    agg = {"w": "sum", "wb": "sum", "n_carriers": n_rule, "n_total": n_rule}
    g = df.groupby(["gene", "mask"], sort=True).agg(agg).reset_index()
    g["beta"] = g["wb"] / g["w"]
    g["se"] = 1.0 / np.sqrt(g["w"])
    g["p"] = clip_p(2.0 * norm.sf(np.abs(g["beta"] / g["se"])))
    return g[["gene", "mask", "beta", "se", "p", "n_carriers", "n_total"]]


def finalize(df: pd.DataFrame, trait: str, cohort: str, source: str) -> pd.DataFrame:
    out = df.copy()
    out["trait"] = trait
    out["cohort"] = cohort
    out["source"] = source
    out["n_carriers"] = out["n_carriers"].astype("int64")
    out["n_total"] = out["n_total"].astype("int64")
    return out[COLUMNS].reset_index(drop=True)
