"""FinnGen R13 LoF gene-burden results (regenie, Finnish biobanks) -> normalized replication rows.

Independent of UK Biobank. pLoF only (frameshift, splice donor/acceptor, stop gained; no LOFTEE, imputed genotypes),
binary endpoints only. Endpoint choices were fixed 2026-09-29 before any result was read (docs/data-sources.md 3.1-3.2).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import common

SOURCE = "finngen_r13"
COHORT = "replication"
MASK = "plof"

# Normalized replication trait -> FinnGen endpoint code(s). The first three are the declared proxies (Amendment 1);
# the rest are the supplementary trade-off outcomes.
ENDPOINTS = {
    "hypercholesterolemia": ["E4_HYPERCHOL"],
    "obesity": ["E4_OBESITY"],
    "hypertension": ["I9_HYPTENS"],
    "type_2_diabetes": ["T2D"],
    "coronary_disease": ["I9_CHD"],
    "cancer_any": ["C3_CANCER"],
    "dementia": ["F5_DEMENTIA"],
    "major_depression": ["F5_DEPRESSIO"],
    "schizophrenia": ["F5_SCHZPHR"],
    "infertility": ["N14_FEMALEINFERT", "N14_MALEINFERT"],
}
# Female and male infertility endpoints are disjoint samples.
N_RULE = {"infertility": "sum"}

ABSENT = {
    "all_cause_mortality": "no all-cause death endpoint in the R13 manifest",
    "fracture": "only site-specific ST19_FRACT_* endpoints; no any-fracture endpoint; not substituted",
    "ldl": "continuous trait; FinnGen LoF burden covers binary endpoints only (proxy: hypercholesterolemia)",
    "bmi": "continuous trait; proxy: obesity",
    "systolic_bp": "continuous trait; proxy: hypertension",
}

# regenie names each burden set "<gene>.<mask>.<maf>", e.g. "A2ML1.Mask1.0.01"; only this one mask is present.
MASK_SUFFIX = r"(?i)\.mask\d+\.[0-9.]+$"

# The LoF carrier allele is A1 (ALLELE1), so a collapsed-genotype frequency above 0.5 means the effect sign may refer to
# the other allele. Such rows are dropped, not flipped; NaN also fails the check.
A1FREQ_MAX = 0.5
log = logging.getLogger(__name__)

USECOLS = ["PHENO", "ID", "A1FREQ", "N", "TEST", "BETA", "SE", "LOG10P"]


def read_lof(path, endpoints=None, chunksize: int = 500_000) -> pd.DataFrame:
    """Stream the (gzipped) merged LoF table and keep only the requested endpoints and the additive test."""
    if endpoints is None:
        endpoints = {e for v in ENDPOINTS.values() for e in v}
    endpoints = set(endpoints)
    keep = []
    reader = pd.read_csv(
        path, sep="\t", usecols=USECOLS, dtype={"PHENO": str, "ID": str, "TEST": str}, chunksize=chunksize
    )
    for chunk in reader:
        sub = chunk[chunk["PHENO"].isin(endpoints) & (chunk["TEST"] == "ADD")]
        if len(sub):
            keep.append(sub)
    if not keep:
        return pd.DataFrame(columns=USECOLS)
    return pd.concat(keep, ignore_index=True)


def guard_a1freq(df: pd.DataFrame):
    """Keep rows with A1FREQ <= 0.5 (NaN fails). Returns (kept, dropped-row counts per PHENO) and logs any drops."""
    ok = df["A1FREQ"] <= A1FREQ_MAX
    dropped = df.loc[~ok].groupby("PHENO").size()
    if len(dropped):
        log.warning("dropped %d rows with A1FREQ > %s or missing: %s", int(dropped.sum()), A1FREQ_MAX, dropped.to_dict())
    return df.loc[ok], dropped


def conversion_summary(raw: pd.DataFrame) -> pd.DataFrame:
    """Per-endpoint row counts: read, dropped by the A1FREQ guard, kept with undefined se (se <= 0 or missing)."""
    kept, dropped = guard_a1freq(raw)
    out = pd.DataFrame({"rows_read": raw.groupby("PHENO").size()})
    out["dropped_a1freq"] = dropped.reindex(out.index).fillna(0).astype(int)
    out["kept_se_undefined"] = (~(kept["SE"] > 0)).groupby(kept["PHENO"]).sum().reindex(out.index).fillna(0).astype(int)
    return out.reset_index().rename(columns={"PHENO": "endpoint"})


def normalize_endpoint(raw: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    """One endpoint -> gene, mask, beta, se, p, n_carriers, n_total.

    The table has no carrier count: n_carriers ~ 2 * A1FREQ * N, where A1FREQ is the frequency of the collapsed
    (max-over-sites) genotype; approximate. Rows failing the A1FREQ <= 0.5 guard are dropped (see guard_a1freq). Rows
    with an undefined se (SE <= 0 or missing) are kept with se = NaN and p intact.
    """
    df = guard_a1freq(raw[raw["PHENO"] == endpoint])[0].rename(columns={"ID": "gene", "BETA": "beta", "SE": "se", "N": "n_total"})
    df["gene"] = df["gene"].str.replace(MASK_SUFFIX, "", regex=True)
    df = df.dropna(subset=["gene", "beta", "LOG10P"]).copy()
    df["se"] = df["se"].where(df["se"] > 0)
    df["p"] = common.clip_p(np.power(10.0, -df["LOG10P"].astype(float)))
    df["n_total"] = df["n_total"].astype("int64")
    df["n_carriers"] = common.approx_carriers(df["A1FREQ"].fillna(0.0), df["n_total"])
    df["mask"] = MASK
    return df[["gene", "mask", "beta", "se", "p", "n_carriers", "n_total"]].reset_index(drop=True)


def build_trait(trait: str, raw: pd.DataFrame) -> pd.DataFrame:
    """Normalized rows for one replication trait; IVW-combines traits defined by two endpoints."""
    per = [normalize_endpoint(raw, e) for e in ENDPOINTS[trait]]
    per = [d for d in per if len(d)]
    if not per:
        return pd.DataFrame(columns=common.COLUMNS)
    df = per[0] if len(per) == 1 else common.ivw_combine(per, N_RULE[trait])
    return common.finalize(df, trait, COHORT, SOURCE)
