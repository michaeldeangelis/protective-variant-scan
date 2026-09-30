"""Per-rule statistics on the normalized table. Every function is a direct reading of one ledger rule."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .schema import Config

CHI2_MEDIAN = float(chi2.ppf(0.5, 1))


def lambda_gc(p) -> float:
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1.0)
    return float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)


def one_sided_p(beta, p, sign):
    """One-sided p for the direction `sign` (+1: beta>0, -1: beta<0) from the two-sided p."""
    p = np.asarray(p, dtype=float)
    return np.where(np.asarray(beta) * sign > 0, p / 2, 1 - p / 2)


def _rows(df, cohort, mask):
    return df[(df["cohort"] == cohort) & (df["mask"] == mask)]


def discovery_hits(df: pd.DataFrame, cfg: Config, mask: str = "plof") -> pd.DataFrame:
    """Gene-level p < discovery threshold in the beneficial direction, panel traits only."""
    r = _rows(df, "discovery", mask)
    r = r[r["trait"].isin(cfg.panel)]
    sign = r["trait"].map(cfg.panel_sign)
    r = r[(r["p"] < cfg.discovery_p) & (r["beta"] * sign > 0)]
    return r[["gene", "trait", "beta", "se", "p", "n_carriers"]].reset_index(drop=True)


def _source_independent(src: pd.Series, cfg: Config) -> pd.Series:
    s = src.astype(str).str.lower()
    allow = pd.Series(False, index=s.index)
    deny = pd.Series(False, index=s.index)
    for tok in cfg.independent_sources:
        allow |= s.str.contains(tok, regex=False)
    for tok in cfg.ukb_sources:
        deny |= s.str.contains(tok, regex=False)
    return allow & ~deny


def independent_replication(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """plof replication rows from an allow-listed source that does not overlap UK Biobank (Amendment 1)."""
    r = _rows(df, "replication", "plof")
    return r[_source_independent(r["source"], cfg)]


def informative(r: pd.DataFrame) -> pd.DataFrame:
    """C11a: rows usable as evidence: finite se > 0 and beta != 0. Other rows serve lambda_GC and coverage only."""
    return r[np.isfinite(r["se"]) & (r["se"] > 0) & (r["beta"] != 0)]


def replication_status(hits: pd.DataFrame, df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Independent-cohort replication. Only traits with a declared proxy can replicate (Amendment 1).

    rep_status: replicated | failed | trait_missing | gene_untested
    """
    rep = informative(independent_replication(df, cfg)).set_index(["gene", "trait"])
    rep_traits = set(rep.index.get_level_values("trait"))
    out = []
    for g, t in zip(hits["gene"], hits["trait"]):
        rec = {"rep_status": "trait_missing", "rep_trait": None, "rep_beta": np.nan, "rep_p_onesided": np.nan}
        if t in cfg.proxies:
            cands = [t, cfg.proxies[t]["trait"]]
            avail = [c for c in cands if c in rep_traits]
            if avail:
                rec["rep_status"] = "gene_untested"
                for c in avail:
                    if (g, c) in rep.index:
                        row = rep.loc[(g, c)]
                        p1 = float(one_sided_p(row["beta"], row["p"], cfg.sign(c)))
                        rec.update(rep_trait=c, rep_beta=float(row["beta"]), rep_p_onesided=p1,
                                   rep_status="replicated" if p1 < cfg.replication_p else "failed")
                        break
        out.append(rec)
    return pd.DataFrame(out, index=hits.index)


def adverse_tradeoffs(genes, df: pd.DataFrame, cfg: Config) -> dict:
    """gene -> list of adverse outcomes: one-sided p in the harmful direction < alpha/N_tradeoff.

    Screens the plof mask in every cohort that reports the outcome.
    """
    r = df[(df["mask"] == "plof") & (df["cohort"].isin(["discovery", "replication"]))
           & (df["trait"].isin(cfg.tradeoff)) & (df["gene"].isin(set(genes)))].copy()
    out = {g: [] for g in genes}
    if r.empty:
        return out
    harmful = -r["trait"].map(cfg.tradeoff_sign)
    r["p_harm"] = one_sided_p(r["beta"], r["p"], harmful.to_numpy())
    for row in r[r["p_harm"] < cfg.tradeoff_p].itertuples():
        out[row.gene].append({"trait": row.trait, "cohort": row.cohort, "beta": float(row.beta),
                              "p_onesided": float(row.p_harm)})
    return out


def tradeoff_tested(genes, df: pd.DataFrame, cfg: Config) -> dict:
    """Distinct trade-off outcomes screened per gene: plof discovery rows plus allow-listed independent replication
    rows only (C10g); untrusted replication rows do not count toward the C1/C9 minimum."""
    genes = set(genes)
    r = pd.concat([_rows(df, "discovery", "plof"), independent_replication(df, cfg)])
    r = r[r["trait"].isin(cfg.tradeoff) & r["gene"].isin(genes)]
    n = r.groupby("gene")["trait"].nunique()
    return {g: int(n.get(g, 0)) for g in genes}


def harmful_panel_hits(genes, df: pd.DataFrame, cfg: Config) -> dict:
    """Informational (C10e): gene -> panel traits with a discovery-threshold association in the HARMFUL direction."""
    r = _rows(df, "discovery", "plof")
    r = r[r["trait"].isin(cfg.panel) & r["gene"].isin(set(genes)) & (r["p"] < cfg.discovery_p)]
    r = r[r["beta"] * r["trait"].map(cfg.panel_sign) < 0]
    out = {g: [] for g in genes}
    for row in r.itertuples():
        out[row.gene].append({"trait": row.trait, "beta": float(row.beta), "p": float(row.p)})
    return out


def mask_status(hits: pd.DataFrame, df: pd.DataFrame, cfg: Config) -> list:
    """Direction of the dmis mask in discovery vs the beneficial direction: consistent | inconsistent | missing."""
    d = _rows(df, "discovery", "dmis").set_index(["gene", "trait"])["beta"]
    out = []
    for g, t in zip(hits["gene"], hits["trait"]):
        if (g, t) not in d.index:
            out.append("missing")
        else:
            out.append("consistent" if d.loc[(g, t)] * cfg.panel_sign[t] > 0 else "inconsistent")
    return out


def eur_status(hits: pd.DataFrame, df: pd.DataFrame, cfg: Config) -> list:
    """EUR-only discovery re-analysis: held (still passes discovery rule) | failed | not_reported."""
    e = _rows(df, "discovery_eur", "plof").set_index(["gene", "trait"])
    out = []
    for g, t in zip(hits["gene"], hits["trait"]):
        if (g, t) not in e.index:
            out.append("not_reported")
        else:
            row = e.loc[(g, t)]
            out.append("held" if row["p"] < cfg.discovery_p and row["beta"] * cfg.panel_sign[t] > 0 else "failed")
    return out
