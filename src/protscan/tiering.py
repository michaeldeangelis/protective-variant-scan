"""Tier A-D assignment and verdict-qualification flags."""
from __future__ import annotations

import pandas as pd

from . import stats
from .schema import Config


def assign_tier(adverse: list, rep_status: str, unscreened: bool = False) -> str:
    """A = discovery + replication + no adverse. B = discovery, not replicable, no adverse.
    C = discovery with an adverse trade-off. D = discovery, replication attempted and failed.
    C1: a would-be A with too few trade-off outcomes screened is capped at B ("trade-off unscreened")."""
    if adverse:
        return "C"
    if rep_status == "replicated":
        return "B" if unscreened else "A"
    if rep_status == "failed":
        return "D"
    return "B"


def build_hits(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """One row per discovery hit (gene x panel trait) with all annotations and its tier."""
    hits = stats.discovery_hits(df, cfg)
    if hits.empty:
        cols = ["gene", "trait", "beta", "se", "p", "n_carriers", "domain", "rep_status", "rep_trait", "rep_beta",
                "rep_p_onesided", "adverse", "n_tradeoff_tested", "tradeoff_unscreened", "mask_status", "eur_status", "tier",
                "lipid_gene", "qualifies"]
        return pd.DataFrame(columns=cols)
    hits = hits.join(stats.replication_status(hits, df, cfg))
    genes = sorted(set(hits["gene"]))
    adv = stats.adverse_tradeoffs(genes, df, cfg)
    tested = stats.tradeoff_tested(genes, df, cfg)
    hits["domain"] = hits["trait"].map(cfg.panel_domain)
    hits["adverse"] = hits["gene"].map(adv)
    hits["n_tradeoff_tested"] = hits["gene"].map(tested)
    hits["tradeoff_unscreened"] = hits["n_tradeoff_tested"] < cfg.min_tradeoffs_screened
    hits["mask_status"] = stats.mask_status(hits, df, cfg)
    hits["eur_status"] = stats.eur_status(hits, df, cfg)
    hits["tier"] = [assign_tier(a, s, u) for a, s, u in zip(hits["adverse"], hits["rep_status"], hits["tradeoff_unscreened"])]
    hits["lipid_gene"] = hits["gene"].isin(cfg.lipid_genes)
    # verdict qualification: non-lipid gene, qualifying trait domain, both masks consistent
    hits["qualifies"] = (~hits["lipid_gene"]) & hits["domain"].isin(cfg.qualifying_domains) & (hits["mask_status"] == "consistent")
    return hits.sort_values(["tier", "p"]).reset_index(drop=True)
