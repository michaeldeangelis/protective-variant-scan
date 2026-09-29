"""Coverage table: which declared trait is present in which source, and why not where absent."""
from __future__ import annotations

import pandas as pd

from . import finngen, genebass

COLUMNS = ["trait", "kind", "discovery_source", "discovery_genes_plof", "discovery_genes_dmis", "discovery_genes_syn",
           "replication_trait", "replication_source", "replication_genes_plof", "replication_role"]


def _genes(df: pd.DataFrame, trait: str, mask: str) -> int:
    if df is None or df.empty:
        return 0
    return int(((df["trait"] == trait) & (df["mask"] == mask)).sum())


def build_coverage(cfg, discovery: pd.DataFrame, replication: pd.DataFrame) -> pd.DataFrame:
    """One row per declared trait (panel, trade-off, control). Absences carry the reason, never a substitute."""
    proxy_of = {k: v["trait"] for k, v in cfg.proxies.items()}
    kinds = [(t, "panel") for t in cfg.panel] + [(t, "tradeoff") for t in cfg.tradeoff] + [(t, "control") for t in cfg.control_sign]
    rows = []
    for trait, kind in kinds:
        row = dict.fromkeys(COLUMNS)
        row.update(trait=trait, kind=kind)
        if trait in genebass.TRAITS:
            row["discovery_source"] = "; ".join(genebass.TRAITS[trait])
        else:
            row["discovery_source"] = "ABSENT: " + genebass.ABSENT[trait]
        for m in ("plof", "dmis", "syn"):
            row[f"discovery_genes_{m}"] = _genes(discovery, trait, m)
        rep_trait = proxy_of.get(trait, trait)
        row["replication_trait"] = rep_trait
        if rep_trait in finngen.ENDPOINTS:
            row["replication_source"] = "; ".join(finngen.ENDPOINTS[rep_trait])
            row["replication_role"] = "declared proxy (Tier A possible)" if trait in proxy_of else (
                "supplementary trade-off rows" if kind == "tradeoff" else "none")
        else:
            row["replication_source"] = "ABSENT: " + finngen.ABSENT.get(trait, "not analysed in FinnGen LoF burden (binary endpoints only)")
            row["replication_role"] = "none (Tier B ceiling)" if kind == "panel" else "none"
        row["replication_genes_plof"] = _genes(replication, rep_trait, "plof")
        rows.append(row)
    return pd.DataFrame(rows, columns=COLUMNS)
