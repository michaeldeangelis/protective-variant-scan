"""gnomAD v4.1 constraint metrics -> constraint table (gene, loeuf, pli) used only to annotate report output."""
from __future__ import annotations

import pandas as pd

URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
COLUMNS = ["gene", "loeuf", "pli"]


def _flag(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().eq("true")


def normalize_constraint(df: pd.DataFrame) -> pd.DataFrame:
    """One row per gene symbol.

    Transcript choice: MANE Select, else canonical. The file lists each transcript under both an Ensembl and a RefSeq
    identifier; the Ensembl-keyed row (gene_id starting ENSG) is preferred. Genes with no LOEUF are dropped.
    """
    d = df.copy()
    d["_mane"] = _flag(d["mane_select"])
    d["_canon"] = _flag(d["canonical"])
    d["_ens"] = d["gene_id"].astype(str).str.startswith("ENSG")
    d = d[d["_mane"] | d["_canon"]]
    d = d.dropna(subset=["gene", "lof.oe_ci.upper"])
    d = d.sort_values(["gene", "_mane", "_ens", "transcript"], ascending=[True, False, False, True])
    d = d.drop_duplicates("gene", keep="first")
    out = d.rename(columns={"lof.oe_ci.upper": "loeuf", "lof.pLI": "pli"})[COLUMNS]
    return out.reset_index(drop=True)
