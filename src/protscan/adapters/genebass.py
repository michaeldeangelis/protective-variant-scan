"""Genebass (UK Biobank exomes) gene-based burden results -> normalized discovery rows.

Data come from the open REST API behind app.genebass.org (no login, no billing project); the Hail tables are
requester-pays and are not used. See docs/data-sources.md for the route, licence and the trait mapping rationale.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import common

SOURCE = "genebass"
COHORT = "discovery"
API = "https://main.genebass.org/api"

# Genebass burdenSet -> normalized mask. Genebass has no REVEL/CADD-filtered set: "missense|LC" (missense plus
# low-confidence pLoF) is the closest available consistency mask and is labelled dmis.
BURDEN_SET_TO_MASK = {"pLoF": "plof", "missense|LC": "dmis", "synonymous": "syn"}

# Normalized trait -> Genebass analysis_id(s). Fixed 2026-09-29 before any result was read (docs/data-sources.md 2.1).
TRAITS = {
    "fluid_intelligence": ["continuous-20016-both_sexes--irnt"],
    "reaction_time": ["continuous-20023-both_sexes--irnt"],
    "hand_grip_strength": ["continuous-hand_grip_strength_custom-both_sexes--custom"],
    "fev1": ["continuous-3063-both_sexes--irnt"],
    "resting_heart_rate": ["continuous-102-both_sexes--irnt"],
    "systolic_bp": ["continuous-4080-both_sexes--irnt"],
    "ldl": ["continuous-30780-both_sexes--irnt"],
    "bmi": ["continuous-21001-both_sexes--irnt"],
    "parental_lifespan": ["continuous-1807-both_sexes--irnt", "continuous-3526-both_sexes--irnt"],
    "triglycerides": ["continuous-30870-both_sexes--irnt"],
    "type_2_diabetes": ["categorical-T2D_custom-both_sexes--custom"],
    "coronary_disease": ["categorical-CAD_custom-both_sexes--custom"],
    "cancer_any": ["categorical-2453-both_sexes--"],
    "dementia": ["icd_first_occurrence-130842-both_sexes--"],
    "major_depression": ["icd_first_occurrence-130894-both_sexes--"],
    "schizophrenia": ["icd_first_occurrence-130874-both_sexes--"],
    "fracture": ["categorical-2463-both_sexes--"],
    "infertility": ["icd_first_occurrence-132156-both_sexes--", "icd_first_occurrence-132084-both_sexes--"],
}
# How multi-analysis traits pool sample sizes: parents' lifespans reuse the same participants, sexes are disjoint.
N_RULE = {"parental_lifespan": "max", "infertility": "sum"}

ABSENT = {
    "numeric_memory": "UKB field 4282 is not among the 4,529 Genebass analyses",
    "pairs_matching_errors": "UKB field 399 is not among the 4,529 Genebass analyses",
    "education_years": "UKB fields 845/6138/22501 are not among the 4,529 Genebass analyses",
    "walking_pace": "UKB field 924 is not among the 4,529 Genebass analyses",
    "all_cause_mortality": "only 'Age at death' (40007, death register, n=11,391) exists; not all-cause mortality, not substituted",
}


def analysis_n_total(meta: dict) -> int:
    """Sample size of one analysis from /phenotypes metadata: cases + controls (controls are null for continuous)."""
    return int(meta["n_cases"]) + int(meta.get("n_controls") or 0)


def qc_pass(qc_records) -> pd.DataFrame:
    """Gene QC for one burden set: keep genes passing Genebass coverage and variant-count flags.

    The lambda-based flags (keep_gene_burden, keep_gene_skat, keep_gene_skato) are deliberately not applied: they
    are derived from synonymous-mask inflation, so filtering on them would make the synonymous negative control circular.
    """
    q = pd.DataFrame(qc_records)
    keep = (q["keep_gene_coverage"] == True) & (q["keep_gene_n_var"] == True)  # noqa: E712 (null counts as fail)
    return q.loc[keep, ["gene_id", "CAF"]].drop_duplicates("gene_id")


def normalize_analysis(records, qc_records, mask: str, n_total: int) -> pd.DataFrame:
    """One analysis x one burden set -> gene, mask, beta, se, p, n_carriers, n_total.

    p is Pvalue_Burden so that p, beta and the reconstructed se describe the same test (SKAT-O p is not used).
    The gene-manhattan endpoint returns no SE and no carrier count: se = |beta| / z(p); n_carriers ~ 2 * CAF * N
    (CAF is gene-level and not phenotype-specific, so the count is approximate). Symbols shared by more than one
    Ensembl gene are ambiguous and dropped.
    """
    df = pd.DataFrame(records)
    df = df.rename(columns={"gene_symbol": "gene", "BETA_Burden": "beta", "Pvalue_Burden": "p"})
    df = df[["gene", "gene_id", "beta", "p"]].dropna(subset=["beta", "p", "gene"])
    df = df[df["gene"].astype(str).str.len() > 0]
    df = df.merge(qc_pass(qc_records), on="gene_id", how="inner")
    df = df[~df["gene"].duplicated(keep=False)].copy()
    df["se"] = common.se_from_beta_p(df["beta"], df["p"])
    df["p"] = common.clip_p(df["p"])
    df = df.dropna(subset=["se"])
    df["n_total"] = int(n_total)
    df["n_carriers"] = common.approx_carriers(df["CAF"].fillna(0.0), df["n_total"])
    df["mask"] = mask
    return df[["gene", "mask", "beta", "se", "p", "n_carriers", "n_total"]].reset_index(drop=True)


def build_trait(trait: str, per_analysis: list) -> pd.DataFrame:
    """Combine the per-analysis frames of one trait (IVW when the trait has two analyses) and finalize."""
    if len(per_analysis) == 1:
        df = per_analysis[0]
    else:
        df = common.ivw_combine(per_analysis, N_RULE[trait])
    return common.finalize(df, trait, COHORT, SOURCE)
