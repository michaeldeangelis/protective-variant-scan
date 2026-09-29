"""Independent ledger-conformance tests.

Source of truth: experiments.md, entry 2026-09-29 protective-variant-scan plus
AMENDMENT 1 (same date). Constants are typed from that text, not read from the
config, except in tests that compare the config TO these constants.
Owned by the reviewer; builders must not edit this file.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "prereg.yaml"

# ---- ledger constants (literal, typed from experiments.md) ----
L_ALPHA = 0.05
L_N_GENES = 20_000
L_N_TARGET = 13
L_N_TRADEOFF = 9
L_DISCOVERY_P = 0.05 / (20_000 * 13)          # ledger prints 1.9e-7 (exactly 1.923e-7)
L_DISCOVERY_P_PRINTED = 1.9e-7
L_ADVERSE_P = 0.05 / 9                        # 5.556e-3
L_REPL_ONE_SIDED_P = 0.05
L_LAMBDA_GC_MAX = 1.10                        # strict: lambda_GC < 1.10
L_SYN_HITS_MAX = 0                            # Amendment 1: zero synonymous genes at threshold

# direction of benefit: +1 higher is better, -1 lower is better
L_PANEL = dict(
    fluid_intelligence=("cognitive", +1),
    reaction_time=("cognitive", -1),              # faster
    numeric_memory=("cognitive", +1),
    pairs_matching_errors=("cognitive", -1),      # fewer
    education_years=("cognitive", +1),            # confounded proxy
    hand_grip_strength=("physical", +1),
    fev1=("physical", +1),
    walking_pace=("physical", +1),                # faster
    resting_heart_rate=("physical", -1),          # lower
    systolic_bp=("physical", -1),                 # lower
    ldl=("metabolic", -1),                        # lower
    bmi=("metabolic", -1),                        # lower
    parental_lifespan=("metabolic", +1),          # longer
)
# adverse direction = higher risk or hazard (beta > 0)
L_TRADEOFF = [
    "type_2_diabetes", "coronary_disease", "cancer_any", "dementia", "major_depression",
    "schizophrenia", "all_cause_mortality", "fracture", "infertility",
]
# Amendment 1: the only declared replication proxies (panel trait -> proxy endpoint).
L_PROXIES = dict(ldl="hypercholesterolemia", bmi="obesity", systolic_bp="hypertension")
L_NO_PROXY = sorted(set(L_PANEL) - set(L_PROXIES))   # Tier B is the ceiling for these
L_UKB_SOURCES = ("genebass", "azphewas", "astrazeneca", "regeneron", "rgc", "ukb", "uk_biobank")
L_VERDICTS = set(["PASS", "LEAD", "KILL"])           # Amendment 1; PROVISIONAL retired


# Tier letters (ledger):
#  A = discovery + independent replication + no adverse trade-off
#  B = discovery, not replicable (trait missing in replication source), no adverse
#  C = discovery with an adverse trade-off
#  D = discovery, replication attempted and failed
def oracle_beneficial(beta, benefit):
    return beta * benefit > 0


def oracle_discovery(beta, p, benefit):
    return p < L_DISCOVERY_P and oracle_beneficial(beta, benefit)


def oracle_replicates(beta, p_two_sided, benefit):
    """One-sided p in the beneficial direction < 0.05 (same direction required)."""
    if not oracle_beneficial(beta, benefit):
        return False
    return (p_two_sided / 2.0) < L_REPL_ONE_SIDED_P


def oracle_adverse(beta, p, adverse_sign=+1):
    return p < L_ADVERSE_P and beta * adverse_sign > 0


# ======================================================================================
# (a) config equals ledger constants
# ======================================================================================
@pytest.fixture(scope="module")
def cfg_yaml():
    assert CONFIG_PATH.exists(), "config/prereg.yaml missing"
    return yaml.safe_load(CONFIG_PATH.read_text())


def test_ledger_constants_self_consistent():
    assert math.isclose(L_DISCOVERY_P, 1.923076923e-7, rel_tol=1e-6)
    assert math.isclose(L_ADVERSE_P, 0.005555555, rel_tol=1e-6)
    assert len(L_PANEL) == L_N_TARGET == 13
    assert len(L_TRADEOFF) == L_N_TRADEOFF == 9
    assert len(L_NO_PROXY) == 10


def test_config_panel_traits_and_direction_of_benefit(cfg_yaml):
    panel = cfg_yaml["panel"]
    assert set(panel) == set(L_PANEL), set(panel) ^ set(L_PANEL)
    assert len(panel) == L_N_TARGET
    for trait, (domain, sign) in L_PANEL.items():
        assert panel[trait]["benefit"] == sign, "%s: benefit sign" % trait
        assert panel[trait]["domain"] == domain, "%s: domain" % trait


def test_config_tradeoff_panel_and_adverse_direction(cfg_yaml):
    t = cfg_yaml["tradeoff"]
    assert set(t) == set(L_TRADEOFF)
    assert len(t) == L_N_TRADEOFF
    # adverse = higher risk, so the benefit sign must be -1 for every trade-off outcome
    assert all(v["benefit"] == -1 for v in t.values())


def test_config_discovery_threshold(cfg_yaml):
    d = cfg_yaml["discovery"]
    assert d["alpha"] == L_ALPHA
    assert d["n_genes"] == L_N_GENES
    assert d["n_traits"] == L_N_TARGET
    thr = d["p_threshold"]
    # both numbers printed in the ledger are acceptable: 1.9e-7 and 0.05/(20000*13)
    assert L_DISCOVERY_P_PRINTED <= thr <= L_DISCOVERY_P * (1 + 1e-9), thr


def test_config_replication_rule(cfg_yaml):
    r = cfg_yaml["replication"]
    assert r["one_sided_p"] == L_REPL_ONE_SIDED_P
    assert dict((k, v["trait"]) for k, v in r["proxies"].items()) == L_PROXIES
    for k, v in r["proxies"].items():
        assert v["benefit"] == L_PANEL[k][1], "proxy %s direction must equal parent trait" % k


def test_config_adverse_threshold(cfg_yaml):
    assert cfg_yaml["tradeoff_screen"]["alpha"] == L_ALPHA  # divided by N_tradeoff = 9 in code


def test_config_controls(cfg_yaml):
    c = cfg_yaml["controls"]
    assert c["lambda_gc_max"] == L_LAMBDA_GC_MAX
    assert c["syn_hits_max"] == L_SYN_HITS_MAX
    pos = c["positive"]
    genesets = [tuple(sorted(p["genes"])) for p in pos]
    assert ("PCSK9",) in genesets
    assert ("ANGPTL4", "APOC3") in genesets
    ldl = [p for p in pos if p["genes"] == ["PCSK9"] and p["trait"] == "ldl"]
    cad = [p for p in pos if p["genes"] == ["PCSK9"] and p["trait"] == "coronary_disease"]
    assert ldl and ldl[0]["at_discovery_threshold"] is True    # LDL lower at the discovery threshold
    assert cad and cad[0]["at_discovery_threshold"] is False   # coronary-disease direction protective
    lip = [p for p in pos if sorted(p["genes"]) == ["ANGPTL4", "APOC3"]]
    assert lip and lip[0]["at_discovery_threshold"] is True


def test_config_masks_and_permutation_retired(cfg_yaml):
    m = cfg_yaml["masks"]
    assert m["primary"] == "plof" and m["negative_control"] == "syn"
    assert "permut" not in yaml.safe_dump(cfg_yaml).lower(), "Amendment 1 retired the permuted run"


def test_config_verdict_domains_and_lipid_genes(cfg_yaml):
    v = cfg_yaml["verdict"]
    assert set(v["qualifying_domains"]) == set(["cognitive", "physical"])
    lipid = set(g.upper() for g in v["lipid_pathway_genes"])
    assert set(["PCSK9", "ANGPTL4", "APOC3"]) <= lipid


def test_config_loads_via_schema_and_derives_thresholds():
    from protscan.schema import load_config
    c = load_config(CONFIG_PATH)
    assert len(c.panel) == 13 and len(c.tradeoff) == 9
    assert math.isclose(c.tradeoff_p, L_ADVERSE_P, rel_tol=1e-12)
    assert L_DISCOVERY_P_PRINTED <= c.discovery_p <= L_DISCOVERY_P * (1 + 1e-9)
    assert c.replication_p == L_REPL_ONE_SIDED_P
    assert c.lambda_gc_max == L_LAMBDA_GC_MAX
    assert c.syn_hits_max == 0
    for t, (_, s) in L_PANEL.items():
        assert c.sign(t) == s
    for t in L_TRADEOFF:
        assert c.sign(t) == -1


# ======================================================================================
# (b) rule behaviour on hand-built tables (unit level, independent oracle)
# ======================================================================================
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from protscan import controls, stats, tiering
from protscan.schema import COLUMNS, load_config, validate

DISC_SRC = "genebass"
REP_SRC = "finngen_r12"
CTRL_TRAITS = ["triglycerides"]


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG_PATH)


def benefit_sign(trait):
    """Direction of benefit from the LEDGER constants (not from the config)."""
    if trait in L_PANEL:
        return L_PANEL[trait][1]
    if trait in L_TRADEOFF or trait in CTRL_TRAITS:
        return -1
    for k, v in L_PROXIES.items():
        if v == trait:
            return L_PANEL[k][1]
    raise KeyError(trait)


def good(trait, mag=0.5):
    return benefit_sign(trait) * mag


def bad(trait, mag=0.5):
    return -benefit_sign(trait) * mag


def R(gene, trait, mask, cohort, beta, p, source=None):
    if source is None:
        source = REP_SRC if cohort == "replication" else DISC_SRC
    return [gene, trait, mask, cohort, beta, 0.1, p, 100, 400000, source]


def T(rows):
    return validate(pd.DataFrame(rows, columns=COLUMNS))


def hit_rows(gene, trait, p=1e-9, dmis_good=True):
    return [
        R(gene, trait, "plof", "discovery", good(trait), p),
        R(gene, trait, "dmis", "discovery", good(trait) if dmis_good else bad(trait), 0.3),
    ]


def null_tradeoffs(gene, cohorts=("discovery", "replication")):
    rows = []
    for c in cohorts:
        for t in L_TRADEOFF:
            rows.append(R(gene, t, "plof", c, 0.01, 0.5))
    return rows


def tier_map(rows, cfg):
    hits = tiering.build_hits(T(rows), cfg)
    return dict(((r.gene, r.trait), r.tier) for r in hits.itertuples())


def drop_rows(rows, gene, trait, cohort):
    return [r for r in rows if not (r[0] == gene and r[1] == trait and r[3] == cohort)]


def sbp_replicated_gene(gene="GA", rep_p=0.02):
    """SBP discovery hit + hypertension replication in the beneficial direction."""
    rows = hit_rows(gene, "systolic_bp") + null_tradeoffs(gene)
    rows.append(R(gene, "hypertension", "plof", "replication", good("hypertension"), rep_p))
    return rows


# ---------------- discovery filter ----------------
@pytest.mark.parametrize("trait", sorted(L_PANEL))
def test_discovery_direction_of_benefit_per_trait(trait, cfg):
    rows = [R("GOOD", trait, "plof", "discovery", good(trait), 1e-9),
            R("HARM", trait, "plof", "discovery", bad(trait), 1e-9)]
    got = stats.discovery_hits(T(rows), cfg)
    assert list(got["gene"]) == ["GOOD"], "trait %s: harmful-direction gene must not be a discovery" % trait


def test_discovery_threshold_is_strict_and_not_looser_than_ledger(cfg):
    rows = [R("G1", "ldl", "plof", "discovery", good("ldl"), 1.5e-7),
            R("G2", "ldl", "plof", "discovery", good("ldl"), L_DISCOVERY_P),
            R("G3", "ldl", "plof", "discovery", good("ldl"), 2.0e-7),
            R("G4", "ldl", "plof", "discovery", good("ldl"), 1e-6)]
    got = set(stats.discovery_hits(T(rows), cfg)["gene"])
    assert "G1" in got
    assert not got.intersection(["G2", "G3", "G4"]), got


def test_discovery_only_plof_discovery_panel_traits(cfg):
    rows = [R("SYNONLY", "ldl", "syn", "discovery", good("ldl"), 1e-12),
            R("DMISONLY", "ldl", "dmis", "discovery", good("ldl"), 1e-12),
            R("REPONLY", "ldl", "plof", "replication", good("ldl"), 1e-12),
            R("EURONLY", "ldl", "plof", "discovery_eur", good("ldl"), 1e-12),
            R("TRADEOFF", "type_2_diabetes", "plof", "discovery", good("type_2_diabetes"), 1e-12),
            R("TG", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-12)]
    assert len(stats.discovery_hits(T(rows), cfg)) == 0
    assert len(tiering.build_hits(T(rows), cfg)) == 0


# ---------------- replication rule ----------------
@pytest.mark.parametrize("p2,expect", [(0.02, True), (0.08, True), (0.0999, True),
                                       (0.10, False), (0.12, False), (0.9, False)])
def test_replication_one_sided_boundary(p2, expect, cfg):
    tm = tier_map(sbp_replicated_gene(rep_p=p2), cfg)
    assert oracle_replicates(good("hypertension"), p2, benefit_sign("hypertension")) == expect
    assert tm[("GA", "systolic_bp")] == ("A" if expect else "D")


def test_replication_opposite_direction_never_replicates(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA")
    rows.append(R("GA", "hypertension", "plof", "replication", bad("hypertension"), 1e-12))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "D"


@pytest.mark.parametrize("trait", L_NO_PROXY)
def test_traits_without_declared_proxy_cap_at_tier_B(trait, cfg):
    rows = hit_rows("GA", trait) + null_tradeoffs("GA")
    # even a strong same-trait replication row must not lift the tier above B (Amendment 1)
    rows.append(R("GA", trait, "plof", "replication", good(trait), 1e-12))
    assert tier_map(rows, cfg)[("GA", trait)] == "B"


@pytest.mark.parametrize("trait", sorted(L_PROXIES))
def test_declared_proxy_can_reach_tier_A(trait, cfg):
    px = L_PROXIES[trait]
    rows = hit_rows("GA", trait) + null_tradeoffs("GA")
    rows.append(R("GA", px, "plof", "replication", good(px), 0.01))
    assert tier_map(rows, cfg)[("GA", trait)] == "A"


def test_replication_proxy_cross_wiring_rejected(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA")
    rows.append(R("GA", "obesity", "plof", "replication", good("obesity"), 1e-12))
    rows.append(R("GA", "hypercholesterolemia", "plof", "replication", good("ldl"), 1e-12))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A"


def test_replication_is_per_gene(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA")
    rows.append(R("GB", "hypertension", "plof", "replication", good("hypertension"), 1e-12))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A"


@pytest.mark.parametrize("src", L_UKB_SOURCES)
def test_ukb_overlapping_source_cannot_yield_tier_A(src, cfg):
    """Replication rows whose source is a UK Biobank exome resource are not independent."""
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(R("GA", "hypertension", "plof", "replication", good("hypertension"), 1e-12, source=src))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A", "source %s counted as independent" % src


@pytest.mark.parametrize("cohort", ["discovery", "discovery_eur"])
def test_within_ukb_cohort_rows_never_count_as_replication(cohort, cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(R("GA", "hypertension", "plof", cohort, good("hypertension"), 1e-12))
    rows.append(R("GA", "systolic_bp", "plof", "discovery_eur", good("systolic_bp"), 1e-12))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A"


# ---------------- adverse trade-off screen ----------------
@pytest.mark.parametrize("outcome", L_TRADEOFF)
def test_adverse_each_outcome_harmful_direction_gives_tier_C(outcome, cfg):
    rows = drop_rows(sbp_replicated_gene(), "GA", outcome, "discovery")
    rows.append(R("GA", outcome, "plof", "discovery", +0.8, 1e-4))   # higher risk, p far below 0.05/9
    assert oracle_adverse(+0.8, 1e-4)
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "C"


@pytest.mark.parametrize("outcome", L_TRADEOFF)
def test_protective_tradeoff_direction_is_not_adverse(outcome, cfg):
    rows = drop_rows(sbp_replicated_gene(), "GA", outcome, "discovery")
    rows.append(R("GA", outcome, "plof", "discovery", -0.8, 1e-9))   # lower risk: beneficial
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "A"


@pytest.mark.parametrize("p", [0.02, 0.05, 0.2])
def test_adverse_threshold_not_looser_than_ledger(p, cfg):
    # harmful-direction signals above 0.05/9 (two-sided) are not adverse
    assert not oracle_adverse(+0.8, p)
    rows = drop_rows(sbp_replicated_gene(), "GA", "coronary_disease", "discovery")
    rows.append(R("GA", "coronary_disease", "plof", "discovery", +0.8, p))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "A", "p=%s wrongly adverse" % p


def test_adverse_and_failed_replication_is_not_silently_A_or_B(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("replication",))
    rows.append(R("GA", "coronary_disease", "plof", "discovery", +0.8, 1e-4))
    rows.append(R("GA", "hypertension", "plof", "replication", bad("hypertension"), 0.5))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] in ("C", "D")   # ledger silent on C vs D


@pytest.mark.xfail(strict=False, reason="ledger silent: a gene with 0 of 9 trade-offs screened is reported "
                                          "Tier A (see reviews/review-1.md)")
def test_tier_A_requires_tradeoff_panel_to_have_been_screened(cfg):
    rows = hit_rows("GA", "systolic_bp")
    rows.append(R("GA", "hypertension", "plof", "replication", good("hypertension"), 0.01))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A"


# ---------------- tier definitions ----------------
def test_tier_definitions_all_four_letters(cfg):
    rows = []
    rows += sbp_replicated_gene("TA")                                     # A
    rows += hit_rows("TB", "fluid_intelligence") + null_tradeoffs("TB")   # B
    rows += drop_rows(sbp_replicated_gene("TC"), "TC", "coronary_disease", "discovery")
    rows.append(R("TC", "coronary_disease", "plof", "discovery", +0.9, 1e-5))          # C
    rows += hit_rows("TD", "systolic_bp") + null_tradeoffs("TD")
    rows.append(R("TD", "hypertension", "plof", "replication", bad("hypertension"), 0.4))   # D
    tm = tier_map(rows, cfg)
    assert tm[("TA", "systolic_bp")] == "A"
    assert tm[("TB", "fluid_intelligence")] == "B"
    assert tm[("TC", "systolic_bp")] == "C"
    assert tm[("TD", "systolic_bp")] == "D"


def test_synonymous_only_gene_never_receives_a_tier(cfg):
    rows = [R("SYNGENE", "systolic_bp", "syn", "discovery", good("systolic_bp"), 1e-15)]
    assert len(tiering.build_hits(T(rows), cfg)) == 0
