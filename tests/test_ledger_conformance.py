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
