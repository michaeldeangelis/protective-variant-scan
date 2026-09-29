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
# Amendment 1 clarifications C1-C4 (experiments.md, written before any gene-level result was opened)
L_MIN_TRADEOFFS_SCREENED = 5                          # C1: Tier A needs >= 5 of 9 trade-off outcomes screened
L_QUALIFYING_DOMAINS = set(["cognitive", "physical"])  # C2: SBP counts as physical; LDL, BMI, lifespan do not
L_LIPID_GENES = set("""PCSK9 ANGPTL3 ANGPTL4 ANGPTL8 APOC3 APOA5 APOB APOE LDLR LDLRAP1 LPL LPA CETP LIPC MTTP
NPC1L1 ABCG5 ABCG8 HMGCR SORT1 LIPG GPIHBP1 LMF1 ANGPTL1""".split())   # C3: frozen at commit 4073b89
L_DISCOVERY_P_ACCEPTED = 1.9e-7                       # C4: accepted as implemented (not 1.923e-7)
L_INDEP_TOKENS = ("finngen", "all_of_us", "allofus", "synthetic_replication")   # C4 allow list (+ fixture token)


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
    """C4: the trade-off screen is one-sided in the harmful direction, alpha/9."""
    if beta * adverse_sign <= 0:
        return False
    return (p / 2.0) < L_ADVERSE_P


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
    # C4 accepts 1.9e-7 (the value printed in the ledger), which is stricter than 0.05/(20000*13) = 1.923e-7
    assert thr == L_DISCOVERY_P_ACCEPTED, thr
    assert thr <= L_DISCOVERY_P


def test_config_replication_rule(cfg_yaml):
    r = cfg_yaml["replication"]
    assert r["one_sided_p"] == L_REPL_ONE_SIDED_P
    assert dict((k, v["trait"]) for k, v in r["proxies"].items()) == L_PROXIES
    for k, v in r["proxies"].items():
        assert v["benefit"] == L_PANEL[k][1], "proxy %s direction must equal parent trait" % k


def test_config_adverse_threshold(cfg_yaml):
    assert cfg_yaml["tradeoff_screen"]["alpha"] == L_ALPHA  # divided by N_tradeoff = 9 in code
    assert cfg_yaml["tradeoff_screen"]["min_tradeoffs_screened"] == L_MIN_TRADEOFFS_SCREENED   # C1


def test_config_replication_source_allow_and_deny_lists(cfg_yaml):
    r = cfg_yaml["replication"]
    allow = [t.lower() for t in r["independent_sources"]]
    deny = [t.lower() for t in r["ukb_overlapping_sources"]]
    assert "finngen" in allow and "all_of_us" in allow                       # C4
    for tok in ("genebass", "azphewas", "astrazeneca", "regeneron", "ukb", "uk_biobank"):
        assert tok in deny, tok
    for a in allow:
        assert not any(d in a for d in deny), "allow token %r contains a UKB token" % a


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
    assert lipid == L_LIPID_GENES, lipid ^ L_LIPID_GENES        # C3: frozen list
    assert len(v["lipid_pathway_genes"]) == len(lipid), "duplicate entries"


def test_config_loads_via_schema_and_derives_thresholds():
    from protscan.schema import load_config
    c = load_config(CONFIG_PATH)
    assert len(c.panel) == 13 and len(c.tradeoff) == 9
    assert math.isclose(c.tradeoff_p, L_ADVERSE_P, rel_tol=1e-12)
    assert c.discovery_p == L_DISCOVERY_P_ACCEPTED
    assert c.min_tradeoffs_screened == L_MIN_TRADEOFFS_SCREENED
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


# ---------------- controls (unit level) ----------------
CHI2_MEDIAN = float(chi2.ppf(0.5, 1))


def baseline_null(n_genes=150, seed=7):
    """Null discovery table (plof/dmis/syn, N(0,1) z-scores) plus the lipid positive controls."""
    rng = np.random.default_rng(seed)
    genes = ["NULLG%04d" % i for i in range(n_genes)]
    traits = sorted(L_PANEL) + list(L_TRADEOFF) + list(CTRL_TRAITS)
    rows = []
    for g in genes:
        for t in traits:
            for m in ("plof", "dmis", "syn"):
                z = float(rng.standard_normal())
                rows.append(R(g, t, m, "discovery", z * 0.1, float(2 * norm.sf(abs(z)))))
    rows.append(R("PCSK9", "ldl", "plof", "discovery", good("ldl"), 1e-30))
    rows.append(R("PCSK9", "coronary_disease", "plof", "discovery", good("coronary_disease"), 1e-6))
    rows.append(R("ANGPTL4", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-20))
    rows.append(R("APOC3", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-20))
    return pd.DataFrame(rows, columns=COLUMNS)


@pytest.fixture(scope="module")
def null_table():
    return baseline_null()


def set_row(df, gene, trait, mask, cohort, beta=None, p=None):
    df = df.copy()
    sel = (df["gene"] == gene) & (df["trait"] == trait) & (df["mask"] == mask) & (df["cohort"] == cohort)
    assert sel.sum() == 1, (gene, trait, mask, cohort)
    if beta is not None:
        df.loc[sel, "beta"] = beta
    if p is not None:
        df.loc[sel, "p"] = p
    return df


def drop_gene(df, gene):
    return df[df["gene"] != gene].copy()


def run_ctrl(df, cfg):
    return controls.run_controls(validate(df), cfg)


def test_controls_valid_on_clean_null_table(null_table, cfg):
    c = run_ctrl(null_table, cfg)
    assert c["negative_synonymous"]["lambda_gc"] < L_LAMBDA_GC_MAX
    assert c["negative_synonymous"]["n_syn_hit_genes"] == 0
    assert c["valid"] is True and c["failed"] == [] and c["not_run"] == []


def test_pcsk9_ldl_wrong_direction_fails(null_table, cfg):
    df = set_row(null_table, "PCSK9", "ldl", "plof", "discovery", beta=bad("ldl"), p=1e-30)
    c = run_ctrl(df, cfg)
    assert c["valid"] is False and c["failed"]


def test_pcsk9_ldl_beneficial_but_below_discovery_threshold_fails(null_table, cfg):
    df = set_row(null_table, "PCSK9", "ldl", "plof", "discovery", p=1e-5)
    c = run_ctrl(df, cfg)
    assert c["valid"] is False and c["failed"]


def test_pcsk9_missing_makes_controls_invalid(null_table, cfg):
    c = run_ctrl(drop_gene(null_table, "PCSK9"), cfg)
    assert c["valid"] is False


def test_pcsk9_coronary_wrong_direction_fails(null_table, cfg):
    df = set_row(null_table, "PCSK9", "coronary_disease", "plof", "discovery", beta=bad("coronary_disease"))
    c = run_ctrl(df, cfg)
    assert c["valid"] is False and c["failed"]


def test_lipid_control_needs_one_of_angptl4_or_apoc3(null_table, cfg):
    one_bad = set_row(null_table, "ANGPTL4", "triglycerides", "plof", "discovery", beta=bad("triglycerides"))
    assert run_ctrl(one_bad, cfg)["valid"] is True            # APOC3 still recovered
    both_bad = set_row(one_bad, "APOC3", "triglycerides", "plof", "discovery", beta=bad("triglycerides"))
    c = run_ctrl(both_bad, cfg)
    assert c["valid"] is False and c["failed"]
    weak = set_row(set_row(null_table, "ANGPTL4", "triglycerides", "plof", "discovery", p=1e-4),
                   "APOC3", "triglycerides", "plof", "discovery", p=1e-4)
    assert run_ctrl(weak, cfg)["valid"] is False              # beneficial but not at discovery threshold


def flat_syn_lambda(df, lam):
    """Set every synonymous p to the value whose chi2 statistic is lam x median(chi2_1): lambda_GC == lam."""
    df = df.copy()
    sel = df["mask"] == "syn"
    df.loc[sel, "p"] = float(chi2.sf(lam * CHI2_MEDIAN, 1))
    df.loc[sel, "beta"] = 0.0
    return df


def test_lambda_gc_boundary(null_table, cfg):
    ok = run_ctrl(flat_syn_lambda(null_table, 1.099), cfg)
    assert ok["negative_synonymous"]["lambda_ok"] is True and ok["valid"] is True
    bad_ = run_ctrl(flat_syn_lambda(null_table, 1.101), cfg)
    assert bad_["negative_synonymous"]["lambda_ok"] is False
    assert bad_["valid"] is False and "negative_synonymous" in bad_["failed"]


def test_synonymous_hit_in_beneficial_direction_fails_control(null_table, cfg):
    df = set_row(null_table, "NULLG0001", "hand_grip_strength", "syn", "discovery",
                 beta=good("hand_grip_strength"), p=1e-9)
    c = run_ctrl(df, cfg)
    assert c["negative_synonymous"]["n_syn_hit_genes"] == 1
    assert c["valid"] is False and "negative_synonymous" in c["failed"]


def test_synonymous_hit_in_harmful_direction_is_not_counted(null_table, cfg):
    df = set_row(null_table, "NULLG0001", "hand_grip_strength", "syn", "discovery",
                 beta=bad("hand_grip_strength"), p=1e-9)
    c = run_ctrl(df, cfg)
    assert c["negative_synonymous"]["n_syn_hit_genes"] == 0


# ---------------- UKB-overlap: name variants a substring denylist can miss ----------------
UKB_NAME_VARIANTS = ["uk-biobank", "az_phewas", "opentargets_gene_burden", "backman_2021_exomes",
                     "UK Biobank 450k exomes", "pan_ukb_burden"]


@pytest.mark.parametrize("src", UKB_NAME_VARIANTS)
def test_ukb_derived_source_name_variants_cannot_yield_tier_A(src, cfg):
    """Ledger: Genebass, AZ and Regeneron/Open Targets are all UKB exomes; none may count as replication."""
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(R("GA", "hypertension", "plof", "replication", good("hypertension"), 1e-12, source=src))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A", "source %r counted as independent" % src


@pytest.mark.parametrize("src", ["finngen_r13", "finngen_r12", "FinnGen_DF13", "all_of_us_aba", "synthetic_replication"])
def test_independent_sources_can_still_yield_tier_A(src, cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(R("GA", "hypertension", "plof", "replication", good("hypertension"), 0.01, source=src))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "A"


# ======================================================================================
# (c) end to end: synthetic fixtures -> python -m protscan run -> independent oracle
# ======================================================================================
import gzip
import importlib.util
import json
import os
import subprocess
import sys

MIN_LIPID = set(["PCSK9", "ANGPTL4", "APOC3", "APOB", "LDLR", "ANGPTL3", "LPL", "LPA"])
N_GENES = 400


def _load_synth():
    spec = importlib.util.spec_from_file_location("make_synthetic_for_conformance", ROOT / "scripts" / "make_synthetic.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def synth():
    return _load_synth()


def write_tables(tables, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(out_dir / (name + ".csv.gz"), index=False)


def run_cli(data_dir, out_json):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    proc = subprocess.run(
        [sys.executable, "-m", "protscan", "run", "--config", str(CONFIG_PATH),
         "--data", str(data_dir), "--out", str(out_json)],
        capture_output=True, text=True, env=env, cwd=str(ROOT), timeout=600)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return json.loads(Path(out_json).read_text())


def all_rows(tables):
    return pd.concat(list(tables.values()), ignore_index=True)


def is_ukb_source(src):
    s = str(src).lower()
    return any(tok in s for tok in L_UKB_SOURCES)


def is_independent_source(src):
    """C4: allow-listed (finngen, all_of_us) and no UK Biobank token."""
    s = str(src).lower()
    return any(tok in s for tok in L_INDEP_TOKENS) and not is_ukb_source(s)


def oracle(df):
    """Independent re-implementation of the ledger + Amendment 1 on a normalized table.

    Returns (controls_valid, tiers, verdict) where tiers maps (gene, trait) -> letter for discovery hits.
    """
    d = df.query("cohort == 'discovery' and mask == 'plof'")
    hits = []
    for r in d.itertuples():
        if r.trait in L_PANEL and oracle_discovery(r.beta, r.p, L_PANEL[r.trait][1]):
            hits.append((r.gene, r.trait))
    rep = df.query("cohort == 'replication' and mask == 'plof'")
    rep = rep[[is_independent_source(s) for s in rep["source"]]]
    repd = dict(((r.gene, r.trait), r) for r in rep.itertuples())
    adverse_genes = set()
    tox = df.query("mask == 'plof' and cohort in ['discovery', 'replication']")
    for r in tox.itertuples():
        if r.trait in L_TRADEOFF and oracle_adverse(r.beta, r.p):
            adverse_genes.add(r.gene)
    screened = dict()
    for r in tox.itertuples():
        if r.trait in L_TRADEOFF:
            screened.setdefault(r.gene, set()).add(r.trait)
    tiers = dict()
    for g, t in hits:
        if g in adverse_genes:
            tiers[(g, t)] = "C"
            continue
        if t not in L_PROXIES:
            tiers[(g, t)] = "B"
            continue
        row = None
        for cand in (t, L_PROXIES[t]):
            if (g, cand) in repd:
                row = repd[(g, cand)]
                break
        if row is None:
            tiers[(g, t)] = "B"
        elif oracle_replicates(row.beta, row.p, benefit_sign(row.trait)):
            tiers[(g, t)] = "A" if len(screened.get(g, ())) >= L_MIN_TRADEOFFS_SCREENED else "B"   # C1
        else:
            tiers[(g, t)] = "D"
    # controls
    syn = df.query("cohort == 'discovery' and mask == 'syn'")
    lam = float(np.median(chi2.isf(np.clip(syn["p"].to_numpy(), 1e-300, 1.0), 1)) / CHI2_MEDIAN) if len(syn) else float("inf")
    syn_hits = [r for r in syn.itertuples() if r.trait in L_PANEL and oracle_discovery(r.beta, r.p, L_PANEL[r.trait][1])]

    def rowof(g, t):
        x = d[(d["gene"] == g) & (d["trait"] == t)]
        return x.iloc[0] if len(x) else None

    ok = lam < L_LAMBDA_GC_MAX and len(syn_hits) == 0
    x = rowof("PCSK9", "ldl")
    ok = ok and x is not None and oracle_discovery(x.beta, x.p, -1)
    x = rowof("PCSK9", "coronary_disease")
    ok = ok and x is not None and x.beta < 0
    tg = [rowof(g, "triglycerides") for g in ("ANGPTL4", "APOC3")]
    ok = ok and any(y is not None and y.beta < 0 and y.p < L_DISCOVERY_P for y in tg)
    # verdict
    dm = df.query("cohort == 'discovery' and mask == 'dmis'")
    dmd = dict(((r.gene, r.trait), r.beta) for r in dm.itertuples())
    qual = []
    for (g, t), letter in tiers.items():
        if g in MIN_LIPID or L_PANEL[t][0] not in ("cognitive", "physical"):
            continue
        if len(screened.get(g, ())) < L_MIN_TRADEOFFS_SCREENED:      # C9: LEAD/PASS need >= 5 outcomes screened
            continue
        if (g, t) in dmd and oracle_beneficial(dmd[(g, t)], L_PANEL[t][1]):
            qual.append(letter)
    if not ok:
        verdict = "KILL"
    elif "A" in qual:
        verdict = "PASS"
    elif "B" in qual:
        verdict = "LEAD"
    else:
        verdict = "KILL"
    return ok, tiers, verdict


def result_tiers(res):
    out = dict()
    for letter in "ABCD":
        for r in res["tiers"][letter]:
            out[(r["gene"], r["trait"])] = letter
    return out




@pytest.fixture(scope="module")
def scenario_runs(synth, tmp_path_factory):
    """Each synthetic scenario written to disk and run through the real CLI once."""
    runs = dict()
    for name in ("pass", "lead", "nolead", "broken_positive", "broken_lambda", "broken_syn_hit"):
        d = tmp_path_factory.mktemp("scn_" + name)
        tables = synth.make_synthetic(name, n_genes=N_GENES)
        write_tables(tables, d)
        res = run_cli(d, d / "out" / "protective-scan.json")
        runs[name] = (tables, res, d)
    return runs


@pytest.mark.parametrize("name,verdict", [("pass", "PASS"), ("lead", "LEAD"), ("nolead", "KILL"),
                                          ("broken_positive", "KILL"), ("broken_lambda", "KILL"),
                                          ("broken_syn_hit", "KILL")])
def test_scenario_verdict_matches_ledger_rules(scenario_runs, name, verdict):
    tables, res, _ = scenario_runs[name]
    assert res["verdict"]["verdict"] == verdict
    assert res["verdict"]["verdict"] in L_VERDICTS
    ok, tiers, o_verdict = oracle(all_rows(tables))
    assert o_verdict == verdict, "oracle and hand expectation disagree (test bug)"
    assert res["controls"]["valid"] == bool(ok)
    assert result_tiers(res) == tiers, "pipeline tiers differ from independent oracle"


def test_results_json_and_report_contain_required_sections(scenario_runs):
    _, res, d = scenario_runs["pass"]
    for k in ("controls", "rungs", "tiers", "verdict"):
        assert k in res
    c = res["controls"]
    assert "lambda_gc" in c["negative_synonymous"] and "n_syn_hit_genes" in c["negative_synonymous"]
    assert len(c["positive"]) >= 3
    assert set(res["rungs"]) == set(["trivial", "simplest", "incumbent", "candidate"])
    assert set(res["tiers"]) == set("ABCD")
    text = (d / "out" / "report.md").read_text()
    for needle in ("Verdict", "lambda_GC", "Tier A", "Tier B", "Tier C", "Tier D", "trivial", "simplest",
                   "incumbent", "candidate", "SYNTHETIC"):
        assert needle in text, needle


def test_trivial_rung_is_the_synonymous_mask(scenario_runs):
    _, res, _ = scenario_runs["broken_syn_hit"]
    assert res["rungs"]["trivial"]["n_genes"] == 1
    assert res["controls"]["negative_synonymous"]["syn_hit_genes"] == ["SYNSYN1"]
    _, res0, _ = scenario_runs["pass"]
    assert res0["rungs"]["trivial"]["n_genes"] == 0


def test_lead_is_never_reported_as_pass(scenario_runs):
    _, res, d = scenario_runs["lead"]
    v = res["verdict"]
    assert v["verdict"] == "LEAD" and v["pass_genes"] == [] and v["lead_genes"] == ["SYNLEAD1"]
    assert result_tiers(res)[("SYNLEAD1", "fluid_intelligence")] == "B"
    text = (d / "out" / "report.md").read_text()
    assert "## Verdict: LEAD" in text and "UNREPLICATED" in text
    assert "## Verdict: PASS" not in text


def test_pass_scenario_tiers_by_hand(scenario_runs):
    _, res, _ = scenario_runs["pass"]
    tm = result_tiers(res)
    assert tm[("SYNPASS1", "systolic_bp")] == "A"
    assert tm[("SYNLEAD1", "fluid_intelligence")] == "B"
    assert tm[("SYNADV1", "systolic_bp")] == "C"
    assert tm[("SYNFAIL1", "systolic_bp")] == "D"
    assert tm[("SYNMASK1", "fev1")] == "B"
    assert res["verdict"]["pass_genes"] == ["SYNPASS1"]      # SYNMASK1 excluded: dmis direction disagrees
    assert "SYNMASK1" not in res["verdict"]["lead_genes"]
    assert "PCSK9" not in res["verdict"]["pass_genes"]       # lipid-pathway gene, Tier A on LDL


# ---------------- mutation tests: break one thing, check verdict and oracle agreement ----------------
def setrow(df, gene, trait, mask, cohort, z=None, beta=None, p=None, source=None):
    sel = (df["gene"] == gene) & (df["trait"] == trait) & (df["mask"] == mask) & (df["cohort"] == cohort)
    assert sel.sum() == 1, (gene, trait, mask, cohort, int(sel.sum()))
    i = df.index[sel][0]
    if z is not None:
        df.loc[i, "beta"] = z * df.loc[i, "se"]
        df.loc[i, "p"] = float(2 * norm.sf(abs(z)))
    if beta is not None:
        df.loc[i, "beta"] = beta
    if p is not None:
        df.loc[i, "p"] = p
    if source is not None:
        df.loc[i, "source"] = source


def zgood(trait, mag):
    return benefit_sign(trait) * abs(mag)


def zbad(trait, mag):
    return -benefit_sign(trait) * abs(mag)


@pytest.fixture(scope="module")
def base_frames(synth):
    return dict((n, all_rows(synth.make_synthetic(n, n_genes=N_GENES)).copy()) for n in ("pass", "lead"))


def mutated_run(base, edit, tmp_path):
    from protscan.run import run_pipeline
    df = base.copy()
    edit(df)
    d = tmp_path / "data"
    d.mkdir()
    df.to_csv(d / "burden_mutated.csv.gz", index=False)
    res = run_pipeline(CONFIG_PATH, d, tmp_path / "out" / "protective-scan.json")
    return res, df


def lam_flat(df, lam):
    sel = df["mask"] == "syn"
    df.loc[sel, "p"] = float(chi2.sf(lam * CHI2_MEDIAN, 1))
    df.loc[sel, "beta"] = 0.0


def edit_flip_pcsk9_ldl(df):
    setrow(df, "PCSK9", "ldl", "plof", "discovery", z=zbad("ldl", 15))


def edit_pcsk9_ldl_weak(df):
    setrow(df, "PCSK9", "ldl", "plof", "discovery", p=1e-5)


def edit_drop_tg_controls(df):
    df.drop(df.index[df["gene"].isin(["ANGPTL4", "APOC3"])], inplace=True)


def edit_flip_both_tg(df):
    setrow(df, "ANGPTL4", "triglycerides", "plof", "discovery", z=zbad("triglycerides", 9))
    setrow(df, "APOC3", "triglycerides", "plof", "discovery", z=zbad("triglycerides", 12))


def edit_pcsk9_cad_flip(df):
    setrow(df, "PCSK9", "coronary_disease", "plof", "discovery", z=zbad("coronary_disease", 4))


def edit_lambda_high(df):
    lam_flat(df, 1.101)


def edit_syn_hit(df):
    setrow(df, "SYNPASS1", "hand_grip_strength", "syn", "discovery", z=zgood("hand_grip_strength", 7))


@pytest.mark.parametrize("edit", [edit_flip_pcsk9_ldl, edit_pcsk9_ldl_weak, edit_drop_tg_controls, edit_flip_both_tg,
                                  edit_pcsk9_cad_flip, edit_lambda_high, edit_syn_hit],
                         ids=lambda f: f.__name__)
def test_breaking_any_control_yields_kill_even_with_tier_A_gene(edit, base_frames, tmp_path):
    res, df = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["controls"]["valid"] is False
    assert res["verdict"]["verdict"] == "KILL", res["verdict"]
    assert res["verdict"]["controls_valid"] is False
    ok, tiers, o_verdict = oracle(df)
    assert not ok and o_verdict == "KILL"
    assert ("SYNPASS1", "systolic_bp") in result_tiers(res)     # a Tier-A style gene existed; KILL still wins
    assert "PASS" != res["verdict"]["verdict"]


def test_lambda_gc_just_below_limit_still_passes(base_frames, tmp_path):
    res, df = mutated_run(base_frames["pass"], lambda x: lam_flat(x, 1.099), tmp_path)
    assert res["controls"]["negative_synonymous"]["lambda_gc"] < L_LAMBDA_GC_MAX
    assert res["verdict"]["verdict"] == "PASS"


def test_synonymous_hit_harmful_direction_does_not_kill(base_frames, tmp_path):
    def edit(df):
        setrow(df, "SYNPASS1", "hand_grip_strength", "syn", "discovery", z=zbad("hand_grip_strength", 7))
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["controls"]["valid"] is True and res["verdict"]["verdict"] == "PASS"


def _rename_gene(old, new):
    def edit(df):
        df.loc[df["gene"] == old, "gene"] = new
    return edit


def _sbp_rep(p):
    def edit(df):
        setrow(df, "SYNPASS1", "hypertension", "plof", "replication", beta=zgood("hypertension", 1) * 0.1, p=p)
    return edit


def edit_ukb_relabel(df):
    sel = df["cohort"] == "replication"
    df.loc[sel, "source"] = "genebass_hypertension_lookup"


def edit_no_replication(df):
    df.drop(df.index[df["cohort"] == "replication"], inplace=True)


def edit_dmis_flip(df):
    setrow(df, "SYNPASS1", "systolic_bp", "dmis", "discovery", z=zbad("systolic_bp", 2))


def edit_dmis_missing(df):
    df.drop(df.index[(df["gene"] == "SYNPASS1") & (df["mask"] == "dmis")], inplace=True)


def edit_adverse(df):
    setrow(df, "SYNPASS1", "coronary_disease", "plof", "discovery", beta=0.9, p=1e-4)


def edit_sbp_p_just_above_threshold(df):
    setrow(df, "SYNPASS1", "systolic_bp", "plof", "discovery", p=2.0e-7)


def edit_sbp_harmful(df):
    setrow(df, "SYNPASS1", "systolic_bp", "plof", "discovery", z=zbad("systolic_bp", 9))


PASS_BLOCKERS = [
    ("ukb_relabelled_replication", edit_ukb_relabel),
    ("no_replication_cohort", edit_no_replication),
    ("dmis_direction_disagrees", edit_dmis_flip),
    ("dmis_missing", edit_dmis_missing),
    ("adverse_tradeoff_coronary", edit_adverse),
    ("replication_one_sided_p_0.055", _sbp_rep(0.11)),
    ("replication_wrong_direction", None),
    ("discovery_p_2e-7", edit_sbp_p_just_above_threshold),
    ("discovery_harmful_direction", edit_sbp_harmful),
    ("lipid_pathway_gene_APOB", _rename_gene("SYNPASS1", "APOB")),
    ("lipid_pathway_gene_LDLR", _rename_gene("SYNPASS1", "LDLR")),
]


def edit_rep_wrong_direction(df):
    setrow(df, "SYNPASS1", "hypertension", "plof", "replication", beta=zbad("hypertension", 1) * 0.1, p=1e-9)


@pytest.mark.parametrize("name,edit", PASS_BLOCKERS, ids=[n for n, _ in PASS_BLOCKERS])
def test_each_single_defect_blocks_pass_and_degrades_to_lead(name, edit, base_frames, tmp_path):
    if edit is None:
        edit = edit_rep_wrong_direction
    res, df = mutated_run(base_frames["pass"], edit, tmp_path)
    ok, tiers, o_verdict = oracle(df)
    assert res["controls"]["valid"] is True, res["controls"]["failed"]
    assert res["verdict"]["verdict"] != "PASS", "defect %s still gave PASS" % name
    assert res["verdict"]["verdict"] == "LEAD"           # SYNLEAD1 (cognitive, Tier B) remains
    assert res["verdict"]["verdict"] == o_verdict
    assert result_tiers(res) == tiers


def test_replication_one_sided_boundary_end_to_end(base_frames, tmp_path):
    res, _ = mutated_run(base_frames["pass"], _sbp_rep(0.0999), tmp_path)
    assert res["verdict"]["verdict"] == "PASS"


def test_same_trait_replication_for_cognitive_trait_still_lead_not_pass(base_frames, tmp_path):
    def edit(df):
        extra = df[(df["gene"] == "SYNLEAD1") & (df["trait"] == "fluid_intelligence")
                   & (df["mask"] == "plof") & (df["cohort"] == "discovery")].copy()
        extra["cohort"] = "replication"
        extra["source"] = "finngen_r13"
        extra["p"] = 1e-12
        df.loc[df.index.max() + 1] = extra.iloc[0]
    res, df = mutated_run(base_frames["lead"], edit, tmp_path)
    assert res["verdict"]["verdict"] == "LEAD" and res["verdict"]["pass_genes"] == []
    assert result_tiers(res)[("SYNLEAD1", "fluid_intelligence")] == "B"


# ---------------- verdict logic: exhaustive over tier/qualification combinations ----------------
def _hits_frame(spec):
    """spec: list of (tier, qualifies). One distinct gene per entry."""
    rows = [dict(gene="G%d" % i, tier=t, qualifies=q) for i, (t, q) in enumerate(spec)]
    return pd.DataFrame(rows, columns=["gene", "tier", "qualifies"])


def _all_specs():
    keys = [(t, q) for t in "ABCD" for q in (True, False)]
    for mask in range(1 << len(keys)):
        yield [keys[i] for i in range(len(keys)) if mask >> i & 1]


CTRL_OK = dict(valid=True, failed=[], not_run=[])
CTRL_FAIL = dict(valid=False, failed=["pcsk9_ldl_lower"], not_run=[])
CTRL_NOTRUN = dict(valid=False, failed=[], not_run=["negative_synonymous"])


def test_decide_is_exhaustively_consistent_with_amendment_1():
    from protscan.run import decide
    n = 0
    for spec in _all_specs():
        h = _hits_frame(spec)
        has_a = ("A", True) in spec
        has_b = ("B", True) in spec
        v = decide(CTRL_OK, h)["verdict"]
        assert v == ("PASS" if has_a else "LEAD" if has_b else "KILL"), spec
        assert decide(CTRL_FAIL, h)["verdict"] == "KILL", spec
        assert decide(CTRL_NOTRUN, h)["verdict"] == "KILL", spec
        n += 1
    assert n == 256


def test_decide_never_lists_tier_B_gene_as_pass_gene():
    from protscan.run import decide
    out = decide(CTRL_OK, _hits_frame([("B", True), ("B", True), ("C", True), ("D", True), ("A", False)]))
    assert out["verdict"] == "LEAD" and out["pass_genes"] == [] and len(out["lead_genes"]) == 2


# ---------------- fuzz: implementation tiers vs independent oracle ----------------
def _fuzz_table(seed, n_genes=60):
    rng = np.random.default_rng(seed)
    panel = sorted(L_PANEL)
    proxies_of = L_PROXIES
    indep = ["finngen_r13", "all_of_us_aba", "FinnGen_R12"]
    ukb = ["genebass", "azphewas_v1", "regeneron_rgc", "UKB_lookup", "opentargets", "mystery_biobank"]
    rows = []
    for i in range(n_genes):
        g = "F%03d" % i
        for t in rng.choice(panel, size=int(rng.integers(1, 4)), replace=False):
            sgn = 1 if rng.random() < 0.7 else -1
            p = float(10 ** rng.uniform(-9, -6))                  # straddles 1.9e-7
            while 1.89e-7 <= p <= 1.93e-7:                          # C4: 1.9e-7 accepted; formula gives 1.923e-7
                p = float(10 ** rng.uniform(-9, -6))
            rows.append(R(g, t, "plof", "discovery", sgn * benefit_sign(t) * 0.3, p))
            rows.append(R(g, t, "dmis", "discovery", float(rng.normal()), 0.4))
            if t in proxies_of and rng.random() < 0.8:
                px = proxies_of[t]
                s2 = 1 if rng.random() < 0.7 else -1
                p2 = float(rng.choice([0.001, 0.03, 0.0999, 0.1, 0.11, 0.4, 0.9]))
                src = str(rng.choice(indep + ukb))
                rows.append(R(g, px, "plof", "replication", s2 * benefit_sign(px) * 0.2, p2, source=src))
        for out in L_TRADEOFF:
            for cohort in ("discovery", "replication"):
                if rng.random() < 0.6:
                    p3 = float(10 ** rng.uniform(-6, 0))
                    rows.append(R(g, out, "plof", cohort, float(rng.normal()), p3))
    return T(rows)


@pytest.mark.parametrize("seed", range(25))
def test_fuzz_tiers_match_independent_oracle(seed, cfg):
    df = _fuzz_table(seed)
    hits = tiering.build_hits(df, cfg)
    got = dict(((r.gene, r.trait), r.tier) for r in hits.itertuples())
    _, want, _ = oracle(df)
    assert got == want


# ---------------- documented gaps (xfail, non-strict): see reviews/review-1.md ----------------
@pytest.mark.xfail(strict=False, reason="R-1: negative control has no minimum coverage; 3 synonymous rows validate the pipeline")
def test_negative_control_must_cover_the_discovery_universe(cfg):
    rows = [R("PCSK9", "ldl", "plof", "discovery", good("ldl"), 1e-30),
            R("PCSK9", "coronary_disease", "plof", "discovery", good("coronary_disease"), 1e-4),
            R("APOC3", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-30),
            R("PCSK9", "ldl", "syn", "discovery", 0.001, 0.6),
            R("G1", "ldl", "syn", "discovery", 0.001, 0.5),
            R("G2", "bmi", "syn", "discovery", 0.001, 0.4)]
    c = controls.run_controls(T(rows), cfg)
    assert c["valid"] is False

@pytest.mark.xfail(strict=False, reason="R-3: any --config is accepted; verdict not tied to the preregistered constants")
def test_edited_config_cannot_silently_produce_a_verdict(synth, tmp_path):
    text = CONFIG_PATH.read_text().replace("syn_hits_max: 0", "syn_hits_max: 3")
    assert "syn_hits_max: 3" in text
    cfg_copy = tmp_path / "loose.yaml"
    cfg_copy.write_text(text)
    d = tmp_path / "data"
    write_tables(synth.make_synthetic("broken_syn_hit", n_genes=N_GENES), d)
    from protscan.run import run_pipeline
    try:
        res = run_pipeline(cfg_copy, d, tmp_path / "o" / "r.json")
    except Exception:
        return
    assert res["verdict"]["verdict"] == "KILL" or res.get("config_matches_ledger") is False


# ======================================================================================
# (d) Amendment 1 clarifications C1-C4
# ======================================================================================
def sbp_gene_with_screen(gene, disc_outcomes, rep_outcomes=()):
    """SBP hit replicated via hypertension; only the listed trade-off outcomes have plof rows."""
    rows = hit_rows(gene, "systolic_bp")
    rows.append(R(gene, "hypertension", "plof", "replication", good("hypertension"), 0.01))
    for t in disc_outcomes:
        rows.append(R(gene, t, "plof", "discovery", 0.01, 0.5))
    for t in rep_outcomes:
        rows.append(R(gene, t, "plof", "replication", 0.01, 0.5))
    return rows


@pytest.mark.parametrize("n,letter", [(0, "B"), (1, "B"), (4, "B"), (5, "A"), (6, "A"), (9, "A")])
def test_c1_tier_A_needs_at_least_five_tradeoffs_screened(n, letter, cfg):
    rows = sbp_gene_with_screen("GA", L_TRADEOFF[:n])
    assert L_MIN_TRADEOFFS_SCREENED == 5
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == letter


def test_c1_screening_is_the_union_of_outcomes_across_cohorts(cfg):
    union = sbp_gene_with_screen("GA", L_TRADEOFF[:3], L_TRADEOFF[3:5])      # 5 distinct outcomes
    assert tier_map(union, cfg)[("GA", "systolic_bp")] == "A"
    same = sbp_gene_with_screen("GB", L_TRADEOFF[:3], L_TRADEOFF[:3])        # 3 distinct outcomes, twice
    assert tier_map(same, cfg)[("GB", "systolic_bp")] == "B"


def test_c1_only_plof_rows_count_as_screened(cfg):
    rows = sbp_gene_with_screen("GA", L_TRADEOFF[:4])
    rows.append(R("GA", L_TRADEOFF[4], "dmis", "discovery", 0.01, 0.5))       # not the plof mask
    rows.append(R("GA", L_TRADEOFF[5], "syn", "discovery", 0.01, 0.5))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "B"


def test_c1_unscreened_would_be_A_is_flagged_and_adverse_still_wins(cfg):
    hits = tiering.build_hits(T(sbp_gene_with_screen("GA", L_TRADEOFF[:4])), cfg)
    assert list(hits["tier"]) == ["B"] and bool(hits["tradeoff_unscreened"].iloc[0]) is True
    rows = sbp_gene_with_screen("GB", L_TRADEOFF[1:4])
    rows.append(R("GB", L_TRADEOFF[0], "plof", "discovery", 0.9, 1e-5))      # harmful and significant
    assert tier_map(rows, cfg)[("GB", "systolic_bp")] == "C"


def test_c1_end_to_end_unscreened_gene_cannot_give_pass(base_frames, tmp_path):
    def edit(df):
        keep = L_TRADEOFF[:4]
        sel = (df["gene"] == "SYNPASS1") & df["trait"].isin(L_TRADEOFF) & (~df["trait"].isin(keep))
        df.drop(df.index[sel], inplace=True)
    res, df = mutated_run(base_frames["pass"], edit, tmp_path)
    assert result_tiers(res)[("SYNPASS1", "systolic_bp")] == "B"
    assert res["verdict"]["verdict"] == "LEAD" and res["verdict"]["pass_genes"] == []
    _, tiers, o_verdict = oracle(df)
    assert tiers[("SYNPASS1", "systolic_bp")] == "B" and o_verdict == "LEAD"


# ---------------- C2: only cognitive/physical traits qualify; PASS is reachable only through SBP ----------------
def qualifies_for(trait, cfg):
    rows = hit_rows("GA", trait) + null_tradeoffs("GA")
    if trait in L_PROXIES:
        px = L_PROXIES[trait]
        rows.append(R("GA", px, "plof", "replication", good(px), 0.01))
    hits = tiering.build_hits(T(rows), cfg)
    assert len(hits) == 1
    return bool(hits["qualifies"].iloc[0]), hits["tier"].iloc[0]


@pytest.mark.parametrize("trait", sorted(L_PANEL))
def test_c2_qualifying_domains(trait, cfg):
    q, tier = qualifies_for(trait, cfg)
    assert q == (L_PANEL[trait][0] in L_QUALIFYING_DOMAINS), trait
    assert tier == ("A" if trait in L_PROXIES else "B")


def test_c2_replicated_ldl_or_bmi_gene_cannot_give_pass(base_frames, tmp_path):
    for trait, px in (("bmi", "obesity"), ("ldl", "hypercholesterolemia")):
        def edit(df, trait=trait, px=px):
            setrow(df, "SYNPASS1", "systolic_bp", "plof", "discovery", z=zbad("systolic_bp", 1))   # remove SBP hit
            setrow(df, "SYNPASS1", trait, "plof", "discovery", z=zgood(trait, 9))
            setrow(df, "SYNPASS1", trait, "dmis", "discovery", z=zgood(trait, 2))
            setrow(df, "SYNPASS1", px, "plof", "replication", z=zgood(px, 4))
        sub = tmp_path / trait
        sub.mkdir()
        res, df = mutated_run(base_frames["pass"], edit, sub)
        assert result_tiers(res)[("SYNPASS1", trait)] == "A"
        assert res["verdict"]["pass_genes"] == [] and res["verdict"]["verdict"] == "LEAD"    # SYNLEAD1 only


# ---------------- C4: accepted behaviours ----------------
@pytest.mark.parametrize("p2,adverse", [(0.004, True), (0.008, True), (0.0109, True), (0.0113, False), (0.05, False)])
def test_c4_tradeoff_screen_is_one_sided(p2, adverse, cfg):
    assert oracle_adverse(+0.8, p2) == adverse
    rows = drop_rows(sbp_replicated_gene(), "GA", "coronary_disease", "discovery")
    rows.append(R("GA", "coronary_disease", "plof", "discovery", +0.8, p2))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == ("C" if adverse else "A")


def test_c4_gene_untested_in_replication_is_tier_B_not_D(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(R("OTHER", "hypertension", "plof", "replication", good("hypertension"), 0.01))   # trait present, gene absent
    hits = tiering.build_hits(T(rows), cfg)
    assert list(hits["tier"]) == ["B"] and list(hits["rep_status"]) == ["gene_untested"]


def test_c4_control_that_cannot_run_gives_kill_controls_not_evaluable(base_frames, tmp_path):
    def edit(df):
        df.drop(df.index[df["mask"] == "syn"], inplace=True)
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["verdict"]["verdict"] == "KILL"
    assert res["verdict"]["reason"].startswith("controls_not_evaluable")
    assert res["controls"]["valid"] is False


def test_c4_eur_only_is_reported_but_does_not_gate(base_frames, tmp_path):
    def edit(df):
        setrow(df, "SYNPASS1", "systolic_bp", "plof", "discovery_eur", beta=0.0, p=0.9)
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    row = [r for r in res["tiers"]["A"] if r["gene"] == "SYNPASS1"][0]
    assert row["eur_status"] == "failed"
    assert res["verdict"]["verdict"] == "PASS"


def test_c4_pcsk9_coronary_control_is_direction_only(null_table, cfg):
    df = set_row(null_table, "PCSK9", "coronary_disease", "plof", "discovery", beta=good("coronary_disease", 0.01), p=0.99)
    assert run_ctrl(df, cfg)["valid"] is True
    df = set_row(null_table, "PCSK9", "coronary_disease", "plof", "discovery", beta=bad("coronary_disease", 0.01), p=0.99)
    assert run_ctrl(df, cfg)["valid"] is False


def test_c4_unknown_replication_source_is_excluded_and_listed(base_frames, tmp_path):
    def edit(df):
        sel = df["cohort"] == "replication"
        df.loc[sel, "source"] = "mystery_biobank_burden"
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["verdict"]["verdict"] == "LEAD"
    assert "mystery_biobank_burden" in res["data"]["replication_sources_excluded"]


# ---------------- C9: the screening minimum also gates LEAD ----------------
@pytest.mark.parametrize("n,counts", [(0, False), (4, False), (5, True), (9, True)])
def test_c9_tier_B_gene_counts_toward_lead_only_if_five_outcomes_screened(n, counts, cfg):
    from protscan.run import decide
    rows = hit_rows("GA", "fluid_intelligence")
    for t in L_TRADEOFF[:n]:
        rows.append(R("GA", t, "plof", "discovery", 0.01, 0.5))
    hits = tiering.build_hits(T(rows), cfg)
    assert list(hits["tier"]) == ["B"]
    assert bool(hits["qualifies"].iloc[0]) is counts
    v = decide(CTRL_OK, hits)
    assert v["verdict"] == ("LEAD" if counts else "KILL")
    assert v["pass_genes"] == []


def test_c9_end_to_end_only_unscreened_lead_gene_gives_kill(base_frames, tmp_path):
    def edit(df):
        keep = L_TRADEOFF[:4]
        sel = (df["gene"] == "SYNLEAD1") & df["trait"].isin(L_TRADEOFF) & (~df["trait"].isin(keep))
        df.drop(df.index[sel], inplace=True)
    res, df = mutated_run(base_frames["lead"], edit, tmp_path)
    assert result_tiers(res)[("SYNLEAD1", "fluid_intelligence")] == "B"
    assert res["controls"]["valid"] is True
    assert res["verdict"]["verdict"] == "KILL" and res["verdict"]["lead_genes"] == []
    _, _, o_verdict = oracle(df)
    assert o_verdict == "KILL"


@pytest.mark.xfail(strict=False, reason="R-2: no control covers the replication-side sign convention (see review-1)")
def test_replication_sign_flip_is_detected_by_a_control_or_flag(synth, tmp_path):
    tables = synth.make_synthetic("pass", n_genes=N_GENES)
    rep = tables["burden_synthetic_replication"].copy()
    rep["beta"] = -rep["beta"]                       # simulate an adapter that reports the wrong allele's effect
    tables["burden_synthetic_replication"] = rep
    d = tmp_path / "data"
    write_tables(tables, d)
    from protscan.run import run_pipeline
    res = run_pipeline(CONFIG_PATH, d, tmp_path / "o" / "r.json")
    flagged = (res["controls"]["valid"] is False) or ("replication_sanity" in res["controls"])
    assert flagged, "PASS silently became %s with controls valid" % res["verdict"]["verdict"]
