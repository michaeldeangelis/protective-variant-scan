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
L_SYN_HITS_MAX = 0                            # Amendment 1 zero-hit gate; SUPERSEDED by Amendment 2 (A2b), kept for the record

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
# C10 (lead decisions on review-1)
L_SYN_MIN_COVERAGE = 0.90                             # C10a: syn rows cover >= 90 percent of plof (gene, trait) pairs
L_SYN_MIN_ROWS = 10000                                # C10a: and number >= 10,000
L_REP_SIGN = (("PCSK9", "hypercholesterolemia", -1), ("LDLR", "hypercholesterolemia", +1))   # C10b: established biology
L_A1FREQ_MAX = 0.5                                    # C10b adapter guard
# Reviewer's own pin of config/prereg.yaml (second, independent copy of schema.PINNED_CONFIG_SHA256, C10c).
# Update only together with a dated ledger entry.
# Pin history: 7507f535... (C10) replaced by the C11 config, which adds only syn_max_p1_fraction and replication_sign_max_p
# (reviewer diff of config/prereg.yaml against the C10 version, and an independent shasum, both checked).
L_CONFIG_SHA256 = "be69827539daac4c78b4c8f6ac7d931723abb41ab0e739c0505c357689310a9f"   # Amendment 2 config (was 7ce18c38... under C11)
# Amendment 2 (post-hoc after run 1; experiments.md), typed from the ledger text
L_SYN_CONTAM_MAX_FRACTION = 0.001                     # A2b: contaminated genes at most 0.1 percent of the genes tested
L_A2_LABEL = "Amendment 2 (post-hoc after run 1)"
# C11 (lead decisions on review-2), typed from experiments.md
L_SYN_MAX_P1_FRACTION = 0.05                          # C11b: more than 5 percent of syn rows at p = 1 makes the control not evaluable
L_REP_SIGN_MAX_P = 0.05                               # C11c: a sign-control arm is evaluable only if its p is below 0.05


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
    assert "syn_hits_max" not in c, "A2b replaced the zero-hit gate"
    assert c["syn_contaminated_max_fraction"] == L_SYN_CONTAM_MAX_FRACTION   # A2b
    assert c["syn_min_coverage"] == L_SYN_MIN_COVERAGE                    # C10a
    assert c["syn_min_rows"] == L_SYN_MIN_ROWS
    assert c["syn_max_p1_fraction"] == L_SYN_MAX_P1_FRACTION               # C11b
    assert c["replication_sign_max_p"] == L_REP_SIGN_MAX_P                 # C11c
    got = tuple((x["gene"], x["trait"], x["expected_beta_sign"]) for x in c["replication_sign"])
    assert got == L_REP_SIGN, got                                          # C10b
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


def test_config_file_bytes_match_reviewer_pin_and_code_pin():
    import hashlib
    from protscan import schema
    assert hashlib.sha256(CONFIG_PATH.read_bytes()).hexdigest() == L_CONFIG_SHA256
    assert schema.PINNED_CONFIG_SHA256 == L_CONFIG_SHA256


def test_pyproject_declares_runtime_dependencies():
    text = (ROOT / "pyproject.toml").read_text()
    dep_line = [ln for ln in text.splitlines() if ln.startswith("dependencies")][0]
    for pkg in ("pandas", "numpy", "scipy", "pyyaml", "requests", "pyarrow"):
        assert '"%s"' % pkg in dep_line, pkg


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
    assert c.syn_contaminated_max_fraction == L_SYN_CONTAM_MAX_FRACTION
    assert not hasattr(c, "syn_hits_max")
    assert c.syn_min_coverage == L_SYN_MIN_COVERAGE and c.syn_min_rows == L_SYN_MIN_ROWS
    assert c.syn_max_p1_fraction == L_SYN_MAX_P1_FRACTION and c.replication_sign_max_p == L_REP_SIGN_MAX_P
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


def baseline_null(n_genes=500, seed=7):
    """Null discovery table (plof/dmis/syn, N(0,1) z-scores) plus the lipid positive controls and the
    replication-sign control rows. 500 genes x 23 traits = 11,500 synonymous rows (C10a floor is 10,000)."""
    rng = np.random.default_rng(seed)
    genes = ["NULLG%04d" % i for i in range(n_genes)]
    traits = sorted(L_PANEL) + list(L_TRADEOFF) + list(CTRL_TRAITS)
    g_col, t_col, m_col = [], [], []
    for g in genes:
        for t in traits:
            for m in ("plof", "dmis", "syn"):
                g_col.append(g)
                t_col.append(t)
                m_col.append(m)
    z = rng.standard_normal(len(g_col))
    df = pd.DataFrame(dict(gene=g_col, trait=t_col, mask=m_col, cohort="discovery", beta=z * 0.1, se=0.1,
                           p=2 * norm.sf(np.abs(z)), n_carriers=100, n_total=400000, source=DISC_SRC))
    extra = [
        R("PCSK9", "ldl", "plof", "discovery", good("ldl"), 1e-30),
        R("PCSK9", "coronary_disease", "plof", "discovery", good("coronary_disease"), 1e-6),
        R("ANGPTL4", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-20),
        R("APOC3", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-20),
        R("PCSK9", "hypercholesterolemia", "plof", "replication", -0.8, 1e-6),      # C10b: PCSK9 lowers risk
        R("LDLR", "hypercholesterolemia", "plof", "replication", +0.9, 1e-8),       # C10b: LDLR raises risk
    ]
    return pd.concat([df[COLUMNS], pd.DataFrame(extra, columns=COLUMNS)], ignore_index=True)


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
    assert c["negative_synonymous"]["n_syn_hit_genes"] == 0 and c["negative_synonymous"]["n_contaminated"] == 0
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
    """Intent kept: one synonymous hit still fails the control when fewer than 1,000 genes are tested (A2b: 1 of 500 is 0.2 percent)."""
    df = set_row(null_table, "NULLG0001", "hand_grip_strength", "syn", "discovery",
                 beta=good("hand_grip_strength"), p=1e-9)
    c = run_ctrl(df, cfg)
    n = c["negative_synonymous"]
    assert n["n_syn_hit_genes"] == 1 and n["n_contaminated"] == 1 and n["n_genes_tested"] == 500
    assert n["lambda_ok"] is True and n["contamination_ok"] is False              # fails on the contamination bound, not on lambda
    assert c["valid"] is False and "negative_synonymous" in c["failed"]


def test_synonymous_hit_in_harmful_direction_is_now_counted_as_contamination(null_table, cfg):
    """Amendment 2 supersedes Amendment 1 here: the direction no longer matters (A2a, either direction)."""
    df = set_row(null_table, "NULLG0001", "hand_grip_strength", "syn", "discovery",
                 beta=bad("hand_grip_strength"), p=1e-9)
    n = run_ctrl(df, cfg)["negative_synonymous"]
    assert n["n_syn_hit_genes"] == 0                       # the old beneficial-direction panel-trait count is informational only
    assert n["n_contaminated"] == 1 and n["contaminated_genes"] == ["NULLG0001"] and n["contamination_ok"] is False


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
N_GENES = 500      # synthetic universe must exceed the C10a floor of 10,000 synonymous rows


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
    syn_all = df.query("cohort == 'discovery' and mask == 'syn'")
    # A2a: any synonymous hit at the discovery threshold, either direction, any trait, contaminates the gene
    contaminated = set(r.gene for r in syn_all.itertuples() if r.p < L_DISCOVERY_P_ACCEPTED)
    hits_all = []
    for r in d.itertuples():
        if r.trait in L_PANEL and oracle_discovery(r.beta, r.p, L_PANEL[r.trait][1]):
            hits_all.append((r.gene, r.trait))
    hits = [(g, t) for g, t in hits_all if g not in contaminated]
    rep = df.query("cohort == 'replication' and mask == 'plof'")
    rep = rep[[is_independent_source(s) for s in rep["source"]]]
    # C11a: a replication row is evidence only if se is finite and above 0 and beta is nonzero
    rep = rep[[(np.isfinite(se) and se > 0 and b != 0) for se, b in zip(rep["se"], rep["beta"])]]
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
    syn_inf = syn[syn["p"] < 1]                                                                    # C11b
    lam = (float(np.median(chi2.isf(np.clip(syn_inf["p"].to_numpy(), 1e-300, 1.0), 1)) / CHI2_MEDIAN)
           if len(syn_inf) else float("inf"))
    p1_fraction = (1.0 - len(syn_inf) / len(syn)) if len(syn) else 1.0
    n_tested = len(set(syn["gene"]))

    def rowof(g, t):
        x = d[(d["gene"] == g) & (d["trait"] == t)]
        return x.iloc[0] if len(x) else None

    ok = lam < L_LAMBDA_GC_MAX and p1_fraction <= L_SYN_MAX_P1_FRACTION + 1e-12
    ok = ok and n_tested > 0 and len(contaminated) * 1000 <= n_tested                          # A2b: at most 0.1 percent, exact integer test
    plof_pairs = set(zip(d["gene"], d["trait"]))
    syn_pairs = set(zip(syn["gene"], syn["trait"]))
    coverage = (len(plof_pairs.intersection(syn_pairs)) / len(plof_pairs)) if plof_pairs else 0.0
    ok = ok and len(syn) >= L_SYN_MIN_ROWS and coverage >= L_SYN_MIN_COVERAGE                    # C10a
    signs = []
    for g, tr, want in L_REP_SIGN:                                                               # C10b
        hit_rows_ = rep[(rep["gene"] == g) & (rep["trait"] == tr)]
        if len(hit_rows_) and hit_rows_.iloc[0].p < L_REP_SIGN_MAX_P:                            # C11c: powered arms only
            signs.append(hit_rows_.iloc[0].beta * want > 0)
    ok = ok and len(signs) >= 1 and all(signs)
    x = rowof("PCSK9", "ldl")
    ok = ok and x is not None and oracle_discovery(x.beta, x.p, -1)
    ok = ok and ("PCSK9", "ldl") in hits_all                                                   # C11d + A2: discovery path, contamination irrelevant
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
    for name in ("pass", "lead", "nolead", "broken_positive", "broken_lambda", "broken_syn_hit", "broken_rep_sign",
                 "unscreened", "contaminated_tierA", "contaminated_systemic", "contaminated_minor"):
        d = tmp_path_factory.mktemp("scn_" + name)
        tables = synth.make_synthetic(name, n_genes=N_GENES)
        write_tables(tables, d)
        res = run_cli(d, d / "out" / "protective-scan.json")
        runs[name] = (tables, res, d)
    return runs


@pytest.mark.parametrize("name,verdict", [("pass", "PASS"), ("lead", "LEAD"), ("nolead", "KILL"),
                                          ("broken_positive", "KILL"), ("broken_lambda", "KILL"),
                                          ("broken_syn_hit", "KILL"), ("broken_rep_sign", "KILL"),
                                          ("unscreened", "KILL"), ("contaminated_tierA", "LEAD"),
                                          ("contaminated_systemic", "KILL"), ("contaminated_minor", "PASS")])
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


def test_trivial_rung_is_the_synonymous_mask(base_big, tmp_path):
    """Trivial rung (Amendment 1): beneficial-direction panel-trait synonymous hits. Kept as such after Amendment 2, which
    adds the wider contamination list next to it (either direction, any trait)."""
    def edit(df):
        setrow(df, "SYNSYN1", "hand_grip_strength", "syn", "discovery", z=zgood("hand_grip_strength", 7))       # counts
        setrow(df, "SYNG00000", "bmi", "syn", "discovery", z=zbad("bmi", 7))                                    # harmful: contamination only
    res, _ = mutated_run(base_big, edit, tmp_path)
    assert res["controls"]["negative_synonymous"]["syn_hit_genes"] == ["SYNSYN1"]
    assert res["rungs"]["trivial"]["n_genes"] == 1
    assert res["contamination"]["n_contaminated"] == 2 and res["controls"]["valid"] is True
    assert res["rungs"]["trivial"]["status"] == "RUN"


def test_trivial_rung_is_zero_on_a_clean_table(scenario_runs):
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
def base_big(synth):
    """Pass scenario with more than 2,000 genes, so that a couple of contaminated genes stay within the 0.1 percent bound."""
    return all_rows(synth.make_synthetic("pass", n_genes=2100)).copy()


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
    """Systemic contamination: three null genes with synonymous hits (either direction, different traits): 3 of ~509 genes."""
    setrow(df, "SYNG00000", "hand_grip_strength", "syn", "discovery", z=zgood("hand_grip_strength", 7))
    setrow(df, "SYNG00001", "bmi", "syn", "discovery", z=zbad("bmi", 7))
    setrow(df, "SYNG00002", "coronary_disease", "syn", "discovery", z=zgood("coronary_disease", 7))


def edit_flip_rep_sign_controls(df):
    setrow(df, "PCSK9", "hypercholesterolemia", "plof", "replication", z=+5)
    setrow(df, "LDLR", "hypercholesterolemia", "plof", "replication", z=-5)


def edit_drop_rep_sign_rows(df):
    sel = (df["cohort"] == "replication") & df["gene"].isin(["PCSK9", "LDLR"]) & (df["trait"] == "hypercholesterolemia")
    df.drop(df.index[sel], inplace=True)


def edit_thin_synonymous(df):
    """Keep only 5 percent of the synonymous rows: coverage and row floor both violated."""
    syn = df.index[df["mask"] == "syn"]
    df.drop(syn[int(len(syn) * 0.05):], inplace=True)


@pytest.mark.parametrize("edit", [edit_flip_pcsk9_ldl, edit_pcsk9_ldl_weak, edit_drop_tg_controls, edit_flip_both_tg,
                                  edit_pcsk9_cad_flip, edit_lambda_high, edit_syn_hit, edit_flip_rep_sign_controls,
                                  edit_drop_rep_sign_rows, edit_thin_synonymous],
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


def test_synonymous_hit_harmful_direction_does_not_kill_but_excludes_the_gene(base_big, tmp_path):
    """Intent kept: one synonymous hit within the bound does not kill. A2a: the gene is excluded, whichever direction."""
    def edit(df):
        setrow(df, "SYNPASS1", "hand_grip_strength", "syn", "discovery", z=zbad("hand_grip_strength", 7))
    res, df = mutated_run(base_big, edit, tmp_path)
    assert res["controls"]["valid"] is True
    assert ("SYNPASS1", "systolic_bp") not in result_tiers(res)               # would have been Tier A
    assert res["verdict"]["verdict"] == "LEAD" and res["verdict"]["pass_genes"] == []
    assert res["verdict"]["contaminated_excluded"] == ["SYNPASS1"]
    ok, tiers, o_verdict = oracle(df)
    assert ok and o_verdict == "LEAD" and ("SYNPASS1", "systolic_bp") not in tiers


def _rename_gene(old, new):
    def edit(df):
        df.loc[df["gene"] == old, "gene"] = new
    return edit


def _sbp_rep(p):
    def edit(df):
        setrow(df, "SYNPASS1", "hypertension", "plof", "replication", beta=zgood("hypertension", 1) * 0.1, p=p)
    return edit


def edit_ukb_relabel(df):
    """Only the candidate gene's replication rows come from a UKB-derived source (controls stay evaluable)."""
    sel = (df["cohort"] == "replication") & (df["gene"] == "SYNPASS1")
    df.loc[sel, "source"] = "synthetic_genebass_hypertension_lookup"   # synthetic prefix: no real/synthetic mix


def edit_no_replication(df):
    """Only the candidate gene is missing from the replication cohort (controls stay evaluable)."""
    df.drop(df.index[(df["cohort"] == "replication") & (df["gene"] == "SYNPASS1")], inplace=True)


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
    ("lipid_pathway_gene_LPL", _rename_gene("SYNPASS1", "LPL")),
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
        extra["source"] = "synthetic_replication"
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
                rrow = R(g, px, "plof", "replication", s2 * benefit_sign(px) * 0.2, p2, source=src)
                u = rng.random()
                if u < 0.12:
                    rrow[5] = float("nan")                          # C11a: undefined se
                elif u < 0.20:
                    rrow[5] = 0.0
                elif u < 0.28:
                    rrow[4] = 0.0                                   # C11a: beta exactly 0
                rows.append(rrow)
        for out in L_TRADEOFF:
            for cohort in ("discovery", "replication"):
                if rng.random() < 0.6:
                    p3 = float(10 ** rng.uniform(-6, 0))
                    rows.append(R(g, out, "plof", cohort, float(rng.normal()), p3))
        if rng.random() < 0.15:                                       # A2a: synonymous hit, any trait, either direction
            t_syn = str(rng.choice(panel + list(L_TRADEOFF)))
            p_syn = float(10 ** rng.uniform(-9, -6))
            while 1.89e-7 <= p_syn <= 1.93e-7:
                p_syn = float(10 ** rng.uniform(-9, -6))
            rows.append(R(g, t_syn, "syn", "discovery", float(rng.choice([-0.4, 0.4])), p_syn))
        if rng.random() < 0.3:                                        # ordinary null synonymous row (no hit), a trait the hit cannot use
            rows.append(R(g, "triglycerides", "syn", "discovery", 0.0, float(rng.uniform(0.01, 1.0))))
    return T(rows)


@pytest.mark.parametrize("seed", range(25))
def test_fuzz_tiers_match_independent_oracle(seed, cfg):
    df = _fuzz_table(seed)
    hits = tiering.build_hits(df, cfg)
    got = dict(((r.gene, r.trait), r.tier) for r in hits.itertuples())
    _, want, _ = oracle(df)
    assert got == want


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
        sel = (df["cohort"] == "replication") & (df["gene"] == "SYNPASS1")
        df.loc[sel, "source"] = "synthetic_mystery_biobank_burden"
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["controls"]["valid"] is True
    assert res["verdict"]["verdict"] == "LEAD"
    assert "synthetic_mystery_biobank_burden" in res["data"]["replication_sources_excluded"]


@pytest.mark.parametrize("label", ["synthetic_mystery_biobank_burden", "synthetic_genebass_lookup", "synthetic_azphewas_v1"])
def test_c10b_whole_replication_cohort_from_untrusted_source_gives_kill_not_evaluable(label, base_frames, tmp_path):
    """No independent replication rows at all: the replication-sign control cannot run, so KILL (C10b)."""
    def edit(df):
        df.loc[df["cohort"] == "replication", "source"] = label
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["verdict"]["verdict"] == "KILL"
    assert res["verdict"]["reason"].startswith("controls_not_evaluable")
    assert "replication_sign" in res["controls"]["not_run"]
    assert res["controls"]["replication_sign"]["status"] == "NOT RUN"
    assert result_tiers(res).get(("SYNPASS1", "systolic_bp")) == "B"        # never Tier A without independent rows
    assert label in res["data"]["replication_sources_excluded"]


def test_c10b_no_replication_cohort_at_all_gives_kill(base_frames, tmp_path):
    def edit(df):
        df.drop(df.index[df["cohort"] == "replication"], inplace=True)
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["verdict"]["verdict"] == "KILL" and "replication_sign" in res["controls"]["not_run"]


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


# ======================================================================================
# (e) C10: closure of review-1 findings R-1 (a), R-2 (b), R-3 (c), R-9 (d), R-7 (e), R-8 (f), R-4/R-6/R-10/R-11 (g)
# ======================================================================================
def neg_of(df, cfg):
    return controls.run_controls(validate(df), cfg)["negative_synonymous"]


# ---------------- C10a: synonymous universe size and coverage ----------------
def test_c10a_three_synonymous_rows_are_not_evaluable_and_give_kill(cfg):
    """R-1 probe: 3 synonymous rows used to validate the pipeline. Now NOT RUN -> KILL controls_not_evaluable."""
    from protscan.run import decide
    rows = [R("PCSK9", "ldl", "plof", "discovery", good("ldl"), 1e-30),
            R("PCSK9", "coronary_disease", "plof", "discovery", good("coronary_disease"), 1e-4),
            R("APOC3", "triglycerides", "plof", "discovery", good("triglycerides"), 1e-30),
            R("PCSK9", "hypercholesterolemia", "plof", "replication", -0.8, 1e-6),
            R("PCSK9", "ldl", "syn", "discovery", 0.001, 0.6),
            R("G1", "ldl", "syn", "discovery", 0.001, 0.5),
            R("G2", "bmi", "syn", "discovery", 0.001, 0.4)]
    df = T(rows)
    c = controls.run_controls(df, cfg)
    n = c["negative_synonymous"]
    assert n["status"] == "NOT RUN" and n["n_rows"] == 3 and "rows" in n["reason"]
    assert "negative_synonymous" in c["not_run"] and c["valid"] is False
    hits = tiering.build_hits(df, cfg)
    v = decide(c, hits)
    assert v["verdict"] == "KILL" and v["reason"].startswith("controls_not_evaluable")


@pytest.mark.parametrize("n_genes,evaluable", [(434, False), (435, True)])
def test_c10a_row_floor_is_10000_synonymous_rows(n_genes, evaluable, cfg):
    df = baseline_null(n_genes=n_genes)
    n = neg_of(df, cfg)
    n_syn = int((df["mask"] == "syn").sum())
    assert (n_syn >= L_SYN_MIN_ROWS) == evaluable, n_syn
    assert (n["status"] != "NOT RUN") == evaluable


def _drop_syn_rows(df, k):
    """Remove the synonymous rows of the first k (gene, trait) pairs of the null block (plof rows stay)."""
    syn_idx = df.index[(df["mask"] == "syn")][:k]
    return df.drop(syn_idx)


def test_c10a_coverage_boundary_is_90_percent_of_plof_pairs(null_table, cfg):
    disc = null_table[null_table["cohort"] == "discovery"]
    n_plof = int((disc["mask"] == "plof").sum())              # includes the 4 control pairs without synonymous rows
    n_syn = int((disc["mask"] == "syn").sum())
    k_ok = n_syn - int(np.ceil(L_SYN_MIN_COVERAGE * n_plof - 1e-9))   # drop so coverage sits exactly at the minimum
    assert k_ok > 0
    ok = neg_of(_drop_syn_rows(null_table, k_ok), cfg)
    assert ok["status"] == "OK", ok.get("reason")
    assert ok["coverage"] >= L_SYN_MIN_COVERAGE - 1e-12
    low = neg_of(_drop_syn_rows(null_table, k_ok + 1), cfg)
    assert low["status"] == "NOT RUN" and "coverage" in low["reason"]
    assert low["coverage"] < L_SYN_MIN_COVERAGE


def test_c10a_coverage_is_measured_against_plof_pairs_not_row_count(null_table, cfg):
    """Plenty of synonymous rows, but for the wrong (gene, trait) pairs: rows above 10,000 yet coverage far below 90 percent."""
    df = null_table.copy()
    plof_null = df.index[(df["mask"] == "plof") & df["gene"].str.startswith("NULLG")]
    half = plof_null[: len(plof_null) // 2]
    df.loc[half, "gene"] = "ORPHAN" + df.loc[half, "gene"]
    n = neg_of(df, cfg)
    assert (df["mask"] == "syn").sum() >= L_SYN_MIN_ROWS
    assert n["status"] == "NOT RUN" and n["coverage"] < L_SYN_MIN_COVERAGE


def test_c10a_no_plof_rows_means_not_evaluable(null_table, cfg):
    df = null_table[null_table["mask"] != "plof"]
    n = neg_of(df, cfg)
    assert n["status"] == "NOT RUN" and n["coverage"] == 0.0


@pytest.mark.xfail(strict=False, reason="R2-1: 90 percent coverage still allows the omitted 10 percent to be exactly the hit genes")
def test_c10a_residual_selective_omission_of_hit_genes_is_detected(null_table, cfg):
    df = null_table.copy()
    hit = df.index[(df["gene"] == "NULLG0001") & (df["mask"] == "plof") & (df["trait"] == "ldl")]
    df.loc[hit, "beta"] = good("ldl")
    df.loc[hit, "p"] = 1e-12                                 # a discovery hit whose synonymous row is then withheld
    df = df.drop(df.index[(df["gene"] == "NULLG0001") & (df["mask"] == "syn")])
    c = controls.run_controls(validate(df), cfg)
    assert c["negative_synonymous"]["status"] != "OK"


# ---------------- C10b: replication-sign control ----------------
def rs_of(df, cfg):
    return controls.run_controls(validate(df), cfg)


def drop_rows_df(df, gene, trait, cohort):
    return df[~((df["gene"] == gene) & (df["trait"] == trait) & (df["cohort"] == cohort))].copy()


def test_c10b_both_evaluable_and_correct_is_ok(null_table, cfg):
    c = rs_of(null_table, cfg)
    assert c["replication_sign"]["status"] == "OK" and c["valid"] is True
    assert [k["status"] for k in c["replication_sign"]["checks"]] == ["OK", "OK"]


@pytest.mark.parametrize("pcsk9_beta,ldlr_beta", [(+0.8, +0.9), (-0.8, -0.9), (+0.8, -0.9)])
def test_c10b_any_evaluable_check_with_wrong_sign_fails_the_run(pcsk9_beta, ldlr_beta, null_table, cfg):
    df = set_row(null_table, "PCSK9", "hypercholesterolemia", "plof", "replication", beta=pcsk9_beta)
    df = set_row(df, "LDLR", "hypercholesterolemia", "plof", "replication", beta=ldlr_beta)
    c = rs_of(df, cfg)
    assert c["replication_sign"]["status"] == "FAIL" and c["valid"] is False
    assert "replication_sign" in c["failed"]


@pytest.mark.parametrize("gene,other", [("PCSK9", "LDLR"), ("LDLR", "PCSK9")])
def test_c10b_one_evaluable_check_is_enough_but_must_hold(gene, other, null_table, cfg):
    df = drop_rows_df(null_table, other, "hypercholesterolemia", "replication")
    assert rs_of(df, cfg)["valid"] is True
    want = dict((g, sgn) for g, _, sgn in L_REP_SIGN)[gene]
    bad_df = set_row(df, gene, "hypercholesterolemia", "plof", "replication", beta=-want * 0.5)
    c = rs_of(bad_df, cfg)
    assert c["valid"] is False and "replication_sign" in c["failed"]


def test_c10b_neither_evaluable_is_not_run_and_kills(null_table, cfg):
    df = drop_rows_df(drop_rows_df(null_table, "PCSK9", "hypercholesterolemia", "replication"),
                      "LDLR", "hypercholesterolemia", "replication")
    c = rs_of(df, cfg)
    assert c["replication_sign"]["status"] == "NOT RUN" and "replication_sign" in c["not_run"] and c["valid"] is False


@pytest.mark.parametrize("src", ["genebass", "azphewas_v1", "mystery_biobank"])
def test_c10b_control_rows_from_untrusted_sources_do_not_count(src, null_table, cfg):
    df = null_table.copy()
    df.loc[df["cohort"] == "replication", "source"] = src
    assert rs_of(df, cfg)["replication_sign"]["status"] == "NOT RUN"


def test_c10b_sign_is_direction_only_so_noise_can_pass_documented(null_table, cfg):
    """Scope of the control (review-2 R2-3): the correct sign at p = 0.99 passes; only the sign convention is tested."""
    df = set_row(null_table, "PCSK9", "hypercholesterolemia", "plof", "replication", beta=-0.001, p=0.99)
    assert rs_of(df, cfg)["replication_sign"]["status"] == "OK"


def test_c10b_end_to_end_flipped_replication_gives_kill_control_failed(synth, tmp_path):
    """R-2 probe: every replication beta sign-flipped. Was PASS -> LEAD with controls valid. Now KILL, control failed."""
    tables = synth.make_synthetic("pass", n_genes=N_GENES)
    rep = tables["burden_synthetic_replication"].copy()
    rep["beta"] = -rep["beta"]
    tables["burden_synthetic_replication"] = rep
    d = tmp_path / "data"
    write_tables(tables, d)
    from protscan.run import run_pipeline
    res = run_pipeline(CONFIG_PATH, d, tmp_path / "o" / "r.json")
    assert res["controls"]["valid"] is False
    assert res["controls"]["replication_sign"]["status"] == "FAIL"
    assert res["controls"]["failed"] == ["replication_sign"]
    assert [k["status"] for k in res["controls"]["replication_sign"]["checks"]] == ["FAIL", "FAIL"]
    assert res["verdict"]["verdict"] == "KILL" and res["verdict"]["reason"].startswith("control_failed")


# ---------------- C10b adapter guard: A1FREQ <= 0.5 ----------------
def _finngen_raw(rows):
    cols = ["PHENO", "ID", "A1FREQ", "N", "TEST", "BETA", "SE", "LOG10P"]
    return pd.DataFrame(rows, columns=cols)


def test_c10b_a1freq_guard_keeps_le_half_drops_gt_half_and_nan():
    from protscan.adapters import finngen
    assert finngen.A1FREQ_MAX == L_A1FREQ_MAX
    raw = _finngen_raw([["E4_HYPERCHOL", "G1", 0.001, 1000, "ADD", -0.5, 0.1, 6.0],
                        ["E4_HYPERCHOL", "G2", 0.5, 1000, "ADD", -0.5, 0.1, 6.0],
                        ["E4_HYPERCHOL", "G3", 0.5001, 1000, "ADD", -0.5, 0.1, 6.0],
                        ["E4_HYPERCHOL", "G4", 0.999, 1000, "ADD", -0.5, 0.1, 6.0],
                        ["E4_HYPERCHOL", "G5", float("nan"), 1000, "ADD", -0.5, 0.1, 6.0]])
    out = finngen.normalize_endpoint(raw, "E4_HYPERCHOL")
    assert sorted(out["gene"]) == ["G1", "G2"]
    assert (out["beta"] == -0.5).all()                       # dropped, never flipped
    summ = finngen.conversion_summary(raw)
    assert int(summ["rows_read"].sum()) == 5 and int(summ["dropped_a1freq"].sum()) == 3


def test_c10b_a1freq_guard_applies_to_every_endpoint_built():
    from protscan.adapters import finngen
    rows = []
    for ep in ("E4_HYPERCHOL", "I9_HYPTENS", "E4_OBESITY", "T2D"):
        rows.append([ep, "GOOD", 0.01, 1000, "ADD", -0.3, 0.1, 4.0])
        rows.append([ep, "FLIPPED", 0.9, 1000, "ADD", 0.3, 0.1, 4.0])
    raw = _finngen_raw(rows)
    for trait in ("hypercholesterolemia", "hypertension", "obesity", "type_2_diabetes"):
        df = finngen.build_trait(trait, raw)
        assert list(df["gene"]) == ["GOOD"], trait


# ---------------- C10c: pinned config hash ----------------
NL = chr(10)


def _run_with_config(synth, tmp_path, text, scenario="pass"):
    cfg_copy = tmp_path / "cfg.yaml"
    cfg_copy.write_bytes(text.encode())
    d = tmp_path / "data"
    write_tables(synth.make_synthetic(scenario, n_genes=N_GENES), d)
    from protscan.run import run_pipeline
    res = run_pipeline(cfg_copy, d, tmp_path / "o" / "protective-scan.json")
    return res, tmp_path / "o"


def test_c10c_edited_config_is_stamped_non_preregistered_everywhere(synth, tmp_path):
    """R-3 probe: a loosened control (here syn_max_p1_fraction 0.05 to 0.5) must be stamped, whatever the verdict."""
    text = CONFIG_PATH.read_text().replace("syn_max_p1_fraction: 0.05", "syn_max_p1_fraction: 0.5")
    assert text != CONFIG_PATH.read_text()
    res, out = _run_with_config(synth, tmp_path, text, scenario="pass")
    assert res["config_matches_ledger"] is False and res["verdict"]["config_matches_ledger"] is False
    assert "NON-PREREGISTERED" in res["verdict"]["label"]
    assert res["config_sha256"] != L_CONFIG_SHA256 and res["config_sha256_pinned"] == L_CONFIG_SHA256
    saved = json.loads((out / "protective-scan.json").read_text())
    assert saved["config_matches_ledger"] is False and "NON-PREREGISTERED" in saved["verdict"]["label"]
    assert "NON-PREREGISTERED" in (out / "report.md").read_text()


@pytest.mark.parametrize("variant", ["trailing_newline", "comment_only", "extra_key", "crlf"])
def test_c10c_cosmetic_edits_cannot_bypass_the_pin(variant, synth, tmp_path):
    text = CONFIG_PATH.read_text()
    if variant == "trailing_newline":
        text = text + NL
    elif variant == "comment_only":
        text = "# harmless comment" + NL + text
    elif variant == "extra_key":
        text = text + "extra_unused: 1" + NL
    else:
        text = text.replace(NL, chr(13) + NL)
    res, _ = _run_with_config(synth, tmp_path, text)
    assert res["config_matches_ledger"] is False, variant


def test_c10c_identical_bytes_at_another_path_match(synth, tmp_path):
    res, _ = _run_with_config(synth, tmp_path, CONFIG_PATH.read_text())
    assert res["config_matches_ledger"] is True and "NON-PREREGISTERED" not in res["verdict"]["label"]


def test_c10c_cli_exit_codes(synth, tmp_path):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    d = tmp_path / "data"
    write_tables(synth.make_synthetic("pass", n_genes=N_GENES), d)
    edited = tmp_path / "edited.yaml"
    edited.write_text(CONFIG_PATH.read_text().replace("syn_max_p1_fraction: 0.05", "syn_max_p1_fraction: 0.5"))
    bad = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(edited), "--data", str(d),
                          "--out", str(tmp_path / "bad" / "r.json")], capture_output=True, text=True, env=env, cwd=str(ROOT))
    assert bad.returncode == 2 and "NON-PREREGISTERED" in bad.stderr
    good_ = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(CONFIG_PATH), "--data", str(d),
                            "--out", str(tmp_path / "good" / "r.json")], capture_output=True, text=True, env=env, cwd=str(ROOT))
    assert good_.returncode == 0 and "NON-PREREGISTERED" not in good_.stdout + good_.stderr
    assert good_.stdout.startswith("verdict: PASS")


# ---------------- C10d: rows with undefined se are kept ----------------
def test_c10d_schema_keeps_rows_with_nan_se_and_still_drops_nan_beta_or_p():
    rows = [R("G1", "ldl", "syn", "discovery", 0.0, 1.0),
            R("G2", "ldl", "syn", "discovery", 0.1, 0.5),
            R("G3", "ldl", "syn", "discovery", float("nan"), 0.5),
            R("G4", "ldl", "syn", "discovery", 0.1, float("nan"))]
    rows[0][5] = float("nan")                                # se undefined (beta 0, p 1)
    df = validate(pd.DataFrame(rows, columns=COLUMNS))
    assert sorted(df["gene"]) == ["G1", "G2"] and df.attrs["dropped_nan_rows"] == 2
    assert df.loc[df["gene"] == "G1", "se"].isna().all()


def test_c10d_lambda_and_coverage_use_p_only_so_nan_se_rows_count(null_table, cfg):
    df = null_table.copy()
    sel = df.index[df["mask"] == "syn"][:2000]
    df.loc[sel, "se"] = float("nan")
    n_all = neg_of(null_table, cfg)
    n_nan = neg_of(df, cfg)
    assert n_nan["n_rows"] == n_all["n_rows"] and n_nan["coverage"] == n_all["coverage"]
    assert n_nan["lambda_gc"] == n_all["lambda_gc"] and n_nan["status"] == "OK"


def test_c10d_pipeline_serializes_results_with_nan_se_rows(base_frames, tmp_path):
    def edit(df):
        sel = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("SYNG")][:500]
        df.loc[sel, "se"] = float("nan")
        df.loc[sel, "beta"] = 0.0
        df.loc[sel, "p"] = 1.0
    res, _ = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["controls"]["valid"] is True and res["verdict"]["verdict"] == "PASS"
    json.loads((tmp_path / "out" / "protective-scan.json").read_text())


def test_c10d_finngen_row_with_degenerate_se_but_significant_p_is_not_used():
    from protscan.adapters import finngen
    raw = _finngen_raw([["I9_HYPTENS", "DEGEN", 0.001, 1000, "ADD", -2.0, 0.0, 6.0],
                        ["I9_HYPTENS", "NANSE", 0.001, 1000, "ADD", -2.0, float("nan"), 6.0],
                        ["I9_HYPTENS", "FINE", 0.001, 1000, "ADD", -0.5, 0.1, 6.0]])
    out = finngen.normalize_endpoint(raw, "I9_HYPTENS")
    assert list(out["gene"]) == ["FINE"]


def test_c10d_replication_row_with_undefined_se_and_tiny_p_cannot_give_tier_A(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    r = R("GA", "hypertension", "plof", "replication", good("hypertension"), 1e-9)
    r[5] = float("nan")
    rows.append(r)
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A"


# ---------------- C10e: informational harmful-direction panel column ----------------
def test_c10e_harmful_panel_column_lists_other_traits_and_does_not_gate(cfg):
    rows = sbp_replicated_gene("GA")
    rows.append(R("GA", "fev1", "plof", "discovery", bad("fev1"), 1e-12))            # significant, harmful
    rows.append(R("GA", "hand_grip_strength", "plof", "discovery", bad("hand_grip_strength"), 1e-3))   # not significant
    hits = tiering.build_hits(T(rows), cfg)
    row = hits[hits["trait"] == "systolic_bp"].iloc[0]
    assert [h["trait"] for h in row["harmful_panel"]] == ["fev1"]
    assert row["tier"] == "A" and bool(row["qualifies"]) is True                      # informational only


def test_c10e_end_to_end_report_shows_the_column(scenario_runs):
    _, res, d = scenario_runs["pass"]
    assert "harmful panel (info)" in (d / "out" / "report.md").read_text()
    assert all("harmful_panel" in r for t in "ABCD" for r in res["tiers"][t])


# ---------------- C10f: caveats printed with any PASS or LEAD ----------------
@pytest.mark.parametrize("name", ["pass", "lead"])
def test_c10f_pass_and_lead_reports_print_c5_and_c7(name, scenario_runs):
    _, res, d = scenario_runs[name]
    assert res["verdict"]["verdict"] in ("PASS", "LEAD")
    text = (d / "out" / "report.md").read_text()
    assert "C5" in text and "missense|LC" in text and "direction-only" in text
    assert "C7" in text and "fluid intelligence and reaction time" in text


def test_c10f_synthetic_flag_reaches_the_verdict_block_and_stdout(scenario_runs):
    _, res, _ = scenario_runs["pass"]
    assert res["verdict"]["synthetic"] is True and "[SYNTHETIC DATA]" in res["verdict"]["label"]


# ---------------- C10g ----------------
def test_c10g_untrusted_replication_rows_do_not_count_as_screening(cfg):
    disc4 = sbp_gene_with_screen("GA", L_TRADEOFF[:4])
    trusted = disc4 + [R("GA", L_TRADEOFF[4], "plof", "replication", 0.01, 0.5, source="finngen_r13")]
    untrusted = disc4 + [R("GA", L_TRADEOFF[4], "plof", "replication", 0.01, 0.5, source="mystery_biobank")]
    assert tier_map(trusted, cfg)[("GA", "systolic_bp")] == "A"
    assert tier_map(untrusted, cfg)[("GA", "systolic_bp")] == "B"


def test_c10g_untrusted_replication_rows_still_count_for_adverse_detection(cfg):
    rows = sbp_gene_with_screen("GA", L_TRADEOFF[:5])
    rows.append(R("GA", L_TRADEOFF[5], "plof", "replication", +0.9, 1e-5, source="mystery_biobank"))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "C"     # conservative: an adverse signal is never ignored by source


def test_c11d_positive_control_is_judged_on_the_discovery_effect_only(null_table, cfg):
    """C11d supersedes the C10g tier-C exclusion: the trade-off screen never changes control status."""
    c = controls.run_controls(validate(null_table), cfg)
    ldl = [p for p in c["positive"] if p["id"] == "pcsk9_ldl_lower"][0]
    assert ldl["via_pipeline_discovery_path"] is True and "via_pipeline_tier_path" not in ldl
    assert ldl["tested"][0]["pipeline_tier"] in ("A", "B", "D")
    # a harmful, significant type 2 diabetes signal on PCSK9 (known biology) makes PCSK9/LDL Tier C but leaves the control OK
    for extra in (R("PCSK9", "type_2_diabetes", "plof", "discovery", +0.9, 1e-5),
                  R("PCSK9", "type_2_diabetes", "plof", "replication", +0.9, 1e-5)):
        df = pd.concat([null_table, pd.DataFrame([extra], columns=COLUMNS)], ignore_index=True)
        c2 = controls.run_controls(validate(df), cfg)
        ldl2 = [p for p in c2["positive"] if p["id"] == "pcsk9_ldl_lower"][0]
        assert ldl2["tested"][0]["pipeline_tier"] == "C"
        assert ldl2["status"] == "OK" and c2["valid"] is True and c2["failed"] == [], extra[3]


def test_c11d_positive_control_still_fails_on_the_discovery_effect(null_table, cfg):
    weak = set_row(null_table, "PCSK9", "ldl", "plof", "discovery", p=1e-5)             # not at the discovery threshold
    wrong = set_row(null_table, "PCSK9", "ldl", "plof", "discovery", beta=bad("ldl"))   # wrong direction
    for df in (weak, wrong):
        c = controls.run_controls(validate(df), cfg)
        assert c["valid"] is False and "pcsk9_ldl_lower" in c["failed"]


def test_c11d_end_to_end_pcsk9_type_2_diabetes_harm_does_not_kill(base_frames, tmp_path):
    def edit(df):
        setrow(df, "PCSK9", "type_2_diabetes", "plof", "discovery", z=6.0)
    res, df = mutated_run(base_frames["pass"], edit, tmp_path)
    assert result_tiers(res)[("PCSK9", "ldl")] == "C"
    assert res["controls"]["valid"] is True and res["verdict"]["verdict"] == "PASS"
    ok, tiers, o_verdict = oracle(df)
    assert ok and o_verdict == "PASS" and tiers[("PCSK9", "ldl")] == "C"


def test_c10g_synthetic_and_real_sources_cannot_be_mixed(tmp_path):
    from protscan.schema import load_burden
    from protscan.run import run_pipeline
    synth_rows = [R("G1", "ldl", "syn", "discovery", 0.1, 0.5, source="synthetic_discovery")]
    real_rows = [R("G2", "ldl", "syn", "discovery", 0.1, 0.5, source="genebass")]
    pd.DataFrame(synth_rows, columns=COLUMNS).to_csv(tmp_path / "burden_a.csv", index=False)
    pd.DataFrame(real_rows, columns=COLUMNS).to_csv(tmp_path / "burden_b.csv", index=False)
    with pytest.raises(ValueError, match="mixes synthetic"):
        load_burden(tmp_path)
    with pytest.raises(ValueError, match="mixes synthetic"):
        run_pipeline(CONFIG_PATH, tmp_path, tmp_path / "o" / "r.json")


# ---------------- review-2 documented residual gaps (non-strict xfail) ----------------
def _deflation_masked_table(null_table, frac_p1=0.45, inflate=2.0):
    """Informative synonymous rows inflated by a factor, diluted with uninformative p = 1 rows (C10d keeps them)."""
    df = null_table.copy()
    syn_idx = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("NULLG")]
    z = np.clip(df.loc[syn_idx, "beta"] / df.loc[syn_idx, "se"] * inflate, -4.0, 4.0)   # no tail hits: isolate lambda
    df.loc[syn_idx, "p"] = 2 * norm.sf(np.abs(z))
    n1 = int(len(syn_idx) * frac_p1)
    df.loc[syn_idx[:n1], "p"] = 1.0
    df.loc[syn_idx[:n1], "beta"] = 0.0
    df.loc[syn_idx[:n1], "se"] = float("nan")
    return df, syn_idx[n1:]


def test_c10d_deflation_masking_setup_is_a_real_inflation(null_table):
    """Sanity for the xfail below: the informative synonymous rows alone are strongly inflated (lambda well above 1.10)."""
    df, informative = _deflation_masked_table(null_table)
    lam = float(np.median(chi2.isf(np.clip(df.loc[informative, "p"].to_numpy(), 1e-300, 1.0), 1)) / CHI2_MEDIAN)
    assert lam > 3.0


def test_c10d_lambda_control_is_not_masked_by_uninformative_rows(null_table, cfg):
    df, _ = _deflation_masked_table(null_table)
    n = neg_of(df, cfg)
    assert n["status"] != "OK"


# ======================================================================================
# (f) C11a-d: closure of review-2 findings R2-2 (a), R2-4 (b), R2-3 (c), R2-5 (d)
# ======================================================================================
def rep_row(gene, trait, beta, p, se=0.1, source="finngen_r13"):
    r = R(gene, trait, "plof", "replication", beta, p, source=source)
    r[5] = se
    return r


# ---------------- C11a: only informative replication rows count ----------------
NAN = float("nan")


@pytest.mark.parametrize("label,beta,se", [("se_nan", -0.3, NAN), ("se_zero", -0.3, 0.0), ("se_negative_like_inf", -0.3, float("inf")),
                                           ("beta_zero", 0.0, 0.1)])
def test_c11a_uninformative_replication_row_cannot_replicate(label, beta, se, cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(rep_row("GA", "hypertension", beta, 1e-9, se=se))
    if se == 0.0:
        rows[-1][5] = 0.0
    hits = tiering.build_hits(T(rows) if se != float("inf") else validate(pd.DataFrame(rows, columns=COLUMNS)), cfg)
    assert list(hits["tier"]) == ["B"], label                                 # never A; no informative row: trait_missing
    assert list(hits["rep_status"]) == ["trait_missing"], label


def test_c11a_smallest_finite_positive_se_and_nonzero_beta_still_count(cfg):
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(rep_row("GA", "hypertension", -1e-9, 1e-9, se=1e-12))
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "A"


def test_c11a_an_uninformative_row_does_not_hide_the_traits_other_genes(cfg):
    """A degenerate row for gene GB must not turn gene GA (no replication row at all) into trait_missing or vice versa."""
    rows = hit_rows("GA", "systolic_bp") + null_tradeoffs("GA", cohorts=("discovery",))
    rows.append(rep_row("GB", "hypertension", -0.3, 1e-9, se=NAN))
    rows.append(rep_row("GC", "hypertension", -0.3, 0.01, se=0.1))
    hits = tiering.build_hits(T(rows), cfg)
    assert list(hits["rep_status"]) == ["gene_untested"]                     # GC supplies the trait; GB row is ignored


def test_c11a_uninformative_rows_still_serve_lambda_and_coverage(null_table, cfg):
    df = null_table.copy()
    sel = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("NULLG")][:300]
    df.loc[sel, ["beta", "p"]] = [0.0, 1.0]
    df.loc[sel, "se"] = NAN
    n = neg_of(df, cfg)
    assert n["status"] == "OK" and n["n_p1_rows"] == 300 and n["coverage"] == neg_of(null_table, cfg)["coverage"]


def test_c11a_finngen_adapter_drops_undefined_se_with_p_below_1_and_counts_it():
    from protscan.adapters import finngen
    raw = _finngen_raw([["I9_HYPTENS", "DEGEN0", 0.001, 1000, "ADD", -2.0, 0.0, 6.0],
                        ["I9_HYPTENS", "DEGENNAN", 0.001, 1000, "ADD", -2.0, NAN, 6.0],
                        ["I9_HYPTENS", "NULLROW", 0.001, 1000, "ADD", 0.0, 0.0, 0.0],       # beta 0, p = 1: kept, se undefined
                        ["I9_HYPTENS", "FINE", 0.001, 1000, "ADD", -0.5, 0.1, 6.0],
                        ["I9_HYPTENS", "FLIPPED", 0.9, 1000, "ADD", 0.5, 0.1, 6.0],
                        ["I9_HYPTENS", "NOBETA", 0.001, 1000, "ADD", NAN, 0.1, 6.0]])
    out = finngen.normalize_endpoint(raw, "I9_HYPTENS").set_index("gene")
    assert sorted(out.index) == ["FINE", "NULLROW"]
    assert np.isnan(out.loc["NULLROW", "se"]) and out.loc["NULLROW", "p"] == 1.0
    summ = finngen.conversion_summary(raw).set_index("endpoint").loc["I9_HYPTENS"]
    assert summ["rows_read"] == 6 and summ["dropped_a1freq"] == 1 and summ["dropped_missing_beta_or_p"] == 1
    assert summ["dropped_undefined_se_p_lt_1"] == 2 and summ["kept_se_undefined"] == 1 and summ["rows_out"] == 2


def test_c11a_genebass_adapter_drops_beta_zero_with_p_below_1_and_keeps_p_1():
    from protscan.adapters import genebass
    recs = [dict(gene_symbol="GA", gene_id="E1", BETA_Burden=-0.5, Pvalue_Burden=1e-6),
            dict(gene_symbol="GB", gene_id="E2", BETA_Burden=0.0, Pvalue_Burden=0.5),      # defect: beta 0 with p < 1
            dict(gene_symbol="GC", gene_id="E3", BETA_Burden=0.0, Pvalue_Burden=1.0),      # uninformative: kept
            dict(gene_symbol="GD", gene_id="E4", BETA_Burden=0.2, Pvalue_Burden=1.0)]      # p = 1 with beta: kept
    qc = [dict(gene_id=e, gene_symbol=g, annotation="pLoF", CAF=0.001, keep_gene_burden=True, keep_gene_coverage=True,
               keep_gene_n_var=True) for e, g in (("E1", "GA"), ("E2", "GB"), ("E3", "GC"), ("E4", "GD"))]
    st = {}
    out = genebass.normalize_analysis(recs, qc, "syn", 1000, stats=st).set_index("gene")
    assert sorted(out.index) == ["GA", "GC", "GD"]
    assert st["rows_in"] == 4 and st["dropped_undefined_se_p_lt_1"] == 1 and st["kept_se_undefined"] == 2
    assert out.loc["GC", "p"] == 1.0 and np.isnan(out.loc["GC", "se"])


def test_c11a_ivw_combine_keeps_only_undefined_se_genes_with_p_1():
    from protscan.adapters import common
    a = pd.DataFrame(dict(gene=["G1", "G2", "G3"], mask="plof", beta=[0.0, 0.0, 0.3], se=[NAN, NAN, 0.1], p=[1.0, 0.4, 0.003],
                          n_carriers=10, n_total=100))
    b = pd.DataFrame(dict(gene=["G1", "G2", "G3"], mask="plof", beta=[0.0, 0.0, 0.5], se=[NAN, NAN, 0.2], p=[1.0, 0.4, 0.01],
                          n_carriers=10, n_total=100))
    out = common.ivw_combine([a, b], "max").set_index("gene")
    assert sorted(out.index) == ["G1", "G3"] and np.isnan(out.loc["G1", "se"]) and out.loc["G1", "p"] == 1.0


# ---------------- C11b: lambda_GC on p < 1 rows, p = 1 fraction gate ----------------
def _with_p1_rows(null_table, k):
    """Set the first k null synonymous rows to p = 1 (beta 0, se undefined)."""
    df = null_table.copy()
    sel = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("NULLG")][:k]
    df.loc[sel, "beta"] = 0.0
    df.loc[sel, "p"] = 1.0
    df.loc[sel, "se"] = NAN
    return df


def test_c11b_p1_fraction_boundary_is_five_percent_of_synonymous_rows(null_table, cfg):
    n_syn = int((null_table["mask"] == "syn").sum() - 0)
    k_ok = int(np.floor(L_SYN_MAX_P1_FRACTION * n_syn + 1e-9))               # exactly at the limit or just below
    ok = neg_of(_with_p1_rows(null_table, k_ok), cfg)
    assert ok["status"] == "OK" and ok["n_p1_rows"] == k_ok and ok["p1_fraction"] <= L_SYN_MAX_P1_FRACTION
    over = neg_of(_with_p1_rows(null_table, k_ok + 1), cfg)
    assert over["status"] == "NOT RUN" and over["p1_fraction"] > L_SYN_MAX_P1_FRACTION and "p = 1" in over["reason"]


def test_c11b_exactly_five_percent_is_evaluable(cfg):
    """A table whose p = 1 share is exactly 5.000 percent (500 of 10,000 rows) is evaluable; 501 is not."""
    def table_with(n_p1):
        rng = np.random.default_rng(11)
        p = rng.uniform(0.0, 0.999, size=10_000)
        p[:n_p1] = 1.0
        rows = []
        for i in range(10_000):
            rows.append(R("S%d" % i, "fev1", "plof", "discovery", 0.0, 0.5))
            rows.append(R("S%d" % i, "fev1", "syn", "discovery", 0.0, float(p[i])))
        return T(rows)
    assert controls.negative_control(table_with(500), cfg)["p1_fraction"] == 0.05
    assert controls.negative_control(table_with(500), cfg)["status"] != "NOT RUN"
    assert controls.negative_control(table_with(501), cfg)["status"] == "NOT RUN"


def test_c11b_lambda_is_computed_on_rows_with_p_below_1_and_both_are_reported(null_table, cfg):
    df = _with_p1_rows(null_table, 400)                                       # 3.5 percent p = 1 rows
    n = neg_of(df, cfg)
    inf_p = df.loc[(df["mask"] == "syn") & (df["p"] < 1), "p"].to_numpy()
    want = float(np.median(chi2.isf(np.clip(inf_p, 1e-300, 1.0), 1)) / CHI2_MEDIAN)
    all_p = df.loc[df["mask"] == "syn", "p"].to_numpy()
    want_all = float(np.median(chi2.isf(np.clip(all_p, 1e-300, 1.0), 1)) / CHI2_MEDIAN)
    assert n["lambda_gc"] == pytest.approx(want, rel=1e-12) and n["lambda_gc_all_rows"] == pytest.approx(want_all, rel=1e-12)
    assert n["lambda_gc_all_rows"] < n["lambda_gc"]                            # the p = 1 rows deflate the all-rows value only


def test_c11b_inflation_hidden_by_a_small_p1_share_is_still_caught(null_table, cfg):
    """4 percent p = 1 rows (evaluable) with inflated informative rows: the gate uses the p < 1 lambda and fails."""
    df = null_table.copy()
    syn_idx = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("NULLG")]
    z = np.clip(df.loc[syn_idx, "beta"] / df.loc[syn_idx, "se"] * 1.5, -4.0, 4.0)
    df.loc[syn_idx, "p"] = 2 * norm.sf(np.abs(z))
    df = _with_p1_rows(df, int(0.04 * len(syn_idx)))
    n = neg_of(df, cfg)
    assert n["status"] == "FAIL" and n["lambda_ok"] is False and n["lambda_gc"] > L_LAMBDA_GC_MAX


def test_c11b_heavy_p1_share_end_to_end_is_kill_not_evaluable(base_frames, tmp_path):
    def edit(df):
        sel = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("SYNG")]
        k = int(0.06 * (df["mask"] == "syn").sum())
        df.loc[sel[:k], ["beta", "p"]] = [0.0, 1.0]
        df.loc[sel[:k], "se"] = NAN
    res, df = mutated_run(base_frames["pass"], edit, tmp_path)
    assert res["controls"]["negative_synonymous"]["status"] == "NOT RUN"
    assert res["verdict"]["verdict"] == "KILL" and res["verdict"]["reason"].startswith("controls_not_evaluable")
    ok, _, o_verdict = oracle(df)
    assert not ok and o_verdict == "KILL"
    assert res["thresholds"]["syn_max_p1_fraction"] == L_SYN_MAX_P1_FRACTION


# ---------------- C11c: sign-control arms are evaluable only when powered ----------------
def _sign_table(null_table, pcsk9=(-0.8, 1e-6), ldlr=(0.9, 1e-8)):
    df = set_row(null_table, "PCSK9", "hypercholesterolemia", "plof", "replication", beta=pcsk9[0], p=pcsk9[1])
    return set_row(df, "LDLR", "hypercholesterolemia", "plof", "replication", beta=ldlr[0], p=ldlr[1])


@pytest.mark.parametrize("p,evaluable", [(0.0499, True), (0.05, False), (0.0501, False), (0.4, False)])
def test_c11c_arm_is_evaluable_only_below_p_005(p, evaluable, null_table, cfg):
    assert L_REP_SIGN_MAX_P == 0.05
    df = _sign_table(null_table, pcsk9=(+0.5, p), ldlr=(0.9, 1e-8))            # PCSK9 arm has the WRONG sign
    c = controls.run_controls(validate(df), cfg)
    arm = [k for k in c["replication_sign"]["checks"] if k["gene"] == "PCSK9"][0]
    if evaluable:
        assert arm["status"] == "FAIL" and c["replication_sign"]["reason"] == "sign_control_failed"
        assert "replication_sign" in c["failed"] and c["valid"] is False
    else:
        assert arm["status"] == "NOT RUN" and "underpowered" in arm["why"]
        assert c["replication_sign"]["status"] == "OK" and c["valid"] is True   # underpowered arms never count as failures


def test_c11c_a_powered_wrong_sign_arm_fails_even_when_the_other_arm_is_right(null_table, cfg):
    df = _sign_table(null_table, pcsk9=(-0.8, 1e-6), ldlr=(-0.9, 1e-8))
    c = controls.run_controls(validate(df), cfg)
    assert c["replication_sign"]["status"] == "FAIL" and c["replication_sign"]["reason"] == "sign_control_failed"


def test_c11c_no_evaluable_arm_is_not_evaluable_and_reason_is_recorded(null_table, cfg):
    from protscan.run import decide
    df = _sign_table(null_table, pcsk9=(-0.8, 0.3), ldlr=(0.9, 0.2))            # right signs but both underpowered
    c = controls.run_controls(validate(df), cfg)
    assert c["replication_sign"]["status"] == "NOT RUN" and c["replication_sign"]["reason"] == "sign_control_not_evaluable"
    assert "replication_sign" in c["not_run"] and c["valid"] is False
    v = decide(c, tiering.build_hits(validate(df), cfg))
    assert v["verdict"] == "KILL" and "sign_control_not_evaluable" in v["reason"] and v["reason"].startswith("controls_not_evaluable")


def test_c11c_failure_reason_reaches_the_verdict_text(null_table, cfg):
    from protscan.run import decide
    df = _sign_table(null_table, pcsk9=(+0.8, 1e-6), ldlr=(0.9, 1e-8))
    v = decide(controls.run_controls(validate(df), cfg), tiering.build_hits(validate(df), cfg))
    assert v["verdict"] == "KILL" and v["reason"].startswith("control_failed") and "sign_control_failed" in v["reason"]


@pytest.mark.parametrize("beta,se", [(0.9, NAN), (0.9, 0.0), (0.0, 0.1)])
def test_c11c_uninformative_rows_cannot_serve_as_an_arm(beta, se, null_table, cfg):
    df = _sign_table(null_table, pcsk9=(-0.8, 1e-6), ldlr=(0.9, 1e-8))
    i = df.index[(df["gene"] == "LDLR") & (df["cohort"] == "replication")][0]
    df.loc[i, ["beta", "se"]] = [beta, se]
    c = controls.run_controls(validate(df), cfg)
    arm = [k for k in c["replication_sign"]["checks"] if k["gene"] == "LDLR"][0]
    assert arm["status"] == "NOT RUN" and "no informative row" in arm["why"]
    assert c["valid"] is True                                                    # the PCSK9 arm still carries the control


def test_c11c_end_to_end_flipped_replication_is_still_control_failed(synth, tmp_path):
    tables = synth.make_synthetic("broken_rep_sign", n_genes=N_GENES)
    d = tmp_path / "data"
    write_tables(tables, d)
    from protscan.run import run_pipeline
    res = run_pipeline(CONFIG_PATH, d, tmp_path / "o" / "r.json")
    assert res["controls"]["replication_sign"]["reason"] == "sign_control_failed"
    assert res["verdict"]["verdict"] == "KILL" and "sign_control_failed" in res["verdict"]["reason"]
    assert "Reading rule" in (tmp_path / "o" / "report.md").read_text()
    assert res["thresholds"]["replication_sign_max_p"] == L_REP_SIGN_MAX_P


# ---------------- review-3 residual gaps (non-strict xfail) and documented conservative behaviour ----------------
@pytest.mark.xfail(strict=False, reason="R3-1: uninformative rows (beta 0, p 1, se undefined) count as trade-off outcomes screened")
def test_uninformative_tradeoff_rows_do_not_count_as_screened(cfg):
    rows = hit_rows("GA", "systolic_bp")
    rows.append(rep_row("GA", "hypertension", -0.3, 0.01))
    for t in L_TRADEOFF[:5]:
        r = R("GA", t, "plof", "discovery", 0.0, 1.0)
        r[5] = NAN
        rows.append(r)
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] != "A"


def _near_one_masked_table(null_table, p_near=0.99999, frac=0.45, inflate=1.5):
    df = null_table.copy()
    syn_idx = df.index[(df["mask"] == "syn") & df["gene"].str.startswith("NULLG")]
    z = np.clip(df.loc[syn_idx, "beta"] / df.loc[syn_idx, "se"] * inflate, -4.0, 4.0)
    df.loc[syn_idx, "p"] = 2 * norm.sf(np.abs(z))
    df.loc[syn_idx[: int(len(syn_idx) * frac)], "p"] = p_near
    return df


@pytest.mark.xfail(strict=False, reason="R3-2: the C11b filter is exact p = 1; uninformative rows stored as p just below 1 still deflate lambda_GC")
def test_lambda_control_is_not_masked_by_p_just_below_one(null_table, cfg):
    n = neg_of(_near_one_masked_table(null_table), cfg)
    assert n["status"] != "OK"


def test_c11b_exact_p_one_masking_is_closed_but_the_near_one_variant_shows_the_scope(null_table, cfg):
    """Scope check that pairs with the xfail above: with exact p = 1 the same construction is not evaluable."""
    df = _near_one_masked_table(null_table, p_near=1.0)
    assert neg_of(df, cfg)["status"] == "NOT RUN"


def test_c11a_adverse_screen_is_deliberately_not_filtered_by_informativeness(cfg):
    """Conservative: a harmful-direction row with undefined se still flags the gene (never hides harm). Documented, review-3 R3-3."""
    rows = sbp_replicated_gene("GA")
    rows = drop_rows(rows, "GA", "coronary_disease", "discovery")
    r = R("GA", "coronary_disease", "plof", "discovery", +0.9, 1e-6)
    r[5] = NAN
    rows.append(r)
    assert tier_map(rows, cfg)[("GA", "systolic_bp")] == "C"


def test_c11d_positive_control_still_requires_the_pipeline_discovery_hit(null_table, cfg, monkeypatch):
    """If the pipeline's own hit path does not return PCSK9/LDL, the control fails even though the direct row test passes."""
    empty = pd.DataFrame(columns=["gene", "trait", "beta", "se", "p", "n_carriers"])
    monkeypatch.setattr(stats, "discovery_hits", lambda d, c, mask="plof": empty)
    c = controls.positive_controls(validate(null_table), cfg)
    ldl = [p for p in c if p["id"] == "pcsk9_ldl_lower"][0]
    assert ldl["status"] == "FAIL"
    # controls not routed through the hit path (coronary direction, triglycerides) are unaffected
    assert [p["status"] for p in c if p["id"] != "pcsk9_ldl_lower"] == ["OK", "OK"]


# ======================================================================================
# (g) Amendment 2 (post-hoc after run 1): A2a contamination filter, A2b contamination bound
# ======================================================================================
def syn_row(gene, trait, beta, p):
    return R(gene, trait, "syn", "discovery", beta, p)


def contam_of(rows, cfg):
    return stats.contaminated_genes(T(rows), cfg)


# ---------------- A2a: definition of a contaminated gene ----------------
@pytest.mark.parametrize("trait", sorted(L_PANEL) + list(L_TRADEOFF) + ["triglycerides"])
@pytest.mark.parametrize("sign", [+1, -1])
def test_a2a_any_trait_and_either_direction_contaminates(trait, sign, cfg):
    got = contam_of([syn_row("GC", trait, sign * 0.5, 1e-9)], cfg)
    assert list(got) == ["GC"] and got["GC"][0]["trait"] == trait


@pytest.mark.parametrize("p,contaminated", [(1.0e-7, True), (1.89e-7, True), (1.9e-7, False), (2.0e-7, False), (1e-5, False)])
def test_a2a_threshold_is_the_discovery_threshold_strict(p, contaminated, cfg):
    assert (len(contam_of([syn_row("GC", "ldl", -0.5, p)], cfg)) == 1) == contaminated


@pytest.mark.parametrize("mask,cohort", [("plof", "discovery"), ("dmis", "discovery"), ("syn", "replication"), ("syn", "discovery_eur")])
def test_a2a_only_discovery_cohort_synonymous_rows_contaminate(mask, cohort, cfg):
    rows = [R("GC", "ldl", mask, cohort, -0.5, 1e-12)]
    assert contam_of(rows, cfg) == dict()


def test_a2a_records_every_hit_trait_and_counts_the_gene_once(cfg):
    got = contam_of([syn_row("GC", "ldl", -0.5, 1e-9), syn_row("GC", "bmi", +0.5, 1e-12), syn_row("GD", "fev1", 0.4, 1e-8)], cfg)
    assert sorted(got) == ["GC", "GD"] and sorted(h["trait"] for h in got["GC"]) == ["bmi", "ldl"]


# ---------------- A2a: excluded from ALL tiers and candidate lists ----------------
def _tier_fixture_genes():
    """One gene per tier letter, each with a clean twin, all with full trade-off screening."""
    rows = []
    rows += sbp_replicated_gene("XA") + sbp_replicated_gene("CA")                                   # A
    rows += hit_rows("XB", "fluid_intelligence") + null_tradeoffs("XB")                             # B
    rows += hit_rows("CB", "fluid_intelligence") + null_tradeoffs("CB")
    for g in ("XC", "CC"):                                                                          # C
        rows += drop_rows(sbp_replicated_gene(g), g, "coronary_disease", "discovery")
        rows.append(R(g, "coronary_disease", "plof", "discovery", +0.9, 1e-5))
    for g in ("XD", "CD"):                                                                          # D
        rows += hit_rows(g, "systolic_bp") + null_tradeoffs(g)
        rows.append(R(g, "hypertension", "plof", "replication", bad("hypertension"), 0.4))
    return rows


def test_a2a_contaminated_gene_is_excluded_in_every_tier_and_clean_twin_is_kept(cfg):
    rows = _tier_fixture_genes()
    before = tier_map(rows, cfg)
    assert [before[("XA", "systolic_bp")], before[("XB", "fluid_intelligence")], before[("XC", "systolic_bp")],
            before[("XD", "systolic_bp")]] == ["A", "B", "C", "D"]
    assert before[("CA", "systolic_bp")] == "A" and before[("CD", "systolic_bp")] == "D"
    for g, tr in (("XA", "fev1"), ("XB", "ldl"), ("XC", "coronary_disease"), ("XD", "hand_grip_strength")):
        rows.append(syn_row(g, tr, (+1 if g in ("XA", "XC") else -1) * 0.5, 1e-9))                 # both directions, various traits
    after = tier_map(rows, cfg)
    assert not [k for k in after if k[0] in ("XA", "XB", "XC", "XD")]
    assert after[("CA", "systolic_bp")] == "A" and after[("CB", "fluid_intelligence")] == "B"
    assert after[("CC", "systolic_bp")] == "C" and after[("CD", "systolic_bp")] == "D"


def test_a2a_exclusion_reaches_verdict_qualification_and_candidate_lists(cfg):
    from protscan.run import decide
    rows = sbp_replicated_gene("XA") + hit_rows("XB", "fluid_intelligence") + null_tradeoffs("XB")
    rows += [syn_row("XA", "bmi", +0.5, 1e-9), syn_row("XB", "bmi", -0.5, 1e-9)]
    hits = tiering.build_hits(T(rows), cfg)
    assert len(hits) == 0
    v = decide(CTRL_OK, hits)
    assert v["verdict"] == "KILL" and v["pass_genes"] == [] and v["lead_genes"] == []


def test_a2a_contaminated_control_genes_still_pass_their_controls_and_leave_the_tiers(null_table, cfg):
    """PCSK9, APOC3, ANGPTL4 and LDLR are not candidates: contamination must not corrupt control status."""
    df = null_table.copy()
    extra = [syn_row("PCSK9", "ldl", 0.3, 1e-12), syn_row("APOC3", "triglycerides", 0.3, 1e-12),
             syn_row("ANGPTL4", "bmi", -0.3, 1e-12)]
    df = pd.concat([df, pd.DataFrame(extra, columns=COLUMNS)], ignore_index=True)
    c = controls.run_controls(validate(df), cfg)
    assert [p["status"] for p in c["positive"]] == ["OK", "OK", "OK"]
    assert c["replication_sign"]["status"] == "OK"
    assert c["negative_synonymous"]["n_contaminated"] == 3
    hits = tiering.build_hits(validate(df), cfg)
    assert "PCSK9" not in set(hits["gene"])                                      # excluded from the candidate tiers


def test_a2a_contaminated_ldlr_does_not_touch_the_replication_sign_control(null_table, cfg):
    df = pd.concat([null_table, pd.DataFrame([syn_row("LDLR", "ldl", 0.3, 1e-12)], columns=COLUMNS)], ignore_index=True)
    assert controls.run_controls(validate(df), cfg)["replication_sign"]["status"] == "OK"


# ---------------- A2b: contamination bound is 0.1 percent of GENES tested ----------------
def _universe_table(n_genes, traits, contaminated):
    """n_genes x len(traits) null plof and syn rows; `contaminated` maps gene index to a list of (trait, sign) synonymous hits."""
    rng = np.random.default_rng(5)
    rows = []
    u = rng.uniform(0.01, 0.99, size=(n_genes, len(traits)))
    for i in range(n_genes):
        for j, t in enumerate(traits):
            rows.append(R("U%d" % i, t, "plof", "discovery", 0.0, float(u[i, j])))
            hit = [x for x in contaminated.get(i, []) if x[0] == t]
            if hit:
                rows.append(syn_row("U%d" % i, t, hit[0][1] * 0.5, 1e-9))
            else:
                rows.append(syn_row("U%d" % i, t, 0.0, float(1.0 - u[i, j] * 0.9)))
    return T(rows)


TEN_TRAITS = ["fev1", "bmi", "ldl", "fluid_intelligence", "reaction_time", "hand_grip_strength", "resting_heart_rate",
              "systolic_bp", "parental_lifespan", "walking_pace"]


@pytest.mark.parametrize("n_genes,k,ok", [(1000, 1, True), (1000, 2, False), (2000, 2, True), (2000, 3, False),
                                          (3000, 3, True), (3000, 4, False)])
def test_a2b_bound_is_exactly_one_in_a_thousand_genes(n_genes, k, ok, cfg):
    contam = dict((i, [(TEN_TRAITS[i % 10], +1 if i % 2 else -1)]) for i in range(k))     # both directions, different traits
    n = controls.negative_control(_universe_table(n_genes, TEN_TRAITS, contam), cfg)
    assert n["n_contaminated"] == k and n["n_genes_tested"] == n_genes
    assert L_SYN_CONTAM_MAX_FRACTION == 0.001 and (k * 1000 <= n_genes) == ok
    assert n["contamination_ok"] is ok
    assert n["status"] == ("OK" if ok else "FAIL")


def test_a2b_denominator_is_genes_tested_not_rows_and_a_gene_counts_once(cfg):
    """1,000 genes x 10 traits = 10,000 rows. One gene with hits on 3 traits is ONE contaminated gene: 1/1000, within the bound.
    Two genes are 0.2 percent of genes: over the bound (a rows denominator would say 0.02 percent)."""
    one = controls.negative_control(_universe_table(1000, TEN_TRAITS, dict([(0, [(TEN_TRAITS[0], 1), (TEN_TRAITS[1], -1), (TEN_TRAITS[2], 1)])])), cfg)
    assert one["n_contaminated"] == 1 and one["contamination_ok"] is True and one["status"] == "OK"
    two = controls.negative_control(_universe_table(1000, TEN_TRAITS, dict([(0, [(TEN_TRAITS[0], 1)]), (1, [(TEN_TRAITS[1], -1)])])), cfg)
    assert two["n_contaminated"] == 2 and two["contamination_ok"] is False and two["status"] == "FAIL"


def test_a2b_lambda_and_other_gates_are_unchanged_by_the_new_bound(null_table, cfg):
    n = neg_of(flat_syn_lambda(null_table, 1.101), cfg)
    assert n["lambda_ok"] is False and n["status"] == "FAIL" and n["n_contaminated"] == 0
    assert neg_of(null_table.iloc[:0], cfg)["status"] == "NOT RUN"


def test_a2b_old_zero_hit_gate_is_gone_but_reported_as_information(base_big, tmp_path):
    def edit(df):
        setrow(df, "SYNSYN1", "hand_grip_strength", "syn", "discovery", z=zgood("hand_grip_strength", 7))     # beneficial panel hit
    res, df = mutated_run(base_big, edit, tmp_path)
    n = res["controls"]["negative_synonymous"]
    assert n["n_syn_hit_genes"] == 1 and n["contamination_ok"] is True and res["controls"]["valid"] is True
    assert "syn_hits_ok" not in n and "syn_hits_max" not in n
    ok, _, o_verdict = oracle(df)
    assert ok and o_verdict == res["verdict"]["verdict"]


def test_a2b_end_to_end_systemic_contamination_is_kill_control_failed(scenario_runs):
    _, res, d = scenario_runs["contaminated_systemic"]
    n = res["controls"]["negative_synonymous"]
    assert n["contamination_ok"] is False and n["n_contaminated"] > n["n_genes_tested"] * 0.001
    assert res["verdict"]["verdict"] == "KILL" and res["verdict"]["reason"].startswith("control_failed")
    assert res["contamination"]["within_bound"] is False
    assert "FAIL" in (d / "out" / "report.md").read_text()


# ---------------- A2: scenarios, oracle, and the required disclosure ----------------
def test_a2_contaminated_tier_a_gene_is_excluded_and_verdict_degrades_pass_to_lead(scenario_runs):
    _, res, _ = scenario_runs["contaminated_tierA"]
    assert res["controls"]["valid"] is True
    assert ("SYNPASS1", "systolic_bp") not in result_tiers(res) and res["verdict"]["pass_genes"] == []
    assert res["verdict"]["verdict"] == "LEAD" and res["verdict"]["contaminated_excluded"] == ["SYNPASS1"]
    row = [g for g in res["contamination"]["genes"] if g["gene"] == "SYNPASS1"][0]
    assert row["excluded_plof_hit_traits"] == ["systolic_bp"] and row["syn_hits"][0]["trait"] == "systolic_bp"


def test_a2_one_contaminated_null_gene_below_the_bound_leaves_the_verdict_alone(scenario_runs):
    _, res, _ = scenario_runs["contaminated_minor"]
    assert res["controls"]["negative_synonymous"]["n_contaminated"] == 1 and res["controls"]["valid"] is True
    assert res["verdict"]["verdict"] == "PASS" and res["verdict"]["pass_genes"] == ["SYNPASS1"]


@pytest.mark.parametrize("name", ["pass", "contaminated_tierA", "contaminated_systemic", "contaminated_minor", "broken_syn_hit"])
def test_a2_verdict_block_label_report_and_cli_carry_the_disclosure(name, scenario_runs):
    _, res, d = scenario_runs[name]
    assert res["verdict"]["amendment"] == L_A2_LABEL and L_A2_LABEL in res["verdict"]["label"]
    assert res["contamination"]["amendment"] == L_A2_LABEL
    text = (d / "out" / "report.md").read_text()
    assert "Rules active: **%s**" % L_A2_LABEL in text and "never a clean preregistered pass" in text
    assert "## Contaminated genes (A2a)" in text
    saved = json.loads((d / "out" / "protective-scan.json").read_text())
    assert saved["verdict"]["amendment"] == L_A2_LABEL


def test_a2_cli_stdout_states_the_amendment(synth, tmp_path):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    d = tmp_path / "data"
    write_tables(synth.make_synthetic("pass", n_genes=N_GENES), d)
    out = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(CONFIG_PATH), "--data", str(d),
                          "--out", str(tmp_path / "o" / "r.json")], capture_output=True, text=True, env=env, cwd=str(ROOT))
    assert out.returncode == 0 and L_A2_LABEL in out.stdout


def test_a2_report_is_not_silent_about_excluded_genes(scenario_runs):
    _, res, d = scenario_runs["contaminated_tierA"]
    text = (d / "out" / "report.md").read_text()
    section = text.split("## Contaminated genes (A2a)")[1].split("## Baseline ladder")[0]
    assert "SYNPASS1" in section and "systolic_bp" in section


# ---------------- review-4 residual gaps (non-strict xfail) ----------------
@pytest.mark.xfail(strict=False, reason="R4-1: a pLoF hit gene with no synonymous row cannot be assessed by A2a and is not flagged")
def test_a2a_hit_gene_without_any_synonymous_row_is_flagged_as_unassessed(cfg):
    rows = sbp_replicated_gene("GA") + [syn_row("OTHER", "ldl", 0.0, 0.5)]
    hits = tiering.build_hits(T(rows), cfg)
    assert list(hits["tier"]) == ["A"]                                       # present in the tiers ...
    assert "syn_assessed" in hits.columns and bool(hits["syn_assessed"].iloc[0]) is False     # ... but flagged


@pytest.mark.xfail(strict=False, reason="R4-3: results ledger_entry header still reads 'with AMENDMENT 1' under Amendment 2 rules")
def test_a2_ledger_entry_header_names_amendment_2(scenario_runs):
    _, res, d = scenario_runs["pass"]
    assert "amendment 2" in res["ledger_entry"].lower()
    assert "amendment 2" in (d / "out" / "report.md").read_text().split("## Verdict")[0].lower()
