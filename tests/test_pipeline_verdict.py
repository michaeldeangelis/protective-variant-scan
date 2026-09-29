"""Amendment 1 clarifications C2 (domain gate on PASS/LEAD) and C3 (pinned lipid-pathway list)."""
from pathlib import Path

import pandas as pd
import pytest

from protscan import tiering
from protscan.run import decide
from protscan.schema import load_config, validate

CONFIG = Path(__file__).resolve().parent.parent / "config" / "prereg.yaml"
cfg = load_config(CONFIG)
VALID = {"valid": True, "failed": [], "not_run": []}


def _r(gene, trait, mask="plof", cohort="discovery", beta=-0.5, p=1e-9):
    src = "finngen_test" if cohort == "replication" else "genebass_test"
    return dict(gene=gene, trait=trait, mask=mask, cohort=cohort, beta=beta, se=0.05, p=p,
                n_carriers=100, n_total=400000, source=src)


def _gene_rows(gene, trait, with_replication):
    s = cfg.panel_sign[trait]
    rows = [_r(gene, trait, beta=0.5 * s), _r(gene, trait, mask="dmis", beta=0.1 * s, p=0.2)]
    rows += [_r(gene, t, beta=0.01, p=0.9) for t in cfg.tradeoff]
    if with_replication:
        rep = cfg.proxies[trait]["trait"] if trait in cfg.proxies else trait   # same-trait row for non-proxy traits
        rows.append(_r(gene, rep, cohort="replication", beta=0.5 * s, p=1e-8))
    return rows


def _verdict(trait, with_replication=True):
    hits = tiering.build_hits(validate(pd.DataFrame(_gene_rows("XG", trait, with_replication))), cfg)
    return decide(VALID, hits), hits.iloc[0]


def test_only_sbp_among_proxy_traits_gives_pass():
    v, h = _verdict("systolic_bp")
    assert (h["tier"], h["qualifies"], v["verdict"]) == ("A", True, "PASS")


@pytest.mark.parametrize("trait", ["ldl", "bmi"])
def test_ldl_and_bmi_reach_tier_A_but_never_qualify(trait):
    v, h = _verdict(trait)
    assert h["tier"] == "A" and not h["qualifies"]
    assert v["verdict"] == "KILL" and v["pass_genes"] == [] and v["lead_genes"] == []


@pytest.mark.parametrize("with_replication", [False, True])
def test_parental_lifespan_never_qualifies(with_replication):
    v, h = _verdict("parental_lifespan", with_replication)
    assert h["tier"] == "B" and not h["qualifies"]          # no declared proxy: capped at B; metabolic: not qualifying
    assert v["verdict"] == "KILL"


@pytest.mark.parametrize("trait", [t for t, d in cfg.panel_domain.items() if d == "cognitive"])
@pytest.mark.parametrize("with_replication", [False, True])
def test_cognitive_traits_cap_at_lead(trait, with_replication):
    v, h = _verdict(trait, with_replication)
    assert h["tier"] == "B" and h["qualifies"]
    assert v["verdict"] == "LEAD" and v["pass_genes"] == []


@pytest.mark.parametrize("trait", ["hand_grip_strength", "fev1", "walking_pace", "resting_heart_rate"])
def test_physical_traits_without_proxy_cap_at_lead(trait):
    v, h = _verdict(trait)
    assert h["tier"] == "B" and v["verdict"] == "LEAD"


def test_no_panel_trait_other_than_sbp_can_pass():
    passing = [t for t in cfg.panel if _verdict(t)[0]["verdict"] == "PASS"]
    assert passing == ["systolic_bp"]


# C3: list fixed at commit 4073b89
LIPID_PATHWAY_GENES_4073B89 = [
    "PCSK9", "ANGPTL3", "ANGPTL4", "ANGPTL8", "APOC3", "APOA5", "APOB", "APOE", "LDLR", "LDLRAP1", "LPL",
    "LPA", "CETP", "LIPC", "MTTP", "NPC1L1", "ABCG5", "ABCG8", "HMGCR", "SORT1", "LIPG", "GPIHBP1", "LMF1",
    "ANGPTL1",
]


def test_lipid_pathway_list_pinned_to_4073b89():
    assert cfg.lipid_genes == frozenset(LIPID_PATHWAY_GENES_4073B89)
    assert len(LIPID_PATHWAY_GENES_4073B89) == 24
