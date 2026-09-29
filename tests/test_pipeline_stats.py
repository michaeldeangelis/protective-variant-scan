from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protscan import controls, stats, tiering
from protscan.schema import load_config, validate

CONFIG = Path(__file__).resolve().parent.parent / "config" / "prereg.yaml"
cfg = load_config(CONFIG)
THR = cfg.discovery_p


def row(gene, trait, mask="plof", cohort="discovery", beta=-0.5, p=1e-9, se=0.05):
    return dict(gene=gene, trait=trait, mask=mask, cohort=cohort, beta=beta, se=se, p=p,
                n_carriers=100, n_total=400000,
                source="finngen_test" if cohort == "replication" else "genebass_test")


def table(*rows):
    return validate(pd.DataFrame(rows))


def test_lambda_gc_uniform_and_inflated():
    rng = np.random.default_rng(1)
    p = rng.uniform(size=200_000)
    assert stats.lambda_gc(p) == pytest.approx(1.0, abs=0.02)
    z = rng.standard_normal(200_000) * 1.3
    from scipy.stats import norm
    assert stats.lambda_gc(2 * norm.sf(np.abs(z))) == pytest.approx(1.69, abs=0.05)


def test_one_sided_p():
    assert stats.one_sided_p(-1.0, 0.04, -1) == pytest.approx(0.02)
    assert stats.one_sided_p(1.0, 0.04, -1) == pytest.approx(0.98)


def test_discovery_direction_of_benefit_and_strict_threshold():
    df = table(
        row("A", "reaction_time", beta=-0.5),          # lower time = beneficial: hit
        row("B", "reaction_time", beta=+0.5),          # slower: not a hit
        row("C", "fluid_intelligence", beta=+0.5),     # higher: hit
        row("D", "fluid_intelligence", beta=-0.5),     # lower: not a hit
        row("E", "fluid_intelligence", beta=+0.5, p=THR),        # p == threshold: not a hit (strict)
        row("F", "coronary_disease", beta=-0.5),       # not a panel trait
        row("G", "fluid_intelligence", mask="dmis", beta=0.5),   # wrong mask
    )
    assert set(stats.discovery_hits(df, cfg).gene) == {"A", "C"}


def test_replication_rules():
    df = table(
        row("G1", "systolic_bp"), row("G2", "systolic_bp"), row("G3", "systolic_bp"), row("G4", "systolic_bp"),
        row("G5", "fluid_intelligence", beta=0.5), row("G6", "systolic_bp"),
        row("G1", "hypertension", cohort="replication", beta=-0.3, p=0.08),   # one-sided 0.04: replicates
        row("G2", "hypertension", cohort="replication", beta=-0.3, p=0.12),   # one-sided 0.06: fails
        row("G3", "hypertension", cohort="replication", beta=+0.3, p=0.001),  # opposite direction: fails
        row("G6", "obesity", cohort="replication", beta=-0.3, p=1e-9),
        row("G5", "fluid_intelligence", cohort="replication", beta=0.5, p=1e-9),  # no declared proxy: ignored
    )
    hits = stats.discovery_hits(df, cfg)
    rs = stats.replication_status(hits, df, cfg).set_index(hits["gene"])["rep_status"]
    assert rs["G1"] == "replicated"
    assert rs["G2"] == "failed"
    assert rs["G3"] == "failed"
    assert rs["G4"] == "gene_untested"      # hypertension is in the replication source, G4 is not
    assert rs["G6"] == "gene_untested"      # a row for a different trait (obesity) is not a proxy for SBP
    assert rs["G5"] == "trait_missing"      # Tier B ceiling even though the same trait is in the replication table


def test_replication_gene_untested_vs_trait_missing():
    df = table(row("G1", "bmi"), row("G2", "ldl"),
               row("Z", "obesity", cohort="replication", beta=-0.3, p=1e-6))
    hits = stats.discovery_hits(df, cfg)
    rs = stats.replication_status(hits, df, cfg).set_index(hits["gene"])["rep_status"]
    assert rs["G1"] == "gene_untested"      # obesity is in the replication source, G1 is not
    assert rs["G2"] == "trait_missing"      # neither ldl nor hypercholesterolemia in the replication source


def test_replication_same_trait_accepted_for_proxy_traits():
    df = table(row("G1", "ldl"), row("G1", "ldl", cohort="replication", beta=-0.3, p=1e-6))
    hits = stats.discovery_hits(df, cfg)
    r = stats.replication_status(hits, df, cfg).iloc[0]
    assert r["rep_status"] == "replicated" and r["rep_trait"] == "ldl"


def test_adverse_threshold_direction_and_mask():
    a = cfg.tradeoff_p                                   # 0.05/9
    df = table(
        row("A", "coronary_disease", beta=+0.3, p=2 * a * 0.99),   # one-sided just under threshold, harmful
        row("B", "coronary_disease", beta=+0.3, p=2 * a * 1.01),   # just over
        row("C", "coronary_disease", beta=-0.3, p=1e-9),           # protective: not adverse
        row("D", "dementia", mask="dmis", beta=+0.3, p=1e-9),      # dmis mask not screened
        row("E", "dementia", cohort="replication", beta=+0.3, p=1e-9),  # replication cohort is screened
    )
    adv = stats.adverse_tradeoffs(list("ABCDE"), df, cfg)
    assert [bool(adv[g]) for g in "ABCDE"] == [True, False, False, False, True]


def test_mask_status():
    df = table(
        row("A", "fev1", beta=0.5), row("A", "fev1", mask="dmis", beta=0.1),
        row("B", "fev1", beta=0.5), row("B", "fev1", mask="dmis", beta=-0.1),
        row("C", "fev1", beta=0.5),
    )
    hits = stats.discovery_hits(df, cfg)
    ms = dict(zip(hits["gene"], stats.mask_status(hits, df, cfg)))
    assert ms == {"A": "consistent", "B": "inconsistent", "C": "missing"}


def test_assign_tier():
    adv = [{"trait": "x"}]
    assert tiering.assign_tier([], "replicated") == "A"
    assert tiering.assign_tier([], "trait_missing") == "B"
    assert tiering.assign_tier([], "gene_untested") == "B"
    assert tiering.assign_tier(adv, "trait_missing") == "C"
    assert tiering.assign_tier(adv, "replicated") == "C"
    assert tiering.assign_tier([], "failed") == "D"
    assert tiering.assign_tier(adv, "failed") == "C"


def test_qualifies_excludes_lipid_domain_and_mask_missing():
    df = table(
        row("LPL", "systolic_bp"), row("LPL", "systolic_bp", mask="dmis"),
        row("X1", "bmi"), row("X1", "bmi", mask="dmis"),               # metabolic domain
        row("X2", "systolic_bp"),                                       # dmis missing
        row("X3", "systolic_bp"), row("X3", "systolic_bp", mask="dmis"),
    )
    h = tiering.build_hits(df, cfg).set_index("gene")
    assert not h.loc["LPL", "qualifies"] and h.loc["LPL", "lipid_gene"]
    assert not h.loc["X1", "qualifies"]
    assert not h.loc["X2", "qualifies"]
    assert h.loc["X3", "qualifies"]


def test_build_hits_empty():
    df = table(row("A", "fev1", beta=0.5, p=0.3))
    assert tiering.build_hits(df, cfg).empty


def test_negative_control_thresholds():
    rng = np.random.default_rng(2)
    n = 5000
    base = [row(f"S{i}", "fev1", mask="syn", beta=0.0, p=float(rng.uniform())) for i in range(n)]
    ok = controls.negative_control(table(*base), cfg)
    assert ok["status"] == "OK"
    hit = controls.negative_control(table(*base, row("H", "fev1", mask="syn", beta=0.5, p=1e-9)), cfg)
    assert hit["status"] == "FAIL" and hit["n_syn_hit_genes"] == 1 and hit["lambda_ok"]
    # a hit in the harmful direction does not count
    wrong = controls.negative_control(table(*base, row("H", "fev1", mask="syn", beta=-0.5, p=1e-9)), cfg)
    assert wrong["status"] == "OK"
    inflated = [row(f"S{i}", "fev1", mask="syn", beta=0.0, p=float(rng.uniform() ** 1.5)) for i in range(n)]
    assert controls.negative_control(table(*inflated), cfg)["status"] == "FAIL"
    assert controls.negative_control(table(row("A", "fev1")), cfg)["status"] == "NOT RUN"


def test_positive_controls():
    good = table(row("PCSK9", "ldl", beta=-0.8), row("PCSK9", "coronary_disease", beta=-0.1, p=0.3),
                 row("APOC3", "triglycerides", beta=-1.0))
    assert [c["status"] for c in controls.positive_controls(good, cfg)] == ["OK", "OK", "OK"]
    weak_ldl = table(row("PCSK9", "ldl", beta=-0.8, p=1e-5), row("PCSK9", "coronary_disease", beta=+0.1, p=0.3),
                     row("ANGPTL4", "triglycerides", beta=+1.0))
    assert [c["status"] for c in controls.positive_controls(weak_ldl, cfg)] == ["FAIL", "FAIL", "FAIL"]
    assert [c["status"] for c in controls.positive_controls(table(row("Q", "ldl")), cfg)] == ["NOT RUN"] * 3
