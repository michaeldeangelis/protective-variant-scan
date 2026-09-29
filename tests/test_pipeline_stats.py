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
        *[row(g, t, beta=0.01, p=0.9) for g in ("LPL", "X1", "X2", "X3") for t in list(cfg.tradeoff)[:5]],
    )
    h = tiering.build_hits(df, cfg).set_index("gene")
    assert not h.loc["LPL", "qualifies"] and h.loc["LPL", "lipid_gene"]
    assert not h.loc["X1", "qualifies"]
    assert not h.loc["X2", "qualifies"]
    assert h.loc["X3", "qualifies"]


def test_build_hits_empty():
    df = table(row("A", "fev1", beta=0.5, p=0.3))
    assert tiering.build_hits(df, cfg).empty


def _syn_universe(n, seed=2, p_power=1.0):
    """n (gene, trait) pairs with null plof and syn rows, so coverage is 1.0 and syn rows == n."""
    rng = np.random.default_rng(seed)
    rows = []
    for mask in ("plof", "syn"):
        u = rng.uniform(size=n) ** (p_power if mask == "syn" else 1.0)
        rows += [row(f"S{i}", "fev1", mask=mask, beta=0.0, p=float(u[i])) for i in range(n)]
    return rows


def test_negative_control_thresholds():
    base = _syn_universe(10_000)
    ok = controls.negative_control(table(*base), cfg)
    assert ok["status"] == "OK" and ok["coverage"] == 1.0
    hit = controls.negative_control(table(*base, row("H", "fev1", mask="syn", beta=0.5, p=1e-9)), cfg)
    assert hit["status"] == "FAIL" and hit["n_syn_hit_genes"] == 1 and hit["lambda_ok"]
    # a hit in the harmful direction does not count
    wrong = controls.negative_control(table(*base, row("H", "fev1", mask="syn", beta=-0.5, p=1e-9)), cfg)
    assert wrong["status"] == "OK"
    assert controls.negative_control(table(*_syn_universe(10_000, p_power=1.5)), cfg)["status"] == "FAIL"
    assert controls.negative_control(table(row("A", "fev1")), cfg)["status"] == "NOT RUN"


def test_positive_controls():
    good = table(row("PCSK9", "ldl", beta=-0.8), row("PCSK9", "coronary_disease", beta=-0.1, p=0.3),
                 row("APOC3", "triglycerides", beta=-1.0))
    assert [c["status"] for c in controls.positive_controls(good, cfg)] == ["OK", "OK", "OK"]
    weak_ldl = table(row("PCSK9", "ldl", beta=-0.8, p=1e-5), row("PCSK9", "coronary_disease", beta=+0.1, p=0.3),
                     row("ANGPTL4", "triglycerides", beta=+1.0))
    assert [c["status"] for c in controls.positive_controls(weak_ldl, cfg)] == ["FAIL", "FAIL", "FAIL"]
    assert [c["status"] for c in controls.positive_controls(table(row("Q", "ldl")), cfg)] == ["NOT RUN"] * 3


# ---- C1: Tier A needs >= min_tradeoffs_screened of 9 outcomes screened (plof) ----
def _sbp_replicated(gene="G", k=9, rep_outcomes=(), disc_mask="plof", disc_cohort="discovery", extra=()):
    outs = list(cfg.tradeoff)
    rows = [row(gene, "systolic_bp"), row(gene, "systolic_bp", mask="dmis", beta=-0.1, p=0.3),
            row(gene, "hypertension", cohort="replication", beta=-0.3, p=0.01)]
    rows += [row(gene, t, mask=disc_mask, cohort=disc_cohort, beta=0.01, p=0.9) for t in outs[:k]]
    rows += [row(gene, t, cohort="replication", beta=0.01, p=0.9) for t in rep_outcomes]
    return table(*rows, *extra)


def test_min_tradeoffs_screened_config():
    assert cfg.min_tradeoffs_screened == 5


@pytest.mark.parametrize("k,tier", [(0, "B"), (3, "B"), (4, "B"), (5, "A"), (6, "A"), (9, "A")])
def test_tier_A_boundary_on_tradeoffs_screened(k, tier):
    h = tiering.build_hits(_sbp_replicated(k=k), cfg).iloc[0]
    assert h["n_tradeoff_tested"] == k and h["tier"] == tier
    assert bool(h["tradeoff_unscreened"]) == (k < 5)


def test_screened_counts_distinct_outcomes_across_cohorts_plof_only():
    outs = list(cfg.tradeoff)
    same = tiering.build_hits(_sbp_replicated(k=3, rep_outcomes=outs[:3]), cfg).iloc[0]
    assert same["n_tradeoff_tested"] == 3 and same["tier"] == "B"
    union = tiering.build_hits(_sbp_replicated(k=3, rep_outcomes=outs[3:5]), cfg).iloc[0]
    assert union["n_tradeoff_tested"] == 5 and union["tier"] == "A"
    dmis_only = tiering.build_hits(_sbp_replicated(k=9, disc_mask="dmis"), cfg).iloc[0]
    assert dmis_only["n_tradeoff_tested"] == 0 and dmis_only["tier"] == "B"
    eur_only = tiering.build_hits(_sbp_replicated(k=9, disc_cohort="discovery_eur"), cfg).iloc[0]
    assert eur_only["n_tradeoff_tested"] == 0 and eur_only["tier"] == "B"


def test_unscreened_cap_does_not_touch_tiers_C_and_D():
    adv = [row("G", "coronary_disease", beta=0.9, p=1e-6)]
    c = tiering.build_hits(_sbp_replicated(k=0, extra=adv), cfg).iloc[0]
    assert c["tier"] == "C"                      # adverse found: stays C although only 1 outcome screened
    failed = _sbp_replicated(k=2)
    failed.loc[failed.trait == "hypertension", ["beta", "p"]] = [0.3, 0.01]
    assert tiering.build_hits(failed, cfg).iloc[0]["tier"] == "D"


def test_report_labels_unscreened():
    from protscan.report import _tier_table
    rec = tiering.build_hits(_sbp_replicated(k=4), cfg).iloc[0].to_dict()
    rec.update(loeuf=None, rep_p_onesided=0.005)
    assert "trade-off unscreened" in "\n".join(_tier_table([rec], 9))
