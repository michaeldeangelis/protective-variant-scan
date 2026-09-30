"""Amendment 2 (post-hoc after run 1): A2a contamination filter, A2b contamination bound. Synthetic data only."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protscan import controls, stats, tiering
from protscan.run import AMENDMENT_2, run_pipeline
from protscan.schema import load_config, validate

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "prereg.yaml"
cfg = load_config(CONFIG)


def row(gene, trait, mask="plof", cohort="discovery", beta=-0.5, p=1e-9, se=0.05, source=None):
    src = source or ("finngen_test" if cohort == "replication" else "genebass_test")
    return dict(gene=gene, trait=trait, mask=mask, cohort=cohort, beta=beta, se=se, p=p,
                n_carriers=100, n_total=400000, source=src)


def table(rows):
    return validate(pd.DataFrame(rows))


def universe(n_genes=10_000, n_contaminated=0, seed=4):
    """n_genes null (gene, fev1) pairs with plof and syn rows; the first n_contaminated syn rows are significant."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(0.05, 1.0, size=2 * n_genes)
    rows = [row(f"S{i}", "fev1", "plof", beta=0.0, p=float(u[i])) for i in range(n_genes)]
    for i in range(n_genes):
        hit = i < n_contaminated
        rows.append(row(f"S{i}", "fev1", "syn", beta=(0.3 if i % 2 else -0.3) if hit else 0.0,
                        p=1e-9 if hit else float(u[n_genes + i])))
    return rows


@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("make_synthetic_a2", ROOT / "scripts" / "make_synthetic.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def make(synth, tmp_path_factory):
    def _make(scenario, name=None):
        d = tmp_path_factory.mktemp(name or scenario)
        synth.write(synth.make_synthetic(scenario, n_genes=500), d)
        return d
    return _make


def run(d):
    return run_pipeline(CONFIG, d, Path(d) / "out" / "protective-scan.json")


def tier_genes(res):
    return {t: {x["gene"] for x in res["tiers"][t]} for t in "ABCD"}


# ---------------- config ----------------
def test_config_replaces_zero_hit_gate_with_contamination_bound():
    assert cfg.syn_contaminated_max_fraction == 0.001
    import yaml
    assert "syn_hits_max" not in yaml.safe_load(CONFIG.read_text())["controls"]
    assert not hasattr(cfg, "syn_hits_max")


# ---------------- A2a ----------------
def test_contaminated_genes_either_direction_any_trait_discovery_syn_only():
    df = table([
        row("A", "fev1", "syn", beta=+0.3, p=1e-9),                    # beneficial direction
        row("B", "fev1", "syn", beta=-0.3, p=1e-9),                    # harmful direction
        row("C", "coronary_disease", "syn", beta=+0.3, p=1e-9),        # trade-off trait
        row("D", "triglycerides", "syn", beta=-0.3, p=1e-9),           # control trait
        row("E", "fev1", "syn", beta=0.3, p=cfg.discovery_p),          # p == threshold: not a hit
        row("F", "fev1", "plof", beta=0.3, p=1e-9),                    # not the synonymous mask
        row("G", "fev1", "dmis", beta=0.3, p=1e-9),
        row("H", "fev1", "syn", cohort="replication", beta=0.3, p=1e-9),
        row("I", "fev1", "syn", cohort="discovery_eur", beta=0.3, p=1e-9),
        row("J", "fev1", "syn", beta=0.3, p=1e-9),
        row("J", "bmi", "syn", beta=-0.3, p=1e-12),                    # two traits for one gene
    ])
    c = stats.contaminated_genes(df, cfg)
    assert sorted(c) == ["A", "B", "C", "D", "J"]
    assert [h["trait"] for h in c["J"]] == ["bmi", "fev1"]             # listed with all traits, most significant first


def test_contaminated_gene_is_excluded_from_every_tier_but_others_remain():
    rows = [row("BAD", "fev1", beta=0.5), row("BAD", "fev1", "syn", beta=0.3, p=1e-9),
            row("OK1", "fev1", beta=0.5), row("OK1", "fev1", "syn", beta=0.3, p=0.5)]
    h = tiering.build_hits(table(rows), cfg)
    assert list(h["gene"]) == ["OK1"]
    # the unfiltered discovery rule (simplest rung) still sees it
    assert set(stats.discovery_hits(table(rows), cfg)["gene"]) == {"BAD", "OK1"}


def test_contaminated_gene_never_reaches_a_verdict_list():
    rows = [row("BAD", "fev1", beta=0.5), row("BAD", "fev1", "dmis", beta=0.1, p=0.2),
            row("BAD", "fev1", "syn", beta=-0.3, p=1e-9)] + [row("BAD", t, beta=0.01, p=0.9) for t in cfg.tradeoff]
    assert tiering.build_hits(table(rows), cfg).empty


def test_positive_control_gene_can_be_contaminated_without_failing_the_control():
    df = table([row("PCSK9", "ldl", beta=-0.8, p=1e-30), row("PCSK9", "ldl", "syn", beta=0.3, p=1e-9),
                row("PCSK9", "coronary_disease", beta=-0.1, p=0.3), row("APOC3", "triglycerides", beta=-1.0, p=1e-20)])
    hits = tiering.build_hits(df, cfg)
    assert "PCSK9" not in set(hits["gene"])                            # excluded from tiers and candidate lists
    assert [c["status"] for c in controls.positive_controls(df, cfg, hits)] == ["OK", "OK", "OK"]


# ---------------- A2b ----------------
@pytest.mark.parametrize("k,status", [(0, "OK"), (1, "OK"), (10, "OK"), (11, "FAIL")])
def test_contamination_bound_is_one_tenth_of_a_percent_of_genes_tested(k, status):
    c = controls.negative_control(table(universe(10_000, k)), cfg)
    assert c["status"] == status and c["n_genes_tested"] == 10_000 and c["n_contaminated"] == k
    assert c["contaminated_fraction"] == pytest.approx(k / 10_000) and c["contaminated_max_fraction"] == 0.001
    assert c["contamination_ok"] == (status == "OK")


def test_contamination_check_is_a_fraction_not_a_count():
    # 11 contaminated genes are fine among 20,000 genes tested, not among 10,000
    assert controls.negative_control(table(universe(20_000, 11)), cfg)["status"] == "OK"
    assert controls.negative_control(table(universe(10_000, 11)), cfg)["status"] == "FAIL"


def test_other_negative_control_gates_are_unchanged_and_independent():
    # small contamination but inflated lambda: still FAIL on lambda
    rng = np.random.default_rng(9)
    z = rng.standard_normal(20_000) * 1.5
    from scipy.stats import norm
    p = 2 * norm.sf(np.abs(z))
    rows = [row(f"S{i}", "fev1", "plof", beta=0.0, p=0.5) for i in range(20_000)]
    rows += [row(f"S{i}", "fev1", "syn", beta=0.0, p=float(min(p[i], 0.999))) for i in range(20_000)]
    c = controls.negative_control(table(rows), cfg)
    assert not c["lambda_ok"] and c["status"] == "FAIL"


def test_old_beneficial_direction_count_is_kept_but_no_longer_gates():
    c = controls.negative_control(table(universe(10_000, 4)), cfg)          # 4 contaminated, mixed directions
    assert c["status"] == "OK" and "n_syn_hit_genes" in c and "syn_hit_genes" in c
    assert c["n_syn_hit_genes"] <= c["n_contaminated"]


# ---------------- scenarios (synthetic only) ----------------
def test_scenario_tierA_gene_contaminated_is_excluded_and_verdict_changes(make):
    base = run(make("pass"))
    r = run(make("contaminated_tierA"))
    assert base["verdict"]["verdict"] == "PASS" and "SYNPASS1" in tier_genes(base)["A"]
    assert r["controls"]["valid"] and r["controls"]["negative_synonymous"]["contamination_ok"]
    assert "SYNPASS1" not in set().union(*tier_genes(r).values())
    assert r["verdict"]["verdict"] == "LEAD" and r["verdict"]["lead_genes"] == ["SYNLEAD1"] and r["verdict"]["pass_genes"] == []
    g = {x["gene"]: x for x in r["contamination"]["genes"]}
    assert list(g) == ["SYNPASS1"]
    assert [h["trait"] for h in g["SYNPASS1"]["syn_hits"]] == ["systolic_bp"]
    assert g["SYNPASS1"]["excluded_plof_hit_traits"] == ["systolic_bp"]
    assert r["verdict"]["contaminated_excluded"] == ["SYNPASS1"]


def test_scenario_contamination_above_bound_kills(make):
    r = run(make("contaminated_systemic"))
    n = r["controls"]["negative_synonymous"]
    assert r["verdict"]["verdict"] == "KILL" and r["verdict"]["reason"] == "control_failed: negative_synonymous"
    assert not n["contamination_ok"] and n["lambda_ok"] and n["contaminated_fraction"] > n["contaminated_max_fraction"]
    assert r["contamination"]["within_bound"] is False and r["contamination"]["n_contaminated"] == n["n_contaminated"]
    dirs = {np.sign(h["beta"]) for g in r["contamination"]["genes"] for h in g["syn_hits"]}
    assert dirs == {-1.0, 1.0}                                    # both directions are counted


def test_scenario_contamination_below_bound_gives_normal_verdict(make):
    r = run(make("contaminated_minor"))
    assert r["controls"]["valid"] and r["contamination"]["within_bound"] is True
    assert [g["gene"] for g in r["contamination"]["genes"]] == ["SYNG00000"]
    assert r["verdict"]["verdict"] == "PASS" and r["verdict"]["pass_genes"] == ["SYNPASS1"]
    assert "SYNG00000" not in set().union(*tier_genes(r).values())


def test_bound_scales_with_the_number_of_genes_tested(synth):
    for n in (500, 3000):
        t = synth.make_synthetic("contaminated_systemic", n_genes=n)
        d = t["burden_synthetic_discovery"]
        s = d[(d["mask"] == "syn") & (d.cohort == "discovery") & (d.p < cfg.discovery_p)]
        assert s.gene.nunique() / d.gene.nunique() > cfg.syn_contaminated_max_fraction, n


# ---------------- reporting ----------------
def test_amendment_2_label_in_verdict_block_json_report_and_cli(make, tmp_path):
    d = make("pass")
    r = run(d)
    assert r["verdict"]["amendment"] == AMENDMENT_2 == "Amendment 2 (post-hoc after run 1)"
    assert f"[{AMENDMENT_2}]" in r["verdict"]["label"]
    md = (d / "out" / "report.md").read_text()
    assert AMENDMENT_2 in md and "never a clean preregistered pass" in md
    assert json.loads((d / "out" / "protective-scan.json").read_text())["verdict"]["amendment"] == AMENDMENT_2
    p = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(CONFIG), "--data", str(d),
                        "--out", str(tmp_path / "o" / "r.json")], capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 0 and AMENDMENT_2 in p.stdout


def test_report_lists_contaminated_genes_count_fraction_limit_and_old_count(make):
    d = make("contaminated_tierA")
    r = run(d)
    md = (d / "out" / "report.md").read_text()
    n = r["controls"]["negative_synonymous"]
    assert "## Contaminated genes (A2a): 1 of" in md
    assert "| SYNPASS1 | systolic_bp (beta=" in md and "excluded" in md
    assert f"1 of {n['n_genes_tested']} genes (fraction {n['contaminated_fraction']:.5f}, limit 0.0010)" in md
    assert "old gate, informational" in md and "replaced by A2b" in md
    assert "A2c (manual" in md


def test_report_flags_contamination_over_bound(make):
    d = make("contaminated_systemic")
    run(d)
    md = (d / "out" / "report.md").read_text()
    assert "exceeds the A2b bound: control FAIL" in md and "| synonymous contaminated genes (A2b)" in md
    assert "| FAIL |" in md


def test_contamination_block_present_when_nothing_is_contaminated(make):
    r = run(make("pass"))
    assert r["contamination"]["n_contaminated"] == 0 and r["contamination"]["genes"] == []
    assert r["verdict"]["contaminated_excluded"] == []
