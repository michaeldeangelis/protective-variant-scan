import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from protscan.run import run_pipeline

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "prereg.yaml"
N_GENES = 500


@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("make_synthetic", ROOT / "scripts" / "make_synthetic.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def make(synth, tmp_path_factory):
    def _make(scenario, name=None, edit=None):
        d = tmp_path_factory.mktemp(name or scenario)
        tables = synth.make_synthetic(scenario, n_genes=N_GENES)
        if edit:
            tables = edit(tables)
        synth.write(tables, d)
        return d
    return _make


def run(data_dir):
    return run_pipeline(CONFIG, data_dir, Path(data_dir) / "out" / "protective-scan.json")


def genes(res, tier):
    return {r["gene"] for r in res["tiers"][tier]}


@pytest.fixture(scope="module")
def pass_res(make):
    return run(make("pass"))


def test_pass_scenario(pass_res):
    r = pass_res
    assert r["verdict"]["verdict"] == "PASS"
    assert r["controls"]["valid"]
    assert r["controls"]["negative_synonymous"]["lambda_gc"] < 1.10
    assert r["verdict"]["pass_genes"] == ["SYNPASS1"]
    assert "SYNPASS1" in genes(r, "A") and "PCSK9" in genes(r, "A")
    assert genes(r, "B") == {"SYNLEAD1", "SYNMASK1"}
    assert genes(r, "C") == {"SYNADV1"}
    assert genes(r, "D") == {"SYNFAIL1"}


def test_mask_inconsistent_gene_is_not_a_lead(pass_res):
    m = [x for x in pass_res["tiers"]["B"] if x["gene"] == "SYNMASK1"][0]
    assert m["mask_status"] == "inconsistent" and not m["qualifies"]
    assert "SYNMASK1" not in pass_res["verdict"]["lead_genes"]


def test_lipid_positive_control_genes_never_qualify(pass_res):
    for t in "ABCD":
        for x in pass_res["tiers"][t]:
            if x["gene"] in ("PCSK9", "ANGPTL4", "APOC3"):
                assert not x["qualifies"]


def test_ladder_counts(pass_res):
    rg = pass_res["rungs"]
    assert rg["trivial"]["n_genes"] == 0
    assert rg["simplest"]["n_genes"] == 6          # PCSK9 + 5 planted (lipid TG traits are not panel traits)
    assert rg["candidate"]["n_genes_by_tier"] == {"A": 2, "B": 2, "C": 1, "D": 1}
    assert rg["candidate"]["n_genes"] == 3         # A/B with dmis consistent: PCSK9, SYNPASS1, SYNLEAD1


def test_incumbent_not_run_when_file_absent(pass_res):
    assert pass_res["rungs"]["incumbent"]["status"] == "NOT RUN"


def test_incumbent_lookup(make):
    d = make("pass", "incumbent")
    pd.DataFrame([
        {"gene": "SYNPASS1", "trait": "systolic_bp"},
        {"gene": "OLDGENE", "trait": "bmi"},
        {"gene": "OTHER", "trait": "not_in_panel"},
    ]).to_csv(d / "incumbent_hits.csv", index=False)
    inc = run(d)["rungs"]["incumbent"]
    assert inc["status"] == "RUN"
    assert inc["n_gene_trait_pairs"] == 2 and inc["n_out_of_panel_dropped"] == 1
    assert inc["n_recovered_by_discovery"] == 1 and inc["n_recovered_tier_A"] == 1
    assert inc["n_tier_A_not_in_incumbent"] == 1    # PCSK9/ldl


def test_constraint_annotation(make):
    d = make("pass", "constraint")
    pd.DataFrame({"gene": ["synpass1", "PCSK9"], "loeuf": [0.41, float("nan")]}).to_csv(d / "constraint.csv", index=False)
    r = run(d)
    a = {x["gene"]: x["loeuf"] for x in r["tiers"]["A"]}
    assert a["SYNPASS1"] == 0.41 and a["PCSK9"] is None
    assert "0.41" in (d / "out" / "report.md").read_text()


def test_ukb_sourced_replication_rows_are_excluded_and_listed(make):
    def relabel(t):
        rep = t["burden_synthetic_replication"]
        t["burden_synthetic_replication"] = rep.assign(source="Genebass_replicate")
        t["burden_synthetic_discovery"] = t["burden_synthetic_discovery"].assign(source="genebass_disc")   # no synthetic/real mix
        return t
    r = run(make("pass", "ukbrep", relabel))
    assert r["data"]["replication_sources_excluded"] == ["Genebass_replicate"]
    assert "SYNPASS1" not in genes(r, "A") and "SYNPASS1" in genes(r, "B")


def test_lead_scenario(make):
    r = run(make("lead"))
    assert r["verdict"]["verdict"] == "LEAD"
    assert r["verdict"]["lead_genes"] == ["SYNLEAD1"]
    assert r["controls"]["valid"]


def test_nolead_scenario_is_kill_with_valid_controls(make):
    r = run(make("nolead"))
    assert r["verdict"]["verdict"] == "KILL"
    assert r["controls"]["valid"]
    assert "neither PASS nor LEAD" in r["verdict"]["reason"]


@pytest.mark.parametrize("scenario,failed", [
    ("broken_positive", "pcsk9_ldl_lower"),
    ("broken_lambda", "negative_synonymous"),
    ("broken_syn_hit", "negative_synonymous"),
])
def test_broken_controls_kill(make, scenario, failed):
    r = run(make(scenario))
    assert r["verdict"]["verdict"] == "KILL"
    assert not r["controls"]["valid"]
    assert failed in r["controls"]["failed"]
    # a planted Tier-A gene must not rescue a run whose controls failed
    if scenario != "broken_positive":
        assert r["verdict"]["pass_genes"] == ["SYNPASS1"]


def test_broken_syn_hit_reports_contamination_beyond_bound(make):
    n = run(make("broken_syn_hit"))["controls"]["negative_synonymous"]
    assert n["lambda_ok"] and not n["contamination_ok"] and "SYNSYN1" in n["contaminated_genes"]
    assert n["contaminated_fraction"] > n["contaminated_max_fraction"]


def test_missing_synonymous_rows_means_controls_not_evaluable(make):
    def drop_syn(t):
        t["burden_synthetic_discovery"] = t["burden_synthetic_discovery"].query("mask != 'syn'")
        return t
    r = run(make("pass", "nosyn", drop_syn))
    assert r["controls"]["not_run"] == ["negative_synonymous"]
    assert r["verdict"]["verdict"] == "KILL" and "not_evaluable" in r["verdict"]["reason"]
    assert r["rungs"]["trivial"]["status"] == "NOT RUN"


def test_lipid_pathway_gene_cannot_pass(make):
    def rename(t):
        for k in t:
            t[k] = t[k].assign(gene=t[k]["gene"].replace({"SYNPASS1": "LPL"}))
        return t
    r = run(make("pass", "lipidrename", rename))
    assert r["verdict"]["verdict"] == "LEAD"          # LPL is Tier A but excluded; SYNLEAD1 remains
    assert "LPL" in genes(r, "A") and r["verdict"]["pass_genes"] == []


def test_non_proxy_trait_capped_at_tier_b_even_with_same_trait_replication(make):
    def add_rep(t):
        rep = t["burden_synthetic_replication"]
        extra = pd.DataFrame([dict(gene="SYNLEAD1", trait="fluid_intelligence", mask="plof", cohort="replication",
                                   beta=0.9, se=0.05, p=1e-30, n_carriers=100, n_total=300000, source="synthetic_replication")])
        t["burden_synthetic_replication"] = pd.concat([rep, extra], ignore_index=True)
        return t
    r = run(make("lead", "capB", add_rep))
    assert genes(r, "B") >= {"SYNLEAD1"} and "SYNLEAD1" not in genes(r, "A")
    assert r["verdict"]["verdict"] == "LEAD"


def test_only_lipid_hits_is_kill_with_valid_controls(make):
    def drop_planted(t):
        d = t["burden_synthetic_discovery"]
        t["burden_synthetic_discovery"] = d[~d["gene"].isin(["SYNPASS1", "SYNLEAD1", "SYNADV1", "SYNFAIL1", "SYNMASK1"])]
        return t
    r = run(make("pass", "onlylipid", drop_planted))
    assert genes(r, "A") == {"PCSK9"} and not (genes(r, "B") | genes(r, "C") | genes(r, "D"))
    assert r["controls"]["valid"]
    assert r["verdict"]["verdict"] == "KILL" and "neither PASS nor LEAD" in r["verdict"]["reason"]


def test_json_and_report_contents(make):
    d = make("pass", "cli")
    out = d / "results" / "protective-scan.json"
    p = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(CONFIG), "--data", str(d),
                        "--out", str(out)], capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 0, p.stderr
    assert "verdict: PASS" in p.stdout
    j = json.loads(out.read_text())
    assert set(j) >= {"controls", "rungs", "tiers", "verdict", "data", "thresholds"}
    assert set(j["rungs"]) == {"trivial", "simplest", "incumbent", "candidate"}
    assert j["data"]["synthetic"] is True
    md = (d / "results" / "report.md").read_text()
    for needle in ("SYNTHETIC DATA", "## Verdict: PASS", "## Controls", "lambda_GC", "pcsk9_ldl_lower",
                   "## Baseline ladder", "trivial", "simplest", "incumbent", "NOT RUN", "candidate",
                   "### Tier A", "### Tier B", "### Tier C", "### Tier D"):
        assert needle in md, needle
    assert j == json.loads(json.dumps(j))


def test_unscreened_scenario_caps_replicated_gene_at_tier_B(make):
    d = make("unscreened")
    r = run(d)
    b = [x for x in r["tiers"]["B"] if x["gene"] == "SYNPASS1"]
    assert b and b[0]["tradeoff_unscreened"] and b[0]["n_tradeoff_tested"] == 4 and b[0]["rep_status"] == "replicated"
    assert "SYNPASS1" not in genes(r, "A")
    assert not b[0]["qualifies"]                      # C9: listed as Tier B, trade-off unscreened; not a lead
    assert r["verdict"]["verdict"] == "KILL" and r["verdict"]["lead_genes"] == [] and r["verdict"]["pass_genes"] == []
    assert r["controls"]["valid"] and "neither PASS nor LEAD" in r["verdict"]["reason"]
    assert r["thresholds"]["min_tradeoffs_screened"] == 5
    assert "trade-off unscreened" in (Path(d) / "out" / "report.md").read_text()
