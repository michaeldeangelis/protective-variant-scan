"""Amendment 1 clarifications C10a-C10g (review-1 fixes)."""
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protscan import controls, stats, tiering
from protscan.run import decide, run_pipeline
from protscan.schema import PINNED_CONFIG_SHA256, load_burden, load_config, validate

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "prereg.yaml"
cfg = load_config(CONFIG)
VALID = {"valid": True, "failed": [], "not_run": []}


def row(gene, trait, mask="plof", cohort="discovery", beta=-0.5, p=1e-9, source=None):
    src = source or ("finngen_test" if cohort == "replication" else "genebass_test")
    return dict(gene=gene, trait=trait, mask=mask, cohort=cohort, beta=beta, se=0.05, p=p,
                n_carriers=100, n_total=400000, source=src)


def table(rows):
    return validate(pd.DataFrame(rows))


def universe(n_pairs, syn_pairs=None, seed=3):
    """Null plof rows for n_pairs (gene, fev1) pairs; syn rows for the first syn_pairs of them."""
    syn_pairs = n_pairs if syn_pairs is None else syn_pairs
    rng = np.random.default_rng(seed)
    u = rng.uniform(size=2 * n_pairs)
    rows = [row(f"S{i}", "fev1", "plof", beta=0.0, p=float(u[i])) for i in range(n_pairs)]
    rows += [row(f"S{i}", "fev1", "syn", beta=0.0, p=float(u[n_pairs + i])) for i in range(syn_pairs)]
    return rows


@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("make_synthetic_c10", ROOT / "scripts" / "make_synthetic.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def make(synth, tmp_path_factory):
    def _make(scenario, name=None, edit=None):
        d = tmp_path_factory.mktemp(name or scenario)
        t = synth.make_synthetic(scenario, n_genes=500)
        if edit:
            t = edit(t)
        synth.write(t, d)
        return d
    return _make


def run(d, config=CONFIG):
    return run_pipeline(config, d, Path(d) / "out" / "protective-scan.json")


# ---------------- C10a ----------------
@pytest.mark.parametrize("n,status", [(9_999, "NOT RUN"), (10_000, "OK")])
def test_syn_min_rows_boundary(n, status):
    c = controls.negative_control(table(universe(n)), cfg)
    assert c["status"] == status and c["n_rows"] == n


@pytest.mark.parametrize("syn_pairs,status", [(8_999, "NOT RUN"), (9_000, "OK")])
def test_syn_min_coverage_boundary(syn_pairs, status):
    # 10,000 plof pairs; syn rows padded to >= 10,000 with pairs that have no plof row so only coverage varies
    rows = universe(10_000, syn_pairs)
    rows += [row(f"X{i}", "fev1", "syn", beta=0.0, p=0.5) for i in range(10_000 - syn_pairs)]
    c = controls.negative_control(table(rows), cfg)
    assert c["status"] == status
    assert c["coverage"] == pytest.approx(syn_pairs / 10_000)
    assert c["n_rows"] >= 10_000


def test_syn_thin_universe_is_not_evaluable_and_kills():
    rows = universe(20, 3) + [row("PCSK9", "ldl", "plof", beta=-0.8, p=1e-30),
                              row("PCSK9", "coronary_disease", beta=-0.1, p=0.3),
                              row("APOC3", "triglycerides", beta=-1.0, p=1e-20)]
    c = controls.run_controls(table(rows), cfg)
    assert "negative_synonymous" in c["not_run"] and not c["valid"]
    assert decide(c, tiering.build_hits(table(rows), cfg))["reason"].startswith("controls_not_evaluable")


def test_report_shows_syn_coverage(make):
    d = make("pass")
    run(d)
    assert re.search(r"coverage = 1\.000 of \d+ pairs", (d / "out" / "report.md").read_text())


# ---------------- C10b ----------------
def sign_rows(pcsk9=None, ldlr=None, source=None, p=0.01):
    rows = []
    if pcsk9 is not None:
        rows.append(row("PCSK9", "hypercholesterolemia", cohort="replication", beta=pcsk9, p=p, source=source))
    if ldlr is not None:
        rows.append(row("LDLR", "hypercholesterolemia", cohort="replication", beta=ldlr, p=p, source=source))
    return rows


@pytest.mark.parametrize("pcsk9,ldlr,status", [
    (-0.1, 0.1, "OK"), (-0.1, None, "OK"), (None, 0.1, "OK"),
    (0.1, 0.1, "FAIL"), (-0.1, -0.1, "FAIL"), (0.1, None, "FAIL"), (None, -0.1, "FAIL"),
    (None, None, "NOT RUN"),
])
def test_replication_sign_control_rules(pcsk9, ldlr, status):
    df = table(sign_rows(pcsk9, ldlr) + [row("Q", "fev1")])
    assert controls.replication_sign_control(df, cfg)["status"] == status


def test_replication_sign_ignores_non_independent_sources():
    df = table(sign_rows(-0.1, 0.1, source="genebass_x") + [row("Q", "fev1")])
    assert controls.replication_sign_control(df, cfg)["status"] == "NOT RUN"


def test_replication_sign_failure_kills_and_not_run_is_not_evaluable(make):
    r = run(make("broken_rep_sign"))
    assert r["verdict"]["verdict"] == "KILL" and r["controls"]["failed"] == ["replication_sign"]
    assert r["controls"]["replication_sign"]["status"] == "FAIL"

    def drop(t):
        rep = t["burden_synthetic_replication"]
        t["burden_synthetic_replication"] = rep[~rep.gene.isin(["PCSK9", "LDLR"])]
        return t
    r2 = run(make("pass", "norepsign", drop))
    assert r2["controls"]["not_run"] == ["replication_sign"]
    assert r2["verdict"]["verdict"] == "KILL" and "not_evaluable" in r2["verdict"]["reason"]


def test_replication_sign_pass_scenario_ok(make):
    r = run(make("pass"))
    rs = r["controls"]["replication_sign"]
    assert rs["status"] == "OK" and [c["status"] for c in rs["checks"]] == ["OK", "OK"]


# ---------------- C10c ----------------
def test_pinned_hash_matches_committed_config():
    assert PINNED_CONFIG_SHA256 == hashlib.sha256(CONFIG.read_bytes()).hexdigest()


def test_committed_config_is_preregistered(make):
    r = run(make("pass"))
    assert r["config_matches_ledger"] is True and r["verdict"]["config_matches_ledger"] is True
    assert "NON-PREREGISTERED" not in (Path(r["data"]["dir"]) / "out" / "report.md").read_text()


def test_edited_config_is_stamped_non_preregistered_and_cli_exits_nonzero(make, tmp_path):
    loose = tmp_path / "loose.yaml"
    loose.write_text(CONFIG.read_text().replace("syn_contaminated_max_fraction: 0.001", "syn_contaminated_max_fraction: 0.5"))
    assert loose.read_text() != CONFIG.read_text()
    d = make("broken_syn_hit", "loose")
    r = run(d, loose)
    assert r["config_matches_ledger"] is False and r["verdict"]["config_matches_ledger"] is False
    assert "[NON-PREREGISTERED]" in r["verdict"]["label"]
    assert "NON-PREREGISTERED" in (d / "out" / "report.md").read_text()
    p = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(loose), "--data", str(d),
                        "--out", str(tmp_path / "o" / "r.json")], capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 2 and "NON-PREREGISTERED" in p.stdout + p.stderr
    ok = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(CONFIG), "--data", str(d),
                         "--out", str(tmp_path / "o2" / "r.json")], capture_output=True, text=True, cwd=ROOT)
    assert ok.returncode == 0


# ---------------- C10d ----------------
def test_undefined_se_rows_are_kept_for_lambda_and_hit_counting():
    rows = universe(10_000)
    rows.append(dict(row("Z1", "fev1", "syn", beta=0.0, p=1.0), se=np.nan))
    rows.append(dict(row("Z2", "fev1", "syn", beta=0.5, p=1e-9), se=np.nan))
    df = table(rows)
    assert df.attrs["dropped_nan_rows"] == 0 and df["se"].isna().sum() == 2
    c = controls.negative_control(df, cfg)
    assert c["n_rows"] == 10_002 and c["status"] == "OK" and c["syn_hit_genes"] == ["Z2"] and c["n_contaminated"] == 1
    assert table([dict(row("A", "fev1"), beta=np.nan)]).empty       # missing beta or p is still dropped


# ---------------- C10e ----------------
def test_harmful_panel_hits_are_reported_but_do_not_gate(make):
    def add(t):
        d = t["burden_synthetic_discovery"]
        i = d.index[(d.gene == "SYNPASS1") & (d.trait == "fev1") & (d["mask"] == "plof") & (d.cohort == "discovery")][0]
        d.loc[i, ["beta", "p"]] = [-0.9, 1e-12]                  # fev1 lower = harmful
        t["burden_synthetic_discovery"] = d
        return t
    d = make("pass", "harm", add)
    r = run(d)
    a = [x for x in r["tiers"]["A"] if x["gene"] == "SYNPASS1"][0]
    assert [h["trait"] for h in a["harmful_panel"]] == ["fev1"]
    assert r["verdict"]["verdict"] == "PASS"                     # informational only
    assert "fev1 (p=1e-12)" in (d / "out" / "report.md").read_text()


def test_harmful_panel_hits_direction_and_threshold():
    df = table([row("G", "systolic_bp"), row("G", "fev1", beta=-0.5, p=1e-9),           # harmful, significant
                row("G", "hand_grip_strength", beta=-0.5, p=1e-5),                       # harmful, not significant
                row("G", "walking_pace", beta=+0.5, p=1e-9),                             # beneficial: also a hit, not harmful
                row("G", "resting_heart_rate", beta=+0.5, p=1e-9)])                      # harmful (higher HR), significant
    hp = stats.harmful_panel_hits(["G"], df, cfg)["G"]
    assert sorted(h["trait"] for h in hp) == ["fev1", "resting_heart_rate"]
    h = tiering.build_hits(df, cfg)
    sbp = h[h.trait == "systolic_bp"].iloc[0]
    assert sorted(x["trait"] for x in sbp["harmful_panel"]) == ["fev1", "resting_heart_rate"]


# ---------------- C10f ----------------
def test_caveats_printed_for_pass_and_lead_only(make):
    for scen, verdict in (("pass", "PASS"), ("lead", "LEAD"), ("nolead", "KILL")):
        d = make(scen, "cav_" + scen)
        r = run(d)
        md = (d / "out" / "report.md").read_text()
        assert r["verdict"]["verdict"] == verdict
        assert ("missense|LC" in md and "fluid intelligence and reaction time" in md) == (verdict != "KILL"), scen
        assert ("C5" in md and "C7" in md) == (verdict != "KILL"), scen


# ---------------- C10g ----------------
def _sbp_gene(k_disc, rep_source, k_rep=0):
    outs = list(cfg.tradeoff)
    rows = [row("G", "systolic_bp"), row("G", "systolic_bp", "dmis", beta=-0.1, p=0.3),
            row("G", "hypertension", cohort="replication", beta=-0.3, p=0.01)]
    rows += [row("G", t, beta=0.01, p=0.9) for t in outs[:k_disc]]
    rows += [row("G", t, cohort="replication", beta=0.01, p=0.9, source=rep_source) for t in outs[k_disc:k_disc + k_rep]]
    return table(rows)


@pytest.mark.parametrize("src,counted", [("finngen_r13", True), ("mystery_biobank", False), ("genebass_x", False)])
def test_screening_counts_only_allow_listed_replication_rows(src, counted):
    h = tiering.build_hits(_sbp_gene(3, src, k_rep=3), cfg).iloc[0]
    assert h["n_tradeoff_tested"] == (6 if counted else 3)
    assert h["tier"] == ("A" if counted else "B")


def test_untrusted_replication_rows_still_flag_adverse():
    df = table([row("G", "systolic_bp"), row("G", "coronary_disease", cohort="replication", beta=0.9, p=1e-6,
                                             source="mystery_biobank")])
    assert stats.adverse_tradeoffs(["G"], df, cfg)["G"]         # conservative: harm evidence is never ignored


def test_requests_and_pyarrow_declared():
    txt = (ROOT / "pyproject.toml").read_text()
    deps = re.search(r"dependencies = \[(.*?)\]", txt, re.S).group(1)
    assert '"requests"' in deps and '"pyarrow"' in deps


def test_positive_control_uses_pipeline_discovery_path_and_ignores_tradeoffs():
    good = [row("PCSK9", "ldl", beta=-0.8, p=1e-30), row("PCSK9", "coronary_disease", beta=-0.1, p=0.3),
            row("APOC3", "triglycerides", beta=-1.0, p=1e-20)]
    ok = controls.positive_controls(table(good), cfg)
    assert [c["status"] for c in ok] == ["OK", "OK", "OK"]
    assert ok[0]["via_pipeline_discovery_path"] and ok[0]["tested"][0]["pipeline_tier"] == "B"
    assert not ok[1]["via_pipeline_discovery_path"] and not ok[2]["via_pipeline_discovery_path"]
    # C11d: a significant adverse trade-off (PCSK9 and type 2 diabetes) makes PCSK9 Tier C but never changes control status
    adverse = good + [row("PCSK9", "type_2_diabetes", beta=0.9, p=1e-6)]
    c = controls.positive_controls(table(adverse), cfg)
    assert [x["status"] for x in c] == ["OK", "OK", "OK"] and c[0]["tested"][0]["pipeline_tier"] == "C"


def test_positive_control_not_in_pipeline_hits_fails(monkeypatch):
    df = table([row("PCSK9", "ldl", beta=-0.8, p=1e-30)])
    monkeypatch.setattr(stats, "discovery_hits", lambda d, c, mask="plof": pd.DataFrame(columns=["gene", "trait"]))
    assert controls.positive_controls(df, cfg)[0]["status"] == "FAIL"


# ---------------- synthetic flag and mixed inputs ----------------
def test_synthetic_flag_in_verdict_block_cli_output_and_label(make, tmp_path):
    d = make("pass")
    r = run(d)
    assert r["verdict"]["synthetic"] is True and "[SYNTHETIC DATA]" in r["verdict"]["label"]
    p = subprocess.run([sys.executable, "-m", "protscan", "run", "--config", str(CONFIG), "--data", str(d),
                        "--out", str(tmp_path / "o" / "r.json")], capture_output=True, text=True, cwd=ROOT)
    assert p.returncode == 0 and "PASS [SYNTHETIC DATA]" in p.stdout
    j = json.loads((tmp_path / "o" / "r.json").read_text())
    assert j["verdict"]["synthetic"] is True


def test_real_source_names_are_not_labelled_synthetic(make):
    def relabel(t):
        return {k: v.assign(source=v["source"].str.replace("synthetic_discovery", "genebass_x").str.replace(
            "synthetic_replication", "finngen_x")) for k, v in t.items()}
    r = run(make("pass", "realnames", relabel))
    assert r["verdict"]["synthetic"] is False and "SYNTHETIC" not in r["verdict"]["label"]


def test_mixed_synthetic_and_real_files_are_refused(make):
    def real_replication(t):
        t["burden_finngen"] = t.pop("burden_synthetic_replication").assign(source="finngen_r13")
        return t
    d = make("pass", "mixed", real_replication)
    with pytest.raises(ValueError, match="mixes synthetic"):
        load_burden(d)
    with pytest.raises(ValueError, match="mixes synthetic"):
        run(d)
