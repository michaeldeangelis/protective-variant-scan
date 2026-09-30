"""Amendment 1 clarifications C11a-C11d (review-2 fixes)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from protscan import controls, stats, tiering
from protscan.run import decide, run_pipeline
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


@pytest.fixture(scope="module")
def synth():
    spec = importlib.util.spec_from_file_location("make_synthetic_c11", ROOT / "scripts" / "make_synthetic.py")
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


def run(d):
    return run_pipeline(CONFIG, d, Path(d) / "out" / "protective-scan.json")


def edit_rows(key, gene, trait, mask, cohort, **vals):
    def _edit(t):
        d = t[key]
        i = d.index[(d.gene == gene) & (d.trait == trait) & (d["mask"] == mask) & (d.cohort == cohort)]
        assert len(i) == 1
        for k, v in vals.items():
            d.loc[i[0], k] = v
        t[key] = d
        return t
    return _edit


# ---------------- C11a ----------------
def _sbp_replicated(beta, se, p=1e-9):
    rows = [row("G", "systolic_bp"), row("G", "systolic_bp", "dmis", beta=-0.1, p=0.3),
            row("G", "hypertension", cohort="replication", beta=beta, se=se, p=p)]
    rows += [row("G", t, beta=0.01, p=0.9) for t in cfg.tradeoff]
    return table(rows)


@pytest.mark.parametrize("beta,se,tier", [
    (-0.3, 0.05, "A"),           # informative
    (-0.3, np.nan, "B"),         # undefined se with a significant p: degenerate fit, not evidence
    (-0.3, 0.0, "B"),            # se = 0
    (0.0, 0.05, "B"),            # beta = 0
    (0.0, np.nan, "B"),
    (-0.3, 1e-12, "A"),          # tiny but positive se is still defined
])
def test_replication_rule_needs_finite_positive_se_and_nonzero_beta(beta, se, tier):
    h = tiering.build_hits(_sbp_replicated(beta, se), cfg).iloc[0]
    assert h["tier"] == tier
    if tier == "B":
        assert h["rep_status"] in ("gene_untested", "trait_missing")    # an unusable row is not a failed replication (no Tier D)


def test_informative_filter():
    df = table([row("A", "obesity", cohort="replication", beta=-0.1, se=0.1),
                row("B", "obesity", cohort="replication", beta=-0.1, se=np.nan),
                row("C", "obesity", cohort="replication", beta=-0.1, se=0.0),
                row("D", "obesity", cohort="replication", beta=0.0, se=0.1)])
    assert list(stats.informative(df).gene) == ["A"]


def test_degenerate_replication_row_cannot_give_tier_A_end_to_end(make):
    # reviewer probe: the failed replication of SYNFAIL1 replaced by a significant row with undefined se
    d = make("pass", "degenerate", edit_rows("burden_synthetic_replication", "SYNFAIL1", "hypertension", "plof", "replication",
                                             beta=-0.3, p=1e-9, se=np.nan))
    r = run(d)
    assert "SYNFAIL1" not in {x["gene"] for x in r["tiers"]["A"]}
    assert {x["gene"] for x in r["tiers"]["B"]} >= {"SYNFAIL1"}
    assert r["verdict"]["pass_genes"] == ["SYNPASS1"]


def test_undefined_se_rows_still_serve_lambda_and_coverage():
    rng = np.random.default_rng(5)
    u = rng.uniform(size=20_000)
    rows = [row(f"S{i}", "fev1", "plof", beta=0.0, p=float(u[i]), se=np.nan) for i in range(10_000)]
    rows += [row(f"S{i}", "fev1", "syn", beta=0.0, p=float(u[10_000 + i]), se=np.nan) for i in range(10_000)]
    c = controls.negative_control(table(rows), cfg)
    assert c["status"] == "OK" and c["coverage"] == 1.0 and c["lambda_gc"] == pytest.approx(1.0, abs=0.05)


def test_sign_arm_ignores_undefined_se_rows():
    df = table([row("PCSK9", "hypercholesterolemia", cohort="replication", beta=0.3, se=np.nan, p=1e-9),
                row("Q", "fev1")])
    rs = controls.replication_sign_control(df, cfg)
    assert rs["status"] == "NOT RUN" and rs["checks"][0]["why"] == "no informative row"


# ---------------- C11b ----------------
def _syn_universe(n, n_p1, seed=11, scale=1.0):
    """n plof + n syn rows; the first n_p1 syn rows have p = 1 (beta 0); the rest are N(0, scale^2) z-scores."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(2 * n) * scale
    p = 2 * norm.sf(np.abs(z))
    rows = [row(f"S{i}", "fev1", "plof", beta=0.0, p=float(p[i]), se=np.nan) for i in range(n)]
    for i in range(n):
        if i < n_p1:
            rows.append(row(f"S{i}", "fev1", "syn", beta=0.0, p=1.0, se=np.nan))
        else:
            rows.append(row(f"S{i}", "fev1", "syn", beta=0.01, p=float(p[n + i]), se=0.05))
    return rows


@pytest.mark.parametrize("n_p1,status", [(0, "OK"), (500, "OK"), (501, "NOT RUN")])
def test_p1_fraction_boundary_is_five_percent(n_p1, status):
    c = controls.negative_control(table(_syn_universe(10_000, n_p1)), cfg)
    assert c["status"] == status
    assert c["n_p1_rows"] == n_p1 and c["p1_fraction"] == pytest.approx(n_p1 / 10_000)
    if status == "NOT RUN":
        assert "p = 1" in c["reason"]


def test_lambda_gate_uses_rows_with_p_below_one_and_both_lambdas_are_reported():
    # informative rows inflated (lambda ~1.15) plus 4.9 percent p = 1 rows: the pooled lambda dips below 1.10 (masking),
    # the p < 1 lambda does not
    rows = _syn_universe(20_000, 980, scale=1.15 ** 0.5)
    c = controls.negative_control(table(rows), cfg)
    assert c["lambda_gc_all_rows"] < cfg.lambda_gc_max <= c["lambda_gc"]
    assert c["status"] == "FAIL" and not c["lambda_ok"]
    assert c["lambda_gc"] == stats.lambda_gc(np.array([r["p"] for r in rows if r["mask"] == "syn" and r["p"] < 1]))


def test_p1_rows_never_count_as_syn_hits_and_lambda_unchanged_without_p1():
    c = controls.negative_control(table(_syn_universe(10_000, 0)), cfg)
    assert c["lambda_gc"] == pytest.approx(c["lambda_gc_all_rows"]) and c["n_syn_hit_genes"] == 0


def test_masked_inflation_end_to_end_is_not_evaluable_or_failed(make, synth):
    def mask_it(t):
        d = t["burden_synthetic_discovery"].copy()
        s = (d["mask"] == "syn") & (d.cohort == "discovery")
        idx = d.index[s]
        rng = np.random.default_rng(1)
        p1 = rng.choice(idx, size=int(0.45 * len(idx)), replace=False)      # 45 percent uninformative rows
        rest = idx.difference(p1)
        z = d.loc[rest, "beta"] / d.loc[rest, "se"] * 1.4                   # informative rows inflated (lambda ~2)
        d.loc[rest, "beta"] = z * d.loc[rest, "se"]
        d.loc[rest, "p"] = 2 * norm.sf(np.abs(z))
        d.loc[p1, ["beta", "p"]] = [0.0, 1.0]
        d.loc[p1, "se"] = np.nan
        t["burden_synthetic_discovery"] = d
        return t
    r = run(make("pass", "masked", mask_it))
    n = r["controls"]["negative_synonymous"]
    assert n["status"] == "NOT RUN" and n["p1_fraction"] > 0.05
    assert r["verdict"]["verdict"] == "KILL" and r["verdict"]["reason"].startswith("controls_not_evaluable")


def test_report_prints_both_lambdas(make):
    d = make("pass")
    run(d)
    md = (d / "out" / "report.md").read_text()
    assert "on rows with p < 1" in md and "on all" in md and "p = 1 rows:" in md


# ---------------- C11c ----------------
def _rep(pcsk9=None, ldlr=None):
    rows = [row("Q", "fev1")]
    for gene, spec in (("PCSK9", pcsk9), ("LDLR", ldlr)):
        if spec is not None:
            beta, p = spec
            rows.append(row(gene, "hypercholesterolemia", cohort="replication", beta=beta, p=p))
    return table(rows)


@pytest.mark.parametrize("p,status", [(0.0499, "FAIL"), (0.05, "NOT RUN"), (0.2, "NOT RUN")])
def test_arm_evaluable_only_below_p_005_wrong_sign(p, status):
    rs = controls.replication_sign_control(_rep(pcsk9=(+0.2, p)), cfg)
    assert rs["status"] == status


@pytest.mark.parametrize("p,status", [(0.0499, "OK"), (0.05, "NOT RUN")])
def test_arm_evaluable_only_below_p_005_right_sign(p, status):
    assert controls.replication_sign_control(_rep(ldlr=(+0.2, p)), cfg)["status"] == status


@pytest.mark.parametrize("pcsk9,ldlr,status,reason", [
    ((+0.2, 0.3), (+0.5, 1e-6), "OK", None),                       # underpowered wrong-sign PCSK9 arm is not a failure
    ((-0.2, 0.3), (+0.5, 1e-6), "OK", None),
    ((+0.2, 0.3), (-0.5, 0.4), "NOT RUN", "sign_control_not_evaluable"),
    (None, (-0.5, 0.4), "NOT RUN", "sign_control_not_evaluable"),
    (None, None, "NOT RUN", "sign_control_not_evaluable"),
    ((+0.2, 0.01), (+0.5, 1e-6), "FAIL", "sign_control_failed"),   # significant wrong sign
    ((-0.2, 0.01), (-0.5, 1e-6), "FAIL", "sign_control_failed"),
    ((+0.2, 0.3), (-0.5, 1e-6), "FAIL", "sign_control_failed"),    # an evaluable arm with the wrong sign fails the control
])
def test_sign_control_arms_and_reasons(pcsk9, ldlr, status, reason):
    rs = controls.replication_sign_control(_rep(pcsk9, ldlr), cfg)
    assert rs["status"] == status and rs["reason"] == reason


def test_decide_records_sign_control_reasons():
    hits = pd.DataFrame(columns=["gene", "trait", "tier", "qualifies"])
    failed = {"valid": False, "failed": ["replication_sign"], "not_run": [], "replication_sign": {"reason": "sign_control_failed"}}
    nrun = {"valid": False, "failed": [], "not_run": ["replication_sign"], "replication_sign": {"reason": "sign_control_not_evaluable"}}
    vf, vn = decide(failed, hits), decide(nrun, hits)
    assert vf["verdict"] == vn["verdict"] == "KILL"
    assert vf["reason"].startswith("control_failed") and "sign_control_failed" in vf["reason"]
    assert vn["reason"].startswith("controls_not_evaluable") and "sign_control_not_evaluable" in vn["reason"]


def test_sign_control_end_to_end(make):
    r = run(make("broken_rep_sign"))
    assert r["verdict"]["verdict"] == "KILL" and "sign_control_failed" in r["verdict"]["reason"]
    # PCSK9 arm underpowered with the wrong sign, LDLR arm intact: control OK, verdict unaffected
    under = edit_rows("burden_synthetic_replication", "PCSK9", "hypercholesterolemia", "plof", "replication", beta=0.1, p=0.3)
    r2 = run(make("pass", "underpowered", under))
    assert r2["controls"]["replication_sign"]["status"] == "OK" and r2["verdict"]["verdict"] == "PASS"
    checks = {c["id"]: c for c in r2["controls"]["replication_sign"]["checks"]}
    assert checks["pcsk9_rep_sign"]["status"] == "NOT RUN" and checks["ldlr_rep_sign"]["status"] == "OK"

    # both arms underpowered: not evaluable
    def both(t):
        return edit_rows("burden_synthetic_replication", "LDLR", "hypercholesterolemia", "plof", "replication",
                         beta=0.1, p=0.3)(under(t))
    r3 = run(make("pass", "bothunder", both))
    assert r3["controls"]["not_run"] == ["replication_sign"]
    assert r3["verdict"]["verdict"] == "KILL" and "sign_control_not_evaluable" in r3["verdict"]["reason"]


def test_report_reading_rule_only_on_control_kills(make):
    d1 = make("broken_rep_sign", "rr1")
    run(d1)
    assert "Reading rule" in (d1 / "out" / "report.md").read_text()
    d2 = make("pass", "rr2")
    run(d2)
    assert "Reading rule" not in (d2 / "out" / "report.md").read_text()
    d3 = make("nolead", "rr3")                       # KILL with valid controls: biology conclusion, no reading rule
    run(d3)
    assert "Reading rule" not in (d3 / "out" / "report.md").read_text()


# ---------------- C11d ----------------
def test_pcsk9_t2d_tradeoff_never_changes_control_status_end_to_end(make):
    adverse = edit_rows("burden_synthetic_discovery", "PCSK9", "type_2_diabetes", "plof", "discovery", beta=1.0, p=1e-6)
    r = run(make("pass", "pcsk9t2d", adverse))
    assert r["controls"]["valid"] and r["verdict"]["verdict"] == "PASS"
    pos = {c["id"]: c for c in r["controls"]["positive"]}
    assert pos["pcsk9_ldl_lower"]["status"] == "OK" and pos["pcsk9_ldl_lower"]["tested"][0]["pipeline_tier"] == "C"
    assert "PCSK9" in {x["gene"] for x in r["tiers"]["C"]}


def test_positive_control_still_needs_discovery_effect():
    weak = table([row("PCSK9", "ldl", beta=-0.8, p=1e-5), row("PCSK9", "coronary_disease", beta=-0.1, p=0.3),
                  row("APOC3", "triglycerides", beta=-1.0, p=1e-20)])
    assert [c["status"] for c in controls.positive_controls(weak, cfg)] == ["FAIL", "OK", "OK"]
    wrong = table([row("PCSK9", "ldl", beta=+0.8, p=1e-30)])
    assert controls.positive_controls(wrong, cfg)[0]["status"] == "FAIL"
