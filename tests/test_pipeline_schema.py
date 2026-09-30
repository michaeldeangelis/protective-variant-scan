import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2

from protscan.schema import COLUMNS, load_burden, load_config, validate

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "prereg.yaml"


@pytest.fixture(scope="session")
def cfg():
    return load_config(CONFIG)


@pytest.fixture(scope="session")
def synth():
    spec = importlib.util.spec_from_file_location("make_synthetic", ROOT / "scripts" / "make_synthetic.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_ledger_constants(cfg):
    assert len(cfg.panel) == 13
    assert len(cfg.tradeoff) == 9
    assert cfg.discovery_p == pytest.approx(1.9e-7)
    assert cfg.replication_p == 0.05
    assert cfg.tradeoff_p == pytest.approx(0.05 / 9)
    assert cfg.lambda_gc_max == 1.10
    assert set(cfg.proxies) == {"ldl", "bmi", "systolic_bp"}


def test_direction_of_benefit(cfg):
    lower_better = {"reaction_time", "pairs_matching_errors", "resting_heart_rate", "systolic_bp", "ldl", "bmi"}
    for t, s in cfg.panel_sign.items():
        assert s == (-1 if t in lower_better else 1), t
    assert all(s == -1 for s in cfg.tradeoff_sign.values())
    assert cfg.control_sign["triglycerides"] == -1


def test_domains(cfg):
    dom = {}
    for t, d in cfg.panel_domain.items():
        dom.setdefault(d, []).append(t)
    assert len(dom["cognitive"]) == 5 and len(dom["physical"]) == 5 and len(dom["metabolic"]) == 3


def _row(**kw):
    r = dict(gene="A", trait="bmi", mask="plof", cohort="discovery", beta=0.1, se=0.1, p=0.3,
             n_carriers=10, n_total=100, source="x")
    r.update(kw)
    return r


def test_validate_ok_and_normalizes():
    df = validate(pd.DataFrame([_row(gene=" pcsk9 ", trait="BMI")]))
    assert df.loc[0, "gene"] == "PCSK9" and df.loc[0, "trait"] == "bmi"


@pytest.mark.parametrize("bad", [
    dict(mask="lof"), dict(cohort="ukb"), dict(p=1.5), dict(se=-1.0),
])
def test_validate_rejects(bad):
    with pytest.raises(ValueError):
        validate(pd.DataFrame([_row(**bad)]))


def test_validate_rejects_missing_column_and_duplicates():
    with pytest.raises(ValueError):
        validate(pd.DataFrame([_row()]).drop(columns=["se"]))
    with pytest.raises(ValueError):
        validate(pd.DataFrame([_row(), _row(source="y")]))


def test_validate_drops_nan_and_counts():
    df = validate(pd.DataFrame([_row(), _row(gene="B", p=np.nan)]))
    assert len(df) == 1 and df.attrs["dropped_nan_rows"] == 1


def test_synthetic_roundtrip(synth, cfg, tmp_path):
    tables = synth.make_synthetic("pass", n_genes=400)
    synth.write(tables, tmp_path)
    df = load_burden(tmp_path)
    assert list(df.columns) == COLUMNS
    assert set(df["cohort"]) == {"discovery", "replication", "discovery_eur"}
    assert set(df["mask"]) == {"plof", "dmis", "syn"}


def test_synthetic_planted_signals(synth, cfg):
    tables = synth.make_synthetic("pass", n_genes=600)
    d = tables["burden_synthetic_discovery"]
    disc = d[d.cohort == "discovery"]
    hits = disc[(disc["mask"] == "plof") & (disc.p < cfg.discovery_p)]
    assert set(hits.gene) == {"PCSK9", "ANGPTL4", "APOC3", "SYNPASS1", "SYNLEAD1", "SYNADV1", "SYNFAIL1", "SYNMASK1"}
    syn = disc[disc["mask"] == "syn"]
    lam = np.median(chi2.isf(syn.p.clip(1e-300), 1)) / chi2.ppf(0.5, 1)
    assert lam < cfg.lambda_gc_max
    assert (syn.p < cfg.discovery_p).sum() == 0


@pytest.mark.parametrize("scenario", ["lead", "nolead", "broken_positive", "broken_lambda", "broken_syn_hit", "contaminated_tierA", "contaminated_systemic", "contaminated_minor"])
def test_synthetic_scenarios_build(synth, scenario):
    tables = synth.make_synthetic(scenario, n_genes=300)
    assert all(list(t.columns) == COLUMNS for t in tables.values())


def test_synthetic_deterministic(synth):
    a = synth.make_synthetic("pass", n_genes=100)["burden_synthetic_discovery"]
    b = synth.make_synthetic("pass", n_genes=100)["burden_synthetic_discovery"]
    pd.testing.assert_frame_equal(a, b)
