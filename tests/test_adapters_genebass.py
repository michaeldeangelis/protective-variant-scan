"""Offline tests for the Genebass adapter on a synthetic snippet in the real API schema."""
import json
from pathlib import Path

import pandas as pd
import pytest
from scipy.stats import norm

from protscan import schema
from protscan.adapters import genebass

FIX = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cfg():
    return schema.load_config(ROOT / "config" / "prereg.yaml")


@pytest.fixture()
def records():
    return json.loads((FIX / "genebass_gene_manhattan_snippet.json").read_text())


@pytest.fixture()
def qc():
    return json.loads((FIX / "genebass_gene_qc_snippet.json").read_text())


def test_qc_pass_uses_coverage_and_nvar_only(qc):
    kept = set(genebass.qc_pass(qc)["gene_id"])
    assert "ENSG00000000002" in kept  # keep_gene_burden=false (lambda-based) is deliberately ignored
    assert "ENSG00000000007" not in kept  # coverage fail
    assert "ENSG00000000009" not in kept  # null n_var flag counts as fail


def test_normalize_analysis_drops_and_derives(records, qc):
    out = genebass.normalize_analysis(records, qc, "plof", 100000).set_index("gene")
    # dropped: null result (GENEC), ambiguous symbol (DUPSYM), beta=0 (GENED), coverage fail (GENEE), null n_var (GENEG)
    assert sorted(out.index) == ["GENEA", "GENEB", "GENEF"]
    a = out.loc["GENEA"]
    assert a["beta"] == -0.5 and a["se"] == pytest.approx(0.1, rel=1e-6)
    assert a["p"] == pytest.approx(2 * norm.sf(5.0))
    assert a["n_carriers"] == 20 and a["n_total"] == 100000 and a["mask"] == "plof"
    assert out.loc["GENEB", "se"] == pytest.approx(0.2 / norm.isf(0.02))
    f = out.loc["GENEF"]
    assert f["p"] == 1e-300 and 0 < f["se"] < 0.05


def test_normalized_rows_pass_the_pipeline_schema(records, qc):
    per = [genebass.normalize_analysis(records, qc, m, 1000) for m in ("plof", "dmis", "syn")]
    df = genebass.build_trait("ldl", [pd.concat(per, ignore_index=True)])
    out = schema.validate(df)
    assert set(out["cohort"]) == {"discovery"} and set(out["source"]) == {"genebass"}
    assert set(out["mask"]) == {"plof", "dmis", "syn"} and len(out) == 9


def test_two_analysis_trait_is_ivw_combined(records, qc):
    father = genebass.normalize_analysis(records, qc, "plof", 291196)
    mother = genebass.normalize_analysis(records, qc, "plof", 234028)
    df = genebass.build_trait("parental_lifespan", [father, mother]).set_index("gene")
    assert df.loc["GENEA", "se"] == pytest.approx(0.1 / 2**0.5, rel=1e-6)
    assert df.loc["GENEA", "n_total"] == 291196  # same participants: max, not sum
    inf = genebass.build_trait("infertility", [father, mother]).set_index("gene")
    assert inf.loc["GENEA", "n_total"] == 291196 + 234028  # disjoint sexes: sum


def test_analysis_n_total():
    assert genebass.analysis_n_total({"n_cases": 100, "n_controls": None}) == 100
    assert genebass.analysis_n_total({"n_cases": 100, "n_controls": 900}) == 1000


def test_trait_coverage_matches_config(cfg):
    declared = set(cfg.panel) | set(cfg.tradeoff) | set(cfg.control_sign)
    assert set(genebass.TRAITS) | set(genebass.ABSENT) == declared
    assert not set(genebass.TRAITS) & set(genebass.ABSENT)
    assert set(genebass.BURDEN_SET_TO_MASK.values()) == set(schema.MASKS)
    assert set(genebass.N_RULE) == {t for t, ids in genebass.TRAITS.items() if len(ids) > 1}
    assert len({a for ids in genebass.TRAITS.values() for a in ids}) == sum(len(v) for v in genebass.TRAITS.values())
