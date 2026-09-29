"""Offline tests for the FinnGen adapter on a synthetic snippet in the real regenie column layout."""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protscan import schema
from protscan.adapters import finngen

FIX = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cfg():
    return schema.load_config(ROOT / "config" / "prereg.yaml")


@pytest.fixture()
def raw(tmp_path):
    gz = tmp_path / "lof.txt.gz"
    gz.write_bytes(gzip.compress((FIX / "finngen_lof_snippet.tsv").read_bytes()))
    return finngen.read_lof(gz, chunksize=3)  # tiny chunks exercise the streaming path


def test_read_lof_filters_endpoints_and_additive_test(raw):
    assert set(raw["PHENO"]) == {"E4_HYPERCHOL", "E4_OBESITY", "I9_HYPTENS", "N14_FEMALEINFERT", "N14_MALEINFERT"}
    assert (raw["TEST"] == "ADD").all() and len(raw) == 10


def test_normalize_endpoint(raw):
    out = finngen.normalize_endpoint(raw, "E4_HYPERCHOL").set_index("gene")
    a = out.loc["GENEA"]
    assert a["beta"] == -0.5 and a["se"] == 0.1 and a["p"] == pytest.approx(10**-5.5)
    assert a["n_total"] == 400000 and a["n_carriers"] == 400 and a["mask"] == "plof"
    assert len(out) == 3  # the DOM-test row is excluded; GENED (se 0) is kept


def test_undefined_se_rows_are_kept_with_nan_se_and_p_intact(raw):
    out = finngen.normalize_endpoint(raw, "I9_HYPTENS")
    assert len(out) == 1 and np.isnan(out.loc[0, "se"]) and out.loc[0, "p"] == 1.0 and out.loc[0, "beta"] == 0.0
    d = finngen.normalize_endpoint(raw, "E4_HYPERCHOL").set_index("gene").loc["GENED"]
    assert np.isnan(d["se"]) and d["p"] == 1.0


def test_a1freq_guard_drops_and_counts_failures(raw, caplog):
    kept, dropped = finngen.guard_a1freq(raw)
    assert (kept["A1FREQ"] <= 0.5).all() and len(kept) == 8
    assert dropped.to_dict() == {"E4_OBESITY": 2}  # A1FREQ 0.7 and A1FREQ missing both fail
    assert "dropped 2 rows" in caplog.text
    ob = finngen.normalize_endpoint(raw, "E4_OBESITY")
    assert list(ob["gene"]) == ["GENEA"]


def test_conversion_summary_counts(raw):
    s = finngen.conversion_summary(raw).set_index("endpoint")
    assert s.loc["E4_OBESITY", "rows_read"] == 3 and s.loc["E4_OBESITY", "dropped_a1freq"] == 2
    assert s.loc["E4_HYPERCHOL", "kept_se_undefined"] == 1 and s.loc["I9_HYPTENS", "kept_se_undefined"] == 1
    assert s["dropped_a1freq"].sum() == 2


def test_kept_rows_never_exceed_half_allele_frequency(raw):
    frames = [finngen.normalize_endpoint(raw, e) for e in ("E4_HYPERCHOL", "E4_OBESITY", "I9_HYPTENS")]
    assert sum(len(f) for f in frames) == 5


def test_proxy_rows_pass_the_pipeline_schema(raw):
    frames = [finngen.build_trait(t, raw) for t in ("hypercholesterolemia", "obesity")]
    out = schema.validate(pd.concat(frames, ignore_index=True))
    assert set(out["cohort"]) == {"replication"} and set(out["source"]) == {"finngen_r13"}
    assert set(out["trait"]) == {"hypercholesterolemia", "obesity"} and set(out["mask"]) == {"plof"}


def test_nan_se_rows_survive_schema_validate(raw):
    df = finngen.build_trait("hypercholesterolemia", raw)
    assert df["se"].isna().sum() == 1
    assert len(schema.validate(df)) == len(df)


def test_infertility_combines_female_and_male(raw):
    df = finngen.build_trait("infertility", raw).set_index("gene")
    assert df.loc["GENEA", "n_total"] == 270000
    w = 1 / 0.2**2 + 1 / 0.3**2
    assert df.loc["GENEA", "beta"] == pytest.approx((0.4 / 0.2**2 + 0.6 / 0.3**2) / w)
    assert df.loc["GENEB", "beta"] == pytest.approx(-0.1)  # male-only gene keeps its own estimate


def test_missing_endpoint_gives_empty_frame(raw):
    assert finngen.build_trait("schizophrenia", raw).empty


def test_endpoints_match_amendment_1_proxies_and_config(cfg):
    proxy_names = {px["trait"] for px in cfg.proxies.values()}
    assert proxy_names == {"hypercholesterolemia", "obesity", "hypertension"}
    assert finngen.ENDPOINTS["hypercholesterolemia"] == ["E4_HYPERCHOL"]
    assert finngen.ENDPOINTS["obesity"] == ["E4_OBESITY"]
    assert finngen.ENDPOINTS["hypertension"] == ["I9_HYPTENS"]
    tradeoff = set(cfg.tradeoff)
    assert set(finngen.ENDPOINTS) - proxy_names <= tradeoff
    assert (tradeoff - set(finngen.ENDPOINTS)) == set(finngen.ABSENT) & tradeoff
    assert set(finngen.ABSENT) - tradeoff <= set(cfg.proxies)  # continuous panel traits with a declared proxy
    assert finngen.MASK in schema.MASKS and finngen.COHORT in schema.COHORTS
