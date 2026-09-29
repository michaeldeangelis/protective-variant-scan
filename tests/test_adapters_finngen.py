"""Offline tests for the FinnGen adapter on a synthetic snippet in the real regenie column layout."""
import gzip
from pathlib import Path

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
    assert (raw["TEST"] == "ADD").all() and len(raw) == 7


def test_normalize_endpoint(raw):
    out = finngen.normalize_endpoint(raw, "E4_HYPERCHOL").set_index("gene")
    a = out.loc["GENEA"]
    assert a["beta"] == -0.5 and a["se"] == 0.1 and a["p"] == pytest.approx(10**-5.5)
    assert a["n_total"] == 400000 and a["n_carriers"] == 400 and a["mask"] == "plof"
    assert len(out) == 2  # the DOM-test row is excluded


def test_zero_se_rows_are_dropped(raw):
    assert finngen.normalize_endpoint(raw, "I9_HYPTENS").empty


def test_proxy_rows_pass_the_pipeline_schema(raw):
    frames = [finngen.build_trait(t, raw) for t in ("hypercholesterolemia", "obesity")]
    out = schema.validate(pd.concat(frames, ignore_index=True))
    assert set(out["cohort"]) == {"replication"} and set(out["source"]) == {"finngen_r13"}
    assert set(out["trait"]) == {"hypercholesterolemia", "obesity"} and set(out["mask"]) == {"plof"}


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
