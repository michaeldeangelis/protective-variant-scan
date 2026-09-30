"""Offline tests for the shared adapter helpers (synthetic numbers only)."""
import base64
import hashlib
import http.server
import threading

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from protscan.adapters import common


def test_se_from_beta_p_roundtrip():
    p = 2 * norm.sf(5.0)
    assert common.se_from_beta_p(-0.5, p)[()] == pytest.approx(0.1, rel=1e-9)
    assert common.se_from_beta_p([0.2, -0.2], [0.04, 0.04]) == pytest.approx([0.2 / norm.isf(0.02)] * 2)


def test_se_undefined_when_beta_zero_and_p_underflow_is_floored():
    assert np.isnan(common.se_from_beta_p(0.0, 1.0)[()])
    se = common.se_from_beta_p(1.0, 0.0)[()]
    assert np.isfinite(se) and 0 < se < 0.05  # p floored at 1e-300, z ~ 37
    assert common.clip_p([0.0, 0.5, 2.0]).tolist() == [common.P_FLOOR, 0.5, 1.0]


def test_approx_carriers_capped_at_n():
    assert common.approx_carriers(0.0001, 100000)[()] == 20
    assert common.approx_carriers(0.9, 1000)[()] == 1000


def _frame(rows):
    return pd.DataFrame(rows, columns=["gene", "mask", "beta", "se", "n_carriers", "n_total"])


def test_ivw_combine_two_components_and_single_component_gene():
    female = _frame([("A", "plof", 0.4, 0.2, 10, 120000), ("B", "plof", 1.0, 0.5, 5, 120000)])
    male = _frame([("A", "plof", 0.6, 0.3, 7, 150000)])
    out = common.ivw_combine([female, male], "sum").set_index("gene")
    w = 1 / 0.2**2 + 1 / 0.3**2
    assert out.loc["A", "beta"] == pytest.approx((0.4 / 0.2**2 + 0.6 / 0.3**2) / w)
    assert out.loc["A", "se"] == pytest.approx(w**-0.5)
    assert out.loc["A", "p"] == pytest.approx(2 * norm.sf(abs(out.loc["A", "beta"] / out.loc["A", "se"])))
    assert out.loc["A", "n_total"] == 270000 and out.loc["A", "n_carriers"] == 17
    assert out.loc["B", "beta"] == pytest.approx(1.0) and out.loc["B", "se"] == pytest.approx(0.5)


def test_ivw_combine_max_rule_and_bad_rule():
    a = _frame([("A", "plof", 0.1, 0.1, 10, 300)])
    b = _frame([("A", "plof", 0.1, 0.1, 12, 200)])
    out = common.ivw_combine([a, b], "max")
    assert out.loc[0, "n_total"] == 300 and out.loc[0, "n_carriers"] == 12
    with pytest.raises(ValueError):
        common.ivw_combine([a, b], "mean")


def test_finalize_columns_and_types():
    df = _frame([("A", "plof", 0.1, 0.1, 10.0, 300.0)])
    df["p"] = 0.3
    out = common.finalize(df, "bmi", "discovery", "genebass")
    assert list(out.columns) == common.COLUMNS
    assert out.loc[0, "n_carriers"] == 10 and out["n_total"].dtype == np.int64


def test_record_manifest_upserts_by_path(tmp_path):
    m = tmp_path / "MANIFEST.tsv"
    rec = {"url": "u1", "path": "raw/a", "bytes": 1, "sha256": "x", "downloaded_utc": "t1"}
    common.record_manifest(m, rec)
    common.record_manifest(m, {**rec, "url": "u2", "bytes": 2})
    common.record_manifest(m, {**rec, "path": "raw/b"})
    df = pd.read_csv(m, sep="\t")
    assert list(df["path"]) == ["raw/a", "raw/b"] and df.loc[0, "url"] == "u2"


def test_download_streams_and_checksums_over_localhost(tmp_path):
    payload = b"synthetic payload " * 1000
    (tmp_path / "srv").mkdir()
    (tmp_path / "srv" / "f.bin").write_bytes(payload)

    class Quiet(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **k):
            super().__init__(*a, directory=str(tmp_path / "srv"), **k)

        def log_message(self, *a):
            pass

    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Quiet)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        rec = common.download(f"http://127.0.0.1:{srv.server_port}/f.bin", tmp_path / "out" / "f.bin")
    finally:
        srv.shutdown()
    assert (tmp_path / "out" / "f.bin").read_bytes() == payload
    assert not list((tmp_path / "out").glob("*.part"))
    assert rec["bytes"] == len(payload)
    assert rec["sha256"] == hashlib.sha256(payload).hexdigest() == common.sha256_file(tmp_path / "out" / "f.bin")
    assert rec["md5_b64"] == base64.b64encode(hashlib.md5(payload).digest()).decode()


def test_ivw_combine_keeps_genes_whose_se_is_undefined():
    a = _frame([("A", "plof", 0.4, 0.2, 10, 100), ("U", "plof", 0.0, np.nan, 5, 100), ("M", "plof", 0.0, np.nan, 5, 100)])
    a["p"] = [0.04, 1.0, 1.0]
    b = _frame([("A", "plof", 0.6, np.nan, 7, 100), ("U", "plof", 0.0, np.nan, 6, 100), ("M", "plof", 0.5, 0.25, 6, 100)])
    b["p"] = [1.0, 1.0, 0.05]
    out = common.ivw_combine([a, b], "sum").set_index("gene")
    assert list(out.index) == ["A", "M", "U"]
    assert out.loc["A", "beta"] == pytest.approx(0.4) and out.loc["A", "se"] == pytest.approx(0.2)  # NaN-se component ignored
    assert out.loc["M", "beta"] == pytest.approx(0.5) and out.loc["M", "se"] == pytest.approx(0.25)
    assert out.loc["U", "beta"] == 0.0 and np.isnan(out.loc["U", "se"]) and out.loc["U", "p"] == 1.0
    assert out.loc["U", "n_carriers"] == 11


def test_drop_undefined_se_p_lt_1_keeps_p1_rows():
    df = pd.DataFrame({"se": [0.1, np.nan, np.nan, np.nan], "p": [0.5, 1.0, 0.2, 0.9999]})
    kept, n = common.drop_undefined_se_p_lt_1(df)
    assert n == 2 and kept["p"].tolist() == [0.5, 1.0]


def test_ivw_combine_drops_all_undefined_se_genes_with_p_below_1():
    a = _frame([("D", "plof", 0.3, np.nan, 5, 100)])
    a["p"] = [0.01]
    b = _frame([("D", "plof", 0.3, np.nan, 6, 100)])
    b["p"] = [0.02]
    assert common.ivw_combine([a, b], "sum").empty
