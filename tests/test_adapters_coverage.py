"""Offline test: the coverage table records absences with reasons and never substitutes."""
from pathlib import Path

import pandas as pd
import pytest

from protscan import schema
from protscan.adapters import coverage

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cfg():
    return schema.load_config(ROOT / "config" / "prereg.yaml")


def _rows(trait, mask, cohort, n):
    return pd.DataFrame({"trait": trait, "mask": mask, "cohort": cohort, "gene": [f"G{i}" for i in range(n)]})


def test_build_coverage(cfg):
    disc = pd.concat([_rows("ldl", "plof", "discovery", 3), _rows("ldl", "syn", "discovery", 2)])
    rep = _rows("hypercholesterolemia", "plof", "replication", 4)
    cov = coverage.build_coverage(cfg, disc, rep).set_index("trait")
    assert len(cov) == 13 + 9 + 1
    assert list(cov.columns) == coverage.COLUMNS[1:]
    ldl = cov.loc["ldl"]
    assert ldl["discovery_genes_plof"] == 3 and ldl["discovery_genes_syn"] == 2 and ldl["discovery_genes_dmis"] == 0
    assert ldl["replication_trait"] == "hypercholesterolemia" and ldl["replication_genes_plof"] == 4
    assert "Tier A possible" in ldl["replication_role"] and ldl["replication_source"] == "E4_HYPERCHOL"
    assert cov.loc["walking_pace", "discovery_source"].startswith("ABSENT")
    assert cov.loc["fev1", "replication_source"].startswith("ABSENT") and "Tier B" in cov.loc["fev1", "replication_role"]
    assert cov.loc["fracture", "replication_source"].startswith("ABSENT")
    assert cov.loc["type_2_diabetes", "replication_role"] == "supplementary trade-off rows"
    assert cov.loc["triglycerides", "replication_role"] == "none"
