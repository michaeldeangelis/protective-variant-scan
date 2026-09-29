"""Offline test for the gnomAD constraint adapter on a synthetic snippet in the real column layout."""
from pathlib import Path

import pandas as pd
import pytest

from protscan.adapters import gnomad

FIX = Path(__file__).parent / "fixtures"


def test_normalize_constraint_picks_one_transcript_per_gene():
    df = pd.read_csv(FIX / "gnomad_constraint_snippet.tsv", sep="\t", na_values=["NA"])
    out = gnomad.normalize_constraint(df).set_index("gene")
    assert list(gnomad.normalize_constraint(df).columns) == ["gene", "loeuf", "pli"]
    assert sorted(out.index) == ["GENEA", "GENEB", "GENEE"]
    assert out.loc["GENEA", "loeuf"] == pytest.approx(0.31)  # MANE Ensembl row, not RefSeq (0.30) or non-canonical (1.50)
    assert out.loc["GENEB", "loeuf"] == pytest.approx(1.20)  # canonical fallback
    assert out.loc["GENEE", "loeuf"] == pytest.approx(0.85)  # MANE beats a different canonical transcript
    assert out.index.is_unique
