#!/usr/bin/env python
"""Build <data>/coverage.csv from the normalized adapter outputs and print it as a markdown table."""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from protscan import schema  # noqa: E402
from protscan.adapters import coverage  # noqa: E402


def markdown(df: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(ROOT / "data"))
    ap.add_argument("--config", default=str(ROOT / "config" / "prereg.yaml"))
    a = ap.parse_args()
    data = Path(a.data)
    cfg = schema.load_config(a.config)
    disc, rep = (pd.read_csv(data / f, dtype={"gene": str, "trait": str}) for f in ("burden_genebass.csv.gz", "burden_finngen.csv.gz"))
    cov = coverage.build_coverage(cfg, disc, rep)
    cov.to_csv(data / "coverage.csv", index=False)
    print(markdown(cov))


if __name__ == "__main__":
    main()
