"""End-to-end run: load, controls, ladder rungs, tiers, verdict, JSON + report."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import controls, report, stats, tiering
from .schema import Config, load_burden, load_config


def _incumbent(data_dir: Path, cfg: Config, hits: pd.DataFrame, disc: pd.DataFrame) -> dict:
    path = data_dir / "incumbent_hits.csv"
    if not path.exists():
        return {"status": "NOT RUN", "reason": f"{path.name} not provided; published top hits are not fabricated"}
    inc = pd.read_csv(path, dtype=str)
    if not {"gene", "trait"} <= set(inc.columns):
        raise ValueError("incumbent_hits.csv needs columns gene, trait")
    inc["gene"] = inc["gene"].str.strip().str.upper()
    inc["trait"] = inc["trait"].str.strip().str.lower()
    in_panel = inc[inc["trait"].isin(cfg.panel)].drop_duplicates(["gene", "trait"])
    pairs = set(zip(in_panel["gene"], in_panel["trait"]))
    disc_pairs = set(zip(disc["gene"], disc["trait"]))
    a = hits[hits["tier"] == "A"]
    a_pairs = set(zip(a["gene"], a["trait"]))
    return {
        "status": "RUN",
        "n_rows_in_file": int(len(inc)),
        "n_out_of_panel_dropped": int((~inc["trait"].isin(cfg.panel)).sum()),
        "n_genes": len({g for g, _ in pairs}),
        "n_gene_trait_pairs": len(pairs),
        "n_recovered_by_discovery": len(pairs & disc_pairs),
        "n_recovered_tier_A": len(pairs & a_pairs),
        "n_tier_A_not_in_incumbent": len(a_pairs - pairs),
    }


def _rungs(cfg, ctrl, disc, hits, incumbent) -> dict:
    neg = ctrl["negative_synonymous"]
    return {
        "trivial": {
            "definition": "synonymous mask, discovery rule (Amendment 1)",
            "status": "RUN" if neg["status"] != "NOT RUN" else "NOT RUN",
            "n_genes": neg.get("n_syn_hit_genes"),
        },
        "simplest": {
            "definition": "pLoF only, one trait at a time, discovery p only; no replication, no trade-off screen",
            "status": "RUN",
            "n_genes": int(disc["gene"].nunique()),
            "n_gene_trait_pairs": int(len(disc)),
            "by_trait": {t: int(n) for t, n in disc["trait"].value_counts().items()},
        },
        "incumbent": incumbent,
        "candidate": {
            "definition": "pLoF + dmis consistency, independent replication, trade-off screen, tiering",
            "status": "RUN",
            "n_genes": int(hits.loc[hits["tier"].isin(["A", "B"]) & (hits["mask_status"] == "consistent"), "gene"].nunique()),
            "n_genes_by_tier": {t: int(hits.loc[hits["tier"] == t, "gene"].nunique()) for t in "ABCD"},
            "n_genes_by_tier_mask_consistent": {
                t: int(hits.loc[(hits["tier"] == t) & (hits["mask_status"] == "consistent"), "gene"].nunique()) for t in "ABCD"},
        },
    }


def decide(ctrl: dict, hits: pd.DataFrame) -> dict:
    """PASS / LEAD / KILL per Amendment 1."""
    q = hits[hits["qualifies"]] if len(hits) else hits
    a = sorted(set(q.loc[q["tier"] == "A", "gene"])) if len(q) else []
    b = sorted(set(q.loc[q["tier"] == "B", "gene"])) if len(q) else []
    base = {"controls_valid": ctrl["valid"], "pass_genes": a, "lead_genes": b}
    if ctrl["failed"]:
        return {**base, "verdict": "KILL", "reason": "control_failed: " + ", ".join(ctrl["failed"])}
    if ctrl["not_run"]:
        return {**base, "verdict": "KILL", "reason": "controls_not_evaluable (pipeline unvalidated): " + ", ".join(ctrl["not_run"])}
    if a:
        return {**base, "verdict": "PASS", "reason": "controls valid; non-lipid Tier-A gene with both masks consistent"}
    if b:
        return {**base, "verdict": "LEAD",
                "reason": "controls valid; no qualifying Tier-A gene; non-lipid Tier-B gene with both masks consistent (UNREPLICATED lead)"}
    return {**base, "verdict": "KILL", "reason": "controls valid; neither PASS nor LEAD"}


def _absent(df: pd.DataFrame, cfg: Config) -> dict:
    disc_expected = set(cfg.panel) | set(cfg.tradeoff) | set(cfg.control_sign)
    rep_expected = {p["trait"] for p in cfg.proxies.values()} | set(cfg.tradeoff)
    have = lambda c: set(df.loc[df["cohort"] == c, "trait"])
    return {"discovery": sorted(disc_expected - have("discovery")), "replication": sorted(rep_expected - have("replication"))}


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


def _clean(rec: dict) -> dict:
    return {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in rec.items()}


def run_pipeline(config_path, data_dir, out_path) -> dict:
    cfg = load_config(config_path)
    data_dir = Path(data_dir)
    df = load_burden(data_dir)
    files, dropped = df.attrs["files"], df.attrs["dropped_nan_rows"]

    hits = tiering.build_hits(df, cfg)
    ctrl = controls.run_controls(df, cfg)
    disc = stats.discovery_hits(df, cfg)
    incumbent = _incumbent(data_dir, cfg, hits, disc)

    constraint = {}
    cpath = data_dir / "constraint.csv"
    if cpath.exists():
        c = pd.read_csv(cpath)
        c["gene"] = c["gene"].astype(str).str.upper()
        constraint = dict(zip(c["gene"], c["loeuf"]))

    def records(t):
        rows = [_clean(r) for r in hits[hits["tier"] == t].to_dict("records")]
        for r in rows:
            r["loeuf"] = constraint.get(r["gene"])
        return rows

    sources = {c: sorted(set(df.loc[df["cohort"] == c, "source"])) for c in sorted(set(df["cohort"]))}
    results = {
        "ledger_entry": cfg.ledger_entry,
        "config_sha256": cfg.sha256,
        "data": {
            "dir": str(data_dir), "files": files, "n_rows": int(len(df)), "dropped_nan_rows": int(dropped),
            "sources": sources,
            "synthetic": any(s.startswith("synthetic") for v in sources.values() for s in v),
            "traits_absent": _absent(df, cfg),
            "constraint_file": bool(constraint),
        },
        "thresholds": {
            "discovery_p": cfg.discovery_p, "replication_one_sided_p": cfg.replication_p,
            "tradeoff_p": cfg.tradeoff_p, "lambda_gc_max": cfg.lambda_gc_max,
        },
        "controls": ctrl,
        "rungs": _rungs(cfg, ctrl, disc, hits, incumbent),
        "tiers": {t: records(t) for t in "ABCD"},
        "verdict": decide(ctrl, hits),
    }

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=_jsonable, allow_nan=False) + "\n")
    (out.parent / "report.md").write_text(report.render(results))
    return results
