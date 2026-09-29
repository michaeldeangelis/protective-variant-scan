"""Negative (synonymous mask) and positive (lipid) controls. Status: OK | FAIL | NOT RUN."""
from __future__ import annotations

import pandas as pd

from . import stats
from .schema import Config


def negative_control(df: pd.DataFrame, cfg: Config) -> dict:
    syn = df[(df["cohort"] == "discovery") & (df["mask"] == "syn")]
    if syn.empty:
        return {"status": "NOT RUN", "reason": "no synonymous-mask rows in the discovery cohort"}
    lam = stats.lambda_gc(syn["p"])
    hits = stats.discovery_hits(df, cfg, mask="syn")
    genes = sorted(set(hits["gene"]))
    lam_ok = lam < cfg.lambda_gc_max
    hits_ok = len(genes) <= cfg.syn_hits_max
    return {
        "status": "OK" if lam_ok and hits_ok else "FAIL",
        "lambda_gc": lam, "lambda_gc_max": cfg.lambda_gc_max, "lambda_ok": lam_ok,
        "n_rows": int(len(syn)),
        "n_syn_hit_genes": len(genes), "syn_hits_max": cfg.syn_hits_max, "syn_hits_ok": hits_ok,
        "syn_hit_genes": genes[:50],
    }


def positive_controls(df: pd.DataFrame, cfg: Config) -> list:
    d = df[(df["cohort"] == "discovery") & (df["mask"] == "plof")].set_index(["gene", "trait"])
    out = []
    for spec in cfg.positive_controls:
        trait, sign = spec["trait"], cfg.sign(spec["trait"])
        tested, passed = [], False
        for g in spec["genes"]:
            if (g.upper(), trait) not in d.index:
                continue
            row = d.loc[(g.upper(), trait)]
            ok = bool(row["beta"] * sign > 0 and (row["p"] < cfg.discovery_p or not spec["at_discovery_threshold"]))
            tested.append({"gene": g, "beta": float(row["beta"]), "p": float(row["p"]), "ok": ok})
            passed = passed or ok
        out.append({
            "id": spec["id"], "trait": trait, "genes": list(spec["genes"]),
            "at_discovery_threshold": bool(spec["at_discovery_threshold"]),
            "status": ("OK" if passed else "FAIL") if tested else "NOT RUN",
            "tested": tested,
        })
    return out


def run_controls(df: pd.DataFrame, cfg: Config) -> dict:
    neg = negative_control(df, cfg)
    pos = positive_controls(df, cfg)
    statuses = [neg["status"]] + [c["status"] for c in pos]
    return {
        "negative_synonymous": neg,
        "positive": pos,
        "valid": all(s == "OK" for s in statuses),
        "failed": [n for n, s in zip(["negative_synonymous"] + [c["id"] for c in pos], statuses) if s == "FAIL"],
        "not_run": [n for n, s in zip(["negative_synonymous"] + [c["id"] for c in pos], statuses) if s == "NOT RUN"],
    }
