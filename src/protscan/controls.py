"""Negative (synonymous mask), positive (lipid) and replication-sign controls. Status: OK | FAIL | NOT RUN."""
from __future__ import annotations

import pandas as pd

from . import stats, tiering
from .schema import Config


def negative_control(df: pd.DataFrame, cfg: Config) -> dict:
    syn = df[(df["cohort"] == "discovery") & (df["mask"] == "syn")]
    if syn.empty:
        return {"status": "NOT RUN", "reason": "no synonymous-mask rows in the discovery cohort"}
    plof = df.loc[(df["cohort"] == "discovery") & (df["mask"] == "plof"), ["gene", "trait"]]
    if plof.empty:
        coverage, n_pairs = 0.0, 0
    else:
        m = plof.merge(syn[["gene", "trait"]], how="left", indicator=True)
        coverage, n_pairs = float((m["_merge"] == "both").mean()), int(len(plof))
    cov = {"n_rows": int(len(syn)), "syn_min_rows": cfg.syn_min_rows, "coverage": coverage,
           "syn_min_coverage": cfg.syn_min_coverage, "n_plof_pairs": n_pairs}
    if len(syn) < cfg.syn_min_rows or coverage < cfg.syn_min_coverage:
        return {"status": "NOT RUN", **cov,
                "reason": f"synonymous universe too small to validate the pipeline: {len(syn)} rows (need >= {cfg.syn_min_rows}), "
                          f"coverage {coverage:.3f} of pLoF (gene, trait) pairs (need >= {cfg.syn_min_coverage})"}
    lam = stats.lambda_gc(syn["p"])
    hits = stats.discovery_hits(df, cfg, mask="syn")
    genes = sorted(set(hits["gene"]))
    lam_ok = lam < cfg.lambda_gc_max
    hits_ok = len(genes) <= cfg.syn_hits_max
    return {
        "status": "OK" if lam_ok and hits_ok else "FAIL",
        "lambda_gc": lam, "lambda_gc_max": cfg.lambda_gc_max, "lambda_ok": lam_ok,
        **cov,
        "n_syn_hit_genes": len(genes), "syn_hits_max": cfg.syn_hits_max, "syn_hits_ok": hits_ok,
        "syn_hit_genes": genes[:50],
    }


def positive_controls(df: pd.DataFrame, cfg: Config, hits: pd.DataFrame | None = None) -> list:
    """Direction and threshold tests on the discovery plof rows. Where the trait is a panel trait and the spec
    is at the discovery threshold, the gene must also come out of the pipeline's own hit and tier path and not be Tier C."""
    if hits is None:
        hits = tiering.build_hits(df, cfg)
    tier = {(g, t): k for g, t, k in zip(hits["gene"], hits["trait"], hits["tier"])}
    d = df[(df["cohort"] == "discovery") & (df["mask"] == "plof")].set_index(["gene", "trait"])
    out = []
    for spec in cfg.positive_controls:
        trait, sign = spec["trait"], cfg.sign(spec["trait"])
        via_pipeline = bool(spec["at_discovery_threshold"] and trait in cfg.panel)
        tested, passed = [], False
        for g in spec["genes"]:
            if (g.upper(), trait) not in d.index:
                continue
            row = d.loc[(g.upper(), trait)]
            ok = bool(row["beta"] * sign > 0 and (row["p"] < cfg.discovery_p or not spec["at_discovery_threshold"]))
            rec = {"gene": g, "beta": float(row["beta"]), "p": float(row["p"])}
            if via_pipeline:
                rec["pipeline_tier"] = tier.get((g.upper(), trait))
                ok = ok and rec["pipeline_tier"] is not None and rec["pipeline_tier"] != "C"
            rec["ok"] = ok
            tested.append(rec)
            passed = passed or ok
        out.append({
            "id": spec["id"], "trait": trait, "genes": list(spec["genes"]),
            "at_discovery_threshold": bool(spec["at_discovery_threshold"]), "via_pipeline_tier_path": via_pipeline,
            "status": ("OK" if passed else "FAIL") if tested else "NOT RUN",
            "tested": tested,
        })
    return out


def replication_sign_control(df: pd.DataFrame, cfg: Config) -> dict:
    """C10b: PCSK9 negative, LDLR positive on the hypercholesterolemia proxy in the independent replication cohort.
    At least one check must be evaluable; every evaluable one must hold."""
    rep = stats.independent_replication(df, cfg).set_index(["gene", "trait"])
    checks = []
    for spec in cfg.replication_sign:
        key = (spec["gene"].upper(), spec["trait"])
        if key not in rep.index:
            checks.append({"id": spec["id"], "gene": spec["gene"], "trait": spec["trait"],
                           "expected_beta_sign": spec["expected_beta_sign"], "status": "NOT RUN"})
            continue
        beta = float(rep.loc[key, "beta"])
        checks.append({"id": spec["id"], "gene": spec["gene"], "trait": spec["trait"],
                       "expected_beta_sign": spec["expected_beta_sign"], "beta": beta,
                       "status": "OK" if beta * spec["expected_beta_sign"] > 0 else "FAIL"})
    sts = [c["status"] for c in checks]
    status = "FAIL" if "FAIL" in sts else ("OK" if "OK" in sts else "NOT RUN")
    return {"status": status, "checks": checks}


def run_controls(df: pd.DataFrame, cfg: Config, hits: pd.DataFrame | None = None) -> dict:
    neg = negative_control(df, cfg)
    pos = positive_controls(df, cfg, hits)
    rsign = replication_sign_control(df, cfg)
    names = ["negative_synonymous"] + [c["id"] for c in pos] + ["replication_sign"]
    statuses = [neg["status"]] + [c["status"] for c in pos] + [rsign["status"]]
    return {
        "negative_synonymous": neg,
        "positive": pos,
        "replication_sign": rsign,
        "valid": all(s == "OK" for s in statuses),
        "failed": [n for n, s in zip(names, statuses) if s == "FAIL"],
        "not_run": [n for n, s in zip(names, statuses) if s == "NOT RUN"],
    }
