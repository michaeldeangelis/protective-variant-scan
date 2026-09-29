#!/usr/bin/env python3
"""Synthetic normalized burden tables with planted signals. SYNTHETIC: no real gene results.

Symbols starting with SYN are invented. PCSK9 / ANGPTL4 / APOC3 carry invented numbers, used only
because the ledger keys the positive controls on those symbols.

Scenarios (expected verdict on the default config):
  pass              PASS   planted Tier-A gene, Tier-B lead, Tier C, Tier D, mask-inconsistent gene
  lead              LEAD   no Tier-A plant; cognitive Tier-B lead remains
  nolead            KILL   controls valid; no Tier-A or Tier-B lead
  broken_positive   KILL   PCSK9 effects removed
  broken_lambda     KILL   synonymous z-scores inflated (lambda_GC ~1.7)
  broken_syn_hit    KILL   one synonymous-mask gene at the discovery threshold, beneficial direction
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from protscan.schema import COLUMNS, load_config

SCENARIOS = ("pass", "lead", "nolead", "broken_positive", "broken_lambda", "broken_syn_hit")
BINARY_SE_SCALE = 3.0
N_TOTAL = {"discovery": 400_000, "replication": 300_000}
CARRIERS = {"plof": 150, "dmis": 1500, "syn": 3000}
DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / "config" / "prereg.yaml"


def _frame(genes, traits, masks, cohort, source, se_scale, rng):
    g, t, m = np.meshgrid(np.arange(len(genes)), np.arange(len(traits)), np.arange(len(masks)), indexing="ij")
    g, t, m = g.ravel(), t.ravel(), m.ravel()
    mask_arr = np.array(masks)[m]
    base = np.array([CARRIERS[x] for x in masks])[m] 
    n_car = np.maximum(3, rng.lognormal(np.log(base), 0.5)).round().astype(int)
    se = np.array([se_scale[x] for x in traits])[t] / np.sqrt(n_car)
    z = rng.standard_normal(len(g))
    return pd.DataFrame({
        "gene": np.array(genes)[g], "trait": np.array(traits)[t], "mask": mask_arr, "cohort": cohort,
        "beta": z * se, "se": se, "p": 2 * norm.sf(np.abs(z)),
        "n_carriers": n_car, "n_total": N_TOTAL[cohort], "source": source,
    })


def _set_z(df, gene, trait, mask, z):
    """Overwrite one row so that beta = z * se (signed z)."""
    i = df.index[(df.gene == gene) & (df.trait == trait) & (df["mask"] == mask)]
    assert len(i) == 1, (gene, trait, mask)
    i = i[0]
    df.loc[i, "beta"] = z * df.loc[i, "se"]
    df.loc[i, "p"] = 2 * norm.sf(abs(z))


def _clip_z(df, genes, traits, zmax=1.5):
    """Keep the planted 'clean' genes free of chance trade-off hits."""
    m = df.gene.isin(genes) & df.trait.isin(traits)
    z = (df.loc[m, "beta"] / df.loc[m, "se"]).clip(-zmax, zmax)
    df.loc[m, "beta"] = z * df.loc[m, "se"]
    df.loc[m, "p"] = 2 * norm.sf(np.abs(z))


def make_synthetic(scenario="pass", seed=20260929, n_genes=3000, config_path=DEFAULT_CONFIG):
    """Return {"burden_synthetic_discovery": df, "burden_synthetic_replication": df}."""
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    cfg = load_config(config_path)
    rng = np.random.default_rng(seed)
    ben = cfg.sign

    def good(trait, z):
        return ben(trait) * abs(z)

    def bad(trait, z):
        return -ben(trait) * abs(z)

    planted = ["PCSK9", "ANGPTL4", "APOC3", "SYNPASS1", "SYNLEAD1", "SYNADV1", "SYNFAIL1", "SYNMASK1", "SYNSYN1"]
    genes = planted + [f"SYNG{i:05d}" for i in range(n_genes)]
    d_traits = list(cfg.panel) + list(cfg.tradeoff) + list(cfg.control_sign)
    se_scale = {t: (BINARY_SE_SCALE if t in cfg.tradeoff else 1.0) for t in d_traits}
    disc = _frame(genes, d_traits, ["plof", "dmis", "syn"], "discovery", "synthetic_discovery", se_scale, rng)

    proxy_traits = [px["trait"] for px in cfg.proxies.values()]
    r_traits = proxy_traits + list(cfg.tradeoff)
    r_scale = {t: BINARY_SE_SCALE for t in r_traits}
    r_genes = planted + [g for g in genes[len(planted):] if rng.random() < 0.5]
    rep = _frame(r_genes, r_traits, ["plof"], "replication", "synthetic_replication", r_scale, rng)

    clean = ["SYNPASS1", "SYNLEAD1", "SYNFAIL1", "SYNMASK1"]
    _clip_z(disc, clean, list(cfg.tradeoff))
    _clip_z(rep, clean, list(cfg.tradeoff))

    # lipid positive controls (invented numbers)
    if scenario != "broken_positive":
        _set_z(disc, "PCSK9", "ldl", "plof", good("ldl", 15))
        _set_z(disc, "PCSK9", "ldl", "dmis", good("ldl", 3))
        _set_z(disc, "PCSK9", "coronary_disease", "plof", good("coronary_disease", 4))
        _set_z(rep, "PCSK9", "hypercholesterolemia", "plof", good("hypercholesterolemia", 6))
    _set_z(disc, "ANGPTL4", "triglycerides", "plof", good("triglycerides", 9))
    _set_z(disc, "APOC3", "triglycerides", "plof", good("triglycerides", 12))

    eur_rows = []

    def hit(gene, trait, z, dmis_z, rep_trait=None, rep_z=None):
        _set_z(disc, gene, trait, "plof", good(trait, z))
        _set_z(disc, gene, trait, "dmis", dmis_z)
        if rep_trait:
            _set_z(rep, gene, rep_trait, "plof", rep_z)
        e = disc[(disc.gene == gene) & (disc.trait == trait) & (disc["mask"] == "plof")].copy()
        e["cohort"] = "discovery_eur"
        e["source"] = "synthetic_discovery_eur"
        e["se"] = e["se"] * 1.03
        e["beta"] = e["beta"] * 0.98
        e["p"] = 2 * norm.sf(abs(e["beta"] / e["se"]))
        eur_rows.append(e)

    if scenario in ("pass", "broken_positive", "broken_lambda", "broken_syn_hit"):
        # beneficial non-lipid gene: SBP lower, replicates via the hypertension proxy, no adverse trade-off
        hit("SYNPASS1", "systolic_bp", 8, good("systolic_bp", 2), "hypertension", good("hypertension", 3))
    if scenario != "nolead":
        # cognitive gene: no declared proxy, so Tier B ceiling; both masks agree
        hit("SYNLEAD1", "fluid_intelligence", 7, good("fluid_intelligence", 2))
    # SBP gene with an adverse trade-off (coronary disease, harmful direction): Tier C
    hit("SYNADV1", "systolic_bp", 6, good("systolic_bp", 2), "hypertension", good("hypertension", 3))
    _set_z(disc, "SYNADV1", "coronary_disease", "plof", bad("coronary_disease", 4))
    # SBP gene that fails replication: Tier D
    hit("SYNFAIL1", "systolic_bp", 6, good("systolic_bp", 2), "hypertension", bad("hypertension", 0.5))
    # FEV1 gene whose damaging-missense mask disagrees in direction: Tier B, mask-inconsistent
    hit("SYNMASK1", "fev1", 7, bad("fev1", 2))

    if scenario == "broken_lambda":
        s = disc["mask"] == "syn"
        z = disc.loc[s, "beta"] / disc.loc[s, "se"] * 1.3
        disc.loc[s, "beta"] = z * disc.loc[s, "se"]
        disc.loc[s, "p"] = 2 * norm.sf(np.abs(z))
    if scenario == "broken_syn_hit":
        _set_z(disc, "SYNSYN1", "hand_grip_strength", "syn", good("hand_grip_strength", 7))

    disc = pd.concat([disc] + eur_rows, ignore_index=True)
    return {
        "burden_synthetic_discovery": disc[COLUMNS],
        "burden_synthetic_replication": rep[COLUMNS],
    }


def write(tables, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(out / f"{name}.csv.gz", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/synthetic")
    ap.add_argument("--scenario", default="pass", choices=SCENARIOS)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--n-genes", type=int, default=3000)
    ap.add_argument("--config", default=str(DEFAULT_CONFIG))
    a = ap.parse_args()
    tables = make_synthetic(a.scenario, a.seed, a.n_genes, a.config)
    write(tables, a.out)
    print(f"wrote {sum(len(t) for t in tables.values())} rows to {a.out} (scenario={a.scenario})")


if __name__ == "__main__":
    main()
