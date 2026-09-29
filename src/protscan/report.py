"""Markdown report from the results dict."""
from __future__ import annotations

MAX_ROWS = 40


def _fmt(x, spec=".3g"):
    return "NA" if x is None else format(x, spec)


def _tier_table(rows: list, n_tradeoff: int) -> list:
    if not rows:
        return ["(none)", ""]
    out = ["| gene | trait | beta | p | replication | dmis | EUR-only | adverse | screened | lipid | qualifies | LOEUF | note |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows[:MAX_ROWS]:
        rep = r["rep_status"] if r["rep_trait"] is None else f"{r['rep_status']} ({r['rep_trait']}, p1={_fmt(r['rep_p_onesided'])})"
        adv = ", ".join(a["trait"] for a in r["adverse"]) or "-"
        out.append(f"| {r['gene']} | {r['trait']} | {_fmt(r['beta'])} | {_fmt(r['p'])} | {rep} | {r['mask_status']} | "
                   f"{r['eur_status']} | {adv} | {r['n_tradeoff_tested']}/{n_tradeoff} | {'yes' if r['lipid_gene'] else 'no'} | {'yes' if r['qualifies'] else 'no'} | "
                   f"{_fmt(r.get('loeuf'))} | {'trade-off unscreened' if r['tradeoff_unscreened'] else ''} |")
    if len(rows) > MAX_ROWS:
        out.append(f"... {len(rows) - MAX_ROWS} more in the JSON")
    return out + [""]


def render(res: dict) -> str:
    v, c, rg, d = res["verdict"], res["controls"], res["rungs"], res["data"]
    L = ["# Protective-variant scan: report", ""]
    if d["synthetic"]:
        L += ["**SYNTHETIC DATA: invented numbers, no biological meaning. Not a result of the preregistered scan.**", ""]
    L += [f"Ledger: {res['ledger_entry']}  ", f"Config sha256: {res['config_sha256'][:16]}  ",
          f"Data: {d['dir']} ({', '.join(d['files'])}; {d['n_rows']} rows; {d['dropped_nan_rows']} NaN rows dropped)  ",
          "Sources: " + "; ".join(f"{k}: {', '.join(s)}" for k, s in d["sources"].items()), ""]

    L += [f"## Verdict: {v['verdict']}", "", v["reason"], ""]
    if v["verdict"] == "LEAD":
        L += ["UNREPLICATED lead only. Not a replicated protective lever.", ""]
    if v["verdict"] == "KILL" and v["reason"].startswith("controls_not_evaluable"):
        L += ["Pipeline unvalidated because a control could not be run; this is not a negative biological result.", ""]
    L += [f"- PASS-qualifying genes (Tier A): {', '.join(v['pass_genes']) or 'none'}",
          f"- LEAD-qualifying genes (Tier B): {', '.join(v['lead_genes']) or 'none'}",
          "- Qualifying = non-lipid-pathway gene, cognitive or physical trait, both masks consistent. "
          "Tier-A genes are candidates, never proven levers.", ""]

    L += ["## Controls", "", "| control | status | detail |", "|---|---|---|"]
    n = c["negative_synonymous"]
    if n["status"] == "NOT RUN":
        L.append(f"| synonymous mask | NOT RUN | {n['reason']} |")
    else:
        L.append(f"| synonymous lambda_GC < {_fmt(n['lambda_gc_max'])} | {'OK' if n['lambda_ok'] else 'FAIL'} | "
                 f"lambda_GC = {n['lambda_gc']:.3f} over {n['n_rows']} rows |")
        L.append(f"| synonymous genes at discovery threshold (beneficial) <= {n['syn_hits_max']} | "
                 f"{'OK' if n['syn_hits_ok'] else 'FAIL'} | {n['n_syn_hit_genes']} genes {', '.join(n['syn_hit_genes'][:10])} |")
    for p in c["positive"]:
        det = "; ".join(f"{t['gene']} beta={t['beta']:.3g} p={t['p']:.2g}" for t in p["tested"]) or "no rows for these genes/trait"
        thr = "at discovery threshold" if p["at_discovery_threshold"] else "beneficial direction"
        L.append(f"| {p['id']} ({'/'.join(p['genes'])}, {p['trait']}, {thr}) | {p['status']} | {det} |")
    L += ["", f"Controls valid: **{c['valid']}**", ""]

    L += ["## Baseline ladder (genes passing at each rung)", "", "| rung | status | genes | note |", "|---|---|---|---|"]
    t, s, i, k = rg["trivial"], rg["simplest"], rg["incumbent"], rg["candidate"]
    L.append(f"| trivial (synonymous mask) | {t['status']} | {t['n_genes'] if t['n_genes'] is not None else 'NA'} | expected 0 |")
    L.append(f"| simplest (pLoF, per trait, discovery p only) | {s['status']} | {s['n_genes']} | {s['n_gene_trait_pairs']} gene-trait pairs |")
    if i["status"] == "NOT RUN":
        L.append(f"| incumbent (published top hits) | NOT RUN | NA | {i['reason']} |")
    else:
        L.append(f"| incumbent (published top hits) | RUN | {i['n_genes']} | {i['n_gene_trait_pairs']} pairs; "
                 f"{i['n_recovered_by_discovery']} recovered by discovery; {i['n_recovered_tier_A']} in Tier A; "
                 f"{i['n_tier_A_not_in_incumbent']} Tier-A pairs not in incumbent list |")
    bt = k["n_genes_by_tier"]
    L.append(f"| candidate (full pipeline, Tier A/B with dmis consistent) | {k['status']} | {k['n_genes']} | "
             f"genes by tier A/B/C/D = {bt['A']}/{bt['B']}/{bt['C']}/{bt['D']} |")
    L.append("")

    names = {"A": "discovery + independent replication + no adverse trade-off",
             "B": "discovery, not replicable, no adverse trade-off",
             "C": "discovery with an adverse trade-off",
             "D": "discovery, replication attempted and failed"}
    L += ["## Tiers", ""]
    for tier in "ABCD":
        rows = res["tiers"][tier]
        L += [f"### Tier {tier}: {names[tier]} (genes: {len({r['gene'] for r in rows})}, gene-trait pairs: {len(rows)})", ""]
        L += _tier_table(rows, res["thresholds"]["n_tradeoff"])

    L += ["## Notes", "",
          "- Replication counts only a cohort independent of UK Biobank; only traits with a declared proxy (LDL, BMI, SBP) can reach Tier A. "
          "Within-UKB consistency (dmis mask) is reported and never counted as replication; AstraZeneca portal lookups are not automated.",
          "- EUR-only column is reported and does not gate the verdict. 'screened' = trade-off outcomes with a pLoF row for the gene; "
          "unscreened outcomes cannot be excluded as adverse.",
          f"- Tier A needs >= {res['thresholds']['min_tradeoffs_screened']} of {res['thresholds']['n_tradeoff']} trade-off outcomes screened "
          "(pLoF) for the gene; otherwise a replicated gene is capped at Tier B and labeled 'trade-off unscreened'.",
          f"- Traits absent from discovery: {', '.join(d['traits_absent']['discovery']) or 'none'}.",
          f"- Traits absent from replication: {', '.join(d['traits_absent']['replication']) or 'none'}.",
          ""]
    return "\n".join(L)
