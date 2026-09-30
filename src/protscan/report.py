"""Markdown report from the results dict."""
from __future__ import annotations

MAX_ROWS = 40


def _fmt(x, spec=".3g"):
    return "NA" if x is None else format(x, spec)


def _tier_table(rows: list, n_tradeoff: int) -> list:
    if not rows:
        return ["(none)", ""]
    out = ["| gene | trait | beta | p | replication | dmis | EUR-only | adverse | screened | harmful panel (info) | lipid | qualifies | LOEUF | note |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows[:MAX_ROWS]:
        rep = r["rep_status"] if r["rep_trait"] is None else f"{r['rep_status']} ({r['rep_trait']}, p1={_fmt(r['rep_p_onesided'])})"
        adv = ", ".join(a["trait"] for a in r["adverse"]) or "-"
        hp = ", ".join(f"{h['trait']} (p={_fmt(h['p'])})" for h in r["harmful_panel"]) or "-"
        out.append(f"| {r['gene']} | {r['trait']} | {_fmt(r['beta'])} | {_fmt(r['p'])} | {rep} | {r['mask_status']} | "
                   f"{r['eur_status']} | {adv} | {r['n_tradeoff_tested']}/{n_tradeoff} | {hp} | {'yes' if r['lipid_gene'] else 'no'} | {'yes' if r['qualifies'] else 'no'} | "
                   f"{_fmt(r.get('loeuf'))} | {'trade-off unscreened' if r['tradeoff_unscreened'] else ''} |")
    if len(rows) > MAX_ROWS:
        out.append(f"... {len(rows) - MAX_ROWS} more in the JSON")
    return out + [""]


def render(res: dict) -> str:
    v, c, rg, d = res["verdict"], res["controls"], res["rungs"], res["data"]
    L = ["# Protective-variant scan: report", ""]
    if d["synthetic"]:
        L += ["**SYNTHETIC DATA: invented numbers, no biological meaning. Not a result of the preregistered scan.**", ""]
    if not res["config_matches_ledger"]:
        L += [f"**NON-PREREGISTERED: config sha256 {res['config_sha256'][:16]} does not match the pinned ledger value "
              f"{res['config_sha256_pinned'][:16]}. Not a result of the preregistered scan.**", ""]
    L += [f"Ledger: {res['ledger_entry']}  ", f"Config sha256: {res['config_sha256'][:16]}  ",
          f"Data: {d['dir']} ({', '.join(d['files'])}; {d['n_rows']} rows; {d['dropped_nan_rows']} NaN rows dropped)  ",
          "Sources: " + "; ".join(f"{k}: {', '.join(s)}" for k, s in d["sources"].items()), ""]

    L += [f"## Verdict: {v['label']}", "", v["reason"], ""]
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
        L.append(f"| synonymous mask (coverage/size gate) | NOT RUN | {n['reason']} |")
    else:
        L.append(f"| synonymous coverage >= {n['syn_min_coverage']} of pLoF pairs, rows >= {n['syn_min_rows']} | OK | "
                 f"coverage = {n['coverage']:.3f} of {n['n_plof_pairs']} pairs; {n['n_rows']} rows |")
        L.append(f"| synonymous lambda_GC < {_fmt(n['lambda_gc_max'])} | {'OK' if n['lambda_ok'] else 'FAIL'} | "
                 f"lambda_GC = {n['lambda_gc']:.3f} on rows with p < 1; {n['lambda_gc_all_rows']:.3f} on all {n['n_rows']} rows; "
                 f"p = 1 rows: {n['n_p1_rows']} ({n['p1_fraction']:.3f}, max {n['syn_max_p1_fraction']}) |")
        L.append(f"| synonymous genes at discovery threshold (beneficial) <= {n['syn_hits_max']} | "
                 f"{'OK' if n['syn_hits_ok'] else 'FAIL'} | {n['n_syn_hit_genes']} genes {', '.join(n['syn_hit_genes'][:10])} |")
    for p in c["positive"]:
        det = "; ".join(f"{t['gene']} beta={t['beta']:.3g} p={t['p']:.2g}"
                        + (f" tier={t['pipeline_tier']}" if "pipeline_tier" in t else "") for t in p["tested"]) or "no rows for these genes/trait"
        thr = "at discovery threshold" if p["at_discovery_threshold"] else "beneficial direction"
        thr += ", also via the pipeline discovery-hit path (trade-offs never change status)" if p["via_pipeline_discovery_path"] else ""
        L.append(f"| {p['id']} ({'/'.join(p['genes'])}, {p['trait']}, {thr}) | {p['status']} | {det} |")
    rs = c["replication_sign"]
    for k in rs["checks"]:
        det = (f"beta={k['beta']:.3g} p={k['p']:.2g}" if "beta" in k else "no informative row in independent replication cohort")
        det += f"; not evaluable: {k['why']}" if k["status"] == "NOT RUN" else ""
        want = "negative" if k["expected_beta_sign"] < 0 else "positive"
        L.append(f"| {k['id']} ({k['gene']} pLoF, {k['trait']}, expect {want}) | {k['status']} | {det} |")
    L += ["", f"Replication-sign control (an arm is evaluable only if p < {_fmt(res['thresholds']['replication_sign_max_p'])}; "
          f"at least one arm evaluable, every evaluable arm must match): **{rs['status']}**"
          + (f" ({rs['reason']})" if rs["reason"] else ""), "",
          f"Controls valid: **{c['valid']}**", ""]
    if v["verdict"] == "KILL" and (c["failed"] or c["not_run"]):
        L += ["Reading rule: a KILL involving the sign control or the positive controls is a pipeline/data conclusion only after checking "
              "that the sign-control arms were adequately powered (C11c) and that the positive controls were judged on the discovery effect "
              "alone (C11d). A biology conclusion requires all controls valid.", ""]

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

    if v["verdict"] in ("PASS", "LEAD"):
        L += ["## Caveats required with any PASS or LEAD (C5, C7)", "",
              "- C5: the `dmis` mask is Genebass missense|LC (missense plus low-confidence pLoF), not a damaging-missense filter. "
              "'Both masks consistent' is direction-only agreement between pLoF and missense|LC, weaker than intended.",
              "- C7: cognitive coverage is limited to fluid intelligence and reaction time; numeric memory, pairs matching, "
              "education years, walking pace and all-cause mortality are absent from the discovery source.",
              ""]
    L += ["## Notes", "",
          "- Replication counts only a cohort independent of UK Biobank; only traits with a declared proxy (LDL, BMI, SBP) can reach Tier A. "
          "Within-UKB consistency (dmis mask) is reported and never counted as replication; AstraZeneca portal lookups are not automated.",
          "- EUR-only column is reported and does not gate the verdict. 'screened' = trade-off outcomes with a pLoF row for the gene; "
          "unscreened outcomes cannot be excluded as adverse.",
          f"- Tier A needs >= {res['thresholds']['min_tradeoffs_screened']} of {res['thresholds']['n_tradeoff']} trade-off outcomes screened "
          "(pLoF) for the gene (discovery plus allow-listed independent replication rows); otherwise a replicated gene is capped at Tier B and labeled 'trade-off unscreened'. "
          "The same minimum applies to LEAD (C9): unscreened Tier-B genes are listed but do not count.",
          f"- Traits absent from discovery: {', '.join(d['traits_absent']['discovery']) or 'none'}.",
          f"- Traits absent from replication: {', '.join(d['traits_absent']['replication']) or 'none'}.",
          ""]
    return "\n".join(L)
