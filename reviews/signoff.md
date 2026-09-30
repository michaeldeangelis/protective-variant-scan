# Reviewer signoff: protective-variant scan v1 (ledger conformance)

Reviewer: independent of the builders. Scope: experiments.md entry 2026-09-29 protective-variant-scan, Amendment 1, clarifications C1-C10.
Reviewed code state: HEAD 6b959b6 (C10 pipeline) with 51b8e05 (data adapters); reviewer tests at 210d3db.

## Verdict

Conforms to the ledger as written, for the synthetic pipeline path and for the code paths the real data will use, with the residual items below. No open item can produce a false PASS by itself; the open items that matter can cause a false KILL, a masked control, or a degenerate replication row counted as evidence (R2-2). Recommendation: resolve R2-2 (small code change) and decide R2-3 and R2-4 (ledger decisions) before the first real-data run. The other open items are Low.

## Verified by running code
- tests/test_ledger_conformance.py: 250 passed, 4 xfailed (documented residual gaps R2-1, R2-2 x2, R2-4). Full repository suite at the same state: 407 passed, 4 xfailed.
- Every constant in the file is typed from the ledger (not read from config/prereg.yaml): thresholds, 13-trait panel with directions and domains, 9 trade-offs, proxies, C1 to C10 numbers, lipid list frozen at 4073b89, and a second copy of the config file sha256.
- An independent oracle (re-implementation of discovery, replication, adverse, C1/C9 screening, tiers, controls including C10a to C10b, verdict) agrees with the pipeline on 8 synthetic scenarios run through python -m protscan run and on 25 fuzz tables.
- Verdict logic exhaustively checked over all 256 (tier, qualifies) combinations: a LEAD is never reported as PASS; any failed or not-run control gives KILL.
- Mutation testing: 24 mutants of the earlier pipeline and 20 mutants of the C10 code, all killed by the reviewer file.
- Review-1 findings R-1 to R-11 each closed by a test that fails when the fix is removed (review-2 section 1).
- Adversarial probes run: config pin bypass variants (comment, trailing newline, extra key, CRLF, identical bytes elsewhere), sign-control gaming, coverage-gate boundaries (90 percent, 10,000 rows), NaN se end to end, A1FREQ boundaries, degenerate FinnGen se, deflation masking of lambda_GC, mixed synthetic and real inputs.

## Not verified
- No real data outcome was opened. Real Genebass and FinnGen sign conventions, real synonymous lambda_GC, real coverage of the synonymous table, real PCSK9 and LDLR replication signs, and the power of the replication-sign arms are all unchecked.
- Nothing was executed over the network (fetch_genebass.py, fetch_finngen.py, fetch_gnomad.py, download helpers).
- A clean-environment install of the declared dependencies was not run (pyproject text checked only).
- Behaviour and memory at real scale (hundreds of thousands of rows) were not measured.
- The Genebass and FinnGen adapters were checked against the builders' fixtures, not against live responses.
- The AstraZeneca portal consistency check is not automated and was not assessed.

## Open items (details in reviews/review-2.md)
- R2-2 Medium: FinnGen rows with undefined se keep their p and can reach Tier A. Fix in adapter (keep only beta 0 or p 1).
- R2-3 Medium: replication-sign control is direction-only; a low-powered PCSK9 arm can produce a false KILL. Ledger decision.
- R2-4 Medium: p = 1 rows can mask synonymous inflation in lambda_GC. Ledger decision.
- R2-1, R2-5 to R2-10 Low: hit-gene omission inside the 90 percent coverage, PCSK9 Tier C via type 2 diabetes, self-referential pin, unreported A1FREQ drops, PASS-qualifying list under KILL, config-relative builder tests, gnomAD transcript choice.

## Conditions on any reported result
- The run must use config/prereg.yaml with the pinned hash (report and JSON show config_matches_ledger true) and must not be labelled NON-PREREGISTERED or SYNTHETIC.
- controls.negative_synonymous must show coverage and row counts, and the fraction of p = 1 rows should be read alongside lambda_GC (R2-4).
- A KILL with a failed replication_sign or pcsk9_ldl_lower control should be read together with R2-3 and R2-5 before it is called a biological or pipeline conclusion.
- PASS or LEAD is a candidate only: cognitive results cap at LEAD, and PASS is reachable only through systolic blood pressure with the FinnGen hypertension proxy (C2).
