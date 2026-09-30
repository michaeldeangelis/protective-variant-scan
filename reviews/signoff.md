# Reviewer signoff: protective-variant scan v1 (ledger conformance)

Reviewer: independent of the builders. Scope: experiments.md entry 2026-09-29 protective-variant-scan, Amendment 1, clarifications C1-C11.
Reviewed code state: pipeline 49c5ee4 (C11) with adapters ec033d2; config sha256 7ce18c381814ba0bdfbd56ff6afa236a413bba5752517ebb901b1de74f834b83; reviewer tests at the commit after f361142.
Supersedes the earlier signoff (C10 state). Details: reviews/review-1.md, review-2.md, review-3.md.

## Verdict
- Suite: green. tests/test_ledger_conformance.py 284 passed, 3 xfailed; full repository 479 passed, 3 xfailed. The 3 xfails are documented residual gaps (R2-1, R3-1, R3-2), not hidden failures.
- Ledger conformance: the pipeline conforms to the ledger as written through C11a-d. No open item creates a false PASS by itself.
- First real run: SIGNED OFF as the first real run, on these conditions. (1) Use config/prereg.yaml unchanged (report and JSON must show config_matches_ledger true, no NON-PREREGISTERED or SYNTHETIC label). (2) For every PASS or LEAD gene, hand-check that at least 5 of its trade-off rows are informative (finite se, beta not 0, p not 1), because R3-1 lets uninformative rows count as screened. (3) Read the synonymous control detail (both lambdas, p = 1 share, coverage) before any KILL or hit list, and read a KILL involving the sign control or the positive controls with the C11c and C11d reading rule. R3-1 is a small code change that should be made before the result is treated as final; it does not block a first run.

## Verified by running code
- Ledger constants typed from the ledger (not from config): thresholds, panel and directions, 9 trade-offs, proxies, C1 to C11 numbers, lipid list frozen at 4073b89, and a literal copy of the config file sha256. The pin was recomputed independently with shasum, and git diff of config/prereg.yaml against the C10 version shows exactly two added keys (syn_max_p1_fraction 0.05, replication_sign_max_p 0.05).
- Independent oracle (discovery, replication with informative rows, adverse, C1/C9 screening, tiers, controls including coverage, p = 1 share, powered sign arms, verdict) agrees with the pipeline on 8 synthetic scenarios through python -m protscan run, on 25 fuzz tables with undefined-se and beta-zero replication rows, and on targeted mutations.
- Verdict logic exhaustive over all 256 (tier, qualifies) combinations: a LEAD is never reported as PASS; failed or not-run controls give KILL.
- Boundaries tested: discovery p, one-sided replication 0.05, adverse alpha/9 one-sided, 5 screened outcomes, synonymous coverage 90 percent and 10,000 rows, p = 1 share 5 percent, sign-arm p 0.05, lambda 1.10, A1FREQ 0.5.
- Mutation testing: 24 mutants (C1-C9 code), 20 mutants (C10 code), 21 mutants (C11 code); all killed (one C11 survivor was found and closed with a new test).
- Fix-removal check: each test converted from xfail to a hard test was shown to fail when its fix is removed.
- Adversarial probes run: config pin bypass variants, sign-control gaming, coverage-gate edges, NaN se end to end, degenerate FinnGen se, deflation masking (exact p = 1 closed, p just below 1 open), A1FREQ boundaries, mixed synthetic and real inputs, PCSK9 trade-off harm.
- Builders' test files read for tautologies (review-2 R2-9): config-relative direction in synthetic fixtures noted; independent direction checks exist in the reviewer file.

## Not verified
- No real data outcome was opened. Real Genebass and FinnGen sign conventions, real synonymous lambda_GC, coverage and p = 1 share, real PCSK9 and LDLR arm power, and real drop counts are unchecked. Counts quoted in docs/data-sources.md (zero rows dropped) are the builders' statements.
- Nothing was executed over the network (fetch scripts, download helpers, gnomAD download).
- A clean-environment install of the declared dependencies was not run (pyproject text only).
- Behaviour and memory at real scale (hundreds of thousands of rows) were not measured.
- Adapters were checked against the builders' fixtures and reviewer-built rows, not live responses.
- The AstraZeneca portal consistency check is not automated and was not assessed.

## Open items
- R3-1 Low-Medium: uninformative trade-off rows count as screened (fix: apply informative() in tradeoff_tested).
- R3-2 Low-Medium: p just below 1 evades the C11b filter (fix: gate the share of p above 0.99).
- R2-1 Low: hit pairs may be omitted inside the 90 percent coverage.
- R3-3 Info: discovery_hits is not informativeness-filtered (adapter-guarded); adverse screen intentionally unfiltered.
- R3-4 / R2-7 Low: adapter drop counts are not in the result record.
- R2-6, R2-8, R2-9, R2-10 Low: self-referential pin (mitigated by the reviewer literal), PASS-gene list printed under KILL, config-relative builder tests, gnomAD transcript order.

## Conditions on any reported result
- PASS or LEAD is a candidate only: cognitive results cap at LEAD, and PASS is reachable only through systolic blood pressure with the FinnGen hypertension proxy (C2). The C5 and C7 caveats must accompany it.
- A KILL from a failed or not-evaluable control is a pipeline or data conclusion, not a biology conclusion (C11 reading rule).
