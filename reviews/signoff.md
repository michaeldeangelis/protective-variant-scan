# Reviewer signoff: protective-variant scan v1 (ledger conformance)

Reviewer: independent of the builders. Scope: experiments.md entry 2026-09-29 protective-variant-scan, Amendment 1 (C1-C11), RUN 1 RESULT and Amendment 2 (A2a-A2d).
Reviewed code state: pipeline fbc0053 (Amendment 2) on top of 49c5ee4 / ec033d2; config sha256 be69827539daac4c78b4c8f6ac7d931723abb41ab0e739c0505c357689310a9f; reviewer tests at the commit after 4d721a5.
Supersedes the signoff for run 1. Details: reviews/review-1.md to review-4.md.

## Decision on run 2
- Suite: green. tests/test_ledger_conformance.py 368 passed, 5 xfailed. Full repository 586 passed, 5 xfailed. The 5 xfails are documented residual gaps (R2-1, R3-1, R3-2, R4-1, R4-3), not hidden failures.
- Run 2: YES, signed off, under the ledger rules as written in Amendment 2, on the conditions below. This is the re-verification A2d asks for before run 2 starts. Nothing found can produce a false PASS by itself; the A2 filter and bound behave exactly as A2a and A2b state.
- Honest framing that must stay attached: Amendment 2 changed the negative control after run 1 failed it, and its 0.1 percent bound (about 18 genes, six times the 3 seen) is post-hoc. Run 1 stands as KILL under the original rule. A PASS or LEAD in run 2 is "passed under Amendment 2, which was written after run 1 failed", never a clean preregistered pass.

## Conditions on run 2 and on anything reported from it
1. Use config/prereg.yaml unchanged: report and JSON must show config_matches_ledger true, the label "Amendment 2 (post-hoc after run 1)" and no NON-PREREGISTERED or SYNTHETIC label. Do not change any rule, bound or list after seeing run 2; a change needs a new dated entry and a new run.
2. Read the contamination table before any hit list: how many genes, which traits, whether they cluster at one locus or one trait. A KILL from the A2b bound is systemic contamination; a single unreliable trait (rare binary outcome) can cause it and is visible in the table.
3. For every PASS or LEAD gene, by hand: (a) A2c, chromosome position and GWAS Catalog lookup of the same trait within 500 kb (the pipeline does not carry coordinates, R4-4); (b) confirm the gene has a synonymous row and note its smallest synonymous p, because a gene without one cannot be assessed by A2a (R4-1) and a sub-threshold signal is not excluded; (c) confirm at least 5 of its trade-off rows are informative (finite se, beta not 0, p not 1), because uninformative rows still count as screened (R3-1).
4. Read a KILL that involves the sign control or the positive controls with the C11 reading rule; a biology conclusion needs all controls valid.
5. PASS or LEAD remains a candidate only: cognitive results cap at LEAD, PASS is reachable only through systolic blood pressure with the FinnGen hypertension proxy (C2), and the C5 and C7 caveats accompany it.
6. Do not open results/run1-*.json tier lists to steer the choice of rules; the reviewer did not open them.

## Verified by running code
- Pin recomputed independently (shasum) and config diff against the previous version: syn_hits_max removed, syn_contaminated_max_fraction 0.001 added, nothing else.
- Ledger constants typed from the ledger (not from config), including the A2 numbers, and a literal copy of the config sha256.
- Independent oracle (discovery, contamination filter, replication with informative rows, adverse, C1/C9 screening, tiers, controls with coverage, p = 1 share, powered sign arms and the exact integer 0.1 percent bound, verdict) agrees with the pipeline on 11 synthetic scenarios through python -m protscan run, on 25 fuzz tables with random synonymous hits, and on targeted mutations.
- A2a and A2b boundaries: strict discovery threshold, either direction, all 23 traits, cohort and mask limits, one contaminated gene per tier with a clean twin, 1/1000, 2/2000, 3/3000 pass and 2/1000, 3/2000, 4/3000 fail, denominator is genes not rows, gene counted once, contaminated control genes still pass their controls, disclosure in verdict block, label, report, JSON and CLI.
- Verdict logic exhaustive over 256 (tier, qualifies) combinations: a LEAD is never reported as PASS.
- Mutation testing: 24 (C1-C9), 20 (C10), 21 (C11) and 18 (A2) mutants, all killed.
- Earlier closures (review-1 to review-3) still hold in the current state.

## Not verified
- No real data outcome was opened, including run 1 tier lists and everything under data/. Real synonymous contamination, its loci, real lambda, coverage, p = 1 share, real sign-arm power and real drop counts are unchecked. The run 1 numbers quoted in the ledger are the lead's statements.
- The claim in the ledger that CD3EAP, ZNF224 and HPR sit in LD with common-variant loci is unverified (the ledger marks it so) and was not assessed.
- Nothing was executed over the network. A clean-environment install of the declared dependencies was not run. Behaviour and memory at real scale were not measured. Adapters were checked against fixtures, not live responses. The AstraZeneca portal check is not automated and was not assessed.

## Open items
- R4-1 Low-Medium: hit gene without a synonymous row is unassessed and unflagged (fix: syn_assessed and min_syn_p columns).
- R3-1 Low-Medium: uninformative trade-off rows count as screened. R3-2 Low-Medium: p just below 1 evades the p = 1 filter.
- R4-2 Low: verdict reason does not say when A2a removed a would-be PASS or LEAD gene. R4-3 Low: header still reads "with AMENDMENT 1".
- R4-4 Info: A2c needs coordinates the pipeline drops. R4-5 Info: post-hoc bound six times the run 1 count.
- R2-1, R3-3, R3-4, R2-6, R2-8, R2-9, R2-10 Low or Info: see reviews 2 and 3.
