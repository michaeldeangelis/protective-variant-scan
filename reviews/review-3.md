# Review 3: diff-scoped re-review of C11 (commits 49c5ee4 pipeline, ec033d2 adapters)

Reviewer: independent (did not write the pipeline). Scope: git diff b30292c..HEAD of src/ and scripts/ only. Ledger basis: experiments.md C11a-d.
No real data outcome was opened. Nothing was fetched from the network. The counts quoted from docs/data-sources.md (rows dropped and kept) are the builders' own statements and were not re-derived.

## 1. Config pin (recomputed independently)
- shasum -a 256 config/prereg.yaml gives 7ce18c381814ba0bdfbd56ff6afa236a413bba5752517ebb901b1de74f834b83, equal to schema.PINNED_CONFIG_SHA256 and to the value in the lead message.
- git diff b30292c HEAD -- config/prereg.yaml: exactly two keys added under controls, syn_max_p1_fraction 0.05 and replication_sign_max_p 0.05, plus their comments. Nothing else changed (thresholds, panel, directions, proxies, lipid list, min screened, coverage and row floors, replication_sign genes and signs are byte-identical).
- The reviewer literal was updated to the new hash; the typed constants L_SYN_MAX_P1_FRACTION and L_REP_SIGN_MAX_P come from the ledger text (C11b, C11c), not from the config, and the config test compares the two.

## 2. Closure of review-2 findings
- R2-2 (undefined-se replication rows counted): CLOSED. stats.informative filters replication_status and the sign control; the FinnGen and Genebass adapters drop undefined-se rows with p below 1 and count them; ivw_combine keeps undefined-se genes only when p is 1. Tests that were xfail now pass as hard tests and were confirmed to fail when the fix is removed (fix-removal run: finngen drop rule, informative(), p = 1 gate plus p<1 lambda each removed, test then fails).
- R2-4 (p = 1 rows mask lambda_GC): CLOSED for exact p = 1. Lambda is computed on p below 1 rows, both lambdas are reported, more than 5 percent p = 1 rows makes the control not evaluable. Residual R3-2 (p just below 1).
- R2-3 (direction-only sign control, low-powered arm): CLOSED. An arm counts only if p is below 0.05; underpowered arms are NOT RUN and never failures; reasons sign_control_failed and sign_control_not_evaluable reach the verdict text. Boundary 0.0499 / 0.05 / 0.0501 tested.
- R2-5 (positive control fails on trade-off biology): CLOSED. The control is judged on the discovery effect only; a Tier C PCSK9 leaves the control OK (unit and end-to-end tests, and an oracle comparison).
- R2-1 (selective omission inside the 90 percent coverage): unchanged, still open, Low (xfail kept).
- R2-6 (self-referential pin): unchanged; mitigated by the reviewer literal, which had to be updated by hand for this change (a working example of the check).
- R2-7 (drop counts not in the result record): still open, see R3-4.
- R2-8 to R2-10 (KILL report lists PASS genes, config-relative builder tests, gnomAD transcript order): unchanged, Low; not re-examined because outside the diff.

## 3. New findings

### R3-1 [Low-Medium] Uninformative rows still count as trade-off outcomes "screened"
File: src/protscan/stats.py:98-105 (tradeoff_tested) counts every plof row of a trade-off outcome. C11a filters replication evidence only. A row with beta 0, p 1 and undefined se is a row where no test was informative (the Genebass adapter keeps such rows on purpose), yet it counts toward the C1/C9 minimum of 5 outcomes screened.
Demonstration (run): a replicated SBP gene whose 5 trade-off rows are all beta 0, p 1, se NaN is Tier A. With real Genebass rare-outcome analyses (schizophrenia, dementia, infertility) a gene can plausibly reach 5 "screened" outcomes with several such rows.
Fix: count an outcome as screened only if its row has finite se above 0 (apply stats.informative in tradeoff_tested); keep adverse detection unfiltered (see R3-3). Ledger decision because C1 does not define "screened".
Test: test_uninformative_tradeoff_rows_do_not_count_as_screened (xfail).

### R3-2 [Low-Medium] The p = 1 filter is exact; uninformative rows stored just below 1 still mask lambda_GC
File: src/protscan/controls.py:26-32. C11b tests p at least 1. If the source writes uninformative results as 0.99999 (or 1 minus rounding), the gate sees informative rows. Demonstration (run): 45 percent of synonymous rows at p 0.99999 and the rest inflated 1.5 times (no tail hits): status OK; with exact p 1 the same table is NOT RUN. Under the null only 1 percent of rows have p above 0.99.
Fix: also report the share of synonymous rows with p above 0.99 and treat more than 5 percent as not evaluable (expected share 1 percent), or compute lambda_GC on rows with p below 0.5 plus a check. Ledger decision.
Test: test_lambda_control_is_not_masked_by_p_just_below_one (xfail), paired with a passing scope test at exact p = 1.

### R3-3 [Info] Adverse screen and discovery hits are deliberately or adapter-guarded, not informativeness-filtered
adverse_tradeoffs (stats.py:80-96) still counts a harmful row with undefined se: conservative (harm is never hidden). Tested as documented behaviour. discovery_hits does not apply informative(): safe today because both adapters guarantee that an undefined se comes only with p = 1, but a future adapter that emits an undefined-se row with a small p would create a discovery hit. Suggest applying informative() in discovery_hits as defence in depth.

### R3-4 [Low] Adapter drop counts are still outside the result record (R2-7 open)
scripts/fetch_finngen.py and scripts/fetch_genebass.py write finngen_conversion_summary.csv and genebass_conversion_summary.csv to data/ (gitignored) and print totals; run.py and report.py read neither. docs/data-sources.md records "0 dropped" for the extracted tables; that statement is not machine-checked at run time. If a re-fetch drops rows, the results JSON and report will not show it.
Fix: have run_pipeline read the two summary files when present and copy their totals into results data, and print them in the report.

## 4. Adversarial checks that found nothing
- Sign-control gaming: an arm needs an informative row (finite se above 0, beta nonzero) and p below 0.05, strict. A right-signed underpowered arm is NOT RUN, so noise cannot pass the control; a wrong-signed underpowered arm cannot fail it; a wrong-signed powered arm fails it even when the other arm is right; a null table passes with probability 1/4 to 1/2 per evaluable arm only when the arm reaches p below 0.05 (about 2.5 percent per arm for a wrong sign by chance). Untrusted-source rows and undefined-se rows are excluded. There is no adapter branch that adjusts sign.
- Lambda with few informative rows: the 10,000-row floor and the 5 percent p = 1 gate leave at least 9,500 rows with p below 1, so the median is never taken over a handful of rows. p = 1 fraction boundary: exactly 5 percent evaluable, one more row not evaluable (unit and 10,000-row exact table).
- se above 0 rule versus adverse screen: see R3-3. Versus tiering: a degenerate replication row leaves the trait as trait_missing (Tier B), never A or D; a row for another gene does not hide it.
- Reason strings: sign_control_failed and sign_control_not_evaluable reach the verdict text, and the KILL reading rule is printed in the report.
- Adapter rules: drop rule applies only to undefined se with p below 1 (p = 1 rows kept), on every FinnGen endpoint and every Genebass analysis and mask; A1FREQ guard unchanged (property test by the builder plus the reviewer tests).
- Positive controls: PCSK9 LDL still requires the discovery threshold, the right direction and presence in the pipeline discovery-hit path; trade-off harm on PCSK9 (discovery or replication cohort) is reported as tier C without changing status.

## 5. Test evidence (run)
- tests/test_ledger_conformance.py: 284 passed, 3 xfailed (R2-1, R3-1, R3-2). Full repository suite at HEAD 866580e plus reviewer tests: 479 passed, 3 xfailed. No test was skipped or weakened; 3 xfail markers (R2-2 twice, R2-4) were removed only after each test was shown to fail with its fix removed.
- New closure tests: C11a (uninformative replication rows, adapter drop rules with counts, ivw_combine, lambda and coverage still use them), C11b (5 percent boundary on a 10,000-row table and on the null table, both lambdas, hidden inflation with 4 percent p = 1, end-to-end KILL at 6 percent), C11c (p boundary 0.0499 / 0.05 / 0.0501, powered wrong sign, both underpowered, reason strings in the verdict, uninformative rows cannot be an arm, end-to-end), C11d (Tier C PCSK9 leaves the control OK, discovery-effect failures still fail it, end-to-end PASS with PCSK9 type 2 diabetes harm, control still requires the pipeline discovery hit).
- Oracle: independent re-implementation now includes C11a (informative replication rows), C11b (lambda on p below 1 and p = 1 share), C11c (powered arms), C11d (any tier). It agrees with the pipeline on 8 synthetic scenarios via python -m protscan run, on 25 fuzz tables that now include undefined-se and beta-zero replication rows, and on the mutation cases.
- Mutants of the new code (21): informative() disabled or weakened (2), sign control and replication_status using uninformative rows, p = 1 gate off / at-or-above / counted at 0.99, lambda on all rows or on p at most 1, sign power gate off or inclusive, underpowered arm as failure, reason strings swapped, positive control tier-C exclusion restored, positive control without the hit path, both adapter drop rules off, drop rule extended to p = 1, ivw_combine keeping p below 1, stale pin. One survived at first (positive control without the hit path, an equivalent-looking expression); a test was added and the mutant is now killed. All 21 killed.
