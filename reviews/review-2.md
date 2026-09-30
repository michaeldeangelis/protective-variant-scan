# Review 2: re-review of C10 fixes (commits 6b959b6 pipeline, 51b8e05 data) and previously skipped files

Reviewer: independent (did not write the pipeline). Reviewed at HEAD 6b959b6 plus reviewer test commits up to 210d3db.
Ledger basis: experiments.md, entry 2026-09-29 + Amendment 1 + C1-C10. No real data outcome was opened; nothing was fetched from the network.

## 1. Closure of review-1 findings (each checked by a test that fails when the fix is removed)

- R-1 vacuous negative control (C10a): CLOSED. Tests: three synonymous rows now NOT RUN and KILL controls_not_evaluable; row floor 434 vs 435 genes; coverage boundary at exactly 90 percent of plof pairs; coverage measured against plof pairs not row count; no plof rows. Mutants (gate AND instead of OR, coverage any(), inner merge) all killed.
- R-2 replication sign unvalidated (C10b): CLOSED for the sign convention, residual R2-3. Tests: three wrong-sign combinations, one-arm-evaluable rules, none evaluable, untrusted sources, end-to-end flipped replication now KILL control_failed (was LEAD with controls valid). A1FREQ guard: 0.5 kept; 0.5001, 0.999 and NaN dropped, never flipped; applied on every endpoint built.
- R-3 any config accepted (C10c): CLOSED, residual R2-6. Edited config stamped in result, label, JSON on disk, report and CLI exit code 2; trailing newline, comment line, extra key and CRLF all change the pin; identical bytes at another path match; reviewer holds an independent literal copy of the file sha256.
- R-4 requests undeclared (C10g): CLOSED by text check (a clean-venv install was not run). pyproject lists requests and pyarrow.
- R-6 untrusted rows supply screening (C10g): CLOSED. Screening counts only allow-listed rows; adverse detection still uses every row (conservative, tested).
- R-7 harmful panel hits invisible (C10e): CLOSED. Column lists significant harmful-direction hits on other panel traits, excludes the row own trait, does not gate.
- R-8 caveats and synthetic flag (C10f): CLOSED. C5 and C7 text in PASS and LEAD reports; synthetic and preregistration flags in verdict block, label and CLI stdout.
- R-9 rows with p equal to 1 or beta equal to 0 dropped (C10d): CLOSED, but see R2-4 (new consequence).
- R-10 controls re-implement the rule (C10g): CLOSED, but see R2-5. PCSK9 LDL must come out of build_hits and not be Tier C; mutant removing the Tier C exclusion killed.
- R-11 mixed synthetic and real inputs (C10g): CLOSED. load_burden and run_pipeline raise ValueError; mutant killed.

## 2. New findings

### R2-1 [Low] Coverage gate allows selective omission of exactly the hit genes
File: src/protscan/controls.py:10-38 (negative_control). 90 percent coverage still permits the omitted 10 percent to be the discovery-hit (gene, trait) pairs, which are the pairs the negative control exists to police. A synonymous table that lacks the rows of hit genes passes.
Fix: also report and gate the fraction of discovery-hit pairs that have a synonymous row (require all, or at least 90 percent), and list the uncovered hit pairs.
Test: test_c10a_residual_selective_omission_of_hit_genes_is_detected (xfail).

### R2-2 [Medium] FinnGen rows with SE at most 0 or missing keep their LOG10P-derived p and can reach Tier A
Files: src/protscan/adapters/finngen.py:101 (se becomes NaN, row kept), src/protscan/adapters/common.py:95-125 (ivw_combine keeps such rows), src/protscan/stats.py replication_status (uses beta and p only).
The C10d rationale is that an undefined se means beta is 0 or p is 1 (uninformative). For Genebass that holds because se is derived from beta and p. For FinnGen se is reported by regenie; SE of 0 or NaN with a large BETA and a significant LOG10P is a degenerate or separated fit. Previously dropped, now kept and used.
Probe (run): synthetic pass scenario, the SYNFAIL1 replication row (a failed replication, p 0.62) replaced by beta -0.3, p 1e-9, se NaN: SYNFAIL1 becomes Tier A and the verdict stays PASS. A unit case at tier level shows the same.
Fix: keep an undefined-se row only if beta is 0 or p is at least 1 (uninformative); otherwise drop it and count it in the conversion summary. Apply the same rule in ivw_combine.
Tests: test_c10d_finngen_row_with_degenerate_se_but_significant_p_is_not_used and test_c10d_replication_row_with_undefined_se_and_tiny_p_cannot_give_tier_A (xfail).

### R2-3 [Medium] Replication-sign control is direction-only, every evaluable arm must hold, and it covers one endpoint
File: src/protscan/controls.py:73-95. The whole test is beta times expected sign above 0.
- No significance requirement: a null table passes with probability 1/4 to 1/2. Conversely a genuine but low-powered arm can fail by chance. PCSK9 pLoF carriers in FinnGen are few and the effect on the hypercholesterolemia endpoint is modest, so the PCSK9 arm has a real chance of a wrong sign from noise; C10b then gives KILL control_failed although the pipeline is fine (a false KILL, not a false PASS). LDLR should be strong. Power could not be checked because no outcome was opened.
- The control validates the sign convention of the hypercholesterolemia endpoint only. PASS depends on the hypertension endpoint. A single regenie table makes a per-endpoint flip implausible, but nothing tests it.
Fix (ledger decision for the lead): count an arm as evaluable only when its one-sided p in either direction is below 0.05 (a significant wrong sign is FAIL, a non-significant arm is NOT RUN), keep the at-least-one rule.
Test: test_c10b_sign_is_direction_only_so_noise_can_pass_documented (passes; documents the scope).

### R2-4 [Medium] Keeping p = 1 rows lets uninformative rows deflate lambda_GC and mask inflation
Files: src/protscan/controls.py (lambda over all synonymous rows), src/protscan/adapters/genebass.py (rows with beta 0 or p 1 now kept). lambda_GC is a median: with 45 percent of synonymous rows at p = 1, informative rows with lambda about 2 give a pooled lambda about 0.05 and the control passes.
Demonstration (run): 45 to 50 percent p = 1 rows, informative rows inflated 1.3 to 1.5 times, no tail hits: pooled lambda 0.045 to 0.052, status OK, while lambda on the informative rows is 1.7 to 2.3. The ledger text (lambda over all genes and traits) does not forbid this and C10d chose to keep the rows, so this is a design consequence, not a code defect. It reverses the direction of the R-9 bias (was a spurious KILL, now a spurious pass).
Fix: report the fraction of synonymous rows with p = 1 and lambda_GC over defined-se rows next to the pooled value, and gate on the larger of the two (dated ledger entry).
Test: test_c10d_lambda_control_is_not_masked_by_uninformative_rows (xfail).

### R2-5 [Low] The pipeline-path positive control can fail on genuine biology
File: src/protscan/controls.py:60. PCSK9 LDL must not be Tier C. PCSK9 loss of function is reported to raise type 2 diabetes risk modestly; a one-sided p below 0.0056 in Genebass or FinnGen makes PCSK9 Tier C, the positive control fails and the verdict is KILL. Decided in C10g; flagged so that a KILL for this reason is read as a control-design effect (see the tier field in the report), not as pipeline invalidity.

### R2-6 [Low] The config pin is self-referential
Files: src/protscan/schema.py:19, tests/test_pipeline_c10.py test_pinned_hash_matches_committed_config. The pin, the config and the builder test live in one repository and change together; the builder test re-derives the hash from the file. The reviewer test carries a second literal copy (L_CONFIG_SHA256) and must be updated with any dated ledger entry that changes the config. On a mismatch verdict["verdict"] still reads PASS; only the label, the flag and the exit code change, so a consumer reading only that field loses the stamp. A mismatching run still writes results/protective-scan.json.

### R2-7 [Low] A1FREQ drops are not part of the result record
Files: src/protscan/adapters/finngen.py:73-80, scripts/fetch_finngen.py. Dropped and undefined-se counts go to a log and to data/finngen_conversion_summary.csv, not to the results JSON or the report. If the guard drops many rows the run cannot show it. The guard tests a frequency, so it detects an allele-coding flip but not other sign errors (see R2-3).
Fix: copy the per-endpoint summary into the data section of the results JSON and print it in the report.

### R2-8 [Low] KILL reports still print a PASS-qualifying gene list
File: src/protscan/report.py:45-47. Under KILL (control failed) the report prints "PASS-qualifying genes (Tier A): SYNPASS1"; tests/test_pipeline_run.py test_broken_controls_kill asserts that list is populated under KILL. The verdict is correct; the list invites misreading. Print it only when the verdict is PASS or LEAD, or label it "would qualify if controls were valid".

### R2-9 [Low] Builders' tests: config-relative and count-only assertions
- tests/test_pipeline_verdict.py:22-23 (_gene_rows) and scripts/make_synthetic.py (good and bad helpers) derive direction from cfg.panel_sign and cfg.sign, so a wrong sign in the config is self-consistent in those tests and in every synthetic scenario. Independent direction checks exist only in tests/test_pipeline_schema.py test_direction_of_benefit (literal set) and in the reviewer file.
- tests/test_pipeline_c10.py test_pinned_hash_matches_committed_config is a consistency check, not an independent pin (R2-6).
- tests/test_adapters_finngen.py test_kept_rows_never_exceed_half_allele_frequency asserts a row count (5), not the property in its name.
- tests/test_pipeline_stats.py test_adverse_threshold_direction_and_mask and THR use cfg.tradeoff_p and cfg.discovery_p, so a drifted config value moves test and code together.
- No skips or xfails in the builders' files (grep). One monkeypatch (test_positive_control_not_in_pipeline_hits_fails) is legitimate.

### R2-10 [Low] Previously skipped files
- src/protscan/adapters/gnomad.py:25 drops rows with missing LOEUF before choosing one transcript per gene. A gene whose MANE transcript has no LOEUF silently takes a canonical or other row, so the transcript rule is not applied uniformly. Annotation only; no tier depends on it. Fix: choose the transcript first, then drop missing LOEUF.
- scripts/fetch_gnomad.py: unlike fetch_finngen.py it does not verify the GCS md5 of the download. Annotation only.
- src/protscan/adapters/coverage.py and scripts/fetch_coverage.py: count rows per trait, mask and cohort; no outcome values are read beyond counts. A trait in neither genebass.TRAITS nor genebass.ABSENT raises KeyError (loud). replication_genes_plof counts rows, which equals genes because keys are unique. No defect.
- tests/test_adapters_gnomad.py and tests/test_adapters_coverage.py assert concrete numbers on fixtures; not tautological.

### R2-11 [Info]
- synthetic_replication is an allow-listed source token in the real config (needed by the fixtures). The refusal to mix synthetic and real sources prevents it from co-existing with real rows, so it cannot smuggle rows into a real run.
- The C8 hold is lifted in the ledger. scripts/fetch_finngen.py and scripts/fetch_genebass.py make network calls and were not executed by the reviewer.

## 3. Adversarial checks that found nothing
- Pin bypass: byte-level sha256 of the file. Cosmetic edits, CRLF, extra keys and comment lines all mismatch; identical bytes elsewhere match; mutants that always match or return exit code 0 are killed. Editing schema.py and the config together is not detectable in code (R2-6) but is by the reviewer literal.
- Sign-control gaming: beta of 0 fails; duplicate keys raise in validate; untrusted-source rows do not count; gene case handled; only independent replication rows are read. The adapter has no branch that adjusts sign to satisfy the control.
- Coverage edge cases: exactly 90 percent passes and one pair less fails; exactly 10,000 rows passes and 9,999 fails; no plof rows gives NOT RUN; synonymous rows for other pairs do not raise coverage; rows with se NaN count.
- NaN se: validate keeps them; results JSON serializes them as null; a hit record with NaN se serializes (probe run).
- A1FREQ: inclusive at 0.5, NaN dropped, applied in every build_trait path.

## 4. Test evidence
- tests/test_ledger_conformance.py at 210d3db: 250 passed, 4 xfailed (R2-1, R2-2 twice, R2-4). Full repository suite: 407 passed, 4 xfailed.
- Mutation runs against the C10 code: 20 new mutants (coverage gate, sign control, pin, CLI exit, mixing refusal, A1FREQ, screening rows, tier-C exclusion, harmful panel, caveats, NaN se, decide) all killed.
