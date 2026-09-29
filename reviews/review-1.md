# Review 1: protective-variant scan, conformance to experiments.md (entry 2026-09-29 + Amendment 1 + C1-C9)

Reviewer: independent (did not write the pipeline). Reviewed at HEAD ff04c22 plus uncommitted test additions.
Tests: tests/test_ledger_conformance.py (constants typed from the ledger, not read from config).

## Verification status

Ran (code executed, results observed):
- 202 conformance tests: 199 pass, 3 xfail (documented gaps R-1, R-2, R-3 below). Full repo suite at ff04c22: 316 passed, 2 xfailed (before the R-2 xfail was added).
- Independent oracle (re-implementation of discovery, replication, adverse, tier, C1/C9, verdict from the ledger text) compared with the pipeline on all 6 synthetic scenarios via python -m protscan run: tiers and verdicts identical. 25 fuzz tables (boundary p-values, UKB/unknown sources, partial trade-off screens): tiers identical.
- Exhaustive check of run.decide over all 256 combinations of (tier, qualifies): PASS iff a qualifying Tier A exists; LEAD iff no such gene and a qualifying Tier B exists; KILL otherwise; any failed or not-run control gives KILL. A LEAD cannot be reported as PASS.
- Mutation check: 24 hand-made mutants of the pipeline (direction flips, threshold changes, sidedness, tier swap, control bypass, LEAD to PASS, lipid/domain/mask filters, proxy cross-wiring, lambda df) were run against the conformance file; all 24 killed.
- Edge probes (script only, not committed): empty table, replication-only table, NaN p on PCSK9, lowercase genes, p=0/se=0 all end in KILL or a normal verdict without a crash.

Not verified / not done:
- No real data was opened (C8 hold). Sign conventions of the real Genebass and FinnGen tables are unchecked; see R-2.
- Not reviewed line by line: adapters/gnomad.py, adapters/coverage.py, scripts/fetch_gnomad.py, scripts/fetch_coverage.py, builders' own tests (tests/test_pipeline_*.py, tests/test_adapters_*.py), fixtures. adapters/common.py, adapters/finngen.py, adapters/genebass.py, scripts/fetch_genebass.py, scripts/fetch_finngen.py were read fully.
- Network fetch scripts were not executed.

## Findings

### R-1 [Medium] Negative control passes vacuously on a tiny or selected synonymous table
File: src/protscan/controls.py:12-25 (negative_control). The only guard is syn.empty. Three synonymous rows give lambda_GC from 3 p-values and zero hits, so controls.valid is True.
Scenario: partial Genebass fetch, or the adapter drops most synonymous rows (see R-9), or synonymous rows exist only for non-hit genes. Probe at HEAD ff04c22: PCSK9/APOC3 control rows + 3 syn rows + one cognitive hit gene with 5 trade-off rows returned LEAD, controls valid True, n_rows 3.
Fix: require the synonymous universe to cover the discovery pLoF universe (for example at least 90 percent of (gene, trait) pairs with a plof row also have a syn row, plus a fixed minimum row count); otherwise status NOT RUN (already maps to KILL controls_not_evaluable). Report the coverage ratio. The ledger does not state this rule, so add it by dated entry.
Test: test_negative_control_must_cover_the_discovery_universe (xfail).

### R-2 [Medium] Nothing validates the replication-side effect-sign convention
Files: src/protscan/controls.py:28-47 reads only the discovery cohort; src/protscan/adapters/finngen.py:66-79 uses BETA as is; A1FREQ is used only for a carrier count.
Scenario: the FinnGen BETA refers to the non-carrier allele, or the sign is otherwise reversed. Probe: flipping every replication beta in the synthetic pass scenario turns PASS into LEAD, controls.valid stays True, and PCSK9 silently becomes Tier D. A sign reversal gives false negatives (true replications fail) and false replications of genes whose true FinnGen effect is harmful. All Tier A/D decisions and the adverse screen depend on that sign.
Fix: (a) adapter guard: assert A1FREQ <= 0.5 for every kept row (LoF carrier allele is A1), drop or flip otherwise; (b) add a diagnostic (report-only, or gating via a dated entry): PCSK9 to hypercholesterolemia must be in the beneficial direction in the replication cohort, listed under Controls. Because C8 holds FinnGen tables closed, this cannot be checked yet.
Test: test_replication_sign_flip_is_detected_by_a_control_or_flag (xfail).

### R-3 [Medium] Any --config is accepted; the verdict is not tied to the preregistered constants
Files: src/protscan/run.py:106-107, src/protscan/schema.py:67-121. The report prints the first 16 hex chars of the config sha256 but nothing compares it with the ledger values. Probe: syn_hits_max set to 3 in a copied config with the broken_syn_hit scenario returns PASS. load_config also tolerates p_threshold within 2 percent of the formula and any n_genes.
Fix: pin the expected sha256 of the committed prereg.yaml (or the dict of ledger constants) in code, and stamp config_matches_ledger false plus a NON-PREREGISTERED line in the report (or refuse) on mismatch. The conformance test is currently the only guard and the CLI does not run it.
Test: test_edited_config_cannot_silently_produce_a_verdict (xfail).

### R-4 [Medium] requests is imported but not declared
Files: src/protscan/adapters/common.py:33, scripts/fetch_genebass.py:18, scripts/fetch_finngen.py:16; pyproject.toml:9 lists only pandas, numpy, scipy, pyyaml. In a clean venv pytest fails tests/test_adapters_common.py::test_download_streams_and_checksums_over_localhost (it passed here only because requests 2.34.2 was installed into the shared .venv). V1.md "Done when" requires pytest to pass.
Fix: add requests to dependencies (or an adapters extra that the test skips without).

### R-5 [was High, FIXED and verified] UK-Biobank-overlapping data could be counted as independent replication
Original defect (commit 3fd92e9): stats.replication_status used every row with cohort equal to replication regardless of source; a Genebass/AZ/Regeneron-derived table labelled replication gave Tier A. Fixed in 1aff7c0 (deny list) and then allow list plus deny list (src/protscan/stats.py:37-52; config/prereg.yaml replication.independent_sources and ukb_overlapping_sources). Verified by 7 named-source tests, 6 name-variant tests (uk-biobank, az_phewas, opentargets, mystery source), 5 independent-source tests, discovery and discovery_eur cohort rows, and the mutation that empties the deny list.
Residual (Low): independence still rests on the free-text source string set by the adapter and on the adapter choosing cohort replication. Substring matching of the allow token means a source called finngen_plus_something is accepted unless it carries a deny token. Keep adapters as the only writers of source.

### R-6 [Low] Adverse screen and coverage count rows from sources the replication filter rejects
File: src/protscan/stats.py:80-102 (adverse_tradeoffs, tradeoff_tested) use every replication-cohort row, including rows from unknown or UKB-overlapping sources. For adverse detection this is conservative. For the C1/C9 minimum of 5 screened outcomes it is lenient: untrusted rows can supply the screening. Fix: count screening only from discovery plus independent replication rows, or state the current behaviour in the ledger.

### R-7 [Low] Trade-off screen structural blind spots (ledger-conformant, worth stating in any report)
- all_cause_mortality is absent from Genebass and FinnGen; fracture is absent from FinnGen (docs/data-sources.md 3.2). With C1 (at least 5 of 9) a no-adverse-trade-off claim covers at most 8 of 9 outcomes.
- Adverse is defined only on the 9 outcomes. A gene that is beneficial for SBP and significantly harmful for a panel trait (FEV1, fluid intelligence, grip) is still Tier A. Add an informational column for harmful-direction panel hits (report-only).

### R-8 [Low] Report omits caveats the ledger requires with any PASS or LEAD
C5 says the dmis = missense|LC limitation (direction-only, weaker than intended) is reported with any LEAD or PASS; C7 says cognitive coverage is fluid intelligence and reaction time only. src/protscan/report.py:86-95 (Notes section) prints neither. The synthetic flag appears only as a report banner (report.py:30-31); the verdict block in the JSON and the CLI stdout line (__main__.py:17) carry no synthetic marker, so a copied verdict line loses it. Fix: add both notes when verdict is PASS or LEAD and data is not synthetic; add a synthetic flag to the verdict block and stdout.

### R-9 [Low] Adapter drops rows with beta equal to 0 or p equal to 1, biasing lambda_GC upward
Files: src/protscan/adapters/common.py:80-86 (se_from_beta_p returns NaN when beta is 0, or z is 0 for p of 1) and genebass.py:86 (dropna on se). Uninformative high-p rows are removed from the synonymous table before lambda_GC is computed. Direction of error: lambda inflated, so a spurious KILL is possible; a false PASS is not. Fix: do not drop; write a small positive placeholder se or make se optional in schema.validate, and log the count of dropped rows per source.

### R-10 [Low] Positive controls re-implement the rule instead of exercising the pipeline path
File: src/protscan/controls.py:38 recomputes the direction and threshold test. A bug in stats.discovery_hits or tiering would not be caught by the controls (here it is caught only by the conformance tests). Fix: additionally assert that PCSK9/LDL is returned by stats.discovery_hits and is not Tier C.

### R-11 [Low] Mixed synthetic and real inputs
File: src/protscan/schema.py:156-165 loads every top-level burden*.csv[.gz]. A leftover synthetic file next to real tables is merged; only the banner marks it (run.py:139). Duplicate PCSK9/APOC3/ANGPTL4 keys currently raise, which is accidental protection. Fix: refuse to run when synthetic and non-synthetic sources coexist.

### R-12 [Info] Process
- C8: scripts/fetch_finngen.py downloads and then normalizes the gene-level FinnGen table in one step (lines 45-69). The ledger holds FinnGen tables closed until the user submits the request form. Do not run it before that.
- Real Genebass synonymous lambda_GC may exceed 1.10 because the lambda-based QC flags are deliberately not applied (C6). If it does, the ledger outcome is KILL control_failed (a data property, not a pipeline bug). Read controls.negative_synonymous before interpreting any hit list.
- PASS is reachable only via SBP + FinnGen hypertension (C2). Most SBP hits will be gene_untested (FinnGen tests about 4.9k genes) and fall to Tier B, i.e. LEAD at best (C4).

## Items checked and found correct (no finding)
- Direction of benefit for all 13 panel traits and all 9 trade-off outcomes; one-sided replication boundary (two-sided 0.0999 replicates, 0.10 does not); wrong-direction replication never replicates; per-gene matching; proxy cross-wiring rejected; traits without a declared proxy capped at B even with a strong same-trait replication row.
- Discovery filter: plof only, discovery cohort only, panel traits only, strict inequality, beneficial sign; synonymous, dmis, EUR and replication rows never create tiers.
- Adverse: 9 outcomes, harmful is beta above 0, one-sided alpha/9 (C4); protective direction never adverse; adverse wins over failed replication (C before D; ledger silent).
- C1/C9 minimum of 5 screened outcomes: union across cohorts, plof rows only, boundary 4/5, applies to the Tier A cap and to LEAD.
- Controls: each break (PCSK9 LDL sign, LDL below threshold, missing gene, CAD sign, both TG genes, lambda 1.101, beneficial synonymous hit) gives KILL even when a Tier-A gene exists; lambda 1.099 does not; a harmful-direction synonymous hit does not count; a control that cannot run gives KILL controls_not_evaluable (C4).
- Config equals the literal ledger constants (panel, directions, domains, 9 trade-offs, proxies, thresholds, lipid list frozen at 4073b89, min screened 5, lambda 1.10, syn hits 0). No permutation run remains (Amendment 1).
