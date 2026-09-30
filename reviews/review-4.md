# Review 4: Amendment 2 (post-hoc after run 1), diff-scoped re-review of fbc0053

Reviewer: independent (did not write the pipeline). Scope: git diff 2b373c3..HEAD of src/ and scripts/ (config, controls, stats, tiering, run, report, schema, make_synthetic).
Ledger basis: experiments.md RUN 1 RESULT and AMENDMENT 2 (A2a-A2d). Not opened: results/run1-protective-scan.json, results/run1-report.md, any file under data/. No real outcome was viewed; no network access.
Framing kept from the ledger: Amendment 2 changes the negative control after run 1 failed it. Run 1 stands as KILL under the original rule. Anything reported from run 2 is "passed under Amendment 2, which was written after run 1 failed".

## 1. Pin and config (independent)
- shasum -a 256 config/prereg.yaml gives be69827539daac4c78b4c8f6ac7d931723abb41ab0e739c0505c357689310a9f, equal to schema.PINNED_CONFIG_SHA256 and to the lead message.
- git diff 2b373c3 HEAD -- config/prereg.yaml: syn_hits_max removed; syn_contaminated_max_fraction 0.001 added, with a comment. Nothing else changed (thresholds, panel, directions, proxies, coverage and row floors, p = 1 share, sign-arm p, lipid list, min screened).
- Reviewer literal updated. A2 constants typed from the ledger text: 0.001; either direction; any trait; excluded from all tiers and candidate lists; label "Amendment 2 (post-hoc after run 1)".

## 2. Conformance of fbc0053 to A2a / A2b (verified by tests that fail when the code is mutated)
- A2a: stats.contaminated_genes takes discovery-cohort synonymous rows with p below the discovery threshold (strict), any trait (13 panel, 9 trade-off, triglycerides), either sign. Plof, dmis, replication and EUR rows never contaminate. A gene with hits on several traits is one gene.
- A2a exclusion: tiering.build_hits removes contaminated genes before replication, adverse and tier assignment, so a contaminated gene appears in none of Tier A to D; tested with one contaminated gene per tier next to a clean twin, and through decide(), the candidate rung, verdict lists and the report.
- A2b: gate is lambda below 1.10, coverage, row floor, p = 1 share (unchanged) and contaminated genes at most 0.1 percent of genes tested. Denominator is the number of distinct genes with a discovery synonymous row (not rows, not plof genes). Boundary tested at 1/1000, 2/2000, 3/3000 (pass) and 2/1000, 3/2000, 4/3000 (fail), with float equality safe because k/(1000k) rounds to the same double as 0.001; a gene with three hit traits counts once; a rows denominator or a plof-gene denominator is caught.
- Disclosure: verdict block (amendment field), label, CLI stdout, report ("Rules active" line and "never a clean preregistered pass"), and results JSON contamination.amendment all carry "Amendment 2 (post-hoc after run 1)" (tested on five scenarios and the CLI).
- Independent oracle updated (filter, exact integer bound, positive controls judged on the unfiltered discovery hit). It agrees with the pipeline on 11 synthetic scenarios (three new: contaminated_tierA gives LEAD, contaminated_systemic gives KILL, contaminated_minor gives PASS), on 25 fuzz tables that now include random synonymous hits, and on targeted mutations.

## 3. Adversarial questions from the coordinator
- Can the filter drop a genuine gene silently? Not silently: every excluded gene is in the report table and in results JSON contamination.genes with its synonymous traits and the pLoF traits that were excluded; verdict.contaminated_excluded lists them; the simplest rung is unfiltered and labelled so. Two gaps remain (R4-1, R4-2 below). The chance of a false exclusion is negligible (about 0.06 expected synonymous hits at 1.9e-7 over 334,000 rows, per the ledger).
- Is the denominator right? Yes: genes tested (distinct genes with a synonymous discovery row), consistent between negative_control and the contamination record. Genes whose synonymous rows are all p = 1 still count as tested (slightly lenient, under 1 percent on run 1 counts).
- Can a contaminated positive-control gene corrupt control status? No. Positive controls are judged on the unfiltered stats.discovery_hits, so a contaminated PCSK9, APOC3 or ANGPTL4 still passes its control; the replication-sign control reads the replication cohort only, so a contaminated LDLR does not touch it. Such genes are excluded from the candidate tiers and counted toward the A2b total (tested). Mutant that feeds the filtered hits into the control is killed.
- Does contamination in replication-only data matter? Not detectable: FinnGen has no synonymous mask, and synonymous rows in replication or EUR cohorts are ignored by A2a (tested). A locus that contaminates only the FinnGen pLoF call cannot be seen by A2a; A2c (manual GWAS Catalog lookup) is the only defence, and it is not automated (R4-4).
- Can A2 be bypassed or over-triggered? Bypass: not through config (pin stamps and CLI exit 2), not through cohort, mask, direction or trait choice. It can be evaded by data shape: a pLoF hit gene with no synonymous row is unassessable and stays in the tiers (R4-1), and a gene with a strong but sub-threshold synonymous signal (for example p 1e-6) stays (by design, ledger threshold). Over-trigger: one trait with unreliable p-values (a rare binary outcome) can add several contaminated genes and trip A2b; the report lists the traits per gene so this can be read (reading rule in section 5).

## 4. Findings
### R4-1 [Low-Medium] A pLoF hit gene with no synonymous row is unassessed and not flagged
Files: src/protscan/stats.py:59-67 (contaminated_genes), src/protscan/tiering.py:25-27. A2a can only exclude genes that have a synonymous row. Coverage is gated at 90 percent of plof pairs, so up to 10 percent of pairs (including hit pairs) may lack a synonymous row; those genes stay in the tiers with no indication. Run 1 coverage was 99.9 percent, so the exposure on real data is small, but the flag is free.
Fix: add syn_assessed and min_syn_p columns to each tier record and print them; require syn_assessed for PASS/LEAD qualification or list unassessed candidates prominently.
Test: test_a2a_hit_gene_without_any_synonymous_row_is_flagged_as_unassessed (xfail).

### R4-2 [Low] A verdict changed by the filter is not said in the verdict reason
File: src/protscan/run.py decide/reason. With SYNPASS1 contaminated the verdict becomes LEAD ("no qualifying Tier-A gene ...") and the reason text does not say that a would-be Tier A gene was excluded. The report and JSON do list it. Fix: append "N candidate gene(s) excluded by A2a" to the reason when the exclusion removed a would-be PASS or LEAD gene.

### R4-3 [Low] Header still says "with AMENDMENT 1"
The config ledger_entry string (config line 3) is unchanged, so results.ledger_entry and the report header line read "2026-09-29 protective-variant-scan (with AMENDMENT 1)". The rules-active line and label carry Amendment 2, so the disclosure is present, but the header is stale. Fix with the next dated config change. Test: test_a2_ledger_entry_header_names_amendment_2 (xfail).

### R4-4 [Info] A2c needs coordinates the pipeline does not carry
A2c requires chromosome position for every PASS or LEAD gene. The Genebass records contain chrom and pos, but normalize_analysis drops them and the constraint table has no coordinates, so the manual step has nothing to read from the result. Suggest an annotation file (gene, chrom, pos) written by the fetch script and added to the tier records; not part of the decision path.

### R4-5 [Info] The 0.1 percent bound is six times the run-1 count
About 18 genes of 18,500 are allowed; run 1 had 3. This is disclosed as post-hoc in the ledger and in the report. No code issue; the disclosure text must accompany any PASS or LEAD.

### Carried open items (unchanged by this diff)
R2-1 (hit pairs inside the 90 percent coverage), R3-1 (uninformative trade-off rows count as screened), R3-2 (p just below 1), R3-4 (adapter drop counts outside the result record), R2-6, R2-8 to R2-10. See review-2 and review-3.

## 5. Reading rule for run 2
- A KILL from negative_synonymous must be read with the contamination table: check whether the excess comes from one trait, one locus or many loci before calling it systemic.
- A PASS or LEAD is "passed under Amendment 2, which was written after run 1 failed". A2c must be completed by hand for each such gene.
- The trivial rung still counts only beneficial-direction panel-trait synonymous hits; the contamination list is the wider set (either direction, any trait).

## 6. Test evidence (run)
- tests/test_ledger_conformance.py: 368 passed, 5 xfailed (R2-1, R3-1, R3-2, R4-1, R4-3). Full repository suite before the last added test: 585 passed, 5 xfailed. No test was weakened. The 9 tests that failed after fbc0053 were updated with intent kept: config controls and pin (new key and hash), config load test, trivial rung (rebuilt on a 2,100-gene frame: beneficial hit counted, harmful hit only in the contamination list), the beneficial-direction syn-hit control test (still FAIL on a 500-gene table, now on the contamination bound with lambda OK), the harmful-direction test (now shows the A2a exclusion and no kill within the bound), the systemic break in the control-break parametrization (three contaminated null genes so the Tier A gene survives), the two config-edit tests (edit syn_max_p1_fraction instead of the removed key), and the hit-path test (monkeypatch of stats.discovery_hits).
- A2 mutants: see section 7.

## 7. Mutants of the A2 code (18, all killed)
Contamination threshold made inclusive; contamination limited to panel traits; limited to one direction; bound denominator changed to rows or to plof genes; bound made strict or always true; lambda dropped from the control status; old zero-hit gate reinstated; contamination read from any cohort; build_hits filter removed; positive control fed the filtered hits; amendment label dropped from the verdict label; rules-active line renamed; verdict.contaminated_excluded emptied; excluded pLoF traits omitted from the record; within_bound forced true; stale pin. No survivor.
