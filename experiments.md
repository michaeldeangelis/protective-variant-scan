# Experiments ledger (append-only; corrections are new lines under the entry)

## 2026-09-29 protective-variant-scan   (rules fixed before the run)

Question: Among human genes, does carrying a predicted loss-of-function (pLoF) variant associate with a BETTER outcome on a fixed cognitive/physical trait panel, replicate in an independent cohort, and show no significant adverse trade-off?
Decision it drives: PASS -> take the Tier-A genes into cell-type and regulatory follow-up (threads 1 and 4). KILL -> write a negative result and make thread 1 (cognition cell-type mapping) the primary project.
Held fixed: public gene-based burden results only (no individual-level data); trait panel and direction of benefit below; thresholds below; discovery/replication cohort assignment (set in Step 0 feasibility, before any outcome is viewed).
Varied: only the pipeline rung (see ladder).

Sources (candidates, confirmed or replaced in Step 0 feasibility; any replacement is a new dated entry):
- Discovery: UK Biobank exome gene-based burden results (e.g. Genebass; AstraZeneca PheWAS portal).
- Replication: an independent cohort with gene-based LoF results (e.g. FinnGen, other biobank burden portals).
- Constraint context: gnomAD v4.

Trait panel, direction of benefit (fixed now):
- Cognitive: fluid intelligence score (higher), reaction time (faster), numeric memory (higher), pairs-matching errors (fewer), education years (higher; flagged as confounded proxy).
- Physical: hand grip strength (higher), FEV1 (higher), walking pace (faster), resting heart rate (lower), systolic BP (lower).
- Metabolic/longevity: LDL (lower), BMI (lower), parental lifespan (longer).
- Panel size N_target = 13. Traits absent from a source are recorded as absent, not substituted.

Trade-off panel (adverse direction = harmful): type 2 diabetes, coronary disease, cancer (any), dementia, major depression, schizophrenia, all-cause mortality, fracture, infertility. N_tradeoff = 9.

Baseline ladder (every rung is run and reported, as counts of genes passing at that rung):
- trivial: phenotype-permuted burden results (carrier status shuffled); expected Tier-A count ~0.
- simplest reasonable: single mask (pLoF only), single trait at a time, discovery p-value only, no replication, no trade-off screen.
- incumbent: published top gene-level hits for the same traits (lookup only, no reanalysis).
- candidate: full pipeline = pLoF mask + damaging-missense mask consistency, replication, trade-off screen, tiering.

Negative controls (must show null):
- Synonymous-variant burden across all genes and traits: genomic inflation lambda_GC must be < 1.10.
- Phenotype-permuted run must yield zero genes below the discovery threshold.

Positive controls (must be recovered, else the pipeline is invalid):
- PCSK9 pLoF: LDL lower at the discovery threshold, and coronary-disease direction protective.
- ANGPTL4 or APOC3 pLoF: triglyceride/lipid-lowering in the beneficial direction.

Thresholds:
- Discovery: gene-level burden p < 0.05 / (20,000 genes x 13 traits) = 1.9e-7, beneficial direction.
- Replication: same direction and one-sided p < 0.05 in the independent cohort, for the same trait or its declared proxy.
- Adverse trade-off: any trade-off outcome with p < 0.05 / 9 in the harmful direction for the gene.
- Tiers: A = discovery + replication + no adverse trade-off. B = discovery, not replicable (trait missing in replication source), no adverse. C = discovery with an adverse trade-off. D = discovery, replication attempted and failed.

Pass: controls valid (positive controls recovered AND lambda_GC < 1.10 AND permuted run null) AND at least one Tier-A gene for a cognitive or physical trait that is not a lipid-pathway gene, with consistent direction under both masks.
Kill: any control fails (pipeline invalid: stop, fix, and start a new dated entry); OR controls valid but zero non-lipid Tier-A genes for cognitive/physical traits (negative result: state plainly that no protective pLoF lever was detectable at this power).
Confirmation: PASS is provisional until (a) both masks agree in direction and (b) the result holds in European-ancestry-only analysis if the source reports it. Tier-A genes are reported as candidates, never as proven levers.
Budget: $0 AWS; local compute and public downloads only; at most 1 working day.
Record: /Users/mike/.dev/genetics-research/results/protective-scan.json

What this can NOT show: causality beyond Mendelian-randomization-style inference, effects of upregulation or drug/editing interventions, effects in non-European ancestries where not reported, effects of variants too rare to be tested at gene level, or anything about polygenic (common-variant) architecture.

### Step 0 feasibility (no outcomes viewed)   (appended 2026-09-29; data access only, no gene or outcome results were opened)

Method: live page fetches and HTTP/bucket probes. VERIFIED = seen on a live page or probe today. UNVERIFIED = search summary only, or not reachable.

**1. Genebass (UKB exomes) -- discovery candidate**
- Bulk data: Hail tables in Google Cloud Storage. Paths quoted by the rivas-lab repo (https://github.com/rivas-lab/phenome-wide-unified-model): `gs://ukbb-exome-public/500k/results/results.mt` (gene burden), `variant_results.mt`, `pheno_results.ht`. Hail 0.2.13 stated as prerequisite. Path-in-code VERIFIED via repo; no official page fetched (genebass.org and the Cell Genomics paper returned SPA shell / 403).
- Bucket is REQUESTER PAYS (VERIFIED by probe: "Bucket is a requester pays bucket but no user project provided", HTTP 400). Any download needs a GCP project with billing enabled. The ledger budget says "$0 AWS", so AWS credit is unaffected, but "public downloads only, $0" is not literally true: a small GCP egress charge is expected. Total size UNVERIFIED (listing blocked without a billing project). Egress at Google list price ~USD 0.12/GB is UNVERIFIED. Mitigation: Hail reads only the needed rows/columns (13 + 9 + control phenotypes, 3 masks), so egress should be a fraction of the full table.
- Cohort: 394,841 exomes in the Cell Genomics paper (https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00110-0, search summary); 426,370 European-ancestry per the Open Targets blog (https://blog.opentargets.org/how-open-targets-integrates-gene-burden/). Which of these the `500k/` folder holds is UNVERIFIED. 4,529 phenotypes (search summary).
- Masks: pLoF, missense|LC, synonymous, pLoF|missense|LC and SKAT-O/burden tests -- UNVERIFIED (recalled; not confirmed on a live page). Must be confirmed at first read of `results.mt` before any outcome is looked at.
- Full gene-by-phenotype table downloadable: YES in principle (matrix table, not only a browser UI). Not tested (needs billing project).
- Local environment: no gcloud, no gsutil, no hail installed; python 3.14.4. Hail compatibility with 3.14 UNVERIFIED; plan a separate Python 3.11 venv.

**2. AstraZeneca PheWAS Portal (UKB exomes) -- consistency check only**
- Portal https://azphewas.com/ is a JS app; landing page has no readable content. ~450k to 484k UKB exomes, 10 to 12 collapsing models, ptv/missense/synonymous mentioned (search summaries; UNVERIFIED on a live page).
- Bulk downloads exist at https://public.cgr.astrazeneca.com/BULK-DOWNLOAD/V01/ (VERIFIED, licence v1.2 effective 26 Nov 2025; research use only, no clinical use, redistribution allowed with the licence, attribution to Wang et al. Nature 2021). The V01 page lists ONLY: CNV gene-level PheWAS (binary + quantitative, per ancestry AFR/ASJ/EAS/NFE/SAS) and variant-level ExWAS files (`exwas_ukb500k_<ANC>_{binary,quantitative}.csv.xz`, ancestries incl. AMR). Gene-level exome collapsing (ptv, missense) is NOT listed there. Direct GET of the files returns 403 without the page's "I agree" step; file sizes UNVERIFIED. Path probes returning HTTP 200 are a catch-all page and prove nothing.
- Consequence: AZ gene-level pLoF results are reachable only through the interactive portal (per-gene/phenotype download UNVERIFIED) -- not a bulk source at present.

**3. Replication candidates (independence from UKB is the constraint)**
- FinnGen (independent of UKB). https://www.finngen.fi/en/access_results: latest DF13 (2 Jun 2026, 500,186 people); gene-based LoF results NOT listed for DF13; earlier releases (DF8, DF10 to DF12) list regenie gene-based burden results (VERIFIED on the results page). Docs: https://finngen.gitbook.io/documentation/methods/lof-variant-burden -- LoF = frameshift, splice donor, stop gained, splice acceptor; MAF <= 0.01; info >= 0.8; 4,793 autosomal genes; 2,751 BINARY endpoints only; regenie max mask. Data path for R11: `/finngen/library-green/finngen_R11/finngen_R11_analysis_data/lof/data` (FinnGen Sandbox naming; whether publicly downloadable UNVERIFIED). Summary-stat download requires an online form; instructions arrive by email. LoF calls come from imputed genotypes (not exome sequencing) and are pLoF-only: no missense or synonymous mask, so the synonymous null cannot be run in FinnGen.
- FinnGen coverage of the panel: none of the 13 target traits is a binary endpoint. Lab-value GWAS (383 measurements, https://labvalues.finngen.fi/) are single-variant, not gene-based (gene-level UNVERIFIED). Trade-off outcomes (T2D, CAD, cancer, dementia, depression, schizophrenia, fracture, infertility) are FinnGen endpoints in principle; exact endpoint names UNVERIFIED. All-cause mortality endpoint UNVERIFIED.
- All of Us "All by All" (independent of UKB): gene-based Hail matrix tables (21 MTs; 3,500+ phenotypes incl. physical and lab measurements) per ancestry and meta-analysis, keyed by gene_id, gene_symbol, annotation, max_MAF (search summary of https://support.researchallofus.org/hc/en-us/articles/27049847988884; the article itself returned 403). Access is inside the Researcher Workbench with permissions -- not an open download; registration/tier requirements UNVERIFIED. Likely covers BMI, SBP, pulse, LDL; unlikely to cover the cognitive tests, grip, reaction time.
- Regeneron (RGC-ME): browser https://rgc-research.regeneron.com/me has variant-level data (983,578 individuals); the burden results integrated in Open Targets are from the same UKB data (454,787), so NOT independent of UKB. Gene-burden bulk download UNVERIFIED.
- Biobank Japan / Taiwan Biobank / Mount Sinai BioMe: no open gene-based burden tables found (array-based GWAS only for TWB/BBJ). UNVERIFIED that none exist.
- UKB overlap: Genebass, AZ and Regeneron-in-Open-Targets all use the same UKB exomes (Open Targets blog: "All three sources utilize the same UK Biobank exome sequencing data"). They are one cohort, not replication of each other.

**4. gnomAD v4.1 constraint -- available, small**
- `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv` -- HTTP 200, 95,546,041 bytes (~91 MiB), last modified 18 Apr 2024 (VERIFIED). Hail table `gnomad.v4.1.constraint_metrics.ht` and `README.txt` (column dictionary incl. `lof.oe_ci.upper` = LOEUF, `lof.pLI`) also present (VERIFIED). No registration. Licence terms not fetched (UNVERIFIED; gnomAD is generally open).

**5. Trait coverage at the UKB level (13 targets)**
Checked against the open Pan-UKB phenotype manifest (https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz, 7,223 rows, VERIFIED). This shows the trait exists in UKB; whether each is included in Genebass or AZ is UNVERIFIED (Genebass's own phenotype list not retrieved).
| Target trait | UKB field / manifest entry | Max n (EUR cases col) | Note |
|---|---|---|---|
| Fluid intelligence | 20016 | 135,088 | subset with the test; low power for rare genes |
| Reaction time | 20023 "Mean time to correctly identify matches" | 417,660 | |
| Numeric memory | 4282 "Maximum digits remembered correctly" | 43,741 | very low power |
| Pairs matching errors | 399 "Number of incorrect matches in round" | 419,951 | |
| Education years | 845 "Age completed full time education" (22501 "Year ended full time education", n 106,229) | 283,575 | not literally years of education; already flagged as confounded proxy |
| Hand grip | 46 (left), 47 (right) | 418,776 / 418,827 | two fields; combine rule needed |
| FEV1 | 3063 (also 20150 best measure) | 383,471 | |
| Walking pace | 924 | 417,933 | ordinal |
| Resting heart rate | 102 | 396,667 | |
| Systolic BP | 4080 (also medication-adjusted `SBP`) | 396,663 | |
| LDL | 30780 (also medication-adjusted `LDLC`) | 400,223 | |
| BMI | 21001 | 419,163 | |
| Parental lifespan | 1807 (father's age at death), 3526 (mother's age at death) | 310,232 / 249,247 | not a single field; combine rule needed |
All 13 exist in UKB. Absent-from-source handling stays as ledger says.

**6. Recommended assignment (keeps cohorts independent)**
- Discovery: Genebass (UKB). Pull only: 13 panel traits + 9 trade-off outcomes + PCSK9/ANGPTL4/APOC3 lipid control phenotypes; pLoF, missense|LC and synonymous masks.
- Within-UKB consistency (NOT replication, not counted in tiers): AZ portal pLoF for the same genes, by hand through the portal or a later bulk file.
- Replication: FinnGen gene-based LoF (binary endpoints; trade-off screen, plus panel-trait proxies) and All of Us All by All (BMI, SBP, pulse, LDL, physical measurements) if access is obtainable. Cognitive traits (fluid intelligence, reaction time, numeric memory, pairs matching, education) have no independent gene-level replication source found; by the ledger's tier rules these can reach Tier B only.
- Constraint: gnomAD v4.1 TSV.
- Estimated download: gnomAD 91 MiB (verified). Genebass filtered pull: UNVERIFIED, expected in the low GB range with Hail row/column filtering (full `results.mt` likely much larger). FinnGen gene-based tables and All of Us: UNVERIFIED.

**7. Issues this raises for the ledger rules (recommendation only; the entry above is NOT edited)**
1. Budget: requester-pays means the public-download step needs a GCP project with billing (small GCP charge; AWS credit unaffected).
2. PASS reachability: PASS requires a Tier-A non-lipid gene for a cognitive or physical trait. Tier A requires replication in the same trait "or its declared proxy". FinnGen has binary endpoints only and All of Us lacks cognitive traits, so unless proxies are declared, cognitive traits cap at Tier B and PASS is nearly unreachable for them. Proxies for physical traits (e.g. FEV1 -> COPD/asthma endpoints; SBP -> hypertension; BMI -> obesity; resting HR -> arrhythmia endpoints) must be declared in a new dated entry BEFORE any outcome is viewed.
3. Permutation controls cannot be run from summary statistics: the trivial rung ("phenotype-permuted burden results") and the negative control ("permuted run yields zero genes") require individual-level data, which the ledger rules out ("public gene-based burden results only"). Replace with summary-level nulls in a new dated entry (synonymous-mask null and lambda_GC already exist; a summary-level substitute for the permuted rung is a decision for the ledger owner).
4. Synonymous null cannot be run in FinnGen (pLoF-only) or in AZ bulk (not listed); only Genebass (mask existence UNVERIFIED) and possibly the AZ portal.
5. Power: numeric memory (n ~44K) and fluid intelligence (n ~135K) make gene-level rare pLoF tests weak at p < 1.9e-7; expect few or no cognitive discoveries regardless of biology. That is a power limit, not evidence against protective variants.
6. Genebass release ambiguity (394,841 vs 426,370; folder `500k`) must be settled from the table metadata at first read, before phenotype rows are opened.

Human actions needed before the run: (a) a GCP project with billing enabled and `gcloud auth application-default login` for requester-pays reads; (b) submit the FinnGen download form (https://elomake.helsinki.fi/lomakkeet/124935/lomake.html) and note which release has gene-based files; (c) decide whether to pursue All of Us Researcher Workbench access; (d) approve a dated amendment entry for items 1 to 3 above.

## 2026-09-29 protective-variant-scan AMENDMENT 1   (rules fixed before any outcome was viewed)

Why: Step 0 feasibility (above) found three facts about data access. No gene-level or trait-level result had been read when this amendment was written; the earlier entry's thresholds, trait panel, direction of benefit, trade-off panel and tier letters are unchanged.
1. No independent replication exists for continuous cognitive/physical traits: FinnGen gene-based LoF covers binary endpoints only; Genebass, AstraZeneca and Regeneron all use the same UK Biobank exomes; All of Us access is unverified.
2. A phenotype-permuted run cannot be produced from summary statistics.
3. The Genebass bucket is requester-pays (needs a GCP billing project and login).

Changes (these replace the corresponding lines of the original entry; all else stands):
- Replication for Tier A requires a cohort independent of UK Biobank. Declared proxies (fixed now): LDL -> FinnGen hypercholesterolemia endpoint; BMI -> FinnGen obesity endpoint; SBP -> FinnGen hypertension endpoint. Exact endpoint codes are recorded in docs/data-sources.md BEFORE any gene-level result table is opened. Cognitive traits, grip strength, FEV1, walking pace, resting heart rate, parental lifespan have no declared proxy: Tier B is their ceiling.
- Within-UKB consistency (AstraZeneca portal, or mask agreement) is reported but never counted as replication.
- Negative control: the phenotype-permuted run is replaced by the synonymous mask: lambda_GC < 1.10 AND zero synonymous-mask genes at the discovery threshold in the beneficial direction. The trivial rung is the synonymous mask.
- Verdicts: PASS = controls valid AND >= 1 non-lipid Tier-A gene (trait with a declared proxy) with both masks consistent. LEAD = controls valid, no non-lipid Tier-A gene, but >= 1 non-lipid Tier-B gene for a cognitive or physical trait (discovery p < 1.9e-7, both masks consistent, no adverse trade-off); reported as an UNREPLICATED lead only. KILL = any control fails, or controls valid and neither PASS nor LEAD. Decision: PASS -> follow-up as originally stated. LEAD -> follow-up flagged unreplicated, plus seek replication access. KILL -> negative result, thread 1 becomes primary.
- Access: using a GCP billing project or any login/payment is NOT authorized. Preferred route is any free, non-requester-pays access to Genebass gene-level results. If none exists, real-data stages are parked pending the user; synthetic and gnomAD stages proceed.
Budget: unchanged ($0 AWS, no GCP spend without the user).

### 2026-09-29 Amendment 1 clarifications (from builder-pipeline's ambiguity list; written before any gene-level result was opened by the lead; data may have been downloaded by builder-data but was not to be analyzed)

- C1: A gene is Tier A only if at least 5 of the 9 trade-off outcomes were screened for it (plof mask). Fewer than 5 screened caps it at Tier B, labeled "trade-off unscreened". Reason: "no adverse trade-off" with no data is not evidence.
- C2: PASS and LEAD require a trait in the cognitive or physical domain. Blood pressure (SBP) counts as physical; LDL, BMI and parental lifespan do not qualify. Consequence stated plainly: PASS is reachable only through SBP (the one physical trait with a declared proxy); every cognitive result is capped at LEAD.
- C3: The lipid-pathway gene list is fixed as `verdict.lipid_pathway_genes` in config/prereg.yaml at commit 4073b89. Additions after this need a new dated entry.
- C4: Accepted as implemented: discovery threshold 1.9e-7 (not 1.923e-7); ANGPTL4/APOC3 positive control on triglycerides at the discovery threshold, PCSK9 coronary direction only; genes untested in the replication cohort are Tier B (not D); trade-off screen is one-sided, plof only; a control that cannot run gives KILL "controls_not_evaluable"; EUR-only is reported and does not gate the verdict; replication sources are allow-listed (finngen, all_of_us) and any UK Biobank token disqualifies.

### 2026-09-29 Amendment 1 data clarifications (from builder-data; written before any pipeline run on real tables)

Data in hand: Genebass discovery via the open REST API (CC BY 4.0), FinnGen R13 LoF burden, gnomAD v4.1 constraint. Details, URLs, checksums: docs/data-sources.md.
- C5: `dmis` mask = Genebass `missense|LC` (missense plus low-confidence pLoF); Genebass has no damaging-missense filter. The "both masks consistent" rule is therefore direction-only agreement between pLoF and missense|LC. This is weaker than intended and is a known limitation, reported with any LEAD or PASS.
- C6: Burden test only (BETA_Burden, Pvalue_Burden). SE is reconstructed as |beta| / z(p); carriers approximated from allele frequency. Only signs are compared across Genebass and FinnGen (different scales). Gene QC: drop genes failing Genebass coverage or variant-count flags; the lambda-based flags are NOT applied (they derive from synonymous inflation and would make the negative control circular); ambiguous gene symbols are dropped. Trait-to-analysis mappings (2.1 and 3.2 in docs/data-sources.md) were fixed by rule before results were read; changing any needs a dated entry.
- C7: Traits absent from discovery, recorded and not substituted: numeric memory, pairs-matching errors, education years, walking pace, all-cause mortality. The discovery threshold stays at 0.05/(20000 x 13) = 1.9e-7 (conservative; 9 of 13 panel traits are present). Cognitive coverage is fluid intelligence and reaction time only.
- C8: FinnGen LoF calls are VEP-based without LOFTEE, so not the same mask as Genebass pLoF; only directional agreement is tested for Tier A. FinnGen documentation asks users to submit its request form; that has NOT been submitted. Hold: FinnGen replication tables are not opened or analyzed until the user has submitted it.
- Disclosure: while debugging the file format, builder-data printed two FinnGen rows (LDLR and DHPS, hypercholesterolemia). No rule was set from them; LDLR is not one of the positive-control genes (PCSK9, ANGPTL4, APOC3).

- C9 (from builder-pipeline's judgment call, decided by the lead before any real-data run): the C1 rule applies to LEAD as well. A Tier B gene counts toward LEAD only if at least 5 of the 9 trade-off outcomes were screened for it; otherwise it is listed as "Tier B, trade-off unscreened" and does not count. Reason: LEAD asserts "no adverse trade-off", which is vacuous with no screening.

### 2026-09-29 Amendment 1 review clarifications C10 (from reviews/review-1.md, decided by the lead before any real-data run)

- C10a (R-1): The synonymous negative control is valid only if the synonymous rows cover at least 90 percent of the (gene, trait) pairs that have a pLoF discovery row AND number at least 10,000. Otherwise status is controls_not_evaluable and the verdict is KILL.
- C10b (R-2): Replication-sign control, from established biology, not from any table: in the replication cohort, PCSK9 pLoF for the hypercholesterolemia proxy must have a negative effect (lower risk), and LDLR pLoF a positive effect (higher risk). At least one of the two must be evaluable; every evaluable one must hold; otherwise KILL. The adapter must assert A1FREQ <= 0.5 for every kept FinnGen row (the LoF carrier allele is A1), and drop-and-log any row failing it. Disclosure: an LDLR row was among the two FinnGen rows printed earlier; this control's direction comes from biology (LDLR loss of function causes familial hypercholesterolemia), not from that row.
- C10c (R-3): The run compares the config's sha256 with a pinned value; on mismatch the report and JSON stamp config_matches_ledger=false and the verdict is labeled NON-PREREGISTERED, and the CLI exits nonzero. The pinned hash changes only with a dated ledger entry.
- C10d (R-9): Rows with p = 1 or beta = 0 are retained for lambda_GC and hit counting (which use p only); se stays undefined for them.
- C10e (R-7): Informational only, no gating: the report adds, per Tier A/B gene, any significant (discovery threshold) harmful-direction association on the other panel traits.
- C10f (R-8): Any PASS or LEAD report prints the C5 (missense|LC) and C7 (cognitive coverage limited to fluid intelligence and reaction time) caveats.
- C10g (R-6, R-4, R-10): screening counts use only allow-listed replication rows; requests and pyarrow are declared in pyproject; positive controls are evaluated through the pipeline's own tier path.

### 2026-09-29 C8 hold lifted
The user registered for FinnGen summary-statistics access (confirmation email from the FinnGen service desk, received 5:54 PM on 2026-09-29 and pasted into the session). The C8 hold on opening the FinnGen replication tables is lifted, effective after the C10 fixes and reviewer re-check. FinnGen asks that publications acknowledge "the participants and investigators of the FinnGen study" and cite Kurki et al., Nature 613:508-518 (2023), doi:10.1038/s41586-022-05473-8; recorded in README.md.

### 2026-09-29 Amendment 1 review-2 clarifications C11 (from reviews/review-2.md and signoff.md; decided by the lead before any real-data run)

- C11a (R2-2): A replication row counts toward the replication rule only if its se is finite and > 0 and its beta is nonzero. Rows with undefined se stay available for lambda_GC and coverage counts only. The FinnGen adapter additionally drops and logs any row with undefined se and p < 1 (a data defect), with counts in its conversion summary.
- C11b (R2-4): lambda_GC for the synonymous control is computed on synonymous rows with p < 1, and both lambdas (with and without p = 1 rows) are reported. If p = 1 rows exceed 5 percent of the synonymous rows, the control is not evaluable (KILL, controls_not_evaluable).
- C11c (R2-3): In the replication-sign control (C10b), an arm (PCSK9 negative, LDLR positive) is evaluable only if its p < 0.05. At least one arm must be evaluable and every evaluable arm must match its expected sign; otherwise KILL, with the reason recorded as "sign_control_not_evaluable" (no evaluable arm) or "sign_control_failed" (an evaluable arm has the wrong sign). Underpowered arms never count as failures.
- C11d (R2-5): The positive controls (PCSK9 LDL-lowering; ANGPTL4 or APOC3 triglyceride-lowering) are judged on the discovery effect only (direction and discovery threshold). The trade-off screen, including PCSK9 and type 2 diabetes, never changes control status.
- Reading rule for any KILL involving the sign control or positive controls: report it as a pipeline/data conclusion only after checking C11c and C11d; a biology conclusion requires all controls valid.

## 2026-09-29 protective-variant-scan RUN 1 RESULT   (appended; the rules above were not changed for this run)

Run: `python -m protscan run --config config/prereg.yaml --data data --out results/protective-scan.json` at commit 2b373c3, config sha256 7ce18c38... (config_matches_ledger = true, not synthetic). Record: results/run1-protective-scan.json and results/run1-report.md (renamed after the run so run 2 cannot overwrite them). Tier and gene lists in those files were NOT opened by the lead before the amendment below was written.

Controls (calibration only):
| Control | Result |
|---|---|
| Positive: PCSK9 LDL-lowering | OK, beta -0.0389, p 3.6e-132 |
| Positive: PCSK9 coronary protective | OK (direction), beta -0.0296, p 1.4e-3 |
| Positive: APOC3 triglyceride-lowering | OK, beta -0.0493, p 1e-300 (floor) |
| Replication sign: LDLR (expected +) | OK, beta +2.885, p 1.7e-19 |
| Replication sign: PCSK9 (expected -) | NOT RUN, no informative row; allowed by C11c since one arm evaluated |
| Synonymous: lambda_GC (p<1 rows) | 1.070 (all rows 1.056), limit 1.10: OK |
| Synonymous: coverage / rows / p=1 fraction | 99.9% / 334,441 / 0.57%: OK |
| Synonymous: beneficial-direction hits at discovery threshold | 3 genes (CD3EAP, HPR, ZNF224); limit 0: FAIL |

Verdict by the pre-written rule: KILL (control_failed: negative_synonymous). The pipeline is invalid as specified; this verdict stands. It is not a biological conclusion (C11 reading rule: all controls must be valid for one).
What this shows: discovery and FinnGen effect-sign conventions reproduce known biology. Global inflation is small. Three genes show a beneficial-direction signal in a mask that cannot reflect loss of function. About 334,000 synonymous rows tested at 1.9e-7 would yield about 0.06 hits by chance, so 3 is not chance.
Interpretation, UNVERIFIED (from memory; the local constraint table has no coordinates): CD3EAP and ZNF224 lie in the APOE/TOMM40 region (chr19q13.32) and HPR beside HP (chr16q22), regions with strong common-variant associations. A synonymous burden there most plausibly picks up common variants in linkage disequilibrium, which would contaminate pLoF burden results in the same regions.

## 2026-09-29 protective-variant-scan AMENDMENT 2   (written after run 1 failed a control, before any tier or gene list was opened)

Honest framing: this changes the negative control after seeing it fail. Run 1 stands as KILL under the original rule. Run 2 is a new run under the rules below; if it passes, it is reported as "passed under Amendment 2, which was written after run 1 failed", never as a clean preregistered pass.
- A2a (contamination filter): any gene with a synonymous-mask hit at the discovery threshold (1.9e-7, either direction, any trait) is CONTAMINATED. Contaminated genes are excluded from all tiers and candidate lists and are listed in the report with the trait(s).
- A2b (control gate): the negative control is valid iff lambda_GC < 1.10, coverage >= 90 percent, >= 10,000 rows, p = 1 fraction <= 5 percent (all unchanged), AND contaminated genes number at most 0.1 percent of the genes tested (about 18 of roughly 18,500). Beyond that the contamination is systemic: KILL. The 0.1 percent bound was chosen knowing run 1 had 3; it is a post-hoc bound and is disclosed as a limitation.
- A2c (locus caveat, manual): every PASS or LEAD gene is reported with its chromosome position and a GWAS Catalog lookup of the same trait within 500 kb, recorded in the ledger. Such a neighbouring signal downgrades the gene to "locus-contaminated lead" and it does not count toward PASS or LEAD.
- A2d: everything else (thresholds, panel, tiers, C1-C11, proxies, verdict definitions) is unchanged. The config pin changes only with the implementing commit, and run 2 is not started until the reviewer re-verifies.

## 2026-09-29 protective-variant-scan RUN 2 RESULT   (under Amendment 2, written after run 1 failed a control)

Run: `python -m protscan run --config config/prereg.yaml --data data --out results/run2-protective-scan.json` at commit be7f43f, config sha256 be698275... (config_matches_ledger = true, not synthetic). Reviewer sign-off for this run: reviews/signoff.md (round 4). Record: results/run2-protective-scan.json, results/run2-report.md. Label on every output: Amendment 2 (post-hoc after run 1).

Controls: all valid. Positive controls (PCSK9 LDL, PCSK9 coronary direction, APOC3 triglycerides) OK; replication sign OK (LDLR arm; PCSK9 arm not run, no informative row, allowed by C11c); synonymous lambda_GC 1.070; contaminated genes 11 of about 18,700 tested (0.059 percent, limit 0.1 percent): CD300LG, CD3EAP, DNM2, HPR, MLXIPL, PHF7, SIDT2, SLC22A3, SLC30A3, TLR9, ZNF224. These 11 (either direction, any trait) are excluded from all tiers. Run 1 counted 3 because it only counted the beneficial direction.

Ladder (gene counts):
| Rung | Result |
|---|---|
| trivial (synonymous mask, discovery rule, beneficial direction) | 3 genes |
| simplest (pLoF only, one trait at a time, discovery p only) | 8 gene-trait pairs: LDL 7, BMI 1; none for any cognitive or physical trait |
| incumbent (published top hits) | NOT RUN (no file provided; not fabricated) |
| candidate (masks consistent, replication, trade-off screen, tiers) | Tier A 1, Tier B 5, Tier C 2, Tier D 0 |

Candidates (all metabolic domain; none qualifies for PASS or LEAD under C2):
- Tier A: APOC3 (LDL; replicated in FinnGen; lipid-pathway gene).
- Tier B (not replicable or not replicated): APOB, PCSK9, ANGPTL3 (LDL, lipid-pathway genes), RRBP1 (LDL), GPR151 (BMI).
- Tier C (beneficial but adverse trade-off): GIGYF1 (LDL lowering; adverse type 2 diabetes and major depression), GCK (LDL lowering; adverse type 2 diabetes).

Verdict by the pre-written rule: KILL, controls valid, neither PASS nor LEAD. This is the negative-result branch. Reported as: no protective pLoF lever for a cognitive or physical trait was detectable in these public tables at this power, under Amendment 2 (post-hoc after run 1).
What this shows: the pipeline recovers known protective lipid genes (APOB, PCSK9, ANGPTL3, APOC3) and a known BMI gene (GPR151), and the trade-off screen separates protective genes from ones with a cost (GIGYF1, GCK). No brain or body-function gene reached the discovery threshold 1.9e-7.
What this does not show: that no such variants exist. Power is limited (fluid intelligence and reaction time only for cognition; four panel traits and mortality absent; gene-based pLoF only, discovery threshold stringent for 13 traits); cognitive traits have no independent replication, so LEAD was the ceiling and was not reached. A2c (manual locus lookup) was not needed because there is no PASS or LEAD gene.
Decision it drives (pre-written): KILL -> negative result stands; thread 1 (cognition cell-type mapping) becomes the primary project.
