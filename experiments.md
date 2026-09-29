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
