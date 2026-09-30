# Protective-variant scan: report

Ledger: 2026-09-29 protective-variant-scan (with AMENDMENT 1)  
Config sha256: 7ce18c381814ba0b  
Data: data (burden_finngen.csv.gz, burden_genebass.csv.gz; 1025318 rows; 0 NaN rows dropped)  
Sources: discovery: genebass; replication: finngen_r13

## Verdict: KILL

control_failed: negative_synonymous

- PASS-qualifying genes (Tier A): none
- LEAD-qualifying genes (Tier B): none
- Qualifying = non-lipid-pathway gene, cognitive or physical trait, both masks consistent. Tier-A genes are candidates, never proven levers.

## Controls

| control | status | detail |
|---|---|---|
| synonymous coverage >= 0.9 of pLoF pairs, rows >= 10000 | OK | coverage = 0.999 of 308172 pairs; 334441 rows |
| synonymous lambda_GC < 1.1 | OK | lambda_GC = 1.070 on rows with p < 1; 1.056 on all 334441 rows; p = 1 rows: 1899 (0.006, max 0.05) |
| synonymous genes at discovery threshold (beneficial) <= 0 | FAIL | 3 genes CD3EAP, HPR, ZNF224 |
| pcsk9_ldl_lower (PCSK9, ldl, at discovery threshold, also via the pipeline discovery-hit path (trade-offs never change status)) | OK | PCSK9 beta=-0.0389 p=3.6e-132 tier=B |
| pcsk9_cad_protective (PCSK9, coronary_disease, beneficial direction) | OK | PCSK9 beta=-0.0296 p=0.0014 |
| lipid_lowering (ANGPTL4/APOC3, triglycerides, at discovery threshold) | OK | APOC3 beta=-0.0493 p=1e-300 |
| pcsk9_rep_sign (PCSK9 pLoF, hypercholesterolemia, expect negative) | NOT RUN | no informative row in independent replication cohort; not evaluable: no informative row |
| ldlr_rep_sign (LDLR pLoF, hypercholesterolemia, expect positive) | OK | beta=2.89 p=1.7e-19 |

Replication-sign control (an arm is evaluable only if p < 0.05; at least one arm evaluable, every evaluable arm must match): **OK**

Controls valid: **False**

Reading rule: a KILL involving the sign control or the positive controls is a pipeline/data conclusion only after checking that the sign-control arms were adequately powered (C11c) and that the positive controls were judged on the discovery effect alone (C11d). A biology conclusion requires all controls valid.

## Baseline ladder (genes passing at each rung)

| rung | status | genes | note |
|---|---|---|---|
| trivial (synonymous mask) | RUN | 3 | expected 0 |
| simplest (pLoF, per trait, discovery p only) | RUN | 8 | 8 gene-trait pairs |
| incumbent (published top hits) | NOT RUN | NA | incumbent_hits.csv not provided; published top hits are not fabricated |
| candidate (full pipeline, Tier A/B with dmis consistent) | RUN | 4 | genes by tier A/B/C/D = 1/5/2/0 |

## Tiers

### Tier A: discovery + independent replication + no adverse trade-off (genes: 1, gene-trait pairs: 1)

| gene | trait | beta | p | replication | dmis | EUR-only | adverse | screened | harmful panel (info) | lipid | qualifies | LOEUF | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| APOC3 | ldl | -0.00718 | 1.14e-14 | replicated (hypercholesterolemia, p1=0.0366) | consistent | not_reported | - | 8/9 | - | yes | no | 1.68 |  |

### Tier B: discovery, not replicable, no adverse trade-off (genes: 5, gene-trait pairs: 5)

| gene | trait | beta | p | replication | dmis | EUR-only | adverse | screened | harmful panel (info) | lipid | qualifies | LOEUF | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| APOB | ldl | -0.0465 | 1.8e-299 | gene_untested | inconsistent | not_reported | - | 8/9 | - | yes | no | 0.537 |  |
| PCSK9 | ldl | -0.0389 | 3.58e-132 | gene_untested | consistent | not_reported | - | 8/9 | - | yes | no | 1.25 |  |
| ANGPTL3 | ldl | -0.0148 | 1.29e-26 | gene_untested | consistent | not_reported | - | 8/9 | - | yes | no | 1.21 |  |
| GPR151 | bmi | -0.0027 | 5.58e-08 | gene_untested | inconsistent | not_reported | - | 8/9 | - | no | no | 1.15 |  |
| RRBP1 | ldl | -0.0184 | 1.14e-07 | gene_untested | consistent | not_reported | - | 8/9 | - | no | no | 0.704 |  |

### Tier C: discovery with an adverse trade-off (genes: 2, gene-trait pairs: 2)

| gene | trait | beta | p | replication | dmis | EUR-only | adverse | screened | harmful panel (info) | lipid | qualifies | LOEUF | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GIGYF1 | ldl | -0.0226 | 6.45e-13 | failed (hypercholesterolemia, p1=0.126) | consistent | not_reported | type_2_diabetes, major_depression | 8/9 | hand_grip_strength (p=6.91e-10) | no | no | 0.735 |  |
| GCK | ldl | -0.0355 | 4.61e-08 | gene_untested | inconsistent | not_reported | type_2_diabetes | 8/9 | - | no | no | 0.866 |  |

### Tier D: discovery, replication attempted and failed (genes: 0, gene-trait pairs: 0)

(none)

## Notes

- Replication counts only a cohort independent of UK Biobank; only traits with a declared proxy (LDL, BMI, SBP) can reach Tier A. Within-UKB consistency (dmis mask) is reported and never counted as replication; AstraZeneca portal lookups are not automated.
- EUR-only column is reported and does not gate the verdict. 'screened' = trade-off outcomes with a pLoF row for the gene; unscreened outcomes cannot be excluded as adverse.
- Tier A needs >= 5 of 9 trade-off outcomes screened (pLoF) for the gene (discovery plus allow-listed independent replication rows); otherwise a replicated gene is capped at Tier B and labeled 'trade-off unscreened'. The same minimum applies to LEAD (C9): unscreened Tier-B genes are listed but do not count.
- Traits absent from discovery: all_cause_mortality, education_years, numeric_memory, pairs_matching_errors, walking_pace.
- Traits absent from replication: all_cause_mortality, fracture.
