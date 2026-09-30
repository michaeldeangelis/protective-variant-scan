# Normalized table interface (fixed)

Adapters emit, and `python -m protscan run` reads, a long table: one row per gene x trait x mask x cohort.

| column | type | meaning |
|---|---|---|
| gene | str | HGNC symbol (upper-cased on load) |
| trait | str | canonical name from `config/prereg.yaml` (lower snake_case; see below) |
| mask | str | `plof`, `dmis` (damaging missense), `syn` (synonymous) |
| cohort | str | role of the row: `discovery`, `replication`, `discovery_eur` |
| beta | float | effect of carrying the mask on the raw trait scale (log-odds/log-hazard for binary outcomes). Sign is NOT flipped to "beneficial"; direction of benefit lives in config |
| se | float | standard error of beta |
| p | float | two-sided p-value as reported by the source |
| n_carriers | int | carriers of the mask (cases + controls) |
| n_total | int | total individuals tested |
| source | str | dataset name (e.g. `genebass`, `finngen_r12`); a name starting with `synthetic` triggers a banner in the report |

Key `(gene, trait, mask, cohort)` must be unique across all files. Rows with NaN beta or p are dropped and counted; se may be NaN (rows with beta = 0 or p = 1 are kept for lambda_GC and hit counting). Invalid mask/cohort, p outside [0,1], negative se, or duplicate keys raise an error.

## Cohorts

- `discovery`: UK Biobank exome burden (Genebass). Masks `plof`, `dmis`, `syn`. Traits: 13 panel traits, 9 trade-off outcomes, `triglycerides` (control).
- `replication`: cohort independent of UK Biobank (FinnGen etc.). Only `plof` is read. Traits: the declared proxy names (below) and/or trade-off outcomes. Same-trait rows for `ldl`, `bmi`, `systolic_bp` are also accepted. Rows for other panel traits are ignored for tiering (Amendment 1: Tier B ceiling).
  Independence is enforced on the `source` string (case-insensitive substring): it must contain an allow token (`finngen`, `all_of_us`, `allofus`, `synthetic_replication`) and no UK Biobank token (`genebass`, `azphewas`, `astrazeneca`, `regeneron`, `rgc`, `ukb`, `uk_biobank`, `opentargets`, `backman`, `pan_ukb`, ...). Rows from any other source are ignored for replication and listed under `data.replication_sources_excluded` in the JSON. Tokens live in `config/prereg.yaml` (`replication.independent_sources`, `replication.ukb_overlapping_sources`); adding one needs a dated ledger entry.
- `discovery_eur`: optional EUR-only re-analysis of discovery, `plof` mask, same traits. Reported, not gating.

## Files

Read from the top level of the `--data` directory only:

- `burden*.csv` or `burden*.csv.gz`: normalized tables (any number; concatenated). Suggested: `burden_genebass.csv.gz`, `burden_finngen.csv.gz`.
- `incumbent_hits.csv` (optional): columns `gene,trait[,source]`; published top gene-level hits for the panel traits. Absent -> incumbent rung reported NOT RUN.
- `constraint.csv` (optional): columns `gene,loeuf[,pli]`; annotates Tier A/B genes in the report.

## Run gates the data must satisfy (Amendment 1 C10)

- Synonymous control: `syn` discovery rows must number >= 10,000 and cover >= 90% of the (gene, trait) pairs that have a `plof` discovery row, and rows with p = 1 must be at most 5% of them, else controls_not_evaluable (KILL). lambda_GC is computed on rows with p < 1 (both lambdas are reported).
- Contamination (Amendment 2, post-hoc after run 1): a gene with a `syn` discovery hit at the discovery threshold (either direction, any trait) is excluded from all tiers and candidate lists and listed in the report and JSON; the negative control fails if such genes exceed 0.1% of the genes tested. The old zero-hit gate is replaced (its count stays in the report).
- Replication rule (C11a): a replication row counts only if se is finite and > 0 and beta is nonzero. Rows with undefined se (beta = 0 or p = 1) are kept in the table for lambda_GC and coverage counts only.
- Replication-sign control: the independent replication cohort should carry `plof` rows for `hypercholesterolemia` for PCSK9 (must be negative) and LDLR (must be positive). An arm is evaluable only with an informative row and p < 0.05; at least one arm must be evaluable and every evaluable arm must match (`sign_control_not_evaluable` / `sign_control_failed`). The adapter must guarantee the LoF carrier-allele sign convention (A1FREQ <= 0.5 guard).
- A data directory must not mix `synthetic*` sources with real ones; the run refuses.
- The config sha256 is pinned in `src/protscan/schema.py`; any other config is stamped NON-PREREGISTERED and the CLI exits 2.

## Trait names

Panel (13): `fluid_intelligence, reaction_time, numeric_memory, pairs_matching_errors, education_years, hand_grip_strength, fev1, walking_pace, resting_heart_rate, systolic_bp, ldl, bmi, parental_lifespan`.
Trade-off (9): `type_2_diabetes, coronary_disease, cancer_any, dementia, major_depression, schizophrenia, all_cause_mortality, fracture, infertility`.
Control: `triglycerides`.
Declared replication proxies (Amendment 1; map the FinnGen endpoint to these names): `ldl -> hypercholesterolemia`, `bmi -> obesity`, `systolic_bp -> hypertension`.

Scale conventions: `reaction_time` in time units (lower = faster); `pairs_matching_errors` counts; `walking_pace` coded so higher = faster; binary outcomes as log-odds. Traits absent from a source are simply not emitted (recorded as absent, never substituted).

## Synthetic data

`python scripts/make_synthetic.py --out data/synthetic --scenario pass` writes `burden_synthetic_{discovery,replication}.csv.gz`. Scenarios: `pass, lead, nolead, broken_positive, broken_lambda, broken_syn_hit`.
