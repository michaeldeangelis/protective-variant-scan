# Data sources (owner: builder-data)

Status log. Written 2026-09-29. Sections 2 and 3 were fixed BEFORE any gene-level result table was opened
(only phenotype metadata, the FinnGen manifest/readme and object listings had been read). Section 6 is appended after downloads.

Analysis only. Public, anonymous HTTPS only: no login, no billing project, no application-gated data, no individual-level data.

## 1. Sources at a glance

| Role | Source | Route | Access verified 2026-09-29 | Size |
|---|---|---|---|---|
| discovery | Genebass (UKB exomes, 394,841 per app home page, 4,529 phenotypes) | public REST API behind app.genebass.org: `https://main.genebass.org/api` | HTTP 200, no auth, no project; data CC BY 4.0 (app /terms page) | ~20 analyses x 3 masks, filtered on the fly (MBs) |
| discovery, NOT USED | Genebass Hail tables `gs://ukbb-exome-public/500k/results/*.mt` | GCS | HTTP 400 "requester pays bucket but no user project provided" (also for the `pheno_results.txt.bgz` link on the app downloads page) | not obtainable without a billing project |
| discovery, consistency only | AstraZeneca PheWAS portal | none in bulk | V01 bulk page lists only CNV gene-level and variant-level ExWAS (Step 0); gene-level exome collapsing not listed | not obtained |
| replication | FinnGen R13 LoF burden (regenie, Finnish biobanks) | public GCS bucket `finngen-public-data-r13`, anonymous HTTPS | listing and HEAD succeed without project or login | `lof/data/finngen_R13_lof.txt.gz` 548,499,199 B (streamed and filtered) |
| constraint | gnomAD v4.1 constraint metrics | public GCS bucket `gcp-public-data--gnomad` | HTTP 200 | 95,546,041 B |

## 2. Genebass: free route found, and what was tried

Tried, in order:
1. `gs://ukbb-exome-public/500k/results/{results.mt,variant_results.mt,pheno_results.ht}` (paths quoted on the app /downloads page): requester-pays. `GET storage/v1/b/ukbb-exome-public/o` and `GET https://storage.googleapis.com/ukbb-exome-public/500k/results/pheno_results.txt.bgz` both return HTTP 400 UserProjectMissing. Not usable without a GCP billing project (NOT authorized, ledger Amendment 1).
2. genebass.org redirects to app.genebass.org (React app). The bundle `main-ed3ab9ea1857efb02221.js` (ukbb-exome-ui 0.13.0, build 2024-02-23) hard-codes an open API base `https://main.genebass.org/api`. Routes seen in the bundle: `/phenotypes`, `/phenotypes/qc`, `/categories`, `/phenotype/{analysis_id}`, `/analysis/{analysis_id}/gene-manhattan?burdenSet={set}`, `/gene-qc/burden-set/{set}`, `/gene/{ensembl_id}`, `/top-associations/{set}`, others.
3. `GET /api/phenotypes` (metadata only): HTTP 200, no auth, CORS `*`, 4,529 analyses, gzip JSON, ~267 KB.
4. Per-analysis, all-gene results: `/api/analysis/{analysis_id}/gene-manhattan?burdenSet={set}`. Used with one request at a time and a delay (polite use of the app backend; no bulk crawl beyond the ~20 analyses needed).

Terms (app /terms page, read): CC BY 4.0, no additional restrictions or embargo, do not attempt re-identification. Attribution: Karczewski et al., Cell Genomics 2022, doi:10.1016/j.xgen.2022.100168, and UK Biobank.

Caveat: the API is an app backend, not a documented bulk service; it can change or rate-limit. Mask names (`burdenSet` values) were taken from the app code and are confirmed at first fetch (section 6).

### 2.1 Genebass analysis mapping (fixed before any result was read)

Rule: a panel trait maps to the Genebass analysis whose description matches it; where several exist, prefer the one with a documented definition, then the larger case/sample count. `pheno_sex=both_sexes` only. Alternatives are listed so that nothing is substituted silently.

| Normalized trait | Genebass analysis_id | Description (Genebass) | n (metadata) | Alternatives not used |
|---|---|---|---|---|
| fluid_intelligence | `continuous-20016-both_sexes--irnt` | Fluid intelligence score | 128,302 | 20191, `Fluid_intelligence_score_custom` (100,900) |
| reaction_time | `continuous-20023-both_sexes--irnt` | Mean time to correctly identify matches | 392,194 | |
| numeric_memory | ABSENT | field 4282 not among the 4,529 analyses | | |
| pairs_matching_errors | ABSENT | field 399 not among the 4,529 analyses | | |
| education_years | ABSENT | fields 845 / 6138 / 22501 not present | | |
| hand_grip_strength | `continuous-hand_grip_strength_custom-both_sexes--custom` | Hand grip strength (dominant hand, Genebass-constructed) | 393,068 | 46 (left), 47 (right), `HGSadjHtWt_custom` |
| fev1 | `continuous-3063-both_sexes--irnt` | FEV1 | 360,076 | 20150 (best measure), 20153, 20154 |
| walking_pace | ABSENT | field 924 not present | | |
| resting_heart_rate | `continuous-102-both_sexes--irnt` | Pulse rate, automated reading | 372,555 | 4194 |
| systolic_bp | `continuous-4080-both_sexes--irnt` | Systolic BP, automated reading | 372,551 | `systolic_blood_pressure_custom` (389,886), `..._adjBMI_custom` |
| ldl | `continuous-30780-both_sexes--irnt` | LDL direct | 376,106 | |
| bmi | `continuous-21001-both_sexes--irnt` | BMI | 393,562 | 23104, `BMI_custom` |
| parental_lifespan | `continuous-1807-both_sexes--irnt` + `continuous-3526-both_sexes--irnt` | Father's / Mother's age at death | 291,196 / 234,028 | no single field; combined per gene by fixed-effect inverse-variance meta-analysis (assumes independence; the two share carriers, so it is mildly anti-conservative; Tier B ceiling applies to this trait anyway) |
| triglycerides (control) | `continuous-30870-both_sexes--irnt` | Triglycerides | 376,509 | |
| type_2_diabetes | `categorical-T2D_custom-both_sexes--custom` | Type 2 Diabetes (T2D) | 20,813 / 371,952 | E11 first occurrence |
| coronary_disease | `categorical-CAD_custom-both_sexes--custom` | Coronary artery disease (CAD) | 13,168 / 363,361 | I25, I21 first occurrence |
| cancer_any | `categorical-2453-both_sexes--` | Cancer diagnosed by doctor (self-report) | 30,948 / 362,566 | no registry-based any-cancer analysis found; Z85 (history of malignant neoplasm) not used |
| dementia | `icd_first_occurrence-130842-both_sexes--` | F03 unspecified dementia | 1,091 / 393,750 | F00, F01, F02, G30, `Alzheimers_custom1/2` |
| major_depression | `icd_first_occurrence-130894-both_sexes--` | F32 depressive episode | 41,680 | F33 (2,893), `Depression_custom` (68,114 cases, no documented definition), 20126 codes |
| schizophrenia | `icd_first_occurrence-130874-both_sexes--` | F20 schizophrenia | 792 | `Schizophrenia_custom` (802, no documented definition) |
| all_cause_mortality | ABSENT | only "Age at death" (40007, n=11,391, death register), not all-cause mortality; not substituted | | |
| fracture | `categorical-2463-both_sexes--` | Fractured/broken bones in last 5 years (self-report) | 37,911 / 354,806 | M80/M81 are osteoporosis, not used |
| infertility | `icd_first_occurrence-132156-both_sexes--` (N97 female) + `icd_first_occurrence-132084-both_sexes--` (N46 male) | female / male infertility | 3,044 / 427 | sexes are disjoint; combined per gene by fixed-effect inverse-variance meta-analysis |

## 3. FinnGen R13 (independent replication)

Independence: FinnGen (Finnish biobanks, 500,186 participants in DF13) shares no participants with UK Biobank.

Access: the FinnGen docs (https://finngen.gitbook.io/documentation/data-download) publish the manifest URL and ask users to fill an online form to receive download instructions by email. The bucket itself is anonymously readable (no login, no requester-pays): object listing, HEAD and range GET succeed without credentials. The form was NOT submitted by this pipeline. Human action for compliance: submit https://elomake.helsinki.fi/lomakkeet/124935/lomake.html.

Files (bucket `finngen-public-data-r13`, prefix `lof/`; listed 2026-09-29):

| Object | Bytes | GCS md5 (base64) |
|---|---|---|
| `lof/data/finngen_R13_lof.txt.gz` | 548,499,199 | `J3vgv3ucUTtPJ5XhsRQgVg==` |
| `lof/data/finngen_R13_lof_variants.txt` | 162,231 | `KHTmo0Sxxt1gh3v32U/WfQ==` |
| `lof/finngen_R13_lof_readme` | 1,435 | `QInxG/+LQPcFS9miTgs7Kw==` |
| `summary_stats/finngen_R13_manifest.tsv` | 820,805 | `pw9Oy1es98I+j11xvhHZNg==` |

Method (readme + docs): regenie step 2, burden mode, max mask over LoF variants; LoF = frameshift, splice donor, splice acceptor, stop gained (VEP, no LOFTEE); MAF <= 0.01; info >= 0.8; 4,909 autosomal genes; core binary endpoints only (PD_DEMENTIA_EXMORE removed); LoF from imputed genotypes. PLoF mask only: no `dmis`, no `syn`.

### 3.1 Declared proxies (Amendment 1), exact endpoint codes (verified in the R13 manifest and Risteys)

| Panel trait | Normalized replication trait | FinnGen endpoint | Long name | Cases / controls (R13 manifest) | Gene-based LoF results downloadable free |
|---|---|---|---|---|---|
| ldl | `hypercholesterolemia` | `E4_HYPERCHOL` | Pure hypercholesterolaemia (ICD-10 E78.0) | 48,649 / 411,530 | yes, if the endpoint is a column/phenotype in `finngen_R13_lof.txt.gz` (confirmed in section 6) |
| bmi | `obesity` | `E4_OBESITY` | Obesity (ICD-10 E66) | 33,933 / 466,093 | as above |
| systolic_bp | `hypertension` | `I9_HYPTENS` | Hypertension (ICD-10 I10-I15, I67.4) | 162,639 / 337,464 | as above |

Same-trait rows for ldl/bmi/systolic_bp: not available (FinnGen LoF burden is binary endpoints only).

### 3.2 Trade-off endpoints in FinnGen (fixed now; supplementary, replication cohort rows)

| Normalized trait | FinnGen endpoint | Note |
|---|---|---|
| type_2_diabetes | `T2D` | definitions combined |
| coronary_disease | `I9_CHD` | major coronary heart disease event |
| cancer_any | `C3_CANCER` | malignant neoplasm |
| dementia | `F5_DEMENTIA` | |
| major_depression | `F5_DEPRESSIO` | Depression (broad; recurrent-only `F5_DEPRESSION_RECURRENT` not used) |
| schizophrenia | `F5_SCHZPHR` | |
| all_cause_mortality | ABSENT | no all-cause death endpoint in the manifest |
| fracture | ABSENT | only site-specific `ST19_FRACT_*` endpoints; none is "any fracture"; not substituted |
| infertility | `N14_FEMALEINFERT` + `N14_MALEINFERT` | sex-specific; fixed-effect inverse-variance meta-analysis, as in Genebass |

## 4. UK Biobank overlap across sources

| Pair | Overlap | Consequence |
|---|---|---|
| Genebass vs AstraZeneca vs Regeneron (incl. as integrated in Open Targets) | same UK Biobank exomes | one cohort; never counted as replication of each other |
| Genebass vs FinnGen R13 | none (Finnish biobanks; no UKB samples) | valid independent replication |
| Genebass vs gnomAD v4.1 constraint | gnomAD v4 exomes include UKB participants (recalled from release notes, not re-verified) | constraint is annotation only (LOEUF context in report), never used for testing |

## 5. Not obtained

- Genebass Hail tables: requester-pays (billing project required; not authorized). Route in section 2 used instead.
- AstraZeneca PheWAS gene-level pLoF: not in bulk; portal only. Consistency check not performed.
- All of Us All by All: Researcher Workbench, gated. Not pursued.
- Numeric memory, pairs matching errors, education years, walking pace: not in Genebass. All-cause mortality: not in Genebass or FinnGen. Fracture: not in FinnGen.
