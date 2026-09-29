# Polygenic architecture: conceptual model, limits, data, tools

Scope: public-data analysis only. Verification level: claims below were checked against publisher or search-result abstracts and site pages (not full-text reads). Anything not confirmed is marked UNVERIFIED.

## 1. Conceptual model: omnigenic / core genes
- Boyle, Li, Pritchard (2017, Cell 169:1177) propose that gene regulatory networks are interconnected enough that most genes expressed in a disease-relevant cell can affect "core" genes; most heritability then comes from peripheral genes acting indirectly, not from core genes. [S1]
- Liu, Li, Pritchard (2019, Cell 177:1022) give a formal model: most heritability driven by weak trans-eQTL SNPs whose effects are mediated through peripheral genes onto core genes. [S2]
- Critiques: Wray et al. (2018, Cell) argue the assumption of a small set of core genes per disease may understate biological complexity. [S3] A later commentary argues for a gradient rather than a sharp core/peripheral boundary (surfaced via search; exact citation UNVERIFIED). Sella and Barton (2019, Annu Rev Genomics Hum Genet 20:461) is the evolutionary-genetics treatment of why effects are spread thin (selection on complex traits); I did not read it in detail. [S4]
- Implication for "single-gene levers": under this model, core genes are the only places where a single-gene change is expected to matter strongly, and they carry a minority of heritability. For most traits, expect many small effects and few levers. Exceptions come from coding/rare variants in core genes (see Sec. 2).

## 2. Missing heritability
- Wainschtein et al. (2022, Nat Genet 54:263): WGS of 25,465 European-ancestry people; estimated heritability 0.68 (SE 0.10) for height and 0.30 (SE 0.10) for BMI. Low-MAF variants in low-LD regions were enriched for heritability, more so for protein-altering variants (consistent with negative selection). Rare variants are a major source of the remaining missing heritability. [S5]
- Read: array-based SNP heritability underestimates because rare variants are poorly tagged; large WGS/exome cohorts are needed to see them.

## 3. Interpretation limits
- Direct vs indirect effects: Howe et al. (2022, Nat Genet 54:581) show within-sibship GWAS estimates are smaller than population estimates for seven phenotypes including height, educational attainment, and cognitive ability, reflecting demographic and indirect genetic effects. [S6] Okbay et al. (2022, Nat Genet 54:437) is the 3-million-person educational attainment GWAS, with within- vs between-family prediction analysis. [S7] Consequence: population GWAS and polygenic scores for cognitive traits overstate direct causal effect.
- Ancestry bias: Martin et al. (2019, Nat Genet 51:584) report roughly 79% of GWAS participants are of European ancestry (figure as given in a search-result summary; verify in paper), and PRS trained in Europeans predict less well elsewhere. [S8]
- Also real but not separately sourced here: LD blurs which variant is causal (fine-mapping addresses this, Sec. 6), assortative mating inflates estimates (partly covered by Howe 2022). Variant-to-gene mapping is probabilistic; Open Targets' L2G is a gradient-boosting model over fine-mapping, QTL colocalisation, and functional genomics features. [S9]

## 4. Regulatory angle
- Most GWAS hits are non-coding; sequence-to-function models are used to rank candidate causal regulatory variants.
- Enformer (Avsec et al. 2021, Nat Methods) integrates interactions up to about 100 kb. [S10] Karollus et al. (2023, Genome Biol) find such models capture promoter determinants but mostly ignore distal enhancers. [S11]
- Personal-genome limit: Huang et al. (2023, bioRxiv) find Enformer captures variation across genes but struggles across individuals; fine-tuning reaches roughly variant-based linear model performance. [S12]
- eQTL sign prediction: AlphaGenome (Nature 2026; bioRxiv 2025) reports mean sign auROC rising from 0.75 (Borzoi) to 0.80; Enformer and Borzoi struggle with direction for distal eQTLs. [S13, S14] Figures are from search-result summaries; UNVERIFIED against paper text.
- Takeaway: use these models to prioritize and generate hypotheses, not as ground truth for a person's variant effect.

## 5. Public data (all summary-level; open unless noted)
| Resource | Contents | Access | URL |
|---|---|---|---|
| GWAS Catalog | Curated published associations plus full summary stats; counts not shown on fetched page | Open; download, REST API | https://www.ebi.ac.uk/gwas/ |
| Open Targets | Fine-mapping credible sets, L2G gene ranking, QTL colocalisation | Open (API/bulk expected; page fetch failed, UNVERIFIED) | https://platform.opentargets.org/ |
| Pan-UK Biobank | 7,228 phenotypes, 6 ancestry groups, 16,131 GWAS | Open summary stats download | https://pan.ukbb.broadinstitute.org/ |
| FinnGen | Release DF13: 500,186 people, 2,755 endpoints, 21.3M variants | Summary stats free; online form for download instructions; individual data via Fingenious | https://www.finngen.fi/en/access_results |
| GTEx v8 | 15,201 RNA-seq samples, 49 tissues, 838 donors; eQTL/sQTL | Portal open for summary; raw individual data via dbGaP controlled access | https://gtexportal.org |
| gnomAD v4 | 807,162 individuals (730,947 exomes, 76,215 genomes); allele frequencies, constraint | Open | https://gnomad.broadinstitute.org/ |
| PGC | Psychiatric GWAS summary stats (ADHD, MDD, SCZ, etc.) | Download with terms: no redistribution, research only, commercial use needs permission | https://pgc.unc.edu/for-researchers/download-results/ |
| PGS Catalog | 6,991 polygenic scores, 813 traits, 836 publications (as of 2026-09-17) | Open; web, FTP, REST, pgsc_calc | https://www.pgscatalog.org/ |
| SSGAC | Educational attainment (Okbay 2022) and cognitive-performance summary stats | Page fetch failed (redirect); access terms UNVERIFIED | https://www.thessgac.org/ |

## 6. Toolchain (cited names only; full-text not read)
- LDSC: Bulik-Sullivan et al. 2015 Nat Genet, separates confounding from polygenicity, estimates heritability and genetic correlation. [S15]
- MAGMA: de Leeuw et al. 2015 PLoS Comput Biol, gene and gene-set analysis. [S15]
- coloc: Giambartolomei et al. 2014 PLoS Genet, colocalisation of GWAS and QTL signals; newer versions use SuSiE. [S15]
- TwoSampleMR / MR-Base: Hemani et al. 2018 eLife, Mendelian randomization. [S15]
- SuSiE (Wang et al. 2020, fine-mapping) and PRS-CS (Ge et al. 2019) / LDpred (Vilhjalmsson 2015, polygenic scoring): named from memory, not confirmed in this pass. UNVERIFIED citations.

## Sources
- S1 https://doi.org/10.1016/j.cell.2017.05.038
- S2 https://www.cell.com/cell/fulltext/S0092-8674(19)30400-3
- S3 https://pubmed.ncbi.nlm.nih.gov/29906445/
- S4 https://www.annualreviews.org/content/journals/10.1146/annurev-genom-083115-022316
- S5 https://www.nature.com/articles/s41588-021-00997-7
- S6 https://www.nature.com/articles/s41588-022-01062-7
- S7 https://www.nature.com/articles/s41588-022-01016-z
- S8 https://www.nature.com/articles/s41588-019-0379-x
- S9 https://www.nature.com/articles/s41588-021-00945-5
- S10 https://www.nature.com/articles/s41592-021-01252-x
- S11 https://link.springer.com/article/10.1186/s13059-023-02899-9
- S12 https://www.biorxiv.org/content/10.1101/2023.03.16.532969v2.full
- S13 https://www.nature.com/articles/s41586-025-10014-0
- S14 https://arxiv.org/html/2411.11158v2
- S15 tool papers surfaced via search: https://pmc.ncbi.nlm.nih.gov/articles/PMC8168525/ (Computational Tools for Causal Inference in Genetics)
- GTEx v8: https://www.science.org/doi/10.1126/science.aaz1776
- gnomAD v4: https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/
