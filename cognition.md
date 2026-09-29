# Polygenic architecture of cognitive traits

Compiled 2026-09-29. Public literature only. Nature.com full text was inaccessible (auth redirect), so figures come from abstracts, search summaries of the papers, and one open-access review. "[2nd]" = taken from a secondary summary, not the primary paper text. UNVERIFIED = not retrieved, not filled from memory.

## 1. Heritability: twin vs SNP vs sequence
| Measure | Value | Source |
|---|---|---|
| Twin, intelligence, early childhood | slightly >20% | Knopik 2017 [2nd, via Genetics of Intelligence review, PMC12416016] |
| Twin, school age | ~40-50% | same |
| Twin, adulthood | up to 60%+ | same |
| SNP h2, general cognitive function | 0.25 | Davies 2018, Nat Commun [2nd, search summary] |
| SNP h2, intelligence (Savage) | UNVERIFIED | - |
| SNP h2, educational attainment (population) | 0.13 (within-sibship 0.04) | Howe 2022 [2nd] |
| WGS-based h2, educational attainment and fluid IQ | 48-61% unadjusted; 40-43% after ancestry PCs; 38-40% after geographic clustering | Wainschtein 2025, Nature [2nd, search summary] |
| Rare variants share of WGS h2, 34 traits | avg 22%; WGS h2 ~ pedigree h2 for 15 traits | Wainschtein 2025 [2nd] |

Gap: twin h2 (~60% adult) vs SNP h2 (~0.13-0.25) is large. Environmental confounding of family-based estimates is part of it; the 2025 WGS result shrinks the apparent heritability substantially once ancestry and geography are adjusted for. Height benchmark: common SNPs ~40%, with rare variants ~68%, twin ~80% (Yengo 2022; Wainschtein 2022 [2nd]).

## 2. Largest GWAS and locus scaling
| Study | N | Loci | Trait |
|---|---|---|---|
| Sniekers 2017 | 78,000 | 18 | intelligence [2nd] |
| Savage 2018, Nat Genet 50:912 | 269,867 | 205 (190 new) | intelligence |
| Davies 2018, Nat Commun | 300,486 | 148 | general cognitive function |
| Lee 2018, Nat Genet | ~1.1M | 1,271 SNPs | educational attainment |
| Okbay 2022, Nat Genet | ~3M | 3,952 SNPs | educational attainment |

Discovery scales roughly with N and has not saturated. Cognition GWAS newer than 2022: none verified. Multi-ancestry MVP work exists (Verma 2024, Science) but no cognition-specific result retrieved.

## 3. Polygenic score performance
- EA: R2 11-13% (Lee 2018 [2nd]); 12-16% depending on validation sample (Okbay 2022 [2nd]); 13.3% (range 7.0-15.8) in US Europeans [2nd]. Cognitive performance score: 7-10% (Lee 2018 [2nd]). Intelligence: ~5.2% (Savage 2018 [2nd, review]).
- Within-family: Howe 2022 (178,086 siblings, 19 cohorts): SNP-effect shrinkage 22% (95% CI 6-37) for cognitive ability, 46% (40-52) for EA. Okbay 2022 within-family PGI R2: UNVERIFIED.
- Indirect effects: 69% of EA genetic variance attributable to indirect effects (Schork 2022 [2nd, review]).
- Cross-ancestry: relative R2 in African vs European ancestry averages 22% and 36% across phenotypes in reports cited by search summary; EA4-specific ratio UNVERIFIED.

## 4. Per-variant effect size
Not retrieved. UNVERIFIED. Qualitatively: thousands of variants and R2 of 5-15% imply each is a tiny fraction of an SD.

## 5. Biology of hits
Savage 2018: genes strongly expressed in brain, specifically striatal medium spiny neurons and hippocampal pyramidal neurons; enrichment in conserved and coding regions (146 nonsynonymous exonic variants); pathways for nervous system development and synaptic structure. Davies 2018: neural and cell developmental pathways. Named individual causal genes: none verified.

## 6. Rare and coding variants
Chen 2023, Nat Genet: exomes of 485,930 adults; 8 genes with large-effect rare coding variants: ADGRB2, KDM5B, GIGYF1, ANKRD12, SLC8A1, RC3H2, CACNA1A, BCAS3. Effect sizes and direction: UNVERIFIED (full text not retrieved; medRxiv returned 403). Related: rare coding variants in SCHEMA schizophrenia genes (29 autosomal genes, FDR<5%) were tested against cognition in 76,783 UK Biobank exomes (medRxiv 2023); result not retrieved. No verified evidence of large-effect variants that raise cognition; the documented large-effect rare variants sit in neurodevelopmental/psychiatric gene sets.

## 7. Confounding
Population estimates for EA and cognition are inflated relative to within-sibship (section 3): consistent with stratification, assortative mating and indirect (parental/environmental) genetic effects. Wainschtein 2025 shows heritability estimates for EA and fluid IQ fall from 48-61% to 38-40% after ancestry and geography adjustment.

## Take-away
Cognition is strongly polygenic: hundreds to thousands of loci, R2 ~5-15%, with a meaningful share of population signal not direct. No single gene to target. Rare large-effect variants are known mainly as deleterious (UNVERIFIED direction for Chen 2023 genes).

## Sources
- Savage et al. 2018, Nat Genet 50:912-919, doi:10.1038/s41588-018-0152-6
- Davies et al. 2018, Nat Commun, https://pmc.ncbi.nlm.nih.gov/articles/PMC5974083/
- Lee et al. 2018, Nat Genet, doi:10.1038/s41588-018-0147-3
- Okbay et al. 2022, Nat Genet, doi:10.1038/s41588-022-01016-z
- Howe et al. 2022, Nat Genet 54:581-592, doi:10.1038/s41588-022-01062-7
- Wainschtein et al. 2022, Nat Genet 54:263-273, doi:10.1038/s41588-021-00997-7
- Wainschtein et al. 2025, Nature, doi:10.1038/s41586-025-09720-6
- Chen et al. 2023, Nat Genet, doi:10.1038/s41588-023-01398-8; preprint https://www.medrxiv.org/content/10.1101/2022.06.24.22276728v1
- Singh et al. 2022 (SCHEMA), Nature, doi:10.1038/s41586-022-04556-w
- Review: The Genetics of Intelligence, https://pmc.ncbi.nlm.nih.gov/articles/PMC12416016/ (cites Knopik 2017, Sniekers 2017, Schork 2022, Yengo 2022)
- Verma et al. 2024, Science (MVP), https://www.science.org/doi/10.1126/science.adj1182
