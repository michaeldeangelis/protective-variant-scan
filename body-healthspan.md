# Polygenic architecture: physical function and healthspan

Scope: public published data, analysis only. Compiled 2026-09-29.
Verification level: figures below come from search-result summaries of the cited papers. Nature.com full text was blocked (login redirect), so no figure was checked against the paper body. Anything not seen in a summary is marked UNVERIFIED.

## 1. Per-trait table

| Trait | Twin / pedigree h2 | Largest GWAS found | Independent GWS loci | Variance explained / PGS R2 |
|---|---|---|---|---|
| Height (benchmark) | 54.3% (Polderman 2015 meta-analysis, 87 studies); ~80% is often quoted, UNVERIFIED here | Yengo 2022, N=5.4M, Nature | 12,111 SNPs in 7,209 loci | ~40% in European ancestry (R2 ~0.40-0.44); ~0.13-0.20 in a Latino/Hispanic group; >54% when combined with parental average height |
| Grip strength / weakness | UNVERIFIED | Jones 2021, N=256,523 (age 60+, 22 cohorts), Nat Commun | 15 (weakness, EWGSOP definition; 46,596 cases, 18.9%) | UNVERIFIED |
| VO2max (measured) | 44-68% range; meta-analysis of twin/sibling studies 59% (95% CI 52-66%) | Klevjer 2022, ~4,500 people | 2 in total sample (CDYL, LOC105371536), plus sex-specific hits; very underpowered | UNVERIFIED |
| Lifespan / longevity | ~20-25% (twin), as low as 7% (pedigree); Ruby 2018 says <10% after assortative-mating correction; Shenhar 2025 (Science) says ~50% for intrinsic lifespan when extrinsic mortality is addressed. Estimates conflict | Timmers 2019, 1,012,240 parental lifespans, eLife; Deelen 2019 (~10,000 long-lived, Nat Commun) | 12 (Timmers, common-effect model); Deelen: APOE (rs7412) plus a locus near GPR78 | Not found |
| BMI | UNVERIFIED | Yengo 2018, ~700,000, Hum Mol Genet | 941 near-independent SNPs | 5% of phenotypic variance |
| Blood pressure | UNVERIFIED | Evangelou 2018, >1M, Nat Genet | 901 loci (535 novel) | SBP 11.2% for all 901 loci (4.6% for the 274 earlier loci) |
| LDL / lipids | UNVERIFIED | Graham 2021, ~1.65M (350,000 non-European), Nature | 941 lipid loci (355 new) across lipid traits | LDL PGS R2 not retrieved |
| Sleep duration | UNVERIFIED | Dashti 2019, N=446,118, Nat Commun | 78 (p<5e-8) | SNP h2 9.8%; 78 loci explain 0.69% |

Correction to my earlier premise: there is no "Jones 2019, 256 loci" grip GWAS. The 256,523 figure is the sample size of Jones 2021, which found 15 loci.

## 2. Few-large-effect versus diffuse

- Diffuse: height (12,111 SNPs, and the paper calls the map "saturated" for European ancestry), BMI (941 SNPs explain only 5%), blood pressure (901 loci, 11% of SBP variance), sleep duration (78 loci explain 0.69%). Typical per-variant effect is tiny; the Timmers lifespan alleles are the only per-variant effects I could quote in natural units.
- Lifespan effects per allele (Timmers 2019): APOE -1.06 years, CHRNA3/5 -0.42, LPA -0.76, HLA-DQA1 +0.56 per minor allele; life-extending alleles range +0.23 to +1.07 years. Only APOE is robust in longevity-style GWAS (Deelen 2019).
- Concentrated: LDL via PCSK9, ANGPTL3/4, APOC3. These are lipid-pathway genes with large effects from rare coding variants; lipids are the exception, not the rule.

## 3. Rare protective variants

| Variant | Effect | Evidence | Trade-off |
|---|---|---|---|
| PCSK9 loss-of-function | Cohen 2006 (ARIC): nonsense variants in 2.6% of 3,363 Black participants, LDL -28%, CHD -88%; sequence variants in 3.2% of 9,524 White participants, LDL -15%, CHD -47%. Dewey 2016: LDL -35 mg/dL (African American) and -13 mg/dL (White) | Multiple cohorts, consistent | Schmidt 2017 (Mendelian randomization): LDL-lowering PCSK9 variants associate with higher fasting glucose, weight, waist-to-hip ratio, and type 2 diabetes risk |
| ANGPTL4 E40K | Dewey 2017 NEJM: CAD OR 0.81 (0.70-0.92, P=0.002); triglycerides -13%, HDL +7%; n=1,661 heterozygotes | Single large exome study | None found in the summary |
| ANGPTL3, APOC3 loss-of-function | Reduced triglycerides and CHD risk; ANGPTL3 heterozygotes ~1 in 300 | Multiple studies | Not retrieved; UNVERIFIED effect sizes |
| DEC2/BHLHE41 P384R | Carriers sleep ~6 h instead of ~8 h, without usual sleep-loss effects (He 2009, Science); mechanism: higher orexin expression | One family (mother and daughter) in the original paper; later fly and mouse work (bioRxiv 2023, iScience 2025 reported) | Human trade-offs not established; evidence is too thin to call it a general lever |

## 4. Cognition versus physical and health traits

- Intelligence and parental longevity: rg = 0.33 (SE 0.08) using the intelligence GWAS, 0.37 (SE 0.07) using their meta-analytic sample (Hill 2019, Mol Psychiatry).
- Davies 2018 (N=300,486, 148 loci): overlap with processing speed and health variables including longevity. Direction: higher cognition tracks better health and longevity. The correlations mix genuine pleiotropy with socioeconomic confounding.
- No evidence retrieved of a cognition versus physical-function trade-off. The trade-offs found were within metabolic biology (PCSK9 lowers LDL but raises diabetes risk), not between brain and body.

## 5. Caveats

- Within-family (sibling) GWAS shrinks estimates for cognitive traits and can reduce genetic correlations (Howe 2022, Nat Genet, 178,076 siblings, 25 phenotypes). Attenuation magnitudes for height, BMI and physical traits were not retrieved: UNVERIFIED.
- Portability: height PGS R2 falls from ~0.40 to ~0.13-0.20 outside European ancestry.
- Lifespan heritability is unsettled, roughly <10% to ~50% depending on method, so the ceiling for polygenic prediction of longevity is unclear.
- Physical-performance GWAS (VO2max, grip) are far smaller than height or BMI, so their locus counts reflect power, not architecture.

## 6. Gaps to close

Twin h2 for grip strength, BMI, BP, sleep, LDL; SNP h2 and PGS R2 for grip, VO2max, lifespan, LDL; within-family attenuation for physical traits; DEC2 replication in larger cohorts; ANGPTL3/APOC3 effect sizes; elite-athlete GWAS (not searched).

## Sources

- Yengo 2022, Nature, "A saturated map of common genetic variants associated with human height": https://www.nature.com/articles/s41586-022-05275-y ; PGS Catalog PGP000382: https://www.pgscatalog.org/publication/PGP000382/
- Jones 2021, Nat Commun, muscle weakness: https://www.nature.com/articles/s41467-021-20918-w
- Klevjer 2022 (VO2max GWAS, via search summary); twin meta-analysis PMC4773888: https://pmc.ncbi.nlm.nih.gov/articles/PMC4773888/
- Timmers 2019, eLife: https://elifesciences.org/articles/39856
- Deelen 2019, Nat Commun: https://www.nature.com/articles/s41467-019-11558-2
- Polderman 2015, Nat Genet: https://gwern.net/doc/genetics/heritable/2015-polderman.pdf
- Ruby 2018, Genetics (assortative mating); Shenhar 2025, Science: https://www.science.org/doi/10.1126/science.adz1187
- Yengo 2018 (BMI, 941 SNPs): https://www.ukbiobank.ac.uk/publications/meta-analysis-of-genome-wide-association-studies-for-height-and-body-mass-index-in-700000-individuals-of-european-ancestry/
- Evangelou 2018, Nat Genet (BP): https://www.nature.com/articles/s41588-018-0205-x
- Graham 2021, Nature (lipids): https://www.nature.com/articles/s41586-021-04064-3
- Dashti 2019, Nat Commun (sleep duration): https://www.nature.com/articles/s41467-019-08917-4
- Cohen 2006, NEJM (PCSK9): https://www.researchgate.net/publication/7223799_Sequence_Variations_in_PCSK9_Low_LDL_and_Protection_against_Coronary_Heart_Disease
- Dewey 2016 PCSK9 (Circ Genom Precis Med): https://pmc.ncbi.nlm.nih.gov/articles/PMC5729040/
- Dewey 2017, NEJM (ANGPTL4): https://www.nejm.org/doi/full/10.1056/NEJMoa1510926
- Schmidt 2017, Lancet Diabetes Endocrinol (PCSK9 and T2D): https://pubmed.ncbi.nlm.nih.gov/27908689/
- He 2009, Science (DEC2): via https://www.ncbi.nlm.nih.gov/clinvar/RCV000004788/ and https://pubmed.ncbi.nlm.nih.gov/40989023/
- Hill 2019, Mol Psychiatry: https://www.nature.com/articles/s41380-017-0001-5
- Davies 2018, Nat Commun: https://pmc.ncbi.nlm.nih.gov/articles/PMC5974083/
- Howe 2022, Nat Genet (within-sibship GWAS): https://www.nature.com/articles/s41588-022-01062-7
