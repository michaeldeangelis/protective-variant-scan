# Polygenic architecture map: cognition and body function

Compiled 2026-09-29 from `cognition.md`, `body-healthspan.md`, `methods-limits.md` (same folder; sources are listed there). Public literature, analysis only. Figures come from abstracts and search summaries, not full-text reads. UNVERIFIED = not retrieved.

## 1. The picture in one table

| Trait | Genetic scale | Variance explained by scores | Notes |
|---|---|---|---|
| Height (benchmark) | 12,111 SNPs, 5.4M people | ~40% (EUR), ~13-20% (Latino/Hispanic group) | Saturated map; the most polygenic trait we understand |
| Educational attainment | 3,952 SNPs, ~3M | 11-16% | Within-sibling effects shrink 46% |
| Intelligence | 205 loci, ~270K (2018) | ~5% | Within-sibling shrinkage 22% (cognitive ability) |
| BMI | 941 SNPs, ~700K | 5% | |
| Blood pressure | 901 loci, >1M | 11.2% (SBP) | |
| Sleep duration | 78 loci, 446K | 0.69% | |
| Lifespan | 12 loci (parental lifespans, 1M) | not found | Heritability disputed: <10% to ~50% |
| Grip strength | 15 loci, 256K | UNVERIFIED | Small GWAS, low power |
| VO2max | 2 loci, ~4.5K | UNVERIFIED | Badly underpowered |
| LDL / lipids | 941 lipid loci, 1.65M | UNVERIFIED | Exception: a few large-effect genes |

## 2. Three tiers of architecture

1. **Diffuse common variants (most of the signal).** Thousands of loci, each a tiny fraction of an SD. Cognition, height, BMI, BP, sleep, longevity. No single gene to target.
2. **Rare coding variants in core genes (large effect, few carriers).**
   - Cognition: Chen 2023 found 8 genes (ADGRB2, KDM5B, GIGYF1, ANKRD12, SLC8A1, RC3H2, CACNA1A, BCAS3), but effect direction is UNVERIFIED. Documented large-effect cognition variants are mostly deleterious, in neurodevelopmental sets. No verified variant raises cognition.
   - Lipids: PCSK9 loss-of-function (LDL -15% to -28%, CHD -47% to -88%), ANGPTL4 E40K (CAD OR 0.81). Lipids are the one domain with clean protective levers.
3. **Trade-offs inside biology.** PCSK9 lowering associates with higher type 2 diabetes risk (Schmidt 2017, MR). DEC2 short sleep rests on one family. No brain-vs-body trade-off found; cognition and parental longevity correlate at rg 0.33-0.37.

## 3. Missing heritability and how to read the numbers

- Twin heritability of adult intelligence is ~60%, against SNP heritability of ~0.13-0.25.
- WGS closes some of the gap: rare variants are a major missing piece (Wainschtein 2022: height 0.68, BMI 0.30). For educational attainment and fluid IQ, the 2025 WGS estimate falls from 48-61% to 38-40% after ancestry and geography adjustment, so part of the twin figure is environmental confounding.
- Population GWAS overstate direct effects for cognition (within-sibling shrinkage of 22-46%, and ~69% of EA genetic variance attributed to indirect effects).
- Scores do not travel: ~79% of GWAS participants are European.

## 4. What this means for the goal ("pull out potential")

- **Genes are not switches for these traits.** Under the omnigenic model, only core genes are expected to work as single levers, and they carry a minority of heritability.
- **The most credible finding paths** are protective rare variants (tier 2 lipids/metabolic) and cell-type and regulatory biology behind the diffuse hits (striatal medium spiny neurons, hippocampal pyramidal neurons for cognition).
- **Turning genes on or off** is a regulation-modeling problem. Sequence models (Enformer, Borzoi, AlphaGenome) can rank candidate causal variants, but they are weak on distal-eQTL direction (AlphaGenome sign auROC 0.80 vs Borzoi 0.75) and across individuals. Use them for hypotheses, not ground truth.

## 5. Data gaps and corrections

- A "Jones 2019, 256 loci" grip GWAS does not exist; the real figures are Jones 2021, N=256,523, 15 loci.
- UNVERIFIED: per-variant effect sizes, Savage 2018 SNP heritability, Okbay within-family R2, Chen 2023 effect sizes and direction, any post-2022 cognition GWAS, twin heritability for most physical traits, grip/VO2max score R2, ANGPTL3/APOC3 effect sizes, SuSiE/PRS-CS/LDpred citations, Open Targets and SSGAC access terms.
- Height twin heritability appears as 54.3% (Polderman 2015) in one file and ~80% in another; unresolved.

## 6. Free data for the next step

GWAS Catalog, Pan-UKB (16,131 GWAS), gnomAD v4, PGS Catalog, GTEx v8 (open summary), FinnGen (form-gated), PGC (terms apply). Toolchain: LDSC, MAGMA, coloc, TwoSampleMR; SuSiE, PRS-CS, LDpred.

## 7. Next steps (user's deferred threads)

1. Cognition and brain route: GWAS to genes to cell types.
2. Body and healthspan route.
3. Protective-variant route (recommended; best signal for least compute).
4. Gene regulation modeling: forward and inverse.
5. Pin down the bio subarea.
