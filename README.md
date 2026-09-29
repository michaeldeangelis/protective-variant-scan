# genetics-research

Public-data analysis of the genetic architecture of cognitive and physical traits, and a preregistered scan for protective loss-of-function variants.

- `architecture-map.md`: sourced map of polygenic architecture (sources in `cognition.md`, `body-healthspan.md`, `methods-limits.md`).
- `experiments.md`: append-only preregistration ledger. Rules are fixed there before any run.
- `V1.md`: build plan and status for the protective-variant scan.

Analysis of public summary data only. No individual-level data, no protocols for modifying human genomes. Figures in the docs come from abstracts and search summaries unless marked otherwise; UNVERIFIED means not retrieved.

## Acknowledgements and data licenses
- We want to acknowledge the participants and investigators of the FinnGen study. Kurki, M.I., Karjalainen, J., Palta, P. et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 613, 508-518 (2023). https://doi.org/10.1038/s41586-022-05473-8
- Genebass (UK Biobank exome gene-based burden results) is used under CC BY 4.0 via its public API; see docs/data-sources.md for URLs, checksums and dates.
- gnomAD v4.1 constraint metrics, used per the gnomAD terms; see docs/data-sources.md.
- No raw data is committed; `data/` is gitignored.
