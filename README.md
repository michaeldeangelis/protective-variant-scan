# Protective variant scan

I am testing whether people who carry a broken copy of one gene do better on brain or body measures without paying for it somewhere else. I wrote the rules down first. The scan has not been run on real data yet, so there is no result.

## Abstract

Most traits, including reaction time, height and blood pressure, are shaped by thousands of small genetic differences. Each one adds a tiny push, so there is rarely one gene to switch. The best exception is a rare broken copy of a single gene that helps. The standard example is PCSK9. People with a broken copy have 15 to 28 percent lower LDL cholesterol and much less heart disease, though a broken copy also goes with a higher risk of diabetes. I mapped what published studies say for cognitive and physical traits in [architecture-map.md](architecture-map.md). Those figures come from abstracts and search summaries, not full texts. For years of education, about 3,950 variants explain 11 to 16 percent of the differences between people. Comparing siblings shrinks the effect by about 46 percent, so part of that signal is family environment and not biology.

Then I wrote a test before opening any results. It uses public summary tables from about 395,000 UK Biobank exomes (Genebass) and a separate Finnish cohort (FinnGen). For each gene it asks three things. Do carriers of a broken copy do better on 13 named traits, such as fluid intelligence, reaction time, grip strength, lung function, resting heart rate, blood pressure, LDL, BMI and parental lifespan? Does the finding hold in the independent cohort? Does the gene avoid harm on nine health outcomes? Genes that pass go into tiers. Three known genes (PCSK9, ANGPTL4, APOC3) and fake-data checks must come out right, or the run stops and reports a failure. A "lead" means a gene that passed in the first cohort but could not be checked in a second one.

Two limits shaped the design. No independent cohort exists for the continuous brain measures, so those results can reach "lead" at most and never a full pass. Four of the 13 planned traits (numeric memory, pairs matching, years of education, walking pace) and one of the nine health outcomes (all-cause mortality) are not in the public tables, and I did not substitute them. Blood pressure is the only body measure that can produce a full pass.

What exists now is the pipeline, more than 470 passing tests, and an independent reviewer who broke the code by hand in more than 40 ways and had every break caught. What does not exist is a real result. The rules were amended (Amendment 1 and clarifications C1 to C11 in [experiments.md](experiments.md)), each dated before any result was opened. One exception is disclosed there: two FinnGen rows were printed while debugging a file format. Nothing here is about causes, treatment or editing people.

The question and what would count as an answer are in [hypothesis.md](hypothesis.md). The decision rules are in [experiments.md](experiments.md). Data sources, checksums and every normalization choice are in [docs/data-sources.md](docs/data-sources.md).

## Reproduce

The pipeline runs on invented data with planted signals. It needs Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
pytest
python scripts/make_synthetic.py --scenario pass --out /tmp/synth
python -m protscan run --config config/prereg.yaml --data /tmp/synth --out /tmp/synth/result.json
```

The last command prints the verdict and writes `result.json` and `report.md` beside it. Scenarios are `pass`, `lead`, `nolead`, `unscreened` and four with a deliberately broken control. A run on synthetic data is labeled SYNTHETIC in the verdict.

The real tables are fetched by `scripts/fetch_genebass.py`, `scripts/fetch_finngen.py` and `scripts/fetch_gnomad.py` into `data/`, which is not committed. About 930 MB. The config hash is pinned, so a run with an edited `config/prereg.yaml` is labeled NON-PREREGISTERED.

## What else is here

`architecture-map.md` is the sourced map of how these traits are built. `cognition.md`, `body-healthspan.md` and `methods-limits.md` hold the figures behind it, with UNVERIFIED marked where a number was not retrieved.

`src/protscan/` is the pipeline: adapters to a common table, the statistics, the controls, the tiering and the report.

`tests/test_ledger_conformance.py` is the reviewer's independent copy of every rule, typed from the ledger and not read from the config. `reviews/` holds the review notes and sign-off.

`V1.md` is the build plan with each decision and how to flip it.

## Acknowledgements and data licenses

- We want to acknowledge the participants and investigators of the FinnGen study. Kurki, M.I., Karjalainen, J., Palta, P. et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature 613, 508-518 (2023). https://doi.org/10.1038/s41586-022-05473-8
- Genebass (UK Biobank exome gene-based burden results) is used under CC BY 4.0 through its public API.
- gnomAD v4.1 constraint metrics are used under the gnomAD terms.
- No raw data is committed.

## Cite

Use the Cite this repository button. It reads `CITATION.cff`.

## License

[MIT](LICENSE). Copyright 2026 Michael DeAngelis.
