# INDAGE — sample size simulation and pilot precision

Code, results and report supporting the sample size and sampling allocation of the proposal
**"Age-related indicators of functional decline in adults with intellectual disability: development and validation of a healthy ageing scale (INDAGE)"**
(Spanish State Research Agency, Knowledge Generation Projects 2026; PIs: Maria Carbó Carreté and Albert Sánchez-Niubó, Universitat de Barcelona).

**Report:** [`report/indage_sample_size_report.html`](report/indage_sample_size_report.html) — download and open it in a browser.

## Contents

| Folder | Content |
|---|---|
| `R/` | Simulation scripts. `00_assumptions.R` holds every planning assumption and the adopted design; `01`–`04` and `06` one criterion each; `05_allocation.R` and `07_design_grid.R` the allocation and the choice of the sample size; `run_all.R` runs everything. |
| `results/` | CSV outputs used by the report, and `sessionInfo.txt`. |
| `data/` | Sampling frame figures (BEDPD 2023), Table 1 of this proposal and the 2025 allocation it is compared with (see `data/README.md`). |
| `report/` | R Markdown report and its rendered HTML. |
| `preregistration/` | Draft of the pilot decision rules, to be registered on OSF Registries before the pilot. |

## Design

1,000 participants: 700 adults with ID of aetiologies other than Down syndrome (DS) and a national minimum of 300 with DS (115 aged 35–44, 115 aged 45–54, 70 aged 55–64). The two strata are modelled separately, so their sizes are set separately. Regional totals are proportional to the certified-disability population by autonomous community with a floor of 15; the age profile keeps the 35–44 reference large, because every onset age is defined relative to its mean, and supports the 65+ end, because natural splines are linear beyond the boundary knots. Pilot: 150 adults in Catalonia, retest 50.

## Criteria

1. **Item parameters** — RMSE of graded response model parameters per domain.
2. **Factor structure** — logic of SENECA (Lorenzo-Seva & Ferrando, 2024) on an assumed population model; to be recomputed with the first 250 nationwide participants.
3. **Measurement invariance** — power to detect uniform DIF across age bands and between adults with and without Down syndrome.
4. **Onset age (binding)** — precision of the age from which domain scores fall 0.5 SD below the 35–44 reference, by the size of the Down syndrome sample. `07_design_grid.R` shows that the two strata move independently and that 700 + 300 is the smallest design meeting both acceptance criteria under the conservative scenario; `05_allocation.R` compares allocations at a fixed size.
5. **Pilot precision** — precision of the pilot decision rules by pilot size (`06`).

## Reproduce

Requirements: R ≥ 4.3 with `mirt`, `lavaan`, `MASS`, `splines`, `parallel` (simulation) and `rmarkdown`, `knitr`, `ggplot2` with pandoc (report). Versions used are in `results/sessionInfo.txt`. Parallel code uses `parallel::mclapply` (Linux or macOS).

```bash
# from the repository root
Rscript R/run_all.R                      # about 35 minutes on 7 cores; writes results/
Rscript -e 'rmarkdown::render("report/indage_sample_size_report.Rmd", knit_root_dir = getwd())'
SMOKE=1 Rscript R/04_onset_age.R         # quick test with few replications
```

Seeds are fixed in each script; rerunning reproduces the CSV files in `results/`.

## License and citation

Code: MIT (see `LICENSE`). Report and data compilation: CC BY 4.0. Please cite the archived version (see `CITATION.cff` and the Zenodo DOI).
