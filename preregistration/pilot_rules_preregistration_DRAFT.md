# INDAGE pilot study — preregistration of decision rules (DRAFT)

*Draft prepared with the proposal (November 2026). To be completed and registered on OSF Registries (template: "OSF Preregistration") before the pilot starts (project month 6). The analysis plan of SO2–SO3 is registered separately before nationwide data collection closes.*

## 1. Study information

- **Title.** INDAGE pilot study: selection of cognitive tests and final version of a two-part healthy ageing scale for adults with intellectual disability.
- **Authors.** Maria Carbó Carreté, Albert Sánchez-Niubó (Universitat de Barcelona).
- **Research questions.** (1) Which NeuronUP Assessment tests are feasible and valid enough in adults with mild or moderate ID to be carried into the nationwide study? (2) Is self-report of functional ability sufficiently consistent with an informant's report in adults with moderate ID? (3) Which items are retained, revised or removed?

## 2. Design plan

- Cross-sectional pilot in Catalonia, recruited through Dincat organisations.
- All participants complete the full scale (intrinsic capacity and functional ability), the context questionnaire and the CAMCOG-DS; a professional who has known them for at least six months completes the functional ability items independently.
- NeuronUP Assessment under a planned missing design: every participant completes a common memory test and a subset of the remaining tests, assigned in balanced combinations so that each test is completed by half of the sample and every pair of tests is observed together.
- Retest of the report items in 50 participants after two to three weeks.

## 3. Sampling plan

- **n = 150** adults aged 35 or over with mild or moderate ID (certified degree ≥ 33%), stratified by Down syndrome, sex and age band (35–49, 50–64, 65+), with **moderate ID oversampled to half** of the sample.
- Justification: each NeuronUP test (completed by ~75 participants) is evaluated with ±.08 for the proportion at floor and a lower 95% limit above .30 for a correlation of .50 with the CAMCOG-DS; 50 retest participants estimate an ICC of .80 within ±.10 (Bonett, 2002). See `R/06_pilot_precision.R`.
- Pilot participants do not take part in the nationwide study.

## 4. Variables

- NeuronUP test scores, completion status and administration time; CAMCOG-DS total and subscale scores.
- Item responses of both parts of the scale; informant responses to functional ability items; level of assistance provided during administration.
- Level of ID, aetiology (DS / other known / unspecified), sex, age.

## 5. Decision rules (fixed before data collection)

**NeuronUP tests.** A test is retained for the nationwide study when all of the following hold:

| Rule | Threshold |
|---|---|
| Completion | ≥ 80% of participants who received it complete it |
| Floor effect | ≤ 15% obtain the minimum score (Terwee et al., 2007) |
| Missing data | ≤ 10% |
| Convergent validity | correlation ≥ .50 with the corresponding CAMCOG-DS score |
| Duration | within the administration time agreed with the participating organisations [to specify, minutes] |

- A floor above 15% **within moderate ID** flags the test for joint review by both PIs (alert, not automatic exclusion).
- If the retained tests do not cover **memory and executive function**, the CAMCOG-DS is administered to all participants in the nationwide study.

**Informant report.** If the upper limit of the 95% confidence interval of person–informant agreement (ICC) on functional ability in adults with moderate ID is below .75, the informant report is retained in the nationwide study.

**Items.** Items are revised or removed on item non-response, floor and ceiling effects by level of ID and test–retest reliability [thresholds to specify]. Content validity criteria from the expert panel (Aiken's V, lower 95% limit ≥ .70 retained; .50–.70 revised) are applied before the pilot.

## 6. Analysis plan (pilot)

- Proportions with 95% Wald intervals; correlations with Fisher-z intervals; ICC (two-way, absolute agreement) with 95% intervals.
- Planned missing data handled by pairwise-available estimation for test-level statistics.
- All decisions documented in a public report archived in Zenodo.

## References

- Bonett, D. G. (2002). *Statistics in Medicine, 21*(9), 1331–1335.
- Graham, J. W., et al. (2006). *Psychological Methods, 11*(4), 323–343.
- Terwee, C. B., et al. (2007). *Journal of Clinical Epidemiology, 60*(1), 34–42.
