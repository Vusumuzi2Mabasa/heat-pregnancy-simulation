# Synthetic data dictionary

All records are generated. Site numbers and neighbourhood labels do not identify actual communities.

## Version 2 extension: `results/extended/`

- `synthetic_longitudinal.csv`: 2,400 rows, four per woman. IDs and places are fictional. `sample_day` and `outcome_day` are offsets from 2010-01-01; outcomes occur seven days after samples. `sample_gestational_week=46` means six weeks postpartum following fixed week-40 delivery, not an ongoing 46-week pregnancy. `temperature_C` is preceding seven-day Tmax; `heat` is standardised exposure. `module` is an expression surrogate. `vascular_z`, `depression_z`, `anxiety_z` are continuous standardised traits, not diagnoses. Remaining columns form X.
- `synthetic_clinical.csv`: later complication, ascertainment day, allele dosage and **generator-only** vascular intercept truth. This truth is not used as an observed predictor; the participant effect is inferred.
- `analysis_inputs.npz`: X (600,4,15), W (X plus module), y (600,4,3), mediator (600,4), Z (600,13) and event (600). Z uses prenatal averages only.
- `generating_truth.json`: seeds, true coefficients/covariances, temperature scaling and event count.
- `{joint,wide_prior,noisy_module}_population_draws.npz`: chain x draw x parameter arrays. B includes X then module; c includes the first 12 X columns then prenatal module. a models the mediator. D, sd, sm and loading describe covariances and clinical link.
- `*_participant_draws.npz`: b (4,2000,600,3), saved as float32; inference uses float64. Domain order: vascular, depression, anxiety.
- `*_diagnostics.csv`: rank R-hat, bulk ESS, tail ESS and mean MCSE for 1,887 targets, including all 1,800 subject effects. Summary JSONs record thresholds/extrema.
- `*_causal_effects.csv`: posterior means/medians, 95% equal-tail intervals, probability positive and generating truth. Continuous effects use standardised units; complications use risk differences. `population` averages empirical dosage; other rows use G=0,1,2. Truth remains the base generator in the noisy-module fit.
- `*_mediation_sensitivity.csv`: indirect effects at fixed rho; these are posterior transformations, not fitted rho. Noisy-module sensitivities describe the naive fitted model, not a corrected latent mechanism.
- `prior_sensitivity_comparison.csv`, `assay_error_comparison.csv`: alternative versus base estimates. `noisy_module_inputs.npz` saves the corrupted module; `assay_error_design.json` gives SD/reliability.
- `posterior_predictive_checks.csv`: 400 new-subject replications at observed covariates, with intervals and one-sided tail proportions (not classical p values). `prior_predictive_checks.csv`: 400 prior replications.
- `transcriptomic_scan.csv`: 500 synthetic features, true/estimated heat beta, clustered SE, p and BH q. `synthetic_expression_scan.npz`: expression (600,4,500), cell proxy and batch, separate from the prespecified mediator.
- `validation.json`, `run_summary.json`: numerical checks and main-run metadata. Assay-error design/diagnostics are separate.

## Legacy files directly under `results/`

The following sections describe a different synthetic cohort retained from version 1; these are not joint-model observations.

## `synthetic_daily_weather.csv.gz`

`date`: fictional daily calendar; `site`: 0-2; `neighbourhood`: 0-5; `tmax_C`: simulated daily maximum air temperature; `year`, `month`: calendar fields. Shared AR(1) weather induces dependence. No latitude or longitude is used.

## `synthetic_cohort.csv`

`synthetic_id`: unique generated ID; `conception`: fictional date; `day_index`: days since 2010-01-01; `heat_C`: mean Tmax on gestational days 98-139 (14 to <20 weeks); `heat_z`: that exposure standardised against a synthetic 18,000-record source pool; `heat_observed`: heat_z plus independent normal error; `season_sin`, `season_cos`: conception-season harmonics; `pc`: observed synthetic ancestry proxy; `G`: one fictitious diploid allele dosage; `Gc`: G-0.6; `age_z`: standardised age surrogate; `deprivation`: standard-normal social covariate with site-specific mean; `risk`: generating outcome probability; `preterm`: Bernoulli outcome; `retained`: independent follow-up indicator. Rows lost to follow-up retain latent outcomes in the synthetic file, but clinical analysis uses retained rows only.

## Molecular outputs

`synthetic_molecular_visits.csv`: participant ID, week (14/20/26), seven-day presampling heat_z, cell proxy and binary batch. `synthetic_expression.npz`: expression array (300 subjects, 3 visits, 500 features), feature labels. Subject order follows first appearance in the visit CSV. Features 1-25 have an injected heat coefficient 0.30; remaining features are null. These are not real genes. No molecular missingness is generated.

`molecular_estimates.csv`: feature, generating beta, estimated heat beta, cluster-robust standard error, unadjusted p and BH q. Visit-level genotype comes from the cohort ID. The molecular exposure z score has its own scale, distinct from the clinical gestational heat_z.

## Model and scenario outputs

`clinical_estimates.csv`: heat OR and interaction ratio of ORs, intervals and p values under oracle and assigned exposures. `trend_estimates.csv`: C/year and HC3 intervals by synthetic site. `design_operating_characteristics.csv`: recruited n, assumed effect ORs, error SD, successful/failed fits, rejection rates at .05 and 5e-8, Monte Carlo SE, log-coefficient bias, interaction coverage and mean retained n. Rejection under a null is a type-I-error diagnostic; under an alternative it is conditional power. Bias/coverage use the true generating coefficient even when noisy exposure changes the fitted estimand.

`summary.json`, `software_versions.json`, `validation.json`: realised counts, runtime versions and independent numerical check outcomes.
