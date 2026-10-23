# Synthetic data dictionary

All records are generated. Site numbers and neighbourhood labels do not identify actual communities.

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
