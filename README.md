# Heat exposure and maternal and mental health

## Research overview

A computational study of joint Bayesian modelling for heat exposure, repeated maternal and mental-health traits, an expression mediator and a later pregnancy complication. The analysis evaluates estimation under a known data-generating mechanism and examines sensitivity to residual confounding, regression priors and molecular assay error.

**All participants, weather observations, genotypes and expression measurements are synthetic.** No empirical heat–health association, clinical diagnosis or biological discovery is established. The saved analysis recovers positive indirect effects under the specified mechanism, while sensitivity analyses demonstrate dependence on the assumptions governing mediation and measurement quality.

[Main results](#main-results) · [Reproduction](#reproduction) · [Model specification](#model-specification) · [Data dictionary](#data-dictionary) · [Full PDF report](output/pdf/heat_pregnancy_simulation_report.pdf)

## Study design

| Component | Specification |
|---|---|
| Cohort | 600 synthetic women, with four complete visits each |
| Visit schedule | Pregnancy weeks 14, 20 and 26, followed by six weeks postpartum |
| Exposure | Seven-day heat window ending at each molecular sample |
| Temporal order | Repeated outcomes measured seven days after the module sample; complication assessed by delivery |
| Repeated outcomes | Standardised vascular, depression and anxiety traits |
| Molecular mediator | Prespecified continuous expression module, generated separately from the exploratory scan |
| Genetic component | One synthetic diploid variant, with heat-by-genotype terms |
| Complication endpoint | Binary probit outcome; 86 events in the recorded cohort |
| Sampling | Four Gibbs chains, each with 1,500 warmup and 2,000 retained draws |
| Sensitivity fits | Wider regression priors and a naive refit after expression-module measurement error |
| Exploratory expression analysis | 500 synthetic features; 25 injected heat effects |

The three repeated outcomes have correlated participant random intercepts. The complication likelihood shares the vascular intercept, allowing the endpoint to inform the longitudinal latent effect. Fixed term delivery and complete observation are design assumptions. The 2,400 molecular measurements represent the simulated design, not a funded cohort or sample-size recommendation. Full generating coefficients and seed 2026102302 are recorded in [generating_truth.json](results/extended/generating_truth.json).

## Statistical methods

The model combines Gaussian repeated-outcome likelihoods, a Gaussian mediator model and a probit complication likelihood. Covariates include heat, centred genotype, their interaction, age, deprivation, an ancestry proxy, baseline distress, site, season and visit. Complication predictors use antenatal averages only. Posterior sampling uses blocked Gibbs updates with exact conditional distributions.

Effects represent a one-standard-deviation increase in heat. Continuous-outcome direct and indirect effects are derived from model coefficients and averaged over observed genotype frequencies. Complication effects are marginal risk differences obtained by integrating over new-participant random effects and prenatal mediator residuals. Direct and indirect contrasts sum to the total contrast under the specified structural model.

Residual-confounding sensitivity fixes rho between −0.6 and 0.6, transforming the base posterior separately for each outcome. Rho is not estimated. The expression scan uses participant demeaning, covariate adjustment, participant-clustered inference and Benjamini–Hochberg adjustment across 500 heat tests. Discovered features are not fed back into mediator selection.

## Main results

### Direct, indirect and total effects

| Outcome | Direct effect (95% CrI) | Indirect effect (95% CrI) | Total effect (95% CrI) |
|---|---:|---:|---:|
| Vascular | 0.145 (0.062, 0.225) | 0.123 (0.096, 0.151) | 0.268 (0.184, 0.350) |
| Depression | 0.191 (0.091, 0.290) | 0.177 (0.141, 0.217) | 0.368 (0.266, 0.470) |
| Anxiety | 0.157 (0.063, 0.250) | 0.142 (0.111, 0.176) | 0.299 (0.205, 0.393) |
| Complication | -0.124 (-7.830, 10.351) | 3.644 (0.831, 7.091) | 3.520 (-5.355, 15.046) |

Values are posterior means with equal-tail 95% credible intervals. Vascular, depression and anxiety effects use standardised outcome units. **Complication effects in this table are percentage points**; the source CSV and figure use unscaled risk differences. Source: [joint effect estimates](results/extended/joint_causal_effects.csv).

![Joint Bayesian effects of a one-SD heat shift](results/extended/06_joint_effects.png)

**Figure 1. Estimated pathways under the synthetic mechanism.** Points and intervals show posterior estimates and 95% credible intervals; crosses identify generating truth. The complication total-effect interval includes zero despite a positive indirect-effect interval. These are single-cohort recovery results, not estimates of repeated-simulation power or coverage.

### Sensitivity to residual confounding

![Indirect-effect sensitivity to assumed residual correlation](results/extended/07_mediation_sensitivity.png)

**Figure 2. Residual-confounding sensitivity.** Increasing positive rho reduces the inferred indirect effect in the recorded analysis. Intervals are conditional on each assumed rho and are induced by reparameterisation of the base posterior. They do not estimate the magnitude of unmeasured confounding. [Download sensitivity estimates](results/extended/joint_mediation_sensitivity.csv).

Wider regression priors produced similar continuous-outcome indirect-effect estimates. Adding independent noise with SD 0.5 to the module reduced its conditional reliability to approximately 0.662 and attenuated the indirect-effect estimates in the naive refit. The latter analysis does not correct measurement error. Detailed comparisons are available for [regression priors](results/extended/prior_sensitivity_comparison.csv) and [assay error](results/extended/assay_error_comparison.csv).

### Exploratory expression scan

Of 25 features with injected heat effects, 22 were detected at BH q < 0.05. One of 475 null features was also selected: realised sensitivity was 88.0% and the realised false-discovery proportion was 4.3%. These describe this generated dataset, not repeated-simulation FDR control. The features are normalised continuous surrogates without real gene identifiers. [Download the complete scan](results/extended/transcriptomic_scan.csv).

### Computational diagnostics

| Model | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS |
|---|---:|---:|---:|
| Base joint model | 1.0065 | 712.0 | 1182.1 |
| Wider regression priors | 1.0070 | 640.8 | 1292.5 |
| Noisy expression module | 1.0062 | 679.6 | 1007.2 |

Each model assessed 1,887 scalar targets, including 1,800 participant effects. All recorded fits met R-hat < 1.01 and bulk/tail ESS > 400. These diagnostics assess sampling behaviour rather than model validity. HMC divergences, treedepth and BFMI do not apply to the Gibbs sampler. [Recorded validation](results/extended/validation.json) also checks conjugate moments, precision sampling, probit integration, effect decomposition, rho=0 equivalence and temporal ordering.

Saved [trace plots](results/extended/08_traceplots.png) and [predictive-check plots](results/extended/09_predictive_checks.png) complement the numerical diagnostics. Predictive checks use 400 new-subject posterior replications and 400 prior replications.

## Interpretation and limitations

This study demonstrates an executable joint model under a transparent synthetic mechanism. It does not establish a causal effect of environmental heat or expression on maternal health. Mediation interpretation depends on structural assumptions, including the base model's independent mediator and outcome errors. The sensitivity analyses show how indirect-effect conclusions change when these assumptions or assay quality change.

The simulation omits early delivery, fetal loss, competing risks, dropout, missing assays, personal-exposure error and time-varying treatment/confounding. Gaussian vascular and mental-health surrogates are not validated diagnostic instruments. Expression is not simulated as RNA-seq counts, and the scan lacks realistic gene–gene correlation and pathway annotation. No GWAS, polygenic-score validation or Genomic SEM is fitted.

Only one generated cohort supports the Bayesian recovery results. Convergence checks and alternative fits do not establish power, interval coverage or an efficiency advantage over separate outcome models. The independent earlier simulation retained under `results/` comprises 1,500 recruits and 7,200 logistic fits across 24 scenarios; its operating characteristics do not apply to this joint Bayesian model.

## Reproduction

Use Python 3.12 with the pinned dependencies in [requirements.txt](requirements.txt):

```bash
python -m venv .venv
# Activate the environment using the command appropriate to your operating system.
python -m pip install -r requirements.txt
python simulate.py
python validate.py
python extend_simulation.py
python extended_sensitivity.py
python validate_extended.py
python build_report.py
```

The extension consumes the earlier synthetic weather file. Scripts resolve paths relative to their source directory. Full execution regenerates datasets, posterior fits, validation outputs and reports; preserve copies before changing assumptions. Numerical results may vary slightly across platforms. This documentation update uses the saved outputs and does not represent a new model run.

For summaries of **unchanged** saved fits, `python extend_simulation.py --reuse` and `python extended_sensitivity.py --reuse` avoid resampling. Reuse is inappropriate after modifying the generator, sampler or priors. Population and participant posterior draws remain available for audit; participant effects are stored as float32, while inference uses float64.

| Resource | Purpose |
|---|---|
| [extend_simulation.py](extend_simulation.py) | Cohort generation, joint sampler, effects, diagnostics, predictive checks and expression scan |
| [extended_sensitivity.py](extended_sensitivity.py) | Expression assay-error refit |
| [validate_extended.py](validate_extended.py) | Numerical and consistency checks |
| [simulate.py](simulate.py), [validate.py](validate.py) | Earlier climate, clinical and design simulations |
| [build_report.py](build_report.py), [build_extended_report.py](build_extended_report.py) | Generate the full report from saved results |
| [Full PDF report](output/pdf/heat_pregnancy_simulation_report.pdf), [report source](report.md) | Detailed results, contextual discussion and references |
| [Software record](results/extended/software_versions.json) | Recorded analysis environment |

The full report retains the wider project-development context. Suggested institutional contributions are proposals rather than confirmed collaborations or endorsement. The contextual fellowship document is not redistributed, and no AWI-Gen participant records are used.

## Model specification

<details>
<summary>Likelihood, priors, estimands and sensitivity transformations</summary>


### Data and order

600 women, visits at nominal pregnancy weeks 14, 20, 26 and six weeks postpartum (fixed term delivery at day 280). A seven-day heat window ends at each module sample; repeated outcomes occur seven days later. A binary complication is ascertained by delivery, after the final antenatal outcome. No early delivery or loss is generated. Clinical predictors use only antenatal measurements. All data are synthetic.

Let X contain intercept, heat H, centred genotype G, H*G, age, deprivation, ancestry proxy, baseline distress, two site indicators, seasonal sine/cosine and three visit indicators. Complete coefficients are in `results/extended/generating_truth.json`.

### Likelihood

```
M_it = X_it a + eM_it,                  eM_it ~ N(0, sm^2)
Y_kit = X_it B_k + betaM_k M_it + u_ki + eY_kit
eY_kit ~ N(0, sigma_k^2),               u_i ~ MVN(0, D)
C_i = I(z_i > 0)
z_i = Z_i c + loading u_vascular,i + eC_i, eC_i ~ N(0,1)
```

The three repeated outcomes are vascular, depression and anxiety. Z includes prenatal-mean heat/module/season predictors and baseline covariates, with redundant averaged visit indicators omitted. The base generator makes eM independent of u and outcome errors. M has no extra random intercept. Conditional independence is a simplifying assumption, not a finding.

The complication observation updates the vascular random intercept in the Gibbs sampler. The correlated intercept covariance links the repeated outcomes. The endpoint and trajectories are fitted in a shared likelihood; this does not imply that mental health causes a maternal outcome.

### Priors and sampling

- Outcome and probit slopes: N(0,1); intercepts: N(0,4 variance).
- Continuous outcome variances: inverse-gamma(shape=3, scale=2).
- Random-intercept covariance: inverse-Wishart(df=5, scale=0.5 I).
- Mediator variance: inverse-gamma(3,2); a given sm^2: MVN(0, sm^2 I).
- Wider prior: regression SDs doubled, mediator conditional covariance multiplied by 4; residual and covariance priors unchanged.

Blocked Gibbs steps sample outcome regressions, latent truncated probit z, probit coefficients/loading, participant effects, covariance D, residual variances and the conjugate mediator block. All model updates use float64. Subject draws are saved in float32. The mediator posterior factors from outcome parameter inference because M is observed and the base likelihood errors/priors are independent.

Four chains, 1,500 warmup + 2,000 retained iterations each. Exact-conditional Gibbs has no NUTS/HMC tuning or divergence diagnostics. ArviZ rank R-hat, bulk ESS, 5%/95% tail ESS and mean MCSE are computed for population and subject parameters. Targets: R-hat<1.01, bulk/tail ESS>400.

### Mediation

For a +1 SD heat shift, at centred genotype g:

```
a_g = a_heat + a_heatG g
NIE_k = a_g betaM_k
NDE_k = B_heat,k + B_heatG,k g
TE_k = NDE_k + NIE_k
```

Population effects average over observed genotype frequencies; G=0,1,2 are also reported. No exposure-mediator interaction occurs in outcome models.

Complication standardisation uses empirical baseline covariates and integrates new-subject u and mean prenatal mediator residual. The probability denominator is:

```
sqrt(1 + loading^2 D_vascular,vascular + c_M^2 sm^2/3)
```

Risks use baseline exposure/module (R00), shifted exposure with baseline module (R10), and shifted exposure/module (R11). NDE=R10-R00, NIE=R11-R10, TE=R11-R00. This is a natural-effect decomposition under the stated structural model, not an odds-ratio product. Participant effects are integrated at their unconditional distribution, valid here because they are independent of X and mediator errors.

### Residual-confounding sensitivity

Each outcome is considered separately. Rho is assumed correlation between structural eM and eY conditional on measured covariates and participant effects; it is not estimated. Grid -0.6,...,0.6. For a Gaussian outcome:

```
kappa = rho * sigma_conditional / (sm sqrt(1-rho^2))
betaM_structural = betaM_conditional - kappa
B_structural = B_conditional + kappa a
```

These transformations leave the observed conditional model and total effect unchanged while reallocating pathways. Distributions are induced by the base posterior, not refitted under independent structural priors. For the probit endpoint replace sm by sm/sqrt(3), and sigma_conditional by 1. The direct contrast gains kappa*deltaM; the indirect contrast loses it. Structural residual variance is 1/(1-rho^2), and its covariance with M must be retained. Algebra gives the same total variance as the base conditional model; this identity is validated. Sensitivities do not define a simultaneous free correlation matrix across outcomes.

### Other sensitivity analyses

1. Wider regression priors on unchanged data.
2. Independent N(0,0.5^2) assay noise added to M. Same outcomes/weather, naive refit without error correction. Conditional reliability=0.49/0.74=0.662.

### Boundaries

No empirical effect, real gene, mechanism, diagnostic performance or national exposure trend is established. No repeated Bayesian power/coverage study, efficiency comparison, survival analysis, RNA-seq count model or Genomic SEM is executed. The separate legacy preterm-birth simulations do not provide power for this joint model.

</details>

## Data dictionary

<details>
<summary>Dataset fields, posterior arrays and output definitions</summary>


All records are generated. Site numbers and neighbourhood labels do not identify actual communities.

### Version 2 extension: `results/extended/`

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

### Legacy files directly under `results/`

The following sections describe a different synthetic cohort retained from version 1; these are not joint-model observations.

### `synthetic_daily_weather.csv.gz`

`date`: fictional daily calendar; `site`: 0-2; `neighbourhood`: 0-5; `tmax_C`: simulated daily maximum air temperature; `year`, `month`: calendar fields. Shared AR(1) weather induces dependence. No latitude or longitude is used.

### `synthetic_cohort.csv`

`synthetic_id`: unique generated ID; `conception`: fictional date; `day_index`: days since 2010-01-01; `heat_C`: mean Tmax on gestational days 98-139 (14 to <20 weeks); `heat_z`: that exposure standardised against a synthetic 18,000-record source pool; `heat_observed`: heat_z plus independent normal error; `season_sin`, `season_cos`: conception-season harmonics; `pc`: observed synthetic ancestry proxy; `G`: one fictitious diploid allele dosage; `Gc`: G-0.6; `age_z`: standardised age surrogate; `deprivation`: standard-normal social covariate with site-specific mean; `risk`: generating outcome probability; `preterm`: Bernoulli outcome; `retained`: independent follow-up indicator. Rows lost to follow-up retain latent outcomes in the synthetic file, but clinical analysis uses retained rows only.

### Molecular outputs

`synthetic_molecular_visits.csv`: participant ID, week (14/20/26), seven-day presampling heat_z, cell proxy and binary batch. `synthetic_expression.npz`: expression array (300 subjects, 3 visits, 500 features), feature labels. Subject order follows first appearance in the visit CSV. Features 1-25 have an injected heat coefficient 0.30; remaining features are null. These are not real genes. No molecular missingness is generated.

`molecular_estimates.csv`: feature, generating beta, estimated heat beta, cluster-robust standard error, unadjusted p and BH q. Visit-level genotype comes from the cohort ID. The molecular exposure z score has its own scale, distinct from the clinical gestational heat_z.

### Model and scenario outputs

`clinical_estimates.csv`: heat OR and interaction ratio of ORs, intervals and p values under oracle and assigned exposures. `trend_estimates.csv`: C/year and HC3 intervals by synthetic site. `design_operating_characteristics.csv`: recruited n, assumed effect ORs, error SD, successful/failed fits, rejection rates at .05 and 5e-8, Monte Carlo SE, log-coefficient bias, interaction coverage and mean retained n. Rejection under a null is a type-I-error diagnostic; under an alternative it is conditional power. Bias/coverage use the true generating coefficient even when noisy exposure changes the fitted estimand.

`summary.json`, `software_versions.json`, `validation.json`: realised counts, runtime versions and independent numerical check outcomes.

</details>
