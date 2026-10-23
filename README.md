# Heat exposure and maternal and mental health

## Research overview

A computational study of joint Bayesian modelling for heat exposure, repeated maternal and mental-health traits, an expression mediator and a later pregnancy complication. The analysis evaluates estimation under a known data-generating mechanism and examines sensitivity to residual confounding, regression priors and molecular assay error.

**All participants, weather observations, genotypes and expression measurements are synthetic.** No empirical heat–health association, clinical diagnosis or biological discovery is established. The saved analysis recovers positive indirect effects under the specified mechanism, while sensitivity analyses demonstrate dependence on the assumptions governing mediation and measurement quality.

[Main results](#main-results) · [Reproduction](#reproduction) · [Model specification](#model-specification) · [Data dictionary](#data-dictionary) · [Full PDF report](output/pdf/heat_pregnancy_simulation_report.pdf)

## Research question

The study examines how a joint model can represent heat exposure, maternal physical-health traits, depression and anxiety trajectories, genetic modification and expression-mediated pathways. It separates three questions: associations with physical health, associations with mental-health trajectories, and molecular susceptibility. The single fictional genotype modifies exposure effects; it is not a mediator of heat exposure.

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

### Exposure and measurement sequence

The synthetic daily climate field spans 2010–2030, with three fictional sites and six neighbourhoods per site. It imposes a 5°C seasonal amplitude and a 0.04°C/year trend, with correlated weather noise. Site baselines are 24, 27 and 29°C. Synthetic conception dates span 2027–2028 and are not forecasts. Visit-level heat has a mean of 27.64°C and SD of 4.36°C; the latter defines the one-SD contrast in this simulation.

Baseline covariates precede seven-day exposure windows, molecular sampling and repeated outcomes one week later. Delivery is fixed at day 280; the postpartum sample occurs at day 322. Only antenatal exposures and modules enter the complication model. This sequence preserves temporal ordering but does not establish causal identification in an observational cohort.

### Generating assumptions

| Component | Assumption |
| --- | --- |
| Cohort | 600 women; balanced fictional sites; 4 visits each. Conception dates span 2027-2028. These are synthetic future dates, not forecasts. |
| Genotype | One fictional diploid allele; logit allele frequency = logit(0.30) + 0.20*ancestry. Dosage centred at 0.6. |
| Mediator | Normal expression module with heat coefficient 0.40, heat x genotype 0.12 and residual SD 0.70. Covariates and visit indicators also enter. |
| Vascular / depression / anxiety | Direct heat coefficients: 0.18 / 0.15 / 0.12. Module coefficients: 0.30 / 0.40 / 0.35. Direct heat x genotype: 0.08 / 0.06 / 0.05. Residual SD: 0.65 / 0.80 / 0.75. |
| Shared variation | Correlated participant random intercepts. Covariance diagonal 0.49, 0.64, 0.49; off-diagonal 0.21, 0.14, 0.32. |
| Complication | Probit intercept -1.35; heat 0.20; module 0.35; heat x genotype 0.10; shared vascular-intercept loading 0.65. Probit residual variance fixed at 1. |
| Measured confounders | Age, deprivation, ancestry proxy, baseline distress, site and seasonal harmonics; visit indicators for repeated measures. Complete coefficient vectors saved in generating_truth.json. |

The same heat series drives module expression and both outcome domains. The module precedes the outcomes. At generation, mediator errors are independent of participant effects and outcome errors: the base case has no residual mediator-outcome confounding. Outcome errors are conditionally independent across time and domains; correlated participant intercepts create dependence.

A separate 500-feature normalised-expression matrix is generated for the exploratory scan. Its features are not used to select the mediator. This avoids outcome-guided mediator selection, but also means the simulation does not validate how a real expression module would be discovered or measured.

## Statistical methods

The model combines Gaussian repeated-outcome likelihoods, a Gaussian mediator model and a probit complication likelihood. Covariates include heat, centred genotype, their interaction, age, deprivation, an ancestry proxy, baseline distress, site, season and visit. Complication predictors use antenatal averages only. Posterior sampling uses blocked Gibbs updates with exact conditional distributions.

Effects represent a one-standard-deviation increase in heat. Continuous-outcome direct and indirect effects are derived from model coefficients and averaged over observed genotype frequencies. Complication effects are marginal risk differences obtained by integrating over new-participant random effects and prenatal mediator residuals. Direct and indirect contrasts sum to the total contrast under the specified structural model.

Residual-confounding sensitivity fixes rho between −0.6 and 0.6, transforming the base posterior separately for each outcome. Rho is not estimated. The expression scan uses participant demeaning, covariate adjustment, participant-clustered inference and Benjamini–Hochberg adjustment across 500 heat tests. Discovered features are not fed back into mediator selection.

## Main results

### Exposure and outcome trajectories

![Synthetic heat exposure by site and repeated outcome trajectories](results/extended/05_joint_trajectories.png)

**Figure 1. Exposure and outcome trajectories.** The left panel shows mean seven-day maximum temperature by fictional site. The right panel shows the three standardised outcomes over four visits, with visit 4 representing postpartum follow-up. Patterns reflect injected time effects, weather and participant heterogeneity, not observed pregnancy symptom trajectories. The cohort contains 86 generated complications among 600 participants (14.3%) and 2,400 complete molecular/outcome pairs.

### Genetic modification

| Response | Heat x genotype: mean (95% CrI) | Injected value |
| --- | --- | --- |
| vascular | 0.128 (0.080, 0.175) | 0.08 |
| depression | 0.033 (-0.026, 0.092) | 0.06 |
| anxiety | 0.066 (0.009, 0.122) | 0.05 |
| Expression module | 0.136 (0.093, 0.178) | 0.12 |

These are coefficients per extra fictional allele, not real loci. A lack of interval separation from zero is retained rather than interpreted as proof of no modification. Allele-specific indirect effects are available in the results CSV; no genome-wide discovery analysis was attempted.

### Direct, indirect and total effects

| Outcome | Direct effect (95% CrI) | Indirect effect (95% CrI) | Total effect (95% CrI) |
|---|---:|---:|---:|
| Vascular | 0.145 (0.062, 0.225) | 0.123 (0.096, 0.151) | 0.268 (0.184, 0.350) |
| Depression | 0.191 (0.091, 0.290) | 0.177 (0.141, 0.217) | 0.368 (0.266, 0.470) |
| Anxiety | 0.157 (0.063, 0.250) | 0.142 (0.111, 0.176) | 0.299 (0.205, 0.393) |
| Complication | -0.124 (-7.830, 10.351) | 3.644 (0.831, 7.091) | 3.520 (-5.355, 15.046) |

Values are posterior means with equal-tail 95% credible intervals. Vascular, depression and anxiety effects use standardised outcome units. **Complication effects in this table are percentage points**; the source CSV and figure use unscaled risk differences. Source: [joint effect estimates](results/extended/joint_causal_effects.csv).

![Joint Bayesian effects of a one-SD heat shift](results/extended/06_joint_effects.png)

**Figure 2. Estimated pathways under the synthetic mechanism.** Points and intervals show posterior estimates and 95% credible intervals; crosses identify generating truth. The complication total-effect interval includes zero despite a positive indirect-effect interval. These are single-cohort recovery results, not estimates of repeated-simulation power or coverage.

### Genotype-specific mediation

| Fictional allele dosage | Depression indirect effect (95% CrI) | Truth |
| --- | --- | --- |
| 0 | 0.141 (0.105, 0.180) | 0.131 |
| 1 | 0.198 (0.159, 0.239) | 0.179 |
| 2 | 0.254 (0.205, 0.305) | 0.227 |

The indirect effect can differ by genotype because heat x genotype is included in the mediator equation. The mediator coefficient in the outcome equation is constant in this generator. This is a specific form of moderated mediation, not an exhaustive model of genetic susceptibility.

Antennatal modules precede the complication endpoint; postpartum modules contribute only to postpartum repeated outcomes. The base model assumes no exposure-induced mediator–outcome confounders, interference or treatment feedback. Correct timing alone would not eliminate confounding by illness, treatment, cell composition or social adversity in a real cohort.

### Sensitivity to residual confounding

![Indirect-effect sensitivity to assumed residual correlation](results/extended/07_mediation_sensitivity.png)

**Figure 3. Residual-confounding sensitivity.** Increasing positive rho reduces the inferred indirect effect in the recorded analysis. Intervals are conditional on each assumed rho and are induced by reparameterisation of the base posterior. They do not estimate the magnitude of unmeasured confounding. [Download sensitivity estimates](results/extended/joint_mediation_sensitivity.csv).

Wider regression priors produced similar continuous-outcome indirect-effect estimates. Adding independent noise with SD 0.5 to the module reduced its conditional reliability to approximately 0.662 and attenuated the indirect-effect estimates in the naive refit. The latter analysis does not correct measurement error. Detailed comparisons are available for [regression priors](results/extended/prior_sensitivity_comparison.csv) and [assay error](results/extended/assay_error_comparison.csv).

### Regression-prior and assay-error comparisons

| Outcome | Base mean | Wider mean | Wider 95% CrI |
| --- | --- | --- | --- |
| vascular | 0.123 | 0.123 | (0.096, 0.152) |
| depression | 0.177 | 0.178 | (0.142, 0.217) |
| anxiety | 0.142 | 0.142 | (0.111, 0.176) |
| complication | 3.644 | 3.516 | (0.868, 6.918) |

The wider-prior fit doubles the regression prior SDs. The inverse-Wishart covariance and residual-variance priors are unchanged, so this is a targeted regression-prior check rather than a complete prior-robustness analysis.

#### Module measurement-error stress test

| Outcome | Base mean | Noisy mean | Noisy 95% CrI |
| --- | --- | --- | --- |
| vascular | 0.123 | 0.074 | (0.053, 0.097) |
| depression | 0.177 | 0.103 | (0.076, 0.134) |
| anxiety | 0.142 | 0.093 | (0.068, 0.121) |
| complication | 3.644 | 2.760 | (0.598, 5.529) |

Independent SD 0.5 noise was added to the observed module while leaving the same weather, outcomes and participants unchanged. Conditional module reliability is 0.49/(0.49+0.25)=0.662. The original model was refitted without correcting the error. This evaluates the consequences of an imperfect assay; it is not a measurement-error-adjusted estimator.

The tables use percentage points for complications and standardised units otherwise. These results motivate assay validation and confounding assessment before interpreting a single indirect-effect estimate.

### Exploratory expression scan

Of 25 features with injected heat effects, 22 were detected at BH q < 0.05. One of 475 null features was also selected: realised sensitivity was 88.0% and the realised false-discovery proportion was 4.3%. These describe this generated dataset, not repeated-simulation FDR control. The features are normalised continuous surrogates without real gene identifiers. [Download the complete scan](results/extended/transcriptomic_scan.csv).

### Computational diagnostics

| Model | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS |
|---|---:|---:|---:|
| Base joint model | 1.0065 | 712.0 | 1182.1 |
| Wider regression priors | 1.0070 | 640.8 | 1292.5 |
| Noisy expression module | 1.0062 | 679.6 | 1007.2 |

Each model assessed 1,887 scalar targets, including 1,800 participant effects. All recorded fits met R-hat < 1.01 and bulk/tail ESS > 400. These diagnostics assess sampling behaviour rather than model validity. HMC divergences, treedepth and BFMI do not apply to the Gibbs sampler. [Recorded validation](results/extended/validation.json) also checks conjugate moments, precision sampling, probit integration, effect decomposition, rho=0 equivalence and temporal ordering.

![Four-chain post-warmup Gibbs traces for selected parameters](results/extended/08_traceplots.png)

**Figure 4. Post-warmup sampling traces.** Four chains are displayed for the heat-to-module coefficient, module-to-depression coefficient and vascular loading in the complication model. Each chain contributes 2,000 retained draws. These selected traces illustrate mixing; numerical diagnostics also cover every population parameter and participant effect. [Full base-model diagnostics](results/extended/joint_diagnostics.csv).

### Predictive checks and numerical validation

![New-subject posterior predictive intervals and observed synthetic summaries](results/extended/09_predictive_checks.png)

**Figure 5. Posterior predictive assessment.** Horizontal intervals describe new-subject posterior replications; crosses mark observed synthetic summaries. The plot displays selected statistics; the complete table includes 15 prespecified summaries. All 15 observed values fall within their central 95% posterior predictive intervals. Replication generates new participant effects, modules, repeated outcomes and complications at the observed covariates.

The analysis uses 400 posterior and 400 prior predictive replications. Prior predictions include extreme complication rates and have not been clinically elicited. Predictive agreement within this generated cohort does not establish external calibration. Download the [posterior checks](results/extended/posterior_predictive_checks.csv) and [prior checks](results/extended/prior_predictive_checks.csv).

The recorded numerical validation passed checks of analytic mediator posterior moments, Gaussian precision sampling, marginal probit integration, direct/indirect decomposition, rho=0 sensitivity equivalence and temporal ordering. These checks assess implementation and internal consistency, not frequentist coverage or clinical prediction.

## Interpretation and limitations

This study demonstrates an executable joint model under a transparent synthetic mechanism. It does not establish a causal effect of environmental heat or expression on maternal health. Mediation interpretation depends on structural assumptions, including the base model's independent mediator and outcome errors. The sensitivity analyses show how indirect-effect conclusions change when these assumptions or assay quality change.

The simulation omits early delivery, fetal loss, competing risks, dropout, missing assays, personal-exposure error and time-varying treatment/confounding. Gaussian vascular and mental-health surrogates are not validated diagnostic instruments. Expression is not simulated as RNA-seq counts, and the scan lacks realistic gene–gene correlation and pathway annotation. No GWAS, polygenic-score validation or Genomic SEM is fitted.

Only one generated cohort supports the Bayesian recovery results. Convergence checks and alternative fits do not establish power, interval coverage or an efficiency advantage over separate outcome models. The independent earlier simulation retained under `results/` comprises 1,500 recruits and 7,200 logistic fits across 24 scenarios; its operating characteristics do not apply to this joint Bayesian model.

## Implications for subsequent study design

A subsequent empirical study would first characterise spatial and seasonal exposure variation, validate sensors and sampling windows, and specify primary clinical and mental-health endpoints. Molecular mediators should be prespecified or independently validated, with exploratory feature discovery kept separate from confirmatory mediation. Assay reliability, missingness and realistic recruitment costs require explicit design work.

The proposed development sequence is exposure and endpoint definition, environmental and specimen-handling pilots, repeated joint-model simulations across cohort and molecular-substudy sizes, and a costed protocol with ethics and data-governance arrangements. Proposed contributions from HE2AT, Bio-HEAT and AWI-Gen/SBIMB investigators remain subject to agreement; no collaboration or institutional endorsement is established by this simulation.

Further simulations should address early delivery, outcome-dependent observation, correlated count-based expression, exposure error and time-varying confounding, and compare joint with separate outcome models. Genomic SEM would require suitable external GWAS summary statistics and ancestry-compatible LD information; it is not evaluated by the present single-variant experiment.

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
| [build_report.py](build_report.py), [build_extended_report.py](build_extended_report.py) | Generate the PDF report from saved results; this README is maintained separately |
| [Full PDF report](output/pdf/heat_pregnancy_simulation_report.pdf) | Detailed results, contextual discussion and references |
| [Software record](results/extended/software_versions.json) | Recorded analysis environment |

This README consolidates the former Markdown report, methods and data dictionary. The PDF remains available as a standalone narrative. Suggested institutional contributions are proposals rather than confirmed collaborations or endorsement. The contextual fellowship document is not redistributed, and no AWI-Gen participant records are used.

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

## References and evidence provenance

The following contextual and methodological sources are retained from the full report. No published effect was used to calibrate the generator. The custom parametric model does not reproduce the semiparametric model cited below; no ERA5-Land observations were analysed.

[1] DELTAS Africa II Career Advancement Fellowship. Draft research programme launch: Modelling of Heat Exposure and Maternal and Mental Health in South Africa. Dr Glory Chidumwa, 22 October 2026. User-supplied PDF, contextual source only; not a funding call or instruction.

[2] Chersich MF et al. Associations between high temperatures in pregnancy and risk of preterm birth, low birth weight, and stillbirths. BMJ 2020;371:m3811. https://doi.org/10.1136/bmj.m3811

[3] A systematic review and meta-analysis of heat exposure impacts on maternal, fetal and neonatal health. Nature Medicine. https://doi.org/10.1038/s41591-024-03395-8

[4] Study protocol for the Bio-HEAT study. Wellcome Open Research, 2025. https://pmc.ncbi.nlm.nih.gov/articles/PMC12209768/

[5] DS-I Africa. HE2AT Center programme description. https://www.dsi-africa.org/project/5

[6] Joint modelling of mental health markers through pregnancy: a Bayesian semi-parametric approach. https://doi.org/10.1080/02664763.2022.2154329

[7] Grotzinger AD et al. Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits. Nature Human Behaviour 2019. https://doi.org/10.1038/s41562-019-0566-x

[8] Vehtari A et al. Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. Bayesian Analysis 2021;16:667-718. https://doi.org/10.1214/20-BA1221

[9] ArviZ Stats. Rank-normalised R-hat and array diagnostics documentation. https://python.arviz.org/projects/stats/en/latest/api/generated/arviz_stats.rhat.html

[10] Tingley D et al. mediation: R Package for Causal Mediation Analysis. Journal of Statistical Software 2014;59(5). https://doi.org/10.18637/jss.v059.i05

[11] Bind MA et al. Causal mediation analysis for longitudinal data with exogenous exposure. Biostatistics 2016;17:122-134. https://doi.org/10.1093/biostatistics/kxv029

[12] ECMWF. ERA5-Land dataset description. https://www.ecmwf.int/en/era5-land

[13] AWI-Gen cohort profile. https://pmc.ncbi.nlm.nih.gov/articles/PMC11790221/
