# Model specification: executed extension

## Data and order

600 women, visits at nominal pregnancy weeks 14, 20, 26 and six weeks postpartum (fixed term delivery at day 280). A seven-day heat window ends at each module sample; repeated outcomes occur seven days later. A binary complication is ascertained by delivery, after the final antenatal outcome. No early delivery or loss is generated. Clinical predictors use only antenatal measurements. All data are synthetic.

Let X contain intercept, heat H, centred genotype G, H*G, age, deprivation, ancestry proxy, baseline distress, two site indicators, seasonal sine/cosine and three visit indicators. Complete coefficients are in `results/extended/generating_truth.json`.

## Likelihood

```
M_it = X_it a + eM_it,                  eM_it ~ N(0, sm^2)
Y_kit = X_it B_k + betaM_k M_it + u_ki + eY_kit
eY_kit ~ N(0, sigma_k^2),               u_i ~ MVN(0, D)
C_i = I(z_i > 0)
z_i = Z_i c + loading u_vascular,i + eC_i, eC_i ~ N(0,1)
```

The three repeated outcomes are vascular, depression and anxiety. Z includes prenatal-mean heat/module/season predictors and baseline covariates, with redundant averaged visit indicators omitted. The base generator makes eM independent of u and outcome errors. M has no extra random intercept. Conditional independence is a simplifying assumption, not a finding.

The complication observation updates the vascular random intercept in the Gibbs sampler. The correlated intercept covariance links the repeated outcomes. The endpoint and trajectories are fitted in a shared likelihood; this does not imply that mental health causes a maternal outcome.

## Priors and sampling

- Outcome and probit slopes: N(0,1); intercepts: N(0,4 variance).
- Continuous outcome variances: inverse-gamma(shape=3, scale=2).
- Random-intercept covariance: inverse-Wishart(df=5, scale=0.5 I).
- Mediator variance: inverse-gamma(3,2); a given sm^2: MVN(0, sm^2 I).
- Wider prior: regression SDs doubled, mediator conditional covariance multiplied by 4; residual and covariance priors unchanged.

Blocked Gibbs steps sample outcome regressions, latent truncated probit z, probit coefficients/loading, participant effects, covariance D, residual variances and the conjugate mediator block. All model updates use float64. Subject draws are saved in float32. The mediator posterior factors from outcome parameter inference because M is observed and the base likelihood errors/priors are independent.

Four chains, 1,500 warmup + 2,000 retained iterations each. Exact-conditional Gibbs has no NUTS/HMC tuning or divergence diagnostics. ArviZ rank R-hat, bulk ESS, 5%/95% tail ESS and mean MCSE are computed for population and subject parameters. Targets: R-hat<1.01, bulk/tail ESS>400.

## Mediation

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

## Residual-confounding sensitivity

Each outcome is considered separately. Rho is assumed correlation between structural eM and eY conditional on measured covariates and participant effects; it is not estimated. Grid -0.6,...,0.6. For a Gaussian outcome:

```
kappa = rho * sigma_conditional / (sm sqrt(1-rho^2))
betaM_structural = betaM_conditional - kappa
B_structural = B_conditional + kappa a
```

These transformations leave the observed conditional model and total effect unchanged while reallocating pathways. Distributions are induced by the base posterior, not refitted under independent structural priors. For the probit endpoint replace sm by sm/sqrt(3), and sigma_conditional by 1. The direct contrast gains kappa*deltaM; the indirect contrast loses it. Structural residual variance is 1/(1-rho^2), and its covariance with M must be retained. Algebra gives the same total variance as the base conditional model; this identity is validated. Sensitivities do not define a simultaneous free correlation matrix across outcomes.

## Other sensitivity analyses

1. Wider regression priors on unchanged data.
2. Independent N(0,0.5^2) assay noise added to M. Same outcomes/weather, naive refit without error correction. Conditional reliability=0.49/0.74=0.662.

## Boundaries

No empirical effect, real gene, mechanism, diagnostic performance or national exposure trend is established. No repeated Bayesian power/coverage study, efficiency comparison, survival analysis, RNA-seq count model or Genomic SEM is executed. The separate legacy preterm-birth simulations do not provide power for this joint model.
