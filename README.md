# Heat exposure, maternal and mental health: joint Bayesian simulation

**SYNTHETIC DESIGN STUDY. No real participants, climate observations, genotypes or transcriptomes.**

Prepared for Vusi, 23 October 2026. Version 2.0 keeps Dr Glory Chidumwa's spatiotemporal heat-exposure focus central. A new cohort links synthetic weather, repeated maternal and mental-health measurements, a prespecified expression mediator, and a later pregnancy complication. AWI-Gen participant records are not used.

## Report and main outputs

- [PDF report](output/pdf/heat_pregnancy_simulation_report.pdf)
- [Editable Markdown report](report.md)
- [Bayesian direct, indirect and total effects](results/extended/joint_causal_effects.csv)
- [Unmeasured-confounding sensitivity](results/extended/joint_mediation_sensitivity.csv)
- [Assay-error sensitivity](results/extended/assay_error_comparison.csv)
- [Prior sensitivity](results/extended/prior_sensitivity_comparison.csv)
- [Model diagnostics](results/extended/joint_diagnostics.csv)
- [Numerical validation](results/extended/validation.json)

## What is actually fitted

- **600 synthetic women, four visits each:** three antenatal and one postpartum.
- Gaussian vascular, depression and anxiety trajectories with correlated participant random intercepts.
- A later binary pregnancy complication fitted jointly using a probit likelihood sharing the vascular random intercept.
- A temporally preceding, prespecified synthetic expression module; heat-by-genotype terms in the module and outcome models.
- Four Gibbs chains per model, 1,500 warmup and 2,000 retained draws per chain.
- Three fits: base joint model, wider regression priors, and naive refit after transcriptomic assay error.
- Residual-confounding sensitivity over rho = -0.6 to 0.6; these are posterior reparameterisations at fixed rho, not estimates of rho.
- A separate 500-feature exploratory expression scan with clustered inference and BH adjustment.

Population and participant posterior draws are included. Participant draws are stored as float32 in separate compressed NPZ files; no inference was performed in float32. Files are relatively large because all participant effects are retained for audit and diagnostics.

## Reproduce

Use Python 3.12 in a virtual environment:

```bash
python -m venv .venv
# Activate with the command appropriate to your operating system.
python -m pip install -r requirements.txt
python simulate.py
python validate.py
python extend_simulation.py
python extended_sensitivity.py
python validate_extended.py
python build_report.py
```

The extension consumes the legacy synthetic weather file. Scripts resolve paths relative to themselves. Seeds and full coefficients are saved in `results/extended/generating_truth.json`. Numerical results may vary slightly across platforms.

For regenerating summaries from **unchanged** saved fits, `extend_simulation.py --reuse` and `extended_sensitivity.py --reuse` avoid resampling. Do not reuse after changing the data generator, sampler or prior. Full reproduction should use the commands above.

## Diagnostics and checks

ArviZ Stats computes rank-normalised folded split R-hat, bulk/tail ESS and mean MCSE for 1,887 scalar targets per fit, including all 1,800 participant effects. Each fit passed R-hat < 1.01 and bulk/tail ESS > 400. Trace plots, 400 new-subject posterior replications and 400 prior-predictive replications are included. HMC divergences, treedepth and BFMI do not apply to this Gibbs sampler.

Validation checks analytic conjugate posterior moments, Gaussian precision sampling, marginal probit integration, direct/indirect decomposition, rho=0 sensitivity equivalence and temporal ordering. This is not repeated Bayesian power/coverage calibration.

## Interpretation and limitations

The effects are injected assumptions. Vascular and mental-health traits are standardised surrogates, not clinical diagnoses. The synthetic complication endpoint is a liability model. Delivery is fixed at term; no delivery-dependent observation process, competing risks, dropout or missing assays are simulated. The complete molecular design entails 2,400 measurements and is not a budget or sample-size recommendation.

Expression is normalised continuous data, not RNA-seq counts. The mediator was prespecified and generated separately from the exploratory gene scan; no real gene or pathway has been discovered. Causal mediation interpretation requires assumptions satisfied by construction in the base case, not established by the fitted model. Sensitivity shows dependence on confounding and assay quality. Genomic SEM is discussed but not fitted.

The original 1,500-person preterm-birth simulation and 7,200 logistic fits remain under `results/`. They are separate from the extension and cannot be used as power calculations for this Bayesian model.

## Files

- `extend_simulation.py`: new cohort, joint Gibbs sampler, effects, diagnostics, predictive checks and expression scan.
- `extended_sensitivity.py`: module measurement-error refit.
- `validate_extended.py`: independent numerical and consistency checks.
- `build_report.py`, `build_extended_report.py`: current report builder.
- `simulate.py`, `validate.py`: retained legacy design simulation.
- `DATA_DICTIONARY.md`: definitions for both datasets.
- `MODEL_SPECIFICATION.md`: equations, priors and sensitivity estimands.

The user-supplied fellowship launch PDF informed context and is not redistributed. No institutional endorsement is implied. This private repository contains synthetic data only and is not a storage plan for identifiable location or genomic records.
