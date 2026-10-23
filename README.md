# Heat exposure, pregnancy and molecular susceptibility

**SYNTHETIC DESIGN STUDY. No real participants, weather observations, genotypes or transcriptomes.**

Prepared for Vusi, 23 October 2026. Incorporates Dr Glory Chidumwa's proposed focus on spatiotemporal heat-exposure trends into a purpose-built prospective pregnancy study. The existing AWI-Gen cohort is not the proposed recruitment population.

Version 1.1 includes a methodological appendix on the expanded maternal physical and mental health question, Bayesian joint modelling, genetic modification, transcriptomic mediation and Genomic SEM. **These extensions are proposed; the completed simulation does not generate mental-health outcomes or fit Bayesian/Genomic SEM models.**

## Read the report

- [PDF report](output/pdf/heat_pregnancy_simulation_report.pdf)
- [Editable Markdown report](report.md)
- [Scenario results](results/design_operating_characteristics.csv)
- [Validation](results/validation.json)

## Reproduce

Use Python 3.12 in a virtual environment:

```bash
python -m venv .venv
# Activate the environment using your operating system's standard command.
python -m pip install -r requirements.txt
python simulate.py
python validate.py
python build_report.py
```

The simulation runs 7,200 clinical fits (24 scenarios, 300 repetitions each), plus a demonstration cohort and an expression substudy. Run from any directory; outputs resolve relative to the scripts. Runtime depends on hardware. Seeds are fixed; exact floating-point results may vary slightly across numerical libraries. The generated configuration is a descriptive summary; edit `simulate.py` to alter generating parameters.

## Contents

- `simulate.py`: weather, cohort, expression and operating-characteristic simulation.
- `validate.py`: independent optimiser comparison and output consistency checks.
- `build_report.py`: report builder using saved numerical outputs.
- `simulation_config.json`: key generating assumptions, not an input file.
- `results/`: synthetic data, estimates, figures, software versions and validation.
- `output/pdf/`: report.
- `DATA_DICTIONARY.md`: data definitions and interpretive boundaries.

## Interpretation

The climate trend and biological effects are injected assumptions, not findings about South Africa. No real susceptibility genes are tested. The molecular matrix represents normalised continuous expression, not RNA-seq counts. The clinical model is a simplified 20-week landmark binary outcome and does not simulate delivery timing, stillbirth, miscarriage or missing molecular visits. Nominal single-variant power is not genome-wide discovery power. The analysis does not demonstrate causal mediation.

The source launch PDF is used only for fellowship context and is not redistributed. No partner endorsement is implied. A private code repository is not a substitute for appropriate governance of real location or genomic data.

## Source references

See the report's references for the Bio-HEAT protocol, HE2AT programme description, AWI-Gen profile, heat-pregnancy reviews and ERA5-Land documentation. No real weather data were downloaded for this exercise.
