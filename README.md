# Risk Measurement with GAS and GARCH Models

## Overview
This project studies Value-at-Risk (VaR) and Expected Shortfall (ES) using
Generalized Autoregressive Score (GAS) models, GARCH models, and hybrid
GAS–GARCH specifications. Model performance is evaluated through extensive
backtesting procedures and portfolio applications.

The project combines theoretical modeling, empirical implementation, and
statistical backtesting.

## Models Implemented
- GAS models:
  - Gaussian GAS
  - Student-t GAS
  - One-factor GAS model
- GARCH models
- Hybrid GAS–GARCH model

## Risk Measures
- Value-at-Risk (VaR)
- Expected Shortfall (ES)

## Backtesting Methods
- VaR backtesting:
  - LRuc test
  - LRcc test
  - Dynamic Quantile (DQ) test
  - AE value
- ES backtesting:
  - Conditional Calibration test
  - ES Regression test
  - e-backtesting

## Portfolio Analysis
- Portfolio construction
- Monthly and daily backtesting
- Comparison between GAS and GARCH models
- Fama–French 3-factor and 5-factor regressions

## Repository Structure
- `src/` — code for model estimation and backtesting
- `data/` — input data (if included)
- `figures/` — generated plots and tables
- `latex/` — LaTeX source for the report

## How to Run
1. Run the scripts in `src/` to estimate models and generate results.
2. Figures and tables are saved in `figures/`.
3. The report can be compiled from the `latex/` directory.

## Notes
This project was developed as an academic research project.

