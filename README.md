# DWSDAR: Dynamic Weighted SDAR for High-Dimensional Survival DTRs

The DWSDAR package implements the Dynamic Weighted SDAR algorithm for estimating optimal dynamic treatment regimes in high-dimensional accelerated failure time models with censored survival outcomes, as described in Zhao & Zhao. This method is specifically designed for multi-stage clinical decision making with high-dimensional covariates.

## Key Features

- High-dimensional capability: Handles covariate dimensions (p) much larger than sample size (n)
- Multi-stage optimization: Supports J-stage treatment decision processes
- Censored survival outcomes: Robust to high censoring rates (up to 45%)
- Accurate variable selection: ℓ₀-regularization via SDAR algorithm
- Clinically interpretable: Direct covariate effect interpretation through AFT models
- Theoretical guarantees: Exponential convergence to O(√(log p/n)) error bounds
