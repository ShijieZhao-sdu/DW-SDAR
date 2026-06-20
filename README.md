# DWSDAR: Dynamic Weighted SDAR for High-Dimensional Survival DTRs

`DWSDAR` is an R package implementing the Dynamic Weighted Sparse Differential Active Set Regularization (Dynamic Weighted SDAR) algorithm. This methodology is specifically engineered to estimate optimal dynamic treatment regimes (DTRs) under high-dimensional accelerated failure time (AFT) models using multi-stage clinical data with right-censored survival outcomes.

By utilizing a **joint matrix fitting approach** coupled with a **backward induction framework** featuring counterfactual survival time scaling, the package accurately conducts simultaneous variable selection for both main effects and treatment interactions in sparse, high-dimensional settings.

## 🌟 Key Features

- **High-Dimensional Capability:** Seamlessly handles sparse datasets where the covariate dimension ($p$) is significantly larger than the sample size ($n$).
- **Joint Matrix Fitting:** Merges main effects and interaction terms into a single design matrix, ensuring highly efficient and simultaneous variable selection.
- **Backward Induction with Counterfactual Scaling:** Accommodates multi-stage ($J \ge 2$) sequential decision-making by dynamically adjusting prior-stage log-survival outcomes based on downstream optimal counterfactual decisions.
- **Adaptive Dynamic Weighting:** Incorporates automated Kaplan-Meier-based inverse censoring weights with an integrated truncation mechanism, rendering the algorithm robust against heavy censoring rates (up to 45%).
- **Clinical Interpretability:** Built directly upon the classical AFT model framework, allowing for straightforward, physical interpretations of the estimated regression coefficients.

## 📦 Core Functions

The package provides a highly streamlined user interface, mapping the entire analytical pipeline into two primary functions:

1. **`generate_data(n, p, k, j, ...)`** A comprehensive simulation engine that generates multi-stage survival data matching the high-dimensional AFT model structure, complete with sequential treatment allocations, stage entry indicators (`xi`), and right-censored survival endpoints.
   
2. **`result_estimate(data, ...)`** **The all-in-one DTR estimation engine.** It accepts the multi-stage data list directly and internally manages the sequential backward pipeline: patient filtering $\rightarrow$ survival sorting $\rightarrow$ dynamic KM weighting $\rightarrow$ joint matrix standardization $\rightarrow$ `Asdar` optimization $\rightarrow$ backward counterfactual propagation. It outputs comprehensive parameter estimates and a sequence of optimal treatment decisions (`pred_opt`) natively mapped back to the original subject index layout.

## 🚀 Quick Start Example

```R
library(DWSDAR)

# 1. Configuration Setup
n <- 500       # Sample size
p <- 1000      # Covariate dimension
k <- 10        # True active signals per stage
j <- 3         # Total treatment decision stages

set.seed(123)

# 2. Step 1: Generate Multi-Stage Simulated Data
sim_data <- generate_data(
  n = n, p = p, k = k, j = j, 
  rho = 0.9, sigma = 1, c.r = 0.3, 
  xi_probs = c(0.9, 0.9)
)

# 3. Step 2: Estimate the Optimal Dynamic Treatment Regimes
estimates <- result_estimate(sim_data)

# 4. Step 3: Evaluate Decision Accuracy Across Stages
for (stage in 1:j) {
  true_opt     <- sim_data$Aopt[[stage]]
  pred_opt     <- estimates$pred_opt[[stage]]
  observed_idx <- which(sim_data$xi[[stage]] == 1)

  if (length(observed_idx) > 0) {
    accuracy <- mean(pred_opt[observed_idx] == true_opt[observed_idx], na.rm = TRUE)
    cat(sprintf("Stage %d Optimal Decision Accuracy: %.4f\n", stage, accuracy))
  }
}
