library(DWSDAR)

n <- 500
p <- 1000
k <- 10
j <- 3
rho <- 0.9
sigma <- 1
c.r <- 0.3
xi_probs <- c(0.9, 0.9)

set.seed(123)

# Generate simulated data with specified parameters
# cens_dependent = FALSE: Censoring is independent of covariates
# same_x = FALSE: Covariates differ between treatment stages
data <- generate_data(n, p, k, j, rho, sigma, c.r, xi_probs,
                      cens_dependent = FALSE, same_x = FALSE)
processed_data <- process_data(data, n_total = n)

# Estimate optimal treatment regimes using DW-SDAR
estimates <- result_estimate(processed_data)

# Evaluate accuracy for each treatment stage
for (stage in 1:j) {
  sort_idx <- processed_data$sort_idx[[stage]]
  true_opt_sorted <- data$Aopt[[stage]][sort_idx]
  pred_opt_sorted <- estimates$pred_opt[[stage]]
  observed_idx <- which(processed_data$xi_sorted[[stage]] == 1)

  if (length(observed_idx) > 0) {
    true_opt_obs <- true_opt_sorted[observed_idx]
    pred_opt_obs <- pred_opt_sorted[observed_idx]

    accuracy <- mean(pred_opt_obs == true_opt_obs, na.rm = TRUE)
    cat(paste("Stage", stage, "Accuracy:", round(accuracy, 4), "\n"))
  } else {
    cat(paste("Stage", stage, "no observed subjects\n"))
  }
}
