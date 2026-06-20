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

data <- generate_data(n, p, k, j, rho, sigma, c.r, xi_probs,
                      cens_dependent = FALSE, same_x = FALSE)

estimates <- result_estimate(data)

for (stage in 1:j) {
  true_opt <- data$Aopt[[stage]]
  pred_opt <- estimates$pred_opt[[stage]]
  xi       <- data$xi[[stage]]

  observed_idx <- which(xi == 1)

  if (length(observed_idx) > 0) {
    true_opt_obs <- true_opt[observed_idx]
    pred_opt_obs <- pred_opt[observed_idx]

    accuracy <- mean(pred_opt_obs == true_opt_obs, na.rm = TRUE)
    cat(sprintf("Stage %d Accuracy: %.4f\n", stage, accuracy))
  } else {
    cat(sprintf("Stage %d: no observed subjects\n", stage))
  }
}
