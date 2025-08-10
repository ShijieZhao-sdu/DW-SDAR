library(DWSDAR)

data <- generate_data(
  n = 500,
  p = 1000,
  k = 10,
  j = 2,
  rho = 0.3,
  sigma = 1,
  c.r = 0.3,
  xi_probs = c(0.9)
)

processed <- process_data(data, n_total = 500)
results <- result_estimate(processed)

for (stage in 1:2) {
  sort_idx <- processed$sort_idx[[stage]]
  true_opt_sorted <- data$Aopt[[stage]][sort_idx]
  pred_opt_sorted <- results$pred_opt[[stage]]
  observed_idx <- which(processed$xi_sorted[[stage]] == 1)

  if (length(observed_idx) > 0) {
    accuracy <- mean(
      pred_opt_sorted[observed_idx] == true_opt_sorted[observed_idx],
      na.rm = TRUE
    )
    cat(sprintf("Stage %d accuracy: %.2f%%\n", stage, accuracy * 100))
  } else {
    cat(sprintf("Stage %d: No observed subjects\n", stage))
  }
}
