#' Exponential Inverse Logit (Expit) Function
#'
#' Computes the inverse logit (expit) transformation: expit(x) = 1 / (1 + exp(-x))
#'
#' @param x Input numeric vector
#' @return Transformed values in (0,1) range
#'
#' @examples
#' expit(c(-1, 0, 1))
#'
#' @export
expit <- function(x) {
  1.0 / (1.0 + exp(-x))
}
#' Generate Simulated Data for Multi-stage Survival Analysis
#'
#' Generates simulation data for evaluating DW-SDAR algorithm under high-dimensional
#' multi-stage accelerated failure time models with censored survival outcomes.
#'
#' @param n Number of subjects
#' @param p Number of covariates
#' @param k Number of non-zero coefficients
#' @param j Number of treatment stages
#' @param rho Covariate correlation parameter (default=0.3)
#' @param sigma Error term standard deviation (default=1.0)
#' @param c.r Censoring rate (default=0.3)
#' @param xi_probs Probability of entering subsequent stages (default=rep(0.9, j-1))
#' @param cens_dependent Whether censoring depends on covariates (default=FALSE)
#' @param same_x Whether all stages share the same covariate matrix (default=FALSE)
#'
#' @return List containing:
#' \itemize{
#'   \item{X: List of covariate matrices per stage}
#'   \item{A: List of treatment assignments per stage}
#'   \item{Aopt: List of optimal treatments per stage}
#'   \item{Y: List of survival times per stage}
#'   \item{delta: List of censoring indicators per stage}
#'   \item{xi: List of stage entry indicators}
#'   \item{Beta1: List of main effect coefficients}
#'   \item{Beta2: List of treatment interaction coefficients}
#' }
#'
#' @examples
#' # Independent censoring, different covariates per stage
#' data1 <- generate_data(n=500, p=1000, k=10, j=2)
#'
#' # Censoring depends on covariates, same covariates all stages
#' data2 <- generate_data(n=500, p=1000, k=10, j=2,
#'                       cens_dependent=TRUE, same_x=TRUE)
#'
#' @export

generate_data <- function(n, p, k, j, rho = 0.3, sigma = 1.0, c.r = 0.3,
                          xi_probs = rep(0.9, j-1), cens_dependent = FALSE,
                          same_x = FALSE) {
  expit <- function(x) { 1.0 / (1.0 + exp(-x)) }

  # Initialize storage lists
  X_list <- vector("list", j)
  A_list <- vector("list", j)
  Aopt_list <- vector("list", j)
  Y_list <- vector("list", j)
  delta_list <- vector("list", j)
  xi_list <- vector("list", j)
  Beta1_list <- vector("list", j)
  Beta2_list <- vector("list", j)

  # ====== 1. Covariate Generation ======
  generate_covariates <- function(n, p, rho) {
    X_b <- matrix(rnorm(n * p), n, p)
    X <- matrix(NA, n, p)
    X[1, ] <- X_b[1, ]
    X[n, ] <- X_b[n, ]
    if (n > 2) {
      for (i in 2:(n-1)) {
        X[i, ] <- X_b[i, ] + rho * (X_b[i-1, ] + X_b[i+1, ])
      }
    }
    return(X)
  }

  # Generate base covariate matrix
  base_X <- generate_covariates(n, p, rho)

  # Replicate or generate new covariates per stage
  for (stage in 1:j) {
    if (same_x) {
      X_list[[stage]] <- base_X  # Reuse same matrix
    } else {
      X_list[[stage]] <- generate_covariates(n, p, rho)  # New matrix
    }
  }

  # ====== 2. Coefficient Generation ======
  for (stage in 1:j) {
    beta1 <- rep(0, p)
    beta2 <- rep(0, p)
    non_zero_indices1 <- sample(1:p, k)
    non_zero_indices2 <- sample(1:p, k)
    beta1[non_zero_indices1] <- runif(k, min = sigma * sqrt(20 * log(p)/n),
                                      max = 50 * sigma * sqrt(2 * log(p)/n))
    beta2[non_zero_indices2] <- runif(k, min = sigma * sqrt(20 * log(p)/n),
                                      max = 50 * sigma * sqrt(2 * log(p)/n))
    Beta1_list[[stage]] <- beta1
    Beta2_list[[stage]] <- beta2
  }

  # ====== 3. Censoring Generation ======
  if (cens_dependent) {
    # Covariate-dependent censoring
    l <- runif(k, min = -1, max = 1)
    non_zero_idx <- sample(1:p, k)
    # Use first stage covariates for censoring generation
    log_odds <- X_list[[1]][, non_zero_idx, drop = FALSE] %*% l
    prob_cens <- expit(log_odds + c.r)
    delta_base <- rbinom(n, 1, prob = 1 - prob_cens)
  } else {
    # Independent censoring (original)
    delta_base <- rbinom(n, 1, prob = 1 - c.r)
  }

  # Apply same censoring indicator to all stages
  for (stage in 1:j) {
    delta_list[[stage]] <- delta_base
  }

  # ====== 4. Stage 1 Data Generation ======
  stage <- 1
  c <- runif(k, min = -1, max = 1)
  non_zero_idx <- sample(1:p, k)
  A_prob <- expit(X_list[[stage]][, non_zero_idx, drop = FALSE] %*% c)
  A_list[[stage]] <- rbinom(n, 1, A_prob)
  Aopt_list[[stage]] <- ifelse(X_list[[stage]] %*% Beta2_list[[stage]] > 0, 1, 0)
  Y_list[[stage]] <- X_list[[stage]] %*% Beta1_list[[stage]] +
    A_list[[stage]] * (X_list[[stage]] %*% Beta2_list[[stage]]) +
    rnorm(n, sd = 0.3)

  # ====== 5. Subsequent Stages ======
  xi_list[[1]] <- rep(1, n)
  current_idx <- 1:n

  if (j >= 2) {
    for (stage in 2:j) {
      xi_vec <- rep(0, n)
      xi_prob <- xi_probs[stage-1]

      if (length(current_idx) > 0) {
        xi_temp <- rbinom(length(current_idx), 1, prob = xi_prob)
        xi_vec[current_idx] <- xi_temp
        valid_idx <- current_idx[xi_temp == 1]
      } else {
        valid_idx <- integer(0)
      }

      A_stage <- rep(NA, n)
      Aopt_stage <- rep(NA, n)
      Y_stage <- rep(NA, n)

      if (length(valid_idx) > 0) {
        d <- runif(k, min = -1, max = 1)
        non_zero_idx <- sample(1:p, k)
        A_prob <- expit(X_list[[stage]][valid_idx, non_zero_idx, drop = FALSE] %*% d)
        A_stage[valid_idx] <- rbinom(length(valid_idx), 1, A_prob)

        Aopt_stage[valid_idx] <- ifelse(
          X_list[[stage]][valid_idx, , drop = FALSE] %*% Beta2_list[[stage]] > 0, 1, 0
        )

        Y_stage[valid_idx] <- X_list[[stage]][valid_idx, , drop = FALSE] %*% Beta1_list[[stage]] +
          A_stage[valid_idx] * (X_list[[stage]][valid_idx, , drop = FALSE] %*% Beta2_list[[stage]]) +
          rnorm(length(valid_idx), sd = 0.3)
      }

      A_list[[stage]] <- A_stage
      Aopt_list[[stage]] <- Aopt_stage
      Y_list[[stage]] <- Y_stage
      xi_list[[stage]] <- xi_vec
      current_idx <- valid_idx
    }
  }

  return(list(
    X = X_list,
    A = A_list,
    Aopt = Aopt_list,
    Y = Y_list,
    delta = delta_list,
    xi = xi_list,
    Beta1 = Beta1_list,
    Beta2 = Beta2_list
  ))
}
#' Process Data for DW-SDAR Algorithm
#'
#' Prepares data for DW-SDAR estimation by sorting, weighting, and standardizing
#' according to the dynamic weighting scheme described in Zhao & Zhao (2025).
#'
#' @param data Output from generate_data()
#' @param n_total Total sample size (same as n in generate_data)
#'
#' @return Processed data list:
#' \itemize{
#'   \item{Yw: Weighted response vectors per stage}
#'   \item{Xstar: Standardized design matrices per stage}
#'   \item{D: Scaling matrices per stage}
#'   \item{A_sorted: Sorted treatment assignments}
#'   \item{xi_sorted: Sorted stage indicators}
#'   \item{X_original_sorted: Sorted covariate matrices}
#'   \item{sort_idx: Sorting indices}
#' }
#'
#' @export
process_data <- function(data, n_total) {
  j <- length(data$X)
  Yw_list <- vector("list", j)
  Xstar_list <- vector("list", j)
  D_list <- vector("list", j)
  A_sorted_list <- vector("list", j)
  xi_sorted_list <- vector("list", j)
  X_original_sorted_list <- vector("list", j)
  sort_idx_list <- vector("list", j)

  for (stage in 1:j) {
    X <- data$X[[stage]]
    Y <- data$Y[[stage]]
    A <- data$A[[stage]]
    xi <- data$xi[[stage]]
    delta <- data$delta[[stage]]

    sort_idx <- order(Y)
    sort_idx_list[[stage]] <- sort_idx

    X_sorted <- X[sort_idx, , drop = FALSE]
    Y_sorted <- Y[sort_idx]
    A_sorted <- A[sort_idx]
    xi_sorted <- xi[sort_idx]
    delta_sorted <- delta[sort_idx]

    n_stage <- length(Y_sorted)
    w <- rep(NA, n_stage)

    w[1] <- delta_sorted[1] / n_total
    if (n_stage > 1) {
      for (i in 2:n_stage) {
        prod_term <- 1
        if (i > 2) {
          for (l in 1:(i-1)) {
            ratio <- (n_total - l) / (n_total - l + 1)
            prod_term <- prod_term * (ratio ^ delta_sorted[l])
          }
        }
        w[i] <- (delta_sorted[i] / (n_total - i + 1)) * prod_term
      }
    }

    X_weighted <- sqrt(w) * X_sorted
    Y_weighted <- sqrt(w) * Y_sorted

    col_norms <- sqrt(colSums(X_weighted^2))
    D <- diag(sqrt(n_total) / col_norms)
    X_star <- X_weighted %*% D

    Yw_list[[stage]] <- Y_weighted
    Xstar_list[[stage]] <- X_star
    D_list[[stage]] <- D
    A_sorted_list[[stage]] <- A_sorted
    xi_sorted_list[[stage]] <- xi_sorted
    X_original_sorted_list[[stage]] <- X_sorted
  }

  return(list(
    Yw = Yw_list,
    Xstar = Xstar_list,
    D = D_list,
    A_sorted = A_sorted_list,
    xi_sorted = xi_sorted_list,
    X_original_sorted = X_original_sorted_list,
    sort_idx = sort_idx_list
  ))
}

#' Estimate Optimal DTRs with DW-SDAR
#'
#' Implements the Dynamic Weighted SDAR algorithm for estimating optimal
#' dynamic treatment regimes in high-dimensional AFT models.
#'
#' @param processed_data Output from process_data()
#' @param varr2 ASDAR tuning parameter (default=0.01)
#' @param tau ASDAR step size (default=10)
#' @param tau1 ASDAR regularization (default=1)
#' @param iter_max Maximum iterations (default=20)
#'
#' @return Estimation results:
#' \itemize{
#'   \item{Beta: Main effect estimates per stage}
#'   \item{Psi: Treatment effect estimates per stage}
#'   \item{pred_opt: Predicted optimal treatments}
#' }
#'
#' @importFrom ASDAR Asdar
#' @export
result_estimate <- function(processed_data, varr2 = 0.01, tau = 10, tau1 = 1, iter_max = 20) {
  j <- length(processed_data$Xstar)
  Beta_list <- vector("list", j)
  Psi_list <- vector("list", j)
  pred_opt_list <- vector("list", j)

  for (stage in 1:j) {
    X_star <- processed_data$Xstar[[stage]]
    Yw <- processed_data$Yw[[stage]]
    D <- processed_data$D[[stage]]
    A_sorted <- processed_data$A_sorted[[stage]]
    xi_sorted <- processed_data$xi_sorted[[stage]]
    X_original_sorted <- processed_data$X_original_sorted[[stage]]

    p <- ncol(X_star)
    ita0 <- rep(0, p)

    idx_control <- which(A_sorted == 0 & xi_sorted == 1)
    if (length(idx_control) == 0) {
      beta_hat <- rep(0, p)
    } else {
      beta_fit <- Asdar(
        X_star[idx_control, , drop = FALSE],
        Yw[idx_control],
        varr2,
        ita0,
        tau,
        tau1,
        D,
        iter_max
      )
      beta_hat <- beta_fit[[1]]
    }

    residuals <- Yw - X_star %*% beta_hat

    idx_treatment <- which(A_sorted == 1 & xi_sorted == 1)
    if (length(idx_treatment) == 0) {
      psi_hat <- rep(0, p)
    } else {
      psi_fit <- Asdar(
        X_star[idx_treatment, , drop = FALSE],
        residuals[idx_treatment],
        varr2,
        ita0,
        tau,
        tau1,
        D,
        iter_max
      )
      psi_hat <- psi_fit[[1]]
    }

    beta_orig <- D %*% beta_hat
    psi_orig <- D %*% psi_hat

    pred_opt <- rep(NA, length(xi_sorted))
    observed_idx <- which(xi_sorted == 1)

    if (length(observed_idx) > 0) {
      pred_opt[observed_idx] <- ifelse(
        X_original_sorted[observed_idx, , drop = FALSE] %*% psi_orig > 0, 1, 0
      )
    }

    Beta_list[[stage]] <- beta_orig
    Psi_list[[stage]] <- psi_orig
    pred_opt_list[[stage]] <- pred_opt
  }

  return(list(
    Beta = Beta_list,
    Psi = Psi_list,
    pred_opt = pred_opt_list
  ))
}
