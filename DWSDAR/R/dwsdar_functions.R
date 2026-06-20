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

#' Compute Kaplan-Meier based weights
#'
#' Internal auxiliary function to calculate dynamic weights based on the Kaplan-Meier
#' estimator for censored survival data.
#'
#' @param y Numeric vector of survival times (logged or raw, sorted inside or passed as ordered)
#' @param delta Numeric vector of censoring indicators (1 = observed, 0 = censored)
#' @return A numeric vector of computed weights corresponding to the original input order
#' @keywords internal
compute_km_weights <- function(y, delta) {
  n <- length(y)
  o <- order(y)
  y_s <- y[o]
  d_s <- delta[o]
  w <- numeric(n)

  w[1] <- d_s[1] / n
  if (n > 1) {
    for (i in 2:n) {
      idx_tmp <- 1:(i-1)
      w[i] <- d_s[i] / (n - i + 1) * prod(((n - idx_tmp) / (n - idx_tmp + 1))^d_s[idx_tmp])
    }
  }
  w_orig <- numeric(n)
  w_orig[o] <- w
  return(w_orig)
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
#'   \item{Y: List of survival times (log scale) per stage}
#'   \item{delta: List of censoring indicators per stage}
#'   \item{xi: List of stage entry indicators}
#'   \item{Beta1: List of main effect coefficients}
#'   \item{Beta2: List of treatment interaction coefficients}
#' }
#'
#' @examples
#' data1 <- generate_data(n=200, p=500, k=10, j=2)
#'
#' @export
generate_data <- function(n, p, k, j, rho = 0.3, sigma = 1.0, c.r = 0.3,
                          xi_probs = rep(0.9, j-1), cens_dependent = FALSE,
                          same_x = FALSE) {

  X_list <- vector("list", j)
  A_list <- vector("list", j)
  Aopt_list <- vector("list", j)
  Y_list <- vector("list", j)
  delta_list <- vector("list", j)
  xi_list <- vector("list", j)
  Beta1_list <- vector("list", j)
  Beta2_list <- vector("list", j)

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

  base_X <- generate_covariates(n, p, rho)

  for (stage in 1:j) {
    if (same_x) {
      X_list[[stage]] <- base_X
    } else {
      X_list[[stage]] <- generate_covariates(n, p, rho)
    }
  }

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

  if (cens_dependent) {
    l <- runif(k, min = -1, max = 1)
    non_zero_idx <- sample(1:p, k)
    log_odds <- X_list[[1]][, non_zero_idx, drop = FALSE] %*% l
    prob_cens <- expit(log_odds + c.r)
    delta_base <- rbinom(n, 1, prob = 1 - prob_cens)
  } else {
    delta_base <- rbinom(n, 1, prob = 1 - c.r)
  }

  for (stage in 1:j) {
    delta_list[[stage]] <- delta_base
  }

  stage <- 1
  c <- runif(k, min = -1, max = 1)
  non_zero_idx <- sample(1:p, k)
  A_prob <- expit(X_list[[stage]][, non_zero_idx, drop = FALSE] %*% c)
  A_list[[stage]] <- rbinom(n, 1, A_prob)
  Aopt_list[[stage]] <- ifelse(X_list[[stage]] %*% Beta2_list[[stage]] > 0, 1, 0)
  Y_list[[stage]] <- X_list[[stage]] %*% Beta1_list[[stage]] +
    A_list[[stage]] * (X_list[[stage]] %*% Beta2_list[[stage]]) +
    rnorm(n, sd = 0.3)

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

#' Estimate Optimal DTRs with DW-SDAR via Joint Fitting Backward Induction
#'
#' Implements the Dynamic Weighted Joint-Matrix SDAR algorithm for estimating optimal
#' dynamic treatment regimes in high-dimensional AFT models using backward induction.
#'
#' @param data A list containing raw components from generate_data() or structured real data (X, A, Y, delta, xi)
#' @param varr2 ASDAR tuning parameter (default=0.01)
#' @param tau ASDAR step size (default=10)
#' @param tau1 ASDAR regularization (default=1)
#' @param iter_max Maximum iterations for ASDAR (default=100)
#'
#' @return A list of estimation results containing per-stage matrices:
#' \itemize{
#'   \item{intercept: Intercept (mean response) estimates per stage}
#'   \item{beta: Main effect parameter estimates per stage}
#'   \item{psi: Treatment interaction effect parameter estimates per stage}
#'   \item{pred_opt: Predicted optimal treatment assignments per stage}
#' }
#'
#' @importFrom ASDAR Asdar
#' @export
result_estimate <- function(data, varr2 = 0.01, tau = 10, tau1 = 1, iter_max = 100) {
  j <- length(data$X)
  n_total <- nrow(data$X[[1]])

  # Initialize storage lists for output
  intercept_list <- vector("list", j)
  beta_list      <- vector("list", j)
  psi_list       <- vector("list", j)
  pred_opt_list  <- vector("list", j)

  # Create a deep copy of working response list to allow sequential backward counterfactual additions
  Y_working <- data$Y

  # ====== Backward Induction: From Stage j down to Stage 1 ======
  for (stage in j:1) {
    X     <- data$X[[stage]]
    A     <- data$A[[stage]]
    Y     <- Y_working[[stage]]
    delta <- data$delta[[stage]]
    xi    <- data$xi[[stage]]

    # 筛选当前阶段存活且进入随访的受试者 (xi == 1)
    idx_active <- which(xi == 1)
    nn <- length(idx_active)
    if (nn < 2) stop(paste("Stage", stage, "sample size is too small for estimation."))

    X_a     <- X[idx_active, , drop = FALSE]
    A_a     <- A[idx_active]
    Y_a     <- Y[idx_active]
    delta_a <- delta[idx_active]

    # 按照当前生存时间进行升序排列
    sort_idx <- order(Y_a)
    X_s     <- X_a[sort_idx, , drop = FALSE]
    A_s     <- A_a[sort_idx]
    Y_s     <- Y_a[sort_idx]
    delta_s <- delta_a[sort_idx]

    # 构建严格符合文献公式(6)的联合设计矩阵 Z = (X, A * X)
    Z_s <- cbind(X_s, A_s * X_s)

    # 计算当前生存状态下的 KM 动态权重
    w <- compute_km_weights(Y_s, delta_s)
    w <- pmin(w, 10 * mean(w, na.rm = TRUE))  # 权重截断防爆

    # 响应变量中心化
    y_mean <- mean(Y_s)
    y_centered <- Y_s - y_mean

    # 应用权重与矩阵标准化
    X_weighted <- diag(sqrt(w)) %*% Z_s
    Y_weighted <- sqrt(w) * y_centered

    col_sums_sq <- colSums(X_weighted^2)
    col_sums_sq[col_sums_sq < 1e-8] <- 1
    dd <- diag(sqrt(nn / col_sums_sq))
    X_star <- X_weighted %*% dd

    # 调用外部 ASDAR 求解联合参数 theta = (beta, psi)
    ita0 <- rep(0, ncol(Z_s))
    fit_theta <- Asdar(X_star, Y_weighted, varr2, ita0, tau, tau1, dd, iter_max)[[1]]
    theta_original <- dd %*% fit_theta

    # 拆解出主效应参数 beta 和 交互效应参数 psi
    p_stage <- ncol(X_s)
    beta_original <- theta_original[1:p_stage]
    psi_original  <- theta_original[(p_stage + 1):(2 * p_stage)]

    # 估计当前阶段的最优决策方案
    estAopt_active <- ifelse(X_a %*% psi_original > 0, 1, 0)

    # ====== 如果不是第一阶段，向前滚动计算反事实调整项 ======
    if (stage > 1) {
      power_term <- (estAopt_active - A_a) * (X_a %*% psi_original)
      power_term <- pmin(pmax(power_term, -5), 5)  # 边界约束，防指数爆炸

      adjust <- rep(0, n_total)
      # 映射加回至上一阶段对应的受试者生存空间内 (注意：输入Y为log尺度，此处复现 exp(y1) + adjust 逻辑)
      adjust[idx_active] <- exp(Y_working[[stage]][idx_active]) * exp(power_term)
      Y_working[[stage-1]] <- log(exp(Y_working[[stage-1]]) + adjust)
    }

    # 填充并保存结果结构
    intercept_list[[stage]] <- y_mean
    beta_list[[stage]]      <- beta_original
    psi_list[[stage]]       <- psi_original

    pred_opt <- rep(NA, n_total)
    pred_opt[idx_active]    <- estAopt_active
    pred_opt_list[[stage]]  <- pred_opt
  }

  return(list(
    intercept = intercept_list,
    beta = beta_list,
    psi = psi_list,
    pred_opt = pred_opt_list
  ))
}
