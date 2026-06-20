#' DWSDAR: Dynamic Weighted SDAR for High-Dimensional Survival DTRs
#'
#' @description
#' The DWSDAR package implements the Dynamic Weighted SDAR algorithm for estimating
#' optimal dynamic treatment regimes (DTRs) in high-dimensional accelerated failure time
#' (AFT) models with censored survival data. This version utilizes a joint matrix fitting
#' approach (Formula 6) and backward induction with counterfactual scaling.
#'
#' @details
#' The package provides a comprehensive multi-stage estimation workflow:
#' \enumerate{
#'   \item \code{\link{generate_data}}: Simulate multi-stage survival data with AFT structure.
#'   \item \code{\link{result_estimate}}: Estimate optimal treatment regimes using backward induction,
#'         automatically handling dynamic KM-weighting, joint matrix construction, and counterfactual adjustments.
#' }
#'
#' Key features include:
#' \itemize{
#'   \item Joint matrix fitting for simultaneous selection of main and interaction effects.
#'   \item Backward induction framework supporting multi-stage (j >= 2) decision making.
#'   \item Dynamic weighting robust to heavy censoring under high-dimensional covariates (p >> n).
#' }
#'
#' @import ASDAR
#' @import survival
#' @import stats
"_PACKAGE"
