#' @rawNamespace import(ASDAR, survival)
#' DWSDAR: Dynamic Weighted SDAR for High-Dimensional Survival DTRs
#'
#' @description
#' The DWSDAR package implements the Dynamic Weighted SDAR algorithm for estimating
#' optimal dynamic treatment regimes in high-dimensional accelerated failure time
#' models with censored survival data. The method is particularly designed for
#' multi-stage clinical decision making with high-dimensional covariates.
#'
#' @details
#' The package provides a complete workflow:
#' \enumerate{
#'   \item \code{\link{generate_data}}: Simulate multi-stage survival data
#'   \item \code{\link{process_data}}: Prepare data for DW-SDAR algorithm
#'   \item \code{\link{result_estimate}}: Estimate optimal treatment regimes
#' }
#'
#' Key features include:
#' \itemize{
#'   \item Handling of high-dimensional covariates (p >> n)
#'   \item Support for multi-stage treatment regimes
#'   \item Robustness to high censoring rates (up to 45\%)
#'   \item ℓ₀-regularization for accurate variable selection
#' }
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{generate_data}}}{Generate simulated survival data}
#'   \item{\code{\link{process_data}}}{Prepare data for analysis}
#'   \item{\code{\link{result_estimate}}}{Estimate optimal DTRs}
#'   \item{\code{\link{expit}}}{Inverse logit transformation}
#' }
#'
#' @docType package
#' @name DWSDAR-package
NULL
