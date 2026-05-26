#' Per-study standardized mean and sampling variance under known L
#'
#' Computes `theta(L) = mean over measures of (mean - scale_min) / (L - 1)`
#' and its sampling variance, assuming within-study i.i.d. measures.
#'
#' @param measures Data.frame with `mean`, `sd`, `n`.
#' @param L Scale length.
#' @param scale_min Lowest scale value.
#' @return List with elements `theta` and `v`.
#' @keywords internal
study_theta_and_var <- function(measures, L, scale_min = 1) {
  if (L <= 1) return(list(theta = NA_real_, v = NA_real_))
  std_means <- (measures$mean - scale_min) / (L - 1)
  std_vars  <- (measures$sd / (L - 1))^2 / measures$n
  theta <- mean(std_means)
  v <- mean(std_vars) / nrow(measures)
  list(theta = theta, v = v)
}
