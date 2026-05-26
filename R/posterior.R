#' Bhatia-Davis variance bound on a bounded scale
#'
#' Maximum admissible variance for a random variable taking values in
#' `[scale_min, scale_min + L - 1]` with mean `mu`.
#'
#' @param mu Numeric, the mean.
#' @param L Integer, the scale length.
#' @param scale_min Numeric, the lowest scale value (default 1).
#' @return The variance bound `(mu - scale_min) * (scale_max - mu)`.
#' @keywords internal
sigma2_max <- function(mu, L, scale_min = 1) {
  b <- scale_min + L - 1
  (mu - scale_min) * (b - mu)
}

#' Profile log-likelihood for a single measure under candidate L
#'
#' Uses the asymptotic sampling distributions `xbar ~ N(mu, sigma^2/n)`
#' and `(n-1) s^2 / sigma^2 ~ chi^2(n-1)` at the MLE values of
#' `(mu, sigma^2)`, and rejects L values structurally incompatible with
#' the observation (variance exceeds the Bhatia-Davis bound, or the
#' implied beta-binomial intraclass correlation falls outside `(0, 1)`).
#'
#' @param xbar Sample mean.
#' @param s2 Sample variance.
#' @param n Sample size.
#' @param L Candidate scale length.
#' @param scale_min Lowest scale value (default 1).
#' @return Scalar log-likelihood; `-Inf` if L is incompatible.
#' @keywords internal
loglik_profile_single <- function(xbar, s2, n, L, scale_min = 1) {
  b <- scale_min + L - 1
  if (xbar <= scale_min || xbar >= b) return(-Inf)
  sm <- sigma2_max(xbar, L, scale_min)
  sigma2_mle <- (n - 1) / n * s2
  if (sigma2_mle >= sm) return(-Inf)
  if (L > 2) {
    vstar <- sigma2_mle / sm
    rho <- (vstar * (L - 1) - 1) / (L - 2)
    if (rho <= 0 || rho >= 1) return(-Inf)
  }
  sigma2 <- s2
  ll_xbar <- dnorm(xbar, mean = xbar, sd = sqrt(sigma2 / n), log = TRUE)
  q <- (n - 1) * s2 / sigma2
  ll_s2 <- dchisq(q, df = n - 1, log = TRUE) + log((n - 1) / sigma2)
  ll_xbar + ll_s2
}

#' Posterior over candidate L for a study
#'
#' Combines the per-measure profile log-likelihoods (assuming within-study
#' uniformity of L) with a prior over candidate L to produce a normalized
#' posterior. When every L is structurally rejected by the likelihood, the
#' posterior falls back to the prior.
#'
#' @param measures A data.frame with columns `mean`, `sd`, `n` (one row
#'   per measure within the study).
#' @param prior Named numeric vector indexed by L (as character).
#' @param L_grid Integer vector of candidate L values.
#' @param scale_min Lowest scale value (default 1).
#' @return Named numeric vector summing to 1, indexed by `as.character(L_grid)`.
#' @examples
#' meas <- data.frame(mean = 3.5, sd = 0.8, n = 100)
#' prior <- resolve_prior("tier1", NULL, 2:11)
#' posterior_L(meas, prior, 2:11)
#' @export
posterior_L <- function(measures, prior, L_grid, scale_min = 1) {
  log_lik_study <- vapply(L_grid, function(L) {
    sum(vapply(seq_len(nrow(measures)), function(i) {
      loglik_profile_single(
        xbar = measures$mean[i],
        s2   = measures$sd[i]^2,
        n    = measures$n[i],
        L    = L,
        scale_min = scale_min
      )
    }, numeric(1)))
  }, numeric(1))

  if (all(!is.finite(log_lik_study))) {
    return(setNames(as.numeric(prior[as.character(L_grid)]) /
                      sum(prior[as.character(L_grid)]),
                    as.character(L_grid)))
  }
  log_prior <- log(as.numeric(prior[as.character(L_grid)]))
  log_post <- log_lik_study + log_prior
  log_post[!is.finite(log_post)] <- -Inf
  m <- max(log_post)
  post <- exp(log_post - m)
  post <- post / sum(post)
  setNames(post, as.character(L_grid))
}
