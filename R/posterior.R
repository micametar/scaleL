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
#' One of three likelihood methods (`profile`, `full`, `empirical`) selected
#' via [get_loglik_fn()]. The profile method is the fast feasibility screen:
#' it uses the asymptotic sampling distributions `xbar ~ N(mu, sigma^2/n)`
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

#' Full 2-D quadrature log-likelihood for a single measure under candidate L
#'
#' The parametric likelihood method: integrates over `(mu, V*)` on an
#' `n_mu x n_v` grid under priors `mu ~ Uniform(scale_min, b)` and
#' `V* ~ Beta(2, 2)` truncated to the beta-binomial-representable region.
#' Slower than the profile screen but produces continuous discrimination
#' among feasible L values.
#'
#' @inheritParams loglik_profile_single
#' @param n_mu,n_v Quadrature grid resolution (default 25 each).
#' @return Scalar log-likelihood; `-Inf` if L is incompatible.
#' @keywords internal
loglik_full_single <- function(xbar, s2, n, L, scale_min = 1,
                               n_mu = 25, n_v = 25) {
  b <- scale_min + L - 1
  if (xbar <= scale_min || xbar >= b) return(-Inf)
  mu_grid <- seq(scale_min + 0.001, b - 0.001, length.out = n_mu)
  vstar_grid <- seq(0.001, 0.999, length.out = n_v)
  log_prior_v <- function(v) dbeta(v, 2, 2, log = TRUE)
  log_grid <- matrix(-Inf, nrow = n_mu, ncol = n_v)
  for (i in seq_along(mu_grid)) {
    mu <- mu_grid[i]
    sm <- sigma2_max(mu, L, scale_min)
    for (j in seq_along(vstar_grid)) {
      v <- vstar_grid[j]
      sigma2 <- v * sm
      if (L > 2) {
        rho <- (v * (L - 1) - 1) / (L - 2)
        if (rho <= 0 || rho >= 1) next
      }
      ll_xbar <- dnorm(xbar, mean = mu, sd = sqrt(sigma2 / n), log = TRUE)
      q <- (n - 1) * s2 / sigma2
      ll_s2 <- dchisq(q, df = n - 1, log = TRUE) + log((n - 1) / sigma2)
      log_grid[i, j] <- ll_xbar + ll_s2 + log_prior_v(v)
    }
  }
  dmu <- diff(mu_grid[1:2])
  dv  <- diff(vstar_grid[1:2])
  m <- max(log_grid)
  if (!is.finite(m)) return(-Inf)
  log(sum(exp(log_grid - m)) * dmu * dv) + m
}

# Package-internal environment caching the empirical reference for the
# duration of a scaleL() call. Avoids a global <<- assignment.
.scaleL_env <- new.env(parent = emptyenv())

# Default settings for the empirical reference. These match Piers's
# simulation study (20,000 synthetic studies per (L, N) cell on a 40x40
# histogram, Laplace smoothing 0.5, rho ~ Beta(2, 5)). R can be shrunk via
# the `emp_R` argument of scaleL() to keep examples/tests fast.
EMP_N_GRID  <- c(50, 100, 200, 500, 1000)
EMP_R       <- 20000L
EMP_NBINS   <- 40L
EMP_LAPLACE <- 0.5
EMP_RHO_A   <- 2
EMP_RHO_B   <- 5
EMP_MU_LO   <- 0.05
EMP_MU_HI   <- 0.95

#' Build the empirical reference distribution for the empirical likelihood
#'
#' Builds, by internal Monte-Carlo simulation (beta-binomial draws via
#' [extraDistr::rbbinom()]), a reference set of 2-D histograms over
#' `(xbar, s^2)` for each candidate L at a grid of anchor sample sizes. The
#' reference does NOT depend on the user's corpus, so the empirical
#' likelihood is always available. Used by [loglik_empirical_single()].
#'
#' @param L_grid Integer vector of candidate L values.
#' @param N_grid Anchor sample sizes for the histograms.
#' @param R Number of synthetic studies per (L, N) cell. Lower values build
#'   faster (use for examples/tests); the default matches the manuscript.
#' @param nbins Histogram resolution per axis.
#' @param laplace Laplace smoothing pseudocount per bin.
#' @param rho_a,rho_b Beta parameters for the simulated intraclass correlation.
#' @param mu_lo,mu_hi Range of the simulated standardized mean.
#' @param scale_min Lowest scale value (default 1).
#' @param verbose Logical; print progress.
#' @return A list keyed by candidate L (as character); each entry holds a
#'   list of per-N log-density histograms.
#' @export
build_empirical_reference <- function(L_grid = 4:12,
                                       N_grid = EMP_N_GRID,
                                       R       = EMP_R,
                                       nbins   = EMP_NBINS,
                                       laplace = EMP_LAPLACE,
                                       rho_a   = EMP_RHO_A,
                                       rho_b   = EMP_RHO_B,
                                       mu_lo   = EMP_MU_LO,
                                       mu_hi   = EMP_MU_HI,
                                       scale_min = 1,
                                       verbose = FALSE) {
  if (verbose) {
    message("Building empirical reference (R = ", R, " per L x N cell) ...")
  }
  REF <- list()
  for (L in L_grid) {
    b <- scale_min + L - 1
    xbar_edges <- seq(scale_min, b, length.out = nbins + 1)
    s2_max_overall <- sigma2_max((scale_min + b) / 2, L, scale_min)
    s2_edges <- seq(0, s2_max_overall, length.out = nbins + 1)
    bin_area <- diff(xbar_edges)[1] * diff(s2_edges)[1]
    cells <- vector("list", length(N_grid))
    for (i_N in seq_along(N_grid)) {
      N <- N_grid[i_N]
      xbars <- numeric(R); s2s <- numeric(R)
      n_valid <- 0L
      attempts <- 0L
      while (n_valid < R) {
        attempts <- attempts + 1L
        if (attempts > 5L * R) break
        mu_tilde <- runif(1, mu_lo, mu_hi)
        rho <- rbeta(1, rho_a, rho_b)
        if (rho <= 0 || rho >= 1) next
        x <- scale_min + extraDistr::rbbinom(
          N, size = L - 1,
          alpha = mu_tilde * (1 - rho) / rho,
          beta  = (1 - mu_tilde) * (1 - rho) / rho)
        n_valid <- n_valid + 1L
        xbars[n_valid] <- mean(x)
        s2s[n_valid]   <- var(x)
      }
      ix <- pmin(pmax(findInterval(xbars[seq_len(n_valid)], xbar_edges), 1L),
                 nbins)
      iy <- pmin(pmax(findInterval(s2s[seq_len(n_valid)],  s2_edges),  1L),
                 nbins)
      counts <- matrix(0, nrow = nbins, ncol = nbins)
      for (r in seq_len(n_valid)) counts[ix[r], iy[r]] <- counts[ix[r], iy[r]] + 1
      smoothed <- (counts + laplace) /
                  ((n_valid + nbins^2 * laplace) * bin_area)
      cells[[i_N]] <- list(
        log_density = log(smoothed),
        xbar_edges  = xbar_edges,
        s2_edges    = s2_edges,
        nbins       = nbins
      )
    }
    REF[[as.character(L)]] <- list(cells = cells, N_grid = N_grid, L = L)
  }
  REF
}

#' Empirical log-likelihood for a single measure under candidate L
#'
#' Looks up the 2-D reference histogram density at the observed
#' `(xbar, s^2)` point, log-linearly interpolating across the two bracketing
#' N anchors. The deterministic Bhatia-Davis feasibility floor is restored on
#' top of the smoothed density. The reference is supplied via `REF` (defaults
#' to the package-internal cache populated by [scaleL()]).
#'
#' @inheritParams loglik_profile_single
#' @param REF Empirical reference list (see [build_empirical_reference()]).
#' @return Scalar log-likelihood; `-Inf` if L is incompatible.
#' @keywords internal
loglik_empirical_single <- function(xbar, s2, n, L, scale_min = 1,
                                     REF = .scaleL_env$empirical_reference) {
  b <- scale_min + L - 1
  if (xbar <= scale_min || xbar >= b) return(-Inf)
  sigma2_mle <- (n - 1) / n * s2
  if (sigma2_mle >= sigma2_max(xbar, L, scale_min)) return(-Inf)
  if (is.null(REF)) {
    stop("Empirical likelihood requested but no reference is available. ",
         "Build one with build_empirical_reference() or call scaleL() with ",
         "lik_method = 'empirical'.")
  }
  Lref <- REF[[as.character(L)]]
  if (is.null(Lref)) return(-Inf)
  cell0 <- Lref$cells[[1]]
  ix <- pmin(pmax(findInterval(xbar, cell0$xbar_edges), 1L), cell0$nbins)
  iy <- pmin(pmax(findInterval(s2,   cell0$s2_edges),   1L), cell0$nbins)
  N_grid <- Lref$N_grid
  if (n <= N_grid[1]) {
    i_N_lo <- i_N_hi <- 1L; w_lo <- 1; w_hi <- 0
  } else if (n >= N_grid[length(N_grid)]) {
    i_N_lo <- i_N_hi <- length(N_grid); w_lo <- 0; w_hi <- 1
  } else {
    i_N_hi <- which(N_grid >= n)[1]
    i_N_lo <- i_N_hi - 1L
    w_hi <- (log(n) - log(N_grid[i_N_lo])) /
            (log(N_grid[i_N_hi]) - log(N_grid[i_N_lo]))
    w_lo <- 1 - w_hi
  }
  ld_lo <- Lref$cells[[i_N_lo]]$log_density[ix, iy]
  ld_hi <- Lref$cells[[i_N_hi]]$log_density[ix, iy]
  w_lo * ld_lo + w_hi * ld_hi
}

#' Dispatch on likelihood method
#'
#' Returns the per-measure log-likelihood function for the chosen method.
#' `empirical` (default in [scaleL()]) is the operational reference-histogram
#' lookup; `full` is the 2-D quadrature parametric likelihood; `profile` is
#' the fast feasibility screen.
#'
#' @param method One of `"empirical"`, `"full"`, `"profile"`.
#' @return A per-measure log-likelihood function.
#' @keywords internal
get_loglik_fn <- function(method = c("empirical", "full", "profile")) {
  method <- match.arg(method)
  switch(method,
         "profile"   = loglik_profile_single,
         "full"      = loglik_full_single,
         "empirical" = loglik_empirical_single)
}

#' Posterior over candidate L for a study
#'
#' Combines the per-measure log-likelihoods (assuming within-study
#' uniformity of L) with a prior over candidate L to produce a normalized
#' posterior. The likelihood method is supplied via `loglik_fn`; by default
#' the empirical reference-histogram likelihood is used (see
#' [get_loglik_fn()]). When every L is structurally rejected by the
#' likelihood, the posterior falls back to the prior.
#'
#' @param measures A data.frame with columns `mean`, `sd`, `n` (one row
#'   per measure within the study).
#' @param prior Named numeric vector indexed by L (as character).
#' @param L_grid Integer vector of candidate L values.
#' @param scale_min Lowest scale value (default 1).
#' @param loglik_fn Per-measure log-likelihood function (default
#'   [get_loglik_fn()], i.e. the empirical likelihood). For the empirical
#'   method a reference must be available; use `lik_method = "profile"` or
#'   `"full"` for a reference-free posterior.
#' @return Named numeric vector summing to 1, indexed by `as.character(L_grid)`.
#' @examples
#' meas <- data.frame(mean = 3.5, sd = 0.8, n = 100)
#' prior <- resolve_prior("tier1", NULL, 4:12)
#' # profile likelihood needs no Monte-Carlo reference:
#' posterior_L(meas, prior, 4:12, loglik_fn = scaleL:::loglik_profile_single)
#' @export
posterior_L <- function(measures, prior, L_grid, scale_min = 1,
                        loglik_fn = get_loglik_fn()) {
  log_lik_study <- vapply(L_grid, function(L) {
    sum(vapply(seq_len(nrow(measures)), function(i) {
      loglik_fn(
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
