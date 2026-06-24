# =============================================================================
# Across-study-means channel ("Estimator B").
#
# A second, independent likelihood over L for a target study. A REML
# random-effects fit to the full-information anchor studies' standardized means
# yields (mu_hat, tau2_hat), which define a predictive N(mu_hat,
# sqrt(tau2_hat + se(L)^2)); the candidate L whose standardized mean theta(L)
# lands in the dense region of that predictive is favoured. Gated at >= 10
# anchors, with a nested moderator fallback (a level is used only if it itself
# carries >= 10 anchors, else the pooled fit).
# =============================================================================

#' REML estimate of (mu, tau2) for a random-effects meta-analysis
#'
#' Iterates the REML estimating equation, truncated at zero. Matches metafor's
#' REML and the across-study channel in Piers's pipeline.
#'
#' @param theta Standardized means (anchors).
#' @param vi Sampling variances of `theta`.
#' @param max_iter,tol Iteration controls.
#' @return List with `mu`, `tau2`, `k`.
#' @export
fit_REML <- function(theta, vi, max_iter = 200L, tol = 1e-8) {
  k <- length(theta)
  if (k < 2L) return(list(mu = if (k == 1L) theta[1] else NA_real_,
                          tau2 = 0, k = k))
  tau2 <- stats::var(theta)
  for (it in seq_len(max_iter)) {
    w  <- 1 / (vi + tau2)
    mu <- sum(w * theta) / sum(w)
    upd <- sum(w^2 * ((theta - mu)^2 - vi)) / sum(w^2) + 1 / sum(w)
    tau2_new <- max(0, upd)
    if (abs(tau2_new - tau2) < tol * (tau2 + 1e-10)) { tau2 <- tau2_new; break }
    tau2 <- tau2_new
  }
  w  <- 1 / (vi + tau2); mu <- sum(w * theta) / sum(w)
  list(mu = mu, tau2 = tau2, k = k)
}

#' Across-study log-likelihood over the L grid for one target study
#'
#' For each candidate L: `theta(L)` and its sampling se come from the target's
#' own moments; the predictive is `N(mu_hat, sqrt(tau2_hat + se(L)^2))`
#' truncated to `(0, 1)`; the contribution is the log truncated-normal density
#' of `theta(L)`. Candidate L that put `theta(L)` outside `(0, 1)` score `-Inf`,
#' so the channel also enforces the feasibility floor.
#'
#' @param measures Data.frame with `mean`, `sd`, `n` for the target study.
#' @param mu_hat,tau2_hat Anchor REML estimates (see [fit_REML()]).
#' @param L_grid Candidate L values.
#' @param scale_min Lower endpoint (default 1).
#' @return Numeric vector of log-likelihoods, one per `L_grid` entry.
#' @export
loglik_across_vec <- function(measures, mu_hat, tau2_hat, L_grid,
                              scale_min = 1) {
  vapply(L_grid, function(L) {
    tv <- study_theta_and_var(measures, L, scale_min)
    theta <- tv$theta
    if (!is.finite(theta) || theta <= 0 || theta >= 1) return(-Inf)
    psd <- sqrt(tau2_hat + tv$v)
    if (!is.finite(psd) || psd <= 0) return(-Inf)
    Z <- pnorm(1, mu_hat, psd) - pnorm(0, mu_hat, psd)
    if (!is.finite(Z) || Z <= 0) return(-Inf)
    dnorm(theta, mu_hat, psd, log = TRUE) - log(Z)
  }, numeric(1))
}

#' Anchor table of standardized means from full-information studies
#'
#' Each full-information (known-L) study contributes its standardized mean and
#' sampling variance at its own known L, tagged with its moderator level.
#'
#' @param studies Named list of per-study measure data.frames (each with a
#'   single known L in column `L`, plus an optional `level` column).
#' @param L_known Integer vector of each study's known L.
#' @param levels Optional character vector of moderator levels per study.
#' @param scale_min Lower endpoint (default 1).
#' @return Data.frame with `theta`, `v`, `level`.
#' @export
anchor_table <- function(studies, L_known, levels = NULL, scale_min = 1) {
  if (is.null(levels)) levels <- rep(NA_character_, length(studies))
  rows <- lapply(seq_along(studies), function(j) {
    if (is.na(L_known[j])) return(NULL)
    tv <- study_theta_and_var(studies[[j]], L_known[j], scale_min)
    if (!is.finite(tv$theta) || !is.finite(tv$v)) return(NULL)
    data.frame(theta = tv$theta, v = tv$v, level = levels[j])
  })
  do.call(rbind, rows)
}

#' Fit the across-study group prior with optional moderator conditioning
#'
#' @param anc Anchor table from [anchor_table()].
#' @param target_level Moderator level of the target study.
#' @param mode `"marginal"` (pool all anchors) or `"conditioned"` (use the
#'   target's level if it has >= `min_k` anchors, else fall back to pooled).
#' @param min_k Minimum anchors in a level to condition (default 10).
#' @return List with `mu`, `tau2`, `used` (`"pooled"`/`"conditioned"`), `k`;
#'   or `NULL` if fewer than two usable anchors.
#' @export
group_prior_fit <- function(anc, target_level = NA,
                            mode = c("marginal", "conditioned"),
                            min_k = 10L) {
  mode <- match.arg(mode)
  use  <- anc
  used <- "pooled"
  if (mode == "conditioned" && !is.na(target_level)) {
    sub <- anc[!is.na(anc$level) & anc$level == target_level, , drop = FALSE]
    if (nrow(sub) >= min_k) { use <- sub; used <- "conditioned" }
  }
  if (is.null(use) || nrow(use) < 2L) return(NULL)
  fit <- fit_REML(use$theta, use$v)
  list(mu = fit$mu, tau2 = fit$tau2, used = used, k = nrow(use))
}
