# =============================================================================
# Bound-normalized standard-deviation imputation for bounded (Likert) scales.
#
# On a bounded response scale [lo, hi] the attainable variance is capped by the
# Bhatia-Davis bound (mu - lo)(hi - mu). The observed SD is recast as a
# dispersion fraction phi = sd / sqrt(bound) in (0, 1), which removes the
# deterministic dependence of the raw SD on mean position and scale length. A
# log-linear model log phi ~ mean_position + log(items) is fit on the studies
# that report an SD, and a missing SD is imputed by predicting phi and
# re-expanding through the target's own ceiling: sd_hat = phi_hat * sqrt(bound).
# Imputation is restricted to known-L anchors (circularity firewall).
#
# NOTE: this module's variance ceiling helper is named `.sd_bound` to avoid a
# name clash with posterior.R's sigma2_max(mu, L, scale_min), which has a
# different signature.
# =============================================================================

# Bhatia-Davis (2000) variance ceiling. Supply hi, or L (then hi = lo + L - 1).
# Internal to the SD module (renamed from the reference's sigma2_max).
.sd_bound <- function(mu, lo = 1, hi = NULL, L = NULL) {
  if (is.null(hi)) {
    if (is.null(L)) stop("supply hi or L")
    hi <- lo + (L - 1)
  }
  out <- (mu - lo) * (hi - mu)
  out[!is.na(out) & out < 0] <- NA_real_
  out
}

#' Standardized mean and its sampling variance on a bounded scale
#'
#' Maps any bounded measure (whatever its `lo`/`hi`/`L`) to a standardized mean
#' `theta` in `(0, 1)` and its sampling variance. The common metric bridging
#' the SD module to the scale-length machinery.
#'
#' @param mean,sd,n Observed mean, SD, and sample size.
#' @param lo Lower endpoint (default 1).
#' @param hi Upper endpoint; if `NULL`, derived from `L`.
#' @param L Scale length (used when `hi` is `NULL`).
#' @return List with `theta` and `v_theta`.
#' @export
std_mean <- function(mean, sd, n, lo = 1, hi = NULL, L = NULL) {
  if (is.null(hi)) hi <- lo + (L - 1)
  span <- hi - lo
  list(theta = (mean - lo) / span,
       v_theta = (sd^2 / n) / span^2)
}

#' Mean inter-item correlation from coefficient alpha (inverted Spearman-Brown)
#'
#' `rbar = alpha / (k - alpha * (k - 1))`, clamped to `[-0.20, 0.95]`.
#'
#' @param alpha Coefficient alpha (vector).
#' @param k Item count (vector).
#' @return Vector of mean inter-item correlations.
#' @export
rbar_from_alpha_k <- function(alpha, k) {
  out <- rep(NA_real_, length(alpha))
  ok <- !is.na(alpha) & !is.na(k) & k >= 2 & alpha > 0 & alpha < 1
  denom <- k - alpha * (k - 1)
  ok <- ok & is.finite(denom) & denom > 0
  out[ok] <- alpha[ok] / denom[ok]
  pmin(pmax(out, -0.20), 0.95)
}

#' Impute mean inter-item correlation with a reliability cascade
#'
#' Fills the per-row mean inter-item correlation `rbar` using a reported ->
#' field -> global -> default cascade. The correction factor
#' `J / (1 + (J - 1) rbar)` is decreasing in `rbar`, so under-correction is the
#' damaging direction; uncertain `rbar` is shifted DOWN (toward more inflation).
#'
#' @param d Data.frame with `alpha`, `items`, `n`, and the `group` column.
#' @param group Column naming the field-level fallback stratum (default
#'   `"dimension"`).
#' @param default_rbar Default `rbar` when reliability is unreported (0.25).
#' @return `d` with added columns `alpha_filled`, `rbar_source`, `rbar`.
#' @export
impute_reliability <- function(d, group = "dimension", default_rbar = 0.25) {
  a <- d$alpha; grp <- as.character(d[[group]])
  ok <- !is.na(a) & a > 0 & a < 1 & !is.na(d$n)
  num <- tapply((a * d$n)[ok], grp[ok], sum)
  den <- tapply(d$n[ok],        grp[ok], sum)
  field <- num / den
  global <- if (any(ok)) sum((a * d$n)[ok]) / sum(d$n[ok]) else NA_real_

  af <- a; na_a <- is.na(af)
  af[na_a] <- field[grp[na_a]]
  field_filled <- na_a & !is.na(af)
  still <- is.na(af); af[still] <- global
  global_filled <- still & !is.na(af)
  d$alpha_filled <- af

  rbar_raw <- rbar_from_alpha_k(af, d$items)
  rbar_raw[is.na(rbar_raw)] <- default_rbar
  src <- ifelse(!na_a & !is.na(d$items), "reported",
         ifelse(field_filled & !is.na(d$items), "field",
         ifelse(global_filled & !is.na(d$items), "global", "default")))
  shift <- c(reported = 0.00, field = 0.05, global = 0.10, default = 0.00)
  d$rbar_source <- src
  d$rbar <- pmin(pmax(rbar_raw - shift[src], default_rbar), rbar_raw)
  d$rbar[src == "default"] <- default_rbar
  d
}

# Composite-variance correction factor J / (1 + (J - 1) rbar).
.composite_factor <- function(items, rbar) {
  J <- ifelse(is.na(items) | items < 1, 1, items)
  J / (1 + (J - 1) * rbar)
}

#' Composite-variance correction (forward / inverse / round-trip check)
#'
#' Inflate an item-level variance to the composite metric and back. Provided as
#' an explicit forward/inverse pair so the round-trip can be asserted.
#'
#' @param sd Observed SD.
#' @param s2c Corrected (composite) variance.
#' @param items Item count `J`.
#' @param rbar Mean inter-item correlation.
#' @param tol Round-trip tolerance.
#' @return `to_corrected_var`: corrected variance. `to_composite_sd`: SD.
#'   `assert_composite_roundtrip`: invisibly `TRUE` (errors otherwise).
#' @name composite_correction
#' @export
to_corrected_var <- function(sd, items, rbar) sd^2 * .composite_factor(items, rbar)

#' @rdname composite_correction
#' @export
to_composite_sd  <- function(s2c, items, rbar) sqrt(s2c / .composite_factor(items, rbar))

#' @rdname composite_correction
#' @export
assert_composite_roundtrip <- function(sd, items, rbar, tol = 1e-8) {
  s2c   <- to_corrected_var(sd, items, rbar)
  sd_rt <- to_composite_sd(s2c, items, rbar)
  stopifnot(all(is.na(sd) | abs(sd_rt - sd) < tol))
  invisible(TRUE)
}

#' Single-imputation feasibility screen and SD imputation
#'
#' Flags impossible SDs (above the Bhatia-Davis ceiling, 2 percent tolerance),
#' fits the log-phi dispersion model on the observed SDs, and single-imputes a
#' point prediction for each missing SD on a known-L anchor (the circularity
#' firewall). This is the version that produced the Study 3 recovery numbers;
#' it does NOT propagate imputation variance (use [sd_impute_mi()] for that).
#'
#' @param d Data.frame with `mean`, `sd` (NA where missing), `n`, `items`,
#'   `response_low`, and `response_high` (or `L`), plus optional `flag_scale`.
#' @param min_obs Minimum observed SDs required to fit the phi model (default 10).
#' @return `d` with `sd` filled where imputed, plus `sd_imputed` and
#'   `flag_scale`.
#' @export
sd_feasibility_and_impute <- function(d, min_obs = 10) {
  if (is.null(d$flag_scale)) d$flag_scale <- NA_character_
  lo <- d$response_low
  hi <- if (!is.null(d$response_high)) d$response_high else lo + (d$L - 1)
  bound <- (d$mean - lo) * (hi - d$mean)
  bound[!is.na(bound) & bound < 0] <- NA_real_

  bad_sd <- !is.na(d$sd) & !is.na(bound) & bound > 0 &
            (d$sd <= 0 | d$sd > sqrt(bound) * 1.02)
  d$flag_scale[is.na(d$flag_scale) & bad_sd] <- "sd_impossible"

  L_avail <- if (!is.null(d$L)) !is.na(d$L) else !is.na(hi)
  L_known <- L_avail & !is.na(d$response_low) & is.na(d$flag_scale)
  obs <- !is.na(d$sd) & d$sd > 0 & !is.na(bound) & bound > 0 & is.na(d$flag_scale)
  d$sd_imputed <- FALSE

  if (sum(obs) >= min_obs) {
    phi <- pmin(pmax(d$sd[obs] / sqrt(bound[obs]), 1e-3), 0.999)
    fitdat <- data.frame(logphi = log(phi),
                         mpos = (d$mean[obs] - lo[obs]) / (hi[obs] - lo[obs]),
                         logJ = log(pmax(d$items[obs], 1)))
    fit  <- tryCatch(lm(logphi ~ mpos + logJ, data = fitdat),
                     error = function(e) NULL)
    gmed <- median(phi)
    need <- which(is.na(d$sd) & !is.na(bound) & bound > 0 & L_known)
    for (i in need) {
      mpos_i <- (d$mean[i] - lo[i]) / (hi[i] - lo[i])
      logJ_i <- log(pmax(d$items[i], 1))
      ph <- if (!is.null(fit) && is.finite(logJ_i))
              exp(predict(fit, newdata = data.frame(mpos = mpos_i,
                                                    logJ = logJ_i))) else gmed
      d$sd[i] <- min(max(ph, 1e-3), 0.999) * sqrt(bound[i])
      d$sd_imputed[i] <- TRUE
    }
    message(sprintf("SD: imputed %d anchor SDs (phi model on %d observed).",
                    length(need), sum(obs)))
  } else {
    message("SD: too few observed SDs to fit the phi model; none imputed.")
  }
  d
}

# Design matrix for the phi model.
.phi_design <- function(mean, lo, hi, items)
  cbind(`(int)` = 1,
        mpos = (mean - lo) / (hi - lo),
        logJ = log(pmax(items, 1)))

# M x n_target matrix of predictive log-phi draws given observed (Xo, y) and
# target design Xt. Incorporates coefficient AND residual uncertainty (Rubin
# 1987 "norm" method).
.phi_mi_logdraws <- function(y, Xo, Xt, M) {
  Xo <- as.matrix(Xo); Xt <- as.matrix(Xt)
  keep <- is.finite(y) & is.finite(rowSums(Xo))
  Xo <- Xo[keep, , drop = FALSE]; y <- y[keep]
  n <- nrow(Xo); p <- ncol(Xo)
  out <- matrix(NA_real_, M, nrow(Xt))
  if (n > 0) {
    cm <- colMeans(Xo)
    for (j in seq_len(ncol(Xt))) {
      bad <- !is.finite(Xt[, j]); if (any(bad)) Xt[bad, j] <- cm[j]
    }
  }
  .fallback <- function() {
    df <- max(n - 1, 1); RSS <- sum((y - mean(y))^2); mu <- mean(y)
    for (m in seq_len(M)) {
      s2 <- RSS / rchisq(1, df)
      out[m, ] <<- (mu + sqrt(s2 / max(n, 1)) * rnorm(1)) +
        rnorm(nrow(Xt), 0, sqrt(s2))
    }
    out
  }
  if (n < p + 2) return(.fallback())
  XtX <- crossprod(Xo); XtX <- XtX + diag(1e-8 * mean(diag(XtX)) + 1e-12, p)
  R <- tryCatch(chol(XtX), error = function(e) NULL)
  if (is.null(R)) return(.fallback())
  bhat <- backsolve(R, backsolve(R, crossprod(Xo, y), transpose = TRUE))
  RSS  <- sum((y - Xo %*% bhat)^2); df <- max(n - p, 1)
  for (m in seq_len(M)) {
    s2   <- RSS / rchisq(1, df)
    beta <- bhat + sqrt(s2) * backsolve(R, rnorm(p))
    out[m, ] <- as.numeric(Xt %*% beta) + rnorm(nrow(Xt), 0, sqrt(s2))
  }
  out
}

#' Multiple imputation of missing SDs on bounded scales (Rubin)
#'
#' Multiply-imputes every missing SD on a known-L bounded scale by Bayesian
#' linear-regression draws of the dispersion fraction (coefficient AND residual
#' uncertainty), so the spread of the M imputed SDs reflects genuine predictive
#' uncertainty. Optionally conditions on a moderator `group` with a nested
#' fallback to the pooled fit (>= `min_group` observed SDs in the level).
#'
#' @param d Data.frame as for [sd_feasibility_and_impute()].
#' @param M Number of imputations (default 50).
#' @param min_obs Minimum observed SDs to fit the model (default 10).
#' @param group Optional moderator column for conditioned fits.
#' @param min_group Minimum observed SDs in a level to condition (default 10).
#' @param seed Optional seed.
#' @param return_draws If `TRUE`, also return the `M x n_need` SD draws.
#' @return List with `d` (filled `sd`, `sd_imputed`, `sd_impute_var` =
#'   between-imputation variance), `need`, and optionally `sd_draws`.
#' @export
sd_impute_mi <- function(d, M = 50, min_obs = 10, group = NULL,
                         min_group = 10, seed = NULL, return_draws = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  lo <- d$response_low
  hi <- if (!is.null(d$response_high)) d$response_high else lo + (d$L - 1)
  bound <- (d$mean - lo) * (hi - d$mean); bound[!is.na(bound) & bound < 0] <- NA_real_
  obs  <- which(!is.na(d$sd) & d$sd > 0 & !is.na(bound) & bound > 0)
  if (length(obs) < min_obs) stop("too few observed SDs to fit the phi model")
  L_avail <- if (!is.null(d$L)) !is.na(d$L) else !is.na(hi)
  need <- which(is.na(d$sd) & !is.na(bound) & bound > 0 & L_avail &
                  !is.na(d$response_low))
  y  <- log(pmin(pmax(d$sd[obs] / sqrt(bound[obs]), 1e-3), 0.999))
  Xo <- .phi_design(d$mean[obs],  lo[obs],  hi[obs],  d$items[obs])
  Xt <- .phi_design(d$mean[need], lo[need], hi[need], d$items[need])
  logdr <- .phi_mi_logdraws(y, Xo, Xt, M)
  if (!is.null(group) && length(need)) {
    go <- as.character(d[[group]])[obs]; gt <- as.character(d[[group]])[need]
    for (gv in unique(gt[!is.na(gt)])) {
      oi <- which(go == gv); ti <- which(gt == gv)
      if (length(oi) >= min_group && length(ti) > 0)
        logdr[, ti] <- .phi_mi_logdraws(y[oi], Xo[oi, , drop = FALSE],
                                        Xt[ti, , drop = FALSE], M)
    }
  }
  phi <- pmin(pmax(exp(logdr), 1e-3), 0.999)
  sd_draws <- sweep(phi, 2, sqrt(bound[need]), `*`)
  d$sd_imputed <- FALSE
  if (length(need)) {
    d$sd[need] <- apply(sd_draws, 2, median); d$sd_imputed[need] <- TRUE
    d$sd_impute_var <- NA_real_
    d$sd_impute_var[need] <- apply(sd_draws, 2, var)
  }
  res <- list(d = d, need = need)
  if (return_draws) res$sd_draws <- sd_draws
  res
}
