#' Random-effects meta-analytic fit (internal)
#'
#' Wraps `metafor::rma` with safe error handling. Requires the `metafor`
#' package to be installed (it's in Suggests).
#'
#' @param yi Vector of effects.
#' @param vi Vector of sampling variances.
#' @param method Heterogeneity estimator (default `"PM"`).
#' @return List with `b`, `se`, `ci_lo`, `ci_hi`, `tau2`, `k`.
#' @keywords internal
fit_re <- function(yi, vi, method = "PM") {
  if (length(yi) < 2) {
    return(list(b = NA, se = NA, ci_lo = NA, ci_hi = NA,
                tau2 = NA, k = length(yi)))
  }
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for meta-analytic pooling. ",
         "Install with install.packages('metafor').")
  }
  safe <- tryCatch(
    metafor::rma(yi = yi, vi = vi, method = method),
    error = function(e) NULL
  )
  if (is.null(safe)) {
    return(list(b = NA, se = NA, ci_lo = NA, ci_hi = NA,
                tau2 = NA, k = length(yi)))
  }
  list(
    b     = as.numeric(safe$b),
    se    = as.numeric(safe$se),
    ci_lo = as.numeric(safe$ci.lb),
    ci_hi = as.numeric(safe$ci.ub),
    tau2  = as.numeric(safe$tau2),
    k     = safe$k
  )
}

#' Run the full multi-method meta-analytic pool
#'
#' Implements four methods: Multiple Imputation with Rubin's rules,
#' Bayes-point (law of total variance), Listwise (observed-L only), and
#' Prior-only modal. Returns a tidy data.frame.
#'
#' @param study_results Per-study results from [run_study_imputation()].
#' @param studies Input data split by study_id.
#' @param base_prior Normalized prior.
#' @param M Number of multiple imputations.
#' @param tau2_method Heterogeneity estimator.
#' @return List with `pooled_df` (data.frame), and per-method named entries.
#' @keywords internal
run_meta_analysis <- function(study_results, studies, base_prior, M,
                              tau2_method = "PM") {
  # A study qualifies for listwise only if every block has observed L.
  is_fully_observed <- vapply(study_results, function(sr) {
    all(vapply(sr$block_data, function(b) !b$is_imputed, logical(1)))
  }, logical(1))
  obs_idx <- which(is_fully_observed)
  listwise <- list()
  for (k in obs_idx) {
    sr <- study_results[[k]]
    L_per_block <- vapply(sr$block_data, function(b) {
      if (b$is_imputed) b$L_modal else b$L_observed
    }, integer(1))
    tv <- study_theta_and_var_blocks(studies[[sr$study_id]],
                                     sr$block_id, L_per_block, sr$scale_min)
    listwise[[length(listwise) + 1]] <- c(theta = tv$theta, v = tv$v)
  }
  yi_lw <- vapply(listwise, function(x) x["theta"], numeric(1))
  vi_lw <- vapply(listwise, function(x) x["v"],     numeric(1))
  res_listwise <- fit_re(yi_lw, vi_lw, tau2_method)

  prior_modal_L <- as.integer(names(base_prior)[which.max(base_prior)])
  prior_only <- list()
  for (k in seq_along(study_results)) {
    sr <- study_results[[k]]
    L_per_block <- vapply(sr$block_data, function(b) {
      if (b$is_imputed) prior_modal_L else b$L_observed
    }, integer(1))
    tv <- study_theta_and_var_blocks(studies[[sr$study_id]],
                                     sr$block_id, L_per_block, sr$scale_min)
    prior_only[[length(prior_only) + 1]] <- c(theta = tv$theta, v = tv$v)
  }
  yi_p <- vapply(prior_only, function(x) x["theta"], numeric(1))
  vi_p <- vapply(prior_only, function(x) x["v"],     numeric(1))
  res_prior_only <- fit_re(yi_p, vi_p, tau2_method)

  # Bayes-point: n-weighted aggregate of per-block Bayes-point thetas.
  yi_b <- numeric(length(study_results))
  vi_b <- numeric(length(study_results))
  for (k in seq_along(study_results)) {
    sr <- study_results[[k]]
    n_b <- length(sr$block_data)
    block_thetas <- vapply(sr$block_data, function(b) b$theta_bayes,
                           numeric(1))
    block_vars   <- vapply(sr$block_data, function(b) b$var_total,
                           numeric(1))
    block_ns <- vapply(seq_len(n_b), function(b) {
      sum(studies[[sr$study_id]]$n[sr$block_id == b])
    }, numeric(1))
    total_n <- sum(block_ns)
    w <- if (total_n > 0) block_ns / total_n else rep(1 / n_b, n_b)
    yi_b[k] <- sum(w * block_thetas)
    vi_b[k] <- sum(w^2 * block_vars)
  }
  res_bayes <- fit_re(yi_b, vi_b, tau2_method)

  mi_fits <- vector("list", M)
  for (m in seq_len(M)) {
    yi_m <- numeric(length(study_results))
    vi_m <- numeric(length(study_results))
    for (k in seq_along(study_results)) {
      sr <- study_results[[k]]
      L_per_block <- vapply(sr$block_data, function(b) b$L_draws[m],
                            integer(1))
      tv <- study_theta_and_var_blocks(studies[[sr$study_id]],
                                       sr$block_id, L_per_block,
                                       sr$scale_min)
      yi_m[k] <- tv$theta
      vi_m[k] <- tv$v
    }
    mi_fits[[m]] <- fit_re(yi_m, vi_m, tau2_method)
  }

  bs   <- vapply(mi_fits, function(x) x$b,   numeric(1))
  ses  <- vapply(mi_fits, function(x) x$se,  numeric(1))
  tau2s <- vapply(mi_fits, function(x) x$tau2, numeric(1))
  ok <- is.finite(bs) & is.finite(ses)
  if (sum(ok) < 2) {
    res_mi <- list(b = NA, se = NA, ci_lo = NA, ci_hi = NA, tau2 = NA,
                   k = length(study_results), imputation_se_share = NA,
                   fmi = NA, df = NA)
  } else {
    bs_ok <- bs[ok]; ses_ok <- ses[ok]; tau2s_ok <- tau2s[ok]
    M_eff <- length(bs_ok)
    Wbar  <- mean(ses_ok^2)
    B     <- var(bs_ok)
    T_pool  <- Wbar + (1 + 1 / M_eff) * B
    se_pool <- sqrt(T_pool)
    b_pool  <- mean(bs_ok)
    imp_share <- sqrt((1 + 1 / M_eff) * B) / se_pool
    r  <- (1 + 1 / M_eff) * B / Wbar
    df <- (M_eff - 1) * (1 + 1 / r)^2
    fmi <- (r + 2 / (df + 3)) / (r + 1)
    crit <- qt(0.975, df = df)
    res_mi <- list(
      b     = b_pool,
      se    = se_pool,
      ci_lo = b_pool - crit * se_pool,
      ci_hi = b_pool + crit * se_pool,
      tau2  = mean(tau2s_ok),
      k     = length(study_results),
      imputation_se_share = imp_share,
      fmi   = fmi,
      df    = df
    )
  }

  pooled_df <- data.frame(
    method   = c("MI", "Bayes-point", "Listwise (observed only)", "Prior-only"),
    estimate = c(res_mi$b, res_bayes$b, res_listwise$b, res_prior_only$b),
    se       = c(res_mi$se, res_bayes$se, res_listwise$se, res_prior_only$se),
    ci_lo    = c(res_mi$ci_lo, res_bayes$ci_lo, res_listwise$ci_lo,
                 res_prior_only$ci_lo),
    ci_hi    = c(res_mi$ci_hi, res_bayes$ci_hi, res_listwise$ci_hi,
                 res_prior_only$ci_hi),
    tau2     = c(res_mi$tau2, res_bayes$tau2, res_listwise$tau2,
                 res_prior_only$tau2),
    k        = c(res_mi$k, res_bayes$k, res_listwise$k, res_prior_only$k),
    imputation_se_share = c(res_mi$imputation_se_share, NA, NA, NA),
    stringsAsFactors = FALSE
  )

  list(
    pooled_df  = pooled_df,
    mi         = res_mi,
    bayes      = res_bayes,
    listwise   = res_listwise,
    prior_only = res_prior_only,
    theta_bar  = res_mi$b,
    rubin_se   = res_mi$se,
    ci_lo      = res_mi$ci_lo,
    ci_hi      = res_mi$ci_hi,
    tau2       = res_mi$tau2,
    FMI        = res_mi$fmi,
    M          = M,
    K          = length(study_results)
  )
}

#' Run an MI pool on a user-supplied external effect column
#'
#' For each MI draw we pool the external effect column (which is invariant
#' across imputations) with its supplied SE. With invariant effects all M
#' fits are identical, so this reduces to a single random-effects fit, but
#' we wrap it through the same Rubin's-rules path for API consistency.
#'
#' @param studies Input data split by study_id.
#' @param study_results Per-study results.
#' @param effect_col Column name for the external effect.
#' @param effect_se_col Column name for the effect SE.
#' @param tau2_method Heterogeneity estimator.
#' @return Result list matching `run_meta_analysis()$mi`.
#' @keywords internal
run_external_effect_meta <- function(studies, study_results, effect_col,
                                     effect_se_col, tau2_method = "PM") {
  yi <- numeric(length(study_results))
  vi <- numeric(length(study_results))
  for (k in seq_along(study_results)) {
    sr <- study_results[[k]]
    meas <- studies[[sr$study_id]]
    if (!effect_col %in% names(meas) || !effect_se_col %in% names(meas)) {
      stop("effect_col / effect_se_col not present in study '", sr$study_id, "'.")
    }
    yi[k] <- meas[[effect_col]][1]
    vi[k] <- meas[[effect_se_col]][1]^2
  }
  fit_re(yi, vi, tau2_method)
}
