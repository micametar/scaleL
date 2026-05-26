#' Per-study imputation loop (internal)
#'
#' Lifted verbatim from the practitioner script: computes the data-driven
#' posterior over L, per-measure theta/SE/ESS, and M multiple-imputation
#' draws for every study in the input.
#'
#' @param d Validated input data.frame.
#' @param base_prior Normalized prior over `L_grid`.
#' @param L_grid Integer vector of candidate L values.
#' @param M Number of multiple imputations.
#' @param instrument_priors Optional named list of Tier 3 instrument priors.
#' @return List of per-study result records.
#' @keywords internal
run_study_imputation <- function(d, base_prior, L_grid, M,
                                 instrument_priors = NULL) {
  studies <- split(d, d$study_id)
  study_results <- vector("list", length(studies))
  names(study_results) <- names(studies)

  for (s_name in names(studies)) {
    measures <- studies[[s_name]]
    scale_min_s <- measures$scale_min[1]

    L_vals <- unique(measures$L[!is.na(measures$L)])
    if (length(L_vals) > 1) {
      stop("Study '", s_name, "' has inconsistent L across measures: ",
           paste(L_vals, collapse = ", "),
           ". Each study must share a single L.")
    }
    L_observed <- if (length(L_vals) == 1) L_vals else NA_integer_
    is_imputed <- is.na(L_observed)

    prior <- base_prior
    if (!is.null(instrument_priors) && !is.na(measures$instrument[1]) &&
        measures$instrument[1] %in% names(instrument_priors)) {
      inst_prior <- instrument_priors[[measures$instrument[1]]]
      L_char <- as.character(L_grid)
      prior <- setNames(rep(0, length(L_char)), L_char)
      prior[names(inst_prior)[names(inst_prior) %in% L_char]] <-
        inst_prior[names(inst_prior)[names(inst_prior) %in% L_char]]
      prior <- prior / sum(prior)
    }

    post <- posterior_L(measures, prior, L_grid, scale_min_s)

    L_keys <- as.integer(names(post))
    per_L <- lapply(L_keys, function(L) {
      study_theta_and_var(measures, L, scale_min_s)
    })
    theta_per_L <- vapply(per_L, function(x) x$theta, numeric(1))
    v_per_L     <- vapply(per_L, function(x) x$v, numeric(1))
    valid <- is.finite(theta_per_L) & is.finite(v_per_L) & post > 0
    if (!any(valid)) {
      stop("Study '", s_name,
           "': no finite (theta, v) for any L with positive posterior mass.")
    }
    post_v <- post
    post_v[!valid] <- 0
    post_v <- post_v / sum(post_v)

    theta_bayes_recovered <- sum(post_v * theta_per_L)
    e_v_given_L           <- sum(post_v * v_per_L)
    var_theta_given_L     <- sum(post_v * (theta_per_L - theta_bayes_recovered)^2)
    var_total_recovered   <- e_v_given_L + var_theta_given_L
    se_recovered          <- sqrt(var_total_recovered)

    L_modal <- L_keys[which.max(post_v)]
    L_bayes <- sum(post_v * L_keys)

    v_oracle <- v_per_L[L_keys == L_modal]
    if (length(v_oracle) == 0 || !is.finite(v_oracle) || v_oracle <= 0) {
      ess_recovered <- NA_real_
    } else {
      ess_recovered <- sum(measures$n) * (v_oracle / var_total_recovered)
    }

    if (is_imputed) {
      theta_bayes <- theta_bayes_recovered
      se_with_imp <- se_recovered
      ess_val     <- ess_recovered
    } else {
      tv_obs <- study_theta_and_var(measures, L_observed, scale_min_s)
      theta_bayes <- tv_obs$theta
      se_with_imp <- sqrt(tv_obs$v)
      ess_val     <- sum(measures$n)
    }

    if (is_imputed) {
      L_used <- L_modal
      per_measure_theta <- numeric(nrow(measures))
      per_measure_se    <- numeric(nrow(measures))
      per_measure_ess   <- numeric(nrow(measures))
      for (i in seq_len(nrow(measures))) {
        theta_i_per_L <- (measures$mean[i] - scale_min_s) / (L_keys - 1)
        v_i_per_L     <- ((measures$sd[i] / (L_keys - 1))^2) / measures$n[i]
        valid_i <- is.finite(theta_i_per_L) & is.finite(v_i_per_L) & post > 0
        pw <- post_v
        pw[!valid_i] <- 0
        if (sum(pw) > 0) {
          pw <- pw / sum(pw)
        } else {
          pw <- post_v
        }
        theta_i <- sum(pw * theta_i_per_L)
        e_v_i   <- sum(pw * v_i_per_L)
        var_t_i <- sum(pw * (theta_i_per_L - theta_i)^2)
        v_i_tot <- e_v_i + var_t_i
        per_measure_theta[i] <- theta_i
        per_measure_se[i]    <- sqrt(v_i_tot)
        v_oracle_i <- v_i_per_L[L_keys == L_modal]
        if (length(v_oracle_i) == 1 && is.finite(v_oracle_i) &&
            v_oracle_i > 0 && is.finite(v_i_tot) && v_i_tot > 0) {
          per_measure_ess[i] <- measures$n[i] * (v_oracle_i / v_i_tot)
        } else {
          per_measure_ess[i] <- NA_real_
        }
      }
    } else {
      L_used <- L_observed
      per_measure_theta <- (measures$mean - scale_min_s) / (L_observed - 1)
      per_measure_se    <- (measures$sd / (L_observed - 1)) /
                             sqrt(measures$n)
      per_measure_ess   <- measures$n
    }

    if (M > 0 && is_imputed) {
      L_draws <- sample(L_keys, size = M, replace = TRUE, prob = post_v)
    } else {
      L_draws <- rep(if (is_imputed) L_modal else L_observed, M)
    }

    true_L_s <- if ("true_L" %in% names(measures)) {
      tl_vals <- unique(measures$true_L[!is.na(measures$true_L)])
      if (length(tl_vals) > 1) {
        warning("Study '", s_name, "': inconsistent true_L; using first value.")
      }
      if (length(tl_vals) >= 1) as.integer(tl_vals[1]) else NA_integer_
    } else NA_integer_

    study_results[[s_name]] <- list(
      study_id           = s_name,
      n_measures         = nrow(measures),
      L_observed         = L_observed,
      L_used             = L_used,
      true_L             = true_L_s,
      is_imputed         = is_imputed,
      post               = post_v,
      L_modal            = L_modal,
      L_modal_prob       = max(post_v),
      L_bayes            = L_bayes,
      L_draws            = L_draws,
      theta_bayes        = theta_bayes,
      se_with_imputation = se_with_imp,
      ess                = ess_val,
      theta_bayes_recovered = theta_bayes_recovered,
      se_recovered          = se_recovered,
      ess_recovered         = ess_recovered,
      per_measure_theta  = per_measure_theta,
      per_measure_se     = per_measure_se,
      per_measure_ess    = per_measure_ess,
      theta_per_L        = setNames(theta_per_L, as.character(L_keys)),
      v_per_L            = setNames(v_per_L, as.character(L_keys)),
      scale_min          = scale_min_s
    )
  }

  study_results
}

#' Build the practical (per-measure) and diagnostic (per-study) data.frames
#'
#' @param study_results Output of [run_study_imputation()].
#' @param studies Splits of the input by study_id.
#' @param L_grid Candidate L grid.
#' @param M Number of multiple imputations.
#' @return Named list with `practical` and `diagnostic` data.frames.
#' @keywords internal
build_output_frames <- function(study_results, studies, L_grid, M) {
  practical_rows <- vector("list", length(study_results))
  diagnostic_rows <- vector("list", length(study_results))

  for (k in seq_along(study_results)) {
    sr <- study_results[[k]]
    meas <- studies[[sr$study_id]]

    practical_rows[[k]] <- data.frame(
      study_id       = meas$study_id,
      measure_id     = meas$measure_id,
      mean           = meas$mean,
      sd             = meas$sd,
      n              = meas$n,
      scale_min      = meas$scale_min,
      L_original     = meas$L,
      L_used         = sr$L_used,
      L_source       = if (sr$is_imputed) "imputed" else "observed",
      L_certainty    = if (sr$is_imputed) sr$L_modal_prob else 1,
      theta          = sr$per_measure_theta,
      se_theta       = sr$per_measure_se,
      ess            = sr$per_measure_ess,
      instrument     = meas$instrument,
      stringsAsFactors = FALSE
    )

    post_cols <- as.list(sr$post)
    names(post_cols) <- paste0("L_posterior_", names(sr$post))
    mi_cols <- as.list(sr$L_draws)
    names(mi_cols) <- paste0("L_imputed_", seq_len(M))

    diag_base <- data.frame(
      study_id                = sr$study_id,
      n_measures              = sr$n_measures,
      L_observed              = sr$L_observed,
      L_was_imputed           = sr$is_imputed,
      L_modal_recovered       = sr$L_modal,
      L_modal_prob            = sr$L_modal_prob,
      L_bayes_recovered       = sr$L_bayes,
      modal_recovers_observed = if (!is.na(sr$L_observed)) {
        sr$L_modal == sr$L_observed
      } else NA,
      posterior_on_observed_L = if (!is.na(sr$L_observed)) {
        key <- as.character(sr$L_observed)
        if (key %in% names(sr$post)) as.numeric(sr$post[[key]]) else 0
      } else NA_real_,
      theta_bayes_recovered   = sr$theta_bayes_recovered,
      se_recovered            = sr$se_recovered,
      ess_recovered           = sr$ess_recovered,
      true_L                  = sr$true_L,
      modal_recovers_truth    = if (!is.na(sr$true_L)) {
        sr$L_modal == sr$true_L
      } else NA,
      stringsAsFactors        = FALSE
    )
    diagnostic_rows[[k]] <- cbind(diag_base,
                                  as.data.frame(post_cols, check.names = FALSE),
                                  as.data.frame(mi_cols,   check.names = FALSE))
  }

  list(
    practical  = do.call(rbind, practical_rows),
    diagnostic = do.call(rbind, diagnostic_rows)
  )
}
