#' Assign within-study uniformity blocks (internal)
#'
#' Given per-measure modal L values and modal probabilities, group eligible
#' measures (modal probability >= `min_cluster_modal_prob`) by their modal L.
#' Each group of at least `min_cluster_size` measures becomes a pooled block.
#' Remaining measures become singleton blocks (one block per measure).
#'
#' For backward compatibility, a study with a single measure always yields
#' exactly one block (a singleton).
#'
#' @param per_meas_modal_L Integer vector, one entry per measure.
#' @param per_meas_modal_prob Numeric vector of per-measure modal posterior
#'   probabilities (observed-L measures have probability 1).
#' @param min_cluster_size Minimum pooled-block size (default 2L).
#' @param min_cluster_modal_prob Minimum per-measure modal probability for
#'   pool eligibility (default 0.5).
#' @return Integer vector of block ids (1..n_blocks), same length as
#'   `per_meas_modal_L`.
#' @keywords internal
assign_blocks_within_study <- function(per_meas_modal_L,
                                       per_meas_modal_prob,
                                       min_cluster_size = 2L,
                                       min_cluster_modal_prob = 0.5) {
  n_meas <- length(per_meas_modal_L)
  block_id <- integer(n_meas)
  next_block <- 1L
  eligible <- per_meas_modal_prob >= min_cluster_modal_prob
  if (any(eligible)) {
    unique_modal_L <- unique(per_meas_modal_L[eligible])
    for (L_candidate in unique_modal_L) {
      members <- which(eligible & per_meas_modal_L == L_candidate)
      if (length(members) >= min_cluster_size) {
        block_id[members] <- next_block
        next_block <- next_block + 1L
      }
    }
  }
  singletons <- which(block_id == 0L)
  for (i in singletons) {
    block_id[i] <- next_block
    next_block <- next_block + 1L
  }
  block_id
}

# Internal: degenerate point-mass posterior on a single observed L.
.point_mass_posterior <- function(L_obs, L_grid) {
  post <- setNames(rep(0, length(L_grid)), as.character(L_grid))
  key <- as.character(as.integer(L_obs))
  if (key %in% names(post)) {
    post[key] <- 1
  } else {
    post <- c(post, setNames(1, key))
  }
  post
}

# Internal: compute per-measure posteriors (one per measure within a study).
.per_measure_posteriors <- function(measures, prior, L_grid, scale_min_s) {
  n_meas <- nrow(measures)
  out <- vector("list", n_meas)
  modal_L    <- integer(n_meas)
  modal_prob <- numeric(n_meas)
  observed   <- !is.na(measures$L)
  for (i in seq_len(n_meas)) {
    if (observed[i]) {
      L_obs_i <- as.integer(measures$L[i])
      pm <- .point_mass_posterior(L_obs_i, L_grid)
      out[[i]] <- pm
      modal_L[i] <- L_obs_i
      modal_prob[i] <- 1
    } else {
      pm <- posterior_L(measures[i, , drop = FALSE], prior, L_grid,
                        scale_min_s)
      out[[i]] <- pm
      idx <- which.max(pm)
      modal_L[i] <- as.integer(names(pm)[idx])
      modal_prob[i] <- as.numeric(pm[idx])
    }
  }
  list(posterior = out, modal_L = modal_L, modal_prob = modal_prob,
       observed = observed)
}

# Internal: per-block posterior, per-L theta/v, Bayes-point quantities, MI draws.
.compute_block <- function(block_measures, idx, prior, L_grid, scale_min_s,
                           M, s_name, b) {
  block_obs_L_vals <- unique(block_measures$L[!is.na(block_measures$L)])
  if (length(block_obs_L_vals) > 1) {
    stop("Internal logic error: block in study '", s_name,
         "' contains measures with inconsistent observed L.")
  }
  block_L_observed <- if (length(block_obs_L_vals) == 1) {
    as.integer(block_obs_L_vals)
  } else NA_integer_
  block_is_imputed <- is.na(block_L_observed)

  if (block_is_imputed) {
    block_post <- posterior_L(block_measures, prior, L_grid, scale_min_s)
  } else {
    block_post <- .point_mass_posterior(block_L_observed, L_grid)
  }

  L_keys_b <- as.integer(names(block_post))
  per_L_b <- lapply(L_keys_b, function(L) {
    study_theta_and_var(block_measures, L, scale_min_s)
  })
  theta_per_L_b <- vapply(per_L_b, function(x) x$theta, numeric(1))
  v_per_L_b     <- vapply(per_L_b, function(x) x$v, numeric(1))
  valid_b <- is.finite(theta_per_L_b) & is.finite(v_per_L_b) &
             block_post > 0
  if (!any(valid_b)) {
    stop("Study '", s_name, "', block ", b,
         ": no finite (theta, v) for any L with positive posterior mass.")
  }
  block_post_v <- block_post
  block_post_v[!valid_b] <- 0
  block_post_v <- block_post_v / sum(block_post_v)

  block_theta_bayes <- sum(block_post_v * theta_per_L_b)
  block_e_v   <- sum(block_post_v * v_per_L_b)
  block_var_t <- sum(block_post_v * (theta_per_L_b - block_theta_bayes)^2)
  block_var_total <- block_e_v + block_var_t
  block_se_pooled <- sqrt(block_var_total)
  block_L_modal <- L_keys_b[which.max(block_post_v)]
  block_L_bayes <- sum(block_post_v * L_keys_b)

  if (M > 0 && block_is_imputed) {
    block_L_draws <- sample(L_keys_b, size = M, replace = TRUE,
                            prob = block_post_v)
  } else {
    L_use_block <- if (block_is_imputed) block_L_modal else block_L_observed
    block_L_draws <- rep(L_use_block, M)
  }

  list(
    idx              = idx,
    size             = length(idx),
    is_imputed       = block_is_imputed,
    L_observed       = block_L_observed,
    post             = block_post_v,
    L_keys           = L_keys_b,
    theta_per_L      = theta_per_L_b,
    v_per_L          = v_per_L_b,
    L_modal          = block_L_modal,
    L_modal_prob     = max(block_post_v),
    L_bayes          = block_L_bayes,
    theta_bayes      = block_theta_bayes,
    se_with_imp      = block_se_pooled,
    var_total        = block_var_total,
    L_draws          = block_L_draws
  )
}

#' Per-study imputation loop (internal)
#'
#' For each study, computes per-measure posteriors over L, groups measures
#' into uniformity blocks (clusters of measures whose modal L agrees), pools
#' each block jointly, and derives per-measure theta/SE/ESS from the
#' measure's block. Backward-compatible with the original single-block
#' (uniform) behaviour when a study has only one measure or all measures
#' agree.
#'
#' @param d Validated input data.frame.
#' @param base_prior Normalized prior over `L_grid`.
#' @param L_grid Integer vector of candidate L values.
#' @param M Number of multiple imputations.
#' @param instrument_priors Optional named list of Tier 3 instrument priors.
#' @param min_cluster_size Minimum size of a pooled block (default 2L).
#' @param min_cluster_modal_prob Minimum per-measure modal probability
#'   required to be eligible for pooling (default 0.5).
#' @return List of per-study result records.
#' @keywords internal
run_study_imputation <- function(d, base_prior, L_grid, M,
                                 instrument_priors = NULL,
                                 min_cluster_size = 2L,
                                 min_cluster_modal_prob = 0.5) {
  studies <- split(d, d$study_id)
  study_results <- vector("list", length(studies))
  names(study_results) <- names(studies)

  for (s_name in names(studies)) {
    measures <- studies[[s_name]]
    scale_min_s <- measures$scale_min[1]
    n_meas <- nrow(measures)

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

    # ---- Per-measure pre-processing ----
    pm <- .per_measure_posteriors(measures, prior, L_grid, scale_min_s)
    per_meas_modal_L    <- pm$modal_L
    per_meas_modal_prob <- pm$modal_prob
    per_meas_observed   <- pm$observed

    # ---- Block assignment via uniformity clustering ----
    block_id <- assign_blocks_within_study(per_meas_modal_L,
                                           per_meas_modal_prob,
                                           min_cluster_size,
                                           min_cluster_modal_prob)
    n_blocks <- max(block_id)

    # ---- Per-block processing ----
    block_data <- vector("list", n_blocks)
    for (b in seq_len(n_blocks)) {
      idx <- which(block_id == b)
      block_data[[b]] <- .compute_block(measures[idx, , drop = FALSE],
                                        idx, prior, L_grid, scale_min_s,
                                        M, s_name, b)
    }

    # ---- Per-measure outputs (derived from each measure's block) ----
    per_measure_theta        <- numeric(n_meas)
    per_measure_se           <- numeric(n_meas)
    per_measure_ess          <- numeric(n_meas)
    per_measure_L_used       <- integer(n_meas)
    per_measure_block_id     <- block_id
    per_measure_block_size   <- integer(n_meas)
    per_measure_pooled_flag  <- logical(n_meas)
    per_measure_L_modal_prob <- numeric(n_meas)

    for (i in seq_len(n_meas)) {
      b <- block_id[i]
      blk <- block_data[[b]]
      L_keys_b <- blk$L_keys
      post_b   <- blk$post

      L_used_i <- blk$L_modal
      per_measure_L_used[i]       <- L_used_i
      per_measure_block_size[i]   <- blk$size
      per_measure_pooled_flag[i]  <- (blk$size >= 2L && blk$is_imputed)
      per_measure_L_modal_prob[i] <- blk$L_modal_prob

      if (!blk$is_imputed) {
        L_obs <- blk$L_observed
        per_measure_theta[i] <- (measures$mean[i] - scale_min_s) / (L_obs - 1)
        per_measure_se[i]    <- (measures$sd[i] / (L_obs - 1)) /
                                  sqrt(measures$n[i])
        per_measure_ess[i]   <- measures$n[i]
      } else {
        theta_i_per_L <- (measures$mean[i] - scale_min_s) / (L_keys_b - 1)
        v_i_per_L     <- ((measures$sd[i] / (L_keys_b - 1))^2) /
                          measures$n[i]
        valid_i <- is.finite(theta_i_per_L) & is.finite(v_i_per_L) &
                   post_b > 0
        pw <- post_b
        pw[!valid_i] <- 0
        if (sum(pw) > 0) {
          pw <- pw / sum(pw)
        } else {
          pw <- post_b
        }
        theta_i <- sum(pw * theta_i_per_L)
        e_v_i   <- sum(pw * v_i_per_L)
        var_t_i <- sum(pw * (theta_i_per_L - theta_i)^2)
        v_i_tot <- e_v_i + var_t_i
        per_measure_theta[i] <- theta_i
        per_measure_se[i]    <- sqrt(v_i_tot)
        v_oracle_i <- v_i_per_L[L_keys_b == L_used_i]
        if (length(v_oracle_i) == 1 && is.finite(v_oracle_i) &&
            v_oracle_i > 0 && is.finite(v_i_tot) && v_i_tot > 0) {
          per_measure_ess[i] <- measures$n[i] * (v_oracle_i / v_i_tot)
        } else {
          per_measure_ess[i] <- NA_real_
        }
      }
    }

    # ---- Study-level data-driven posterior (diagnostic only) ----
    data_driven_post <- posterior_L(measures, prior, L_grid, scale_min_s)
    dd_L_keys <- as.integer(names(data_driven_post))
    dd_modal_idx <- which.max(data_driven_post)
    dd_L_modal      <- dd_L_keys[dd_modal_idx]
    dd_L_modal_prob <- as.numeric(data_driven_post[dd_modal_idx])
    dd_L_bayes      <- sum(data_driven_post * dd_L_keys)

    # ---- Study-level summaries (dominant block) ----
    n_imputed_blocks   <- sum(vapply(block_data,
                                     function(x) x$is_imputed, logical(1)))
    n_clustered_blocks <- sum(vapply(block_data,
                                     function(x) x$size >= 2L && x$is_imputed,
                                     logical(1)))
    study_is_imputed    <- n_imputed_blocks > 0
    study_is_nonuniform <- n_blocks > 1
    block_sizes <- vapply(block_data, function(x) x$size, integer(1))
    dom_b   <- which.max(block_sizes)
    dom_blk <- block_data[[dom_b]]

    if (!study_is_imputed && !study_is_nonuniform) {
      L_observed <- dom_blk$L_observed
      is_imputed <- FALSE
    } else {
      L_observed <- NA_integer_
      is_imputed <- TRUE
    }

    post_v        <- dom_blk$post
    L_keys        <- dom_blk$L_keys
    theta_per_L   <- dom_blk$theta_per_L
    v_per_L       <- dom_blk$v_per_L
    theta_bayes_recovered <- dom_blk$theta_bayes
    var_total_recovered   <- dom_blk$var_total
    se_recovered          <- sqrt(var_total_recovered)
    L_modal       <- dom_blk$L_modal
    L_bayes       <- dom_blk$L_bayes
    v_oracle_dom  <- dom_blk$v_per_L[dom_blk$L_keys == L_modal]
    if (length(v_oracle_dom) == 0 || !is.finite(v_oracle_dom) ||
        v_oracle_dom <= 0) {
      ess_recovered <- NA_real_
    } else {
      total_n_dom <- sum(measures$n[block_id == dom_b])
      ess_recovered <- total_n_dom * (v_oracle_dom / var_total_recovered)
    }
    theta_bayes  <- dom_blk$theta_bayes
    se_with_imp  <- dom_blk$se_with_imp
    ess_val      <- ess_recovered
    L_used_study <- dom_blk$L_modal
    L_draws      <- dom_blk$L_draws

    true_L_s <- if ("true_L" %in% names(measures)) {
      tl_vals <- unique(measures$true_L[!is.na(measures$true_L)])
      if (length(tl_vals) > 1) {
        warning("Study '", s_name, "': inconsistent true_L; using first value.")
      }
      if (length(tl_vals) >= 1) as.integer(tl_vals[1]) else NA_integer_
    } else NA_integer_

    study_results[[s_name]] <- list(
      study_id           = s_name,
      n_measures         = n_meas,
      L_observed         = L_observed,
      L_used             = L_used_study,
      true_L             = true_L_s,
      is_imputed         = is_imputed,
      is_nonuniform      = study_is_nonuniform,
      n_blocks           = n_blocks,
      n_clustered_blocks = n_clustered_blocks,
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
      data_driven_post      = data_driven_post,
      dd_L_modal            = dd_L_modal,
      dd_L_modal_prob       = dd_L_modal_prob,
      dd_L_bayes            = dd_L_bayes,
      per_measure_theta        = per_measure_theta,
      per_measure_se           = per_measure_se,
      per_measure_ess          = per_measure_ess,
      per_measure_L_used       = per_measure_L_used,
      per_measure_block_id     = per_measure_block_id,
      per_measure_block_size   = per_measure_block_size,
      per_measure_pooled_flag  = per_measure_pooled_flag,
      per_measure_L_modal_prob = per_measure_L_modal_prob,
      per_meas_observed        = per_meas_observed,
      block_data         = block_data,
      block_id           = block_id,
      theta_per_L        = setNames(theta_per_L, as.character(L_keys)),
      v_per_L            = setNames(v_per_L, as.character(L_keys)),
      scale_min          = scale_min_s
    )
  }

  study_results
}

# Internal: aggregate per-block (theta, v) into a single study-level pair
# using n-weighted averaging across blocks (independent within a study).
study_theta_and_var_blocks <- function(meas, block_id, L_per_block,
                                       scale_min) {
  n_b <- max(block_id)
  if (length(L_per_block) != n_b) {
    stop("L_per_block length mismatch in study_theta_and_var_blocks")
  }
  block_thetas <- numeric(n_b)
  block_vars   <- numeric(n_b)
  block_ns     <- numeric(n_b)
  for (b in seq_len(n_b)) {
    idx <- which(block_id == b)
    block_meas <- meas[idx, , drop = FALSE]
    tv <- study_theta_and_var(block_meas, L_per_block[b], scale_min)
    block_thetas[b] <- tv$theta
    block_vars[b]   <- tv$v
    block_ns[b]     <- sum(block_meas$n)
  }
  if (any(!is.finite(block_thetas)) || any(!is.finite(block_vars))) {
    return(list(theta = NA_real_, v = NA_real_))
  }
  total_n <- sum(block_ns)
  if (total_n <= 0) return(list(theta = NA_real_, v = NA_real_))
  w <- block_ns / total_n
  list(theta = sum(w * block_thetas), v = sum(w^2 * block_vars))
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

    L_used_vec      <- sr$per_measure_L_used
    L_source_vec    <- ifelse(sr$per_meas_observed, "observed", "imputed")
    L_certainty_vec <- ifelse(sr$per_meas_observed, 1,
                              sr$per_measure_L_modal_prob)

    practical_rows[[k]] <- data.frame(
      study_id       = meas$study_id,
      measure_id     = meas$measure_id,
      mean           = meas$mean,
      sd             = meas$sd,
      n              = meas$n,
      scale_min      = meas$scale_min,
      L_original     = meas$L,
      L_used         = L_used_vec,
      L_source       = L_source_vec,
      L_certainty    = L_certainty_vec,
      theta          = sr$per_measure_theta,
      se_theta       = sr$per_measure_se,
      ess            = sr$per_measure_ess,
      block_id       = sr$per_measure_block_id,
      block_size     = sr$per_measure_block_size,
      pooled_flag    = sr$per_measure_pooled_flag,
      instrument     = meas$instrument,
      stringsAsFactors = FALSE
    )

    dd_post <- sr$data_driven_post
    post_cols <- as.list(dd_post)
    names(post_cols) <- paste0("L_posterior_", names(dd_post))
    mi_cols <- as.list(sr$L_draws)
    names(mi_cols) <- paste0("L_imputed_", seq_len(M))

    # Per-block modal-L probabilities (named "block<b>_modal_prob")
    block_mp <- vapply(sr$block_data, function(b) b$L_modal_prob,
                       numeric(1))
    block_modal_cols <- as.list(block_mp)
    names(block_modal_cols) <- paste0("block", seq_along(block_mp),
                                      "_modal_prob")

    diag_base <- data.frame(
      study_id                = sr$study_id,
      n_measures              = sr$n_measures,
      n_blocks                = sr$n_blocks,
      n_clustered_blocks      = sr$n_clustered_blocks,
      study_nonuniform        = sr$is_nonuniform,
      L_observed              = sr$L_observed,
      L_was_imputed           = sr$is_imputed,
      L_modal_recovered       = sr$dd_L_modal,
      L_modal_prob            = sr$dd_L_modal_prob,
      L_bayes_recovered       = sr$dd_L_bayes,
      modal_recovers_observed = if (!is.na(sr$L_observed)) {
        sr$dd_L_modal == sr$L_observed
      } else NA,
      posterior_on_observed_L = if (!is.na(sr$L_observed)) {
        key <- as.character(sr$L_observed)
        if (key %in% names(dd_post)) as.numeric(dd_post[[key]]) else 0
      } else NA_real_,
      theta_bayes_recovered   = sr$theta_bayes_recovered,
      se_recovered            = sr$se_recovered,
      ess_recovered           = sr$ess_recovered,
      true_L                  = sr$true_L,
      modal_recovers_truth    = if (!is.na(sr$true_L)) {
        sr$dd_L_modal == sr$true_L
      } else NA,
      stringsAsFactors        = FALSE
    )
    diagnostic_rows[[k]] <- cbind(
      diag_base,
      as.data.frame(post_cols, check.names = FALSE),
      as.data.frame(mi_cols,   check.names = FALSE),
      as.data.frame(block_modal_cols, check.names = FALSE)
    )
  }

  # rbind tolerant of variable block-column counts across studies
  pad_cols <- function(rows) {
    all_cols <- unique(unlist(lapply(rows, names)))
    lapply(rows, function(df) {
      missing <- setdiff(all_cols, names(df))
      for (m in missing) df[[m]] <- NA
      df[, all_cols, drop = FALSE]
    })
  }

  list(
    practical  = do.call(rbind, practical_rows),
    diagnostic = do.call(rbind, pad_cols(diagnostic_rows))
  )
}
