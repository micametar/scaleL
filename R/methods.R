#' Print a `scaleL` fit
#'
#' @param x A `scaleL` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.scaleL <- function(x, ...) {
  K <- length(x$studies)
  n_imp <- sum(vapply(x$studies, function(s) s$is_imputed, logical(1)))
  pct_imp <- 100 * n_imp / max(K, 1)
  modal_tab <- table(vapply(x$studies, function(s) s$L_modal, integer(1)))

  cat("scaleL fit\n")
  cat(sprintf("  Studies:       %d (%d imputed, %.1f%%)\n", K, n_imp, pct_imp))
  cat(sprintf("  Prior:         %s\n", x$prior_type))
  cat(sprintf("  L grid:        %s\n", paste(x$L_grid, collapse = ", ")))
  cat(sprintf("  M imputations: %d\n", x$M))
  cat("  Modal-L distribution:\n")
  for (nm in names(modal_tab)) {
    cat(sprintf("    L=%s: %d study(ies)\n", nm, as.integer(modal_tab[[nm]])))
  }
  if (!is.null(x$meta)) {
    cat(sprintf("  Pooled (MI):   theta = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",
                x$meta$theta_bar, x$meta$rubin_se,
                x$meta$ci_lo, x$meta$ci_hi))
    if (!is.null(x$meta$tau2) && is.finite(x$meta$tau2)) {
      cat(sprintf("                 tau^2 = %.4f, FMI = %.3f\n",
                  x$meta$tau2, x$meta$FMI))
    }
  }
  invisible(x)
}

#' Summarize a `scaleL` fit
#'
#' @param object A `scaleL` object.
#' @param ... Ignored.
#' @return A `summary.scaleL` object.
#' @export
summary.scaleL <- function(object, ...) {
  K <- length(object$studies)
  is_imp <- vapply(object$studies, function(s) s$is_imputed, logical(1))
  n_imp <- sum(is_imp); n_obs <- K - n_imp

  mean_cert <- if (n_imp > 0) {
    mean(vapply(object$studies,
                function(s) if (s$is_imputed) s$L_modal_prob else NA_real_,
                numeric(1)), na.rm = TRUE)
  } else NA_real_
  mean_ess_ratio <- if (n_imp > 0) {
    mean(vapply(object$studies, function(s) {
      if (s$is_imputed) s$ess / sum(object$data$n[object$data$study_id == s$study_id])
      else NA_real_
    }, numeric(1)), na.rm = TRUE)
  } else NA_real_

  obs_vals <- vapply(object$studies, function(s) as.integer(s$L_observed), integer(1))
  modal_vals <- vapply(object$studies, function(s) as.integer(s$L_modal), integer(1))
  obs_mask <- !is.na(obs_vals)
  recovery <- if (any(obs_mask)) mean(modal_vals[obs_mask] == obs_vals[obs_mask]) else NA_real_

  res <- list(
    K = K, n_imputed = n_imp, n_observed = n_obs,
    mean_modal_cert = mean_cert, mean_ess_ratio = mean_ess_ratio,
    observed_recovery = recovery,
    prior_type = object$prior_type,
    M = object$M,
    meta = object$meta,
    modal_table = table(modal_vals)
  )
  class(res) <- "summary.scaleL"
  res
}

#' @rdname summary.scaleL
#' @param x A `summary.scaleL` object.
#' @export
print.summary.scaleL <- function(x, ...) {
  cat("Bayesian imputation of Likert scale length: summary\n")
  cat("====================================================\n\n")
  cat(sprintf("Studies total:    %d\n", x$K))
  cat(sprintf("L observed:       %d\n", x$n_observed))
  cat(sprintf("L imputed:        %d\n", x$n_imputed))
  cat(sprintf("Prior used:       %s\n", x$prior_type))
  cat(sprintf("M imputations:    %d\n\n", x$M))

  cat("Imputation diagnostics (over imputed studies):\n")
  cat(sprintf("  Mean modal-L posterior probability:  %.3f\n", x$mean_modal_cert))
  cat(sprintf("  Mean ESS as fraction of raw n:       %.3f\n", x$mean_ess_ratio))
  if (!is.na(x$observed_recovery)) {
    cat(sprintf("\nRecovery on observed-L studies (built-in validation):\n"))
    cat(sprintf("  Modal-L recovery: %.1f%%\n", 100 * x$observed_recovery))
  }
  cat("\nModal-L distribution across studies:\n")
  for (nm in names(x$modal_table)) {
    cat(sprintf("  L=%s: %d\n", nm, as.integer(x$modal_table[[nm]])))
  }
  if (!is.null(x$meta)) {
    cat("\nPooled meta-analytic results:\n")
    print(x$meta$pooled_df, row.names = FALSE)
    if (!is.null(x$meta$FMI) && is.finite(x$meta$FMI)) {
      cat(sprintf("\nFraction of missing information (FMI): %.3f\n", x$meta$FMI))
    }
  }
  invisible(x)
}

#' Coerce a `scaleL` fit to a data.frame (the practical per-measure file)
#'
#' Returns the drop-in per-measure data.frame with `L_used`, `theta`,
#' `se_theta`, and `ess` columns, suitable for passing directly to
#' downstream meta-analytic tools such as `metafor::rma`.
#'
#' @param x A `scaleL` object.
#' @param row.names Unused.
#' @param optional Unused.
#' @param ... Unused.
#' @return The per-measure data.frame in `x$practical`.
#' @export
as.data.frame.scaleL <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$practical
}

#' Simple base-R plot for a `scaleL` fit
#'
#' Produces either a posterior heatmap (studies x L) or a modal-L bar chart.
#'
#' @param x A `scaleL` object.
#' @param type One of `"heatmap"` (default) or `"modal"`.
#' @param ... Passed to underlying plot functions.
#' @return Invisibly returns `NULL`.
#' @export
plot.scaleL <- function(x, type = c("heatmap", "modal"), ...) {
  type <- match.arg(type)
  if (type == "modal") {
    modal_vals <- vapply(x$studies, function(s) as.integer(s$L_modal), integer(1))
    barplot(table(modal_vals), xlab = "Modal L",
            ylab = "Number of studies",
            main = "Modal-L distribution", ...)
    return(invisible(NULL))
  }
  K <- length(x$studies)
  L_keys <- as.integer(names(x$studies[[1]]$post))
  M_mat <- matrix(0, nrow = K, ncol = length(L_keys))
  for (i in seq_len(K)) {
    M_mat[i, ] <- as.numeric(x$studies[[i]]$post)
  }
  image(x = seq_len(K), y = L_keys, z = M_mat,
        xlab = "Study index", ylab = "L",
        main = "Posterior over L (studies x L)",
        col = grDevices::hcl.colors(64, "YlOrRd", rev = TRUE), ...)
  invisible(NULL)
}
