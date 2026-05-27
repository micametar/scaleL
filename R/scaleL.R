#' Bayesian imputation of Likert scale length for meta-analysis
#'
#' Main user-facing entry point. Takes a per-measure data.frame, computes
#' the posterior over candidate scale lengths L for each study, produces
#' multiple imputations and per-measure standardized means with calibrated
#' SEs, and (optionally) pools the result with random-effects meta-analysis
#' using Rubin's rules.
#'
#' @param data Data.frame with one row per measure. Required columns:
#'   `study_id`, `mean`, `sd`, `n`. Optional: `L` (NA where missing),
#'   `scale_min` (default 1), `measure_id`, `instrument`, `true_L`.
#' @param prior Prior specification: `"tier1"`, `"tier1_shifted"`,
#'   `"field:<name>"` (where `<name>` is one of the entries in
#'   [FIELD_PRIORS]), or `"custom"`. Default `"tier1"`.
#' @param custom_prior Named numeric vector when `prior = "custom"`.
#' @param instrument_priors Optional named list of Tier 3 instrument-specific
#'   priors. Studies whose `instrument` matches a name in this list use the
#'   matching prior; others fall back to `prior`.
#' @param L_grid Integer vector of candidate scale lengths. Default `2:11`.
#' @param M Number of multiple imputations. Default 50.
#' @param method Reserved for future use. Currently the returned object
#'   contains MI, modal, and bayes-point quantities; the method argument
#'   only controls which is highlighted by `print()` (default `"MI"`).
#' @param compute_meta Logical; if `TRUE`, also fits a random-effects
#'   meta-analysis (requires the `metafor` package).
#' @param effect_col,effect_se_col Optional column names: when supplied,
#'   meta-analytic pooling is done on this external effect column instead of
#'   the standardized mean. Imputation of L is still performed for use
#'   downstream.
#' @param tau2_method Heterogeneity estimator passed to `metafor::rma`.
#'   Default `"PM"` (Paule-Mandel).
#' @param min_cluster_size Minimum number of measures whose per-measure
#'   modal-L posteriors agree before they are pooled into a joint block.
#'   Default `2L`.
#' @param min_cluster_modal_prob Minimum per-measure modal posterior
#'   probability required for a measure to be eligible for pooling with
#'   other measures sharing its modal L. Default `0.5`.
#' @param seed Optional integer seed.
#' @param verbose Logical; emit progress messages.
#' @section Within-study L clustering:
#' Each measure's per-measure posterior over L is computed independently.
#' Measures whose modal L agrees (with modal probability at least
#' `min_cluster_modal_prob`) are grouped; groups of at least
#' `min_cluster_size` measures form a pooled block and are imputed
#' jointly. Remaining measures are imputed as singletons. The per-measure
#' output `pooled_flag` indicates which measures were pooled.
#' @return An S3 object of class `"scaleL"` with elements `data`, `studies`,
#'   `practical`, `diagnostic`, `meta` (if `compute_meta = TRUE`), `call`,
#'   `prior_used`, `L_grid`, `M`, `method`.
#' @examples
#' d <- read.csv(system.file("extdata", "sample_data.csv", package = "scaleL"))
#' fit <- scaleL(d, prior = "tier1", M = 10, compute_meta = FALSE)
#' print(fit)
#' @export
scaleL <- function(data,
                   prior = "tier1",
                   custom_prior = NULL,
                   instrument_priors = NULL,
                   L_grid = 2:11,
                   M = 50L,
                   method = c("MI", "modal", "bayes_point"),
                   compute_meta = FALSE,
                   effect_col = NULL,
                   effect_se_col = NULL,
                   tau2_method = "PM",
                   min_cluster_size = 2L,
                   min_cluster_modal_prob = 0.5,
                   seed = NULL,
                   verbose = FALSE) {
  call <- match.call()
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  d <- validate_input(as.data.frame(data), verbose = verbose)

  base_prior <- resolve_prior(prior, custom_prior, L_grid)
  if (isTRUE(verbose)) {
    message(sprintf("Using prior '%s' over L in {%s}.",
                    prior, paste(L_grid, collapse = ", ")))
  }

  studies <- split(d, d$study_id)
  if (isTRUE(verbose)) {
    message(sprintf("Computing posteriors for %d studies ...", length(studies)))
  }
  study_results <- run_study_imputation(d, base_prior, L_grid, M,
                                        instrument_priors,
                                        min_cluster_size = min_cluster_size,
                                        min_cluster_modal_prob =
                                          min_cluster_modal_prob)
  frames <- build_output_frames(study_results, studies, L_grid, M)

  meta <- NULL
  if (isTRUE(compute_meta)) {
    if (isTRUE(verbose)) message("Fitting meta-analysis ...")
    if (!is.null(effect_col)) {
      ext <- run_external_effect_meta(studies, study_results, effect_col,
                                      effect_se_col, tau2_method)
      meta <- list(theta_bar = ext$b, rubin_se = ext$se,
                   ci_lo = ext$ci_lo, ci_hi = ext$ci_hi,
                   tau2 = ext$tau2, FMI = NA_real_, M = M,
                   K = length(study_results),
                   pooled_df = data.frame(
                     method = "External effect (random-effects)",
                     estimate = ext$b, se = ext$se,
                     ci_lo = ext$ci_lo, ci_hi = ext$ci_hi,
                     tau2 = ext$tau2, k = ext$k,
                     imputation_se_share = NA,
                     stringsAsFactors = FALSE
                   ),
                   external = TRUE)
    } else {
      meta <- run_meta_analysis(study_results, studies, base_prior, M,
                                tau2_method)
      meta$external <- FALSE
    }
  }

  out <- list(
    data        = d,
    studies     = study_results,
    practical   = frames$practical,
    diagnostic  = frames$diagnostic,
    meta        = meta,
    call        = call,
    prior_used  = base_prior,
    prior_type  = prior,
    L_grid      = L_grid,
    M           = M,
    method      = method
  )
  class(out) <- "scaleL"
  out
}

#' Convenience wrapper: run [scaleL()] from a CSV path
#'
#' @param path Path to a CSV file in the input format documented in [scaleL()].
#' @param ... Further arguments passed to [scaleL()].
#' @return A `scaleL` object.
#' @examples
#' \dontrun{
#' fit <- scaleL_from_csv("my_data.csv", prior = "tier1")
#' }
#' @export
scaleL_from_csv <- function(path, ...) {
  d <- as.data.frame(readr::read_csv(path, show_col_types = FALSE))
  scaleL(d, ...)
}

#' Write a `scaleL` fit to disk as CSVs and a summary text file
#'
#' Produces `imputed_data.csv` (practical), `imputation_diagnostics.csv`
#' (per-study), `mi_meta_results.csv` (if a meta-analysis was fit), and
#' `pooled_results.txt` (human-readable summary).
#'
#' @param fit A `scaleL` object.
#' @param dir Directory to write to; created if missing.
#' @return Invisibly, the vector of file paths written.
#' @export
write_scaleL <- function(fit, dir) {
  stopifnot(inherits(fit, "scaleL"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  files <- character(0)
  prac_path <- file.path(dir, "imputed_data.csv")
  diag_path <- file.path(dir, "imputation_diagnostics.csv")
  utils::write.csv(fit$practical, prac_path, row.names = FALSE)
  utils::write.csv(fit$diagnostic, diag_path, row.names = FALSE)
  files <- c(files, prac_path, diag_path)
  if (!is.null(fit$meta)) {
    meta_path <- file.path(dir, "mi_meta_results.csv")
    utils::write.csv(fit$meta$pooled_df, meta_path, row.names = FALSE)
    files <- c(files, meta_path)
  }
  summary_path <- file.path(dir, "pooled_results.txt")
  writeLines(capture_summary_text(fit), summary_path)
  files <- c(files, summary_path)
  invisible(files)
}

# Internal helper: render summary as character vector
capture_summary_text <- function(fit) {
  utils::capture.output(print(summary(fit)))
}
