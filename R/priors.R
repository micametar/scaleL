#' Tier 1 generic prior over Likert scale length
#'
#' Generic prior over L for published Likert scales, concentrated on L = 5
#' and L = 7. Sources: South et al. (2022), Simms et al. (2019), Preston &
#' Colman (2000), Taherdoost (2019), Willits et al. (2016).
#'
#' @format Named numeric vector; names are candidate L (as character),
#'   values are prior probabilities summing to 1.
#' @export
TIER1_PRIOR <- c("4" = 0.035, "5" = 0.650, "6" = 0.020, "7" = 0.270,
                 "8" = 0.002, "9" = 0.004, "10" = 0.012, "11" = 0.005,
                 "12" = 0.002)

#' Tier 1 shifted prior (sensitivity arm)
#'
#' Variant of [TIER1_PRIOR] with more mass on L = 5 and L = 7, for use as
#' a sensitivity-analysis arm.
#'
#' @format Named numeric vector.
#' @export
TIER1_SHIFTED_PRIOR <- c("4" = 0.020, "5" = 0.720, "6" = 0.010, "7" = 0.222,
                         "8" = 0.001, "9" = 0.002, "10" = 0.015, "11" = 0.008,
                         "12" = 0.002)

#' Field-specific priors (Tier 2)
#'
#' Built-in field-specific priors. Replace or extend with field knowledge
#' when available.
#'
#' @format Named list of named numeric vectors. Available fields:
#'   `well_being`, `clinical`, `marketing`, `organizational`.
#' @export
FIELD_PRIORS <- list(
  well_being = c("4" = 0.30, "5" = 0.25, "7" = 0.10, "10" = 0.15, "11" = 0.20),
  clinical   = c("3" = 0.20, "4" = 0.30, "5" = 0.45, "7" = 0.05),
  marketing  = c("5" = 0.45, "7" = 0.40, "10" = 0.05, "11" = 0.10),
  organizational = c("5" = 0.55, "7" = 0.40, "6" = 0.02, "4" = 0.03)
)

#' Estimate a Tier 2 corpus-empirical prior over L
#'
#' Tabulates observed L values in the analyst's corpus over the candidate
#' grid and applies Laplace smoothing with pseudocount `alpha` so that L
#' values absent from the corpus still receive small positive mass:
#' `pi(L) = (n_L + alpha) / (K_obs + |grid| * alpha)`.
#'
#' @param observed_L Integer vector of observed (study-level) L values; NAs
#'   are dropped.
#' @param L_grid Integer vector of candidate L values.
#' @param alpha Laplace pseudocount (default 0.5).
#' @return Named numeric vector indexed by L (as character), summing to 1.
#' @export
estimate_prior_from_corpus <- function(observed_L, L_grid, alpha = 0.5) {
  observed_L <- observed_L[!is.na(observed_L)]
  if (length(observed_L) == 0L) {
    stop("prior = 'corpus' requires at least one observed L value, but no ",
         "studies have L reported.")
  }
  L_char <- as.character(L_grid)
  counts <- setNames(rep(0, length(L_char)), L_char)
  for (L_val in observed_L) {
    key <- as.character(L_val)
    if (key %in% L_char) counts[key] <- counts[key] + 1
  }
  out_of_grid <- !as.character(observed_L) %in% L_char
  if (any(out_of_grid)) {
    bad_vals <- sort(unique(observed_L[out_of_grid]))
    warning("prior = 'corpus': ", sum(out_of_grid), " observed L value(s) ",
            "fall outside L_grid and were ignored: ",
            paste(bad_vals, collapse = ", "),
            ". Consider extending L_grid.")
  }
  smoothed <- counts + alpha
  smoothed / sum(smoothed)
}

#' Resolve a prior specification into a normalized prior vector
#'
#' Implements the tier system: tier 1 generic, tier 1 shifted, tier 2
#' corpus-empirical (Laplace-smoothed observed-L frequencies), field-flavoured
#' tier 1 variants, or user-supplied custom prior. The resulting vector is
#' defined over the full candidate L grid, with zero mass on L values not in
#' the named input prior.
#'
#' @param prior_type One of `"tier1"`, `"tier1_shifted"`, `"corpus"`,
#'   `"field:<name>"`, or `"custom"`.
#' @param custom_prior Named numeric vector. Required if
#'   `prior_type = "custom"`.
#' @param L_grid Integer vector of candidate L values.
#' @param field_priors Named list of field-specific priors (defaults to
#'   [FIELD_PRIORS]).
#' @param observed_L Integer vector of observed study-level L values; required
#'   for `prior_type = "corpus"`.
#' @param laplace_alpha Laplace pseudocount for the corpus prior (default 0.5).
#' @param corpus_min_n Minimum observed-L study count below which the corpus
#'   prior triggers a warning (default 10).
#' @return Named numeric vector summing to 1, indexed by `as.character(L_grid)`.
#' @examples
#' resolve_prior("tier1", NULL, 4:12)
#' resolve_prior("corpus", NULL, 4:12, observed_L = c(5, 5, 7, 5, 7))
#' @export
resolve_prior <- function(prior_type, custom_prior, L_grid,
                          field_priors = FIELD_PRIORS,
                          observed_L = NULL,
                          laplace_alpha = 0.5,
                          corpus_min_n = 10L) {
  L_char <- as.character(L_grid)
  if (identical(prior_type, "tier1")) {
    p <- TIER1_PRIOR
  } else if (identical(prior_type, "tier1_shifted")) {
    p <- TIER1_SHIFTED_PRIOR
  } else if (startsWith(prior_type, "field:")) {
    field_name <- sub("^field:", "", prior_type)
    if (!field_name %in% names(field_priors)) {
      stop("Unknown field prior '", field_name,
           "'. Available: ", paste(names(field_priors), collapse = ", "))
    }
    p <- field_priors[[field_name]]
  } else if (identical(prior_type, "corpus")) {
    if (is.null(observed_L)) {
      stop("prior = 'corpus' requires observed_L; this is supplied ",
           "automatically when called via scaleL().")
    }
    p <- estimate_prior_from_corpus(observed_L, L_grid, alpha = laplace_alpha)
    n_corpus <- sum(!is.na(observed_L))
    if (n_corpus < corpus_min_n) {
      warning("prior = 'corpus' was estimated from only ", n_corpus,
              " observed-L studies (below corpus_min_n = ", corpus_min_n,
              "). The empirical frequencies are noisy; consider ",
              "prior = 'tier1' if the observed-L subset is small.")
    }
  } else if (identical(prior_type, "custom")) {
    if (is.null(custom_prior)) {
      stop("prior_type = 'custom' requires custom_prior to be supplied.")
    }
    p <- custom_prior
  } else {
    stop("Unknown prior_type '", prior_type,
         "'. Use tier1, tier1_shifted, corpus, field:<name>, or custom.")
  }

  prior_vec <- setNames(rep(0, length(L_char)), L_char)
  prior_vec[names(p)[names(p) %in% L_char]] <-
    p[names(p)[names(p) %in% L_char]]
  if (sum(prior_vec) <= 0) {
    stop("Prior places zero mass on every candidate L. Check prior_type.")
  }
  prior_vec / sum(prior_vec)
}

#' Estimate a corpus prior over ranges (origin, L)
#'
#' Mirrors [estimate_prior_from_corpus()] but the unit is the full range
#' `(a, L)` (lower endpoint and length). Laplace-smoothed over the candidate
#' origin x L grid. If no study reports an origin, falls back to a flat origin
#' prior times the supplied L prior so the machinery still runs.
#'
#' @param obs_scale_min Observed lower endpoints (NA where missing).
#' @param obs_L Observed L values (NA where missing).
#' @param origin_grid Candidate origin values.
#' @param L_grid Candidate L values.
#' @param alpha Laplace pseudocount (default 0.5).
#' @param base_L_prior Optional named L prior used for the flat-origin
#'   fallback.
#' @return Named numeric vector indexed by `"origin|L"`, summing to 1.
#' @export
estimate_range_prior_from_corpus <- function(obs_scale_min, obs_L,
                                             origin_grid, L_grid,
                                             alpha = 0.5,
                                             base_L_prior = NULL) {
  keys <- character(0)
  for (a in origin_grid) for (L in L_grid) keys <- c(keys, paste(a, L, sep = "|"))
  counts <- setNames(rep(0, length(keys)), keys)
  ok <- !is.na(obs_scale_min) & !is.na(obs_L)
  if (any(ok)) {
    for (i in which(ok)) {
      key <- paste(obs_scale_min[i], obs_L[i], sep = "|")
      if (key %in% keys) counts[key] <- counts[key] + 1
    }
    smoothed <- counts + alpha
    return(smoothed / sum(smoothed))
  }
  Lp <- if (is.null(base_L_prior))
          setNames(rep(1, length(L_grid)), as.character(L_grid)) else base_L_prior
  for (a in origin_grid) for (L in L_grid)
    counts[paste(a, L, sep = "|")] <- as.numeric(Lp[as.character(L)])
  counts <- counts + alpha
  counts / sum(counts)
}

#' Joint posterior over (origin, L) for a study with unknown lower endpoint
#'
#' Generalizes [posterior_L()] from lengths to ranges: each candidate
#' `(a, L)` is scored by the same per-measure likelihood evaluated at
#' `scale_min = a`, weighted by the corpus range prior.
#'
#' @param measures Data.frame with `mean`, `sd`, `n`.
#' @param range_prior Named numeric vector from
#'   [estimate_range_prior_from_corpus()].
#' @param L_grid Candidate L values.
#' @param origin_grid Candidate origin values.
#' @param loglik_fn Per-measure log-likelihood (default [get_loglik_fn()]).
#'   Origin imputation requires the empirical likelihood in the manuscript.
#' @return Data.frame with columns `scale_min`, `L`, `post` (joint, sums to 1).
#' @export
posterior_range <- function(measures, range_prior, L_grid, origin_grid,
                            loglik_fn = get_loglik_fn()) {
  rows <- list()
  for (a in origin_grid) {
    ll <- vapply(L_grid, function(L) {
      sum(vapply(seq_len(nrow(measures)), function(i) {
        loglik_fn(measures$mean[i], measures$sd[i]^2, measures$n[i],
                  L = L, scale_min = a)
      }, numeric(1)))
    }, numeric(1))
    lp <- vapply(L_grid, function(L) {
      key <- paste(a, L, sep = "|")
      if (key %in% names(range_prior)) log(range_prior[[key]]) else -Inf
    }, numeric(1))
    rows[[as.character(a)]] <- data.frame(scale_min = a, L = L_grid,
                                          logpost = ll + lp)
  }
  res <- do.call(rbind, rows)
  finite <- is.finite(res$logpost)
  if (!any(finite)) {
    res$post <- vapply(seq_len(nrow(res)), function(j) {
      key <- paste(res$scale_min[j], res$L[j], sep = "|")
      if (key %in% names(range_prior)) range_prior[[key]] else 0
    }, numeric(1))
    res$post <- res$post / sum(res$post)
    return(res[, c("scale_min", "L", "post")])
  }
  m <- max(res$logpost[finite])
  res$post <- exp(res$logpost - m)
  res$post[!is.finite(res$post)] <- 0
  res$post <- res$post / sum(res$post)
  res[, c("scale_min", "L", "post")]
}

#' Impute a missing scale origin to its posterior mode
#'
#' Marginalizes the joint range posterior over L and returns the modal origin
#' and its posterior probability. Resolving the origin to its mode (rather than
#' propagating it through multiple imputation) keeps the downstream L pipeline
#' unchanged: the study simply enters it with an imputed `scale_min`.
#'
#' @inheritParams posterior_range
#' @return List with `scale_min` (modal origin) and `prob`.
#' @export
impute_origin <- function(measures, range_prior, L_grid, origin_grid,
                          loglik_fn = get_loglik_fn()) {
  pr  <- posterior_range(measures, range_prior, L_grid, origin_grid, loglik_fn)
  agg <- tapply(pr$post, pr$scale_min, sum)
  a_hat <- as.numeric(names(agg))[which.max(agg)]
  list(scale_min = a_hat, prob = as.numeric(max(agg)))
}
