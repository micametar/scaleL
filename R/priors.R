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

#' Resolve a prior specification into a normalized prior vector
#'
#' Implements the tier system: tier 1 generic, tier 1 shifted, tier 2
#' (field), or user-supplied custom prior. The resulting vector is defined
#' over the full candidate L grid, with zero mass on L values not in the
#' named input prior.
#'
#' @param prior_type One of `"tier1"`, `"tier1_shifted"`, `"field:<name>"`,
#'   or `"custom"`.
#' @param custom_prior Named numeric vector. Required if
#'   `prior_type = "custom"`.
#' @param L_grid Integer vector of candidate L values.
#' @param field_priors Named list of field-specific priors (defaults to
#'   [FIELD_PRIORS]).
#' @return Named numeric vector summing to 1, indexed by `as.character(L_grid)`.
#' @examples
#' resolve_prior("tier1", NULL, 2:11)
#' @export
resolve_prior <- function(prior_type, custom_prior, L_grid,
                          field_priors = FIELD_PRIORS) {
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
  } else if (identical(prior_type, "custom")) {
    if (is.null(custom_prior)) {
      stop("prior_type = 'custom' requires custom_prior to be supplied.")
    }
    p <- custom_prior
  } else {
    stop("Unknown prior_type '", prior_type,
         "'. Use tier1, tier1_shifted, field:<name>, or custom.")
  }

  prior_vec <- setNames(rep(0, length(L_char)), L_char)
  prior_vec[names(p)[names(p) %in% L_char]] <-
    p[names(p)[names(p) %in% L_char]]
  if (sum(prior_vec) <= 0) {
    stop("Prior places zero mass on every candidate L. Check prior_type.")
  }
  prior_vec / sum(prior_vec)
}
