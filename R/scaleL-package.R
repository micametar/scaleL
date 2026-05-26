#' scaleL: Bayesian Imputation of Likert Scale Length for Meta-Analysis
#'
#' Implements the Steel & Fariborzi (2026) procedure for imputing missing
#' Likert scale length L in meta-analytic data sets, with calibrated
#' uncertainty propagation via multiple imputation and Rubin's rules.
#'
#' The main user-facing function is [scaleL()]. See the package README for
#' usage examples and scope-of-applicability notes.
#'
#' @keywords internal
#' @importFrom stats dnorm dchisq qt setNames var
#' @importFrom utils write.csv
"_PACKAGE"
