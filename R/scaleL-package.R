#' scaleL: Bayesian Imputation of Likert Scale Length for Meta-Analysis
#'
#' Bayesian imputation of missing Likert scale length L in meta-analytic
#' data sets, with calibrated uncertainty propagation via multiple
#' imputation and Rubin's rules.
#'
#' The main user-facing function is [scaleL()]. Key capabilities:
#' \itemize{
#'   \item Three likelihood methods selected by `lik_method`: `"empirical"`
#'     (default; internal Monte-Carlo reference histograms), `"full"` (2-D
#'     quadrature), `"profile"` (fast feasibility screen).
#'   \item A structured prior over L: `"tier1"`, `"tier1_shifted"`,
#'     `"corpus"` (Laplace-smoothed observed-L frequencies),
#'     `"field:<name>"`, instrument-specific, or `"custom"`.
#'   \item Joint origin + L imputation when `scale_min` is `NA`.
#'   \item SD imputation on bounded scales ([sd_feasibility_and_impute()],
#'     [sd_impute_mi()]) with composite (J) correction.
#'   \item An across-study-means channel ([fit_REML()],
#'     [loglik_across_vec()]) for a second likelihood over L.
#' }
#' The default candidate grid is `L_grid = 4:12`. See the README for usage
#' examples and scope-of-applicability notes.
#'
#' @keywords internal
#' @importFrom stats dnorm dchisq dbeta qt setNames var lm median predict
#'   rbeta rchisq rnorm runif pnorm
#' @importFrom utils write.csv
#' @importFrom graphics barplot image
#' @importFrom extraDistr rbbinom
"_PACKAGE"
