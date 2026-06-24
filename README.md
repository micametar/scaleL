# scaleL

Bayesian imputation of Likert scale length L for meta-analysis.

`scaleL` recovers the missing scale length L of primary studies that
report only a mean, standard deviation, and sample size. It produces a
per-study posterior over L, multiply imputed L draws, per-measure
standardized means with calibrated standard errors, and (optionally) a
pooled random-effects estimate via Rubin's rules.

## Installation

```r
# install.packages("remotes")
remotes::install_local("scaleL")          # from a local clone
# or
remotes::install_github("micametar/scaleL")
```

Required: R >= 4.1, `readr`, `extraDistr`. Suggested (for meta-pooling):
`metafor`.

## Quick start

```r
library(scaleL)

data <- read.csv(system.file("extdata", "sample_data.csv", package = "scaleL"))

# Empirical likelihood is the default. The Monte-Carlo reference is built once
# per call; shrink it with emp_R for quick runs.
fit  <- scaleL(data, prior = "tier1", M = 50, emp_R = 5000, compute_meta = TRUE)

# Reference-free fast screen:
fit_fast <- scaleL(data, prior = "tier1", M = 50, lik_method = "profile")

print(fit)
summary(fit)

# Drop-in per-measure file for downstream meta-analysis
imputed <- as.data.frame(fit)
# e.g. metafor::rma(yi = imputed$theta, sei = imputed$se_theta)
```

## Likelihood methods

Select with `lik_method`:

- `"empirical"` (default): histogram lookup against an internal Monte-Carlo
  reference built with `extraDistr::rbbinom`. Always available (does not depend
  on your corpus). Use `emp_R` to trade accuracy for speed.
- `"full"`: 2-D quadrature parametric likelihood.
- `"profile"`: fast feasibility screen (no Monte-Carlo reference).

## Input format

A data.frame (or CSV) with one row per measure:

| column       | required | description                                        |
|--------------|----------|----------------------------------------------------|
| study_id     | yes      | groups rows that share L within a study            |
| mean         | yes      | reported sample mean                               |
| sd           | yes      | reported sample SD                                 |
| n            | yes      | sample size                                        |
| L            | optional | scale length; leave NA for studies to impute       |
| scale_min    | optional | lowest scale value (default 1; **NA triggers joint origin + L imputation**) |
| measure_id   | optional | label only                                         |
| instrument   | optional | used by `instrument_priors` if supplied            |
| n_items      | optional | item count J for the composite (J) correction      |
| alpha        | optional | coefficient alpha (recovers mean inter-item r)     |
| inter_item_r | optional | mean inter-item correlation, if known              |
| true_L       | optional | reference truth, validation only (not used)        |

## Origin imputation

If a study's `scale_min` is `NA`, `scaleL()` imputes the full range
`(origin, L)` jointly and resolves the origin to its posterior mode (using the
empirical likelihood, as in the manuscript), then runs the unchanged L
pipeline. Candidate origins default to the observed `scale_min` values among
known-L studies (override with `candidate_origins`). The outputs gain
`origin_imputed`, `scale_min_used`, and `origin_prob`.

## SD imputation

For meta-analyses missing standard deviations on bounded scales, the
`sd_feasibility_and_impute()` (single) and `sd_impute_mi()` (multiple
imputation, Rubin) functions impute SDs by modelling the dispersion fraction
`phi = sd / sqrt(Bhatia-Davis bound)` and re-expanding through each target's
ceiling. `impute_reliability()` and the composite-correction helpers apply the
`J / (1 + (J - 1) rbar)` variance inflation (default `rbar = 0.25`).

## Across-study-means channel

`fit_REML()`, `loglik_across_vec()`, `anchor_table()`, and `group_prior_fit()`
build a second likelihood over L from a REML fit to full-information anchor
studies, gated at >= 10 anchors with nested moderator fallback. Exported as a
standalone API (see `NEWS.md`).

## Outputs

A `scaleL` S3 object with:

- `$practical`  per-measure data.frame with `L_used`, `theta`, `se_theta`, `ess`
- `$diagnostic` per-study data.frame with the full posterior, MI draws, and
  recovery diagnostics
- `$studies`    list of per-study records (posterior, modal/Bayes L, draws,
  per-measure theta and SE)
- `$meta`       pooled meta-analytic results if `compute_meta = TRUE`,
  including `theta_bar`, `rubin_se`, `ci_lo`, `ci_hi`, `tau2`, `FMI`

Writing to disk is opt-in:

```r
write_scaleL(fit, "output_dir")   # writes the four CSV/text files
```

## Within-study uniformity clustering

Within-study L uniformity is now tested rather than assumed. Measures whose
per-measure modal-L posteriors agree are pooled into a block and imputed
jointly; measures that disagree are imputed independently. Configurable via
`min_cluster_size` (default 2) and `min_cluster_modal_prob` (default 0.5).
The `$practical` frame includes `block_id`, `block_size`, and `pooled_flag`
columns indicating each measure's block assignment.

## Priors

Tiered prior system, all supplied as `prior =` to `scaleL()`:

- `"tier1"`         generic prior over published Likert scales (default)
- `"tier1_shifted"` sensitivity arm (more mass on L = 5, 7)
- `"corpus"`        Tier 2 corpus-empirical prior: Laplace-smoothed
                      (alpha = 0.5) frequencies of the observed-L studies in
                      your data (warns below 10 observed-L studies)
- `"field:<name>"`  built-ins: `well_being`, `clinical`, `marketing`,
                      `organizational`
- `"custom"`        with `custom_prior = c("5" = 0.6, "7" = 0.4)`

For instrument-specific priors (Tier 3), pass `instrument_priors` as a
named list and ensure your input has an `instrument` column.

## Scope of applicability

`scaleL` is most useful when:

- Your meta-analysis has K less than about 200 studies missing scale length;
- Reported Likert distributions are approximately symmetric.

It is **not** appropriate when:

- The data exhibits strong floor or ceiling effects;
- Responses are bimodal or polarized (e.g., political-attitude items
  near the extremes);
- Strong acquiescence response styles dominate the reported moments;
- K is large (around 200 studies or more): listwise deletion may then
  outperform imputation. In that regime, report both as a sensitivity
  analysis.

The variance-bound likelihood used here is parametric and relies on
within-study unimodality of the latent response distribution. Violations
of those assumptions degrade both the per-study posterior and the
pooled FMI estimate.

## License

GPL-3.
