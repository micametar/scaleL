# scaleL

Bayesian imputation of Likert scale length L for meta-analysis.
Companion R package to:

> Steel, P., & Fariborzi, H. (2026). *Bayesian Imputation of Likert Scale
> Length in Meta-Analysis.*

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

Required: R >= 4.1, `readr`. Suggested (for meta-pooling): `metafor`.

## Quick start

```r
library(scaleL)

data <- read.csv(system.file("extdata", "sample_data.csv", package = "scaleL"))
fit  <- scaleL(data, prior = "tier1", M = 50, compute_meta = TRUE)

print(fit)
summary(fit)

# Drop-in per-measure file for downstream meta-analysis
imputed <- as.data.frame(fit)
# e.g. metafor::rma(yi = imputed$theta, sei = imputed$se_theta)
```

## Input format

A data.frame (or CSV) with one row per measure:

| column     | required | description                                        |
|------------|----------|----------------------------------------------------|
| study_id   | yes      | groups rows that share L within a study            |
| mean       | yes      | reported sample mean                               |
| sd         | yes      | reported sample SD                                 |
| n          | yes      | sample size                                        |
| L          | optional | scale length; leave NA for studies to impute       |
| scale_min  | optional | lowest scale value (default 1; use 0 for 0-10)     |
| measure_id | optional | label only                                         |
| instrument | optional | used by `instrument_priors` if supplied            |
| true_L     | optional | reference truth, validation only (not used)        |

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

## Priors

Tiered prior system, all supplied as `prior =` to `scaleL()`:

- `"tier1"`         generic prior over published Likert scales (default)
- `"tier1_shifted"` sensitivity arm (more mass on L = 5, 7)
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

## Citation

```r
citation("scaleL")
```

Or cite the paper directly:

> Steel, P., & Fariborzi, H. (2026). Bayesian Imputation of Likert Scale
> Length in Meta-Analysis.

## License

GPL-3.
