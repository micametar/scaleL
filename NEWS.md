# scaleL 0.2.0

Major methodological expansion bringing the package in line with Piers
Steel's June 2026 revision.

## New features

* **Three likelihood methods.** `scaleL()` gains a `lik_method` argument:
  `"empirical"` (new default; Monte-Carlo reference-histogram lookup built
  internally via `extraDistr::rbbinom`, always available), `"full"` (2-D
  quadrature), and `"profile"` (fast feasibility screen, the previous
  behaviour). New functions `build_empirical_reference()`,
  `loglik_full_single()`, `loglik_empirical_single()`, `get_loglik_fn()`.
  The empirical reference is built once per `scaleL()` call and cached in a
  package-internal environment (no global assignment). Use `emp_R` to shrink
  the Monte-Carlo size for fast examples/tests.
* **Corpus (Tier 2) prior.** `prior = "corpus"` estimates a Laplace-smoothed
  (alpha = 0.5) prior from the observed-L studies in your data
  (`estimate_prior_from_corpus()`), warning below 10 observed-L studies.
* **Joint origin + L imputation.** When `scale_min` is `NA` for a study, the
  origin is imputed to its posterior mode (`posterior_range()`,
  `impute_origin()`, `estimate_range_prior_from_corpus()`) and the L pipeline
  runs unchanged. New `practical`/`diagnostic` columns `origin_imputed`,
  `scale_min_used`, `origin_prob`. Per the manuscript this step uses the
  empirical likelihood; if `lik_method != "empirical"` and an origin is
  missing, the empirical likelihood is built for the origin step (with a
  warning).
* **SD imputation module** (`R/sd_impute.R`): `sd_feasibility_and_impute()`
  (single), `sd_impute_mi()` (multiple imputation with between-imputation
  variance), `impute_reliability()`, `std_mean()`, `rbar_from_alpha_k()`, and
  the composite-variance correction helpers `to_corrected_var()`,
  `to_composite_sd()`, `assert_composite_roundtrip()` (J-correction
  `J / (1 + (J - 1) rbar)`; default rbar = 0.25; Spearman-Brown inversion).
* **Across-study-means channel** (`R/across_study.R`): `fit_REML()`,
  `loglik_across_vec()`, `anchor_table()`, `group_prior_fit()` provide a
  second likelihood over L from a REML fit to full-information anchors, gated
  at >= 10 anchors with nested moderator fallback.
* **Default `L_grid` is now `4:12`** (was `2:11`).

## Dependencies

* Added `extraDistr` to Imports (beta-binomial RNG for the empirical
  reference).

## Deviations from Piers's scripts (documented)

* The across-study-means channel is shipped as a clearly-documented, exported
  standalone API (`fit_REML`, `loglik_across_vec`, `anchor_table`,
  `group_prior_fit`) rather than being auto-wired into the main `scaleL()`
  posterior. The reference's wiring lives in a simulation harness
  (`study2_across_study.R`) that depends on tidyverse/furrr and a separate
  simulation engine; the building blocks are faithful and can be combined by
  the user, but enabling the channel on by default inside `scaleL()` was
  judged too entangled to do safely for this release.
* The empirical reference defaults match the manuscript (R = 20000), but tests
  and examples use a smaller `emp_R` (1500-2000) to stay fast. Results from a
  small reference are noisier than the manuscript's.
* The origin is resolved to its posterior **mode** (not propagated through
  multiple imputation), matching Piers's `impute_origin()`.
