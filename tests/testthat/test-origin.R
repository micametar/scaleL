test_that("origin imputation recovers a known origin on synthetic NA case", {
  # Corpus of known-origin (=1), L=5 studies plus one study with NA scale_min
  # generated at origin 1. The empirical likelihood should recover origin 1.
  d <- data.frame(
    study_id = paste0("S", 1:8),
    mean = c(3.0, 3.2, 3.1, 2.9, 3.3, 3.0, 3.1, 3.0),
    sd   = c(.9, .8, .85, .9, .8, .9, .85, .88),
    n    = rep(200, 8),
    L    = rep(5L, 8),
    scale_min = c(1, 1, 1, 1, 1, 1, 1, NA)
  )
  fit <- scaleL(d, prior = "tier1", M = 5, seed = 1, emp_R = 1500,
                lik_method = "empirical")
  expect_true(fit$diagnostic$origin_imputed[8])
  expect_equal(fit$diagnostic$scale_min_used[8], 1)
})

test_that("posterior_range returns a normalized joint posterior", {
  meas <- data.frame(mean = 3.0, sd = 0.8, n = 200)
  rp <- estimate_range_prior_from_corpus(
    obs_scale_min = c(1, 1, 1), obs_L = c(5, 5, 7),
    origin_grid = c(0, 1), L_grid = 4:12)
  pr <- posterior_range(meas, rp, 4:12, c(0, 1),
                        loglik_fn = scaleL:::loglik_profile_single)
  expect_equal(sum(pr$post), 1, tolerance = 1e-10)
})
