test_that("scaleL runs end-to-end on the bundled sample data", {
  path <- system.file("extdata", "sample_data.csv", package = "scaleL")
  skip_if(path == "", "sample_data.csv not installed")
  d <- read.csv(path)
  fit <- scaleL(d, prior = "tier1", M = 10, seed = 1, compute_meta = FALSE,
                lik_method = "profile")
  expect_s3_class(fit, "scaleL")
  expect_equal(nrow(fit$practical), nrow(d))
  expect_equal(nrow(fit$diagnostic), length(unique(d$study_id)))
  expect_true(all(c("L_used", "theta", "se_theta", "ess") %in%
                    names(fit$practical)))
  # Per-measure SEs must be positive and finite
  expect_true(all(is.finite(fit$practical$se_theta) &
                    fit$practical$se_theta > 0))
})

test_that("scaleL accepts synthetic fixture with mix of observed/imputed", {
  d <- make_synthetic_data()
  fit <- scaleL(d, prior = "tier1", M = 5, seed = 42, lik_method = "profile")
  imp_studies <- vapply(fit$studies, function(s) s$is_imputed, logical(1))
  expect_true(sum(imp_studies) == 1)
})

test_that("print and summary work without error", {
  d <- make_synthetic_data()
  fit <- scaleL(d, prior = "tier1", M = 5, seed = 42, lik_method = "profile")
  expect_output(print(fit), "scaleL fit")
  expect_output(print(summary(fit)), "summary")
})

test_that("as.data.frame returns the practical file", {
  d <- make_synthetic_data()
  fit <- scaleL(d, M = 5, seed = 1, lik_method = "profile")
  out <- as.data.frame(fit)
  expect_equal(out, fit$practical)
})
