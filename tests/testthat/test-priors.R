test_that("resolve_prior normalizes built-in priors", {
  p <- resolve_prior("tier1", NULL, 4:12)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0))
})

test_that("resolve_prior accepts field and custom specs", {
  p <- resolve_prior("field:organizational", NULL, 4:12)
  expect_equal(sum(p), 1)
  pc <- resolve_prior("custom", c("5" = 0.6, "7" = 0.4), 4:12)
  expect_equal(sum(pc), 1)
  expect_equal(unname(pc["5"]), 0.6)
})

test_that("resolve_prior errors on unknown spec or zero mass", {
  expect_error(resolve_prior("bogus", NULL, 4:12))
  expect_error(resolve_prior("custom", c("99" = 1), 4:12))
  expect_error(resolve_prior("custom", NULL, 4:12))
})

test_that("corpus prior is Laplace-smoothed and concentrates on observed L", {
  obs <- c(5, 5, 5, 7, 5, 7, 5, 5, 7, 5, 5, 7)  # mostly 5
  p <- resolve_prior("corpus", NULL, 4:12, observed_L = obs)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_true(p["5"] == max(p))
  expect_true(all(p > 0))  # Laplace smoothing -> no zeros
})

test_that("corpus prior warns when too few observed-L studies", {
  expect_warning(resolve_prior("corpus", NULL, 4:12, observed_L = c(5, 7)),
                 "observed-L studies")
})
