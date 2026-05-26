test_that("sigma2_max returns the Bhatia-Davis bound", {
  # mu = 3 on [1, 5] -> (3-1)*(5-3) = 4
  expect_equal(scaleL:::sigma2_max(3, L = 5, scale_min = 1), 4)
  expect_equal(scaleL:::sigma2_max(5, L = 11, scale_min = 0), 25)
})

test_that("posterior_L returns a normalized vector over L_grid", {
  meas <- data.frame(mean = 3.5, sd = 0.8, n = 100)
  prior <- resolve_prior("tier1", NULL, 2:11)
  post <- posterior_L(meas, prior, 2:11)
  expect_equal(sum(post), 1, tolerance = 1e-10)
  expect_named(post, as.character(2:11))
})

test_that("modal-L recovery on a clearly-1-to-5 synthetic case", {
  # mean=3, sd=0.8, n=300 is compatible with L=5 but not L=2/3
  meas <- data.frame(mean = 3.0, sd = 0.8, n = 300)
  prior <- resolve_prior("tier1", NULL, 2:11)
  post <- posterior_L(meas, prior, 2:11)
  modal <- as.integer(names(post)[which.max(post)])
  expect_true(modal == 5)
})
