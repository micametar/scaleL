test_that("fit_REML matches a known random-effects fit", {
  theta <- c(0.5, 0.55, 0.45, 0.6, 0.4)
  vi    <- rep(0.01, 5)
  fit <- fit_REML(theta, vi)
  expect_true(fit$mu > 0.4 && fit$mu < 0.6)
  expect_true(fit$tau2 >= 0)
  expect_equal(fit$k, 5)
})

test_that("loglik_across_vec favours the feasible, well-located L", {
  meas <- data.frame(mean = 3.0, sd = 0.8, n = 200)
  # anchors centred near theta = 0.5 with modest heterogeneity
  ll <- loglik_across_vec(meas, mu_hat = 0.5, tau2_hat = 0.01, L_grid = 4:12)
  expect_equal(length(ll), length(4:12))
  modal <- (4:12)[which.max(ll)]
  expect_true(modal == 5)  # (3-1)/(5-1) = 0.5 lands on mu_hat
})

test_that("group_prior_fit falls back to pooled below the anchor gate", {
  anc <- data.frame(theta = runif(20, 0.4, 0.6), v = rep(0.01, 20),
                    level = rep(c("x", "y"), each = 10))
  gp_marg <- group_prior_fit(anc, "x", mode = "marginal")
  expect_equal(gp_marg$used, "pooled")
  gp_cond <- group_prior_fit(anc, "x", mode = "conditioned", min_k = 10)
  expect_equal(gp_cond$used, "conditioned")  # x has 10 anchors
  gp_cond2 <- group_prior_fit(anc, "x", mode = "conditioned", min_k = 11)
  expect_equal(gp_cond2$used, "pooled")      # below the gate -> fall back
})
