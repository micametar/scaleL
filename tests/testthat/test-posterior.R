test_that("sigma2_max returns the Bhatia-Davis bound", {
  # mu = 3 on [1, 5] -> (3-1)*(5-3) = 4
  expect_equal(scaleL:::sigma2_max(3, L = 5, scale_min = 1), 4)
  expect_equal(scaleL:::sigma2_max(5, L = 11, scale_min = 0), 25)
})

test_that("posterior_L returns a normalized vector over L_grid (profile)", {
  meas <- data.frame(mean = 3.5, sd = 0.8, n = 100)
  prior <- resolve_prior("tier1", NULL, 4:12)
  post <- posterior_L(meas, prior, 4:12,
                      loglik_fn = scaleL:::loglik_profile_single)
  expect_equal(sum(post), 1, tolerance = 1e-10)
  expect_named(post, as.character(4:12))
})

test_that("all three likelihood methods agree on feasibility and rank true L", {
  # Clean 1..5 data: mean 4.5 (infeasible for L = 4, whose max is 4), sd 1.0,
  # large n. True L = 5.
  meas <- data.frame(mean = 4.5, sd = 1.0, n = 400)
  prior <- resolve_prior("tier1", NULL, 4:12)
  REF <- build_empirical_reference(L_grid = 4:12, R = 1500, verbose = FALSE)

  fns <- list(
    profile   = scaleL:::loglik_profile_single,
    full      = scaleL:::loglik_full_single,
    empirical = function(...) scaleL:::loglik_empirical_single(..., REF = REF)
  )
  for (nm in names(fns)) {
    post <- posterior_L(meas, prior, 4:12, loglik_fn = fns[[nm]])
    modal <- as.integer(names(post)[which.max(post)])
    # feasibility agreement: L = 4 is structurally impossible (mean 4.5 > max 4)
    expect_equal(unname(post["4"]), 0, info = nm)
    # true L = 5 ranked highest by all three methods
    expect_true(modal == 5, info = paste("method", nm, "modal was", modal))
  }
})

test_that("empirical likelihood errors clearly without a reference", {
  expect_error(
    scaleL:::loglik_empirical_single(3, 0.64, 100, L = 5, REF = NULL),
    "reference")
})
