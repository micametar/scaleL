make_sd_frame <- function(seed = 20260520L, K = 60, hide = 18) {
  set.seed(seed)
  L  <- sample(c(5, 7), K, replace = TRUE)
  mu <- runif(K, 2, L - 1)
  s  <- sqrt(runif(K, 0.05, 0.6) * (mu - 1) * (L - mu))  # true SD below ceiling
  d <- data.frame(mean = mu, sd = s, n = sample(80:400, K, TRUE),
                  items = sample(3:8, K, TRUE),
                  response_low = 1, response_high = L, L = L,
                  alpha = runif(K, 0.6, 0.9),
                  dimension = sample(letters[1:4], K, TRUE),
                  stringsAsFactors = FALSE)
  d$sd[sample(K, hide)] <- NA
  list(d = d, L = L)
}

test_that("composite-correction round-trips", {
  expect_true(assert_composite_roundtrip(c(0.8, 1.2, NA), c(5, 7, 4),
                                         c(0.25, 0.3, 0.25)))
})

test_that("reliability/J-correction inverts Spearman-Brown", {
  # rbar = alpha / (k - alpha (k-1)); alpha = 0.8, k = 5 -> 0.444...
  rb <- rbar_from_alpha_k(0.8, 5)
  expect_equal(rb, 0.8 / (5 - 0.8 * 4), tolerance = 1e-8)
  # Row 2 has neither alpha nor a field/global fallback (its dimension never
  # reports alpha) -> default rbar.
  d <- data.frame(alpha = c(0.8, NA), items = c(5, 4), n = c(100, 100),
                  dimension = c("a", "b"))
  d <- impute_reliability(d)
  expect_true(all(d$rbar >= 0.25 - 1e-9))      # default floor
  expect_equal(d$rbar_source[2], "global")     # filled from global alpha
})

test_that("sd_feasibility_and_impute fills missing SDs within the bound", {
  fr <- make_sd_frame()
  out <- suppressMessages(sd_feasibility_and_impute(fr$d))
  imp <- which(out$sd_imputed)
  expect_true(length(imp) > 0)
  bound <- (out$mean - out$response_low) * (out$response_high - out$mean)
  expect_true(all(out$sd[imp] <= sqrt(bound[imp]) + 1e-8))  # respects ceiling
})

test_that("sd_impute_mi recovers SDs, respects the bound, returns between-imp var", {
  fr <- make_sd_frame()
  res <- sd_impute_mi(fr$d, M = 30, seed = 7)
  need <- res$need
  expect_true(length(need) > 0)
  bound <- (res$d$mean - res$d$response_low) *
    (res$d$response_high - res$d$mean)
  expect_true(all(res$d$sd[need] <= sqrt(bound[need]) + 1e-8))
  # between-imputation variance is recorded and positive
  expect_true(all(is.finite(res$d$sd_impute_var[need])))
  expect_true(all(res$d$sd_impute_var[need] > 0))
})
