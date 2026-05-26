test_that("resolve_prior normalizes built-in priors", {
  p <- resolve_prior("tier1", NULL, 2:11)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0))
})

test_that("resolve_prior accepts field and custom specs", {
  p <- resolve_prior("field:organizational", NULL, 2:11)
  expect_equal(sum(p), 1)
  pc <- resolve_prior("custom", c("5" = 0.6, "7" = 0.4), 2:11)
  expect_equal(sum(pc), 1)
  expect_equal(unname(pc["5"]), 0.6)
})

test_that("resolve_prior errors on unknown spec or zero mass", {
  expect_error(resolve_prior("bogus", NULL, 2:11))
  expect_error(resolve_prior("custom", c("99" = 1), 2:11))
  expect_error(resolve_prior("custom", NULL, 2:11))
})
