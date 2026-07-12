# ================================================================================================ #
# Family functions and validation
# ================================================================================================ #

test_that("family functions and strings normalize to internal families", {
  expect_equal(validate_family(stats::gaussian())$family, "Gaussian")
  expect_equal(validate_family(stats::binomial())$family, "Binomial")
  expect_equal(validate_family(bernoulli())$family, "Binomial")
  expect_equal(validate_family(weibull())$family, "Weibull")
  expect_equal(validate_family(cox())$family, "Cox")

  expect_equal(validate_family("gaussian")$family, "Gaussian")
  expect_equal(validate_family("bernoulli")$family, "Binomial")
  expect_equal(validate_family("binomial")$family, "Binomial")
  expect_equal(validate_family("weibull")$family, "Weibull")
  expect_equal(validate_family("cox")$family, "Cox")
  # capitalized legacy spellings still resolve
  expect_equal(validate_family("Gaussian")$family, "Gaussian")
})

test_that("cox() carries the baseline-hazard intervals", {
  expect_null(cox()$intervals)
  expect_equal(cox(intervals = 10)$intervals, 10L)
  expect_error(cox(intervals = -1), "positive")
})

test_that("unsupported input errors", {
  expect_error(validate_family("poisson"), "Invalid family")
  expect_error(validate_family(stats::poisson()), "Unsupported family")
  expect_error(validate_family(stats::gaussian(link = "log")), "identity")
  expect_error(validate_family(42), "family function")
})
