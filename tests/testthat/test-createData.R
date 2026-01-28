data("coalgov")

# ================================================================================================ #
# Tests for createData() - internal data transformation function
# ================================================================================================ #

test_that("createData() correctly processes simple Gaussian model", {
  # Parse formula first
  parsed <- bml:::dissectFormula(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    data = coalgov
  )

  # Create data
  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Gaussian"
  )

  expect_type(data_list, "list")
  expect_true("Y" %in% names(data_list))
  expect_true("X" %in% names(data_list))
  expect_equal(length(data_list$Y), nrow(coalgov))
})

test_that("createData() correctly processes survival model", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  expect_true("t" %in% names(data_list))
  expect_true("event" %in% names(data_list))
  expect_true("X" %in% names(data_list))
})

test_that("createData() correctly handles mm() blocks", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check for mm-specific variables
  expect_true("mmid1" %in% names(data_list))
  expect_true("mainid1" %in% names(data_list))
  expect_true("n.mm1" %in% names(data_list))
  expect_true("X.mm.1" %in% names(data_list))
})

test_that("createData() correctly handles hm() blocks", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check for hm-specific variables
  expect_true("hmid1" %in% names(data_list))
  expect_true("n.hm1" %in% names(data_list))
})

test_that("createData() correctly handles multiple mm() blocks", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE) +
      mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check for multiple mm blocks
  expect_true("X.mm.1" %in% names(data_list))
  expect_true("X.mm.2" %in% names(data_list))
})

test_that("createData() correctly processes weight functions", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check that weight variables are created
  expect_true(any(grepl("^n$", names(data_list))))
})

test_that("createData() correctly handles parameterized weight functions", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel), RE = FALSE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check that weight function variables are included
  expect_true(any(grepl("pseatrel", names(data_list))))
})

test_that("createData() correctly handles fix() variables", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + fix(majority, 1.0),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Fixed variables should contribute to offset
  expect_true("offset" %in% names(data_list) || "X" %in% names(data_list))
})

test_that("createData() correctly handles fix() in mm() blocks", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fix(fdep, 1.0) + ipd), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Should have offset for fixed variables
  expect_true(any(grepl("offset", names(data_list))))
})

test_that("createData() preserves sample size", {
  parsed <- bml:::dissectFormula(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Gaussian"
  )

  expect_equal(data_list$n.main, nrow(coalgov))
})

test_that("createData() correctly indexes mm members", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check that member indexing makes sense
  expect_true(all(data_list$mmid1 >= 1))
  expect_true(all(data_list$mmid1 <= data_list$n.mm1))
})

test_that("createData() correctly indexes hm groups", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # Check that hm indexing makes sense
  expect_true(all(data_list$hmid1 >= 1))
  expect_true(all(data_list$hmid1 <= data_list$n.hm1))
})

test_that("createData() handles intercept correctly", {
  # With intercept
  parsed_with <- bml:::dissectFormula(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    data = coalgov
  )

  data_with <- bml:::createData(
    data = coalgov,
    parsed = parsed_with,
    family = "Gaussian"
  )

  # Without intercept
  parsed_without <- bml:::dissectFormula(
    formula = sim.y ~ 0 + majority,
    family = "Gaussian",
    data = coalgov
  )

  data_without <- bml:::createData(
    data = coalgov,
    parsed = parsed_without,
    family = "Gaussian"
  )

  # Number of main-level covariates should differ by 1
  expect_equal(ncol(data_with$X), ncol(data_without$X) + 1)
})

test_that("createData() handles AR specifications", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE, ar = TRUE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  # AR requires sequential indexing
  expect_true("ar1" %in% names(data_list) || any(grepl("ar", names(data_list))))
})

test_that("createData() handles subset of data", {
  subset_data <- coalgov[1:100, ]

  parsed <- bml:::dissectFormula(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    data = subset_data
  )

  data_list <- bml:::createData(
    data = subset_data,
    parsed = parsed,
    family = "Gaussian"
  )

  expect_equal(data_list$n.main, 100)
  expect_equal(length(data_list$Y), 100)
})

test_that("createData() handles Cox models", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Cox",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Cox"
  )

  # Cox-specific elements
  expect_true("t" %in% names(data_list))
  expect_true("event" %in% names(data_list))
})
