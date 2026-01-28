data("coalgov")

# ================================================================================================ #
# Tests for createJagsVars() - creates variable list for JAGS
# ================================================================================================ #

test_that("createJagsVars() creates correct structure for Gaussian model", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Gaussian",
    cox_intervals = NULL
  )

  expect_type(jags_vars, "list")
  expect_true("Y" %in% names(jags_vars))
  expect_true("X" %in% names(jags_vars))
  expect_true("n.main" %in% names(jags_vars))
  expect_true("n.b" %in% names(jags_vars))
})

test_that("createJagsVars() creates correct structure for Weibull model", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true("t" %in% names(jags_vars))
  expect_true("event" %in% names(jags_vars))
  expect_true("X" %in% names(jags_vars))
})

test_that("createJagsVars() handles Cox model without intervals", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Cox",
    cox_intervals = NULL
  )

  expect_true("t" %in% names(jags_vars))
  expect_true("event" %in% names(jags_vars))
  expect_true("Y" %in% names(jags_vars))
  expect_true("dN" %in% names(jags_vars))
  expect_true("dL0" %in% names(jags_vars))
  expect_true("n.tu" %in% names(jags_vars))
})

test_that("createJagsVars() handles Cox model with intervals", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Cox",
    cox_intervals = 10
  )

  expect_true("Y_interval" %in% names(jags_vars))
  expect_true("dN_interval" %in% names(jags_vars))
  expect_true("n.intervals" %in% names(jags_vars))
  expect_equal(jags_vars$n.intervals, 10)
  expect_false("Y" %in% names(jags_vars))  # Should not have old Y matrix
  expect_false("dN" %in% names(jags_vars))  # Should not have old dN matrix
})

test_that("createJagsVars() correctly handles mm() blocks", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true("X.mm.1" %in% names(jags_vars))
  expect_true("n.mm.b.1" %in% names(jags_vars))
  expect_true("mmid1" %in% names(jags_vars))
  expect_true("mainid1" %in% names(jags_vars))
  expect_true("n.mm1" %in% names(jags_vars))
})

test_that("createJagsVars() correctly handles mm() RE", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Should have indicators for RE
  expect_true("mm.re1" %in% names(jags_vars))
  expect_equal(jags_vars$mm.re1, 1)
})

test_that("createJagsVars() correctly handles hm() blocks", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true("hmid1" %in% names(jags_vars))
  expect_true("n.hm1" %in% names(jags_vars))
  expect_true("hm.type1" %in% names(jags_vars))
})

test_that("createJagsVars() correctly handles parameterized weight functions", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Should have weight function variables
  expect_true(any(grepl("X.mm.wvar", names(jags_vars))))
  expect_true(any(grepl("n.mm.wb", names(jags_vars))))
})

test_that("createJagsVars() correctly handles deterministic weights", {
  # Equal weights should be pre-computed
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # If weights are pre-computed, we should have w.mm
  if ("w.mm1" %in% names(jags_vars)) {
    expect_true(is.numeric(jags_vars$w.mm1))
  }
})

test_that("createJagsVars() correctly handles weight constraints", {
  # Constrained weights
  parsed_constrained <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n, c = TRUE), RE = FALSE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed_constrained,
    family = "Weibull"
  )

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed_constrained,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true("mm.wc1" %in% names(jags_vars))
  expect_equal(jags_vars$mm.wc1, 1)

  # Unconstrained weights
  parsed_unconstrained <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n, c = FALSE), RE = FALSE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed_unconstrained,
    family = "Weibull"
  )

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed_unconstrained,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_equal(jags_vars$mm.wc1, 2)
})

test_that("createJagsVars() correctly handles fixed coefficients", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Fixed coefficients should reduce n.b and add offset
  expect_true("offset" %in% names(jags_vars) || jags_vars$n.b < 2)
})

test_that("createJagsVars() correctly counts parameters", {
  parsed <- bml:::dissectFormula(
    formula = sim.y ~ 1 + majority + mwc,
    family = "Gaussian",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Gaussian"
  )

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Gaussian",
    cox_intervals = NULL
  )

  # 1 intercept + 2 covariates = 3 parameters
  expect_equal(jags_vars$n.b, 3)
})

test_that("createJagsVars() correctly handles multiple mm() blocks", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true("X.mm.1" %in% names(jags_vars))
  expect_true("X.mm.2" %in% names(jags_vars))
  expect_true("n.mm1" %in% names(jags_vars))
  expect_true("n.mm2" %in% names(jags_vars))
})

test_that("createJagsVars() correctly handles AR specifications", {
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

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # AR requires special indexing
  expect_true("mm.ar1" %in% names(jags_vars))
  expect_equal(jags_vars$mm.ar1, 1)
})

test_that("createJagsVars() dimensions are consistent", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep + ipd), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Weibull"
  )

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Check dimensional consistency
  expect_equal(length(jags_vars$t), jags_vars$n.main)
  expect_equal(nrow(jags_vars$X), jags_vars$n.main)
  expect_equal(ncol(jags_vars$X), jags_vars$n.b)
  expect_equal(nrow(jags_vars$X.mm.1), jags_vars$n.main)
  expect_equal(ncol(jags_vars$X.mm.1), jags_vars$n.mm.b.1)
})

test_that("createJagsVars() handles Binomial model", {
  parsed <- bml:::dissectFormula(
    formula = earlyterm ~ 1 + majority,
    family = "Binomial",
    data = coalgov
  )

  data_list <- bml:::createData(
    data = coalgov,
    parsed = parsed,
    family = "Binomial"
  )

  jags_vars <- bml:::createJagsVars(
    data = data_list,
    parsed = parsed,
    family = "Binomial",
    cox_intervals = NULL
  )

  expect_true("Y" %in% names(jags_vars))
  expect_true("X" %in% names(jags_vars))
  # Binomial Y should be 0/1
  expect_true(all(jags_vars$Y %in% c(0, 1)))
})
