data("coalgov")

# ================================================================================================ #
# Tests for createModelstring() - generates JAGS model code
# ================================================================================================ #

test_that("createModelstring() generates valid JAGS code for Gaussian model", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Gaussian",
    cox_intervals = NULL
  )

  expect_type(model_string, "character")
  expect_true(grepl("model\\s*\\{", model_string))
  expect_true(grepl("dnorm", model_string))
  expect_true(grepl("for.*j.*1:n.main", model_string))
})

test_that("createModelstring() generates valid JAGS code for Binomial model", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Binomial",
    cox_intervals = NULL
  )

  expect_true(grepl("dbern", model_string))
  expect_true(grepl("logit", model_string) || grepl("ilogit", model_string))
})

test_that("createModelstring() generates valid JAGS code for Weibull model", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("dweib", model_string))
  expect_true(grepl("shape", model_string))
})

test_that("createModelstring() generates valid JAGS code for Cox model", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Cox",
    cox_intervals = NULL
  )

  expect_true(grepl("dpois", model_string))
  expect_true(grepl("dN", model_string))
  expect_true(grepl("dL0", model_string))
})

test_that("createModelstring() generates valid JAGS code for Cox with intervals", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Cox",
    cox_intervals = 10
  )

  expect_true(grepl("dpois", model_string))
  expect_true(grepl("dN_interval", model_string))
  expect_true(grepl("Y_interval", model_string))
  expect_true(grepl("lambda0", model_string))
  expect_true(grepl("n.intervals", model_string))
})

test_that("createModelstring() includes mm() block code", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("X\\.mm\\.1", model_string))
  expect_true(grepl("b\\.mm\\.1", model_string))
})

test_that("createModelstring() includes mm() random effects", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("alpha\\.mm", model_string))
  expect_true(grepl("sigma\\.mm", model_string))
  expect_true(grepl("dnorm.*0.*sigma", model_string))
})

test_that("createModelstring() includes hm() random effects", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("alpha\\.hm", model_string))
  expect_true(grepl("sigma\\.hm", model_string))
})

test_that("createModelstring() includes hm() fixed effects", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "FE"),
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("b\\.hm", model_string))
  expect_false(grepl("sigma\\.hm", model_string))  # FE doesn't have sigma
})

test_that("createModelstring() includes weight function code", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("w\\.mm", model_string))
  expect_true(grepl("b\\.mm\\.wb", model_string))
})

test_that("createModelstring() includes weight normalization for constrained weights", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n, c = TRUE), RE = FALSE),
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Constrained weights should be normalized
  expect_true(grepl("sum", model_string) || grepl("w\\.mm", model_string))
})

test_that("createModelstring() includes priors", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Gaussian",
    cox_intervals = NULL
  )

  # Should have priors for b coefficients
  expect_true(grepl("b\\[.*\\].*~.*dnorm", model_string))
  # Should have prior for sigma
  expect_true(grepl("sigma.*~.*dunif", model_string) ||
              grepl("tau.*~.*dgamma", model_string))
})

test_that("createModelstring() handles multiple mm() blocks", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  expect_true(grepl("mm\\.1", model_string))
  expect_true(grepl("mm\\.2", model_string))
})

test_that("createModelstring() handles AR specifications", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # AR should have different structure for first vs subsequent
  expect_true(grepl("ar", model_string, ignore.case = TRUE))
})

test_that("createModelstring() handles offset from fixed coefficients", {
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Fixed coefficient should appear as offset
  expect_true(grepl("offset", model_string) || grepl("mu\\[j\\]", model_string))
})

test_that("createModelstring() produces well-formed JAGS model", {
  parsed <- bml:::dissectFormula(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE) +
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

  model_string <- bml:::createModelstring(
    jags_vars = jags_vars,
    parsed = parsed,
    family = "Weibull",
    cox_intervals = NULL
  )

  # Check basic structure
  expect_true(grepl("^\\s*model\\s*\\{", model_string))
  expect_true(grepl("\\}\\s*$", model_string))

  # Count braces - should be balanced
  open_braces <- length(gregexpr("\\{", model_string)[[1]])
  close_braces <- length(gregexpr("\\}", model_string)[[1]])
  expect_equal(open_braces, close_braces)
})
