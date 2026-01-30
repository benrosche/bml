data("coalgov")

# ================================================================================================ #
# Tests for createJagsVars() - creates variable list for JAGS
# ================================================================================================ #

# Helper function to set up jags variables (mirrors bml.R flow)
setup_jags_vars <- function(formula, family, data = coalgov, cox_intervals = NULL) {
  # 1. Parse formula
  formula_parts <- bml:::dissectFormula(formula, family, data)
  mm <- formula_parts$mm
  hm <- formula_parts$hm

  # 2. Create data structures
  data_parts <- bml:::createData(data, formula_parts)
  data      <- data_parts$data
  mm_blocks <- data_parts$mm_blocks
  main      <- data_parts$main
  hm_blocks <- data_parts$hm_blocks

  # 3. Create jags variables
  jags_vars <- bml:::createJagsVars(
    data = data,
    family = family,
    mm_blocks = mm_blocks,
    main = main,
    hm_blocks = hm_blocks,
    mm = mm,
    hm = hm,
    monitor = TRUE,
    modelfile = FALSE,
    n.chains = 2,
    inits = NULL,
    cox_intervals = cox_intervals
  )

  list(
    jags_vars = jags_vars,
    formula_parts = formula_parts,
    data_parts = data_parts
  )
}

test_that("createJagsVars() creates correct structure for Gaussian model", {
  result <- setup_jags_vars(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian"
  )

  jv <- result$jags_vars

  expect_type(jv, "list")
  expect_true("jags.data" %in% names(jv))
  expect_true("jags.params" %in% names(jv))
  expect_true("Ns" %in% names(jv))

  # Check JAGS data contains expected elements
  expect_true("Y" %in% names(jv$jags.data))
  expect_true("X.main" %in% names(jv$jags.data))
  expect_true("n.main" %in% names(jv$jags.data))
  # n.Xmain is in Ns, not jags.data (not needed by JAGS model)
  expect_true("n.Xmain" %in% names(jv$Ns))
})

test_that("createJagsVars() creates correct structure for Weibull model", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Weibull"
  )

  jv <- result$jags_vars

  expect_true("t" %in% names(jv$jags.data))
  expect_true("ct.lb" %in% names(jv$jags.data))
  expect_true("censored" %in% names(jv$jags.data))
  expect_true("X.main" %in% names(jv$jags.data))
})

test_that("createJagsVars() handles Cox model without intervals", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Cox"
  )

  jv <- result$jags_vars

  expect_true("Y" %in% names(jv$jags.data))
  expect_true("dN" %in% names(jv$jags.data))
  expect_true("t.unique" %in% names(jv$jags.data))
  expect_true("n.tu" %in% names(jv$jags.data))
})

test_that("createJagsVars() handles Cox model with intervals", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Cox",
    cox_intervals = 10
  )

  jv <- result$jags_vars

  expect_true("Y_interval" %in% names(jv$jags.data))
  expect_true("dN_interval" %in% names(jv$jags.data))
  expect_true("n.intervals" %in% names(jv$jags.data))
  expect_equal(jv$jags.data$n.intervals, 10)
})

test_that("createJagsVars() correctly handles mm() blocks", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  expect_true("X.mm.1" %in% names(jv$jags.data))
  expect_true("n.Xmm.1" %in% names(jv$jags.data))
  # mmid is only in jags.data when RE = TRUE
  expect_true("n.mm" %in% names(jv$jags.data))
  expect_true("mmi1" %in% names(jv$jags.data))
  expect_true("mmi2" %in% names(jv$jags.data))
})

test_that("createJagsVars() correctly handles mm() RE", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # RE requires unique mm count
  expect_true("n.umm" %in% names(jv$jags.data))
  expect_true(jv$jags.data$n.umm > 0)
})

test_that("createJagsVars() correctly handles hm() blocks", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull"
  )

  jv <- result$jags_vars

  expect_true("hmid" %in% names(jv$jags.data))
  expect_true("n.hm" %in% names(jv$jags.data))
})

test_that("createJagsVars() correctly handles parameterized weight functions", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel), RE = FALSE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # Weight function variables should be in X.w.1
  expect_true("X.w.1" %in% names(jv$jags.data))
})

test_that("createJagsVars() correctly handles deterministic weights (Phase 2 optimization)", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # Phase 2: When fn has no parameters, weights should be pre-computed
  # Check that w.1 is passed as data (pre-computed)
  expect_true("w.1" %in% names(jv$jags.data))
  expect_true(is.numeric(jv$jags.data$w.1))
})

test_that("createJagsVars() correctly handles weight constraints with accumulator", {
  # Constrained weights with parameters use accumulator pattern
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel, c = TRUE), RE = FALSE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # Accumulator pattern requires grp.mm
  expect_true("grp.mm" %in% names(jv$jags.data))
})

test_that("createJagsVars() correctly handles fixed coefficients (Phase 1 optimization)", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + fix(majority, 1.0),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # Fixed coefficients should produce an offset
  expect_true("offset.main" %in% names(jv$jags.data))
})

test_that("createJagsVars() correctly counts main-level parameters", {
  result <- setup_jags_vars(
    formula = sim.y ~ 1 + majority + mwc,
    family = "Gaussian"
  )

  jv <- result$jags_vars

  # 1 intercept + 2 covariates = 3 parameters
  # n.Xmain is in Ns, not jags.data
  expect_equal(jv$Ns$n.Xmain, 3)
  expect_equal(ncol(jv$jags.data$X.main), 3)
})

test_that("createJagsVars() correctly handles multiple mm() blocks", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE) +
      mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  expect_true("X.mm.1" %in% names(jv$jags.data))
  expect_true("X.mm.2" %in% names(jv$jags.data))
  expect_true("n.Xmm.1" %in% names(jv$jags.data))
  expect_true("n.Xmm.2" %in% names(jv$jags.data))
})

test_that("createJagsVars() correctly handles AR specifications", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE, ar = TRUE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # AR requires special indexing
  expect_true("n.GPn" %in% names(jv$jags.data))
  expect_true("n.GPNi" %in% names(jv$jags.data))
})

test_that("createJagsVars() dimensions are consistent", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep + ipd), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  jv <- result$jags_vars

  # Check dimensional consistency
  expect_equal(length(jv$jags.data$t), jv$jags.data$n.main)
  expect_equal(nrow(jv$jags.data$X.main), jv$jags.data$n.main)
  expect_equal(ncol(jv$jags.data$X.main), jv$Ns$n.Xmain)
})

test_that("createJagsVars() handles Binomial model", {
  result <- setup_jags_vars(
    formula = earlyterm ~ 1 + majority,
    family = "Binomial"
  )

  jv <- result$jags_vars

  expect_true("Y" %in% names(jv$jags.data))
  expect_true("X.main" %in% names(jv$jags.data))
  # Binomial Y should be 0/1
  expect_true(all(jv$jags.data$Y %in% c(0, 1)))
})

test_that("createJagsVars() returns jags.params for monitoring", {
  result <- setup_jags_vars(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian"
  )

  jv <- result$jags_vars

  expect_true(length(jv$jags.params) > 0)
  expect_true("b" %in% jv$jags.params)
  expect_true("sigma" %in% jv$jags.params)
})

test_that("createJagsVars() returns jags.inits", {
  result <- setup_jags_vars(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Weibull"
  )

  jv <- result$jags_vars

  expect_type(jv$jags.inits, "list")
})
