data("coalgov")

# ================================================================================================ #
# Tests for createModelstring() - generates JAGS model code
# ================================================================================================ #

# Helper function to set up model components (mirrors bml.R flow)
setup_model <- function(formula, family, data = coalgov, priors = NULL, cox_intervals = NULL) {
  DIR <- system.file(package = "bml")

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

  # 3. Create modelstring
  model_string <- bml:::createModelstring(
    family, priors, mm_blocks, main, hm_blocks, mm, hm,
    DIR, monitor = TRUE, modelfile = FALSE, cox_intervals
  )

  list(
    model_string = model_string,
    formula_parts = formula_parts,
    data_parts = data_parts,
    mm = mm,
    hm = hm
  )
}

test_that("createModelstring() generates valid JAGS code for Gaussian model", {
  result <- setup_model(
    formula = event_wkb ~ 1 + majority,
    family = "Gaussian"
  )

  expect_type(result$model_string, "character")
  expect_true(grepl("model\\s*\\{", result$model_string))
  expect_true(grepl("dnorm", result$model_string))
  expect_true(grepl("for.*j.*1:n.main", result$model_string))
})

test_that("createModelstring() generates valid JAGS code for Binomial model", {
  result <- setup_model(
    formula = event_wkb ~ 1 + majority,
    family = "Binomial"
  )

  expect_true(grepl("dbern", result$model_string))
  expect_true(grepl("logit", result$model_string) || grepl("ilogit", result$model_string))
})

test_that("createModelstring() generates valid JAGS code for Weibull model", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority,
    family = "Weibull"
  )

  expect_true(grepl("dweib", result$model_string))
  expect_true(grepl("shape", result$model_string))
})

test_that("createModelstring() generates valid JAGS code for Cox model", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority,
    family = "Cox"
  )

  expect_true(grepl("dpois", result$model_string))
  expect_true(grepl("dN", result$model_string))
  expect_true(grepl("dL0", result$model_string))
})

test_that("createModelstring() generates valid JAGS code for Cox with intervals", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority,
    family = "Cox",
    cox_intervals = 10
  )

  expect_true(grepl("dpois", result$model_string))
  expect_true(grepl("dN_interval", result$model_string))
  expect_true(grepl("Y_interval", result$model_string))
  expect_true(grepl("lambda0", result$model_string))
  expect_true(grepl("n.intervals", result$model_string))
})

test_that("createModelstring() includes mm() block code", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull"
  )

  expect_true(grepl("X\\.mm\\.1", result$model_string))
  expect_true(grepl("b\\.mm\\.1", result$model_string))
})

test_that("createModelstring() includes mm() random effects", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  expect_true(grepl("re\\.mm", result$model_string))
  expect_true(grepl("sigma\\.mm", result$model_string))
  expect_true(grepl("tau\\.mm", result$model_string))
})

test_that("createModelstring() includes hm() random effects", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull"
  )

  expect_true(grepl("re\\.hm\\.1", result$model_string))
  expect_true(grepl("sigma\\.hm\\.1", result$model_string))
  expect_true(grepl("tau\\.hm\\.1", result$model_string))
})

test_that("createModelstring() includes hm() fixed effects", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      hm(id = id(cid), type = "FE"),
    family = "Weibull"
  )

  expect_true(grepl("b\\.hm", result$model_string))
  expect_false(grepl("sigma\\.hm", result$model_string))  # FE doesn't have sigma
})

test_that("createModelstring() includes weight function code", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b0 + b1 * pseat), RE = FALSE),
    family = "Weibull"
  )

  # Weight variable: w.1[i]
  expect_true(grepl("w\\.1\\[", result$model_string))
  # Weight function parameters: b.w.1[1], b.w.1[2]
  expect_true(grepl("b\\.w\\.1", result$model_string))
})

test_that("createModelstring() includes weight normalization for constrained weights", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n, c = TRUE), RE = FALSE),
    family = "Weibull"
  )

  # Constrained weights should be normalized
  expect_true(grepl("sum", result$model_string) || grepl("w\\.mm", result$model_string))
})

test_that("createModelstring() includes priors", {
  result <- setup_model(
    formula = event_wkb ~ 1 + majority,
    family = "Gaussian"
  )

  # Should have priors for b coefficients
  expect_true(grepl("b\\[.*\\].*~.*dnorm", result$model_string))
  # Should have prior for tau (scaled gamma)
  expect_true(grepl("tau.*~.*dscaled\\.gamma", result$model_string))
})

test_that("createModelstring() handles multiple mm() blocks", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE) +
      mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  expect_true(grepl("mm\\.1", result$model_string))
  expect_true(grepl("mm\\.2", result$model_string))
})

test_that("createModelstring() handles AR specifications", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE, ar = TRUE),
    family = "Weibull"
  )

  # AR should have different structure for first vs subsequent
  expect_true(grepl("ar", result$model_string, ignore.case = TRUE))
})

test_that("createModelstring() handles offset from fixed coefficients", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + fix(majority, 1.0),
    family = "Weibull"
  )

  # Fixed coefficient should appear as offset
  expect_true(grepl("offset", result$model_string) || grepl("mu\\[j\\]", result$model_string))
})

test_that("createModelstring() produces well-formed JAGS model", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE) +
      hm(id = id(cid), type = "RE"),
    family = "Weibull"
  )

  # Check basic structure
  expect_true(grepl("^\\s*model\\s*\\{", result$model_string))
  expect_true(grepl("\\}\\s*$", result$model_string))

  # Count braces - should be balanced
  open_braces <- length(gregexpr("\\{", result$model_string)[[1]])
  close_braces <- length(gregexpr("\\}", result$model_string)[[1]])
  expect_equal(open_braces, close_braces)
})

# ================================================================================================ #
# Tests for accumulator pattern optimization
# ================================================================================================ #

test_that("createModelstring() uses accumulator pattern for constrained weights with parameters", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b0 + b1 * pseat, c = TRUE), RE = FALSE),
    family = "Weibull"
  )

  # Should have cumulative sum pattern (accumulator)
  expect_true(grepl("cum\\.uw\\.1\\[1\\]\\s*<-\\s*0", result$model_string),
              info = "Should initialize cumulative sum to 0")

  expect_true(grepl("cum\\.uw\\.1\\[i\\+1\\]\\s*<-\\s*cum\\.uw\\.1\\[i\\]\\s*\\+\\s*uw\\.1\\[i\\]", result$model_string),
              info = "Should have accumulator loop pattern")

  # Should have group sums computed once per group (using per-group indices)
  expect_true(grepl("sum\\.uw\\.1\\[j\\]\\s*<-\\s*cum\\.uw\\.1\\[mmi2\\.1\\[j\\]\\+1\\]\\s*-\\s*cum\\.uw\\.1\\[mmi1\\.1\\[j\\]\\]", result$model_string),
              info = "Should compute group sums using cumulative sum differences")

  # Should use grp.mm for normalization (per-group naming)
  expect_true(grepl("w\\.1\\[i\\]\\s*<-\\s*uw\\.1\\[i\\]\\s*/\\s*sum\\.uw\\.1\\[grp\\.mm\\.1\\[i\\]\\]", result$model_string),
              info = "Should normalize using grp.mm.1 indexing")
})

test_that("createModelstring() does NOT use accumulator for unconstrained weights", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b0 + b1 * pseat, c = FALSE), RE = FALSE),
    family = "Weibull"
  )

  # Should NOT have cumulative sum pattern
  expect_false(grepl("cum\\.uw", result$model_string),
               info = "Should not have cumulative sums for unconstrained weights")

  # Should NOT have group sums
  expect_false(grepl("sum\\.uw", result$model_string),
               info = "Should not have group sums for unconstrained weights")

  # Should NOT reference grp.mm
  expect_false(grepl("grp\\.mm", result$model_string),
               info = "Should not use grp.mm for unconstrained weights")
})

test_that("createModelstring() does NOT use accumulator when weights are pre-computed", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n, c = TRUE), RE = FALSE),
    family = "Weibull"
  )

  # When fn has no parameters (like 1/n), weights are pre-computed in R
  # So the accumulator pattern should NOT be in the JAGS model
  expect_false(grepl("cum\\.uw", result$model_string),
               info = "Should not have accumulator when weights are pre-computed")

  # Should indicate weights are pre-computed
  expect_true(grepl("pre-computed", result$model_string, ignore.case = TRUE),
              info = "Should indicate weights are pre-computed")
})

test_that("accumulator pattern produces balanced braces", {
  result <- setup_model(
    formula = Surv(dur_wkb, event_wkb) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b0 + b1 * pseat, c = TRUE), RE = TRUE),
    family = "Weibull"
  )

  # Count braces - should still be balanced with accumulator pattern
  open_braces <- length(gregexpr("\\{", result$model_string)[[1]])
  close_braces <- length(gregexpr("\\}", result$model_string)[[1]])
  expect_equal(open_braces, close_braces,
               info = "Accumulator pattern should maintain balanced braces")
})
