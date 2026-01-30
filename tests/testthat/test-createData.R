data("coalgov")

# ================================================================================================ #
# Tests for createData() - internal data transformation function
# ================================================================================================ #

# Helper to parse formula and create data
setup_data <- function(formula, family, data = coalgov) {
  formula_parts <- bml:::dissectFormula(formula, family, data)
  data_parts <- bml:::createData(data, formula_parts)
  list(
    data_parts = data_parts,
    formula_parts = formula_parts
  )
}

test_that("createData() returns correct structure", {
  result <- setup_data(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian"
  )

  expect_type(result$data_parts, "list")
  expect_true("data" %in% names(result$data_parts))
  expect_true("mm_blocks" %in% names(result$data_parts))
  expect_true("main" %in% names(result$data_parts))
  expect_true("hm_blocks" %in% names(result$data_parts))
})

test_that("createData() correctly processes simple Gaussian model", {
  result <- setup_data(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian"
  )

  # Main level should have lhs and vars
  expect_true("lhs" %in% names(result$data_parts$main))
  expect_true("vars" %in% names(result$data_parts$main))
  expect_true("dat" %in% names(result$data_parts$main))

  # Should include the outcome variable name
  expect_true("sim.y" %in% result$data_parts$main$lhs)
})

test_that("createData() correctly processes survival model", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Weibull"
  )

  # Survival models have time and event in lhs
  expect_true(length(result$data_parts$main$lhs) == 2)
  expect_true("govdur" %in% result$data_parts$main$lhs)
  expect_true("earlyterm" %in% result$data_parts$main$lhs)
})

test_that("createData() correctly handles mm() blocks", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull"
  )

  # Should have mm_blocks
  expect_false(is.null(result$data_parts$mm_blocks))
  expect_equal(length(result$data_parts$mm_blocks), 1)

  # mm block should have expected structure
  block <- result$data_parts$mm_blocks[[1]]
  expect_true("vars" %in% names(block))
  expect_true("dat" %in% names(block))
  expect_true("fn" %in% names(block))
  expect_true("RE" %in% names(block))
  expect_true("ar" %in% names(block))

  # Variables should be extracted

  expect_true("fdep" %in% block$vars)
})

test_that("createData() correctly handles hm() blocks", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull"
  )

  # Should have hm_blocks
  expect_false(is.null(result$data_parts$hm_blocks))
  expect_equal(length(result$data_parts$hm_blocks), 1)

  # hm block should have expected structure
  block <- result$data_parts$hm_blocks[[1]]
  expect_true("id" %in% names(block))
  expect_true("type" %in% names(block))
  expect_equal(block$type, "RE")
})

test_that("createData() correctly handles multiple mm() blocks", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE) +
      mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  # Should have 2 mm blocks
  expect_equal(length(result$data_parts$mm_blocks), 2)

  # Each block should have its own variables
  expect_true("fdep" %in% result$data_parts$mm_blocks[[1]]$vars)
  expect_true("ipd" %in% result$data_parts$mm_blocks[[2]]$vars)
})

test_that("createData() correctly stores weight function info", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull"
  )

  block <- result$data_parts$mm_blocks[[1]]

  # fn should be preserved
  expect_false(is.null(block$fn))
  expect_true("string" %in% names(block$fn))
  expect_true("vars" %in% names(block$fn))
})

test_that("createData() correctly handles parameterized weight functions", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel), RE = FALSE),
    family = "Weibull"
  )

  block <- result$data_parts$mm_blocks[[1]]

  # fn should have parameters
  expect_true("params" %in% names(block$fn))
  expect_true("b0" %in% block$fn$params || "b1" %in% block$fn$params)

  # fn should have variables
  expect_true("pseatrel" %in% block$fn$vars)
})

test_that("createData() correctly handles fix() variables at main level", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + fix(majority, 1.0),
    family = "Weibull"
  )

  # Should have fixed variables in main
  expect_true("vars_fixed" %in% names(result$data_parts$main))
})

test_that("createData() correctly handles fix() in mm() blocks", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fix(fdep, 1.0) + ipd), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull"
  )

  block <- result$data_parts$mm_blocks[[1]]

  # Should have fixed variables
  expect_true("vars_fixed" %in% names(block))
  expect_false(is.null(block$vars_fixed))

  # Free variable should be ipd
  expect_true("ipd" %in% block$vars)
})

test_that("createData() preserves sample size in main data", {
  result <- setup_data(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian"
  )

  # Main dat should have same number of rows as original data
  expect_equal(nrow(result$data_parts$main$dat), nrow(coalgov))
})

test_that("createData() creates sequential IDs", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  # Data should have sequential mmid and mainid
  expect_true("mmid" %in% names(result$data_parts$data))
  expect_true("mainid" %in% names(result$data_parts$data))

  # IDs should start at 1
  expect_equal(min(result$data_parts$data$mmid), 1)
  expect_equal(min(result$data_parts$data$mainid), 1)
})

test_that("createData() creates sequential hm IDs", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull"
  )

  # Data should have sequential hmid
  expect_true("hmid" %in% names(result$data_parts$data))
  expect_equal(min(result$data_parts$data$hmid), 1)
})

test_that("createData() handles intercept correctly", {
  # With intercept
  result_with <- setup_data(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian"
  )

  # Without intercept
  result_without <- setup_data(
    formula = sim.y ~ 0 + majority,
    family = "Gaussian"
  )

  # With intercept should have X0 in vars
  expect_true("X0" %in% result_with$data_parts$main$vars)

  # Without intercept should not have X0
  expect_false("X0" %in% result_without$data_parts$main$vars)
})

test_that("createData() handles AR specifications", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE, ar = TRUE),
    family = "Weibull"
  )

  block <- result$data_parts$mm_blocks[[1]]

  # AR flag should be preserved
  expect_true(block$ar)
})

test_that("createData() handles subset of data", {
  subset_data <- coalgov[1:100, ]

  result <- setup_data(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    data = subset_data
  )

  expect_equal(nrow(result$data_parts$main$dat), 100)
})

test_that("createData() handles Cox models", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Cox"
  )

  # Should have survival lhs variables
  expect_true("govdur" %in% result$data_parts$main$lhs)
  expect_true("earlyterm" %in% result$data_parts$main$lhs)
})

test_that("createData() stores mm_blocks attributes", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull"
  )

  # mm_blocks should have attributes
  expect_true(!is.null(attr(result$data_parts$mm_blocks, "has_RE")))
  expect_true(!is.null(attr(result$data_parts$mm_blocks, "has_vars")))

  # has_RE should be TRUE since we specified RE = TRUE
  expect_true(attr(result$data_parts$mm_blocks, "has_RE"))
})

test_that("createData() handles hm with fixed effects type", {
  result <- setup_data(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "FE"),
    family = "Weibull"
  )

  block <- result$data_parts$hm_blocks[[1]]
  expect_equal(block$type, "FE")
})
