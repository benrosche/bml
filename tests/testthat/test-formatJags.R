data("coalgov")

# ================================================================================================ #
# Tests for formatJags() - formats JAGS output for display
# ================================================================================================ #

# Helper function to create mock JAGS output
create_mock_jags_output <- function(params) {
  # Create mock BUGS output structure
  bugs_output <- list()

  for (param in names(params)) {
    bugs_output[[param]] <- array(
      rnorm(1000, mean = params[[param]], sd = 0.1),
      dim = c(250, 2, 2)  # 250 iterations, 2 chains, 2 parameters (for indexed)
    )
  }

  bugs_output$mean <- lapply(params, function(x) x)
  bugs_output$sd <- lapply(params, function(x) 0.1)
  bugs_output$`2.5%` <- lapply(params, function(x) x - 0.2)
  bugs_output$`97.5%` <- lapply(params, function(x) x + 0.2)
  bugs_output$DIC <- 1000

  list(
    BUGSoutput = bugs_output,
    model = "mock model"
  )
}

test_that("formatJags() creates correct output structure", {
  # Create minimal parsed structure
  parsed <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = NULL,
    l2 = NULL
  )

  # Mock JAGS output
  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  expect_type(result, "list")
  expect_true("reg.table" %in% names(result))
  expect_true("jags.out" %in% names(result))
  expect_true("modelstring" %in% names(result))
})

test_that("formatJags() creates regression table with correct columns", {
  parsed <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = NULL,
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  reg_table <- result$reg.table

  expect_true("Parameter" %in% names(reg_table))
  expect_true("mean" %in% names(reg_table))
  expect_true("sd" %in% names(reg_table))
  expect_true("lb" %in% names(reg_table))
  expect_true("ub" %in% names(reg_table))
})

test_that("formatJags() correctly labels parameters", {
  parsed <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = NULL,
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  reg_table <- result$reg.table

  # Should have intercept and covariate labels
  expect_true("Intercept" %in% reg_table$Parameter || "b[1]" %in% reg_table$Parameter)
  expect_true("majority" %in% reg_table$Parameter || "b[2]" %in% reg_table$Parameter)
})

test_that("formatJags() includes family-specific parameters", {
  # Weibull should have shape parameter
  parsed_weibull <- list(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Weibull",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = NULL,
    l2 = NULL
  )

  jags_out_weibull <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "shape" = 1.5
  ))

  result_weibull <- bml:::formatJags(
    jags_out = jags_out_weibull,
    parsed = parsed_weibull,
    data = coalgov,
    family = "Weibull",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  expect_true("shape" %in% result_weibull$reg.table$Parameter)
})

test_that("formatJags() includes mm() block parameters", {
  parsed <- list(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
    family = "Weibull",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = list(
      list(
        id = c("pid", "gid"),
        vars = c("fdep"),
        fn = list(vars = "n", params = character(0)),
        RE = TRUE
      )
    ),
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "b.mm.1[1]" = -0.3,
    "sigma.mm1" = 0.5,
    "shape" = 1.5
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Weibull",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  reg_table <- result$reg.table

  # Should have mm-specific parameters
  expect_true(any(grepl("mm", reg_table$Parameter, ignore.case = TRUE)))
})

test_that("formatJags() includes hm() block parameters", {
  parsed <- list(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      hm(id = id(cid), type = "RE"),
    family = "Weibull",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = NULL,
    l2 = list(
      list(
        id = "cid",
        vars = NULL,
        type = "RE"
      )
    )
  )

  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma.hm1" = 0.3,
    "shape" = 1.5
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Weibull",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  reg_table <- result$reg.table

  # Should have hm-specific parameters
  expect_true(any(grepl("hm", reg_table$Parameter, ignore.case = TRUE)))
})

test_that("formatJags() adds metadata attributes", {
  parsed <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(
      vars = c("majority"),
      fixed = NULL,
      intercept = TRUE
    ),
    l1 = NULL,
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "sigma" = 0.8
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  reg_table <- result$reg.table

  # Check for metadata attributes
  expect_false(is.null(attr(reg_table, "estimate_type")))
  expect_false(is.null(attr(reg_table, "credible_interval")))
  expect_false(is.null(attr(reg_table, "DIC")))
  expect_false(is.null(attr(reg_table, "outcome_family")))
})

test_that("formatJags() correctly reports outcome family", {
  # Gaussian
  parsed_gaussian <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(vars = c("majority"), fixed = NULL, intercept = TRUE),
    l1 = NULL,
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list("b[1]" = 0.5, "b[2]" = 1.2, "sigma" = 0.8))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed_gaussian,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  expect_true(grepl("Gaussian", attr(result$reg.table, "outcome_family")))
})

test_that("formatJags() reports Cox intervals correctly", {
  parsed_cox <- list(
    formula = Surv(govdur, earlyterm) ~ 1 + majority,
    family = "Cox",
    l0 = list(vars = c("majority"), fixed = NULL, intercept = TRUE),
    l1 = NULL,
    l2 = NULL,
    lhs = c("govdur", "earlyterm")
  )

  jags_out <- create_mock_jags_output(list("b[1]" = 0.5, "b[2]" = 1.2))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed_cox,
    data = coalgov,
    family = "Cox",
    cox_intervals = 10,
    modelstring = "model { }"
  )

  # Should mention intervals in outcome description
  expect_true(grepl("10.*interval", attr(result$reg.table, "outcome_family"), ignore.case = TRUE))
})

test_that("formatJags() includes level specification", {
  parsed <- list(
    formula = Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE) +
      hm(id = id(cid), type = "RE"),
    family = "Weibull",
    l0 = list(vars = c("majority"), fixed = NULL, intercept = TRUE),
    l1 = list(
      list(
        id = c("pid", "gid"),
        vars = c("fdep"),
        fn = list(vars = "n", params = character(0)),
        RE = TRUE
      )
    ),
    l2 = list(
      list(
        id = "cid",
        vars = NULL,
        type = "RE"
      )
    )
  )

  jags_out <- create_mock_jags_output(list(
    "b[1]" = 0.5,
    "b[2]" = 1.2,
    "b.mm.1[1]" = -0.3,
    "sigma.mm1" = 0.5,
    "sigma.hm1" = 0.3,
    "shape" = 1.5
  ))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Weibull",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  expect_false(is.null(attr(result$reg.table, "level_spec")))
  expect_true(grepl("mm", attr(result$reg.table, "level_spec"), ignore.case = TRUE))
  expect_true(grepl("hm", attr(result$reg.table, "level_spec"), ignore.case = TRUE))
})

test_that("formatJags() stores modelstring", {
  parsed <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(vars = c("majority"), fixed = NULL, intercept = TRUE),
    l1 = NULL,
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list("b[1]" = 0.5, "b[2]" = 1.2, "sigma" = 0.8))

  test_modelstring <- "model { Y[j] ~ dnorm(mu[j], tau) }"

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = test_modelstring
  )

  expect_equal(result$modelstring, test_modelstring)
})

test_that("formatJags() stores JAGS output object", {
  parsed <- list(
    formula = sim.y ~ 1 + majority,
    family = "Gaussian",
    l0 = list(vars = c("majority"), fixed = NULL, intercept = TRUE),
    l1 = NULL,
    l2 = NULL
  )

  jags_out <- create_mock_jags_output(list("b[1]" = 0.5, "b[2]" = 1.2, "sigma" = 0.8))

  result <- bml:::formatJags(
    jags_out = jags_out,
    parsed = parsed,
    data = coalgov,
    family = "Gaussian",
    cox_intervals = NULL,
    modelstring = "model { }"
  )

  expect_false(is.null(result$jags.out))
  expect_true("BUGSoutput" %in% names(result$jags.out))
})
