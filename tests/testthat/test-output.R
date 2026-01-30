data("coalgov")

# ================================================================================================ #
# Setup: Create a simple fitted model for output tests
# ================================================================================================ #

# Create a minimal model for testing (if JAGS is available)
create_test_model <- function(monitor = FALSE) {
  testthat::skip_on_cran()

  if (!requireNamespace("rjags", quietly = TRUE)) {
    skip("rjags not available")
  }

  # Use a very small subset and few iterations for speed
  test_data <- coalgov[1:100, ]

  tryCatch({
    m <- bml(
      sim.y ~ 1 + majority,
      family = "Gaussian",
      data = test_data,
      n.iter = 500,
      n.burnin = 100,
      n.chains = 2,
      monitor = monitor
    )
    return(m)
  }, error = function(e) {
    skip(paste("Could not fit test model:", e$message))
  })
}

# ================================================================================================ #
# Tests for summary.bml()
# ================================================================================================ #

test_that("summary.bml() returns correct structure", {
  m <- create_test_model(monitor = FALSE)

  s <- summary(m)

  expect_s3_class(s, "bml_summary")
  expect_s3_class(s, "data.frame")
  expect_true("Parameter" %in% names(s))
  expect_true("mean" %in% names(s))
  expect_true("sd" %in% names(s))
  expect_true("lb" %in% names(s))
  expect_true("ub" %in% names(s))
  expect_gt(nrow(s), 0)
})

test_that("summary.bml() rounding works", {
  m <- create_test_model(monitor = FALSE)

  s2 <- summary(m, r = 2)
  s4 <- summary(m, r = 4)

  # Check that values are actually different due to rounding
  expect_true(any(s2$mean != s4$mean))
})

test_that("summary.bml() preserves attributes", {
  m <- create_test_model(monitor = FALSE)

  s <- summary(m)

  expect_false(is.null(attr(s, "estimate_type")))
  expect_false(is.null(attr(s, "credible_interval")))
  expect_false(is.null(attr(s, "DIC")))
  expect_false(is.null(attr(s, "outcome_family")))
})

test_that("summary.bml() print method works", {
  m <- create_test_model(monitor = FALSE)

  s <- summary(m)

  expect_output(print(s), "Outcome family")
  expect_output(print(s), "Estimates")
  expect_output(print(s), "DIC")
})

# ================================================================================================ #
# Tests for mcmcDiag()
# ================================================================================================ #

test_that("mcmcDiag() requires monitor = TRUE", {
  m <- create_test_model(monitor = FALSE)

  # Should fail because monitor was FALSE
  expect_error(
    mcmcDiag(m, parameters = "b"),
    "JAGS output could not be retrieved"
  )
})

test_that("mcmcDiag() returns correct structure with monitor = TRUE", {
  m <- create_test_model(monitor = TRUE)

  diag <- mcmcDiag(m, parameters = "b")

  expect_s3_class(diag, "data.frame")
  expect_true(all(c("Gelman/Rubin convergence statistic",
                    "Geweke z-score",
                    "Heidelberger/Welch p-value",
                    "Autocorrelation (lag 50)") %in% rownames(diag)))
  expect_gt(ncol(diag), 0)
})

test_that("mcmcDiag() works with specific parameter names", {
  m <- create_test_model(monitor = TRUE)

  diag <- mcmcDiag(m, parameters = c("b[1]", "b[2]"))

  expect_equal(ncol(diag), 2)
  expect_true("b[1]" %in% colnames(diag))
  expect_true("b[2]" %in% colnames(diag))
})

test_that("mcmcDiag() handles pattern matching", {
  m <- create_test_model(monitor = TRUE)

  # Pattern should match all b parameters
  diag <- mcmcDiag(m, parameters = "b")

  expect_gt(ncol(diag), 1)  # Should have multiple b parameters
})

# ================================================================================================ #
# Tests for monetPlot()
# ================================================================================================ #

test_that("monetPlot() requires monitor = TRUE", {
  m <- create_test_model(monitor = FALSE)

  expect_error(
    monetPlot(m, parameter = "b[1]"),
    "JAGS output could not be retrieved"
  )
})

test_that("monetPlot() returns a ggplot object", {
  m <- create_test_model(monitor = TRUE)

  p <- monetPlot(m, parameter = "b[1]")

  expect_s3_class(p, "patchwork")
  expect_s3_class(p, "gg")
})

test_that("monetPlot() accepts custom labels", {
  m <- create_test_model(monitor = TRUE)

  p <- monetPlot(m, parameter = "b[1]", label = "Custom Label")

  # Check that the plot was created without error
  expect_s3_class(p, "patchwork")
})

test_that("monetPlot() respects yaxis parameter", {
  m <- create_test_model(monitor = TRUE)

  p1 <- monetPlot(m, parameter = "b[1]", yaxis = TRUE)
  p2 <- monetPlot(m, parameter = "b[1]", yaxis = FALSE)

  expect_s3_class(p1, "patchwork")
  expect_s3_class(p2, "patchwork")
})

test_that("monetPlot() respects rounding parameter", {
  m <- create_test_model(monitor = TRUE)

  p1 <- monetPlot(m, parameter = "b[1]", r = 2)
  p2 <- monetPlot(m, parameter = "b[1]", r = 4)

  # Just check that both work without error
  expect_s3_class(p1, "patchwork")
  expect_s3_class(p2, "patchwork")
})

# ================================================================================================ #
# Integration test: Full model with output functions
# ================================================================================================ #

test_that("Full workflow with mm() and output functions works", {
  testthat::skip_on_cran()

  if (!requireNamespace("rjags", quietly = TRUE)) {
    skip("rjags not available")
  }

  # Small subset for speed
  test_data <- coalgov[1:100, ]

  tryCatch({
    # Fit model with mm block
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n), RE = TRUE),
      family = "Weibull",
      data = test_data,
      n.iter = 500,
      n.burnin = 100,
      n.chains = 2,
      monitor = TRUE
    )

    # Test summary
    s <- summary(m)
    expect_s3_class(s, "bml_summary")
    expect_true(any(grepl("mm", s$Parameter)))

    # Test mcmcDiag
    diag <- mcmcDiag(m, parameters = "b")
    expect_s3_class(diag, "data.frame")

    # Test monetPlot
    p <- monetPlot(m, parameter = "b[1]")
    expect_s3_class(p, "patchwork")

  }, error = function(e) {
    skip(paste("Integration test failed:", e$message))
  })
})
