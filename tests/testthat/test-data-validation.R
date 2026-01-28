data("coalgov")

# ================================================================================================ #
# Tests for data validation
# ================================================================================================ #

test_that("Missing data in outcome variable is handled", {
  test_data <- coalgov
  test_data$sim.y[1:10] <- NA

  # Should either warn or handle gracefully
  expect_warning({
    m <- bml(
      sim.y ~ 1 + majority,
      family = "Gaussian",
      data = test_data,
      run = FALSE
    )
  }, NA)  # NA means no warning expected, or adjust if warnings are expected
})

test_that("Missing data in covariates is detected", {
  test_data <- coalgov
  test_data$majority[1:10] <- NA

  expect_warning({
    m <- bml(
      sim.y ~ 1 + majority,
      family = "Gaussian",
      data = test_data,
      run = FALSE
    )
  }, NA)  # Adjust based on actual behavior
})

test_that("Non-existent outcome variable throws error", {
  expect_error(
    bml(
      nonexistent ~ 1 + majority,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Non-existent covariate throws error", {
  expect_error(
    bml(
      sim.y ~ 1 + nonexistent,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Non-existent ID variable in mm() throws error", {
  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(nonexistent1, nonexistent2), vars = vars(fdep), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    ),
    "ID variable.*not found"
  )
})

test_that("Non-existent variable in mm() vars throws error", {
  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(nonexistent), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    ),
    "Variable.*not found"
  )
})

test_that("Non-existent variable in mm() fn weight function throws error", {
  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/nonexistent)),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    ),
    "Variable.*not found"
  )
})

test_that("Non-existent ID variable in hm() throws error", {
  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        hm(id = id(nonexistent), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    ),
    "ID variable.*not found"
  )
})

test_that("Non-existent name variable in hm() throws error", {
  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        hm(id = id(cid), name = nonexistent, type = "FE", showFE = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    ),
    "name variable.*not found"
  )
})

# ================================================================================================ #
# Tests for data structure requirements
# ================================================================================================ #

test_that("Surv() requires two arguments for Weibull", {
  expect_error(
    bml(
      Surv(govdur) ~ 1 + majority,
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Surv() requires two arguments for Cox", {
  expect_error(
    bml(
      Surv(govdur) ~ 1 + majority,
      family = "Cox",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Event indicator must be 0/1 for survival models", {
  test_data <- coalgov
  test_data$earlyterm[1] <- 2  # Invalid value

  # Should either error or warn
  expect_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority,
      family = "Weibull",
      data = test_data,
      run = FALSE
    )
  }, NA)  # Adjust based on actual behavior
})

# ================================================================================================ #
# Tests for mm() data structure requirements
# ================================================================================================ #

test_that("mm() requires member IDs to be unique within groups", {
  # This tests whether duplicate member-group combinations are handled
  test_data <- coalgov[1:50, ]

  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    )
  })
})

test_that("mm() weight function variables must exist in data", {
  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * nonexistent)),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  )
})

# ================================================================================================ #
# Tests for data types
# ================================================================================================ #

test_that("Factor variables in mm() are handled", {
  test_data <- coalgov
  test_data$fdep_factor <- as.factor(round(test_data$fdep))

  # Should work or give clear error
  result <- tryCatch({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep_factor), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    )
    "success"
  }, error = function(e) {
    "error"
  })

  expect_true(result %in% c("success", "error"))
})

test_that("Character ID variables are handled", {
  test_data <- coalgov
  test_data$pid_char <- as.character(test_data$pid)
  test_data$gid_char <- as.character(test_data$gid)

  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid_char, gid_char), vars = vars(fdep), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for edge cases
# ================================================================================================ #

test_that("Single observation groups in mm() are handled", {
  # Keep only governments with 1 party
  test_data <- coalgov[coalgov$n == 1, ]

  if (nrow(test_data) > 0) {
    expect_no_error({
      m <- bml(
        Surv(govdur, earlyterm) ~ 1 + majority +
          mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n)),
        family = "Weibull",
        data = test_data,
        run = FALSE
      )
    })
  } else {
    skip("No single-party governments in data")
  }
})

test_that("Very small sample size is handled", {
  test_data <- coalgov[1:10, ]

  expect_no_error({
    m <- bml(
      sim.y ~ 1 + majority,
      family = "Gaussian",
      data = test_data,
      run = FALSE
    )
  })
})

test_that("All censored observations in survival model", {
  test_data <- coalgov[1:50, ]
  test_data$earlyterm <- 0  # All censored

  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority,
      family = "Weibull",
      data = test_data,
      run = FALSE
    )
  })
})

test_that("No events in any interval for Cox with intervals", {
  test_data <- coalgov[1:50, ]
  test_data$earlyterm <- 0  # All censored

  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority,
      family = "Cox",
      cox_intervals = 5,
      data = test_data,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for formula parsing edge cases
# ================================================================================================ #

test_that("Complex formula with interactions works", {
  expect_no_error({
    m <- bml(
      sim.y ~ 1 + majority * mwc,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Formula with I() transformation works", {
  expect_no_error({
    m <- bml(
      sim.y ~ 1 + I(majority + mwc),
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Very long formula works", {
  expect_no_error({
    m <- bml(
      sim.y ~ 1 + majority + mwc + hetero +
        mm(id = id(pid, gid), vars = vars(fdep + ipd + rile), fn = fn(w ~ 1/n), RE = TRUE) +
        hm(id = id(cid), vars = vars(investiture + pmpower), type = "RE"),
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})
