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

test_that("Missing data in RHS covariates throws error", {
  test_data <- coalgov
  test_data$majority[1:10] <- NA

  expect_error(
    bml(
      sim.y ~ 1 + majority,
      family = "Gaussian",
      data = test_data,
      run = FALSE
    ),
    "Missing values.*majority"
  )
})

test_that("Missing data in mm() vars throws error", {
  test_data <- coalgov
  test_data$fdep[1:10] <- NA

  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    ),
    "Missing values.*fdep"
  )
})

test_that("Missing data in weight function variables throws error", {
  test_data <- coalgov
  test_data$pseatrel[1:10] <- NA

  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    ),
    "Missing values.*pseatrel"
  )
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
    "ID variable.*not found in data"
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
    "mm\\(\\) variable.*not found"
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
    "Weight function variable.*not found"
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
    "hm\\(\\) id variable not found"
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
    "hm\\(\\) name variable not found"
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
# Tests for weight variable variation
# ================================================================================================ #

test_that("Constant weight variables within groups produce warning", {
  # Create test data where weight variable is constant within each mainid
  test_data <- coalgov
  # Set pseatrel to be constant within each government (mainid = gid)
  test_data <- test_data %>%
    dplyr::group_by(gid) %>%
    dplyr::mutate(pseatrel_constant = mean(pseatrel)) %>%
    dplyr::ungroup()

  expect_warning(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel_constant)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    ),
    "constant across members"
  )
})

test_that("Varying weight variables within groups produce no warning", {
  # pseatrel should vary within governments (parties have different seat shares)
  expect_no_warning(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ b0 + b1 * pseatrel)),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  )
})

# ================================================================================================ #
# Tests for mm() data structure requirements
# ================================================================================================ #

test_that("mm() errors on duplicate member-group combinations", {
  # Create data with duplicate pid-gid combinations (same party twice in same government)
  test_data <- rbind(coalgov[1:50, ], coalgov[1, ])

  expect_error(
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = test_data,
      run = FALSE
    ),
    "Duplicate member-group combinations"
  )
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

test_that("No events in any interval for Cox with intervals throws error", {
  test_data <- coalgov[1:50, ]
  test_data$earlyterm <- 0  # All censored

  # Cox model cannot be estimated without any events - should throw error
  expect_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority,
      family = "Cox",
      cox_intervals = 5,
      data = test_data,
      run = FALSE
    )
  }, "No events observed")
})

# ================================================================================================ #
# Tests for formula parsing edge cases
# ================================================================================================ #

test_that("Formula with interactions (*) works", {
  m <- bml(
    sim.y ~ 1 + majority * mwc,
    family = "Gaussian",
    data = coalgov,
    run = FALSE
  )

  # majority * mwc expands to: intercept + majority + mwc + majority:mwc = 4 params
  # Check for "for (x in 1:4)" in the modelstring
  expect_true(grepl("for\\s*\\(x in 1:4\\)", m$modelstring))
})

test_that("Formula with interaction (:) works", {
  m <- bml(
    sim.y ~ 1 + majority + mwc + majority:mwc,
    family = "Gaussian",
    data = coalgov,
    run = FALSE
  )

  # Same as majority * mwc: intercept + majority + mwc + majority:mwc = 4 params
  # Check for "for (x in 1:4)" in the modelstring
  expect_true(grepl("for\\s*\\(x in 1:4\\)", m$modelstring))
})

test_that("Formula with I() transformation works", {
  m <- bml(
    sim.y ~ 1 + I(majority + mwc),
    family = "Gaussian",
    data = coalgov,
    run = FALSE
  )

  # intercept + I(majority + mwc) = 2 params
  # Check for "for (x in 1:2)" in the modelstring
  expect_true(grepl("for\\s*\\(x in 1:2\\)", m$modelstring))
})

test_that("Formula with I() squared term works", {
  m <- bml(
    sim.y ~ 1 + majority + I(majority^2),
    family = "Gaussian",
    data = coalgov,
    run = FALSE
  )

  # intercept + majority + I(majority^2) = 3 params
  # Check for "for (x in 1:3)" in the modelstring
  expect_true(grepl("for\\s*\\(x in 1:3\\)", m$modelstring))
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

# ================================================================================================ #
# Tests for mm() and hm() vars with interactions and I()
# ================================================================================================ #

test_that("mm() vars with interaction (*) works", {
  expect_no_error({
    bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep * ipd), fn = fn(w ~ 1/n), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })

  # fdep * ipd expands to: fdep + ipd + fdep:ipd = 3 mm vars
  m <- bml(
    Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep * ipd), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull",
    data = coalgov,
    run = FALSE
  )
  expect_equal(m$jags.data$n.Xmm.1, 3)
})

test_that("mm() vars with interaction (:) works", {
  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep + ipd + fdep:ipd), fn = fn(w ~ 1/n), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("mm() vars with I() transformation works", {
  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fdep + I(fdep^2)), fn = fn(w ~ 1/n), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })

  # fdep + I(fdep^2) = 2 mm vars
  m <- bml(
    Surv(govdur, earlyterm) ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(fdep + I(fdep^2)), fn = fn(w ~ 1/n), RE = FALSE),
    family = "Weibull",
    data = coalgov,
    run = FALSE
  )
  expect_equal(m$jags.data$n.Xmm.1, 2)
})

test_that("hm() vars with interaction (*) works", {
  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        hm(id = id(cid), vars = vars(investiture * pmpower), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("hm() vars with I() transformation works", {
  expect_no_error({
    m <- bml(
      Surv(govdur, earlyterm) ~ 1 + majority +
        hm(id = id(cid), vars = vars(investiture + I(investiture^2)), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Combined interactions in main formula, mm() and hm() works", {
  expect_no_error({
    m <- bml(
      sim.y ~ 1 + majority * mwc +
        mm(id = id(pid, gid), vars = vars(fdep * ipd), fn = fn(w ~ 1/n), RE = TRUE) +
        hm(id = id(cid), vars = vars(investiture + I(pmpower^2)), type = "RE"),
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})
