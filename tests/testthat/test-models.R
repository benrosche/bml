data("coalgov")

# ================================================================================================ #
# Tests for different model families
# ================================================================================================ #

test_that("Gaussian model setup works", {
  expect_no_error({
    m <- bml(
      event_wkb ~ 1 + majority,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Binomial model setup works", {
  expect_no_error({
    m <- bml(
      event_wkb ~ 1 + majority,
      family = "Binomial",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Weibull model setup works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority,
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Cox model setup works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority,
      family = "Cox",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Cox model with intervals works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority,
      family = "Cox",
      cox_intervals = 10,
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for mm() blocks
# ================================================================================================ #

test_that("mm() with variables and equal weights works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("mm() with variables and RE works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("mm() with RE only (no variables) works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("mm() with multiple variables works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance + cohesion), fn = fn(w ~ 1/n), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("mm() with parameterized weight function works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 +
        majority +
        mm(
          id = id(pid, gid),
          vars = vars(finance),
          fn = fn(w ~ 1 / (1 + (n - 1) * exp(-(b1 * cohesion)))),
          RE = FALSE
        ),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("mm() with aggregation functions in weight function works", {
  # Create bounded version of pseat for valid weights (0 to 1)
  test_data <- coalgov
  test_data$pseatrel01 <- (test_data$pseat - min(test_data$pseat, na.rm = TRUE)) /
    (max(test_data$pseat, na.rm = TRUE) - min(test_data$pseat, na.rm = TRUE))

  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 +
        majority +
        mm(
          id = id(pid, gid),
          vars = vars(finance),
          fn = fn(w ~ b1 * min(pseatrel01) + (1 - b1) * mean(pseatrel01)),
          RE = FALSE
        ),
      family = "Weibull",
      data = test_data,
      run = FALSE
    )
  })

  # When fn has parameters, X.w.1 contains the aggregated columns
  m <- bml(
    event_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b1 * min(pseatrel01) + (1 - b1) * max(pseatrel01)), RE = FALSE),
    family = "Gaussian",
    data = test_data,
    run = FALSE
  )
  expect_true("pseatrel01_min" %in% colnames(m$jags.data$X.w.1))
  expect_true("pseatrel01_max" %in% colnames(m$jags.data$X.w.1))
})

test_that("mm() with quantile aggregation in weight function works", {
  # Create bounded version of pseat for valid weights (0 to 1)
  test_data <- coalgov
  test_data$pseatrel01 <- (test_data$pseat - min(test_data$pseat, na.rm = TRUE)) /
    (max(test_data$pseat, na.rm = TRUE) - min(test_data$pseat, na.rm = TRUE))

  expect_no_error({
    m <- bml(
      event_wkb ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b1 * quantile(pseatrel01, 0.75) / max(pseatrel01)), RE = FALSE),
      family = "Gaussian",
      data = test_data,
      run = FALSE
    )
  })

  # Verify quantile column naming (with parameter so X.w.1 exists)
  m <- bml(
    event_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ b1 * quantile(pseatrel01, 0.25)), RE = FALSE),
    family = "Gaussian",
    data = test_data,
    run = FALSE
  )
  expect_true("pseatrel01_q25" %in% colnames(m$jags.data$X.w.1))
})

test_that("mm() with unconstrained weights works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n, c = FALSE), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Multiple mm() blocks work", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE) +
        mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for hm() blocks
# ================================================================================================ #

test_that("hm() with random effects works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        hm(id = id(cid), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("hm() with fixed effects works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        hm(id = id(cid), type = "FE", showFE = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("hm() with variables works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        hm(id = id(cid), vars = vars(investiture), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("hm() with name specification works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        hm(id = id(cid), name = country, type = "FE", showFE = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Multiple hm() blocks work", {
  expect_no_error({
    m <- bml(
      event_wkb ~ 1 + majority +
        hm(id = id(cid), type = "RE") +
        hm(id = id(gid), type = "RE"),
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for combined mm() and hm() blocks
# ================================================================================================ #

test_that("mm() and hm() blocks work together", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE) +
        hm(id = id(cid), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("Complex model with multiple blocks works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority + mwc +
        mm(id = id(pid, gid), vars = vars(finance + cohesion), fn = fn(w ~ 1/n), RE = TRUE) +
        hm(id = id(cid), vars = vars(investiture), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for fix() functionality
# ================================================================================================ #

test_that("fix() with main-level variables works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + fix(majority, 1.0),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("fix() within mm() works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = vars(fix(finance, 1.0) + cohesion), fn = fn(w ~ 1/n), RE = FALSE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("fix() within hm() works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        hm(id = id(cid), vars = vars(fix(investiture, 0.5)), type = "RE"),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for AR specifications
# ================================================================================================ #

test_that("AR in mm() works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE, ar = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

test_that("AR in hm() works", {
  expect_no_error({
    m <- bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        hm(id = id(cid), type = "RE", ar = TRUE),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for model without intercept
# ================================================================================================ #

test_that("Model without intercept works", {
  expect_no_error({
    m <- bml(
      event_wkb ~ 0 + majority,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Tests for intercept-only model
# ================================================================================================ #

test_that("Intercept-only model works", {
  expect_no_error({
    m <- bml(
      event_wkb ~ 1,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  })
})

# ================================================================================================ #
# Error handling tests
# ================================================================================================ #

test_that("Invalid family throws error", {
  expect_error(
    bml(
      event_wkb ~ 1 + majority,
      family = "InvalidFamily",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Surv() with non-survival family throws error", {
  expect_error(
    bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority,
      family = "Gaussian",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Non-Surv() with survival family throws error", {
  expect_error(
    bml(
      dur_wkb ~ 1 + majority,
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Missing ID variables throw error", {
  expect_error(
    bml(
      Surv(dur_wkb, event_wkb) ~ 1 + majority +
        mm(id = id(nonexistent1, nonexistent2), vars = vars(finance), fn = fn(w ~ 1/n)),
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  )
})

test_that("Missing covariate variables throw error", {
  expect_error(
    bml(
      Surv(dur_wkb, event_wkb) ~ 1 + nonexistent_var,
      family = "Weibull",
      data = coalgov,
      run = FALSE
    )
  )
})
