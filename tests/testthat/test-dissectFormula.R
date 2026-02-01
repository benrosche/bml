data("coalgov")
data <- coalgov

test_that("dissectFormula() works with all model families and with and without intercept", {

  # Gaussian
  expect_no_error(
    dissectFormula(
      formula(
        event_wkb ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Gaussian",
      data
    ),
    message = "LHS CHECK (Gaussian)"
  )

  expect_error(
    dissectFormula(
      formula(
        Surv(dur_wkb, event_wkb) ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Gaussian",
      data
    ),
    label = "LHS CHECK (Gaussian) - Surv() should fail"
  )

  # Binomial
  expect_no_error(
    dissectFormula(
      formula(
        event_wkb ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Binomial",
      data
    ),
    message = "LHS CHECK (Binomial)"
  )

  expect_error(
    dissectFormula(
      formula(
        Surv(dur_wkb, event_wkb) ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Binomial",
      data
    ),
    label = "LHS CHECK (Binomial) - Surv() should fail"
  )

  # Weibull
  expect_no_error(
    dissectFormula(
      formula(
        Surv(dur_wkb, event_wkb) ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Weibull",
      data
    ),
    message = "LHS CHECK (Weibull)"
  )

  expect_error(
    dissectFormula(
      formula(
        event_wkb ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Weibull",
      data
    ),
    label = "LHS CHECK (Weibull) - Surv() should fail"
  )

  expect_no_error(
    dissectFormula(
      formula(
        event_wkb ~ majority +
          mwc +
          hm(id = id(cid), name = country, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Gaussian",
      data
    ),
    message = "INTERCEPT MISSING"
  )

})

test_that("dissectFormula() works with different ways of specifying mm()", {

  # RHS

  expect_error(
    dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(vars = vars(cohesion), fn = fn(w ~ 1/exp(cohesion), c = FALSE))),
      family="Weibull",
      data
    ),
    label="mm() -> id() missing"
  )

  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), fn = fn(w ~ 1/exp(cohesion), c = FALSE))),
      family="Weibull",
      data
    ),
    message="mm() -> vars missing should auto-set RE = TRUE"
  )

  expect_error(
    dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion))),
      family="Weibull",
      data
    ),
    label="mm() -> fn() missing"
  )

  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE), RE = TRUE)),
      family="Weibull",
      data
    ),
    message="mm() -> no vars (NULL) should be possible"
  )

  expect_true(
    {dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE), RE = TRUE)),
      family="Weibull",
      data)}$mm[[1]]$fn$constraint,
    label="mm() -> c = TRUE (constrained) should return TRUE"
  )

  expect_false(
    {dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = FALSE), RE = TRUE)),
      family="Weibull",
      data)}$mm[[1]]$fn$constraint,
    label="mm() -> c = FALSE (unconstrained) should return FALSE"
  )

  expect_true(
    {dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)),
      family="Weibull",
      data)}$mm[[1]]$fn$constraint,
    label="mm() -> default c (TRUE) should return TRUE"
  )

  expect_no_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(
            id = id(pid, gid),
            vars = vars(finance),
            fn = fn(w ~ 1 / (1 + (n - 1) * exp(-(b1 * pseat))), c = FALSE),
            RE = TRUE
          )
      ),
      family = "Weibull",
      data
    ),
    message = "mm() -> complex weight function should be possible"
  )

  expect_no_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(
            id = id(pid, gid),
            vars = vars(finance),
            fn = fn(w ~ b1 * pseat + (1-b1) * pseat, c = FALSE),
            RE = TRUE
          )
      ),
      family = "Weibull",
      data
    ),
    message = "mm() -> complex weight function should be possible"
  )

  expect_no_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(
            id = id(pid, gid),
            vars = vars(finance),
            fn = fn(w ~ b1*pseat^n, c = FALSE),
            RE = TRUE
          )
      ),
      family = "Weibull",
      data
    ),
    message = "mm() -> complex weight function should be possible"
  )

})

test_that("Main formula intercept handling follows lm() behavior", {
  # Default: intercept included (like lm)
  result_default <- dissectFormula(
    formula(event_wkb ~ majority),
    family = "Gaussian",
    data
  )
  expect_true(
    "X0" %in% result_default$mainvars,
    label = "Default formula should include intercept (X0)"
  )

  # Explicit 1+: intercept included

  result_explicit <- dissectFormula(
    formula(event_wkb ~ 1 + majority),
    family = "Gaussian",
    data
  )
  expect_true(
    "X0" %in% result_explicit$mainvars,
    label = "Explicit 1+ should include intercept (X0)"
  )

  # Using 0+: no intercept
  result_no_int_0 <- dissectFormula(
    formula(event_wkb ~ 0 + majority),
    family = "Gaussian",
    data
  )
  expect_false(
    "X0" %in% result_no_int_0$mainvars,
    label = "0+ should exclude intercept"
  )

  # Using -1: no intercept
  result_no_int_minus <- dissectFormula(
    formula(event_wkb ~ -1 + majority),
    family = "Gaussian",
    data
  )
  expect_false(
    "X0" %in% result_no_int_minus$mainvars,
    label = "-1 should exclude intercept"
  )
})

test_that("vars() in mm() ignores numeric literals", {

  # vars(1 + finance) should only extract finance, not 1
  result_with_1 <- dissectFormula(
    formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = vars(1 + finance), fn = fn(w ~ 1/n))),
    family = "Weibull",
    data
  )
  expect_equal(result_with_1$mm[[1]]$vars$free, "finance",
    label = "vars(1 + finance) should only extract finance")

  # vars(-1 + finance) should only extract finance, not -1
  result_with_minus1 <- dissectFormula(
    formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = vars(-1 + finance), fn = fn(w ~ 1/n))),
    family = "Weibull",
    data
  )
  expect_equal(result_with_minus1$mm[[1]]$vars$free, "finance",
    label = "vars(-1 + finance) should only extract finance")

  # vars(0 + finance) should only extract finance
  result_with_0 <- dissectFormula(
    formula(Surv(govdur, event_wkb) ~ 1 + mm(id = id(pid, gid), vars = vars(0 + finance), fn = fn(w ~ 1/n))),
    family = "Weibull",
    data
  )
  expect_equal(result_with_0$mm[[1]]$vars$free, "finance",
    label = "vars(0 + finance) should only extract finance")

})

test_that("vars() requires formula-style syntax with +", {

  # Comma-separated should throw an error
  expect_error(
    vars(finance, cohesion),
    "must be specified using formula-style",
    label = "vars(x, y) should throw error"
  )

  # Formula-style should work
  expect_no_error(
    vars(finance + cohesion),
    message = "vars(x + y) should work"
  )

  # Single variable should work
  expect_no_error(
    vars(finance),
    message = "vars(x) should work"
  )

})

test_that("dissectFormula() works with different ways of specifying hm()", {

  expect_error(
    dissectFormula(
      formula(event_wkb ~ 1 + hm(name = country, type = "RE", showFE = F)),
      family="Gaussian",
      data
    ),
    label="hm() -> id() missing"
  )

  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, event_wkb) ~ 1 + hm(id = id(cid), name = country)),
      family="Weibull",
      data
    ),
    message="hm() -> no type or showFE should be possible"
  )

})

test_that("mm() blocks with different mmid but same mainid are allowed", {
  # Create data with two different mmid columns
  data_multi_mmid <- coalgov %>%
    dplyr::mutate(pid2 = pid + 1000)  # Create second mmid column

  # Two mm blocks with different mmid (pid and pid2) but same mainid (gid)
  # This should be allowed
  expect_no_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE) +
          mm(id = id(pid2, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = FALSE)
      ),
      family = "Weibull",
      data_multi_mmid
    ),
    message = "mm() blocks with different mmid but same mainid should be allowed"
  )

  # Verify mm_groups is correctly populated
  result <- dissectFormula(
    formula(
      Surv(govdur, event_wkb) ~ 1 +
        mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE) +
        mm(id = id(pid2, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = FALSE)
    ),
    family = "Weibull",
    data_multi_mmid
  )

  expect_true(!is.null(result$mm_groups), label = "mm_groups should be populated")
  expect_equal(length(result$mm_groups), 2, label = "Should have 2 mmid groups")
  expect_true("pid" %in% names(result$mm_groups), label = "mm_groups should contain 'pid'")
  expect_true("pid2" %in% names(result$mm_groups), label = "mm_groups should contain 'pid2'")
})

test_that("mm() blocks with different mainid are not allowed", {
  # Create data with two different mainid columns
  data_multi_mainid <- coalgov %>%
    dplyr::mutate(gid2 = gid + 1000)  # Create second mainid column

  # Two mm blocks with different mainid should fail
  expect_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = FALSE) +
          mm(id = id(pid, gid2), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = FALSE)
      ),
      family = "Weibull",
      data_multi_mainid
    ),
    regexp = "same mainid",
    label = "mm() blocks with different mainid should fail"
  )
})

test_that("RE=TRUE for multiple mm() blocks with same mmid is not allowed", {
  # Two mm blocks with same mmid (pid) and both with RE=TRUE should fail
  expect_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(id = id(pid, gid), vars = vars(finance), fn = fn(w ~ 1/n), RE = TRUE) +
          mm(id = id(pid, gid), vars = vars(cohesion), fn = fn(w ~ 1/n), RE = TRUE)
      ),
      family = "Weibull",
      data
    ),
    regexp = "RE = TRUE can only be specified for one",
    label = "Multiple RE=TRUE with same mmid should fail"
  )
})

test_that("RE=TRUE for multiple mm() blocks with different mmid is allowed", {
  # Create data with two different mmid columns
  data_multi_mmid <- coalgov %>%
    dplyr::mutate(pid2 = pid + 1000)  # Create second mmid column

  # Two mm blocks with different mmid and both with RE=TRUE should be allowed
  expect_no_error(
    dissectFormula(
      formula(
        Surv(govdur, event_wkb) ~ 1 +
          mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE) +
          mm(id = id(pid2, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)
      ),
      family = "Weibull",
      data_multi_mmid
    ),
    message = "Multiple RE=TRUE with different mmid should be allowed"
  )
})
