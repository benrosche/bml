data("coalgov")
data <- coalgov

test_that("dissectFormula() works with all model families and with and without intercept", {

  # Gaussian
  expect_no_error(
    dissectFormula(
      formula(
        sim.y ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Gaussian",
      data
    ),
    message = "LHS CHECK (Gaussian)"
  )

  expect_error(
    dissectFormula(
      formula(
        Surv(sim.st, sim.e) ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
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
        earlyterm ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Binomial",
      data
    ),
    message = "LHS CHECK (Binomial)"
  )

  expect_error(
    dissectFormula(
      formula(
        Surv(sim.st, sim.e) ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
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
        Surv(sim.st, sim.e) ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Weibull",
      data
    ),
    message = "LHS CHECK (Weibull)"
  )

  expect_error(
    dissectFormula(
      formula(
        sim.y ~ 1 +
          majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
      ),
      family = "Weibull",
      data
    ),
    label = "LHS CHECK (Weibull) - Surv() should fail"
  )

  expect_no_error(
    dissectFormula(
      formula(
        sim.y ~ majority +
          mwc +
          hm(id = id(cid), name = cname, type = "RE", showFE = F) +
          mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1 / n, c = TRUE))
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
      formula(Surv(govdur, earlyterm) ~ 1 + mm(vars = vars(ipd), fn = fn(w ~ 1/exp(ipd), c = FALSE))),
      family="Weibull",
      data
    ),
    label="mm() -> id() missing"
  )

  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), fn = fn(w ~ 1/exp(ipd), c = FALSE))),
      family="Weibull",
      data
    ),
    message="mm() -> vars missing should auto-set RE = TRUE"
  )

  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = vars(ipd))),
      family="Weibull",
      data
    ),
    label="mm() -> fn() missing"
  )

  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE), RE = TRUE)),
      family="Weibull",
      data
    ),
    message="mm() -> no vars (NULL) should be possible"
  )

  expect_true(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE), RE = TRUE)),
      family="Weibull",
      data)}$mm[[1]]$fn$constraint,
    label="mm() -> c = TRUE (constrained) should return TRUE"
  )

  expect_false(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = FALSE), RE = TRUE)),
      family="Weibull",
      data)}$mm[[1]]$fn$constraint,
    label="mm() -> c = FALSE (unconstrained) should return FALSE"
  )

  expect_true(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)),
      family="Weibull",
      data)}$mm[[1]]$fn$constraint,
    label="mm() -> default c (TRUE) should return TRUE"
  )

  expect_no_error(
    dissectFormula(
      formula(
        Surv(govdur, earlyterm) ~ 1 +
          mm(
            id = id(pid, gid),
            vars = vars(fdep),
            fn = fn(w ~ 1 / (1 + (n - 1) * exp(-(b1 * pseatrel))), c = FALSE),
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
        Surv(govdur, earlyterm) ~ 1 +
          mm(
            id = id(pid, gid),
            vars = vars(fdep),
            fn = fn(w ~ b1 * pseatrel + (1-b1) * pseatrel, c = FALSE),
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
        Surv(govdur, earlyterm) ~ 1 +
          mm(
            id = id(pid, gid),
            vars = vars(fdep),
            fn = fn(w ~ b1*pseatrel^n, c = FALSE),
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
    formula(sim.y ~ majority),
    family = "Gaussian",
    data
  )
  expect_true(
    "X0" %in% result_default$mainvars,
    label = "Default formula should include intercept (X0)"
  )

  # Explicit 1+: intercept included

  result_explicit <- dissectFormula(
    formula(sim.y ~ 1 + majority),
    family = "Gaussian",
    data
  )
  expect_true(
    "X0" %in% result_explicit$mainvars,
    label = "Explicit 1+ should include intercept (X0)"
  )

  # Using 0+: no intercept
  result_no_int_0 <- dissectFormula(
    formula(sim.y ~ 0 + majority),
    family = "Gaussian",
    data
  )
  expect_false(
    "X0" %in% result_no_int_0$mainvars,
    label = "0+ should exclude intercept"
  )

  # Using -1: no intercept
  result_no_int_minus <- dissectFormula(
    formula(sim.y ~ -1 + majority),
    family = "Gaussian",
    data
  )
  expect_false(
    "X0" %in% result_no_int_minus$mainvars,
    label = "-1 should exclude intercept"
  )
})

test_that("vars() in mm() ignores numeric literals", {

  # vars(1 + fdep) should only extract fdep, not 1
  result_with_1 <- dissectFormula(
    formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = vars(1 + fdep), fn = fn(w ~ 1/n))),
    family = "Weibull",
    data
  )
  expect_equal(result_with_1$mm[[1]]$vars$free, "fdep",
    label = "vars(1 + fdep) should only extract fdep")

  # vars(-1 + fdep) should only extract fdep, not -1
  result_with_minus1 <- dissectFormula(
    formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = vars(-1 + fdep), fn = fn(w ~ 1/n))),
    family = "Weibull",
    data
  )
  expect_equal(result_with_minus1$mm[[1]]$vars$free, "fdep",
    label = "vars(-1 + fdep) should only extract fdep")

  # vars(0 + fdep) should only extract fdep
  result_with_0 <- dissectFormula(
    formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = vars(0 + fdep), fn = fn(w ~ 1/n))),
    family = "Weibull",
    data
  )
  expect_equal(result_with_0$mm[[1]]$vars$free, "fdep",
    label = "vars(0 + fdep) should only extract fdep")

})

test_that("vars() requires formula-style syntax with +", {

  # Comma-separated should throw an error
  expect_error(
    vars(fdep, ipd),
    "must be specified using formula-style",
    label = "vars(x, y) should throw error"
  )

  # Formula-style should work
  expect_no_error(
    vars(fdep + ipd),
    message = "vars(x + y) should work"
  )

  # Single variable should work
  expect_no_error(
    vars(fdep),
    message = "vars(x) should work"
  )

})

test_that("dissectFormula() works with different ways of specifying hm()", {

  expect_error(
    dissectFormula(
      formula(sim.y ~ 1 + hm(name = cname, type = "RE", showFE = F)),
      family="Gaussian",
      data
    ),
    label="hm() -> id() missing"
  )

  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + hm(id = id(cid), name = cname)),
      family="Weibull",
      data
    ),
    message="hm() -> no type or showFE should be possible"
  )

})
