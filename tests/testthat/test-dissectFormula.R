data("coalgov")
data <- coalgov

test_that("dissectFormula() works with all model families and with and without intercept", {
  
  # Gaussian
  expect_no_error(
    dissectFormula(
      formula(sim.y ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Gaussian",
      data
    ),
    message="LHS CHECK (Gaussian)"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(sim.st, sim.e) ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Gaussian",
      data
    ),
    label="LHS CHECK (Gaussian)"
  )
  
  
  # Binomial 
  expect_no_error(
    dissectFormula(
      formula(sim.y ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Binomial",
      data
    ),
    message="LHS CHECK (Binomial)"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(sim.st, sim.e) ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Binomial",
      data
    ),
    label="LHS CHECK (Binomial)"
  )

  # Weibull
  expect_no_error(
    dissectFormula(
      formula(Surv(sim.st, sim.e) ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Weibull",
      data
    ),
    message="LHS CHECK (Weibull)"
  )
  
  expect_error(
    dissectFormula(
      formula(sim.y ~ 1 + majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Weibull",
      data
    ),
    label="LHS CHECK (Weibull)"
  )
  
  expect_no_error(
    dissectFormula(
      formula(sim.y ~ majority + mwc + hm(id = cid, name = cname, type = RE, showFE = F) + mm(id = id(pid, gid), vars = vars(ipd), fn = fn(w ~ 1/n, c = TRUE))),
      family="Gaussian",
      data
    ),
    message="INTERCEPT MISSING"
  )
  
  # 2do: add Cox
  
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
  
  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), fn = fn(w ~ 1/exp(ipd), c = FALSE))),
      family="Weibull",
      data
    ),
    label="mm() -> vars missing"
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
  
  expect_equal(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = TRUE), RE = TRUE)),
      family="Weibull",
      data)}[["l1"]][["mmwconstraint"]],
    1,
    label="mm() -> c = TRUE (constrained) should return 1"
  )

  expect_equal(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n, c = FALSE), RE = TRUE)),
      family="Weibull",
      data)}[["l1"]][["mmwconstraint"]],
    2,
    label="mm() -> c = FALSE (unconstrained) should return 2"
  )

  expect_equal(
    {dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + mm(id = id(pid, gid), vars = NULL, fn = fn(w ~ 1/n), RE = TRUE)),
      family="Weibull",
      data)}[["l1"]][["mmwconstraint"]],
    1,
    label="mm() -> default c (TRUE) should return 1"
  )
  
})

test_that("dissectFormula() works with different ways of specifying hm()", {
  
  expect_error(
    dissectFormula(
      formula(sim.y ~ 1 + hm(name = cname, type = RE, showFE = F)),
      family="Gaussian",
      data
    ),
    label="hm() -> id() missing"
  )
  
  expect_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + hm(id = cid, name =, type = RE, showFE = F)),
      family="Weibull",
      data
    ),
    label="hm() -> name improperly specified"
  )
  
  expect_no_error(
    dissectFormula(
      formula(Surv(govdur, earlyterm) ~ 1 + hm(id = cid, name = cname)),
      family="Weibull",
      data
    ),
    message="hm() -> no type or showFE should be possible"
  )
  
 
  
})

