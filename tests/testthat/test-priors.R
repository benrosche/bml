# ================================================================================================ #
# prior(): construction, translation, get_prior, resolution
# ================================================================================================ #

data("coalgov")

test_that("prior() builds and combines", {
  p <- prior(normal(0, 10), class = "b")
  expect_s3_class(p, "bml_prior")
  expect_equal(p$prior, "normal(0, 10)")

  ps <- c(prior(normal(0, 10), class = "b"),
          prior(exponential(1), class = "sd"),
          "b.w.1[1] ~ dnorm(0, 1)")
  expect_s3_class(ps, "bml_prior")
  expect_equal(nrow(ps), 3)
  expect_equal(ps$class[3], "raw")
})

test_that("translate_dist: SD -> precision and friends", {
  expect_equal(translate_dist("normal(0, 10)"), "dnorm(0, 0.01)")
  expect_equal(translate_dist("normal(0, 10)", positive = TRUE), "dnorm(0, 0.01)T(0,)")
  expect_equal(translate_dist("student_t(3, 0, 2)"), "dt(0, 0.25, 3)")
  expect_equal(translate_dist("cauchy(0, 5)"), "dt(0, 0.04, 1)")
  expect_equal(translate_dist("exponential(1)"), "dexp(1)")
  expect_equal(translate_dist("gamma(2, 0.5)"), "dgamma(2, 0.5)")
  expect_equal(translate_dist("uniform(0, 1)"), "dunif(0, 1)")
  expect_equal(translate_dist("lkj(2)"), "__lkj__(2)")
  expect_equal(translate_dist("dnorm(0, 0.5)"), "dnorm(0, 0.5)")  # raw JAGS passthrough
  expect_error(translate_dist("normal(0, -1)"), "sigma")
  expect_error(translate_dist("beta(1, 1)"), "Unsupported")
})

test_that("get_prior enumerates settable priors", {
  gp <- get_prior(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = re(1 + cohesion)),
    data = coalgov, family = gaussian()
  )
  expect_true(all(c("class", "coef", "group", "node", "default") %in% names(gp)))
  expect_true("Intercept" %in% gp$class)
  expect_true(any(gp$class == "b" & gp$coef == "majority"))
  expect_true(any(gp$class == "b" & gp$coef == "A_cohesion" & gp$group == "mm.1"))
  expect_true(any(gp$class == "sd" & gp$coef == "Intercept"))
  expect_true(any(gp$class == "sd" & gp$coef == "cohesion"))
  expect_true(any(gp$class == "sigma"))
})

test_that("get_prior includes w, fn, and cor classes when present", {
  gp <- suppressMessages(suppressWarnings(get_prior(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion),
         w = w(~ ilogit(b0 + b1 * pseat)),
         fn = fn("threshold", c = est(), kappa = 10)) +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
         RE = re(1 + cohesion, cor = TRUE)),
    data = coalgov, family = gaussian()
  )))
  expect_true(any(gp$class == "w" & gp$coef == "b0"))
  expect_true(any(gp$class == "fn" & gp$coef == "c"))
  expect_true(any(gp$class == "cor"))
})

test_that("resolution: narrowing by class/coef/group; unmatched errors; raw passes through", {
  fp <- suppressMessages(dissectFormula(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    "Gaussian", coalgov
  ))
  tab <- build_prior_table(fp, "Gaussian")

  res <- resolve_priors(c(
    prior(normal(0, 5), class = "b"),
    prior(normal(0, 1), class = "b", coef = "majority"),
    prior(exponential(1), class = "sd", group = "mm"),
    "tau ~ dgamma(1, 1)"
  ), tab)

  # coef-specific override wins because it is applied later
  expect_equal(res$overrides[["b[2]"]], "dnorm(0, 1)")
  expect_equal(res$overrides[["b.fn.1[1]"]], "dnorm(0, 0.04)")
  expect_equal(res$overrides[["sigma.mm.1"]], "dexp(1)")
  expect_equal(res$raw, "tau ~ dgamma(1, 1)")

  expect_error(resolve_priors(prior(normal(0, 1), class = "b", coef = "nope"), tab),
               "matches no parameter")
})

test_that("legacy raw list input still passes through", {
  fp <- suppressMessages(dissectFormula(
    dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    "Gaussian", coalgov
  ))
  tab <- build_prior_table(fp, "Gaussian")
  res <- resolve_priors(list("b[1] ~ dnorm(0, 0.5)"), tab)
  expect_equal(res$raw, "b[1] ~ dnorm(0, 0.5)")
})
