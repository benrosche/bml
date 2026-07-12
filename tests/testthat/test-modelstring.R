# ================================================================================================ #
# Generated JAGS code (run = FALSE): fn types, DSL, re/fe, interactions, priors
# ================================================================================================ #

data("coalgov")

code_of <- function(formula, prior = NULL, family = stats::gaussian(), monitor = TRUE) {
  out <- suppressWarnings(suppressMessages(
    bml(formula, data = coalgov, family = family, prior = prior,
        run = FALSE, monitor = monitor)
  ))
  out
}

test_that("sum block: precomputed feature, coefficient outside, log_lik + mu monitored", {
  out <- code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE)
  )
  expect_match(out$modelstring, "b.fn.1\\[1\\] \\* F.1\\[j\\]")
  expect_match(out$modelstring, "log_lik\\[j\\] <- logdensity.norm")
  expect_true(all(c("F.1", "w.1") %in% names(out$jags.data)))
  expect_true(all(c("b.fn.1", "mu", "log_lik", "pred") %in% out$jags.params))
})

test_that("emergent types emit the expected reductions", {
  out <- code_of(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion),
         w = w(~ ilogit(b0 + b1 * pseat)), fn = fn("var")) +
      mm(id = id(pid, gid), vars = vars(cohesion),
         w = w(~ ilogit(b0 + b1 * pseat)), fn = fn("smax", kappa = 5)) +
      mm(id = id(pid, gid), vars = vars(cohesion),
         w = w(~ ilogit(b0 + b1 * pseat)), fn = fn("threshold", c = est(), kappa = 10))
  )
  ms <- out$modelstring
  # var: two-pass accumulator
  expect_match(ms, "A.1\\[j\\] <- inprod")
  expect_match(ms, "pow\\(X.mm.1\\[i,1\\] - A.1\\[grp.mm.1\\[i\\]\\], 2\\)")
  # smax: log-sum-exp
  expect_match(ms, "\\(1 / 5\\) \\* log\\(sum\\(e.2\\[")
  # threshold with estimated cutpoint
  expect_match(ms, "ilogit\\(10 \\* \\(X.mm.3\\[i,1\\] - fn.c.3\\)\\)")
  expect_match(ms, "fn.c.3 ~ dunif\\(")
  expect_true("fn.c.3" %in% out$jags.params)
})

test_that("hhi/effn/entropy read the weights alone", {
  out <- code_of(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), w = w(~ ilogit(b0 + b1 * pseat)), fn = fn("hhi"), RE = FALSE)
  )
  expect_match(out$modelstring, "F.1\\[j\\] <- inprod\\(w.1\\[")
})

test_that("cov emits the two-attribute reduction", {
  out <- code_of(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion, pseat),
         w = w(~ ilogit(b0 + b1 * finance)), fn = fn("cov"))
  )
  expect_match(out$modelstring, "A1.1\\[j\\]")
  expect_match(out$modelstring, "A2.1\\[j\\]")
  expect_match(out$modelstring,
               "\\(X.mm.1\\[i,1\\] - A1.1\\[grp.mm.1\\[i\\]\\]\\) \\* \\(X.mm.1\\[i,2\\] - A2.1\\[grp.mm.1\\[i\\]\\]\\)")
})

test_that("DSL compiles to one accumulator pass per E() node", {
  out <- code_of(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion),
         w = w(~ ilogit(b0 + b1 * pseat)),
         fn = fn(~ E((cohesion - E(cohesion))^2)))
  )
  ms <- out$modelstring
  expect_match(ms, "E1.1\\[j\\] <- sum\\(t1.1\\[")            # inner mean
  expect_match(ms, "E1.1\\[grp.mm.1\\[i\\]\\]")               # referenced at member level
  expect_match(ms, "E2.1\\[j\\] <- sum\\(t2.1\\[")            # outer reduction
  expect_match(ms, "F.1\\[j\\] <- E2.1\\[j\\]")
})

test_that("re slopes: independent by default, correlated 2x2 on opt-in", {
  out_i <- code_of(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = re(1 + cohesion))
  )
  expect_match(out_i$modelstring, "re.mm.1.s1\\[i\\] ~ dnorm\\(0, tau.mm.1.s1\\)")
  expect_no_match(out_i$modelstring, "rho.mm.1")

  out_c <- code_of(
    dur_wkb ~ 1 +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
         RE = re(1 + cohesion, cor = TRUE))
  )
  expect_match(out_c$modelstring, "dmnorm.vcov")
  expect_match(out_c$modelstring, "rho.mm.1 <- 2 \\* z.rho.mm.1 - 1")
  expect_true("rho.mm.1" %in% out_c$jags.params)
})

test_that("hm: re() and fe() grammar, reference coding for fe", {
  out <- code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
      hm(id = id(cid), FE = fe(1, showFE = TRUE))
  )
  expect_match(out$modelstring, "alpha.hm.1\\[1\\] <- 0")
  expect_match(out$modelstring, "alpha.hm.1\\[c\\] ~ dnorm\\(0, 0.0001\\)")
  expect_true("alpha.hm.1" %in% out$jags.params)

  out2 <- suppressWarnings(code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
      hm(id = id(cid), RE = re(1 + majority))
  ))
  expect_match(out2$modelstring, "re.hm.1.s1\\[hmid\\[j\\]\\] \\* X.re.hm.1\\[j,1\\]")
})

test_that("interactions enter mu with b.int", {
  out <- code_of(
    dur_wkb ~ 1 + majority + Ax:majority +
      mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE)
  )
  expect_match(out$modelstring, "b.int\\[1\\] \\* F.1\\[j\\] \\* X.main\\[j,2\\]")

  out2 <- code_of(
    dur_wkb ~ 1 + Ax:Vx +
      mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
      mm(name = Vx, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("var"))
  )
  expect_match(out2$modelstring, "b.int\\[1\\] \\* F.1\\[j\\] \\* F.2\\[j\\]")
})

test_that("prior() overrides land on the right nodes; sd overrides reparameterize", {
  out <- code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    prior = c(prior(normal(0, 5), class = "b"),
              prior(student_t(3, 0, 2), class = "Intercept"),
              prior(exponential(1), class = "sd"),
              prior(normal(0, 10), class = "sigma"))
  )
  ms <- out$modelstring
  expect_match(ms, "b\\[1\\] ~ dt\\(0, 0.25, 3\\)")
  expect_match(ms, "b\\[2\\] ~ dnorm\\(0, 0.04\\)")
  expect_match(ms, "b.fn.1\\[1\\] ~ dnorm\\(0, 0.04\\)")
  expect_match(ms, "sigma.mm.1 ~ dexp\\(1\\)")
  expect_match(ms, "tau.mm.1 <- pow\\(sigma.mm.1, -2\\)")
  expect_match(ms, "sigma ~ dnorm\\(0, 0.01\\)T\\(0,\\)")
  expect_match(ms, "tau <- pow\\(sigma, -2\\)")
})

test_that("raw JAGS strings replace prior lines (escape hatch)", {
  out <- code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    prior = c(prior(normal(0, 5), class = "b"), "b.fn.1[1] ~ dnorm(0, 0.5)")
  )
  expect_match(out$modelstring, "b.fn.1\\[1\\] ~ dnorm\\(0, 0.5\\)", fixed = FALSE)
})

test_that("monitor = FALSE drops pred/mu/log_lik", {
  out <- code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    monitor = FALSE
  )
  expect_no_match(out$modelstring, "log_lik")
  expect_false(any(c("mu", "pred", "log_lik") %in% out$jags.params))
})

test_that("removed legacy bml() arguments give a migration error", {
  expect_error(
    bml(dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n)),
        data = coalgov, n.iter = 100),
    "Removed argument"
  )
  expect_error(
    bml(dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n)),
        data = coalgov, priors = list("b[1] ~ dnorm(0, 1)")),
    "Removed argument"
  )
})

test_that("make_jagscode / make_jagsdata return the pieces", {
  code <- make_jagscode(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    data = coalgov, family = gaussian()
  )
  expect_s3_class(code, "bml_jagscode")
  expect_match(as.character(code), "model \\{")

  dat <- make_jagsdata(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = TRUE),
    data = coalgov, family = gaussian()
  )
  expect_true(all(c("Y", "X.main", "F.1") %in% names(dat)))
})

test_that("AR random walks: participation-order vs time-indexed (gap-scaled)", {
  out_o <- code_of(
    dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
                     RE = re(1, ar = TRUE))
  )
  expect_match(out_o$modelstring,
               "re.mm.1\\[i,t\\] ~ dnorm\\(re.mm.1\\[i,t-1\\], tau.mm.1\\)")
  expect_no_match(out_o$modelstring, "gap.mm.1")

  out_t <- code_of(
    dur_wkb ~ 1 + mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
                     RE = re(1, ar = gid))
  )
  expect_match(out_t$modelstring,
               "re.mm.1\\[i,t\\] ~ dnorm\\(re.mm.1\\[i,t-1\\], tau.mm.1 / gap.mm.1\\[i,t\\]\\)")
  expect_true("gap.mm.1" %in% names(out_t$jags.data))
  gm <- out_t$jags.data[["gap.mm.1"]]
  expect_true(all(gm > 0))
  expect_true(all(gm[, 1] == 1))  # first walk step has no gap

  # hm level
  out_h <- code_of(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE) +
      hm(id = id(cid), RE = re(1, ar = gid))
  )
  expect_match(out_h$modelstring,
               "re.hm.1\\[c,t\\] ~ dnorm\\(re.hm.1\\[c,t-1\\], tau.hm.1 / gap.hm.1\\[c,t\\]\\)")
  expect_true("gap.hm.1" %in% names(out_h$jags.data))
})
