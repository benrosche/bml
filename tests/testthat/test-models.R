# ================================================================================================ #
# Fitted models (JAGS): output structure, methods, equivalences
# ================================================================================================ #

test_that("gaussian sum + RE fit: structure and labels", {
  m <- create_test_model()

  expect_s3_class(m, "bml")
  expect_true(all(c("(Intercept)", "majority", "A_cohesion") %in% m$reg.table$Parameter))
  expect_true(any(grepl("^sd\\(Intercept\\)", m$reg.table$Parameter)))
  expect_equal(m$input$iter, 800)
  expect_equal(m$input$chains, 2)

  s <- summary(m)
  expect_s3_class(s, "bml_summary")
  expect_output(print(s), "Coefficients \\(class \"b\"\\)")
  expect_output(print(s), "Variance parameters")

  td <- tidy(m)
  expect_true(all(c("term", "estimate", "std.error", "conf.low", "conf.high") %in% names(td)))
  expect_true("(Intercept)" %in% td$term)

  g <- glance(m)
  expect_equal(g$nobs, m$input$n.main)
})

test_that("posterior methods: fixef, ranef, as_draws, fitted, predict, loo", {
  m <- create_test_model()

  fx <- fixef(m)
  expect_true(all(c("(Intercept)", "majority", "A_cohesion") %in% rownames(fx)))
  expect_equal(colnames(fx), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))

  rf <- ranef(m)
  expect_length(rf$mm[[1]], m$input$n.umm)

  dd <- as_draws_df(m)
  expect_true(posterior::ndraws(dd) > 0)

  ft <- fitted(m)
  expect_equal(nrow(ft), m$input$n.main)

  pp <- posterior_predict(m, ndraws = 10)
  expect_equal(dim(pp), c(10, m$input$n.main))

  pr <- predict(m)
  expect_equal(nrow(pr), m$input$n.main)

  ll <- log_lik(m)
  expect_equal(ncol(ll), m$input$n.main)

  l <- suppressWarnings(loo(m))
  expect_s3_class(l, "loo")
  wc <- suppressWarnings(waic(m))
  expect_s3_class(wc, "waic")

  expect_error(posterior_predict(m, newdata = coalgov), "not supported")
})

test_that("varDecomp still discovers the components", {
  m <- create_test_model()
  vd <- varDecomp(m)
  expect_s3_class(vd, "bml_varDecomp")
  expect_true(any(grepl("Residual", vd$Component)))
  expect_true(any(grepl("MM", vd$Component)))
})

test_that("fixed weights + no RE reproduces lm() on the precomputed feature", {
  skip_if_no_jags()

  ax <- coalgov |>
    dplyr::group_by(gid) |>
    dplyr::summarise(Ax = mean(cohesion), dur = dplyr::first(dur_wkb),
                     majority = dplyr::first(majority))
  lmfit <- stats::lm(dur ~ majority + Ax, data = ax)

  m <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE),
    data = coalgov, family = gaussian(),
    iter = 3000, warmup = 1000, chains = 2, seed = 2
  )))

  b_bml <- fixef(m)["A_cohesion", "Estimate"]
  b_lm <- stats::coef(lmfit)["Ax"]
  se_lm <- summary(lmfit)$coefficients["Ax", 2]
  expect_lt(abs(b_bml - b_lm), 3 * se_lm)
})

test_that("fn(\"var\") reproduces lm() on the precomputed weighted variance", {
  skip_if_no_jags()

  vx <- coalgov |>
    dplyr::group_by(gid) |>
    dplyr::summarise(Vx = mean((cohesion - mean(cohesion))^2),
                     dur = dplyr::first(dur_wkb), majority = dplyr::first(majority))
  lmv <- stats::lm(dur ~ majority + Vx, data = vx)

  m <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), fn = fn("var")),
    data = coalgov, family = gaussian(),
    iter = 3000, warmup = 1000, chains = 2, seed = 3
  )))

  expect_true("V_cohesion" %in% rownames(fixef(m)))
  expect_lt(abs(fixef(m)["V_cohesion", "Estimate"] - stats::coef(lmv)["Vx"]),
            3 * summary(lmv)$coefficients["Vx", 2])
})

test_that("estimated weights fit with class-w prior; weight matrix returned", {
  skip_if_no_jags()

  m <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion),
         w = w(~ ilogit(b0 + b1 * pseat), scale = TRUE), RE = FALSE),
    data = coalgov, family = gaussian(),
    prior = prior(normal(0, 1), class = "w"),
    monitor = "parameters",
    iter = 1000, warmup = 400, chains = 2, seed = 4
  )))

  expect_true(any(grepl("^w\\[", m$reg.table$Parameter)))
  expect_equal(nrow(m$w[[1]]), m$input$n.main)
  # weights are normalized: rows sum to ~1
  rs <- rowSums(m$w[[1]], na.rm = TRUE)
  expect_true(all(abs(rs - 1) < 1e-6))
})

test_that("random slopes fit (independent and correlated)", {
  skip_if_no_jags()

  m <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = re(1 + cohesion)),
    data = coalgov, family = gaussian(),
    monitor = "random_effects",
    iter = 800, warmup = 300, chains = 2, seed = 5
  )))
  expect_true(any(grepl("sd\\(cohesion\\)", m$reg.table$Parameter)))
  expect_true(is.list(ranef(m)$mm[[1]]))

  mc <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
         RE = re(1 + cohesion, cor = TRUE)),
    data = coalgov, family = gaussian(),
    monitor = "random_effects",
    iter = 800, warmup = 300, chains = 2, seed = 6
  )))
  expect_true(any(grepl("^cor\\(Intercept,cohesion\\)", mc$reg.table$Parameter)))
})

test_that("cross-level interaction matches lm()", {
  skip_if_no_jags()

  ax <- coalgov |>
    dplyr::group_by(gid) |>
    dplyr::summarise(Ax = mean(cohesion), dur = dplyr::first(dur_wkb),
                     majority = dplyr::first(majority))
  lmx <- stats::lm(dur ~ majority * Ax, data = ax)

  m <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority + Ax:majority +
      mm(name = Ax, id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE),
    data = coalgov, family = gaussian(),
    iter = 3000, warmup = 1000, chains = 2, seed = 7
  )))

  expect_true("Ax:majority" %in% rownames(fixef(m)))
  expect_lt(abs(fixef(m)["Ax:majority", "Estimate"] - stats::coef(lmx)["majority:Ax"]),
            4 * summary(lmx)$coefficients["majority:Ax", 2])
})

test_that("bernoulli family fits and supports waic", {
  skip_if_no_jags()

  m <- suppressWarnings(suppressMessages(bml(
    event_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n), RE = FALSE),
    data = coalgov, family = bernoulli(),
    monitor = "log_lik",
    iter = 800, warmup = 300, chains = 2, seed = 8
  )))
  expect_true("majority" %in% rownames(fixef(m)))
  expect_s3_class(suppressWarnings(waic(m)), "waic")
})

test_that("time-indexed AR random walk fits", {
  skip_if_no_jags()

  m <- suppressWarnings(suppressMessages(bml(
    dur_wkb ~ 1 + majority +
      mm(id = id(pid, gid), vars = vars(cohesion), w = w(~ 1/n),
         RE = re(1, ar = gid)),
    data = coalgov, family = gaussian(),
    monitor = "random_effects",
    iter = 800, warmup = 300, chains = 2, seed = 12
  )))
  expect_true(any(grepl("^sd\\(Intercept\\)", m$reg.table$Parameter)))
  # AR REs come back as a member x walk-step matrix
  expect_true(is.matrix(ranef(m)$mm[[1]]))
})
